// Fault network: tip resolution, throw profiles, kinematics, hanging-wall fields, restoration.
#include "stratonix/engines/faults.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "stratonix/core/error.hpp"
#include "stratonix/core/roots.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/structure.hpp"
#include "stratonix/engines/wells.hpp"

// Plugin anchors called below (C linkage, defined by STRATONIX_REGISTER in each plugin file).
STRATONIX_DECLARE_PLUGIN_ANCHOR(kin_dip_slip);
STRATONIX_DECLARE_PLUGIN_ANCHOR(kin_vertical);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_rigid);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_layer_rigid);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_decay);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_linear);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_vertical_shear);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_inclined_shear);

namespace stratonix {

constexpr double kInf = std::numeric_limits<double>::infinity();

// ---- profiles ------------------------------------------------------------------------------
ProfileType profile_from(const std::string& s) {
  if (s == "sine") return ProfileType::Sine;
  if (s == "elliptical") return ProfileType::Elliptical;
  if (s == "linear") return ProfileType::Linear;
  return ProfileType::FlatTopped;
}
double profile_value(ProfileType p, double taper, double t) {
  if (t <= 0 || t >= 1) return 0;
  switch (p) {
    case ProfileType::Sine: return std::sin(kPi * t);
    case ProfileType::Elliptical: return std::sqrt(std::max(0.0, 1 - (2 * t - 1) * (2 * t - 1)));
    case ProfileType::Linear: return 1 - std::fabs(2 * t - 1);
    case ProfileType::FlatTopped:
      if (taper <= 0) return 1;
      if (t < taper) return 0.5 * (1 - std::cos(kPi * t / taper));
      if (t > 1 - taper) return 0.5 * (1 - std::cos(kPi * (1 - t) / taper));
      return 1;
  }
  return 0;
}

// ---- curved faults (STX-TN-05 section 2) ---------------------------------------------------
double CurvedFault::zf(double u) const {
  if (shape == Shape::Listric) return zd - (zd - z0) * std::exp(-u * tan0 / (zd - z0));
  const std::size_t n = pu.size();
  std::size_t k;
  if (u <= pu[0]) k = 0;
  else if (u >= pu[n - 1]) k = n - 2;
  else k = static_cast<std::size_t>(std::upper_bound(pu.begin(), pu.end(), u) - pu.begin()) - 1;
  return pz[k] + (pz[k + 1] - pz[k]) * (u - pu[k]) / (pu[k + 1] - pu[k]);
}
double CurvedFault::dzdu(double u) const {
  if (shape == Shape::Listric) return tan0 * std::exp(-u * tan0 / (zd - z0));
  const std::size_t n = pu.size();
  std::size_t k;
  if (u <= pu[0]) k = 0;
  else if (u >= pu[n - 1]) k = n - 2;
  else k = static_cast<std::size_t>(std::upper_bound(pu.begin(), pu.end(), u) - pu.begin()) - 1;
  return (pz[k + 1] - pz[k]) / (pu[k + 1] - pu[k]);
}
double CurvedFault::u_at(double z) const {
  if (shape == Shape::Listric) {
    if (z >= zd) return std::numeric_limits<double>::infinity();
    return -(zd - z0) / tan0 * std::log((zd - z) / (zd - z0));
  }
  // the first segment dips (validated); flats are skipped by taking the first point reaching z
  const std::size_t n = pu.size();
  if (z <= pz[0]) return pu[0] + (z - pz[0]) * (pu[1] - pu[0]) / (pz[1] - pz[0]);
  for (std::size_t k = 0; k + 1 < n; ++k)
    if (z <= pz[k + 1] && pz[k + 1] > pz[k]) return pu[k] + (z - pz[k]) * (pu[k + 1] - pu[k]) / (pz[k + 1] - pz[k]);
  double dz = pz[n - 1] - pz[n - 2];
  if (dz <= 0) return std::numeric_limits<double>::infinity();
  return pu[n - 1] + (z - pz[n - 1]) * (pu[n - 1] - pu[n - 2]) / dz;
}
double CurvedFault::G(double xi) const {
  if (ta == 0) return zf(xi);
  // The shear line meets the fault once (rule F12): phi decreases strictly with z.
  auto phi = [&](double z) { return zf(xi + s * z * ta) - z; };
  double z = expand_and_solve(phi, zf(xi), 50.0, 1.0e7, 1e-10);
  return z;
}
void CurvedFault::forward(double u, double z, double h, double& u1, double& z1) const {
  double xi = u - s * z * ta, xi1 = xi + h;
  z1 = z + G(xi1) - G(xi);
  u1 = xi1 + s * z1 * ta;
}
void CurvedFault::inverse(double u1, double z1, double h, double& u, double& z) const {
  double xi1 = u1 - s * z1 * ta, xi = xi1 - h;
  z = z1 - G(xi1) + G(xi);
  u = xi + s * z * ta;
}

// ---- kinematics and hanging-wall plugins ---------------------------------------------------
namespace {
class DipSlip final : public Kinematics {
 public:
  Vec3 slip_vector(const Vec3& dd, double dip) const override { return dd / std::sin(dip); }
  std::string name() const override { return "dip_slip"; }
};
class VerticalOffset final : public Kinematics {
 public:
  Vec3 slip_vector(const Vec3&, double) const override { return {0, 0, -1}; }
  std::string name() const override { return "vertical_offset"; }
};
class RigidHW final : public HangingWallField {
 public:
  double scale(double, double, double) const override { return 1.0; }
  std::string name() const override { return "rigid"; }
};
class LayerRigidHW final : public HangingWallField {
 public:
  double scale(double, double, double) const override { return 1.0; }
  std::string name() const override { return "layer_rigid"; }
};
class DecayHW final : public HangingWallField {
 public:
  double scale(double sigma, double, double infl) const override {
    return infl > 0 ? std::max(0.0, 1 - std::fabs(sigma) / infl) : 1.0;
  }
  std::string name() const override { return "decay"; }
};
/// Vertical and inclined shear: the displacement is the shear construction of the curved fault
/// (CurvedFault), not a scale; scale() is not used for them.
class ShearHW final : public HangingWallField {
 public:
  explicit ShearHW(std::string n) : name_(std::move(n)) {}
  double scale(double, double, double) const override { return 1.0; }
  std::string name() const override { return name_; }
 private:
  std::string name_;
};
class LinearBetweenHW final : public HangingWallField {
 public:
  double scale(double sigma, double next, double) const override {
    if (!std::isfinite(next) || next <= 0) return 1.0;
    return std::max(0.0, 1 - std::fabs(sigma) / next);
  }
  bool needs_next_fault() const override { return true; }
  std::string name() const override { return "linear_between_faults"; }
};
}  // namespace

STRATONIX_REGISTER(kin_dip_slip, registry::Kind::Kinematics, "dip_slip", "1.0.0",
                   "Hanging wall moves down the fault dip; vertical component equals the throw",
                   [] { return std::make_unique<DipSlip>(); })
STRATONIX_REGISTER(kin_vertical, registry::Kind::Kinematics, "vertical_offset", "1.0.0",
                   "Hanging wall moves vertically by the throw (prototype mode)",
                   [] { return std::make_unique<VerticalOffset>(); })
STRATONIX_REGISTER(hw_rigid, registry::Kind::HangingWall, "rigid", "1.0.0",
                   "Full displacement throughout the moving block", [] { return std::make_unique<RigidHW>(); })
STRATONIX_REGISTER(hw_layer_rigid, registry::Kind::HangingWall, "layer_rigid", "1.0.0",
                   "Each layer moves as a rigid slab by the throw where that layer meets the fault",
                   [] { return std::make_unique<LayerRigidHW>(); })
STRATONIX_REGISTER(hw_decay, registry::Kind::HangingWall, "decay", "1.0.0",
                   "Displacement decays linearly to zero at an influence distance",
                   [] { return std::make_unique<DecayHW>(); })
STRATONIX_REGISTER(hw_linear, registry::Kind::HangingWall, "linear_between_faults", "1.0.0",
                   "Displacement decays linearly to zero at the next fault plane",
                   [] { return std::make_unique<LinearBetweenHW>(); })
STRATONIX_REGISTER(hw_vertical_shear, registry::Kind::HangingWall, "vertical_shear", "1.0.0",
                   "Curved faults: hanging wall collapses onto the fault by vertical simple shear (Chevron construction)",
                   [] { return std::make_unique<ShearHW>("vertical_shear"); })
STRATONIX_REGISTER(hw_inclined_shear, registry::Kind::HangingWall, "inclined_shear", "1.0.0",
                   "Curved faults: hanging wall collapses onto the fault by simple shear on inclined lines",
                   [] { return std::make_unique<ShearHW>("inclined_shear"); })

void register_fault_plugins() {
  STRATONIX_PLUGIN_ANCHOR(kin_dip_slip);
  STRATONIX_PLUGIN_ANCHOR(kin_vertical);
  STRATONIX_PLUGIN_ANCHOR(hw_rigid);
  STRATONIX_PLUGIN_ANCHOR(hw_layer_rigid);
  STRATONIX_PLUGIN_ANCHOR(hw_decay);
  STRATONIX_PLUGIN_ANCHOR(hw_linear);
  STRATONIX_PLUGIN_ANCHOR(hw_vertical_shear);
  STRATONIX_PLUGIN_ANCHOR(hw_inclined_shear);
}

// ---- fault ---------------------------------------------------------------------------------
void Fault::tips(double a, double& bu, double& bl) const {
  double L2 = 0.5 * length;
  if (a_tab.empty() || std::fabs(a) > L2) { bu = bl = 0; return; }
  double t = (a + L2) / length * (a_tab.size() - 1);
  std::size_t i = std::min<std::size_t>(static_cast<std::size_t>(t), a_tab.size() - 2);
  double f = t - i;
  bu = bu_tab[i] + f * (bu_tab[i + 1] - bu_tab[i]);
  bl = bl_tab[i] + f * (bl_tab[i + 1] - bl_tab[i]);
  if (elliptical) {
    double e = std::sqrt(std::max(0.0, 1 - (a / L2) * (a / L2)));
    bool ub = !upper_free(), lb = !lower_free();
    if (ub && lb) {
      double c = 0.5 * (bu + bl), h = 0.5 * (bl - bu) * e;
      bu = c - h;
      bl = c + h;
    } else if (lb) {
      bl = bu + (bl - bu) * e;
    } else if (ub) {
      bu = bl - (bl - bu) * e;
    }
  }
}

double Fault::throw_at(double a, double b) const {
  double L2 = 0.5 * length;
  if (std::fabs(a) >= L2) return 0;
  double bu, bl;
  tips(a, bu, bl);
  if (!(bl > bu)) return 0;
  double eta = (b - bu) / (bl - bu);
  if (eta <= 0 && !upper_free()) return 0;
  if (eta >= 1 && !lower_free()) return 0;
  if (upper_free() && eta < 0.5) eta = 0.5;
  if (lower_free() && eta > 0.5) eta = 0.5;
  double xi = (a + L2) / length;
  return tmax * profile_value(pdip, taper_dip, eta) * profile_value(pstrike, taper_strike, xi);
}

double Fault::heave_at(double a) const {
  if (!curved) return 0;
  double L2 = 0.5 * length;
  if (std::fabs(a) >= L2) return 0;
  return curved->heave * profile_value(pstrike, taper_strike, (a + L2) / length);
}

double Fault::layer_b(double s, double a) const {
  const std::size_t na = la_tab.size();
  if (na < 2) return 0;
  auto column = [&](std::size_t i) {
    const double* ls = ls_tab.data() + l_start[i];
    const double* lb = lb_tab.data() + l_start[i];
    std::size_t n = l_start[i + 1] - l_start[i];
    if (s <= ls[0]) return lb[0];
    if (s >= ls[n - 1]) return lb[n - 1];
    std::size_t j = static_cast<std::size_t>(std::upper_bound(ls, ls + n, s) - ls);
    double w = (s - ls[j - 1]) / (ls[j] - ls[j - 1]);
    return lb[j - 1] + w * (lb[j] - lb[j - 1]);
  };
  double t = std::clamp((a - la_tab.front()) / (la_tab.back() - la_tab.front()), 0.0, 1.0) * static_cast<double>(na - 1);
  std::size_t i0 = std::min(static_cast<std::size_t>(t), na - 2);
  double fa = t - static_cast<double>(i0);
  double b0 = column(i0), b1 = column(i0 + 1);
  return b0 + fa * (b1 - b0);
}

double Fault::extent_coord(double a, double b) const {
  double L2 = 0.5 * length;
  double ra = std::fabs(a) / L2;
  if (ra >= 1) return ra;
  double bu, bl;
  tips(a, bu, bl);
  if (!(bl > bu)) return 1.0 + 1e-9 + ra;
  double eta = (b - bu) / (bl - bu);
  if (upper_free() && eta < 0.5) eta = 0.5;
  if (lower_free() && eta > 0.5) eta = 0.5;
  return std::max(ra, std::fabs(2 * eta - 1));
}

// ---- network -------------------------------------------------------------------------------
static TipRule tip_rule(const desc::TipDesc& t, const desc::Description& d) {
  TipRule r;
  if (t.rule == "horizon") {
    r.kind = TipRule::Kind::Horizon;
    r.horizon = d.horizon_index(t.horizon);
    r.offset = t.offset.si();
  } else if (t.rule == "within_unit") {
    r.kind = TipRule::Kind::WithinUnit;
    auto& u = d.units[static_cast<std::size_t>(d.unit_index(t.unit))];
    r.top = d.horizon_index(u.top);
    r.base = d.horizon_index(u.base);
    r.fraction = t.fraction.si();
  } else if (t.rule == "depth") {
    r.kind = TipRule::Kind::Depth;
    r.depth = t.depth.si();
  }
  r.unfaulted = t.stratigraphy == "unfaulted";
  return r;
}

/// Elevation of the tip rule level at (x, y) on the restored structure.
static double rule_level(const TipRule& r, const Structure& st, double x, double y, const Domain& dom, bool upper) {
  switch (r.kind) {
    case TipRule::Kind::Horizon: return st.z(static_cast<std::size_t>(r.horizon), x, y) - r.offset;
    case TipRule::Kind::WithinUnit: {
      double zt = st.z(static_cast<std::size_t>(r.top), x, y), zb = st.z(static_cast<std::size_t>(r.base), x, y);
      return zt + r.fraction * (zb - zt);
    }
    case TipRule::Kind::Depth: return -r.depth;
    case TipRule::Kind::Free: return upper ? dom.z_top : dom.z_base;
  }
  return kNaN;
}

FaultNetwork::FaultNetwork(const desc::Description& d, const Structure& st, BuildLog& log) {
  auto& reg = registry::Registry::global();
  for (std::size_t k = 0; k < d.faults.size(); ++k) {
    const auto& fd = d.faults[k];
    Fault f;
    f.name = fd.name;
    f.index = static_cast<int>(k);
    f.age = fd.age.si();
    f.dip = fd.dip.present ? fd.dip.si() : 0.0;
    if (fd.curved()) {
      auto c = std::make_shared<CurvedFault>();
      if (fd.geom_type == "listric") {
        c->shape = CurvedFault::Shape::Listric;
        c->z0 = fd.reference_depth.si();
        c->tan0 = std::tan(fd.dip.si());
        c->zd = fd.detachment_depth.si();
      } else {
        c->shape = CurvedFault::Shape::RampFlat;
        for (std::size_t i = 0; i < fd.profile.rows(); ++i) {
          c->pu.push_back(fd.profile.si(2 * i));
          c->pz.push_back(fd.profile.si(2 * i + 1));
        }
        f.dip = std::atan2(c->pz[1] - c->pz[0], c->pu[1] - c->pu[0]);
      }
      c->s = fd.shear_sense == "synthetic" ? 1.0 : -1.0;
      c->ta = fd.hw_type == "inclined_shear" ? std::tan(fd.shear_angle.si()) : 0.0;
      c->heave = fd.heave.si();
      f.curved = c;
    }
    if (fd.trace.present) {
      Vec3 a{fd.trace.si(0), fd.trace.si(1), -fd.trace_depth.si()};
      Vec3 b{fd.trace.si(2), fd.trace.si(3), -fd.trace_depth.si()};
      Vec3 sv = normalized(b - a);
      double strike = std::atan2(sv.x, sv.y);
      // dip_side right: dips to the right of the trace direction (right-hand rule), left: opposite
      if (fd.dip_side == "left") strike += kPi;
      f.strike = strike;
      f.p0 = (a + b) * 0.5;
    } else {
      f.strike = fd.strike.si();
      double z0 = fd.reference_depth.present ? fd.reference_depth.si() : (f.curved ? f.curved->zf(0) : 0.0);
      f.p0 = {fd.reference_point.si(0), fd.reference_point.si(1), -z0};
    }
    f.s = azimuth_vector(f.strike);
    f.dh = azimuth_vector(f.strike + kPi / 2);
    f.dd = f.dh * std::cos(f.dip) + Vec3{0, 0, -std::sin(f.dip)};
    f.n = f.dh * std::sin(f.dip) + Vec3{0, 0, std::cos(f.dip)};
    f.length = fd.length.si();
    f.elliptical = fd.tip_line_shape == "elliptical";
    f.upper = tip_rule(fd.upper, d);
    f.lower = tip_rule(fd.lower, d);
    f.tmax = f.curved ? f.curved->heave : fd.throw_max.si();
    f.pdip = profile_from(fd.profile_dip.type);
    f.pstrike = profile_from(fd.profile_strike.type);
    f.taper_dip = fd.profile_dip.taper.present ? fd.profile_dip.taper.si() : 0.25;
    f.taper_strike = fd.profile_strike.taper.present ? fd.profile_strike.taper.si() : 0.25;
    f.kin = reg.create<Kinematics>(registry::Kind::Kinematics, fd.slip_mode);
    f.slip = f.kin->slip_vector(f.dd, f.dip);
    f.hw = reg.create<HangingWallField>(registry::Kind::HangingWall, fd.hw_type);
    f.hw_name = f.hw->name();
    f.layered = f.hw->name() == "layer_rigid";
    f.influence = fd.influence_distance.present ? fd.influence_distance.si() : 0.0;
    if (!fd.abuts_fault.empty()) {
      f.abuts = d.fault_index(fd.abuts_fault);
      f.abut_hw_side = fd.abuts_side == "hanging_wall";
    }
    faults_.push_back(std::move(f));
  }
  st_ = &st;
  {
    const Domain dom = d.domain();
    z_top_ = dom.z_top;
    z_base_ = dom.z_base;
    std::size_t nl = 0;
    for (auto& f : faults_) nl += f.layered ? 1 : 0;
    if (nl != 0 && nl != faults_.size())
      throw Error(ErrorCode::Build, "hanging wall 'layer_rigid' must be used by all faults of a model or by none", "faults");
    layered_ = nl != 0 && nl == faults_.size();
    vertical_ = true;
    for (auto& f : faults_) vertical_ = vertical_ && std::fabs(f.slip.x) < 1e-15 && std::fabs(f.slip.y) < 1e-15;
  }
  // Ages: youngest first (stable for equal ages: later in the list is treated as younger).
  order_.resize(faults_.size());
  std::iota(order_.begin(), order_.end(), 0);
  std::stable_sort(order_.begin(), order_.end(), [&](int a, int b) {
    if (faults_[a].age != faults_[b].age) return faults_[a].age > faults_[b].age;
    return a > b;
  });
  rank_.assign(faults_.size(), 0);
  for (std::size_t i = 0; i < order_.size(); ++i) rank_[static_cast<std::size_t>(order_[i])] = static_cast<int>(i);
  // Frame rule: fixed point (only needed when there are faults).
  if (faults_.empty()) return;
  if (d.frame_point.present) {
    fixed_ = {d.frame_point.si(0), d.frame_point.si(1), d.frame_point.si(2)};
  } else {
    std::string wn = !d.frame_well.empty() ? d.frame_well : d.ref_well;
    const desc::WellDesc* w = wn.empty() ? nullptr : d.well(wn);
    if (!w) throw Error(ErrorCode::Build, "frame rule: no fixed well or point", "frame");
    // Fix the block at the deepest top of the well (on the well path).
    double zdeep = -d.depth_range.si(0);
    for (auto& [h, p] : w->tops) zdeep = std::min(zdeep, -p.si());
    fixed_ = make_well_path(*w, d).at_elevation(zdeep);
  }
  fixed_block_ = 0;
  for (int f : order_) {
    Vec3 q = fixed_;  // the fixed point never moves
    if (faults_[static_cast<std::size_t>(f)].sigma(q) > 0) fixed_block_ |= (BlockCode(1) << f);
  }
  for (auto& f : faults_)
    if (f.curved && ((fixed_block_ >> f.index) & 1u))
      throw Error(ErrorCode::Build,
                  "fault " + f.name + " (rule F11): the frame's fixed point lies in the hanging wall of a curved fault; the "
                  "shear construction moves the hanging wall over a fixed footwall, so fix a point or well in the footwall",
                  "frame.fixed");
  resolve_tips(d, st, log);
}

/// Elevation of horizon h at (x, y) in the frame of fault `rank` (stratigraphy displaced by the
/// faults older than it, which are already resolved).
double FaultNetwork::level_in_frame(const Structure& st, int h, double x, double y, int f, const Domain& dom) const {
  const std::size_t r_older = static_cast<std::size_t>(rank_[static_cast<std::size_t>(f)]) + 1;
  auto s_at = [&](double z) {
    Vec3 q{x, y, z};
    if (layered_) {
      // material on horizon h is at level h: restore it through the older faults directly
      q = chain(q, r_older, order_.size(), static_cast<double>(h), nullptr, nullptr, 0);
      return q.z - st.z(static_cast<std::size_t>(h), q.x, q.y);
    }
    for (std::size_t k = r_older; k < order_.size(); ++k) {
      int fi = order_[k];
      const Fault& F = faults_[static_cast<std::size_t>(fi)];
      bool side = F.sigma(q) > 0;
      bool fixed_side = ((fixed_block_ >> fi) & 1u) != 0;
      if (side != fixed_side) q = undo(fi, q, side ? 1 : -1);
    }
    return q.z - st.z(static_cast<std::size_t>(h), q.x, q.y);
  };
  // Start at the restored horizon and expand a bracket; older faults move it by at most their throws.
  double z0 = st.z(static_cast<std::size_t>(h), x, y);
  double z = expand_and_solve(s_at, z0, 25.0, 1.0e5, 1e-7);
  (void)dom;
  return std::isnan(z) ? z0 : z;
}

bool FaultNetwork::cannot_move(int fi, const Vec3& q) const {
  const Fault& f = faults_[static_cast<std::size_t>(fi)];
  double a, b;
  f.profile_coords(q, a, b);
  if (std::fabs(a) >= 0.5 * f.length) return true;  // slip has no strike component
  if (f.layered || f.curved) return false;  // a layer can meet the fault far from the point
  if (b + f.slip_b < f.b_lo || b - f.slip_b > f.b_hi) return true;
  if (f.decays && std::fabs(f.sigma(q)) - f.slip_sigma >= f.influence) return true;
  return false;
}

void FaultNetwork::resolve_tips(const desc::Description& d, const Structure& st, BuildLog& log) {
  const Domain dom = d.domain();
  // Oldest first: each fault's tips refer to the stratigraphy as displaced by older faults.
  for (std::size_t r = order_.size(); r-- > 0;) {
    const int k = order_[r];
    Fault& f = faults_[static_cast<std::size_t>(k)];
    if (f.curved) {
      resolve_curved(d, st, log, k, r);
      continue;
    }
    bool older_exist = r + 1 < order_.size();
    auto level = [&](const TipRule& rule, double x, double y, bool upper) -> double {
      if (!older_exist || rule.unfaulted) return rule_level(rule, st, x, y, dom, upper);
      switch (rule.kind) {
        case TipRule::Kind::Horizon: return level_in_frame(st, rule.horizon, x, y, k, dom) - rule.offset;
        case TipRule::Kind::WithinUnit: {
          double zt = level_in_frame(st, rule.top, x, y, k, dom), zb = level_in_frame(st, rule.base, x, y, k, dom);
          return zt + rule.fraction * (zb - zt);
        }
        default: return rule_level(rule, st, x, y, dom, upper);
      }
    };
    // Resolve tip tables by bisection along b for each strike position.
    const std::size_t K = 513;
    double sd = std::sin(f.dip);
    double b_top = (f.p0.z - dom.z_top) / sd, b_base = (f.p0.z - dom.z_base) / sd;
    double margin = 0.5 * (b_base - b_top) + 1000;
    f.a_tab.assign(K, 0.0);
    f.bu_tab.assign(K, 0.0);
    f.bl_tab.assign(K, 0.0);
    f.clamped_tab.assign(K, 0);
    // Strike positions are independent: solve them in parallel (each entry written once).
    parallel_for(K, 16, [&](std::size_t i0, std::size_t i1) {
      for (std::size_t i = i0; i < i1; ++i) {
        double a = -0.5 * f.length + f.length * static_cast<double>(i) / (K - 1);
        auto solve = [&](const TipRule& r, bool upper) -> double {
          if (r.kind == TipRule::Kind::Free) return upper ? b_top : b_base;
          auto g = [&](double b) {
            Vec3 p = f.at(a, b);
            return p.z - level(r, p.x, p.y, upper);
          };
          double lo = b_top - margin, hi = b_base + margin;
          double glo = g(lo), ghi = g(hi);
          if (!(glo > 0 && ghi < 0)) {
            f.clamped_tab[i] = 1;
            return glo <= 0 ? lo : hi;
          }
          return illinois(g, lo, hi, glo, ghi, 1e-7);
        };
        f.a_tab[i] = a;
        f.bu_tab[i] = solve(f.upper, true);
        f.bl_tab[i] = solve(f.lower, false);
      }
    });
    int clamped = 0;
    for (std::size_t i = 0; i < K; ++i) {
      clamped += f.clamped_tab[i];
      if (f.bl_tab[i] <= f.bu_tab[i])
        throw Error(ErrorCode::Build, "fault " + f.name + ": the upper tip lies below the lower tip at strike position " +
                                          fmt_fixed(f.a_tab[i], 1) + " m", "faults[" + std::to_string(k) + "].extent");
    }
    if (clamped) log.note("fault " + f.name + ": tip rule level not met inside the search range at " + std::to_string(clamped) +
                          " strike positions; tips clamped (see checks)");
    // Cached bounds for exact culling of points the fault cannot move.
    f.b_lo = 1e300;
    f.b_hi = -1e300;
    for (std::size_t i = 0; i < K; ++i) {
      double bu, bl;
      f.tips(f.a_tab[i], bu, bl);
      if (bl > bu) { f.b_lo = std::min(f.b_lo, bu); f.b_hi = std::max(f.b_hi, bl); }
    }
    if (f.upper_free()) f.b_lo = -1e300;
    if (f.lower_free()) f.b_hi = 1e300;
    double tabs = std::fabs(f.tmax);
    f.slip_b = tabs * std::fabs(f.slip.z) / std::sin(f.dip) + 1e-9;
    f.slip_sigma = tabs * std::fabs(dot(f.slip, f.n)) + 1e-9;
    f.decays = f.hw->name() == "decay" && f.influence > 0;
    // Throw gradient along the slip: material moved by tau samples the profile tau * slip_b / tmax
    // further along dip. If the throw changes faster than one metre per metre of slip, two
    // restored positions map to one present position (the wall would overlap itself) and the
    // model is not unique, so the build is refused.
    {
      double gmax = 0;
      const double per_tau = std::fabs(f.slip.z) / std::sin(f.dip);
      for (std::size_t i = 0; i < K; i += 4) {
        double bu, bl;
        f.tips(f.a_tab[i], bu, bl);
        if (!(bl > bu) || !std::isfinite(bu) || !std::isfinite(bl)) continue;
        const int NB = 800;
        double db = (bl - bu) / NB;
        double prev = f.throw_at(f.a_tab[i], bu);
        for (int k = 1; k <= NB; ++k) {
          double t = f.throw_at(f.a_tab[i], bu + k * db);
          gmax = std::max(gmax, std::fabs(t - prev) / db * per_tau);
          prev = t;
        }
      }
      f.max_slip_gradient = gmax;
      log.value("fault." + f.name + ".max_throw_gradient_along_slip", fmt_fixed(gmax, 3));
      if (gmax >= 1.0 && !f.layered)
        throw Error(ErrorCode::Build,
                    "fault " + f.name + ": the throw changes by " + fmt_fixed(gmax, 2) +
                        " m per metre of slip near a tip (must be below 1): the hanging wall would overlap itself. "
                        "Lengthen the down-dip taper, move the tips apart or reduce the throw.",
                    "faults." + f.name + ".throw");
      if (gmax >= 0.8 && !f.layered)
        log.note("fault " + f.name + ": throw gradient along slip " + fmt_fixed(gmax, 2) + " is close to the limit of 1");
    }
    double bu0, bl0;
    f.tips(0, bu0, bl0);
    log.value("fault." + f.name + ".tips_depth_at_reference_m",
              fmt_fixed(-f.at(0, bu0).z, 2) + " / " + fmt_fixed(-f.at(0, bl0).z, 2));
    if (layered_) {
      // Where every level meets the plane in this fault's frame (older faults already resolved):
      // exactly at the horizons, and from the level of plane points every 10 m down dip between.
      const std::size_t nh = st.horizons(), KL = 129;
      const std::size_t r_older = r + 1;
      std::vector<std::vector<std::pair<double, double>>> cols(KL);
      f.la_tab.assign(KL, 0.0);
      parallel_for(KL, 4, [&](std::size_t i0, std::size_t i1) {
        for (std::size_t i = i0; i < i1; ++i) {
          double a = -0.5 * f.length + f.length * static_cast<double>(i) / (KL - 1);
          f.la_tab[i] = a;
          auto& c = cols[i];
          for (std::size_t h = 0; h < nh; ++h) {
            auto g = [&](double b) {
              Vec3 p = f.at(a, b);
              return p.z - level_in_frame(st, static_cast<int>(h), p.x, p.y, k, dom);
            };
            double lo = b_top - margin, hi = b_base + margin;
            double glo = g(lo), ghi = g(hi);
            if (glo > 0 && ghi < 0) c.emplace_back(illinois(g, lo, hi, glo, ghi, 1e-7), static_cast<double>(h));
          }
          const double db = 10.0;
          for (double b = b_top; b <= b_base + 1e-9; b += db)
            c.emplace_back(b, solve_layer(f.at(a, b), r_older, nullptr, 0));
          std::sort(c.begin(), c.end());
          // keep the level increasing down dip
          std::vector<std::pair<double, double>> m;
          for (auto& e : c)
            if (m.empty() || e.second > m.back().second + 1e-12) m.push_back(e);
          c.swap(m);
        }
      });
      f.l_start.assign(KL + 1, 0);
      f.ls_tab.clear();
      f.lb_tab.clear();
      for (std::size_t i = 0; i < KL; ++i) {
        f.l_start[i] = f.ls_tab.size();
        for (auto& [b, sl] : cols[i]) { f.lb_tab.push_back(b); f.ls_tab.push_back(sl); }
        if (cols[i].empty()) { f.lb_tab.push_back(b_top); f.ls_tab.push_back(-1.0); }
      }
      f.l_start[KL] = f.ls_tab.size();
      // Layers must keep their order after faulting: between neighbouring levels the throw may
      // change by less than the layer is thick at the fault.
      double worst = 0;
      for (std::size_t i = 0; i < KL; i += 8)
        for (std::size_t j = f.l_start[i] + 1; j < f.l_start[i + 1]; ++j) {
          double ba = f.lb_tab[j - 1], bb = f.lb_tab[j];
          double dz = (bb - ba) * std::sin(f.dip);
          double dt = std::fabs(f.throw_at(f.la_tab[i], bb) - f.throw_at(f.la_tab[i], ba));
          if (dz > 1e-6) worst = std::max(worst, dt / dz);
        }
      log.value("fault." + f.name + ".max_layer_throw_change_per_thickness", fmt_fixed(worst, 3));
      for (double e : {-0.5, 0.5}) {
        Vec3 end = f.p0 + f.s * (e * f.length);
        if (dom.inside_xy(end.x, end.y)) {
          const bool steps = f.pstrike == ProfileType::FlatTopped && f.taper_strike <= 0;
          if (steps)
            log.note("fault " + f.name + ": its strike end lies inside the domain; with layer-rigid hanging walls and no "
                     "taper along strike the block offset ends there along a line that is not a fault (give the strike "
                     "profile a taper, or make the fault longer than the model)");
          else
            log.note("fault " + f.name + ": its strike end lies inside the domain; the throw tapers to zero there along "
                     "strike, so the layer-rigid block offset dies out towards the end");
          break;
        }
      }
      if (worst >= 1.0)
        throw Error(ErrorCode::Build, "fault " + f.name + ": the throw changes faster than the layers thicken at the fault (" +
                                          fmt_fixed(worst, 2) + " m per metre); layers would overlap. Reduce the throw or move the tips apart.",
                    "faults." + f.name + ".throw");
    }
  }
}

void FaultNetwork::resolve_curved(const desc::Description& d, const Structure& st, BuildLog& log, int k, std::size_t r) {
  Fault& f = faults_[static_cast<std::size_t>(k)];
  const CurvedFault& c = *f.curved;
  const Domain dom = d.domain();
  // Tips are free (rule F10): the fault spans the model from its top to beyond the far side, in
  // the horizontal coordinate u along the dip direction.
  double u_top = c.u_at(-dom.z_top), u_far = u_top;
  for (double cx : {dom.x0, dom.x0 + dom.sx})
    for (double cy : {dom.y0, dom.y0 + dom.sy}) u_far = std::max(u_far, dot(Vec3{cx, cy, f.p0.z} - f.p0, f.dh));
  u_far += 1000 + std::fabs(c.heave);
  const std::size_t K = 513;
  f.a_tab.assign(K, 0.0);
  f.bu_tab.assign(K, u_top);
  f.bl_tab.assign(K, u_far);
  f.clamped_tab.assign(K, 0);
  for (std::size_t i = 0; i < K; ++i) f.a_tab[i] = -0.5 * f.length + f.length * static_cast<double>(i) / (K - 1);
  f.b_lo = -1e300;
  f.b_hi = 1e300;
  f.slip_b = f.slip_sigma = 0;
  f.max_slip_gradient = 0;
  log.value("fault." + f.name + ".geometry", d.faults[static_cast<std::size_t>(k)].geom_type + ", " + f.hw_name +
                                                (c.ta > 0 ? (c.s < 0 ? " antithetic " : " synthetic ") + fmt_fixed(rad2deg(std::atan(c.ta)), 2) + " deg" : ""));

  // Report at the reference strike position (a = 0): per horizon heave, throw, rollover width and
  // maximum rollover dip (STX-TN-05 section 2.4).
  const double h = f.heave_at(0);
  if (h <= 0) return;
  const double tol = d.tol_rollover.present ? d.tol_rollover.si() : 1.0;
  const bool older = r + 1 < order_.size();
  auto at_u = [&](double u) { return f.p0 + f.dh * u; };
  auto hdepth = [&](std::size_t hz, double u) {  // horizon depth in the fault's frame at u (a = 0)
    Vec3 p = at_u(u);
    return older ? -level_in_frame(st, static_cast<int>(hz), p.x, p.y, k, dom) : -st.z(hz, p.x, p.y);
  };
  for (std::size_t hz = 0; hz < st.horizons(); ++hz) {
    const std::string key = f.name + "." + d.horizons[hz].name;
    // Footwall cut-off: where the fault reaches the horizon.
    auto g = [&](double u) { return c.zf(u) - hdepth(hz, u); };
    double lo = u_top, glo = g(lo), uc = kNaN;
    for (double u = u_top + 50; u <= u_far; u += 50) {
      double gu = g(u);
      if ((gu >= 0) != (glo >= 0)) { uc = illinois(g, lo, u, glo, gu, 1e-9); break; }
      lo = u;
      glo = gu;
    }
    if (std::isnan(uc)) {
      log.note("fault " + f.name + ": horizon " + d.horizons[hz].name + " is not cut by the curved fault at the reference strike position");
      continue;
    }
    double zc = c.zf(uc), u1c, z1c;
    c.forward(uc, zc, h, u1c, z1c);
    log.value(key + ".heave", fmt_fixed(h, 3));
    log.value(key + ".throw", fmt_fixed(z1c - zc, 3));
    // Rollover: the horizon's material subsides by delta(xi) = G(xi + h) - G(xi); it ends where that
    // falls to the tolerance.
    const double xic = uc - c.s * zc * c.ta;
    auto delta = [&](double xi) { return c.G(xi + h) - c.G(xi) - tol; };
    double xr = kNaN, a0 = xic, d0 = delta(a0);
    for (double xi = xic + 5; xi <= xic + 3 * (u_far - u_top); xi += 5) {
      double dx = delta(xi);
      if ((dx > 0) != (d0 > 0)) { xr = illinois(delta, a0, xi, d0, dx, 1e-9); break; }
      a0 = xi;
      d0 = dx;
    }
    // Present position of the horizon's material on shear line xi.
    auto present = [&](double xi, double& u1, double& z1) {
      double z = expand_and_solve([&](double zz) { return hdepth(hz, xi + c.s * zz * c.ta) - zz; }, zc, 25.0, 1.0e6, 1e-9);
      c.forward(xi + c.s * z * c.ta, z, h, u1, z1);
    };
    if (!std::isnan(xr)) {
      double u1r, z1r;
      present(xr, u1r, z1r);
      log.value(key + ".rollover_width", fmt_fixed(u1r - u1c, 2));
      double dipmax = 0, pu = u1c, pz = z1c;
      for (int i = 1; i <= 200; ++i) {
        double xi = xic + (xr - xic) * i / 200.0, u1, z1;
        present(xi, u1, z1);
        if (u1 - pu > 1e-9) dipmax = std::max(dipmax, std::fabs(std::atan2(z1 - pz, u1 - pu)));
        pu = u1;
        pz = z1;
      }
      log.value(key + ".rollover_dip_max", fmt_fixed(rad2deg(dipmax), 3));
    }
  }
}

double FaultNetwork::next_fault_distance(int f, const Vec3& q, double sigma) const {
  const Fault& F = faults_[static_cast<std::size_t>(f)];
  Vec3 base = q - F.n * sigma;
  Vec3 dir = F.n * (sigma >= 0 ? 1.0 : -1.0);
  double best = kInf;
  for (std::size_t g = 0; g < faults_.size(); ++g) {
    if (static_cast<int>(g) == f) continue;
    const Fault& G = faults_[g];
    double den = dot(G.n, dir);
    if (std::fabs(den) < 1e-9) continue;
    double t = -G.sigma(base) / den;
    if (t <= std::fabs(sigma) || t >= best) continue;
    Vec3 x = base + dir * t;
    double a, b;
    G.ab(x, a, b);
    if (G.extent_coord(a, b) <= 1) best = t;
  }
  return best;
}

Vec3 FaultNetwork::relative_displacement(int fi, const Vec3& q, bool* active) const {
  const Fault& f = faults_[static_cast<std::size_t>(fi)];
  if (active) *active = true;
  if (f.abuts >= 0) {
    bool hw = faults_[static_cast<std::size_t>(f.abuts)].sigma(q) > 0;
    if (hw != f.abut_hw_side) {
      if (active) *active = false;
      return {0, 0, 0};
    }
  }
  if (f.curved) {
    if (f.sigma(q) <= 0) return {0, 0, 0};
    double a = dot(q - f.p0, f.s), u = dot(q - f.p0, f.dh), h = f.heave_at(a);
    if (h == 0) return {0, 0, 0};
    double u1, z1;
    f.curved->forward(u, -q.z, h, u1, z1);
    return f.dh * (u1 - u) + Vec3{0, 0, -z1 - q.z};
  }
  if (layered_) {
    double s = solve_layer(q, static_cast<std::size_t>(rank_[static_cast<std::size_t>(fi)]) + 1, nullptr, 0);
    return f.slip * layer_throw(fi, s, q);
  }
  double a, b;
  f.profile_coords(q, a, b);
  double t = f.throw_at(a, b);
  if (t == 0) return {0, 0, 0};
  double sig = f.sigma(q);
  double next = f.hw->needs_next_fault() ? next_fault_distance(fi, q, sig) : kInf;
  double sc = f.hw->scale(sig, next, f.influence);
  return f.slip * (t * sc);
}

Vec3 FaultNetwork::undo(int fi, const Vec3& q, int sgn, double* tau_out) const {
  const Fault& f = faults_[static_cast<std::size_t>(fi)];
  if (f.curved) {
    // Only the hanging wall moves (rule F11). tau: vertical displacement (present minus restored depth).
    if (tau_out) *tau_out = 0;
    if (sgn != 1) return q;
    if (f.abuts >= 0 && (faults_[static_cast<std::size_t>(f.abuts)].sigma(q) > 0) != f.abut_hw_side) return q;
    double a = dot(q - f.p0, f.s), u1 = dot(q - f.p0, f.dh), h = f.heave_at(a);
    if (h == 0) return q;
    double u, z;
    f.curved->inverse(u1, -q.z, h, u, z);
    if (tau_out) *tau_out = -q.z - z;
    Vec3 r = q + f.dh * (u - u1);
    r.z = -z;
    return r;
  }
  if (layered_) {
    const std::size_t r = static_cast<std::size_t>(rank_[static_cast<std::size_t>(fi)]);
    double s = solve_layer(q, r, nullptr, sgn);
    bool active = true;
    if (f.abuts >= 0) active = (faults_[static_cast<std::size_t>(f.abuts)].sigma(q) > 0) == f.abut_hw_side;
    double tau = active ? layer_throw(fi, s, q) : 0.0;
    if (tau_out) *tau_out = tau;
    return q - f.slip * (sgn * tau);
  }
  const Vec3 slip = f.slip;
  // tau(q0) along the line q0 = q - sgn * slip * tau; g(tau) = tau - amount(q0(tau)).
  auto amount = [&](double tau) {
    Vec3 q0 = q - slip * (sgn * tau);
    Vec3 r = relative_displacement(fi, q0);
    return dot(r, slip) / dot(slip, slip);
  };
  double tmax = f.tmax;
  double lo = std::min(0.0, tmax), hi = std::max(0.0, tmax);
  double glo = lo - amount(lo), ghi = hi - amount(hi);
  double tau;
  // |amount| <= |tmax|, so g(lo) <= 0 <= g(hi) in exact arithmetic. When the root sits at an
  // end of the bracket (the profile's plateau at tmax), rounding can give g a few ulps of the
  // wrong sign there; such an end is the root.
  const double gtol = 1e-9 * std::max(1.0, std::fabs(tmax));
  if (std::fabs(glo) <= gtol) tau = lo;
  else if (std::fabs(ghi) <= gtol) tau = hi;
  else if ((glo < 0) == (ghi < 0)) tau = std::fabs(glo) < std::fabs(ghi) ? lo : hi;  // not reached
  else {
    tau = illinois([&](double t) { return t - amount(t); }, lo, hi, glo, ghi, 1e-11, 1e-11, 200);
  }
  if (tau_out) *tau_out = tau;
  return q - slip * (sgn * tau);
}

Vec3 FaultNetwork::restore(const Vec3& p, BlockCode* block, const BlockCode* forced) const {
  if (layered_) {
    double s = solve_layer(p, 0, forced, 0);
    return chain(p, 0, order_.size(), s, block, forced, 0);
  }
  const bool fast = optimised_kernels();
  Vec3 q = p;
  BlockCode code = 0;
  for (int fi : order_) {
    const Fault& f = faults_[static_cast<std::size_t>(fi)];
    bool side = forced ? ((*forced >> fi) & 1u) != 0 : f.sigma(q) > 0;
    if (side) code |= (BlockCode(1) << fi);
    bool fixed_side = ((fixed_block_ >> fi) & 1u) != 0;
    if (side == fixed_side) continue;
    if (fast && cannot_move(fi, q)) continue;
    q = undo(fi, q, side ? 1 : -1);
  }
  if (block) *block = code;
  return q;
}

Vec3 FaultNetwork::restore_from_frame(const Vec3& p, int f, bool hw) const {
  if (layered_) {
    const std::size_t r = static_cast<std::size_t>(rank_[static_cast<std::size_t>(f)]);
    double s = solve_layer(p, r, nullptr, hw ? 2 : -2);
    return chain(p, r, order_.size(), s, nullptr, nullptr, hw ? 2 : -2);
  }
  Vec3 q = p;
  bool fixed_side = ((fixed_block_ >> f) & 1u) != 0;
  if (hw != fixed_side) q = undo(f, q, hw ? 1 : -1);
  for (std::size_t k = static_cast<std::size_t>(rank_[static_cast<std::size_t>(f)]) + 1; k < order_.size(); ++k) {
    int fi = order_[k];
    const Fault& F = faults_[static_cast<std::size_t>(fi)];
    bool side = F.sigma(q) > 0;
    bool fs = ((fixed_block_ >> fi) & 1u) != 0;
    if (side != fs) q = undo(fi, q, side ? 1 : -1);
  }
  return q;
}

Vec3 FaultNetwork::to_fault_frame(const Vec3& p, int f) const {
  if (layered_) {
    double s = solve_layer(p, 0, nullptr, 0);
    return chain(p, 0, static_cast<std::size_t>(rank_[static_cast<std::size_t>(f)]), s, nullptr, nullptr, 0);
  }
  const bool fast = optimised_kernels();
  Vec3 q = p;
  for (int fi : order_) {
    if (fi == f) break;
    const Fault& F = faults_[static_cast<std::size_t>(fi)];
    bool side = F.sigma(q) > 0;
    bool fixed_side = ((fixed_block_ >> fi) & 1u) != 0;
    if (side == fixed_side) continue;
    if (fast && cannot_move(fi, q)) continue;
    q = undo(fi, q, side ? 1 : -1);
  }
  return q;
}

// ---- layer-rigid hanging walls -------------------------------------------------------------
namespace {
/// Level of elevation zr in a column whose horizon elevations are z[0..nh).
double level_in_column(double zr, const double* z, std::size_t nh, double z_top, double z_base) {
  const double top = std::max(z_top, z[0]), base = std::min(z_base, z[nh - 1]);
  if (zr >= top) return -1.0 - (zr - top) / 1000.0;
  if (zr <= base) return static_cast<double>(nh) + (base - zr) / 1000.0;
  double up = top;
  for (std::size_t l = 0; l <= nh; ++l) {
    double dn = l < nh ? z[l] : base;
    if (zr >= dn) {
      double th = up - dn;
      return static_cast<double>(l) - 1.0 + (th > 1e-12 ? (up - zr) / th : 0.0);
    }
    up = std::min(up, dn);
  }
  return static_cast<double>(nh);
}
}  // namespace

double FaultNetwork::strat_coord(const Vec3& r) const {
  const std::size_t nh = st_->horizons();
  double zs[256];
  std::vector<double> big;
  double* z = zs;
  if (nh > 256) { big.resize(nh); z = big.data(); }
  st_->z_all(r.x, r.y, z);
  return level_in_column(r.z, z, nh, z_top_, z_base_);
}

double FaultNetwork::layer_displacement(int fi, const Vec3& p, bool* active) const {
  const Fault& f = faults_[static_cast<std::size_t>(fi)];
  Vec3 q = to_fault_frame(p, fi);
  if (active) *active = true;
  if (f.abuts >= 0 && (faults_[static_cast<std::size_t>(f.abuts)].sigma(q) > 0) != f.abut_hw_side) {
    if (active) *active = false;
    return 0.0;
  }
  double s = solve_layer(p, 0, nullptr, 0);
  return std::fabs(layer_throw(fi, s, q)) * norm(f.slip);
}

double FaultNetwork::layer_throw(int fi, double s, const Vec3& q) const {
  const Fault& f = faults_[static_cast<std::size_t>(fi)];
  double a = dot(q - f.p0, f.s);
  if (std::fabs(a) >= 0.5 * f.length) return 0.0;
  return f.throw_at(a, f.layer_b(s, a));
}

Vec3 FaultNetwork::chain(const Vec3& q0, std::size_t r0, std::size_t r1, double s, BlockCode* code, const BlockCode* forced,
                         int first) const {
  Vec3 q = q0;
  BlockCode c = 0;
  for (std::size_t r = r0; r < r1 && r < order_.size(); ++r) {
    const int fi = order_[r];
    const Fault& F = faults_[static_cast<std::size_t>(fi)];
    const bool fixed_side = ((fixed_block_ >> fi) & 1u) != 0;
    bool side;
    int sgn = 0;
    if (r == r0 && (first == 1 || first == -1)) {
      side = first > 0;
      sgn = first;
    } else {
      if (r == r0 && (first == 2 || first == -2)) side = first > 0;
      else side = forced ? ((*forced >> fi) & 1u) != 0 : F.sigma(q) > 0;
      if (side != fixed_side) sgn = side ? 1 : -1;
    }
    if (side) c |= (BlockCode(1) << fi);
    if (sgn == 0) continue;
    if (F.abuts >= 0 && (faults_[static_cast<std::size_t>(F.abuts)].sigma(q) > 0) != F.abut_hw_side) continue;
    q = q - F.slip * (sgn * layer_throw(fi, s, q));
  }
  if (code) *code = c;
  return q;
}

double FaultNetwork::solve_layer(const Vec3& q, std::size_t r0, const BlockCode* forced, int first) const {
  // With vertical slips the column never moves: read its horizons once.
  const std::size_t nh = st_->horizons();
  double zs[256];
  std::vector<double> big;
  double* zc = zs;
  if (nh > 256) { big.resize(nh); zc = big.data(); }
  if (vertical_) st_->z_all(q.x, q.y, zc);
  auto g = [&](double s) {
    Vec3 r = chain(q, r0, order_.size(), s, nullptr, forced, first);
    return (vertical_ ? level_in_column(r.z, zc, nh, z_top_, z_base_) : strat_coord(r)) - s;
  };
  // Bracket from the largest possible displacement; widen if a dipping slip moved out of it.
  double d = 0;
  for (std::size_t r = r0; r < order_.size(); ++r) {
    const Fault& F = faults_[static_cast<std::size_t>(order_[r])];
    d += std::fabs(F.tmax) * norm(F.slip);
  }
  double lo = strat_coord(q + Vec3{0, 0, d + 1e-6}), hi = strat_coord(q - Vec3{0, 0, d + 1e-6});
  double glo = g(lo), ghi = g(hi);
  if (!(glo >= 0 && ghi <= 0)) {
    lo = -1.0 - 1e4;
    hi = static_cast<double>(st_->horizons()) + 1e4;
    glo = g(lo);
    ghi = g(hi);
  }
  if (glo == 0) return lo;
  if (ghi == 0) return hi;
  if ((glo > 0) == (ghi > 0)) return std::fabs(glo) < std::fabs(ghi) ? lo : hi;
  return illinois(g, lo, hi, glo, ghi, 1e-12, 1e-13, 200);
}

}  // namespace stratonix
