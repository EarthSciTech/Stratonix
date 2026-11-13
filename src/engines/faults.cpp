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
class DecayHW final : public HangingWallField {
 public:
  double scale(double sigma, double, double infl) const override {
    return infl > 0 ? std::max(0.0, 1 - std::fabs(sigma) / infl) : 1.0;
  }
  std::string name() const override { return "decay"; }
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
STRATONIX_REGISTER(hw_decay, registry::Kind::HangingWall, "decay", "1.0.0",
                   "Displacement decays linearly to zero at an influence distance",
                   [] { return std::make_unique<DecayHW>(); })
STRATONIX_REGISTER(hw_linear, registry::Kind::HangingWall, "linear_between_faults", "1.0.0",
                   "Displacement decays linearly to zero at the next fault plane",
                   [] { return std::make_unique<LinearBetweenHW>(); })

void register_fault_plugins() {
  STRATONIX_PLUGIN_ANCHOR(kin_dip_slip);
  STRATONIX_PLUGIN_ANCHOR(kin_vertical);
  STRATONIX_PLUGIN_ANCHOR(hw_rigid);
  STRATONIX_PLUGIN_ANCHOR(hw_decay);
  STRATONIX_PLUGIN_ANCHOR(hw_linear);
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
    f.dip = fd.dip.si();
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
      f.p0 = {fd.reference_point.si(0), fd.reference_point.si(1), -fd.reference_depth.si()};
    }
    f.s = azimuth_vector(f.strike);
    f.dh = azimuth_vector(f.strike + kPi / 2);
    f.dd = f.dh * std::cos(f.dip) + Vec3{0, 0, -std::sin(f.dip)};
    f.n = f.dh * std::sin(f.dip) + Vec3{0, 0, std::cos(f.dip)};
    f.length = fd.length.si();
    f.elliptical = fd.tip_line_shape == "elliptical";
    f.upper = tip_rule(fd.upper, d);
    f.lower = tip_rule(fd.lower, d);
    f.tmax = fd.throw_max.si();
    f.pdip = profile_from(fd.profile_dip.type);
    f.pstrike = profile_from(fd.profile_strike.type);
    f.taper_dip = fd.profile_dip.taper.present ? fd.profile_dip.taper.si() : 0.25;
    f.taper_strike = fd.profile_strike.taper.present ? fd.profile_strike.taper.si() : 0.25;
    f.kin = reg.create<Kinematics>(registry::Kind::Kinematics, fd.slip_mode);
    f.slip = f.kin->slip_vector(f.dd, f.dip);
    f.hw = reg.create<HangingWallField>(registry::Kind::HangingWall, fd.hw_type);
    f.influence = fd.influence_distance.present ? fd.influence_distance.si() : 0.0;
    if (!fd.abuts_fault.empty()) {
      f.abuts = d.fault_index(fd.abuts_fault);
      f.abut_hw_side = fd.abuts_side == "hanging_wall";
    }
    faults_.push_back(std::move(f));
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
  resolve_tips(d, st, log);
}

/// Elevation of horizon h at (x, y) in the frame of fault `rank` (stratigraphy displaced by the
/// faults older than it, which are already resolved).
double FaultNetwork::level_in_frame(const Structure& st, int h, double x, double y, int f, const Domain& dom) const {
  auto s_at = [&](double z) {
    Vec3 q{x, y, z};
    for (std::size_t k = static_cast<std::size_t>(rank_[static_cast<std::size_t>(f)]) + 1; k < order_.size(); ++k) {
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
    bool older_exist = r + 1 < order_.size();
    auto level = [&](const TipRule& rule, double x, double y, bool upper) -> double {
      if (!older_exist) return rule_level(rule, st, x, y, dom, upper);
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
    double bu0, bl0;
    f.tips(0, bu0, bl0);
    log.value("fault." + f.name + ".tips_depth_at_reference_m",
              fmt_fixed(-f.at(0, bu0).z, 2) + " / " + fmt_fixed(-f.at(0, bl0).z, 2));
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
  if (glo == 0) tau = lo;
  else if (ghi == 0) tau = hi;
  else if ((glo < 0) == (ghi < 0)) tau = amount(0.0);  // no bracket (cannot happen for |amount| <= |tmax|)
  else {
    int side = 0;
    tau = 0;
    for (int it = 0; it < 100; ++it) {
      tau = (lo * ghi - hi * glo) / (ghi - glo);
      double gt = tau - amount(tau);
      if (std::fabs(gt) < 1e-11 || hi - lo < 1e-11) break;
      if ((gt < 0) == (glo < 0)) {
        lo = tau; glo = gt;
        if (side == -1) ghi *= 0.5;
        side = -1;
      } else {
        hi = tau; ghi = gt;
        if (side == 1) glo *= 0.5;
        side = 1;
      }
    }
  }
  if (tau_out) *tau_out = tau;
  return q - slip * (sgn * tau);
}

Vec3 FaultNetwork::restore(const Vec3& p, BlockCode* block, const BlockCode* forced) const {
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

}  // namespace stratonix
