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
STRATONIX_DECLARE_PLUGIN_ANCHOR(kin_oblique);
STRATONIX_DECLARE_PLUGIN_ANCHOR(kin_strike_slip);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_rigid);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_layer_rigid);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_decay);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_linear);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_vertical_shear);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_inclined_shear);
STRATONIX_DECLARE_PLUGIN_ANCHOR(hw_fpf);

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
void CurvedFault::setup_fpf() {
  const std::size_t n = pu.size(), ns = n - 1;
  fe_u.assign(ns, 0); fe_z.assign(ns, 0); fn_u.assign(ns, 0); fn_z.assign(ns, 0);
  for (std::size_t k = 0; k < ns; ++k) {
    double du = pu[k + 1] - pu[k], dz = pz[k + 1] - pz[k], l = std::hypot(du, dz);
    fe_u[k] = du / l; fe_z[k] = dz / l;
    fn_u[k] = fe_z[k]; fn_z[k] = -fe_u[k];  // towards the hanging wall (shallower)
  }
  fm_u.assign(n, 0); fm_z.assign(n, 0);
  for (std::size_t k = 0; k < n; ++k) {
    if (k == 0) { fm_u[k] = fn_u[0]; fm_z[k] = fn_z[0]; continue; }
    if (k == ns) { fm_u[k] = fn_u[ns - 1]; fm_z[k] = fn_z[ns - 1]; continue; }
    double c = 1 + fn_u[k - 1] * fn_u[k] + fn_z[k - 1] * fn_z[k];
    fm_u[k] = (fn_u[k - 1] + fn_u[k]) / c;
    fm_z[k] = (fn_z[k - 1] + fn_z[k]) / c;
  }
  fl0.assign(ns, 0); fl1.assign(ns, 0);
  fpf_dmax = std::numeric_limits<double>::infinity();
  for (std::size_t k = 0; k < ns; ++k) {
    fl0[k] = (pu[k + 1] - pu[k]) * fe_u[k] + (pz[k + 1] - pz[k]) * fe_z[k];
    fl1[k] = (fm_u[k + 1] - fm_u[k]) * fe_u[k] + (fm_z[k + 1] - fm_z[k]) * fe_z[k];
    // the end segments are unbounded; inner ones vanish where the bisectors cross
    if (k > 0 && k + 1 < ns && fl1[k] < 0) fpf_dmax = std::min(fpf_dmax, -fl0[k] / fl1[k]);
  }
  fpf = true;
}
void CurvedFault::fpf_coords(double u, double z, double& s, double& d) const {
  const std::size_t ns = pu.size() - 1;
  std::size_t k0;
  if (u <= pu[0]) k0 = 0;
  else if (u >= pu[ns]) k0 = ns - 1;
  else k0 = static_cast<std::size_t>(std::upper_bound(pu.begin(), pu.end(), u) - pu.begin()) - 1;
  auto test = [&](std::size_t k, double& t, double& dd) {
    dd = (u - pu[k]) * fn_u[k] + (z - pz[k]) * fn_z[k];
    t = (u - pu[k]) * fe_u[k] + (z - pz[k]) * fe_z[k] - dd * (fm_u[k] * fe_u[k] + fm_z[k] * fe_z[k]);
    return (k == 0 || t >= -1e-9) && (k + 1 == ns || t <= fl0[k] + dd * fl1[k] + 1e-9);
  };
  double t = 0, dd = 0;
  std::size_t k = k0;
  bool found = test(k0, t, dd);
  for (std::size_t r = 1; !found && r <= ns; ++r) {
    if (k0 + r < ns && test(k0 + r, t, dd)) { k = k0 + r; found = true; break; }
    if (r <= k0 && test(k0 - r, t, dd)) { k = k0 - r; found = true; break; }
  }
  if (!found) test(k = k0, t, dd);
  s = t;
  for (std::size_t j = 0; j < k; ++j) s += fl0[j] + dd * fl1[j];
  d = dd;
}
void CurvedFault::fpf_at(double s, double d, double& u, double& z) const {
  const std::size_t ns = pu.size() - 1;
  double acc = 0;
  for (std::size_t k = 0; k < ns; ++k) {
    double len = fl0[k] + d * fl1[k];
    if (s < acc + len || k + 1 == ns) {
      if (k == 0 && s < 0) acc = 0;
      u = pu[k] + fm_u[k] * d + fe_u[k] * (s - acc);
      z = pz[k] + fm_z[k] * d + fe_z[k] * (s - acc);
      return;
    }
    acc += len;
  }
}
void CurvedFault::forward(double u, double z, double h, double& u1, double& z1) const {
  if (fpf) {
    double s, d;
    fpf_coords(u, z, s, d);
    fpf_at(s + h, d, u1, z1);
    return;
  }
  double xi = u - s * z * ta, xi1 = xi + h;
  z1 = z + G(xi1) - G(xi);
  u1 = xi1 + s * z1 * ta;
}
void CurvedFault::inverse(double u1, double z1, double h, double& u, double& z) const {
  if (fpf) {
    double s, d;
    fpf_coords(u1, z1, s, d);
    fpf_at(s - h, d, u, z);
    return;
  }
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
/// Oblique slip: unit vector at the rake in the fault plane (Aki and Richards: from the strike
/// direction, positive up dip); the displacement magnitude is the net slip.
class Oblique final : public Kinematics {
 public:
  Vec3 slip_vector(const Vec3& dd, double dip) const override { return dd / std::sin(dip); }
  Vec3 slip_vector(const Vec3& s, const Vec3& dd, double, double rake) const override {
    return s * std::cos(rake) - dd * std::sin(rake);
  }
  std::string name() const override { return "oblique"; }
};
class StrikeSlip final : public Kinematics {
 public:
  Vec3 slip_vector(const Vec3& dd, double dip) const override { return dd / std::sin(dip); }
  Vec3 slip_vector(const Vec3& s, const Vec3&, double, double rake) const override { return s * std::cos(rake); }
  std::string name() const override { return "strike_slip"; }
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
STRATONIX_REGISTER(kin_oblique, registry::Kind::Kinematics, "oblique", "1.0.0",
                   "Hanging wall moves at a rake in the fault plane by the net slip",
                   [] { return std::make_unique<Oblique>(); })
STRATONIX_REGISTER(kin_strike_slip, registry::Kind::Kinematics, "strike_slip", "1.0.0",
                   "Hanging wall moves along strike by the net slip (left- or right-lateral)",
                   [] { return std::make_unique<StrikeSlip>(); })
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
STRATONIX_REGISTER(hw_fpf, registry::Kind::HangingWall, "fault_parallel_flow", "1.0.0",
                   "Curved faults: hanging wall flows parallel to the fault (fault-bend folds)",
                   [] { return std::make_unique<ShearHW>("fault_parallel_flow"); })
STRATONIX_REGISTER(hw_inclined_shear, registry::Kind::HangingWall, "inclined_shear", "1.0.0",
                   "Curved faults: hanging wall collapses onto the fault by simple shear on inclined lines",
                   [] { return std::make_unique<ShearHW>("inclined_shear"); })

void register_fault_plugins() {
  STRATONIX_PLUGIN_ANCHOR(kin_dip_slip);
  STRATONIX_PLUGIN_ANCHOR(kin_vertical);
  STRATONIX_PLUGIN_ANCHOR(kin_oblique);
  STRATONIX_PLUGIN_ANCHOR(kin_strike_slip);
  STRATONIX_PLUGIN_ANCHOR(hw_rigid);
  STRATONIX_PLUGIN_ANCHOR(hw_layer_rigid);
  STRATONIX_PLUGIN_ANCHOR(hw_decay);
  STRATONIX_PLUGIN_ANCHOR(hw_linear);
  STRATONIX_PLUGIN_ANCHOR(hw_vertical_shear);
  STRATONIX_PLUGIN_ANCHOR(hw_inclined_shear);
  STRATONIX_PLUGIN_ANCHOR(hw_fpf);
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
  double fs;
  if (relay_with >= 0 && a >= std::min(ra_one, ra_zero) && a <= std::max(ra_one, ra_zero))
    fs = (a - ra_zero) / (ra_one - ra_zero);
  else
    fs = profile_value(pstrike, taper_strike, xi);
  return tmax * profile_value(pdip, taper_dip, eta) * fs;
}

double Fault::slip_at(double a) const {
  double L2 = 0.5 * length;
  if (std::fabs(a) >= L2) return 0;
  return std::fabs(tmax) / std::sin(dip) * profile_value(pstrike, taper_strike, (a + L2) / length);
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
      if (fd.hw_type == "fault_parallel_flow") {
        if (c->shape == CurvedFault::Shape::Listric) {
          // 600-segment polyline from the model top to near the detachment (STX-TN-06 section 2)
          const double ztop = d.domain().z_top;
          double ua = c->u_at(std::min(-ztop, c->z0)) - 1000, ub = -(c->zd - c->z0) / c->tan0 * std::log(1e-3);
          for (int i = 0; i <= 600; ++i) {
            double u = ua + (ub - ua) * i / 600.0;
            c->pu.push_back(u);
            c->pz.push_back(c->zf(u));
          }
          c->shape = CurvedFault::Shape::RampFlat;
        }
        c->heave = (fd.sense == "reverse" ? -1.0 : 1.0) * fd.net_slip.si();
        c->setup_fpf();
      } else {
        c->heave = fd.heave.si();
      }
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
    if (f.curved) f.tmax = f.curved->heave;
    else if (fd.strike_component()) f.tmax = fd.net_slip.si();
    else f.tmax = fd.sense == "reverse" ? -fd.throw_max.si() : fd.throw_max.si();
    if (fd.slip_mode == "oblique") f.rake = fd.rake.si();
    if (fd.slip_mode == "strike_slip") f.rake = fd.sense == "right_lateral" ? kPi : 0.0;
    f.pdip = profile_from(fd.profile_dip.type);
    f.pstrike = profile_from(fd.profile_strike.type);
    f.taper_dip = fd.profile_dip.taper.present ? fd.profile_dip.taper.si() : 0.25;
    f.taper_strike = fd.profile_strike.taper.present ? fd.profile_strike.taper.si() : 0.25;
    f.kin = reg.create<Kinematics>(registry::Kind::Kinematics, fd.slip_mode);
    f.slip = f.kin->slip_vector(f.s, f.dd, f.dip, f.rake);
    f.hw = reg.create<HangingWallField>(registry::Kind::HangingWall, fd.hw_type);
    f.hw_name = f.hw->name();
    f.layered = f.hw->name() == "layer_rigid";
    f.influence = fd.influence_distance.present ? fd.influence_distance.si() : 0.0;
    if (!fd.abuts_fault.empty()) {
      f.abuts = d.fault_index(fd.abuts_fault);
      f.abut_hw_side = fd.abuts_side == "hanging_wall";
      f.transfer = fd.abuts_transfer;
    }
    if (fd.fold_type == "trishear") {
      auto t = std::make_shared<TrishearFold>();
      t->tanphi = std::tan(0.5 * fd.fold_apical_angle.si());
      t->ps = fd.fold_ps_ratio.si();
      t->v0 = f.tmax < 0 ? 1.0 : -1.0;  // reverse: towards the tip
      f.tri = t;
    }
    if (!fd.drag_type.empty()) {
      f.drag = fd.drag_type == "reverse" ? -1 : 1;
      f.drag_fraction = fd.drag_fraction.si();
      f.drag_width = fd.drag_width.si();
    }
    faults_.push_back(std::move(f));
  }
  // Relays: complementary linear ramps of the strike factor over the overlap (STX-TN-05 s.4.2).
  for (std::size_t k = 0; k < d.faults.size(); ++k) {
    if (d.faults[k].relay_fault.empty()) continue;
    Fault& F = faults_[k];
    Fault& G = faults_[static_cast<std::size_t>(d.fault_index(d.faults[k].relay_fault))];
    double c2 = dot(G.p0 - F.p0, F.s);
    double lo = std::max(-0.5 * F.length, c2 - 0.5 * G.length), hi = std::min(0.5 * F.length, c2 + 0.5 * G.length);
    // F's inner tip is the end of F inside the overlap; F is full at the other end of the overlap.
    bool f_ends_high = hi >= 0.5 * F.length - 1e-6;
    F.relay_with = G.index;
    G.relay_with = F.index;
    F.ra_one = f_ends_high ? lo : hi;
    F.ra_zero = f_ends_high ? hi : lo;
    auto g_coord = [&](double aF) { return dot(F.p0 + F.s * aF - G.p0, G.s); };
    G.ra_one = g_coord(F.ra_zero);
    G.ra_zero = g_coord(F.ra_one);
    log.value("fault." + F.name + ".relay." + G.name + ".overlap_m", fmt_fixed(hi - lo, 2));
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
  // Linked splays (abutting with transfer): the block beyond the branch line slides on both faults,
  // so as a rigid block it can only move parallel to the branch line; the splay's slip vector is
  // the branch-line direction, scaled to a vertical component of one (STX-TN-05 section 4.1).
  for (auto& f : faults_) {
    if (!f.transfer || f.abuts < 0 || f.curved) continue;
    const Fault& m = faults_[static_cast<std::size_t>(f.abuts)];
    if (m.curved) continue;
    Vec3 dir = cross(m.n, f.n);
    double ln = norm(dir);
    if (ln < 1e-9)
      throw Error(ErrorCode::Build, "fault " + f.name + " (rule F15): a linked splay must not be parallel to its master", "faults");
    dir = dir / ln;
    if (dir.z > 0) dir = dir * -1.0;
    double plunge = std::asin(-dir.z);
    if (plunge < deg2rad(10))
      throw Error(ErrorCode::Build, "fault " + f.name + " (rule F15): the branch line with " + m.name + " plunges " +
                                        fmt_fixed(rad2deg(plunge), 1) + " degrees; a linked splay slips along it, which needs a plunge of "
                                        "at least 10 degrees. Declare the splay without transfer, or change its strike or dip",
                  "faults");
    // hanging wall moves down along the branch line for a normal splay
    f.slip = dir / (-dir.z) * (f.kin->slip_vector(f.dd, f.dip).z < 0 ? 1.0 : -1.0);
    log.value("fault." + f.name + ".slip_vector_per_m_throw", fmt_fixed(f.slip.x, 4) + ", " + fmt_fixed(f.slip.y, 4) + ", " + fmt_fixed(f.slip.z, 4));
  }
  for (auto& f : faults_)
    if (f.tri && ((fixed_block_ >> f.index) & 1u))
      throw Error(ErrorCode::Build, "fault " + f.name + " (rule F17): the frame's fixed point lies on the hanging-wall side of a "
                                                         "trishear fault; fix a point or well in the footwall",
                  "frame.fixed");
  for (auto& f : faults_)
    if (f.curved && ((fixed_block_ >> f.index) & 1u))
      throw Error(ErrorCode::Build,
                  "fault " + f.name + " (rule F11): the frame's fixed point lies in the hanging wall of a curved fault; the "
                  "shear construction moves the hanging wall over a fixed footwall, so fix a point or well in the footwall",
                  "frame.fixed");
  resolve_tips(d, st, log);
  for (auto& f : faults_) {
    if (!f.tri) continue;
    if (norm(trishear_map(f.index, fixed_, false) - fixed_) > 1e-9)
      throw Error(ErrorCode::Build, "fault " + f.name + " (rule F17): the frame's fixed point lies in the trishear zone's sweep; "
                                                         "fix a point or well outside it", "frame.fixed");
    double bu, bl, S = f.slip_at(0);
    f.tips(0, bu, bl);
    log.value(f.name + ".trishear.slip", fmt_fixed(S, 3));
    log.value(f.name + ".trishear.tip_present_depth", fmt_fixed(-f.at(0, bu).z, 3));
    log.value(f.name + ".trishear.tip_initial_depth", fmt_fixed(-f.at(0, bu + f.tri->ps * S).z, 3));
    // forelimb dip of each horizon at the reference strike position: march the restored horizon forward
    const Vec3 tip = f.at(0, bu);
    for (std::size_t hz = 0; hz < st.horizons(); ++hz) {
      if (st.z(hz, tip.x, tip.y) < tip.z) {  // below the present tip: cut by the fault, not folded over it
        log.value(f.name + "." + d.horizons[hz].name + ".trishear", "cut by the fault");
        continue;
      }
      double dipmax = 0;
      Vec3 prev{};
      bool have = false;
      for (double b = bu - 3000; b <= bu + f.tri->ps * S + 3000; b += 10) {
        Vec3 p = f.at(0, b);
        // the horizon in this fault's frame at the map position of the plane point (no older faults assumed)
        Vec3 h0{p.x, p.y, st.z(hz, p.x, p.y)};
        Vec3 pr = trishear_map(f.index, h0, true);
        if (have) {
          double du = std::hypot(pr.x - prev.x, pr.y - prev.y);
          if (du > 1e-6) dipmax = std::max(dipmax, std::atan(std::fabs(pr.z - prev.z) / du));
        }
        prev = pr;
        have = true;
      }
      log.value(f.name + "." + d.horizons[hz].name + ".forelimb_dip_max", fmt_fixed(rad2deg(dipmax), 3));
    }
  }
  // Branch lines of splays (STX-TN-05 section 4.1): the line along the cross product of the normals.
  for (auto& f : faults_) {
    if (f.abuts < 0 || f.curved) continue;
    const Fault& m = faults_[static_cast<std::size_t>(f.abuts)];
    if (m.curved) continue;
    Vec3 dir = cross(m.n, f.n);
    double ln = norm(dir);
    if (ln < 1e-9) continue;  // parallel planes: no branch line
    dir = dir / ln;
    if (dir.z > 0) dir = dir * -1.0;
    // the point of the line in the plane through the splay's reference point normal to the line
    // solves n_m.x = n_m.p0_m, n_f.x = n_f.p0_f, dir.x = dir.p0_f (Cramer's rule)
    Vec3 r0{dot(m.n, m.p0), dot(f.n, f.p0), dot(dir, f.p0)};
    auto det3 = [](const Vec3& a, const Vec3& b, const Vec3& c) { return dot(a, cross(b, c)); };
    Vec3 c0{m.n.x, f.n.x, dir.x}, c1{m.n.y, f.n.y, dir.y}, c2{m.n.z, f.n.z, dir.z};
    double D = det3(c0, c1, c2);
    Vec3 pt{det3(r0, c1, c2) / D, det3(c0, r0, c2) / D, det3(c0, c1, r0) / D};
    std::string key = f.name + ".branch." + m.name;
    double trend = std::fmod(rad2deg(std::atan2(dir.x, dir.y)) + 360.0, 360.0);
    log.value(key + ".trend", fmt_fixed(trend, 4));
    log.value(key + ".plunge", fmt_fixed(rad2deg(std::asin(-dir.z)), 4));
    log.value(key + ".point", fmt_fixed(pt.x, 2) + ", " + fmt_fixed(pt.y, 2) + ", " + fmt_fixed(-pt.z, 2));
    if (f.transfer) log.note("fault " + f.name + " is a linked splay of " + m.name + " (throw transfer across the branch line)");
  }
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
      q = step(fi, q, faults_[static_cast<std::size_t>(fi)].sigma(q) > 0, false);
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
  if (f.transfer || f.tri) return false;  // strike component of slip; trishear moves both sides
  // a point beyond the strike ends stays there unless the slip carries it along strike
  if (std::fabs(a) >= 0.5 * f.length + std::fabs(dot(f.slip, f.s)) * std::fabs(f.tmax)) return true;
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
      // along strike: the slip's strike component samples the strike profile elsewhere
      const double per_tau_a = std::fabs(dot(f.slip, f.s));
      if (per_tau_a > 1e-12)
        for (std::size_t i = 1; i < K; ++i) {
          double bu, bl;
          f.tips(f.a_tab[i], bu, bl);
          if (!(bl > bu) || !std::isfinite(bu) || !std::isfinite(bl)) continue;
          double bm = 0.5 * (bu + bl), da = f.a_tab[i] - f.a_tab[i - 1];
          double t1 = f.throw_at(f.a_tab[i], bm), t0 = f.throw_at(f.a_tab[i - 1], bm);
          gmax = std::max(gmax, std::fabs(t1 - t0) / da * per_tau_a);
        }
      f.max_slip_gradient = gmax;
      log.value("fault." + f.name + ".max_throw_gradient_along_slip", fmt_fixed(gmax, 3));
      if (f.tri) gmax = 0;  // constant slip behind the tip; the zone conserves area
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
    if (std::fabs(dot(f.slip, f.s)) > 1e-12 || f.slip.z > 0 || f.tmax < 0) {
      Vec3 u = normalized(f.slip) * (f.tmax < 0 ? -1.0 : 1.0);  // hanging-wall motion per metre of slip
      log.value(f.name + ".slip_vector", fmt_fixed(u.x, 4) + ", " + fmt_fixed(u.y, 4) + ", " + fmt_fixed(u.z, 4));
    }
    if (f.tmax < 0 && !f.layered) {
      // Reverse faults: repeated thickness of each horizon at the reference strike position
      // (STX-TN-05 section 5.2); for a plane it equals the throw where the horizon meets the fault.
      for (std::size_t hz = 0; hz < st.horizons(); ++hz) {
        auto g = [&](double b) {
          Vec3 p = f.at(0, b);
          double zl = older_exist ? level_in_frame(st, static_cast<int>(hz), p.x, p.y, k, dom) : st.z(hz, p.x, p.y);
          return p.z - zl;
        };
        double lo = b_top - margin, hi = b_base + margin, glo = g(lo), ghi = g(hi);
        if (!(glo > 0 && ghi < 0)) continue;
        double bh = illinois(g, lo, hi, glo, ghi, 1e-7);
        if (bh < bu0 || bh > bl0) continue;
        log.value(f.name + "." + d.horizons[hz].name + ".repeat", fmt_fixed(std::fabs(f.throw_at(0, bh)), 3));
      }
    }
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
  double u_top = c.u_at(-dom.z_top), u_far = u_top, u_near = 1e300;
  for (double cx : {dom.x0, dom.x0 + dom.sx})
    for (double cy : {dom.y0, dom.y0 + dom.sy}) {
      double u = dot(Vec3{cx, cy, f.p0.z} - f.p0, f.dh);
      u_far = std::max(u_far, u);
      u_near = std::min(u_near, u);
    }
  if (!std::isfinite(u_top)) u_top = std::min(c.pu.front(), u_near) - 1000 - std::fabs(c.heave);  // upper flat below the top
  if (!std::isfinite(u_far) || u_far < u_top) u_far = u_top;
  for (double cx : {dom.x0, dom.x0 + dom.sx})
    for (double cy : {dom.y0, dom.y0 + dom.sy}) u_far = std::max(u_far, dot(Vec3{cx, cy, f.p0.z} - f.p0, f.dh));
  u_far += 1000 + std::fabs(c.heave);
  if (c.fpf) {
    // Rule F18: flow lines of neighbouring bends cross above the bisector-crossing offset.
    double need = 0;
    for (double u = u_near; u <= u_far; u += 25) need = std::max(need, c.zf(u) + dom.z_top);
    log.value("fault." + f.name + ".fpf.bisector_crossing_offset",
              std::isfinite(c.fpf_dmax) ? fmt_fixed(c.fpf_dmax, 2) : std::string("none"));
    if (c.fpf_dmax < need)
      throw Error(ErrorCode::Build, "fault " + f.name + " (rule F18): the bisectors of neighbouring bends cross " +
                                        fmt_fixed(c.fpf_dmax, 1) + " m above the fault, below the model top (" + fmt_fixed(need, 1) +
                                        " m above it): fault-parallel flow lines would overlap. Lengthen the short segments "
                                        "between bends or soften the bends", "faults." + f.name + ".geometry.profile");
  }
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
  if (c.fpf) {
    log.value(f.name + ".fpf.slip", fmt_fixed(-h, 3));
    // per horizon: its offset above the fault (constant along its flow line)
    for (std::size_t hz = 0; hz < st.horizons(); ++hz) {
      Vec3 p = f.p0 + f.dh * c.pu.back();
      double z = -st.z(hz, p.x, p.y), s0, d0;
      c.fpf_coords(c.pu.back(), z, s0, d0);
      if (d0 <= 0) continue;
      log.value(f.name + "." + d.horizons[hz].name + ".fpf_offset_above_fault", fmt_fixed(d0, 3));
    }
    return;
  }
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
  if (f.tri) return trishear_map(fi, q, true) - q;
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
  if (f.tri) {
    Vec3 r = trishear_map(fi, q, false);
    if (tau_out) *tau_out = q.z - r.z;
    (void)sgn;
    return r;
  }
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

Vec3 FaultNetwork::drag_forward(int fi, const Vec3& p1, bool hw) const {
  const Fault& f = faults_[static_cast<std::size_t>(fi)];
  if (!f.drag) return {0, 0, 0};
  double w = f.drag_shape(f.sigma(p1));
  if (w == 0) return {0, 0, 0};
  if (f.abuts >= 0 && (faults_[static_cast<std::size_t>(f.abuts)].sigma(p1) > 0) != f.abut_hw_side) return {0, 0, 0};
  double a, b;
  f.profile_coords(p1, a, b);
  // normal drag: the hanging wall moves against its slip, the footwall with it
  double sg = (hw ? -1.0 : 1.0) * f.drag;
  return f.slip * (sg * w * f.throw_at(a, b));
}

Vec3 FaultNetwork::drag_inverse(int fi, const Vec3& q, bool hw) const {
  const Fault& f = faults_[static_cast<std::size_t>(fi)];
  if (!f.drag) return q;
  double w = f.drag_shape(f.sigma(q));  // sigma is unchanged by a displacement along the slip
  if (w == 0) return q;
  if (f.abuts >= 0 && (faults_[static_cast<std::size_t>(f.abuts)].sigma(q) > 0) != f.abut_hw_side) return q;
  const double c = (hw ? -1.0 : 1.0) * f.drag * w;
  // q = q1 + slip * c * T(q1), q1 = q - slip * t: solve t = c * T(q - slip t); the throw gradient
  // along the slip is below 1 and |c| <= 1, so the root is unique.
  auto g = [&](double t) {
    double a, b;
    f.profile_coords(q - f.slip * t, a, b);
    return t - c * f.throw_at(a, b);
  };
  double lim = std::fabs(c * f.tmax) + 1e-9;
  double glo = g(-lim), ghi = g(lim);
  double t = (glo == 0) ? -lim : (ghi == 0) ? lim : illinois(g, -lim, lim, glo, ghi, 1e-11, 1e-11, 200);
  return q - f.slip * t;
}

Vec3 FaultNetwork::trishear_map(int fi, const Vec3& q, bool fwd) const {
  const Fault& f = faults_[static_cast<std::size_t>(fi)];
  const TrishearFold& t = *f.tri;
  const Vec3 d = q - f.p0;
  const double a = dot(d, f.s), S = f.slip_at(a);
  if (S <= 0) return q;
  double bu, bl;
  f.tips(a, bu, bl);
  const double b0 = dot(d, f.dd), y0 = dot(d, f.n), tp = t.tanphi, ps = t.ps, v0 = t.v0;
  auto tipb = [&](double s) { return bu + ps * (S - s); };
  auto in_zone = [&](double b, double y, double s) { double x = tipb(s) - b; return x > 0 && std::fabs(y) < x * tp; };
  // Closed form when the relative path (a straight line along the fault) never meets the zone.
  {
    const double sa = fwd ? 0.0 : S, sb = fwd ? S : 0.0;
    const double rate = y0 > 0 ? v0 : 0.0;  // db/ds = -rate outside the zone
    auto b_at = [&](double s) { return b0 - rate * (s - sa); };
    if (!in_zone(b_at(sa), y0, sa) && !in_zone(b_at(sb), y0, sb)) return q + f.dd * (b_at(sb) - b0);
  }
  // Integrate the velocity field (fourth-order Runge-Kutta, steps shrinking towards the apex).
  auto vel = [&](double b, double y, double s, double& db, double& dy) {
    double x = tipb(s) - b;
    if (x <= 0 || y >= x * tp) {
      db = y > 0 || x > 0 ? -v0 : 0.0;
      if (x <= 0 && y <= 0) db = 0;
      dy = 0;
      return;
    }
    if (y <= -x * tp) { db = 0; dy = 0; return; }
    double r = y / (x * tp);
    db = -0.5 * v0 * (r + 1);
    dy = 0.25 * v0 * tp * (r * r - 1);
  };
  double b = b0, y = y0, s = fwd ? 0.0 : S;
  const double send = fwd ? S : 0.0, dir = fwd ? 1.0 : -1.0;
  while (dir * (send - s) > 1e-12) {
    double x = std::fabs(tipb(s) - b);
    double h = std::clamp(0.02 * x * tp + 0.05, 0.05, 5.0);
    h = std::min(h, dir * (send - s));
    double hs = dir * h;
    double k1b, k1y, k2b, k2y, k3b, k3y, k4b, k4y;
    vel(b, y, s, k1b, k1y);
    vel(b + 0.5 * hs * k1b, y + 0.5 * hs * k1y, s + 0.5 * hs, k2b, k2y);
    vel(b + 0.5 * hs * k2b, y + 0.5 * hs * k2y, s + 0.5 * hs, k3b, k3y);
    vel(b + hs * k3b, y + hs * k3y, s + hs, k4b, k4y);
    b += hs / 6 * (k1b + 2 * k2b + 2 * k3b + k4b);
    y += hs / 6 * (k1y + 2 * k2y + 2 * k3y + k4y);
    s += hs;
  }
  return q + f.dd * (b - b0) + f.n * (y - y0);
}

Vec3 FaultNetwork::step(int fi, const Vec3& q, bool side, bool fast) const {
  const Fault& f = faults_[static_cast<std::size_t>(fi)];
  if (f.tri) return trishear_map(fi, q, false);  // every side of the zone moves
  Vec3 r = f.drag ? drag_inverse(fi, q, side) : q;
  bool fixed_side = ((fixed_block_ >> fi) & 1u) != 0;
  if (side == fixed_side) return r;
  if (fast && cannot_move(fi, r)) return r;
  return undo(fi, r, side ? 1 : -1);
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
    q = step(fi, q, side, fast);
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
  Vec3 q = step(f, p, hw, false);
  for (std::size_t k = static_cast<std::size_t>(rank_[static_cast<std::size_t>(f)]) + 1; k < order_.size(); ++k) {
    int fi = order_[k];
    q = step(fi, q, faults_[static_cast<std::size_t>(fi)].sigma(q) > 0, false);
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
    q = step(fi, q, faults_[static_cast<std::size_t>(fi)].sigma(q) > 0, fast);
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
