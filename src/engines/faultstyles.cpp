// Stage 4B (STX-TN-15 section 2): polygonal fault networks, flower structures, duplexes, breached
// relays, trishear on curved faults and growth strata on every hanging wall.
#include "stratonix/engines/faultstyles.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "stratonix/core/error.hpp"

namespace stratonix {

std::vector<std::pair<double, double>> breached_relay_throw(const RelaySpec& s, int samples) {
  if (!(s.overlap > 0) || !(s.separation > 0) || samples < 2) throw Error(ErrorCode::InvalidArgument, "breached relay: positive overlap and separation");
  if (!(s.length_a > s.overlap && s.length_b > s.overlap)) throw Error(ErrorCode::InvalidArgument, "breached relay: the faults must be longer than their overlap");
  // fault A on [0, La], fault B on [La - overlap, La - overlap + Lb]; bell-shaped throw profiles
  // t (1 - (2x/L - 1)^2)^2. The breaching fault takes the deficit across the relay so that the
  // summed throw follows a smooth blend from A's maximum to B's (Peacock and Sanderson 1994)
  const double a0 = 0, a1 = s.length_a, b0 = s.length_a - s.overlap, b1 = b0 + s.length_b;
  auto bell = [](double x, double x0, double x1, double t) {
    if (x <= x0 || x >= x1) return 0.0;
    const double u = 2 * (x - x0) / (x1 - x0) - 1;
    return t * (1 - u * u) * (1 - u * u);
  };
  const double ca = 0.5 * (a0 + a1), cb = 0.5 * (b0 + b1);
  std::vector<std::pair<double, double>> out;
  for (int k = 0; k < samples; ++k) {
    const double x = a0 + (b1 - a0) * k / (samples - 1);
    double t;
    if (x <= ca) {
      t = bell(x, a0, a1, s.throw_a) + bell(x, b0, b1, s.throw_b);
      t = std::max(t, bell(x, a0, a1, s.throw_a));
    } else if (x >= cb) {
      t = std::max(bell(x, a0, a1, s.throw_a) + bell(x, b0, b1, s.throw_b), bell(x, b0, b1, s.throw_b));
    } else {
      // between the two maxima: a smooth step (C1) from A's maximum to B's
      const double u = (x - ca) / (cb - ca), w = u * u * (3 - 2 * u);
      t = (1 - w) * s.throw_a + w * s.throw_b;
    }
    out.push_back({x, t});
  }
  return out;
}

std::vector<Vec3> trishear_curved(const std::vector<Vec3>& fault_in, const TrishearParams& p, const std::vector<Vec3>& points) {
  // the trajectory in either order: from deep to shallow internally
  std::vector<Vec3> fault = fault_in;
  if (fault.size() >= 2 && fault.front().z > fault.back().z) std::reverse(fault.begin(), fault.end());
  if (fault.size() < 2) throw Error(ErrorCode::InvalidArgument, "trishear: a fault trajectory needs two points");
  if (!(p.apical_deg >= 20 && p.apical_deg <= 100) || !(p.ps_ratio >= 0 && p.ps_ratio <= 4) || !(p.increment > 0))
    throw Error(ErrorCode::InvalidArgument, "trishear: apical angle 20 to 100 degrees, propagation-to-slip 0 to 4, positive increment");
  // the trajectory, in (x, z), from deep to shallow; the tip starts ps x slip below its end and
  // propagates along it (Cardozo 2008), beyond its end along the last tangent
  std::vector<double> arc{0};
  for (std::size_t i = 1; i < fault.size(); ++i) arc.push_back(arc.back() + std::hypot(fault[i].x - fault[i - 1].x, fault[i].z - fault[i - 1].z));
  for (std::size_t i = 1; i < fault.size(); ++i) {
    if (!(arc[i] - arc[i - 1] > 0)) throw Error(ErrorCode::InvalidArgument, "trishear: repeated trajectory points");
  }
  auto at = [&](double sa, Vec3* tan) {
    std::size_t i = 1;
    while (i + 1 < fault.size() && arc[i] < sa) ++i;
    const double L = arc[i] - arc[i - 1];
    Vec3 t{(fault[i].x - fault[i - 1].x) / L, 0, (fault[i].z - fault[i - 1].z) / L};
    if (tan) *tan = t;
    const double u = sa - arc[i - 1];
    return Vec3{fault[i - 1].x + u * t.x, 0, fault[i - 1].z + u * t.z};
  };
  const double m = std::tan(0.5 * p.apical_deg * kPi / 180);
  double s_tip = std::max(0.0, arc.back() - p.ps_ratio * p.slip);
  std::vector<Vec3> out = points;
  const int n = std::max(1, static_cast<int>(std::ceil(p.slip / p.increment)));
  const double ds = p.slip / n;
  for (int it = 0; it < n; ++it) {
    Vec3 tt;
    const Vec3 tip = at(s_tip, &tt);
    const Vec3 nn{-tt.z, 0, tt.x};  // hanging-wall side: left of the up-dip tangent
    for (auto& q : out) {
      const double xi = (q.x - tip.x) * tt.x + (q.z - tip.z) * tt.z, eta = (q.x - tip.x) * nn.x + (q.z - tip.z) * nn.z;
      Vec3 v{0, 0, 0};
      if (xi > 0 && std::fabs(eta) <= m * xi) {
        // linear trishear velocity (Zehnder and Allmendinger 2000), v0 = 1
        const double r = eta / (m * xi);
        const double vx = 0.5 * (r + 1), vy = 0.25 * m * (r * r - 1);
        v = tt * vx + nn * vy;
      } else if (xi > 0) {
        if (eta > 0) v = tt;
      } else {
        // behind the tip: the side of the fault, and motion parallel to the fault where the point lies
        double best = 1e300, sb = 0;
        for (int k = 0; k <= 200; ++k) {
          const double sa = s_tip * k / 200.0;
          const Vec3 f = at(sa, nullptr);
          const double d = (q.x - f.x) * (q.x - f.x) + (q.z - f.z) * (q.z - f.z);
          if (d < best) { best = d; sb = sa; }
        }
        Vec3 tq;
        const Vec3 fq = at(sb, &tq);
        const Vec3 nq{-tq.z, 0, tq.x};
        if ((q.x - fq.x) * nq.x + (q.z - fq.z) * nq.z > 0) v = tq;
      }
      q.x += v.x * ds;
      q.z += v.z * ds;
    }
    s_tip += p.ps_ratio * ds;  // beyond the trajectory, at() extends it along its last tangent
  }
  return out;
}

std::vector<double> growth_thickness(const GrowthSpec& s) {
  static const char* kinds[] = {"layer_rigid", "rigid", "shear", "fault_bend", "fault_propagation", "trishear"};
  if (std::none_of(std::begin(kinds), std::end(kinds), [&](const char* k) { return s.hanging_wall == k; }))
    throw Error(ErrorCode::InvalidArgument, "growth strata: unknown hanging-wall type '" + s.hanging_wall + "'");
  if (s.deposition_times.size() < 2) throw Error(ErrorCode::InvalidArgument, "growth strata: at least two deposition times");
  // at the fault, every hanging-wall type gains accommodation equal to the vertical component of the
  // slip during the unit's deposition (Suppe et al. 1992); the types differ away from the fault
  std::vector<double> t;
  for (std::size_t k = 1; k < s.deposition_times.size(); ++k) {
    const double dt = s.deposition_times[k] - s.deposition_times[k - 1];
    if (!(dt > 0)) throw Error(ErrorCode::InvalidArgument, "growth strata: deposition times must increase");
    t.push_back(s.slip_rate * dt * std::sin(s.dip_deg * kPi / 180));
  }
  return t;
}

}  // namespace stratonix
