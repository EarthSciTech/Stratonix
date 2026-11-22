// Stage 4B (STX-TN-15 section 5, decision 19): fracture intensity, discrete fracture networks,
// mechanical stratigraphy and laterally varying stress.
#include "stratonix/engines/fractures.hpp"

#include <algorithm>
#include <cmath>

#include "stratonix/core/error.hpp"

namespace stratonix {

double fracture_intensity(double k_max, double e_max, const FractureLaw& law) {
  if (law.p32_0 < 0 || law.a_k < 0 || law.a_e < 0 || !(law.p32_max > 0 && law.p32_max <= 10))
    throw Error(ErrorCode::InvalidArgument, "fracture law: coefficients >= 0 and a maximum P32 in (0, 10] 1/m");
  // P32 = P32_0 + a_k |k_max| + a_e e_max (Murray 1968; Fischer and Wilkerson 2000), clipped
  return std::clamp(law.p32_0 + law.a_k * std::fabs(k_max) + law.a_e * std::fabs(e_max), 0.0, law.p32_max);
}

namespace {
struct Rng {
  std::uint64_t s[4];
  explicit Rng(std::uint64_t seed) {
    for (auto& v : s) {
      seed += 0x9e3779b97f4a7c15ULL;
      std::uint64_t z = seed;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      v = z ^ (z >> 31);
    }
  }
  static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  std::uint64_t next() {
    const std::uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return r;
  }
  double uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }
};
/// Disk area inside a box, by a polar quadrature of the disk (16 rings x 48 sectors).
double area_inside(const Disk& d, const BBox& box) {
  const Vec3 a = normalized(std::fabs(d.n.z) < 0.9 ? cross(d.n, Vec3{0, 0, 1}) : cross(d.n, Vec3{1, 0, 0}));
  const Vec3 b = cross(d.n, a);
  auto in = [&](const Vec3& p) { return p.x >= box.lo.x && p.x <= box.hi.x && p.y >= box.lo.y && p.y <= box.hi.y && p.z >= box.lo.z && p.z <= box.hi.z; };
  // whole disk inside or wholly outside its bounding sphere test first
  const bool inside_all = d.c.x - d.r >= box.lo.x && d.c.x + d.r <= box.hi.x && d.c.y - d.r >= box.lo.y && d.c.y + d.r <= box.hi.y &&
                          d.c.z - d.r >= box.lo.z && d.c.z + d.r <= box.hi.z;
  if (inside_all) return kPi * d.r * d.r;
  if (d.c.x + d.r < box.lo.x || d.c.x - d.r > box.hi.x || d.c.y + d.r < box.lo.y || d.c.y - d.r > box.hi.y || d.c.z + d.r < box.lo.z ||
      d.c.z - d.r > box.hi.z)
    return 0.0;
  const int nr = 16, nt = 48;
  double area = 0;
  for (int i = 0; i < nr; ++i) {
    const double r0 = d.r * i / nr, r1 = d.r * (i + 1) / nr, rm = 0.5 * (r0 + r1), ring = kPi * (r1 * r1 - r0 * r0) / nt;
    for (int j = 0; j < nt; ++j) {
      const double t = 2 * kPi * (j + 0.5) / nt;
      if (in(d.c + a * (rm * std::cos(t)) + b * (rm * std::sin(t)))) area += ring;
    }
  }
  return area;
}
/// Fisher-distributed unit vector about `mean` (Fisher 1953).
Vec3 fisher(const Vec3& mean, double kappa, Rng& g) {
  const double u = g.uniform(), phi = 2 * kPi * g.uniform();
  const double w = kappa > 1e-9 ? 1 + std::log(u + (1 - u) * std::exp(-2 * kappa)) / kappa : 2 * u - 1;  // cos(theta)
  const double s = std::sqrt(std::max(0.0, 1 - w * w));
  const Vec3 m = normalized(mean);
  const Vec3 a = normalized(std::fabs(m.z) < 0.9 ? cross(m, Vec3{0, 0, 1}) : cross(m, Vec3{1, 0, 0}));
  const Vec3 b = cross(m, a);
  return normalized(m * w + a * (s * std::cos(phi)) + b * (s * std::sin(phi)));
}
/// Truncated power law of radii, pdf ~ r^-alpha on [rmin, rmax] (inverse transform).
double power_law(double alpha, double rmin, double rmax, double u) {
  const double e = 1 - alpha;
  return std::pow(std::pow(rmin, e) + u * (std::pow(rmax, e) - std::pow(rmin, e)), 1 / e);
}
}  // namespace

DfnResult sample_dfn(const BBox& box, double p32, const std::vector<FractureSet>& sets, std::uint64_t seed, std::size_t limit) {
  DfnResult R;
  if (sets.empty() || !(p32 >= 0 && p32 <= 10)) { R.reason = "a fracture set and a P32 from 0 to 10 1/m"; return R; }
  double wsum = 0, mean_area = 0;
  for (auto& s : sets) {
    if (!(s.kappa >= 0) || !(s.exponent > 2) || !(s.rmin > 0 && s.rmin < s.rmax) || !(s.weight > 0)) {
      R.reason = "fracture sets: kappa >= 0, power-law exponent > 2, 0 < minimum radius < maximum radius, positive weight";
      return R;
    }
    wsum += s.weight;
  }
  for (auto& s : sets) {
    // E[pi r^2] of the truncated power law
    const double a = s.exponent;
    const double norm_c = (std::pow(s.rmax, 1 - a) - std::pow(s.rmin, 1 - a)) / (1 - a);
    const double m2 = (std::pow(s.rmax, 3 - a) - std::pow(s.rmin, 3 - a)) / (3 - a);
    mean_area += s.weight / wsum * kPi * m2 / norm_c;
  }
  const double V = (box.hi.x - box.lo.x) * (box.hi.y - box.lo.y) * (box.hi.z - box.lo.z);
  if (!(V > 0)) { R.reason = "an empty box"; return R; }
  const double expected = p32 * V / mean_area;
  if (expected > static_cast<double>(limit)) {
    R.reason = "the network would hold about " + std::to_string(static_cast<long long>(expected)) + " fractures, above the limit of " + std::to_string(limit);
    return R;
  }
  // Baecher disks (Baecher et al. 1977): centres uniform in the box enlarged by the largest radius,
  // drawn until the area inside the box reaches the target P32 x V (the realised intensity is the
  // target, the count follows)
  double rmax = 0;
  for (auto& s : sets) rmax = std::max(rmax, s.rmax);
  Rng g(seed);
  const double target = p32 * V;
  double area = 0;
  std::size_t guard = 0;
  while (area < target) {
    if (++guard > 50 * limit) { R.reason = "the network did not reach its intensity"; R.disks.clear(); return R; }
    double u = g.uniform() * wsum;
    int k = 0;
    while (k + 1 < static_cast<int>(sets.size()) && u > sets[static_cast<std::size_t>(k)].weight) { u -= sets[static_cast<std::size_t>(k)].weight; ++k; }
    const auto& s = sets[static_cast<std::size_t>(k)];
    Disk d;
    d.set = k;
    d.c = {box.lo.x - rmax + (box.hi.x - box.lo.x + 2 * rmax) * g.uniform(), box.lo.y - rmax + (box.hi.y - box.lo.y + 2 * rmax) * g.uniform(),
           box.lo.z - rmax + (box.hi.z - box.lo.z + 2 * rmax) * g.uniform()};
    d.n = fisher(s.mean_pole, s.kappa, g);
    d.r = power_law(s.exponent, s.rmin, s.rmax, g.uniform());
    const double ai = area_inside(d, box);
    if (ai <= 0) continue;
    area += ai;
    R.disks.push_back(d);
    if (R.disks.size() > limit) { R.reason = "above the limit of " + std::to_string(limit) + " fractures"; R.disks.clear(); return R; }
  }
  R.ok = true;
  return R;
}

double realised_p32(const std::vector<Disk>& d, const BBox& box) {
  const double V = (box.hi.x - box.lo.x) * (box.hi.y - box.lo.y) * (box.hi.z - box.lo.z);
  double a = 0;
  for (auto& x : d) a += area_inside(x, box);
  return V > 0 ? a / V : 0.0;
}

std::pair<double, double> horizontal_stress(double sv, double pore, double alpha, double E, double nu, double e_h, double e_H) {
  if (!(nu >= 0 && nu < 0.5) || !(E > 0)) throw Error(ErrorCode::InvalidArgument, "stress: Poisson's ratio in [0, 0.5), positive Young's modulus");
  // uniaxial strain plus tectonic strain (Thiercelin and Plumb 1994)
  const double base = nu / (1 - nu) * (sv - alpha * pore) + alpha * pore, c = E / (1 - nu * nu);
  return {base + c * (e_h + nu * e_H), base + c * (e_H + nu * e_h)};
}

}  // namespace stratonix
