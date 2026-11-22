// Stage 4B (STX-TN-15 sections 2.1 to 2.3): planners of polygonal fault networks, flower
// structures and duplexes.
#include "stratonix/core/faultsystems.hpp"

#include <algorithm>
#include <cmath>

#include "stratonix/core/error.hpp"

namespace stratonix {

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
  double normal() {
    const double u1 = std::max(uniform(), 1e-300), u2 = uniform();
    return std::sqrt(-2 * std::log(u1)) * std::cos(2 * kPi * u2);
  }
};
}  // namespace

PolygonalNetwork polygonal_network(const PolygonalSpec& s, std::uint64_t seed) {
  if (!(s.spacing > 0) || !(s.x1 > s.x0) || !(s.y1 > s.y0)) throw Error(ErrorCode::InvalidArgument, "polygonal network: a positive spacing and a non-empty area");
  if (!(s.unit_top > s.unit_base)) throw Error(ErrorCode::InvalidArgument, "polygonal network: the unit's top must lie above its base");
  // seeds on a jittered lattice of the labelled spacing; one fault on the perpendicular bisector of
  // each pair of lattice neighbours (a quasi-Voronoi tessellation; Cartwright and Lonergan 1996)
  Rng g(seed);
  const int nx = std::max(2, static_cast<int>((s.x1 - s.x0) / s.spacing)), ny = std::max(2, static_cast<int>((s.y1 - s.y0) / s.spacing));
  std::vector<Vec3> seeds(static_cast<std::size_t>(nx * ny));
  const double jit = 0.25 * s.spacing;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      seeds[static_cast<std::size_t>(j * nx + i)] = {s.x0 + (i + 0.5) * s.spacing + jit * (2 * g.uniform() - 1), s.y0 + (j + 0.5) * s.spacing + jit * (2 * g.uniform() - 1), s.unit_top};
  PolygonalNetwork R;
  double sum = 0;
  auto add = [&](int a, int b) {
    const Vec3 pa = seeds[static_cast<std::size_t>(a)], pb = seeds[static_cast<std::size_t>(b)];
    const Vec3 m = (pa + pb) * 0.5, d = pb - pa;
    const double dist = std::hypot(d.x, d.y);
    const Vec3 t = normalized(Vec3{-d.y, d.x, 0});
    PolygonalFault f;
    f.length = s.spacing;
    f.a = m - t * (0.5 * f.length);
    f.b = m + t * (0.5 * f.length);
    f.dip_deg = std::clamp(s.dip_mean + s.dip_sd * g.normal(), 30.0, 70.0);
    f.throw_m = s.throw_ratio * f.length;  // displacement-length scaling (C10)
    R.faults.push_back(f);
    sum += dist;
  };
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      if (i + 1 < nx) add(j * nx + i, j * nx + i + 1);
      if (j + 1 < ny) add(j * nx + i, (j + 1) * nx + i);
    }
  R.mean_spacing = R.faults.empty() ? 0 : sum / static_cast<double>(R.faults.size());
  return R;
}

FlowerResult flower_splays(const FlowerSpec& s) {
  if (s.splays < 2 || s.splays > 9) throw Error(ErrorCode::InvalidArgument, "flower structure: from 2 to 9 splays");
  if (!(s.spread_deg >= 10 && s.spread_deg <= 60)) throw Error(ErrorCode::InvalidArgument, "flower structure: spread from 10 to 60 degrees");
  FlowerResult R;
  // splays fan symmetrically about the vertical master fault; each carries an equal share of the
  // vertical component and of the horizontal slip (Harding 1985; Sylvester 1988)
  for (int k = 0; k < s.splays; ++k) {
    const double f = s.splays == 1 ? 0.0 : static_cast<double>(k) / (s.splays - 1) - 0.5;
    R.dips_deg.push_back(90.0 - std::fabs(f) * s.spread_deg);
    R.vertical_throw.push_back(s.vertical_component / s.splays);
    R.horizontal_slip.push_back(s.horizontal_slip / s.splays);
  }
  return R;
}

DuplexResult duplex(const DuplexSpec& s) {
  if (s.horses < 1 || s.horses > 12) throw Error(ErrorCode::InvalidArgument, "duplex: from 1 to 12 horses");
  if (!(s.ramp_angle_deg >= 10 && s.ramp_angle_deg <= 45) || !(s.slip > 0) || !(s.horse_length > 0))
    throw Error(ErrorCode::InvalidArgument, "duplex: ramp angle from 10 to 45 degrees, positive slip and horse length");
  DuplexResult R;
  // each horse is carried over the next ramp by its slip (Boyer and Elliott 1982; Mitra 1986)
  R.shortening = s.horses * s.slip;
  const double r = s.slip / s.horse_length;
  if (r < 0.95) {
    R.style = "hinterland_dipping";
    R.stack_thickness = 2 * s.thickness;  // each horse overlaps the next one: two thicknesses
  } else if (r <= 1.05) {
    R.style = "antiformal_stack";
    R.stack_thickness = s.horses * s.thickness;  // the horses stack on one another
  } else {
    R.style = "foreland_dipping";
    R.stack_thickness = 2 * s.thickness;
  }
  return R;
}

}  // namespace stratonix
