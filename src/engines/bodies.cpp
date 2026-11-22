// Stage 4B (STX-TN-15 section 1): salt sheets, sills, saucer sills, dykes and sand injectites.
#include "stratonix/engines/bodies.hpp"

#include <algorithm>
#include <cmath>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"

namespace stratonix {

namespace {
bool inside_polygon(const std::vector<Vec3>& poly, double x, double y) {
  bool in = false;
  for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
    const auto& a = poly[i];
    const auto& b = poly[j];
    if ((a.y > y) != (b.y > y) && x < (b.x - a.x) * (y - a.y) / (b.y - a.y) + a.x) in = !in;
  }
  return in;
}
double distance_to_edges(const std::vector<Vec3>& poly, double x, double y) {
  double best = 1e300;
  for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
    const double ex = poly[i].x - poly[j].x, ey = poly[i].y - poly[j].y, l2 = ex * ex + ey * ey;
    const double t = l2 > 0 ? std::clamp(((x - poly[j].x) * ex + (y - poly[j].y) * ey) / l2, 0.0, 1.0) : 0.0;
    best = std::min(best, std::hypot(poly[j].x + t * ex - x, poly[j].y + t * ey - y));
  }
  return best;
}

// SplitMix64 / xoshiro256** (engines do not use the package generator)
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

double SheetBody::thickness_at(double x, double y) const {
  if (footprint.size() < 3 || !inside_polygon(footprint, x, y)) return 0;
  if (!(taper > 0)) return thickness;
  // linear taper with the distance to the footprint edge
  return thickness * std::min(1.0, distance_to_edges(footprint, x, y) / taper);
}

double sheet_volume(const SheetBody& s, double spacing) {
  if (!(spacing > 0)) throw Error(ErrorCode::InvalidArgument, "sheet volume: the spacing must be positive");
  BBox b;
  for (auto& p : s.footprint) b.add(p);
  // cell-centred midpoint integration; sub-sampled 4 x 4 per cell along the taper for accuracy
  double v = 0;
  const std::size_t nx = static_cast<std::size_t>(std::ceil((b.hi.x - b.lo.x) / spacing)), ny = static_cast<std::size_t>(std::ceil((b.hi.y - b.lo.y) / spacing));
  for (std::size_t j = 0; j < ny; ++j)
    for (std::size_t i = 0; i < nx; ++i) {
      const double x0 = b.lo.x + static_cast<double>(i) * spacing, y0 = b.lo.y + static_cast<double>(j) * spacing;
      double sum = 0;
      for (int a = 0; a < 4; ++a)
        for (int c = 0; c < 4; ++c) sum += s.thickness_at(x0 + (a + 0.5) * spacing / 4, y0 + (c + 0.5) * spacing / 4);
      v += sum / 16 * spacing * spacing;
    }
  return v;
}

double SillBody::thickness_at(double r) const {
  if (!(radius > 0) || r >= radius) return 0;
  return thickness * std::sqrt(1 - r * r / (radius * radius));
}

double sill_volume(const SillBody& s) { return 2.0 / 3.0 * kPi * s.radius * s.radius * s.thickness; }

double forced_fold_uplift(const SillBody& s, double r) {
  if (!(s.radius > 0) || r >= s.radius) return 0;
  // fold volume (1/3) pi a^2 w0 = f (2/3) pi a^2 h0, so w0 = 2 f h0 (Pollard and Johnson 1973)
  const double w0 = 2 * std::clamp(s.accommodation, 0.0, 1.0) * s.thickness, q = 1 - r * r / (s.radius * s.radius);
  return w0 * q * q;
}

double forced_fold_volume(const SillBody& s, double spacing) {
  if (!(spacing > 0)) throw Error(ErrorCode::InvalidArgument, "forced fold volume: the spacing must be positive");
  const int n = static_cast<int>(std::ceil(s.radius / spacing));
  double v = 0;
  for (int j = -n; j < n; ++j)
    for (int i = -n; i < n; ++i) {
      double sum = 0;
      for (int a = 0; a < 2; ++a)
        for (int c = 0; c < 2; ++c)
          sum += forced_fold_uplift(s, std::hypot((i + (a + 0.5) / 2) * spacing, (j + (c + 0.5) / 2) * spacing));
      v += sum / 4 * spacing * spacing;
    }
  return v;
}

InjectiteComplex build_injectites(const InjectiteSpec& s, std::uint64_t seed) {
  InjectiteComplex R;
  if (!(s.ceiling > s.parent_top)) { R.reason = "the ceiling must lie above (be younger than) the parent unit's top"; return R; }
  if (s.sheets < 1 || !(s.parent_sand_volume > 0) || !(s.max_fraction > 0 && s.max_fraction <= 1)) {
    R.reason = "at least one sheet, a positive parent sand volume and a fraction in (0, 1]";
    return R;
  }
  if (!(s.area.hi.x > s.area.lo.x && s.area.hi.y > s.area.lo.y)) { R.reason = "the placement area is empty"; return R; }
  Rng g(seed);
  const double climb = s.climb_deg * kPi / 180, rise = s.ceiling - s.parent_top;
  for (int k = 0; k < s.sheets; ++k) {
    const double x = s.area.lo.x + g.uniform() * (s.area.hi.x - s.area.lo.x), y = s.area.lo.y + g.uniform() * (s.area.hi.y - s.area.lo.y);
    const double az = g.uniform() * 2 * kPi;
    const double len = std::max(1.0, s.length_mean + s.length_sd * g.normal());
    const double th = std::max(0.1, s.thickness_mean + s.thickness_sd * g.normal());
    // a wing from the parent top, climbing at the labelled angle, capped by the ceiling
    const double run = std::min(len * std::cos(climb), rise / std::max(1e-9, std::tan(climb)));
    const double up = std::min(rise, run * std::tan(climb));
    InjectiteSheet sh;
    sh.start = {x, y, s.parent_top};
    sh.end = {x + run * std::sin(az), y + run * std::cos(az), s.parent_top + up};
    sh.thickness = th;
    sh.width = 0.5 * len;  // labelled aspect: width half the length
    R.sheets.push_back(sh);
  }
  double vol = 0;
  for (auto& sh : R.sheets) vol += norm(sh.end - sh.start) * sh.width * sh.thickness;
  // the parent cannot give more than the labelled fraction: scale the thicknesses down to it
  const double cap = s.max_fraction * s.parent_sand_volume;
  if (vol > cap) {
    const double f = cap / vol;
    for (auto& sh : R.sheets) sh.thickness *= f;
    vol = cap;
  }
  R.injected_volume = vol;
  R.parent_volume_after = s.parent_sand_volume - vol;
  R.ok = true;
  return R;
}

// ---- bodies in a model ------------------------------------------------------------------------

namespace {
struct Placed {
  BodyReport rep;
  std::string type;
  // salt sheet
  SheetBody sheet;
  GridSpec base_grid;
  std::vector<double> base;  ///< base elevation on base_grid (NaN outside)
  // sill
  SillBody sill;
  double inner_radius = 0, outer_radius = 0, climb = 0, decay = 0;
  // dyke
  std::vector<Vec3> trace;
  double width = 0, top = 0, bottom = 0;
  // injectite
  std::vector<InjectiteSheet> sheets;

  double base_at(double x, double y) const {
    if (base.empty()) return sheet.base_z;
    const double fx = std::clamp((x - base_grid.x0) / base_grid.dx, 0.0, static_cast<double>(base_grid.nx - 1));
    const double fy = std::clamp((y - base_grid.y0) / base_grid.dy, 0.0, static_cast<double>(base_grid.ny - 1));
    const std::size_t i = std::min<std::size_t>(static_cast<std::size_t>(fx), base_grid.nx - 2), j = std::min<std::size_t>(static_cast<std::size_t>(fy), base_grid.ny - 2);
    const double u = fx - static_cast<double>(i), v = fy - static_cast<double>(j);
    auto at = [&](std::size_t a, std::size_t b) { return base[b * base_grid.nx + a]; };
    return (1 - u) * (1 - v) * at(i, j) + u * (1 - v) * at(i + 1, j) + (1 - u) * v * at(i, j + 1) + u * v * at(i + 1, j + 1);
  }
  // mid-plane, thickness of a (saucer) sill at radius r
  double sill_mid(double r) const { return sill.centre.z + (r > inner_radius ? (r - inner_radius) * std::tan(climb) : 0.0); }
  double sill_thick(double r) const {
    if (r >= outer_radius) return 0;
    return sill.thickness * std::sqrt(1 - r * r / (outer_radius * outer_radius));
  }
  bool contains(const Vec3& p) const {
    if (type == "salt_sheet") {
      const double t = sheet.thickness_at(p.x, p.y);
      if (!(t > 0)) return false;
      const double b = base_at(p.x, p.y);
      return p.z >= b && p.z <= b + t;
    }
    if (type == "sill" || type == "saucer_sill") {
      const double r = std::hypot(p.x - sill.centre.x, p.y - sill.centre.y), h = sill_thick(r);
      return h > 0 && std::fabs(p.z - sill_mid(r)) <= 0.5 * h;
    }
    if (type == "dyke") {
      if (p.z > top || p.z < bottom) return false;
      for (std::size_t i = 1; i < trace.size(); ++i) {
        const double ex = trace[i].x - trace[i - 1].x, ey = trace[i].y - trace[i - 1].y, l2 = ex * ex + ey * ey;
        const double t = l2 > 0 ? std::clamp(((p.x - trace[i - 1].x) * ex + (p.y - trace[i - 1].y) * ey) / l2, 0.0, 1.0) : 0.0;
        if (std::hypot(trace[i - 1].x + t * ex - p.x, trace[i - 1].y + t * ey - p.y) <= 0.5 * width) return true;
      }
      return false;
    }
    if (type == "injectite") {
      for (auto& sh : sheets) {
        const Vec3 d = sh.end - sh.start;
        const double L = norm(d);
        if (!(L > 0)) continue;
        const Vec3 u = d / L;
        Vec3 lat = normalized(Vec3{-u.y, u.x, 0});
        if (!(norm(lat) > 0)) lat = {1, 0, 0};
        const Vec3 n = normalized(cross(u, lat));
        const Vec3 q = p - sh.start;
        const double a = dot(q, u), b = dot(q, lat), c = dot(q, n);
        if (a >= 0 && a <= L && std::fabs(b) <= 0.5 * sh.width && std::fabs(c) <= 0.5 * sh.thickness) return true;
      }
      return false;
    }
    return false;
  }
  /// uplift of the host at p (roof of a sill, overburden of a salt sheet)
  double uplift(const Vec3& p) const {
    if (type == "salt_sheet") {
      const double t = sheet.thickness_at(p.x, p.y);
      return (t > 0 && p.z > base_at(p.x, p.y) + t) ? t : 0.0;
    }
    if (type == "sill" || type == "saucer_sill") {
      const double r = std::hypot(p.x - sill.centre.x, p.y - sill.centre.y), h = sill_thick(r);
      if (!(h > 0)) return 0;
      const double roof = sill_mid(r) + 0.5 * h;
      if (p.z <= roof) return 0;
      // the roof is lifted by f h(r): no gap and no overlap; the lift decays upward per km
      const double w = std::clamp(sill.accommodation, 0.0, 1.0) * h;
      return w * std::max(0.0, 1 - decay * (p.z - roof) / 1000.0);
    }
    return 0;
  }
};

/// Elevation of horizon h at (x, y): where the unit changes from h to h + 1 downward (NaN if none).
double horizon_z(const Model& m, int h, double x, double y) {
  const auto& dom = m.info().domain;
  const int n = 600;
  std::vector<Vec3> pts(n);
  for (int k = 0; k < n; ++k) pts[static_cast<std::size_t>(k)] = {x, y, dom.z_top - (dom.z_top - dom.z_base) * (k + 0.5) / n};
  EvalResult er;
  m.evaluate(pts, kEvalUnit, er);
  for (int k = 1; k < n; ++k) {
    if (er.unit[static_cast<std::size_t>(k - 1)] <= h && er.unit[static_cast<std::size_t>(k)] > h) {
      double a = pts[static_cast<std::size_t>(k - 1)].z, b = pts[static_cast<std::size_t>(k)].z;
      for (int it = 0; it < 30; ++it) {
        const double zm = 0.5 * (a + b);
        (m.unit_at({x, y, zm}) <= h ? a : b) = zm;
      }
      return 0.5 * (a + b);
    }
  }
  return kNaN;
}

class BodyModel final : public Model {
 public:
  BodyModel(ModelPtr inner, std::vector<Placed> bodies) : inner_(std::move(inner)), b_(std::move(bodies)), info_(inner_->info()) {
    for (auto& b : b_) {
      UnitInfo u;
      u.name = b.rep.name;
      info_.units.push_back(u);
      b.rep.unit = static_cast<int>(info_.units.size()) - 1;
    }
    // youngest first: a younger body occupies a point before an older one
    std::stable_sort(b_.begin(), b_.end(), [](const Placed& a, const Placed& c) { return a.rep.age > c.rep.age; });
    for (auto& b : b_) reports_.push_back(b.rep);
  }
  const ModelInfo& info() const override { return info_; }
  void evaluate(std::span<const Vec3> pts, unsigned fields, EvalResult& out) const override {
    std::vector<Vec3> q(pts.begin(), pts.end());
    for (auto& p : q) p.z -= total_uplift(p);
    inner_->evaluate(q, fields, out);
    if (fields & kEvalUnit)
      for (std::size_t i = 0; i < pts.size(); ++i)
        for (auto& b : b_)
          if (b.contains(pts[i])) { out.unit[i] = b.rep.unit; break; }
  }
  void scalar_in_block(std::span<const Vec3> pts, BlockCode block, int horizon, std::span<double> out) const override {
    std::vector<Vec3> q(pts.begin(), pts.end());
    for (auto& p : q) p.z -= total_uplift(p);
    inner_->scalar_in_block(q, block, horizon, out);
  }
  void evaluate_fault(int fault, std::span<const Vec3> pts, FaultEval& out) const override { inner_->evaluate_fault(fault, pts, out); }
  void restore(std::span<const Vec3> pts, std::span<Vec3> out) const override {
    // restoration removes the bodies first: the host at the undeformed position
    std::vector<Vec3> q(pts.begin(), pts.end());
    for (auto& p : q) p.z -= total_uplift(p);
    inner_->restore(q, out);
  }
  const void* capability(std::string_view name) const override {
    if (name == "bodies") return &reports_;
    return inner_->capability(name);
  }

 private:
  double total_uplift(const Vec3& p) const {
    double u = 0;
    for (auto& b : b_) u += b.uplift(p);
    return u;
  }
  ModelPtr inner_;
  std::vector<Placed> b_;
  ModelInfo info_;
  std::vector<BodyReport> reports_;
};

std::vector<Vec3> polygon_of(const desc::Param& p) {
  std::vector<Vec3> v;
  for (std::size_t r = 0; r < p.rows(); ++r) v.push_back({p.si(2 * r), p.si(2 * r + 1), 0});
  return v;
}
/// sample points of a body, to test the faults that cross it
std::vector<Vec3> samples_of(const Placed& b) {
  std::vector<Vec3> s;
  if (b.type == "salt_sheet") {
    BBox bb;
    for (auto& p : b.sheet.footprint) bb.add(p);
    for (int j = 0; j <= 20; ++j)
      for (int i = 0; i <= 20; ++i) {
        const double x = bb.lo.x + (bb.hi.x - bb.lo.x) * i / 20, y = bb.lo.y + (bb.hi.y - bb.lo.y) * j / 20, t = b.sheet.thickness_at(x, y);
        if (t > 0) s.push_back({x, y, b.base_at(x, y) + 0.5 * t});
      }
  } else if (b.type == "sill" || b.type == "saucer_sill") {
    for (int k = 0; k < 24; ++k)
      for (double f : {0.0, 0.5, 0.95}) {
        const double r = f * b.outer_radius, a = 2 * kPi * k / 24;
        s.push_back({b.sill.centre.x + r * std::cos(a), b.sill.centre.y + r * std::sin(a), b.sill_mid(r)});
      }
  } else if (b.type == "dyke") {
    for (std::size_t i = 1; i < b.trace.size(); ++i)
      for (int k = 0; k <= 10; ++k)
        for (double z : {b.bottom, 0.5 * (b.top + b.bottom), b.top}) {
          const Vec3 p = b.trace[i - 1] + (b.trace[i] - b.trace[i - 1]) * (k / 10.0);
          s.push_back({p.x, p.y, z});
        }
  } else {
    for (auto& sh : b.sheets) {
      s.push_back(sh.start);
      s.push_back(sh.end);
    }
  }
  return s;
}
}  // namespace

ModelPtr wrap_bodies(const desc::Description& d, ModelPtr base, BuildLog& log) {
  if (d.bodies.empty()) return base;
  const auto& info = base->info();
  const auto& dom = info.domain;
  std::vector<Placed> placed;
  for (auto& bd : d.bodies) {
    Placed b;
    b.type = bd.type;
    b.rep.name = bd.name;
    b.rep.type = bd.type;
    b.rep.age = bd.age.si();
    if (bd.type == "salt_sheet") {
      b.sheet.footprint = polygon_of(bd.footprint);
      b.sheet.thickness = bd.thickness.si();
      b.sheet.taper = bd.taper.present ? bd.taper.si() : 0.0;
      if (bd.base_depth.present) {
        b.sheet.base_z = -bd.base_depth.si();
      } else {
        // base at a horizon of the host plus a labelled offset (positive up), sampled on a grid
        const int h = d.horizon_index(bd.base_horizon);
        BBox bb;
        for (auto& p : b.sheet.footprint) bb.add(p);
        const double sp = std::max(25.0, std::max(bb.hi.x - bb.lo.x, bb.hi.y - bb.lo.y) / 60.0);
        b.base_grid.x0 = bb.lo.x;
        b.base_grid.y0 = bb.lo.y;
        b.base_grid.dx = b.base_grid.dy = sp;
        b.base_grid.nx = static_cast<std::size_t>(std::ceil((bb.hi.x - bb.lo.x) / sp)) + 1;
        b.base_grid.ny = static_cast<std::size_t>(std::ceil((bb.hi.y - bb.lo.y) / sp)) + 1;
        b.base.assign(b.base_grid.nx * b.base_grid.ny, kNaN);
        const double off = bd.base_offset.present ? bd.base_offset.si() : 0.0;
        parallel_for(b.base.size(), 16, [&](std::size_t a0, std::size_t a1) {
          for (std::size_t n = a0; n < a1; ++n)
            b.base[n] = horizon_z(*base, h, b.base_grid.x0 + static_cast<double>(n % b.base_grid.nx) * sp, b.base_grid.y0 + static_cast<double>(n / b.base_grid.nx) * sp) + off;
        });
        // fill nodes where the horizon was not found by the mean of the found ones
        double mean = 0;
        std::size_t nf = 0;
        for (double z : b.base)
          if (std::isfinite(z)) { mean += z; ++nf; }
        if (nf == 0) throw Error(ErrorCode::Build, "body " + bd.name + ": horizon " + bd.base_horizon + " was not found under the sheet's footprint");
        mean /= static_cast<double>(nf);
        for (double& z : b.base)
          if (!std::isfinite(z)) z = mean;
        b.sheet.base_z = mean;
      }
      b.rep.volume = sheet_volume(b.sheet, std::max(10.0, std::sqrt(std::fabs(dom.sx * dom.sy)) / 200.0));
      b.rep.min_thickness = b.sheet.thickness;
    } else if (bd.type == "sill" || bd.type == "saucer_sill") {
      b.sill.centre = {bd.centre.si(0), bd.centre.si(1), -bd.depth.si()};
      b.sill.radius = bd.radius.si();
      b.sill.thickness = bd.thickness.si();
      b.sill.accommodation = bd.accommodation.present ? bd.accommodation.si() : 1.0;
      b.inner_radius = bd.radius.si();
      b.outer_radius = bd.type == "saucer_sill" ? bd.outer_radius.si() : bd.radius.si();
      b.climb = bd.type == "saucer_sill" ? bd.climb.si() : 0.0;
      b.decay = bd.uplift_decay.present ? bd.uplift_decay.si() : 0.0;
      SillBody whole = b.sill;
      whole.radius = b.outer_radius;
      b.rep.volume = sill_volume(whole);  // elliptical profile over the full (outer) radius
      b.rep.min_thickness = b.sill.thickness;
    } else if (bd.type == "dyke") {
      b.trace = polygon_of(bd.trace);
      b.width = bd.width.si();
      b.top = -bd.top_depth.si();
      b.bottom = -bd.base_depth.si();
      double len = 0;
      for (std::size_t i = 1; i < b.trace.size(); ++i) len += norm(b.trace[i] - b.trace[i - 1]);
      b.rep.volume = len * b.width * (b.top - b.bottom);
      b.rep.min_thickness = b.width;
    } else if (bd.type == "injectite") {
      // parent sand volume: the parent unit's thickness under the area times its net-to-gross
      int pu = -1, ch = d.horizon_index(bd.ceiling);
      for (std::size_t u = 0; u < info.units.size(); ++u)
        if (info.units[u].name == bd.parent) pu = static_cast<int>(u);
      if (pu < 0) throw Error(ErrorCode::Build, "body " + bd.name + ": parent unit " + bd.parent + " is not in the model");
      const auto area = polygon_of(bd.footprint);
      BBox bb;
      for (auto& p : area) bb.add(p);
      const int n = 20;
      double vol = 0, top = 0, ceil = 0;
      int cnt = 0;
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
          const double x = bb.lo.x + (bb.hi.x - bb.lo.x) * (i + 0.5) / n, y = bb.lo.y + (bb.hi.y - bb.lo.y) * (j + 0.5) / n;
          bool in = false;
          for (std::size_t a = 0, c = area.size() - 1; a < area.size(); c = a++)
            if ((area[a].y > y) != (area[c].y > y) && x < (area[c].x - area[a].x) * (y - area[a].y) / (area[c].y - area[a].y) + area[a].x) in = !in;
          if (!in) continue;
          const double zt = horizon_z(*base, pu - 1, x, y), zb = horizon_z(*base, pu, x, y), zc = horizon_z(*base, ch, x, y);
          if (!std::isfinite(zt) || !std::isfinite(zb) || !std::isfinite(zc)) continue;
          vol += (zt - zb) * (bb.hi.x - bb.lo.x) * (bb.hi.y - bb.lo.y) / (n * n);
          top += zt;
          ceil += zc;
          ++cnt;
        }
      if (cnt == 0) throw Error(ErrorCode::Build, "body " + bd.name + ": the parent unit or the ceiling is absent under the area");
      InjectiteSpec sp;
      sp.parent_top = top / cnt;
      sp.ceiling = ceil / cnt;
      sp.parent_sand_volume = vol * bd.net_to_gross.si();
      sp.max_fraction = bd.max_fraction.present ? bd.max_fraction.si() : 0.3;
      sp.sheets = static_cast<int>(std::llround(bd.sheets.si()));
      sp.length_mean = bd.length.si();
      sp.length_sd = bd.length_sd.present ? bd.length_sd.si() : 0.0;
      sp.thickness_mean = bd.thickness.si();
      sp.thickness_sd = bd.thickness_sd.present ? bd.thickness_sd.si() : 0.0;
      sp.climb_deg = bd.climb.si() * 180 / kPi;
      sp.area = bb;
      auto cx = build_injectites(sp, static_cast<std::uint64_t>(std::llround(bd.seed.si())));
      if (!cx.ok) throw Error(ErrorCode::Build, "body " + bd.name + ": " + cx.reason);
      b.sheets = cx.sheets;
      b.rep.volume = cx.injected_volume;
      b.rep.min_thickness = 1e300;
      for (auto& sh : b.sheets) b.rep.min_thickness = std::min(b.rep.min_thickness, sh.thickness);
      log.value("bodies." + bd.name + ".mass_balance", "parent sand " + fmt_shortest(sp.parent_sand_volume / 1e6) + " Mm3, injected " +
                                                          fmt_shortest(cx.injected_volume / 1e6) + " Mm3, parent after " + fmt_shortest(cx.parent_volume_after / 1e6) + " Mm3");
    }
    // a fault younger than the body that crosses it would have to displace it: not built
    const auto smp = samples_of(b);
    for (std::size_t f = 0; f < info.faults.size(); ++f) {
      const auto& fd = d.faults;
      double fage = kNaN;
      for (auto& x : fd)
        if (x.name == info.faults[f].name && x.age.present) fage = x.age.si();
      if (!(fage > b.rep.age) || smp.empty()) continue;
      FaultEval fe;
      base->evaluate_fault(static_cast<int>(f), smp, fe);
      bool pos = false, neg = false;
      for (std::size_t i = 0; i < smp.size(); ++i) {
        if (!(fe.extent[i] <= 1.0)) continue;
        (fe.distance[i] > 0 ? pos : neg) = true;
      }
      if (pos && neg)
        throw Error(ErrorCode::Build, "body " + bd.name + " (age " + fmt_shortest(b.rep.age) + ") is crossed by the younger fault " + info.faults[f].name +
                                          " (age " + fmt_shortest(fage) + "): displacing a body by a younger fault is not built; give the body an age after the fault's");
    }
    const double cell = d.grid_spacing.present ? d.grid_spacing.si() : 50.0;
    if (b.rep.min_thickness < 0.5 * cell) b.rep.note = "thinner than half the grid spacing (" + fmt_shortest(cell) + " m): coarse grids miss it";
    log.value("bodies." + bd.name, bd.type + ", age " + fmt_shortest(b.rep.age) + ", volume " + fmt_shortest(b.rep.volume / 1e6) + " Mm3" +
                                       (b.rep.note.empty() ? "" : "; " + b.rep.note));
    placed.push_back(std::move(b));
  }
  return std::make_shared<BodyModel>(base, std::move(placed));
}

}  // namespace stratonix
