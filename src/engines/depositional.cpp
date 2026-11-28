// Stage 4B (STX-TN-15 section 3): stratal terminations, clinoform sets, channel belts and lobes.
#include "stratonix/engines/depositional.hpp"

#include <algorithm>
#include <cmath>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"

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
};
/// Plan distance from (x, y) to a polyline.
double plan_distance(const std::vector<Vec3>& line, double x, double y) {
  double best = 1e300;
  for (std::size_t i = 1; i < line.size(); ++i) {
    const double ax = line[i - 1].x, ay = line[i - 1].y, dx = line[i].x - ax, dy = line[i].y - ay;
    const double L2 = dx * dx + dy * dy;
    const double t = L2 > 0 ? std::clamp(((x - ax) * dx + (y - ay) * dy) / L2, 0.0, 1.0) : 0.0;
    best = std::min(best, std::hypot(x - ax - t * dx, y - ay - t * dy));
  }
  return best;
}
}  // namespace

std::vector<Vec3> onlap_line(const std::function<double(double, double)>& horizon, const std::function<double(double, double)>& surface,
                             const GridSpec& g) {
  // where h - s changes sign between neighbouring nodes, the crossing by linear interpolation
  std::vector<Vec3> out;
  auto d = [&](double x, double y) { return horizon(x, y) - surface(x, y); };
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i) {
      const double x = g.x0 + g.dx * static_cast<double>(i), y = g.y0 + g.dy * static_cast<double>(j), d0 = d(x, y);
      auto edge = [&](double x1, double y1) {
        const double d1 = d(x1, y1);
        if ((d0 > 0) == (d1 > 0) || d0 == d1) return;
        const double t = d0 / (d0 - d1), xc = x + t * (x1 - x), yc = y + t * (y1 - y);
        out.push_back({xc, yc, horizon(xc, yc)});
      };
      if (i + 1 < g.nx) edge(x + g.dx, y);
      if (j + 1 < g.ny) edge(x, y + g.dy);
    }
  return out;
}

std::vector<Clinoform> clinoform_set(const ClinoformSpec& s) {
  if (s.count < 2) throw Error(ErrorCode::InvalidArgument, "clinoform set: at least two clinoforms");
  if (!(s.foreset_dip_deg >= 0.5 && s.foreset_dip_deg <= 35)) throw Error(ErrorCode::InvalidArgument, "clinoform set: foreset dip from 0.5 to 35 degrees");
  if (!(s.trajectory_deg >= -5 && s.trajectory_deg <= 15)) throw Error(ErrorCode::InvalidArgument, "clinoform set: trajectory angle from -5 to 15 degrees");
  if (!(s.spacing > 0) || !(s.shelf_z > s.base_z)) throw Error(ErrorCode::InvalidArgument, "clinoform set: positive spacing, shelf above base");
  if (s.shape != "oblique" && s.shape != "sigmoidal" && s.shape != "tangential")
    throw Error(ErrorCode::InvalidArgument, "clinoform set: shape oblique, sigmoidal or tangential");
  const double tdip = std::tan(s.foreset_dip_deg * kPi / 180), ttr = std::tan(s.trajectory_deg * kPi / 180);
  std::vector<Clinoform> out;
  for (int i = 0; i < s.count; ++i) {
    Clinoform c;
    c.shelf_break_x = s.x_start + i * s.spacing;
    c.shelf_break_z = s.shelf_z + i * s.spacing * ttr;
    const double H = c.shelf_break_z - s.base_z;
    // the profile z(x) (Steel and Olsen 2002): topset to the break, foreset, bottomset
    auto z = [&](double x) {
      const double u = x - c.shelf_break_x;
      if (s.shape == "oblique") return std::max(s.base_z, c.shelf_break_z - std::max(0.0, u) * tdip);
      if (s.shape == "tangential") return u <= 0 ? c.shelf_break_z : s.base_z + H * std::exp(-u * tdip / H);
      // sigmoidal: tanh with its steepest slope tan(dip), the break at the upper bend (u = -w)
      const double w = H / (2 * tdip);
      return s.base_z + 0.5 * H * (1 - std::tanh((u - w) / w));
    };
    const double x0 = c.shelf_break_x - 2 * s.spacing, x1 = c.shelf_break_x + 6 * H / tdip;
    double maxdip = 0, zp = z(x0);
    for (int k = 0; k <= 2000; ++k) {
      const double x = x0 + (x1 - x0) * k / 2000.0, zz = z(x);
      c.profile.push_back({x, 0, zz});
      if (k > 0) maxdip = std::max(maxdip, std::atan((zp - zz) / ((x1 - x0) / 2000.0)));
      zp = zz;
    }
    c.foreset_dip_deg = maxdip * 180 / kPi;
    out.push_back(std::move(c));
  }
  return out;
}

bool Geobody::contains(const Vec3& p) const {
  if (p.x < x0 || p.x > x1 || p.y < y0 || p.y > y1 || p.z < z_base || p.z >= z_base + thickness) return false;
  if (kind == "lobe") {
    const double ca = std::cos(orientation), sa = std::sin(orientation);
    const double dx = p.x - centre.x, dy = p.y - centre.y;
    const double along = dx * sa + dy * ca, across = dx * ca - dy * sa;  // azimuth convention
    const double b = radius / aspect;
    const double r2 = (along * along) / (radius * radius) + (across * across) / (b * b);
    if (r2 >= 1) return false;
    return p.z < z_base + thickness * (1 - r2);
  }
  const double d = plan_distance(centreline, p.x, p.y);
  if (!(d < 0.5 * width)) return false;
  if (shape == "parabolic") {
    const double q = 2 * d / width;
    return p.z >= z_base + thickness * q * q;  // the channel's top is flat at z_base + thickness
  }
  return true;
}

double Geobody::volume() const {
  if (kind == "lobe") return 0.5 * kPi * radius * (radius / aspect) * thickness;  // paraboloid lens
  double L = 0;
  for (std::size_t i = 1; i < centreline.size(); ++i) L += std::hypot(centreline[i].x - centreline[i - 1].x, centreline[i].y - centreline[i - 1].y);
  return L * width * thickness * (shape == "parabolic" ? 2.0 / 3.0 : 1.0);
}

std::vector<Geobody> place_channels(const ChannelSpec& s, const UnitBox& u, std::uint64_t seed) {
  if (s.count < 1 || !(s.width > 0) || !(s.thickness > 0)) throw Error(ErrorCode::InvalidArgument, "channels: a count, positive width and thickness");
  const double wt = s.width / s.thickness;
  if (wt < 5 || wt > 500) throw Error(ErrorCode::InvalidArgument, "channels: width-to-thickness ratio from 5 to 500");
  if (s.thickness > u.z_top - u.z_base) throw Error(ErrorCode::InvalidArgument, "channels: thicker than their unit");
  if (s.shape != "rectangular" && s.shape != "parabolic") throw Error(ErrorCode::InvalidArgument, "channels: shape rectangular or parabolic");
  if (s.amplitude > 0 && !(s.wavelength > 0)) throw Error(ErrorCode::InvalidArgument, "channels: a sinuous channel needs a wavelength");
  if (s.amplitude > 0) {
    // sinuosity of a sine centreline: arc length per wavelength / wavelength (at most 3)
    const double k = 2 * kPi * s.amplitude / s.wavelength;
    double L = 0;
    for (int i = 0; i < 400; ++i) L += std::sqrt(1 + k * k * std::pow(std::cos(2 * kPi * (i + 0.5) / 400), 2)) / 400;
    if (L > 3) throw Error(ErrorCode::InvalidArgument, "channels: sinuosity above 3");
  }
  Rng g(seed);
  const double az = s.orientation_deg * kPi / 180;
  const Vec3 along{std::sin(az), std::cos(az), 0}, across{std::cos(az), -std::sin(az), 0};
  const Vec3 c0{0.5 * (u.x0 + u.x1), 0.5 * (u.y0 + u.y1), 0};
  // extent of the box along and across the channels' direction
  double amin = 1e300, amax = -1e300, cmin = 1e300, cmax = -1e300;
  for (double x : {u.x0, u.x1})
    for (double y : {u.y0, u.y1}) {
      const Vec3 d{x - c0.x, y - c0.y, 0};
      amin = std::min(amin, dot(d, along)); amax = std::max(amax, dot(d, along));
      cmin = std::min(cmin, dot(d, across)); cmax = std::max(cmax, dot(d, across));
    }
  std::vector<Geobody> out;
  for (int k = 0; k < s.count; ++k) {
    Geobody b;
    b.kind = "channel";
    b.width = s.width;
    b.thickness = s.thickness;
    b.shape = s.shape;
    double off, phase = 0;
    if (s.random_placement) {
      off = cmin + (cmax - cmin) * g.uniform();
      b.z_base = u.z_base + (u.z_top - u.z_base - s.thickness) * g.uniform();
      phase = 2 * kPi * g.uniform();
    } else {
      off = cmin + (cmax - cmin) * (k + 0.5) / s.count;
      b.z_base = u.z_base + 0.5 * (u.z_top - u.z_base - s.thickness);
    }
    const int n = s.amplitude > 0 ? std::max(2, static_cast<int>(std::ceil((amax - amin) / s.wavelength * 48))) : 1;
    for (int i = 0; i <= n; ++i) {
      const double a = amin + (amax - amin) * i / n;
      const double c = off + (s.amplitude > 0 ? s.amplitude * std::sin(2 * kPi * a / s.wavelength + phase) : 0.0);
      b.centreline.push_back(c0 + along * a + across * c);
    }
    b.x0 = b.y0 = 1e300;
    b.x1 = b.y1 = -1e300;
    for (auto& q : b.centreline) {
      b.x0 = std::min(b.x0, q.x - 0.5 * s.width); b.x1 = std::max(b.x1, q.x + 0.5 * s.width);
      b.y0 = std::min(b.y0, q.y - 0.5 * s.width); b.y1 = std::max(b.y1, q.y + 0.5 * s.width);
    }
    out.push_back(std::move(b));
  }
  return out;
}

std::vector<Geobody> place_lobes(const LobeSpec& s, const UnitBox& u, std::uint64_t seed) {
  if (s.count < 1 || !(s.radius > 0) || !(s.thickness > 0) || !(s.aspect >= 1))
    throw Error(ErrorCode::InvalidArgument, "lobes: a count, positive radius and thickness, aspect of at least 1");
  const double wt = 2 * s.radius / s.aspect / s.thickness;
  if (wt < 5 || wt > 500) throw Error(ErrorCode::InvalidArgument, "lobes: width-to-thickness ratio from 5 to 500");
  if (s.thickness > u.z_top - u.z_base) throw Error(ErrorCode::InvalidArgument, "lobes: thicker than their unit");
  Rng g(seed);
  std::vector<Geobody> out;
  for (int k = 0; k < s.count; ++k) {
    Geobody b;
    b.kind = "lobe";
    b.radius = s.radius;
    b.aspect = s.aspect;
    b.thickness = s.thickness;
    b.orientation = s.orientation_deg * kPi / 180;
    b.centre = {u.x0 + (u.x1 - u.x0) * g.uniform(), u.y0 + (u.y1 - u.y0) * g.uniform(), 0};
    b.z_base = u.z_base + (u.z_top - u.z_base - s.thickness) * g.uniform();
    b.x0 = b.centre.x - s.radius; b.x1 = b.centre.x + s.radius;
    b.y0 = b.centre.y - s.radius; b.y1 = b.centre.y + s.radius;
    out.push_back(std::move(b));
  }
  return out;
}

double net_to_gross(const std::vector<Geobody>& b, const UnitBox& u, double spacing) {
  if (!(spacing > 0)) throw Error(ErrorCode::InvalidArgument, "net-to-gross: positive lattice spacing");
  const auto nx = static_cast<std::size_t>(std::max(1.0, std::round((u.x1 - u.x0) / spacing)));
  const auto ny = static_cast<std::size_t>(std::max(1.0, std::round((u.y1 - u.y0) / spacing)));
  const auto nz = static_cast<std::size_t>(std::max(1.0, std::round((u.z_top - u.z_base) / spacing)));
  const double dx = (u.x1 - u.x0) / nx, dy = (u.y1 - u.y0) / ny, dz = (u.z_top - u.z_base) / nz;
  std::size_t in = 0;
  for (std::size_t j = 0; j < ny; ++j)
    for (std::size_t i = 0; i < nx; ++i) {
      const double x = u.x0 + (i + 0.5) * dx, y = u.y0 + (j + 0.5) * dy;
      std::vector<const Geobody*> here;
      for (auto& g : b)
        if (x >= g.x0 && x <= g.x1 && y >= g.y0 && y <= g.y1) here.push_back(&g);
      if (here.empty()) continue;
      for (std::size_t k = 0; k < nz; ++k) {
        const Vec3 p{x, y, u.z_base + (k + 0.5) * dz};
        for (auto* g : here)
          if (g->contains(p)) { ++in; break; }
      }
    }
  return static_cast<double>(in) / static_cast<double>(nx * ny * nz);
}

// ---- geobodies in models ------------------------------------------------------------------
namespace {
struct GeoSet {
  GeobodyReport rep;
  int host = -1;                 ///< host unit index in the model
  std::vector<Geobody> bodies;   ///< in box coordinates: z from 0 (unit base) to box_thickness (unit top)
};

class GeobodyModel final : public Model {
 public:
  GeobodyModel(ModelPtr inner, std::vector<GeoSet> sets) : inner_(std::move(inner)), sets_(std::move(sets)), info_(inner_->info()) {
    col_ = static_cast<const ColumnQuery*>(inner_->capability("column"));
    nh_ = info_.horizons.size();
    for (auto& s : sets_) reports_.push_back(s.rep);
  }
  const ModelInfo& info() const override { return info_; }
  void set_units(std::vector<UnitInfo> extra) {
    for (std::size_t k = 0; k < extra.size(); ++k) {
      info_.units.push_back(extra[k]);
      sets_[k].rep.unit = static_cast<int>(info_.units.size()) - 1;
      reports_[k].unit = sets_[k].rep.unit;
    }
  }
  void evaluate(std::span<const Vec3> pts, unsigned fields, EvalResult& out) const override {
    const unsigned f = fields | kEvalUnit;
    inner_->evaluate(pts, f, out);
    if (!(fields & kEvalUnit)) return;
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < pts.size(); ++i)
      for (auto& s : sets_)
        if (out.unit[i] == s.host) { idx.push_back(i); break; }
    if (idx.empty()) return;
    std::vector<Vec3> q(idx.size()), r(idx.size());
    for (std::size_t k = 0; k < idx.size(); ++k) q[k] = pts[idx[k]];
    inner_->restore(q, r);
    std::vector<double> e(nh_);
    for (std::size_t k = 0; k < idx.size(); ++k) {
      const std::size_t i = idx[k];
      col_->column(r[k].x, r[k].y, e.data());
      // later sets first (erosional): they take a point before earlier ones
      for (std::size_t s = sets_.size(); s-- > 0;) {
        const GeoSet& S = sets_[s];
        if (out.unit[i] != S.host) continue;
        const double top = e[static_cast<std::size_t>(S.host) - 1], base = e[static_cast<std::size_t>(S.host)];
        if (!(top > base)) continue;
        const double fz = std::clamp((r[k].z - base) / (top - base), 0.0, 1.0);
        const Vec3 b{r[k].x, r[k].y, fz * S.rep.box_thickness};
        bool in = false;
        for (auto& g : S.bodies)
          if (g.contains(b)) { in = true; break; }
        if (in) { out.unit[i] = S.rep.unit; break; }
      }
    }
  }
  void scalar_in_block(std::span<const Vec3> pts, BlockCode block, int horizon, std::span<double> out) const override {
    inner_->scalar_in_block(pts, block, horizon, out);
  }
  void evaluate_fault(int fault, std::span<const Vec3> pts, FaultEval& out) const override { inner_->evaluate_fault(fault, pts, out); }
  void evaluate_fault_geometry(int fault, std::span<const Vec3> pts, FaultEval& out) const override { inner_->evaluate_fault_geometry(fault, pts, out); }
  void restore(std::span<const Vec3> pts, std::span<Vec3> out) const override { inner_->restore(pts, out); }
  const void* capability(std::string_view name) const override {
    if (name == "geobodies") return &reports_;
    return inner_->capability(name);
  }

 private:
  ModelPtr inner_;
  std::vector<GeoSet> sets_;
  std::vector<GeobodyReport> reports_;
  ModelInfo info_;
  const ColumnQuery* col_ = nullptr;
  std::size_t nh_ = 0;
};
}  // namespace

ModelPtr wrap_geobodies(const desc::Description& d, ModelPtr base, BuildLog& log) {
  if (d.geobodies.empty()) return base;
  const auto* col = static_cast<const ColumnQuery*>(base->capability("column"));
  if (!col) throw Error(ErrorCode::Build, "geobodies need a model with a restored stratigraphic column (rule F39)", "geobodies");
  const ModelInfo& info = base->info();
  const Domain& dom = info.domain;
  const std::size_t nh = info.horizons.size();
  std::vector<GeoSet> sets;
  std::vector<UnitInfo> units;
  std::vector<double> e(nh);
  col->column(dom.x0 + 0.5 * dom.sx, dom.y0 + 0.5 * dom.sy, e.data());
  for (auto& g : d.geobodies) {
    GeoSet S;
    S.rep.name = g.name;
    S.rep.type = g.type;
    S.rep.host = g.unit;
    for (std::size_t k = 0; k < info.units.size(); ++k)
      if (info.units[k].name == g.unit) S.host = static_cast<int>(k);
    if (S.host < 1 || static_cast<std::size_t>(S.host) >= nh)
      throw Error(ErrorCode::Build, "geobodies " + g.name + ": unit " + g.unit + " is not between two horizons of the model (rule F39)", "geobodies");
    const double T = e[static_cast<std::size_t>(S.host) - 1] - e[static_cast<std::size_t>(S.host)];
    if (!(T > 0)) throw Error(ErrorCode::Build, "geobodies " + g.name + ": unit " + g.unit + " has no thickness at the domain centre", "geobodies");
    S.rep.box_thickness = T;
    const UnitBox box{dom.x0, dom.y0, dom.x0 + dom.sx, dom.y0 + dom.sy, 0.0, T};
    const auto seed = static_cast<std::uint64_t>(g.seed.si());
    try {
      if (g.type == "channels") {
        ChannelSpec cs;
        cs.count = static_cast<int>(g.count.si());
        cs.width = g.width.si();
        cs.thickness = g.thickness.si();
        cs.wavelength = g.wavelength.present ? g.wavelength.si() : 0.0;
        cs.amplitude = g.amplitude.present ? g.amplitude.si() : 0.0;
        cs.orientation_deg = rad2deg(g.orientation.si());
        cs.shape = g.shape;
        S.bodies = place_channels(cs, box, seed);
      } else {
        LobeSpec ls;
        ls.count = static_cast<int>(g.count.si());
        ls.radius = g.radius.si();
        ls.aspect = g.aspect.si();
        ls.thickness = g.thickness.si();
        ls.orientation_deg = rad2deg(g.orientation.si());
        S.bodies = place_lobes(ls, box, seed);
      }
    } catch (const Error& er) {
      throw Error(ErrorCode::Build, "geobodies " + g.name + " (rule F39): " + er.what(), "geobodies");
    }
    S.rep.count = static_cast<int>(S.bodies.size());
    double vol = 0;
    for (auto& b : S.bodies) vol += b.volume();
    S.rep.ntg_closed = vol / (dom.sx * dom.sy * T);
    // lattice estimate in the box: 200 x 200 columns, 40 levels
    {
      const std::size_t n = 200, nz = 40;
      std::size_t in = 0;
      for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i)
          for (std::size_t k = 0; k < nz; ++k) {
            const Vec3 p{box.x0 + (i + 0.5) * (box.x1 - box.x0) / n, box.y0 + (j + 0.5) * (box.y1 - box.y0) / n, (k + 0.5) * T / nz};
            for (auto& b : S.bodies)
              if (b.contains(p)) { ++in; break; }
          }
      S.rep.ntg = static_cast<double>(in) / static_cast<double>(n * n * nz);
    }
    log.value("geobodies." + g.name, g.type + " in " + g.unit + ": " + std::to_string(S.rep.count) + ", net-to-gross " + fmt_fixed(S.rep.ntg, 4) +
                                         " (closed form " + fmt_fixed(S.rep.ntg_closed, 4) + "), unit thickness at the domain centre " + fmt_fixed(T, 1) + " m");
    UnitInfo u;
    u.name = g.name;
    if (!g.color.empty()) u.color = g.color;
    units.push_back(u);
    sets.push_back(std::move(S));
  }
  auto m = std::make_shared<GeobodyModel>(std::move(base), std::move(sets));
  m->set_units(std::move(units));
  return m;
}

}  // namespace stratonix
