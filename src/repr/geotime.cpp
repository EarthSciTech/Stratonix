// Stage 4A (STX-TN-14 section 5, decision D-047): geological-time coordinates. (u, v) is the
// restored horizontal position (the model's own restoration); t is interpolated in elevation in
// the restored-frame column (capability "column") between horizon ages, or - for a model without
// a column - between the horizon scalar fields at the present-day point (relative time of the
// implicit scalar field, Lajaunie et al. 1997). No volume parameterisation is solved.
#include "stratonix/repr/geotime.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"

namespace stratonix {

namespace {
const ColumnQuery* column_of(const Model& m) { return static_cast<const ColumnQuery*>(m.capability("column")); }
const StratigraphyQuery* strat_of(const Model& m) { return static_cast<const StratigraphyQuery*>(m.capability("stratigraphy")); }

/// t at elevation z in a column of horizon elevations e (top first) with ages a (young first).
/// Above the top and below the base: the nearest interval's rate. Zero-thickness intervals (absent
/// horizons) are jumped: the hiatus.
double interpolate_time(const std::vector<double>& e, const std::vector<double>& a, double z) {
  const std::size_t n = e.size();
  if (n == 1) return a[0];
  // first and last intervals of positive thickness, for extrapolation
  auto rate = [&](std::size_t i, std::size_t j) { return (a[j] - a[i]) / (e[i] - e[j]); };
  std::size_t top = 0;
  while (top + 1 < n && !(e[top] - e[top + 1] > 1e-9)) ++top;
  if (top + 1 >= n) return z >= e[0] ? a[0] : a[n - 1];
  if (z >= e[0]) return a[top] - (z - e[top]) * rate(top, top + 1);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    if (!(e[i] - e[i + 1] > 1e-9)) continue;
    if (z <= e[i] && z >= e[i + 1]) return a[i] + (e[i] - z) / (e[i] - e[i + 1]) * (a[i + 1] - a[i]);
  }
  std::size_t bot = n - 1;
  while (bot > 0 && !(e[bot - 1] - e[bot] > 1e-9)) --bot;
  if (bot == 0) return a[n - 1];
  return a[bot] + (e[bot] - z) * rate(bot - 1, bot);
}

/// Restored elevation of time t in a column (inverse of interpolate_time).
double interpolate_elevation(const std::vector<double>& e, const std::vector<double>& a, double t) {
  const std::size_t n = e.size();
  std::size_t top = 0;
  while (top + 1 < n && !(e[top] - e[top + 1] > 1e-9)) ++top;
  if (top + 1 >= n) return e[0];
  auto rate = [&](std::size_t i, std::size_t j) { return (a[j] - a[i]) / (e[i] - e[j]); };
  if (t <= a[top]) return e[top] - (t - a[top]) / rate(top, top + 1);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    if (!(e[i] - e[i + 1] > 1e-9)) continue;
    if (t >= a[i] && t <= a[i + 1]) return e[i] - (t - a[i]) / (a[i + 1] - a[i]) * (e[i] - e[i + 1]);
  }
  std::size_t bot = n - 1;
  while (bot > 0 && !(e[bot - 1] - e[bot] > 1e-9)) --bot;
  if (bot == 0) return e[n - 1];
  return e[bot] - (t - a[bot]) / rate(bot - 1, bot);
}
}  // namespace

GeologicalTime::GeologicalTime(const Model& m, std::vector<double> ages) : m_(m), ages_(std::move(ages)) {
  const std::size_t nh = m.info().horizons.size();
  if (nh == 0) throw Error(ErrorCode::InvalidArgument, "geological time needs at least one horizon");
  if (ages_.empty()) {
    ordinal_ = true;
    for (std::size_t h = 0; h < nh; ++h) ages_.push_back(static_cast<double>(h + 1));
  }
  if (ages_.size() != nh) throw Error(ErrorCode::InvalidArgument, "geological time: one age per horizon (" + std::to_string(nh) + "), got " + std::to_string(ages_.size()));
  for (std::size_t h = 0; h + 1 < nh; ++h)
    if (!(ages_[h + 1] > ages_[h]))
      throw Error(ErrorCode::InvalidArgument, "geological time: horizon ages must increase downwards (" + m.info().horizons[h].name + " " +
                                                  std::to_string(ages_[h]) + " Ma, " + m.info().horizons[h + 1].name + " " + std::to_string(ages_[h + 1]) + " Ma)");
}

GeoTimePoint GeologicalTime::forward(const Vec3& x) const {
  GeoTimePoint g;
  const std::size_t nh = ages_.size();
  Vec3 r;
  m_.restore(std::span<const Vec3>(&x, 1), std::span<Vec3>(&r, 1));
  if (!std::isfinite(r.x) || !std::isfinite(r.y) || !std::isfinite(r.z)) return g;
  EvalResult er;
  m_.evaluate(std::span<const Vec3>(&x, 1), kEvalAll, er);
  g.unit = er.unit.empty() ? -1 : er.unit[0];
  g.block = er.block.empty() ? 0 : er.block[0];
  g.u = r.x;
  g.v = r.y;
  if (const ColumnQuery* cq = column_of(m_)) {
    std::vector<double> e(nh);
    cq->column(r.x, r.y, e.data());
    for (double v : e)
      if (!std::isfinite(v)) return g;
    g.t = interpolate_time(e, ages_, r.z);
  } else {
    // relative time from the scalar fields: s_h is about the height above horizon h
    std::vector<double> e(nh);
    for (std::size_t h = 0; h < nh; ++h) {
      double s = er.s(h, 0);
      if (!std::isfinite(s)) return g;
      e[h] = x.z - s;  // local elevation of horizon h under (or over) the point
    }
    for (std::size_t h = 1; h < nh; ++h) e[h] = std::min(e[h], e[h - 1]);
    g.t = interpolate_time(e, ages_, x.z);
  }
  g.defined = std::isfinite(g.t);
  return g;
}

std::optional<Vec3> GeologicalTime::inverse(double u, double v, double t, double tolerance) const {
  const std::size_t nh = ages_.size();
  const auto& dom = m_.info().domain;
  if (const ColumnQuery* cq = column_of(m_)) {
    std::vector<double> e(nh);
    cq->column(u, v, e.data());
    for (double z : e)
      if (!std::isfinite(z)) return std::nullopt;
    const Vec3 R{u, v, interpolate_elevation(e, ages_, t)};
    // fixed point on the restoration: x <- x + R - restore(x)
    Vec3 x = R;
    for (int it = 0; it < 200; ++it) {
      Vec3 r;
      m_.restore(std::span<const Vec3>(&x, 1), std::span<Vec3>(&r, 1));
      if (!std::isfinite(r.x) || !std::isfinite(r.z)) return std::nullopt;
      Vec3 d = R - r;
      x = x + d;
      if (norm(d) <= tolerance) return x;
    }
    return std::nullopt;
  }
  // scalar-field models (no restoration of their own: restore is the identity or a translation)
  Vec3 x{u, v, 0.5 * (dom.z_top + dom.z_base)};
  for (int it = 0; it < 200; ++it) {
    Vec3 r;
    m_.restore(std::span<const Vec3>(&x, 1), std::span<Vec3>(&r, 1));
    x.x += u - r.x;
    x.y += v - r.y;
    auto g0 = forward(x), g1 = forward(x + Vec3{0, 0, 1});
    if (!g0.defined || !g1.defined) return std::nullopt;
    const double dtdz = g1.t - g0.t;
    if (!(std::fabs(dtdz) > 1e-12)) return std::nullopt;
    const double dz = std::clamp((t - g0.t) / dtdz, -0.25 * (dom.z_top - dom.z_base), 0.25 * (dom.z_top - dom.z_base));
    x.z += dz;
    if (std::fabs(dz) <= tolerance && std::fabs(u - r.x) <= tolerance && std::fabs(v - r.y) <= tolerance) return x;
  }
  return std::nullopt;
}

std::vector<Hiatus> GeologicalTime::hiatuses(double u, double v) const {
  std::vector<Hiatus> out;
  const std::size_t nh = ages_.size();
  const StratigraphyQuery* sq = strat_of(m_);
  const ColumnQuery* cq = column_of(m_);
  std::vector<bool> absent(nh, false);
  std::vector<double> e(nh, kNaN);
  if (cq) cq->column(u, v, e.data());
  for (std::size_t h = 0; h < nh; ++h) {
    if (sq && sq->horizon_absent(static_cast<int>(h), u, v)) absent[h] = true;
    // a horizon that coincides with the one below it in the restored column has no record here
    // (eroded or onlapped against it)
    if (h + 1 < nh && std::isfinite(e[h]) && std::isfinite(e[h + 1]) && !(e[h] - e[h + 1] > 1e-6) && h > 0) absent[h] = true;
  }
  for (std::size_t h = 0; h < nh; ++h) {
    if (!absent[h]) continue;
    Hiatus hi;
    hi.horizon_missing = static_cast<int>(h);
    for (std::size_t k = h; k-- > 0;)
      if (!absent[k]) { hi.t_young = ages_[k]; break; }
    for (std::size_t k = h + 1; k < nh; ++k)
      if (!absent[k]) { hi.t_old = ages_[k]; break; }
    out.push_back(hi);
  }
  return out;
}

WheelerDiagram wheeler_diagram(const GeologicalTime& gt, const std::string& name, const Polyline& line, std::size_t columns, std::size_t rows) {
  WheelerDiagram w;
  w.name = name;
  w.ages = gt.ages();
  w.ordinal = gt.ordinal();
  if (line.size() < 2 || columns < 2 || rows < 2) return w;
  const auto& dom = gt.model().info().domain;
  std::vector<double> cum{0};
  for (std::size_t i = 1; i < line.size(); ++i) cum.push_back(cum.back() + std::hypot(line[i].x - line[i - 1].x, line[i].y - line[i - 1].y));
  const double L = cum.back();
  // the time axis is clipped half a mean interval beyond the first and last horizons
  const double half = w.ages.size() > 1 ? 0.5 * (w.ages.back() - w.ages.front()) / static_cast<double>(w.ages.size() - 1) : 0.5;
  const double tmin = w.ages.front() - half, tmax = w.ages.back() + half;
  w.columns.resize(columns);
  std::vector<std::size_t> undefined(columns, 0);
  parallel_for(columns, 4, [&](std::size_t b, std::size_t e) {
    for (std::size_t c = b; c < e; ++c) {
      const double d = L * static_cast<double>(c) / static_cast<double>(columns - 1);
      std::size_t seg = 1;
      while (seg + 1 < line.size() && cum[seg] < d) ++seg;
      const double f = cum[seg] > cum[seg - 1] ? (d - cum[seg - 1]) / (cum[seg] - cum[seg - 1]) : 0;
      const double x = line[seg - 1].x + f * (line[seg].x - line[seg - 1].x), y = line[seg - 1].y + f * (line[seg].y - line[seg - 1].y);
      auto& col = w.columns[c];
      col.station = d;
      auto record = [&](int unit, double ta, double tb) {
        double a = std::clamp(std::min(ta, tb), tmin, tmax), bb = std::clamp(std::max(ta, tb), tmin, tmax);
        if (!(bb > a)) return;
        // merge with the column's last interval of the same unit when they touch
        if (!col.intervals.empty() && col.intervals.back().unit == unit && std::fabs(col.intervals.back().t_old - a) < 1e-6)
          col.intervals.back().t_old = std::max(col.intervals.back().t_old, bb);
        else
          col.intervals.push_back({a, bb, unit});
      };
      GeoTimePoint prev;
      double zprev = 0;
      for (std::size_t r = 0; r < rows; ++r) {
        const double z = dom.z_top - (dom.z_top - dom.z_base) * (static_cast<double>(r) + 0.5) / static_cast<double>(rows);
        auto g = gt.forward({x, y, z});
        if (!g.defined) { ++undefined[c]; prev = g; zprev = z; continue; }
        if (prev.defined && prev.unit == g.unit && prev.block == g.block) {
          record(g.unit, prev.t, g.t);
        } else if (prev.defined) {
          // a boundary between the samples: bisect it, and extend both sides to it. A time jump
          // across a horizon is a hiatus (blank); across a fault it is a fault cut-out.
          double za = zprev, zb = z;
          GeoTimePoint ga = prev, gb = g;
          for (int it = 0; it < 20; ++it) {
            const double zm = 0.5 * (za + zb);
            auto gm = gt.forward({x, y, zm});
            if (gm.defined && gm.unit == prev.unit && gm.block == prev.block) { za = zm; ga = gm; } else { zb = zm; gb = gm; }
          }
          record(prev.unit, prev.t, ga.t);
          if (gb.defined && ga.block != gb.block) record(kWheelerFault, ga.t, gb.t);
          if (gb.defined && gb.unit == g.unit && gb.block == g.block) record(g.unit, gb.t, g.t);
        }
        prev = g;
        zprev = z;
      }
    }
  });
  for (auto u : undefined) w.undefined += u;
  return w;
}

}  // namespace stratonix
