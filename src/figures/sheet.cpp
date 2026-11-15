// Section sheet: a presentation page for one section with a main panel, a detail panel, a
// legend by role, computed key numbers and notes. Unit fills are polygons per fault block
// between horizon lines found by root finding in that block, so offsets at faults are sharp.
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>

#include "stratonix/core/roots.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/figures/figures.hpp"

namespace stratonix::fig {
namespace {

constexpr double kInfD = std::numeric_limits<double>::infinity();
const char* kFont = "DejaVu Sans";
const Color kAxis = rgb(68, 68, 68), kTick = rgb(85, 85, 85), kTitle = rgb(34, 34, 34), kLine = rgb(21, 21, 21),
            kBlack = rgb(0, 0, 0), kCompletion = rgb(232, 17, 45), kNote = rgb(176, 0, 31), kGrey = rgb(85, 85, 85);

// ---- section geometry -----------------------------------------------------------------------
struct Interval {
  double top = 0, bot = 0;      ///< elevations, top > bot
  BlockCode b = 0;
  std::vector<double> h;        ///< per horizon: root (elevation), +inf above, -inf below the search range
};
struct Geometry {
  Polyline line;
  double length = 0, z_top = 0, z_bot = 0;
  std::vector<double> st;
  std::vector<Vec3> xy;
  std::vector<std::vector<double>> fz;    ///< [col][fault] elevation of the fault surface or NaN
  std::vector<std::vector<double>> fext;  ///< [col][fault] extent coordinate there (<= 1 inside the tips)
  std::vector<std::vector<Interval>> iv;  ///< [col] intervals top to bottom
};

Vec3 point_at(const Polyline& l, double s) {
  double acc = 0;
  for (std::size_t i = 1; i < l.size(); ++i) {
    double d = std::hypot(l[i].x - l[i - 1].x, l[i].y - l[i - 1].y);
    if (s <= acc + d || i + 1 == l.size()) {
      double t = d > 0 ? std::clamp((s - acc) / d, 0.0, 1.0) : 0.0;
      return {l[i - 1].x + t * (l[i].x - l[i - 1].x), l[i - 1].y + t * (l[i].y - l[i - 1].y), 0};
    }
    acc += d;
  }
  return l.front();
}
double polyline_length(const Polyline& l) {
  double acc = 0;
  for (std::size_t i = 1; i < l.size(); ++i) acc += std::hypot(l[i].x - l[i - 1].x, l[i].y - l[i - 1].y);
  return acc;
}

Geometry build_geometry(const Model& m, const Polyline& line, double s0, double s1, double ds, double z_top, double z_bot) {
  Geometry g;
  g.line = line;
  g.length = polyline_length(line);
  g.z_top = z_top;
  g.z_bot = z_bot;
  const std::size_t nf = m.info().faults.size(), nh = m.info().horizons.size();
  std::size_t n = static_cast<std::size_t>(std::ceil((s1 - s0) / ds - 1e-9)) + 1;
  for (std::size_t i = 0; i < n; ++i) {
    double s = std::min(s1, s0 + ds * static_cast<double>(i));
    g.st.push_back(s);
    g.xy.push_back(point_at(line, s));
  }
  g.fz.assign(n, std::vector<double>(nf, kNaN));
  g.fext.assign(n, std::vector<double>(nf, kNaN));
  g.iv.assign(n, {});
  const double margin = 150.0;
  parallel_for(n, 16, [&](std::size_t a, std::size_t b) {
    for (std::size_t i = a; i < b; ++i) {
      const double x = g.xy[i].x, y = g.xy[i].y;
      std::vector<double> cuts;
      for (std::size_t f = 0; f < nf; ++f) {
        const int NP = 24;
        std::vector<Vec3> p(NP + 1);
        for (int k = 0; k <= NP; ++k) p[static_cast<std::size_t>(k)] = {x, y, z_top + (z_bot - z_top) * k / NP};
        FaultEval fe;
        m.evaluate_fault(static_cast<int>(f), p, fe);
        auto dist = [&](double z) {
          FaultEval e;
          Vec3 q{x, y, z};
          m.evaluate_fault(static_cast<int>(f), std::span<const Vec3>(&q, 1), e);
          return e.distance[0];
        };
        for (int k = 0; k < NP; ++k) {
          double da = fe.distance[static_cast<std::size_t>(k)], db = fe.distance[static_cast<std::size_t>(k + 1)];
          if ((da > 0) == (db > 0)) continue;
          double z = illinois(dist, p[static_cast<std::size_t>(k)].z, p[static_cast<std::size_t>(k + 1)].z, da, db, 1e-6);
          Vec3 q{x, y, z};
          FaultEval e;
          m.evaluate_fault(static_cast<int>(f), std::span<const Vec3>(&q, 1), e);
          g.fz[i][f] = z;
          g.fext[i][f] = e.extent[0];
          cuts.push_back(z);
          break;
        }
      }
      std::sort(cuts.begin(), cuts.end(), std::greater<double>());
      std::vector<double> edges{z_top};
      for (double c : cuts)
        if (c < edges.back() - 1e-6 && c > z_bot + 1e-6) edges.push_back(c);
      edges.push_back(z_bot);
      for (std::size_t k = 0; k + 1 < edges.size(); ++k) {
        Interval iv;
        iv.top = edges[k];
        iv.bot = edges[k + 1];
        Vec3 mid{x, y, 0.5 * (iv.top + iv.bot)};
        EvalResult r;
        m.evaluate(std::span<const Vec3>(&mid, 1), kEvalBlock, r);
        iv.b = r.block[0];
        iv.h.assign(nh, kNaN);
        double hi = iv.top + margin, lo = iv.bot - margin;
        for (std::size_t h = 0; h < nh; ++h) {
          auto sc = [&](double z) {
            Vec3 q{x, y, z};
            double out = 0;
            m.scalar_in_block(std::span<const Vec3>(&q, 1), iv.b, static_cast<int>(h), std::span<double>(&out, 1));
            return out;
          };
          double fh = sc(hi), fl = sc(lo);
          if (fh < 0) iv.h[h] = kInfD;          // the horizon lies above the range
          else if (fl > 0) iv.h[h] = -kInfD;    // below it
          else iv.h[h] = illinois(sc, lo, hi, fl, fh, 1e-5);
        }
        g.iv[i].push_back(std::move(iv));
      }
    }
  });
  return g;
}

const Interval* block_interval(const Geometry& g, std::size_t i, BlockCode b) {
  for (auto& iv : g.iv[i])
    if (iv.b == b) return &iv;
  return nullptr;
}
double clampv(double v, double lo, double hi) { return std::min(hi, std::max(lo, v)); }

/// Cut-off of horizon h on fault f, on the side whose block bit equals `hw` (1 hanging wall).
bool cutoff(const Geometry& g, int h, int f, bool hw, double& s_out, double& z_out) {
  const std::size_t n = g.st.size();
  auto side_root = [&](std::size_t i, double& v) {
    double zf = g.fz[i][static_cast<std::size_t>(f)];
    if (std::isnan(zf)) return false;
    for (auto& iv : g.iv[i]) {
      bool bit = ((iv.b >> f) & 1u) != 0;
      if (bit != hw) continue;
      if (std::fabs(iv.top - zf) < 1e-6 || std::fabs(iv.bot - zf) < 1e-6) {
        double r = iv.h[static_cast<std::size_t>(h)];
        if (!std::isfinite(r)) return false;
        v = r - zf;
        return true;
      }
    }
    return false;
  };
  for (std::size_t i = 0; i + 1 < n; ++i) {
    double a, b;
    if (!side_root(i, a) || !side_root(i + 1, b)) continue;
    if ((a > 0) == (b > 0) && a != 0) continue;
    double t = a == b ? 0 : a / (a - b);
    s_out = g.st[i] + t * (g.st[i + 1] - g.st[i]);
    z_out = g.fz[i][static_cast<std::size_t>(f)] + t * (g.fz[i + 1][static_cast<std::size_t>(f)] - g.fz[i][static_cast<std::size_t>(f)]);
    return true;
  }
  return false;
}

/// Fault trace inside the tip line: points (s, z), shallow end first.
std::vector<std::pair<double, double>> fault_trace(const Geometry& g, int f) {
  std::vector<std::pair<double, double>> pts;
  const std::size_t fi = static_cast<std::size_t>(f);
  for (std::size_t i = 0; i < g.st.size(); ++i) {
    double z = g.fz[i][fi], e = g.fext[i][fi];
    if (std::isnan(z)) continue;
    bool in = e <= 1.0;
    // ends: interpolate where the extent coordinate passes 1 between neighbouring columns
    if (i > 0 && !std::isnan(g.fz[i - 1][fi])) {
      double e0 = g.fext[i - 1][fi];
      if ((e0 <= 1.0) != in) {
        double t = (1.0 - e0) / (e - e0);
        pts.push_back({g.st[i - 1] + t * (g.st[i] - g.st[i - 1]), g.fz[i - 1][fi] + t * (z - g.fz[i - 1][fi])});
      }
    }
    if (in) pts.push_back({g.st[i], z});
  }
  std::sort(pts.begin(), pts.end(), [](auto& a, auto& b) { return a.second > b.second; });
  return pts;
}

std::string fmt_dec(double v, int d) {
  std::string s = fmt_fixed(v, d < 0 ? 1 : d);
  if (s == "-0" || s == "-0.0" || s == "-0.00") s = s.substr(1);
  return s;
}
std::string minus(const std::string& s) { return !s.empty() && s[0] == '-' ? "−" + s.substr(1) : s; }
std::string num_word(std::size_t n) {
  static const char* w[] = {"no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve"};
  return n < 13 ? w[n] : std::to_string(n);
}
std::string fmt_g(double v) {
  std::ostringstream o;
  o << v;
  return o.str();
}
/// Replace {key} and {key:N} placeholders.
std::string fill(std::string t, const std::map<std::string, double>& vals, const std::map<std::string, std::string>& words = {}) {
  for (auto& [k, v] : words) {
    std::string tag = "{" + k + "}";
    for (std::size_t p; (p = t.find(tag)) != std::string::npos;) t.replace(p, tag.size(), v);
  }
  for (auto& [k, v] : vals) {
    for (std::size_t p; (p = t.find("{" + k)) != std::string::npos;) {
      std::size_t e = t.find('}', p);
      if (e == std::string::npos) break;
      std::string spec = t.substr(p + 1 + k.size(), e - p - 1 - k.size());
      std::string out = fmt_g(v);
      if (!spec.empty() && spec[0] == ':') out = fmt_dec(v, std::atoi(spec.c_str() + 1));
      t.replace(p, e - p + 1, out);
    }
  }
  return t;
}
std::vector<std::string> split_lines(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (s[i] == '\\' && i + 1 < s.size() && s[i + 1] == 'n') { out.push_back(cur); cur.clear(); ++i; continue; }
    if (s[i] == '\n') { out.push_back(cur); cur.clear(); continue; }
    cur += s[i];
  }
  out.push_back(cur);
  return out;
}

// ---- panel ------------------------------------------------------------------------------------
struct Panel {
  double x0, y0, x1, y1;   // pt
  double s0, s1, z0, z1;   // data window: distance, elevation (z0 top)
  Pt p(double s, double z) const { return {x0 + (s - s0) / (s1 - s0) * (x1 - x0), y0 + (z0 - z) / (z0 - z1) * (y1 - y0)}; }
};

struct Ctx {
  const FigureContext& c;
  const SheetSpec& sp;
  Figure& f;
  double tw(const std::string& s, double size, bool bold = false) const { return text_width(s, size, bold, kFont); }
};

std::vector<double> ticks(double lo, double hi, double step) {
  std::vector<double> t;
  double a = std::ceil(lo / step - 1e-9) * step;
  for (double v = a; v <= hi + 1e-9 * step; v += step) t.push_back(std::fabs(v) < 1e-9 ? 0.0 : v);
  return t;
}
double nice_step(double span, int target) {
  double raw = span / target, mag = std::pow(10.0, std::floor(std::log10(raw)));
  for (double k : {1.0, 2.0, 2.5, 5.0, 10.0})
    if (k * mag >= raw) return k * mag;
  return 10 * mag;
}

void text_lines(Ctx& x, double px, double py, const std::string& t, double size, Color col, bool bold = false) {
  auto lines = split_lines(t);
  for (std::size_t k = 0; k < lines.size(); ++k)
    x.f.text({px, py + 1.2 * size * static_cast<double>(k)}, lines[k], size, Anchor::Start, bold, col);
}
struct Box { double x, y, w, h; };
Box text_box(Ctx& x, double px, double py, const std::string& t, double size) {
  auto lines = split_lines(t);
  double w = 0;
  for (auto& l : lines) w = std::max(w, x.tw(l, size));
  double h = size * (0.93 + 1.2 * static_cast<double>(lines.size() - 1)) + 0.24 * size;
  return {px, py - 0.93 * size, w, h};
}
/// Line from `target` to the edge of `b` along the direction to its centre.
Pt box_edge(const Box& b, Pt target) {
  Pt c{b.x + b.w / 2, b.y + b.h / 2};
  double dx = target.x - c.x, dy = target.y - c.y;
  double t = 1e30;
  if (std::fabs(dx) > 1e-9) t = std::min(t, (b.w / 2) / std::fabs(dx));
  if (std::fabs(dy) > 1e-9) t = std::min(t, (b.h / 2) / std::fabs(dy));
  if (t > 1) t = 1;
  return {c.x + dx * t, c.y + dy * t};
}

void draw_panel(Ctx& x, const Geometry& g, const Panel& pn, bool detail) {
  const Model& m = *x.c.model;
  const auto& inf = m.info();
  const std::size_t nh = inf.horizons.size(), n = g.st.size();
  Figure& f = x.f;
  f.clip(pn.x0, pn.y0, pn.x1 - pn.x0, pn.y1 - pn.y0);
  // blocks present
  std::vector<BlockCode> blocks;
  for (auto& col : g.iv)
    for (auto& iv : col)
      if (std::find(blocks.begin(), blocks.end(), iv.b) == blocks.end()) blocks.push_back(iv.b);
  std::sort(blocks.begin(), blocks.end());
  // unit fills, block by block, run by run
  for (BlockCode b : blocks) {
    std::size_t i = 0;
    while (i < n) {
      if (!block_interval(g, i, b)) { ++i; continue; }
      std::size_t j = i;
      while (j + 1 < n && block_interval(g, j + 1, b)) ++j;
      for (std::size_t u = 0; u <= nh; ++u) {
        // one polygon per stretch where the unit has thickness (a zero-thickness stretch would
        // leave a stroked hairline of this unit's colour along a block edge)
        std::vector<double> tz(j - i + 1), bz(j - i + 1);
        for (std::size_t k = i; k <= j; ++k) {
          const Interval& iv = *block_interval(g, k, b);
          double t = u == 0 ? iv.top : clampv(iv.h[u - 1], iv.bot, iv.top);
          double bt = u == nh ? iv.bot : clampv(iv.h[u], iv.bot, iv.top);
          tz[k - i] = t;
          bz[k - i] = std::min(bt, t);
        }
        auto thick = [&](std::size_t q) { return tz[q] - bz[q] > 1e-6; };
        std::size_t q = 0;
        const std::size_t nq = tz.size();
        while (q < nq) {
          if (!thick(q)) { ++q; continue; }
          std::size_t a0 = q > 0 ? q - 1 : q, a1 = q;
          while (a1 + 1 < nq && thick(a1 + 1)) ++a1;
          std::size_t e1 = a1 + 1 < nq ? a1 + 1 : a1;
          std::vector<Pt> poly;
          for (std::size_t r = a0; r <= e1; ++r) poly.push_back(pn.p(g.st[i + r], tz[r]));
          for (std::size_t r = e1 + 1; r-- > a0;) poly.push_back(pn.p(g.st[i + r], bz[r]));
          Color cc = unit_color(x.c, u);
          f.path(poly, true, cc, cc, 0.3);
          q = a1 + 1;
        }
      }
      i = j + 1;
    }
  }
  // horizon lines, clipped to each block's interval
  for (BlockCode b : blocks)
    for (std::size_t h = 0; h < nh; ++h) {
      std::vector<Pt> cur;
      auto flush = [&]() {
        if (cur.size() > 1) f.path(cur, false, none(), kLine, 0.45);
        cur.clear();
      };
      for (std::size_t k = 0; k + 1 < n; ++k) {
        const Interval* a = block_interval(g, k, b);
        const Interval* c = block_interval(g, k + 1, b);
        if (!a || !c || !std::isfinite(a->h[h]) || !std::isfinite(c->h[h])) { flush(); continue; }
        // visible part of the segment: bot(t) <= h(t) <= top(t), all linear in t
        double ha = a->h[h], hc = c->h[h];
        double t0 = 0, t1 = 1;
        auto keep = [&](double fa, double fc) {  // keep where fa + t (fc - fa) >= 0
          if (fa < 0 && fc < 0) { t0 = 1; t1 = 0; return; }
          if (fa >= 0 && fc >= 0) return;
          double t = fa / (fa - fc);
          if (fa < 0) t0 = std::max(t0, t); else t1 = std::min(t1, t);
        };
        keep(a->top - ha, c->top - hc);
        keep(ha - a->bot, hc - c->bot);
        if (t1 - t0 <= 1e-9) { flush(); continue; }
        Pt pa = pn.p(g.st[k] + t0 * (g.st[k + 1] - g.st[k]), ha + t0 * (hc - ha));
        Pt pc = pn.p(g.st[k] + t1 * (g.st[k + 1] - g.st[k]), ha + t1 * (hc - ha));
        if (cur.empty() || t0 > 1e-9) { flush(); cur.push_back(pa); }
        cur.push_back(pc);
        if (t1 < 1 - 1e-9) flush();
      }
      flush();
    }
  // faults between their tips
  for (std::size_t ff = 0; ff < inf.faults.size(); ++ff) {
    auto tr = fault_trace(g, static_cast<int>(ff));
    if (tr.size() < 2) continue;
    std::vector<Pt> pl;
    for (auto& [s, z] : tr) pl.push_back(pn.p(s, z));
    f.path(pl, false, none(), kBlack, 1.6);
  }
  // wells within 300 m of the line: path, then completions
  for (auto& w : x.c.wells) {
    double best = 1e30, bs = 0;
    for (std::size_t k = 0; k < n; ++k) {
      double d = std::hypot(g.xy[k].x - w.x(), g.xy[k].y - w.y());
      if (d < best) { best = d; bs = g.st[k]; }
    }
    if (best > 300) continue;
    std::vector<Pt> pl;
    for (auto& p : w.points()) pl.push_back(pn.p(bs, std::min(p.z, pn.z0 + 1000)));
    f.path(pl, false, none(), kBlack, detail ? 1.4 : 1.4);
    for (auto& c : x.sp.completions)
      if (c.well == w.name()) f.line(pn.p(bs, -c.top), pn.p(bs, -c.base), kCompletion, detail ? 3.4 : 3.0);
  }
  f.unclip();
}

void draw_axes(Ctx& x, const Panel& pn, double xstep, double zstep, double font) {
  Figure& f = x.f;
  for (double t : ticks(pn.s0, pn.s1, xstep)) {
    Pt p = pn.p(t, pn.z1);
    f.line(p, {p.x, p.y + 3.5}, kTick, 0.8);
    f.text({p.x, p.y + 3.5 + 3.5 + 0.73 * font}, fmt_dec(t, 0), font, Anchor::Middle, false, kTick);
  }
  for (double t : ticks(pn.z1, pn.z0, zstep)) {
    Pt p = pn.p(pn.s0, t);
    f.line({p.x - 3.5, p.y}, p, kTick, 0.8);
    f.text({p.x - 3.5 - 3.5, p.y + 0.36 * font}, minus(fmt_dec(t, 0)), font, Anchor::End, false, kTick);
  }
  f.rect(pn.x0, pn.y0, pn.x1 - pn.x0, pn.y1 - pn.y0, none(), kAxis, 0.8);
}

/// Mid-depth of unit u at column k (the interval holding most of it).
bool unit_mid(const Geometry& g, std::size_t k, std::size_t u, std::size_t nh, double& z) {
  double best = 0;
  bool ok = false;
  for (auto& iv : g.iv[k]) {
    double t = u == 0 ? iv.top : clampv(iv.h[u - 1], iv.bot, iv.top);
    double b = u == nh ? iv.bot : clampv(iv.h[u], iv.bot, iv.top);
    if (t - b > best) { best = t - b; z = 0.5 * (t + b); ok = true; }
  }
  return ok;
}
std::size_t col_at(const Geometry& g, double s) {
  std::size_t best = 0;
  for (std::size_t k = 1; k < g.st.size(); ++k)
    if (std::fabs(g.st[k] - s) < std::fabs(g.st[best] - s)) best = k;
  return best;
}
/// Elevation of horizon h at distance s (first interval that holds it).
bool horizon_at(const Geometry& g, double s, std::size_t h, double& z) {
  std::size_t k = col_at(g, s);
  for (auto& iv : g.iv[k]) {
    double r = iv.h[h];
    if (std::isfinite(r) && r <= iv.top && r >= iv.bot) { z = r; return true; }
  }
  return false;
}

}  // namespace

SheetResult section_sheet(const FigureContext& c, const SheetSpec& sp) {
  SheetResult res;
  Figure& f = res.figure;
  f.set_font_family(kFont);
  const Model& m = *c.model;
  const auto& inf = m.info();
  const std::size_t nh = inf.horizons.size();
  Ctx x{c, sp, f};
  const double L = polyline_length(sp.line);

  // ---- main panel geometry ----
  const double d_top = sp.main.d1 > sp.main.d0 ? sp.main.d0 : -inf.domain.z_top;
  const double d_bot = sp.main.d1 > sp.main.d0 ? sp.main.d1 : -inf.domain.z_base;
  const double ms0 = sp.main.s1 > sp.main.s0 ? sp.main.s0 : 0.0, ms1 = sp.main.s1 > sp.main.s0 ? sp.main.s1 : L;
  Geometry gm = build_geometry(m, sp.line, ms0, ms1, std::max(2.0, (ms1 - ms0) / 4000.0), -d_top, -d_bot);
  Panel pm;
  pm.x0 = 61.4;
  pm.x1 = 887.3;
  pm.y0 = 62.9;
  pm.s0 = ms0; pm.s1 = ms1; pm.z0 = -d_top; pm.z1 = -d_bot;
  pm.y1 = pm.y0 + (pm.x1 - pm.x0) * (d_bot - d_top) / (ms1 - ms0) * sp.main.ve;
  {  // keep the stated exaggeration: a panel too tall for the page is made narrower instead
    const double max_h = (sp.detail.present ? 393.7 : 600.0) - pm.y0;
    if (pm.y1 - pm.y0 > max_h) {
      pm.x1 = pm.x0 + (pm.x1 - pm.x0) * max_h / (pm.y1 - pm.y0);
      pm.y1 = pm.y0 + max_h;
    }
  }

  // ---- key numbers (computed before drawing so annotations can use them) ----
  auto well_path = [&](const std::string& n) -> const WellPath* {
    for (auto& w : c.wells)
      if (w.name() == n) return &w;
    return nullptr;
  };
  auto well_station = [&](const std::string& n, double& s) {
    const WellPath* w = well_path(n);
    if (!w) return false;
    double best = 1e30;
    for (std::size_t k = 0; k < gm.st.size(); ++k) {
      double d = std::hypot(gm.xy[k].x - w->x(), gm.xy[k].y - w->y());
      if (d < best) { best = d; s = gm.st[k]; }
    }
    return true;
  };
  auto horizon_at_well = [&](const std::string& wn, int h, Vec3& p) {
    const WellPath* w = well_path(wn);
    if (!w) return false;
    for (auto& t : repr::well_intersections(m, *w, 1.0))
      if (t.horizon_index == h) {
        // refine: the well is vertical or nearly so; take the exact root in its block
        p = t.p;
        return true;
      }
    return false;
  };
  auto block_of = [&](const Vec3& p) {
    EvalResult r;
    m.evaluate(std::span<const Vec3>(&p, 1), kEvalBlock, r);
    return r.block[0];
  };
  // Fault blocks lettered A, B, ... in the order they hold horizon h along the section.
  std::map<int, std::vector<BlockCode>> block_order;
  auto letter = [&](int hh, BlockCode b) -> std::string {
    auto it = block_order.find(hh);
    if (it == block_order.end()) {
      std::vector<BlockCode> o;
      for (std::size_t i = 0; i < gm.st.size(); ++i)
        for (auto& iv : gm.iv[i]) {
          double r = iv.h[static_cast<std::size_t>(hh)];
          if (std::isfinite(r) && r <= iv.top && r >= iv.bot && std::find(o.begin(), o.end(), iv.b) == o.end()) o.push_back(iv.b);
        }
      it = block_order.emplace(hh, std::move(o)).first;
    }
    auto p = std::find(it->second.begin(), it->second.end(), b);
    if (p == it->second.end()) return "?";
    std::size_t n = static_cast<std::size_t>(p - it->second.begin());
    return n < 26 ? std::string(1, static_cast<char>('A' + n)) : std::to_string(n + 1);
  };
  // Block on side hw of fault ff next to the fault surface at column i.
  auto side_block = [&](std::size_t i, int ff, bool hw, BlockCode& out) {
    double zf = gm.fz[i][static_cast<std::size_t>(ff)];
    if (std::isnan(zf)) return false;
    for (auto& iv : gm.iv[i]) {
      if ((((iv.b >> ff) & 1u) != 0) != hw) continue;
      if (std::fabs(iv.top - zf) < 1e-6 || std::fabs(iv.bot - zf) < 1e-6) { out = iv.b; return true; }
    }
    return false;
  };
  // Vertical throw of fault ff at horizon hh (from the footwall cut-off).
  auto throw_at_h = [&](int ff, int hh, double& out) {
    double s, z;
    if (!cutoff(gm, hh, ff, false, s, z)) return false;
    Vec3 xy = point_at(sp.line, s);
    const Vec3& nrm = inf.faults[static_cast<std::size_t>(ff)].normal;
    Vec3 q{xy.x - nrm.x * 0.01, xy.y - nrm.y * 0.01, z - nrm.z * 0.01};  // just on the footwall side
    FaultEval e;
    m.evaluate_fault(ff, std::span<const Vec3>(&q, 1), e);
    double vert = inf.faults[static_cast<std::size_t>(ff)].slip_mode == "dip_slip" ? std::sqrt(std::max(0.0, 1 - nrm.z * nrm.z)) : 1.0;
    out = e.displacement[0] * vert;
    return true;
  };
  // Thickness between horizons ha and hb where they meet fault ff on its footwall: the vertical
  // distance between their two footwall cut-offs (the interval the fault cuts).
  auto footwall_thickness = [&](int ff, int ha, int hb, double& out) {
    double sa, za, sb, zb;
    if (!cutoff(gm, ha, ff, false, sa, za) || !cutoff(gm, hb, ff, false, sb, zb)) return false;
    out = za - zb;
    return true;
  };
  auto unit_label = [&](std::size_t u) {
    return u < c.unit_labels.size() && !c.unit_labels[u].empty() ? c.unit_labels[u] : inf.units[u].name;
  };
  auto unit_at = [&](const Vec3& p) {
    BlockCode b = block_of(p);
    std::size_t u = 0;
    for (std::size_t hh = 0; hh < nh; ++hh) {
      double o = 0;
      m.scalar_in_block(std::span<const Vec3>(&p, 1), b, static_cast<int>(hh), std::span<double>(&o, 1));
      if (o < 0) u = hh + 1;  // below horizon hh
    }
    return u;
  };
  // Unfaulted elevation of horizon hh along the section columns.
  auto unfaulted_profile = [&](int hh) {
    std::vector<double> z(gm.st.size(), kNaN);
    if (sp.unfaulted)
      for (std::size_t i = 0; i < gm.st.size(); ++i) z[i] = sp.unfaulted(hh, gm.xy[i].x, gm.xy[i].y);
    return z;
  };
  auto compass = [](double az) {
    static const char* w[] = {"north", "north-east", "east", "south-east", "south", "south-west", "west", "north-west"};
    int k = static_cast<int>(std::floor(std::fmod(az + 360.0 + 22.5, 360.0) / 45.0)) % 8;
    return std::string(w[k]);
  };
  for (auto& k : sp.keys) {
    SheetValue v;
    v.label = k.label;
    v.on_sheet = k.on_sheet;
    std::map<std::string, double> bv;
    std::map<std::string, std::string> bw;
    const int h = k.horizon;
    if (k.kind == "horizon_at_well") {
      Vec3 p;
      if (horizon_at_well(k.well, h, p)) { v.value = -p.z; v.text = fmt_dec(v.value, k.decimals); }
    } else if (k.kind == "horizon_highest") {
      double best = -kInfD, best_s = 0;
      BlockCode best_b = 0;
      for (std::size_t i = 0; i < gm.st.size(); ++i)
        for (auto& iv : gm.iv[i]) {
          double r = iv.h[static_cast<std::size_t>(h)];
          if (std::isfinite(r) && r <= iv.top && r >= iv.bot && r > best) { best = r; best_s = gm.st[i]; best_b = iv.b; }
        }
      for (std::size_t ff = 0; ff < inf.faults.size(); ++ff)
        for (bool hw : {false, true}) {
          double s, z;
          BlockCode b;
          if (cutoff(gm, h, static_cast<int>(ff), hw, s, z) && z > best && side_block(col_at(gm, s), static_cast<int>(ff), hw, b)) {
            best = z; best_s = s; best_b = b;
          }
        }
      v.value = -best;
      v.text = fmt_dec(v.value, k.decimals);
      Vec3 xy = point_at(sp.line, best_s);
      bv = {{"x", xy.x}, {"y", xy.y}, {"s", best_s}};
      bw["block"] = letter(h, best_b);
    } else if (k.kind == "block_offset" && k.faults.size() == 2) {
      double sa = 0, za = 0, sb = 0, zb = 0;
      bool oka = cutoff(gm, h, k.faults[0], true, sa, za) || cutoff(gm, h, k.faults[0], false, sa, za);
      bool okb = cutoff(gm, h, k.faults[1], true, sb, zb) || cutoff(gm, h, k.faults[1], false, sb, zb);
      Vec3 pw;
      double zmid;
      if (oka && okb && horizon_at_well(k.well, h, pw) && horizon_at(gm, 0.5 * (sa + sb), static_cast<std::size_t>(h), zmid)) {
        Vec3 xy = point_at(sp.line, 0.5 * (sa + sb));
        Vec3 pa{xy.x, xy.y, zmid};
        Vec3 ra, rw;
        m.restore(std::span<const Vec3>(&pa, 1), std::span<Vec3>(&ra, 1));
        m.restore(std::span<const Vec3>(&pw, 1), std::span<Vec3>(&rw, 1));
        v.value = (pw.z - rw.z) - (pa.z - ra.z);
        v.text = fmt_dec(v.value, k.decimals);
      }
    } else if (k.kind == "fault_tips") {
      auto tr = fault_trace(gm, k.fault);
      if (tr.size() >= 2) {
        v.value = -tr.front().second;
        v.text = fmt_dec(-tr.front().second, k.decimals) + " / " + fmt_dec(-tr.back().second, k.decimals);
      }
    } else if (k.kind == "fault_throw") {
      if (throw_at_h(k.fault, h, v.value)) v.text = fmt_dec(v.value, k.decimals);
    } else if (k.kind == "well_fault_distance") {
      double sw = 0;
      std::vector<std::string> parts;
      if (well_station(k.well, sw)) {
        int ht = -1, hb = -1;
        if (k.unit >= 0 && static_cast<std::size_t>(k.unit) < inf.units.size()) {
          ht = inf.units[static_cast<std::size_t>(k.unit)].top_horizon;
          hb = inf.units[static_cast<std::size_t>(k.unit)].base_horizon;
        }
        for (int ff : k.faults) {
          double sc = kNaN;
          double prev = kNaN, prev_s = 0;
          for (std::size_t i = 0; i < gm.st.size() && ht >= 0 && hb >= 0; ++i) {
            double zf = gm.fz[i][static_cast<std::size_t>(ff)];
            if (std::isnan(zf)) { prev = kNaN; continue; }
            double zt, zb;
            if (sp.unfaulted) {
              zt = sp.unfaulted(ht, gm.xy[i].x, gm.xy[i].y);
              zb = sp.unfaulted(hb, gm.xy[i].x, gm.xy[i].y);
            } else {
              Vec3 q{gm.xy[i].x, gm.xy[i].y, zf};
              BlockCode b = block_of(q);
              const Interval* iv = block_interval(gm, i, b);
              if (!iv) { prev = kNaN; continue; }
              zt = iv->h[static_cast<std::size_t>(ht)];
              zb = iv->h[static_cast<std::size_t>(hb)];
            }
            double d = zf - (zt + k.fraction * (zb - zt));
            if (!std::isnan(prev) && ((prev > 0) != (d > 0))) {
              double t = prev / (prev - d);
              sc = prev_s + t * (gm.st[i] - prev_s);
              if (sp.unfaulted) {
                // refine on the exact fault surface and unfaulted level between the two columns
                auto gfun = [&](double sv) {
                  Vec3 xy = point_at(sp.line, sv);
                  double zg = zf;
                  auto dist = [&](double z) {
                    Vec3 q{xy.x, xy.y, z};
                    FaultEval e;
                    m.evaluate_fault(ff, std::span<const Vec3>(&q, 1), e);
                    return e.distance[0];
                  };
                  double lo = zg - 200, hi = zg + 200, dl = dist(lo), dh = dist(hi);
                  if ((dl > 0) == (dh > 0)) return kNaN;
                  double zfault = illinois(dist, lo, hi, dl, dh, 1e-7);
                  double a = sp.unfaulted(ht, xy.x, xy.y), bz = sp.unfaulted(hb, xy.x, xy.y);
                  return zfault - (a + k.fraction * (bz - a));
                };
                double ga = gfun(prev_s), gb = gfun(gm.st[i]);
                if (!std::isnan(ga) && !std::isnan(gb) && (ga > 0) != (gb > 0)) sc = illinois(gfun, prev_s, gm.st[i], ga, gb, 1e-6);
              }
              break;
            }
            prev = d;
            prev_s = gm.st[i];
          }
          parts.push_back(std::isnan(sc) ? "-" : fmt_dec(std::fabs(sw - sc), k.decimals));
          if (!std::isnan(sc) && v.value == 0) v.value = std::fabs(sw - sc);
        }
        std::string t;
        for (std::size_t i = 0; i < parts.size(); ++i) t += (i ? " / " : "") + parts[i];
        v.text = t;
      }
    } else if (k.kind == "block_highest_at_fault" || k.kind == "relief_to_end") {
      Vec3 pw;
      if (horizon_at_well(k.well, h, pw)) {
        BlockCode wb = block_of(pw);
        bool hw = ((wb >> k.fault) & 1u) != 0;
        double s, z;
        if (cutoff(gm, h, k.fault, hw, s, z)) {
          if (k.kind == "block_highest_at_fault") {
            v.value = -z;
          } else {
            std::size_t ce = k.end == "first" ? 0 : gm.st.size() - 1;
            const Interval* iv = block_interval(gm, ce, wb);
            double ze = iv ? iv->h[static_cast<std::size_t>(h)] : kNaN;
            v.value = z - ze;
          }
          v.text = fmt_dec(v.value, k.decimals);
        }
      }
    } else if (k.kind == "text") {
      v.text = k.value;
    } else if (k.kind == "completion") {
      for (auto& cp : sp.completions)
        if (cp.well == k.well && cp.name == k.completion) {
          v.value = cp.top;
          v.text = fmt_dec(cp.top, k.decimals >= 0 ? k.decimals : 0) + "-" + fmt_dec(cp.base, k.decimals >= 0 ? k.decimals : 0);
        }
    } else if (k.kind == "vertical_exaggeration") {
      v.value = sp.main.ve;
      v.text = fmt_dec(v.value, k.decimals);
    } else if (k.kind == "fold_relief" || k.kind == "limb_dip") {
      auto z = unfaulted_profile(h);
      std::size_t ic = gm.st.size();
      for (std::size_t i = 0; i < z.size(); ++i)
        if (!std::isnan(z[i]) && (ic == gm.st.size() || z[i] > z[ic])) ic = i;
      if (ic < gm.st.size()) {
        Vec3 xy = point_at(sp.line, gm.st[ic]);
        bv = {{"x", xy.x}, {"y", xy.y}, {"s", gm.st[ic]}};
        std::size_t ie = k.end == "first" ? 0 : z.size() - 1;
        if (k.kind == "fold_relief") {
          v.value = z[ic] - z[ie];
        } else {
          std::size_t a = std::min(ic, ie), b = std::max(ic, ie);
          double mx = 0;
          for (std::size_t i = a; i < b; ++i) {
            double ds = gm.st[i + 1] - gm.st[i];
            if (ds > 0 && !std::isnan(z[i]) && !std::isnan(z[i + 1])) mx = std::max(mx, std::fabs(z[i + 1] - z[i]) / ds);
          }
          v.value = std::atan(mx) * 180.0 / kPi;
        }
        v.text = fmt_dec(v.value, k.decimals);
      }
    } else if (k.kind == "fault_dip") {
      const Vec3& dd = inf.faults[static_cast<std::size_t>(k.fault)].dip_dir;
      double hz = std::hypot(dd.x, dd.y);
      v.value = std::atan2(std::fabs(dd.z), hz) * 180.0 / kPi;
      double az = std::atan2(dd.x, dd.y) * 180.0 / kPi;
      v.text = fmt_dec(v.value, k.decimals) + " deg towards " + compass(az);
    } else if (k.kind == "fault_tip_units") {
      auto tr = fault_trace(gm, k.fault);
      if (tr.size() >= 2) {
        Vec3 a = point_at(sp.line, tr.front().first), b = point_at(sp.line, tr.back().first);
        a.z = tr.front().second - 0.05;  // just below the upper tip
        b.z = tr.back().second + 0.05;   // just above the lower tip
        v.text = unit_label(unit_at(a)) + " / " + unit_label(unit_at(b));
      }
    } else if (k.kind == "fault_cutoff_position") {
      double s, z;
      if (cutoff(gm, h, k.fault, k.side == "hanging_wall", s, z)) {
        Vec3 xy = point_at(sp.line, s);
        v.value = k.coordinate == "x" ? xy.x : k.coordinate == "y" ? xy.y : s;
        v.text = fmt_dec(v.value, k.decimals);
        bv = {{"x", xy.x}, {"y", xy.y}, {"s", s}};
      }
    } else if (k.kind == "fault_cutoff_depths") {
      double sf, zf, sh, zh;
      BlockCode bf, bh;
      if (cutoff(gm, h, k.fault, false, sf, zf) && cutoff(gm, h, k.fault, true, sh, zh)) {
        v.value = -zf;
        v.text = fmt_dec(-zf, k.decimals) + " / " + fmt_dec(-zh, k.decimals);
        if (side_block(col_at(gm, sf), k.fault, false, bf)) bw["fw_block"] = letter(h, bf);
        if (side_block(col_at(gm, sh), k.fault, true, bh)) bw["hw_block"] = letter(h, bh);
      }
    } else if (k.kind == "fault_overlap" || k.kind == "throw_ratio") {
      double t, thk;
      int ha = k.from, hb = k.to;
      if (k.kind == "fault_overlap" && k.unit >= 0) {
        ha = inf.units[static_cast<std::size_t>(k.unit)].top_horizon;
        hb = inf.units[static_cast<std::size_t>(k.unit)].base_horizon;
      }
      if (ha >= 0 && hb >= 0 && throw_at_h(k.fault, h, t) && footwall_thickness(k.fault, ha, hb, thk) && std::fabs(thk) > 1e-9) {
        v.value = k.kind == "fault_overlap" ? 1.0 - std::fabs(t) / std::fabs(thk) : std::fabs(t) / std::fabs(thk);
        v.text = fmt_dec(v.value, k.decimals);
        bv = {{"throw", std::fabs(t)}, {"thickness", std::fabs(thk)}};
      }
    } else if (k.kind == "growth_index" && k.unit >= 0) {
      // Thorsen (1963): hanging-wall over footwall thickness of the unit at the fault (stage 2C)
      const int ha = inf.units[static_cast<std::size_t>(k.unit)].top_horizon, hb = inf.units[static_cast<std::size_t>(k.unit)].base_horizon;
      double ta, tb, thk;
      if (ha >= 0 && hb >= 0 && throw_at_h(k.fault, ha, ta) && throw_at_h(k.fault, hb, tb) && footwall_thickness(k.fault, ha, hb, thk) &&
          std::fabs(thk) > 1e-9) {
        v.value = (std::fabs(thk) + std::fabs(tb) - std::fabs(ta)) / std::fabs(thk);
        v.text = fmt_dec(v.value, k.decimals);
        bv = {{"throw_top", std::fabs(ta)}, {"throw_base", std::fabs(tb)}, {"footwall_thickness", std::fabs(thk)}};
      }
    } else if (k.kind == "fault_gap") {
      auto tr = fault_trace(gm, k.faults[0]);
      if (tr.size() >= 2) {
        const double za = tr.back().second, sa = tr.back().first;
        const std::size_t fb = static_cast<std::size_t>(k.faults[1]);
        double best = kNaN;
        for (std::size_t i = 0; i + 1 < gm.st.size(); ++i) {
          double a = gm.fz[i][fb] - za, b = gm.fz[i + 1][fb] - za;
          if (std::isnan(a) || std::isnan(b) || ((a > 0) == (b > 0) && a != 0)) continue;
          double t = a == b ? 0 : a / (a - b);
          double sb = gm.st[i] + t * (gm.st[i + 1] - gm.st[i]);
          if (std::isnan(best) || std::fabs(sb - sa) < best) best = std::fabs(sb - sa);
        }
        if (!std::isnan(best)) { v.value = best; v.text = fmt_dec(best, k.decimals); }
      }
    }
    if (v.text.empty()) v.text = "n/a";
    v.basis = fill(k.basis, bv, bw);
    res.keys.push_back(v);
  }

  // ---- title ----
  std::size_t units_shown = 0;
  for (std::size_t u = 0; u <= nh; ++u) {
    bool present = false;
    for (std::size_t k = 0; k < gm.st.size() && !present; k += 50) {
      double z;
      present = unit_mid(gm, k, u, nh, z);
    }
    units_shown += present ? 1 : 0;
  }
  std::map<std::string, double> tv{{"units", static_cast<double>(units_shown)},
                                   {"surfaces", static_cast<double>(units_shown + 1)},
                                   {"faults", static_cast<double>(inf.faults.size())},
                                   {"ve", sp.main.ve}};
  std::map<std::string, std::string> tw{{"faults_word", num_word(inf.faults.size())}};
  f.text({pm.x0, 57.1}, fill(sp.title, tv, tw), 11, Anchor::Start, false, kTitle);

  // ---- main panel ----
  draw_panel(x, gm, pm, false);
  draw_axes(x, pm, nice_step(ms1 - ms0, 8), nice_step(d_bot - d_top, 8), 8);
  f.text({0.5 * (pm.x0 + pm.x1), pm.y1 + 24.5}, sp.along_label, 8, Anchor::Middle, false, kTitle);
  f.text({21.4, 0.5 * (pm.y0 + pm.y1)}, sp.depth_label, 8, Anchor::Middle, false, kTitle, -90);
  if (!sp.end_left.empty()) f.text({pm.x0 + 3.2, pm.y0 + 10.6}, sp.end_left, 11, Anchor::Start, true, kTitle);
  if (!sp.end_right.empty()) f.text({pm.x1 - 3.2, pm.y0 + 10.6}, sp.end_right, 11, Anchor::End, true, kTitle);
  // unit labels
  auto labels = [&](const Geometry& g, const Panel& pn, const SheetPanel& spn, double size) {
    if (spn.label_at < 0) return;
    std::size_t k = col_at(g, spn.label_at);
    for (int u : spn.label_units) {
      double z;
      if (!unit_mid(g, k, static_cast<std::size_t>(u), nh, z)) continue;
      std::string t = static_cast<std::size_t>(u) < c.unit_labels.size() && !c.unit_labels[static_cast<std::size_t>(u)].empty()
                          ? c.unit_labels[static_cast<std::size_t>(u)]
                          : inf.units[static_cast<std::size_t>(u)].name;
      Pt p = pn.p(g.st[k], z);
      f.text({p.x, p.y + 0.36 * size}, t, size, Anchor::Middle, false, kTitle);
    }
  };
  labels(gm, pm, sp.main, 8);
  // fault labels at the shallow tip
  auto fault_labels = [&](const Geometry& g, const Panel& pn, bool boxed_default, double size) {
    for (std::size_t ff = 0; ff < inf.faults.size(); ++ff) {
      auto tr = fault_trace(g, static_cast<int>(ff));
      if (tr.size() < 2) continue;
      Pt p = pn.p(tr.front().first, tr.front().second);
      if (p.y > pn.y1 || p.x < pn.x0 || p.x > pn.x1) continue;
      bool boxed = boxed_default;
      const std::string& nm = inf.faults[ff].name;
      double w = x.tw(nm, size, true);
      double by = p.y - 4.0;  // baseline above the tip
      if (std::fabs(tr.front().second - g.z_top) < 1e-6 || by - 0.76 * size < pn.y0) {
        // the tip lies above the panel: label the fault where it enters, inside the frame
        by = pn.y0 + 0.76 * size + 4.0;
        boxed = true;
      }
      if (boxed) f.rect(p.x - w / 2 - 1.6, by - 0.76 * size - 1.4, w + 3.2, 0.76 * size + 2.8, Color{255, 255, 255, 178}, none(), 0);
      f.text({p.x, by}, nm, size, Anchor::Middle, true, kTitle);
    }
  };
  fault_labels(gm, pm, false, 9);
  // well labels
  for (auto& w : c.wells) {
    double s = 0;
    if (!well_station(w.name(), s)) continue;
    Pt p = pm.p(s, pm.z0);
    f.text({p.x + 4.0, pm.y0 + 8.9}, w.name(), 9, Anchor::Start, true, kTitle);
  }
  // annotations
  auto annotate = [&](const Geometry& g, const Panel& pn, const SheetPanel& spn, double size) {
    for (auto& a : spn.annotations) {
      Pt tp = pn.p(a.text_s, -a.text_depth);
      if (a.type == "leader") {
        double z;
        if (!horizon_at(g, a.at, static_cast<std::size_t>(a.horizon), z)) continue;
        Pt target = pn.p(a.at, z);
        Box b = text_box(x, tp.x, tp.y, a.text, size);
        f.line(box_edge(b, target), target, kTitle, 0.6);
        text_lines(x, tp.x, tp.y, a.text, size, kTitle);
      } else if (a.type == "horizon_at_fault") {
        double s, z;
        if (!cutoff(g, a.horizon, a.fault, !a.footwall, s, z)) continue;
        std::string t = fill(a.text, {{"value", -z}});
        Pt q = pn.p(s, z);
        Box b = text_box(x, tp.x, tp.y, t, size);
        f.line(q, {b.x, q.y}, kBlack, 1.0, {3.7, 1.6});
        f.rect(b.x - 2.0, b.y - 1.6, b.w + 4.0, b.h + 3.2, Color{255, 255, 255, 204}, none(), 0);
        text_lines(x, tp.x, tp.y, t, size, kTitle);
      } else if (a.type == "completion") {
        const WellPath* w = well_path(a.well);
        double s = 0;
        if (!w || !well_station(a.well, s)) continue;
        for (auto& cpl : sp.completions) {
          if (cpl.well != a.well || (!a.completion.empty() && cpl.name != a.completion)) continue;
          std::string t = fill(a.text, {{"top", cpl.top}, {"base", cpl.base}});
          Pt target = pn.p(s, -0.5 * (cpl.top + cpl.base));
          Box b = text_box(x, tp.x, tp.y, t, size);
          f.line(target, box_edge(b, target), kNote, 0.9);
          text_lines(x, tp.x, tp.y, t, size, kNote);
        }
      }
    }
  };
  annotate(gm, pm, sp.main, 8.2);

  // ---- detail panel ----
  if (sp.detail.present) {
    const auto& dp = sp.detail;
    Geometry gd = build_geometry(m, sp.line, dp.s0, dp.s1, std::max(1.0, (dp.s1 - dp.s0) / 3000.0), -dp.d0, -dp.d1);
    Panel pd;
    pd.x0 = 61.4;
    pd.x1 = 337.9;
    pd.y0 = 449.3;
    pd.s0 = dp.s0; pd.s1 = dp.s1; pd.z0 = -dp.d0; pd.z1 = -dp.d1;
    pd.y1 = pd.y0 + (pd.x1 - pd.x0) * (dp.d1 - dp.d0) / (dp.s1 - dp.s0) * dp.ve;
    if (pd.y1 > 650.0) {  // too tall: narrower, same exaggeration
      pd.x1 = pd.x0 + (pd.x1 - pd.x0) * (650.0 - pd.y0) / (pd.y1 - pd.y0);
      pd.y1 = 650.0;
    }
    std::map<std::string, double> dv = tv;
    dv["ve"] = dp.ve;
    dv["s0_km"] = dp.s0 / 1000.0;
    dv["s1_km"] = dp.s1 / 1000.0;
    f.text({pd.x0, 443.0}, fill(dp.title, dv, tw), 9.5, Anchor::Start, false, kTitle);
    draw_panel(x, gd, pd, true);
    draw_axes(x, pd, nice_step(dp.s1 - dp.s0, 6), nice_step(dp.d1 - dp.d0, 9), 8);
    f.text({0.5 * (pd.x0 + pd.x1), pd.y1 + 24.5}, sp.along_label, 8, Anchor::Middle, false, kTitle);
    f.text({21.4, 0.5 * (pd.y0 + pd.y1)}, sp.depth_label, 8, Anchor::Middle, false, kTitle, -90);
    labels(gd, pd, dp, 8);
    fault_labels(gd, pd, true, 9);
    annotate(gd, pd, dp, 7.5);
    // window outline on the main panel
    Pt a = pm.p(dp.s0, -dp.d0), b = pm.p(dp.s1, -dp.d1);
    f.path({a, {b.x, a.y}, b, {a.x, b.y}}, true, none(), kBlack, 1.1, {1.1, 1.8});
  }

  // ---- legend ----
  double legend_bottom = 61.2;
  {
    const double sx = 916.3, tx = 934.8, pitch = 10.3, size = 7.45;
    double y = 61.2;
    std::string group = "\x01";
    for (std::size_t u = 0; u < inf.units.size(); ++u) {
      std::string role = u < c.unit_roles.size() ? c.unit_roles[u] : std::string();
      if (role.empty()) continue;
      bool present = false;
      for (std::size_t k = 0; k < gm.st.size() && !present; k += 25) {
        double z;
        present = unit_mid(gm, k, u, nh, z);
      }
      if (!present) continue;
      if (role != group) {
        group = role;
        std::string up = role;
        for (auto& ch : up) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        f.text({tx, y}, up, size, Anchor::Start, true, kBlack);
        y += pitch;
      }
      std::string t = u < c.unit_labels.size() && !c.unit_labels[u].empty() ? c.unit_labels[u] : inf.units[u].name;
      f.rect(sx, y - 2.5 - 3.35, 12.0, 6.7, unit_color(c, u), rgb(30, 30, 30), 0.5);
      f.text({tx, y}, t, size, Anchor::Start, false, kBlack);
      y += pitch;
    }
    if (!sp.well_legend.empty() || !sp.completion_legend.empty()) {
      y += pitch;
      f.text({tx, y}, "WELL", size, Anchor::Start, true, kBlack);
      y += pitch;
      if (!sp.well_legend.empty()) {
        f.line({sx - 0.5, y - 2.5}, {sx + 12.5, y - 2.5}, kBlack, 1.4);
        f.text({tx, y}, sp.well_legend, size, Anchor::Start, false, kBlack);
        y += pitch;
      }
      if (!sp.completion_legend.empty()) {
        f.line({sx - 1.5, y - 2.5}, {sx + 13.5, y - 2.5}, kCompletion, 3.0);
        f.text({tx, y}, sp.completion_legend, size, Anchor::Start, false, kBlack);
        y += pitch;
      }
    }
    if (!sp.fault_legend.empty()) {
      y += pitch;
      f.text({tx, y}, "FAULTS", size, Anchor::Start, true, kBlack);
      y += pitch;
      f.line({sx - 0.5, y - 2.5}, {sx + 12.5, y - 2.5}, kBlack, 1.6);
      for (auto& l : sp.fault_legend) {
        f.text({tx, y}, l, size, Anchor::Start, false, kBlack);
        y += pitch;
      }
    }
    legend_bottom = y;
  }

  // ---- key numbers: bottom right, right-aligned to the legend column ----
  double key_top = f.height();
  {
    std::vector<const SheetValue*> shown;
    for (auto& v : res.keys)
      if (v.on_sheet) shown.push_back(&v);
    if (!shown.empty()) {
      const double size = 7.2, pitch = 10.4, head = 15.1, right = f.width() - 22.0, bottom = f.height() - 16.0;
      double w = x.tw(sp.key_title, 8, true);
      for (auto* v : shown) w = std::max(w, x.tw(v->label + ": " + v->text, size));
      double kx = right - w;
      double ky = bottom - pitch * static_cast<double>(shown.size() - 1) - head;
      key_top = ky - 8;
      f.text({kx, ky}, sp.key_title, 8, Anchor::Start, true, kTitle);
      ky += head;
      for (auto* v : shown) {
        f.text({kx, ky}, v->label + ": " + v->text, size, Anchor::Start, false, kTitle);
        ky += pitch;
      }
    }
  }
  // ---- notes: under the legend, above the key numbers ----
  if (!sp.notes.empty()) {
    double ny = std::max(legend_bottom + 14.0, 0.0);
    for (auto& l : sp.notes) {
      if (ny > key_top - 4) break;
      f.text({909.6, ny}, fill(l, {}, {{"figure", sp.name}}), 7.36, Anchor::Start, false, kGrey);
      ny += 8.8;
    }
  }
  return res;
}

}  // namespace stratonix::fig
