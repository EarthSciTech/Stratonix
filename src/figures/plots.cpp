// Figure set: sections, maps, 3D views, fence diagram; writing in four formats.
#include <algorithm>
#include <set>
#include <cmath>
#include <map>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/figures/figures.hpp"

namespace stratonix::fig {

namespace {
const Color kInk = rgb(30, 30, 30), kGrid = rgb(225, 225, 225), kFault = rgb(170, 20, 20), kWell = rgb(0, 0, 0);
const char* kPalette[] = {"#e8eef4", "#f2d39b", "#a9cfa4", "#e5a98c", "#9ec5e0", "#d9c2e6", "#f0e68c",
                          "#c8b98f", "#b5d6d0", "#e8b4c8", "#cfd8a0", "#b9b4aa"};

std::vector<double> nice_ticks(double lo, double hi, int target) {
  double span = hi - lo;
  if (!(span > 0)) return {lo};
  double raw = span / target, mag = std::pow(10, std::floor(std::log10(raw)));
  double step = mag;
  for (double m : {1.0, 2.0, 2.5, 5.0, 10.0})
    if (mag * m >= raw) { step = mag * m; break; }
  std::vector<double> t;
  for (double v = std::ceil(lo / step) * step; v <= hi + 1e-9 * span; v += step) t.push_back(std::fabs(v) < 1e-9 * span ? 0 : v);
  return t;
}
std::string km_or_m(double v) { return fmt_fixed(v, 0); }

void footer(Figure& f, const FigureContext& c) {
  const auto& inf = c.model->info();
  double y = f.height() - 12;
  f.line({18, y - 11}, {f.width() - 18, y - 11}, rgb(180, 180, 180), 0.5);
  std::string s = "Coordinate system: " + inf.crs.label() + "   |   Vertical datum: " + inf.crs.vertical_datum +
                  "   |   Evidence basis: " + basis_label(inf.basis) + "   |   " + c.model_name + " v" + c.model_version +
                  "   |   Stratonix " + build_info().version;
  double size = 6.5;
  while (text_width(s, size) > f.width() - 36 && size > 4.5) size -= 0.25;
  f.text({18, y}, s, size, Anchor::Start, false, rgb(80, 80, 80));
  if (inf.basis == EvidenceBasis::Design)
    f.text({f.width() - 18, 20}, "DESIGN MODEL - not an interpretation", 7.5, Anchor::End, true, rgb(170, 20, 20));
}

void title(Figure& f, const std::string& t, const std::string& sub) {
  double size = 13;
  while (text_width(t, size, true) > f.width() - 260 && size > 9) size -= 0.5;
  f.text({18, 24}, t, size, Anchor::Start, true, kInk);
  f.text({18, 38}, sub, 8, Anchor::Start, false, rgb(70, 70, 70));
}

/// Legend entries drawn top-down from (x, y); returns the next y.
struct LegendEntry { std::string group, label; Color fill; bool line = false; Color line_color; double width = 1; int marker = 0; };  // marker 1 diamond, 2 triangle, 3 cross
void legend_symbol(Figure& f, double x, double y, const LegendEntry& e) {
  if (e.marker == 1) f.path({{x + 7, y + 4}, {x + 10, y + 1}, {x + 13, y + 4}, {x + 10, y + 7}}, true, e.fill, kWell, 0.9);
  else if (e.marker == 2) f.path({{x + 10, y}, {x + 14, y + 7}, {x + 6, y + 7}}, true, e.fill, none());
  else if (e.marker == 3) { f.line({x + 7, y + 1}, {x + 13, y + 7}, e.line_color, 1.4); f.line({x + 7, y + 7}, {x + 13, y + 1}, e.line_color, 1.4); }
  else if (e.line) f.line({x + 2, y + 4}, {x + 18, y + 4}, e.line_color, e.width);
  else f.rect(x + 2, y, 16, 8, e.fill, rgb(90, 90, 90), 0.3);
}

/// Legend from (x, y) down, grouped; when one column would pass max_y, entries are set two to a
/// row in a smaller font. Returns the next y.
double legend(Figure& f, double x, double y, double w, const std::vector<LegendEntry>& es, double max_y = 1e9) {
  std::size_t groups = 0;
  std::string g0;
  for (auto& e : es)
    if (e.group != g0) { g0 = e.group; ++groups; }
  const bool two = y + 16.0 * static_cast<double>(groups) + 11.0 * static_cast<double>(es.size()) > max_y;
  const double cw = two ? w / 2 : w, row = two ? 9.5 : 11, font = two ? 6.2 : 7;
  std::string group;
  int col = 0;
  for (auto& e : es) {
    if (e.group != group) {
      if (col == 1) { y += row; col = 0; }
      group = e.group;
      y += 4;
      f.text({x, y + 8}, group, 7.5, Anchor::Start, true, kInk);
      y += 12;
    }
    double xx = x + col * cw;
    legend_symbol(f, xx, y, e);
    std::string lab = e.label;
    while (text_width(lab, font) > cw - 26 && lab.size() > 4) lab = lab.substr(0, lab.size() - 2);
    f.text({xx + 24, y + 7}, lab, font, Anchor::Start, false, kInk);
    if (two && col == 0) col = 1;
    else { col = 0; y += row; }
  }
  if (col == 1) y += row;
  return y;
}

double key_box_height(const FigureContext& c) { return c.key_numbers.empty() ? 0 : 16 + 10.5 * static_cast<double>(c.key_numbers.size()); }

double key_box(Figure& f, const FigureContext& c, double x, double y, double w) {
  if (c.key_numbers.empty()) return y;
  double h = 16 + 10.5 * static_cast<double>(c.key_numbers.size());
  f.rect(x, y, w, h, rgb(248, 248, 244), rgb(150, 150, 150), 0.5);
  f.text({x + 6, y + 11}, "Key numbers", 7.5, Anchor::Start, true, kInk);
  double yy = y + 23;
  for (auto& [k, v] : c.key_numbers) {
    std::string lab = k;
    while (text_width(lab + ": " + v, 6.5) > w - 10 && lab.size() > 6) lab = lab.substr(0, lab.size() - 2);
    f.text({x + 6, yy}, lab + ": " + v, 6.5, Anchor::Start, false, kInk);
    yy += 10.5;
  }
  return y + h;
}

std::vector<LegendEntry> unit_legend(const FigureContext& c) {
  std::vector<LegendEntry> es;
  const auto& inf = c.model->info();
  for (std::size_t u = 0; u < inf.units.size(); ++u) {
    std::string role = u < c.unit_roles.size() && !c.unit_roles[u].empty() ? c.unit_roles[u] : "Units";
    es.push_back({role, inf.units[u].name, unit_color(c, u)});
  }
  // groups in order of first appearance (stratigraphic order); units without a role last
  std::vector<std::string> order;
  for (auto& e : es)
    if (e.group != "Units" && std::find(order.begin(), order.end(), e.group) == order.end()) order.push_back(e.group);
  order.push_back("Units");
  auto rank = [&](const std::string& g) { return std::find(order.begin(), order.end(), g) - order.begin(); };
  std::stable_sort(es.begin(), es.end(), [&](auto& a, auto& b) { return rank(a.group) < rank(b.group); });
  return es;
}
}  // namespace

std::string figure_basename(const FigureContext& c, const std::string& kind, const std::string& name) {
  std::string s = "STX_" + c.model_name + "_v" + c.model_version + "_" + kind + (name.empty() ? "" : "_" + name);
  for (auto& ch : s)
    if (!(std::isalnum(static_cast<unsigned char>(ch)) || ch == '-' || ch == '_' || ch == '.')) ch = '-';
  return s;
}

Color unit_color(const FigureContext& c, std::size_t u) {
  if (u < c.unit_colors.size() && !c.unit_colors[u].empty()) return Color::hex(c.unit_colors[u]);
  const auto& inf = c.model->info();
  if (u == inf.horizons.size()) return Color::hex(kPalette[11]);  // below the base horizon (bodies come after it)
  return Color::hex(kPalette[u % 11]);
}

// =========================================================================================
// Section
// =========================================================================================
namespace {
struct Frame2 {
  double x0, y0, x1, y1;   // plot rectangle (pt)
  double s0, s1, z0, z1;   // data window: distance, elevation (z0 top > z1 bottom)
  Pt p(double s, double z) const {
    return {x0 + (s - s0) / (s1 - s0) * (x1 - x0), y0 + (z0 - z) / (z0 - z1) * (y1 - y0)};
  }
};

void draw_section_panel(Figure& f, const FigureContext& c, const repr::Section& s, const Frame2& fr, bool axes_labels, double font) {
  const auto& inf = c.model->info();
  const std::size_t nh = inf.horizons.size();
  f.clip(fr.x0, fr.y0, fr.x1 - fr.x0, fr.y1 - fr.y0);
  // Unit fills, strip by strip.
  for (std::size_t col = 0; col + 1 < s.ns; ++col) {
    double sa = s.station[col], sb = s.station[col + 1];
    if (sb < fr.s0 || sa > fr.s1) continue;
    bool smooth = true;
    for (std::size_t h = 0; h < nh && smooth; ++h)
      smooth = !std::isnan(s.horizon_z[h][col]) && !std::isnan(s.horizon_z[h][col + 1]);
    for (std::size_t ff = 0; ff < s.fault_z.size() && smooth; ++ff)
      smooth = std::isnan(s.fault_z[ff][col]) && std::isnan(s.fault_z[ff][col + 1]);
    if (smooth) {
      for (std::size_t u = 0; u <= nh; ++u) {
        double ta = u == 0 ? fr.z0 : s.horizon_z[u - 1][col], tb = u == 0 ? fr.z0 : s.horizon_z[u - 1][col + 1];
        double ba = u == nh ? fr.z1 : s.horizon_z[u][col], bb = u == nh ? fr.z1 : s.horizon_z[u][col + 1];
        ta = std::clamp(ta, fr.z1, fr.z0); tb = std::clamp(tb, fr.z1, fr.z0);
        ba = std::clamp(ba, fr.z1, fr.z0); bb = std::clamp(bb, fr.z1, fr.z0);
        if (ta - ba < 1e-6 && tb - bb < 1e-6) continue;
        Color cc = unit_color(c, u);
        f.path({fr.p(sa, ta), fr.p(sb, tb), fr.p(sb, bb), fr.p(sa, ba)}, true, cc, cc, 0.25);
      }
    } else {
      // run-length rectangles from the section's unit grid (column col)
      std::size_t r = 0;
      while (r < s.nz) {
        std::int32_t u = s.unit[r * s.ns + col];
        std::size_t r2 = r;
        while (r2 + 1 < s.nz && s.unit[(r2 + 1) * s.ns + col] == u) ++r2;
        double zt = s.z_top - (static_cast<double>(r) - 0.5) * s.dz, zb = s.z_top - (static_cast<double>(r2) + 0.5) * s.dz;
        zt = std::clamp(zt, fr.z1, fr.z0);
        zb = std::clamp(zb, fr.z1, fr.z0);
        if (zt > zb) {
          Color cc = unit_color(c, static_cast<std::size_t>(u));
          f.path({fr.p(sa, zt), fr.p(sb, zt), fr.p(sb, zb), fr.p(sa, zb)}, true, cc, cc, 0.25);
        }
        r = r2 + 1;
      }
    }
  }
  // stage 4B: bodies (units after the stratigraphic column) over the fills, from the unit grid
  if (inf.units.size() > nh + 1)
    for (std::size_t col = 0; col + 1 < s.ns; ++col) {
      const double sa = s.station[col], sb = s.station[col + 1];
      if (sb < fr.s0 || sa > fr.s1) continue;
      std::size_t r = 0;
      while (r < s.nz) {
        const std::int32_t u = s.unit[r * s.ns + col];
        std::size_t r2 = r;
        while (r2 + 1 < s.nz && s.unit[(r2 + 1) * s.ns + col] == u) ++r2;
        if (u > static_cast<std::int32_t>(nh)) {
          double zt = std::clamp(s.z_top - (static_cast<double>(r) - 0.5) * s.dz, fr.z1, fr.z0);
          double zb = std::clamp(s.z_top - (static_cast<double>(r2) + 0.5) * s.dz, fr.z1, fr.z0);
          if (zt > zb) {
            const Color cc = unit_color(c, static_cast<std::size_t>(u));
            f.path({fr.p(sa, zt), fr.p(sb, zt), fr.p(sb, zb), fr.p(sa, zb)}, true, cc, cc, 0.25);
          }
        }
        r = r2 + 1;
      }
    }
  // Horizon lines (per fault block where blocks overlap, so both sheets are drawn).
  auto draw_line = [&](const std::vector<double>& zz) {
    std::vector<Pt> cur;
    for (std::size_t col = 0; col < s.ns; ++col) {
      double z = zz[col];
      if (std::isnan(z)) { if (cur.size() > 1) f.path(cur, false, none(), kInk, 0.6); cur.clear(); continue; }
      cur.push_back(fr.p(s.station[col], z));
    }
    if (cur.size() > 1) f.path(cur, false, none(), kInk, 0.6);
  };
  for (std::size_t h = 0; h < nh; ++h) {
    if (h < s.horizon_block_z.size() && !s.horizon_block_z[h].empty())
      for (auto& [bc, zz] : s.horizon_block_z[h]) draw_line(zz);
    else
      draw_line(s.horizon_z[h]);
  }
  // Fault traces (between tips).
  for (std::size_t ff = 0; ff < s.fault_z.size(); ++ff) {
    std::vector<Pt> cur;
    for (std::size_t col = 0; col < s.ns; ++col) {
      double z = s.fault_z[ff][col];
      if (std::isnan(z)) { if (cur.size() > 1) f.path(cur, false, none(), kFault, 1.3); cur.clear(); continue; }
      cur.push_back(fr.p(s.station[col], z));
    }
    if (cur.size() > 1) f.path(cur, false, none(), kFault, 1.3);
  }
  // Wells within 300 m of the line, projected.
  for (auto& w : c.wells) {
    double best = 1e30, bs = 0;
    for (std::size_t col = 0; col < s.ns; ++col) {
      double d = std::hypot(s.xy[col].x - w.x(), s.xy[col].y - w.y());
      if (d < best) { best = d; bs = s.station[col]; }
    }
    if (best > 300) continue;
    const auto& pts = w.points();
    std::vector<Pt> pl;
    for (auto& p : pts) pl.push_back(fr.p(bs, p.z));
    f.path(pl, false, none(), kWell, 1.0);
    auto tops = repr::well_intersections(*c.model, w, 2.0);
    for (auto& t : tops) {
      Pt q = fr.p(bs, t.p.z);
      f.path({{q.x - 2.2, q.y}, {q.x + 2.2, q.y}}, false, none(), kWell, 1.4);
    }
    Pt top = fr.p(bs, std::min(fr.z0, pts.front().z));
    f.text({top.x, top.y - 3}, w.name(), font, Anchor::Middle, true, kWell);
  }
  f.unclip();
  // Frame and ticks.
  f.rect(fr.x0, fr.y0, fr.x1 - fr.x0, fr.y1 - fr.y0, none(), kInk, 0.6);
  for (double t : nice_ticks(fr.s0, fr.s1, 6)) {
    Pt p = fr.p(t, fr.z1);
    f.line(p, {p.x, p.y + 3}, kInk, 0.5);
    f.text({p.x, p.y + 11}, km_or_m(t), font, Anchor::Middle, false, kInk);
  }
  for (double t : nice_ticks(-fr.z0, -fr.z1, 6)) {
    Pt p = fr.p(fr.s0, -t);
    f.line({p.x - 3, p.y}, p, kInk, 0.5);
    f.text({p.x - 5, p.y + 2.5}, km_or_m(t), font, Anchor::End, false, kInk);
  }
  if (axes_labels) {
    f.text({0.5 * (fr.x0 + fr.x1), fr.y1 + 22}, "Distance along section (m)", font + 0.5, Anchor::Middle, false, kInk);
    f.text({fr.x0 - 42, 0.5 * (fr.y0 + fr.y1)}, "Depth below datum (m)", font + 0.5, Anchor::Middle, false, kInk, -90);
  }
}
}  // namespace

Figure section_figure(const FigureContext& c, const repr::Section& s, const std::vector<double>* dw) {
  Figure f(842, 540);
  const auto& inf = c.model->info();
  // Elevation window: horizons +/- a margin, inside the domain.
  double zmax = -1e30, zmin = 1e30;
  for (auto& hz : s.horizon_z)
    for (double z : hz)
      if (!std::isnan(z)) { zmax = std::max(zmax, z); zmin = std::min(zmin, z); }
  double span = std::max(100.0, zmax - zmin);
  double z0 = std::min(inf.domain.z_top, zmax + 0.35 * span), z1 = std::max(inf.domain.z_base, zmin - 0.35 * span);
  double x0 = 70, x1 = 842 - 232, y0 = 58;
  // Height from the requested vertical exaggeration, limited by the page.
  double h_ve = (x1 - x0) * (z0 - z1) / s.length * c.vertical_exaggeration;
  double y1 = y0 + std::clamp(h_ve, 120.0, 540.0 - 58 - 92);
  double ve = (y1 - y0) / (x1 - x0) * s.length / (z0 - z1);
  Frame2 fr{x0, y0, x1, y1, 0, s.length, z0, z1};
  std::string sub = "Section " + s.name + " from (" + fmt_fixed(s.line.front().x, 0) + ", " + fmt_fixed(s.line.front().y, 0) + ") to (" +
                    fmt_fixed(s.line.back().x, 0) + ", " + fmt_fixed(s.line.back().y, 0) + "); vertical exaggeration " + fmt_fixed(ve, 1) + "x";
  title(f, (c.title.empty() ? c.model_name : c.title) + " - section " + s.name, sub);
  draw_section_panel(f, c, s, fr, true, 7);
  // Bends of the section line: dashed verticals at interior vertices.
  {
    double acc = 0;
    for (std::size_t v = 1; v + 1 < s.line.size(); ++v) {
      acc += std::hypot(s.line[v].x - s.line[v - 1].x, s.line[v].y - s.line[v - 1].y);
      Pt a = fr.p(acc, z0), b = fr.p(acc, z1);
      f.line(a, b, rgb(90, 90, 90), 0.5, {3, 2});
      f.text({a.x, a.y - 4}, "bend", 6, Anchor::Middle, false, rgb(90, 90, 90));
    }
  }
  f.text({x0, y0 - 4}, "A", 9, Anchor::Middle, true, kInk);
  f.text({x1, y0 - 4}, "A'", 9, Anchor::Middle, true, kInk);
  // Fault labels at the shallow end of each trace.
  for (std::size_t ff = 0; ff < s.fault_z.size(); ++ff) {
    double bz = -1e30;
    std::size_t bc = 0;
    for (std::size_t col = 0; col < s.ns; ++col)
      if (!std::isnan(s.fault_z[ff][col]) && s.fault_z[ff][col] > bz && s.fault_z[ff][col] <= z0) { bz = s.fault_z[ff][col]; bc = col; }
    if (bz > -1e29) {
      Pt p = fr.p(s.station[bc], bz);
      f.text({p.x + 3, p.y - 3}, inf.faults[ff].name, 7, Anchor::Start, true, kFault);
    }
  }
  // Right column: legend, key numbers, detail panel.
  double rx = 842 - 214, rw = 196, ry = 50;
  auto es = unit_legend(c);
  es.push_back({"Lines", "Horizon", none(), true, kInk, 0.6});
  es.push_back({"Lines", "Fault (between tips)", none(), true, kFault, 1.3});
  if (!c.wells.empty()) es.push_back({"Lines", "Well and tops", none(), true, kWell, 1.0});
  ry = legend(f, rx, ry, rw, es, 540.0 - 30 - key_box_height(c) - 8 - (dw && dw->size() == 4 ? 80 : 0)) + 6;
  ry = key_box(f, c, rx, ry, rw) + 8;
  if (dw && dw->size() == 4) {
    double ds0 = (*dw)[0], ds1 = (*dw)[1], dz0 = -(*dw)[2], dz1 = -(*dw)[3];
    Frame2 det{rx + 30, ry + 12, rx + rw, std::min(540.0 - 60, ry + 12 + 150), ds0, ds1, dz0, dz1};
    if (det.y1 - det.y0 > 60) {
      f.text({rx, ry + 7}, "Detail", 7.5, Anchor::Start, true, kInk);
      draw_section_panel(f, c, s, det, false, 5.5);
      // outline of the detail window on the main panel
      Pt a = fr.p(ds0, dz0), b = fr.p(ds1, dz1);
      f.rect(a.x, a.y, b.x - a.x, b.y - a.y, none(), rgb(20, 60, 160), 0.8);
    }
  }
  footer(f, c);
  return f;
}

// =========================================================================================
// Map
// =========================================================================================
namespace {
Color ramp(double t) {  // shallow (t = 0) warm, deep (t = 1) cool
  static const double stops[][3] = {{254, 240, 190}, {240, 200, 120}, {180, 210, 150}, {110, 180, 190}, {60, 110, 170}, {40, 60, 120}};
  t = std::clamp(t, 0.0, 1.0) * 5;
  int i = std::min(4, static_cast<int>(t));
  double u = t - i;
  return rgb(static_cast<int>(stops[i][0] + u * (stops[i + 1][0] - stops[i][0])), static_cast<int>(stops[i][1] + u * (stops[i + 1][1] - stops[i][1])),
             static_cast<int>(stops[i][2] + u * (stops[i + 1][2] - stops[i][2])));
}
}  // namespace

Color termination_color(const std::string& kind) {
  if (kind == "onlap") return rgb(0, 90, 200);
  if (kind == "truncation") return rgb(205, 105, 0);
  if (kind == "pinch-out") return rgb(130, 40, 160);
  if (kind == "salt") return rgb(200, 30, 120);
  return rgb(0, 130, 130);  // clinoform
}
std::string termination_label(const std::string& kind) {
  if (kind == "truncation") return "Subcrop (erosional truncation) line";
  if (kind == "onlap") return "Onlap line";
  if (kind == "pinch-out") return "Pinch-out line";
  if (kind == "salt") return "Salt piercement line";
  return "Clinoform termination line";
}

Figure map_figure(const FigureContext& c, int h, const repr::FaultedSurface& surf,
                  const std::vector<std::vector<repr::CutoffLine>>& cutoffs, const repr::Closure* cl,
                  const std::vector<repr::TerminationLine>* terms) {
  Figure f(842, 595);
  const auto& inf = c.model->info();
  const auto& d = inf.domain;
  const Grid2D& g = surf.z;
  double side = std::min(595.0 - 110, 842.0 - 280);
  double sc = side / std::max(d.sx, d.sy);
  double mx0 = 60, my0 = 60, mw = d.sx * sc, mh = d.sy * sc;
  auto P = [&](double x, double y) { return Pt{mx0 + (x - d.x0) * sc, my0 + mh - (y - d.y0) * sc}; };
  double zmin = 1e30, zmax = -1e30;
  for (double v : g.v)
    if (!std::isnan(v)) { zmin = std::min(zmin, v); zmax = std::max(zmax, v); }
  const bool flat = !(zmax - zmin >= 0.5);  // relief below 0.5 m: one colour, no contours
  title(f, (c.title.empty() ? c.model_name : c.title) + " - " + inf.horizons[static_cast<std::size_t>(h)].name + " depth map",
        "Depth below datum (m), contour interval shown in the colour bar; fault polygons from hanging-wall and footwall cut-offs");
  // Shaded depth image.
  Image im;
  // Image resolution independent of the grid (about 3 pixels per point at the drawn size).
  im.w = std::max(1, static_cast<int>(std::lround(mw * 3)));
  im.h = std::max(1, static_cast<int>(std::lround(mh * 3)));
  im.rgba.resize(static_cast<std::size_t>(im.w) * im.h * 4);
  // Pixel depths: bilinear within one fault block. Where the four nodes are not all defined in
  // one block (near a fault), each block present among them is extrapolated to the pixel and
  // kept only if the model puts that point in that block; a pixel no block claims is in the
  // fault gap. This draws the gap where the model has it, not where the grid cell is.
  std::vector<double> pz(static_cast<std::size_t>(im.w) * im.h, std::nan(""));
  std::vector<Vec3> cand;
  std::vector<std::pair<std::size_t, BlockCode>> cand_of;
  const auto& gs = g.spec;
  for (int yy = 0; yy < im.h; ++yy)
    for (int xx = 0; xx < im.w; ++xx) {
      double x = d.x0 + (xx + 0.5) / im.w * d.sx, y = d.y0 + d.sy - (yy + 0.5) / im.h * d.sy;
      double fx = std::clamp((x - gs.x0) / gs.dx, 0.0, static_cast<double>(gs.nx - 1));
      double fy = std::clamp((y - gs.y0) / gs.dy, 0.0, static_cast<double>(gs.ny - 1));
      std::size_t i0 = std::min<std::size_t>(static_cast<std::size_t>(fx), gs.nx > 1 ? gs.nx - 2 : 0);
      std::size_t j0 = std::min<std::size_t>(static_cast<std::size_t>(fy), gs.ny > 1 ? gs.ny - 2 : 0);
      double u = fx - i0, v = fy - j0;
      std::size_t k[4] = {gs.index(i0, j0), gs.index(i0 + 1, j0), gs.index(i0, j0 + 1), gs.index(i0 + 1, j0 + 1)};
      double wt[4] = {(1 - u) * (1 - v), u * (1 - v), (1 - u) * v, u * v};
      std::size_t pix = static_cast<std::size_t>(yy) * im.w + xx;
      bool all = true;
      for (int c4 = 0; c4 < 4 && all; ++c4) all = !std::isnan(g.v[k[c4]]) && surf.block[k[c4]] == surf.block[k[0]];
      if (all) {
        double z = 0;
        for (int c4 = 0; c4 < 4; ++c4) z += wt[c4] * g.v[k[c4]];
        pz[pix] = z;
        continue;
      }
      std::vector<BlockCode> seen;
      for (int c4 = 0; c4 < 4; ++c4) {
        if (std::isnan(g.v[k[c4]])) continue;
        BlockCode bc = surf.block[k[c4]];
        if (std::find(seen.begin(), seen.end(), bc) != seen.end()) continue;
        seen.push_back(bc);
        // weighted mean of this block's nodes plus a linear extension along the cell edges
        double zs = 0, ws = 0;
        for (int c5 = 0; c5 < 4; ++c5)
          if (!std::isnan(g.v[k[c5]]) && surf.block[k[c5]] == bc) { zs += wt[c5] * g.v[k[c5]]; ws += wt[c5]; }
        double z = ws > 1e-12 ? zs / ws : g.v[k[c4]];
        cand.push_back({x, y, z});
        cand_of.push_back({pix, bc});
      }
    }
  if (!cand.empty()) {
    EvalResult er;
    c.model->evaluate(cand, kEvalBlock, er);
    for (std::size_t q = 0; q < cand.size(); ++q) {
      if (er.block[q] != cand_of[q].second) continue;
      double& t = pz[cand_of[q].first];
      if (std::isnan(t) || cand[q].z > t) t = cand[q].z;  // overlap (reverse fault): upper kept
    }
  }
  for (int yy = 0; yy < im.h; ++yy)
    for (int xx = 0; xx < im.w; ++xx) {
      double v = pz[static_cast<std::size_t>(yy) * im.w + xx];
      std::uint8_t* p = &im.rgba[(static_cast<std::size_t>(yy) * im.w + xx) * 4];
      if (std::isnan(v)) { p[0] = p[1] = p[2] = 235; p[3] = 255; continue; }
      Color cc = flat ? ramp(0.5) : ramp((zmax - v) / std::max(1.0, zmax - zmin));
      p[0] = cc.r; p[1] = cc.g; p[2] = cc.b; p[3] = 255;
    }
  const Image* img = f.own(std::move(im));
  f.clip(mx0, my0, mw, mh);
  f.image(mx0, my0, mw, mh, img);
  // Contours.
  auto ticks = flat ? std::vector<double>{} : nice_ticks(-zmax, -zmin, 12);
  double ci = ticks.size() > 1 ? ticks[1] - ticks[0] : 10;
  int nlab = 0;
  std::vector<Pt> placed;
  for (std::size_t k = 0; k < ticks.size(); ++k) {
    double level = -ticks[k];
    bool major = static_cast<long long>(std::llround(ticks[k] / ci)) % 2 == 0;
    for (auto& line : repr::contour_lines(g, level)) {
      if (line.size() < 5) continue;  // fragments at fault gaps
      std::vector<Pt> pl;
      for (auto& p : line) pl.push_back(P(p.x, p.y));
      f.path(pl, false, none(), major ? rgb(60, 60, 60) : rgb(120, 120, 120), major ? 0.55 : 0.3);
      if (major && line.size() > 40 && nlab < 14) {
        Pt q = pl[pl.size() / 2];
        bool clash = false;
        for (auto& w : c.wells) {
          Pt wp = P(w.x(), w.y());
          if (std::fabs(wp.x - q.x) < 30 && std::fabs(wp.y - q.y) < 12) clash = true;
        }
        for (auto& o : placed)
          if (std::fabs(o.x - q.x) < 30 && std::fabs(o.y - q.y) < 10) clash = true;
        if (clash) continue;
        placed.push_back(q);
        std::string lab = fmt_fixed(ticks[k], 0);
        double tw = text_width(lab, 5.5);
        f.rect(q.x - tw / 2 - 1, q.y - 4.5, tw + 2, 6, rgb(255, 255, 255), none());
        f.text({q.x, q.y}, lab, 5.5, Anchor::Middle, false, rgb(40, 40, 40));
        ++nlab;
      }
    }
  }
  // Fault polygons.
  for (std::size_t ff = 0; ff < cutoffs.size(); ++ff) {
    const repr::CutoffLine* hw = nullptr;
    const repr::CutoffLine* fw = nullptr;
    for (auto& l : cutoffs[ff]) (l.hanging_wall ? hw : fw) = &l;
    if (!hw || !fw) continue;
    // join consecutive columns where both walls exist into polygons
    std::vector<Pt> a, b;
    auto flush = [&] {
      if (a.size() >= 2) {
        std::vector<Pt> poly = a;
        poly.insert(poly.end(), b.rbegin(), b.rend());
        f.path(poly, true, rgb(235, 235, 235), kFault, 0.7);
      }
      a.clear();
      b.clear();
    };
    for (std::size_t i = 0; i < hw->by_column.size(); ++i) {
      const Vec3 &p = hw->by_column[i], &q = fw->by_column[i];
      if (std::isnan(p.x) || std::isnan(q.x)) { flush(); continue; }
      a.push_back(P(q.x, q.y));
      b.push_back(P(p.x, p.y));
    }
    flush();
  }
  // Closure.
  if (cl && cl->found) {
    for (auto& line : repr::contour_lines(cl->mask, 0.5)) {
      std::vector<Pt> pl;
      for (auto& p : line) pl.push_back(P(p.x, p.y));
      f.path(pl, false, none(), rgb(0, 140, 60), 1.6);
    }
    Pt sp = P(cl->spill_point.x, cl->spill_point.y);
    f.line({sp.x - 4, sp.y - 4}, {sp.x + 4, sp.y + 4}, rgb(0, 110, 50), 1.4);
    f.line({sp.x - 4, sp.y + 4}, {sp.x + 4, sp.y - 4}, rgb(0, 110, 50), 1.4);
    Pt cp = P(cl->crest.x, cl->crest.y);
    f.path({{cp.x, cp.y - 4}, {cp.x + 3.5, cp.y + 2.5}, {cp.x - 3.5, cp.y + 2.5}}, true, rgb(0, 140, 60), none());
  }
  // Stratal termination lines of the units above and below this horizon (stage 2C).
  std::vector<std::string> term_kinds;
  if (terms)
    for (auto& t : *terms) {
      std::vector<Pt> pl;
      for (auto& p : t.line) pl.push_back(P(p.x, p.y));
      f.path(pl, false, none(), termination_color(t.kind), 1.3, {4, 2});
      if (std::find(term_kinds.begin(), term_kinds.end(), t.kind) == term_kinds.end()) term_kinds.push_back(t.kind);
    }
  // Wells.
  for (auto& w : c.wells) {
    Pt p = P(w.x(), w.y());
    f.path({{p.x - 2.5, p.y}, {p.x, p.y - 2.5}, {p.x + 2.5, p.y}, {p.x, p.y + 2.5}}, true, rgb(255, 255, 255), kWell, 0.9);
    f.text({p.x + 4, p.y - 3}, w.name(), 7, Anchor::Start, true, kWell);
  }
  f.unclip();
  f.rect(mx0, my0, mw, mh, none(), kInk, 0.6);
  for (double t : nice_ticks(d.x0, d.x0 + d.sx, 6)) {
    Pt p = P(t, d.y0);
    f.line(p, {p.x, p.y + 3}, kInk, 0.5);
    f.text({p.x, p.y + 11}, km_or_m(t), 6.5, Anchor::Middle, false, kInk);
  }
  for (double t : nice_ticks(d.y0, d.y0 + d.sy, 6)) {
    Pt p = P(d.x0, t);
    f.line({p.x - 3, p.y}, p, kInk, 0.5);
    f.text({p.x - 5, p.y + 2.5}, km_or_m(t), 6.5, Anchor::End, false, kInk);
  }
  f.text({mx0 + mw / 2, my0 + mh + 24}, "Easting (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({mx0 - 40, my0 + mh / 2}, "Northing (m)", 7.5, Anchor::Middle, false, kInk, -90);
  // North arrow and scale bar.
  double nx = mx0 + mw - 18, ny = my0 + 14;
  f.path({{nx, ny - 8}, {nx + 5, ny + 6}, {nx, ny + 3}, {nx - 5, ny + 6}}, true, kInk, none());
  f.text({nx, ny + 16}, "N", 7, Anchor::Middle, true, kInk);
  auto sb = nice_ticks(0, d.sx / 4, 1);
  double sbl = sb.size() > 1 ? sb[1] : d.sx / 5;
  f.rect(mx0 + 10, my0 + mh - 16, sbl * sc, 4, kInk, kInk, 0.3);
  f.text({mx0 + 10 + sbl * sc / 2, my0 + mh - 20}, fmt_fixed(sbl, 0) + " m", 6.5, Anchor::Middle, false, kInk);
  // Colour bar and legend.
  double rx = 842 - 200, ry = 60;
  f.text({rx, ry}, "Depth (m)", 7.5, Anchor::Start, true, kInk);
  if (flat) {
    f.rect(rx, ry + 6, 14, 14, ramp(0.5), kInk, 0.4);
    f.text({rx + 20, ry + 16}, fmt_fixed(-zmax, 1) + " (flat: relief below 0.5 m, no contours)", 6.5, Anchor::Start, false, kInk);
  } else {
    for (int k = 0; k < 100; ++k) f.rect(rx, ry + 6 + k * 1.6, 14, 1.7, ramp(k / 99.0), none());
    f.rect(rx, ry + 6, 14, 160, none(), kInk, 0.4);
    for (double t : nice_ticks(-zmax, -zmin, 6)) {
      double yy = ry + 6 + (t + zmax) / std::max(1.0, zmax - zmin) * 160;
      f.line({rx + 14, yy}, {rx + 17, yy}, kInk, 0.4);
      f.text({rx + 20, yy + 2.5}, fmt_fixed(t, 0), 6.5, Anchor::Start, false, kInk);
    }
    f.text({rx, ry + 178}, "Contours every " + fmt_shortest(ci) + " m (labelled every " + fmt_shortest(2 * ci) + " m)", 6.5, Anchor::Start,
           false, kInk);
  }
  std::vector<LegendEntry> es;
  es.push_back({"Symbols", "Fault polygon (gap)", rgb(235, 235, 235)});
  es.push_back({"Symbols", "Fault cut-off lines", none(), true, kFault, 0.7});
  if (cl && cl->found) {
    es.push_back({"Symbols", "Closure outline", none(), true, rgb(0, 140, 60), 1.6});
    es.push_back({"Symbols", "Crest", rgb(0, 140, 60), false, none(), 1, 2});
    es.push_back({"Symbols", "Spill point", none(), false, rgb(0, 110, 50), 1, 3});
  }
  for (auto& k : term_kinds) es.push_back({"Symbols", termination_label(k), none(), true, termination_color(k), 1.3});
  if (!c.wells.empty()) es.push_back({"Symbols", "Well", rgb(255, 255, 255), false, none(), 1, 1});
  double ly = legend(f, rx, ry + 190, 190, es, 595.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, rx, ly, 190);
  footer(f, c);
  return f;
}

Figure unfolded_figure(const FigureContext& c, const std::string& section, const std::vector<repr::UnfoldedHorizon>& beds) {
  Figure f(842, 595);
  const auto& inf = c.model->info();
  title(f, (c.title.empty() ? c.model_name : c.title) + " - section " + section + " unfolded by flexural slip",
        "Each bed straightened along its own length from a pin at the section's start, fault gaps closed, at its depth at the pin "
        "(Dahlstrom 1969); bed lengths in the legend");
  const double px = 70, py = 70, pw = 842 - 300, ph = 595 - 140;
  double lmax = 1, zmin = 1e300, zmax = -1e300;
  for (auto& b : beds) {
    lmax = std::max(lmax, b.traced_length);
    for (auto& p : b.restored) { zmin = std::min(zmin, p.z); zmax = std::max(zmax, p.z); }
  }
  if (!(zmax > zmin)) { zmin -= 1; zmax += 1; }
  const double pad = 0.08 * (zmax - zmin);
  zmin -= pad;
  zmax += pad;
  auto X = [&](double s) { return px + s / lmax * pw; };
  auto Y = [&](double z) { return py + (zmax - z) / (zmax - zmin) * ph; };
  f.rect(px, py, pw, ph, none(), kInk, 0.6);
  for (double t : nice_ticks(0, lmax, 6)) {
    f.line({X(t), py + ph}, {X(t), py + ph + 3}, kInk, 0.5);
    f.text({X(t), py + ph + 12}, fmt_fixed(t, 0), 6.5, Anchor::Middle, false, kInk);
  }
  for (double t : nice_ticks(-zmax, -zmin, 6)) {
    f.line({px - 3, Y(-t)}, {px, Y(-t)}, kInk, 0.5);
    f.text({px - 5, Y(-t) + 2.5}, fmt_fixed(t, 0), 6.5, Anchor::End, false, kInk);
  }
  f.text({px + pw / 2, py + ph + 28}, "Restored bed length from the pin (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({px - 45, py + ph / 2}, "Depth at the pin (m)", 7.5, Anchor::Middle, false, kInk, -90);
  std::vector<LegendEntry> es;
  for (auto& b : beds) {
    if (b.restored.size() < 2) continue;
    Color col = unit_color(c, static_cast<std::size_t>(b.horizon) + 1);
    Color dark = rgb(static_cast<std::uint8_t>(col.r * 0.6), static_cast<std::uint8_t>(col.g * 0.6), static_cast<std::uint8_t>(col.b * 0.6));
    f.line({X(0), Y(b.restored.front().z)}, {X(b.traced_length), Y(b.restored.back().z)}, dark, 1.4);
    es.push_back({"Beds", inf.horizons[static_cast<std::size_t>(b.horizon)].name + ": " + fmt_fixed(b.traced_length, 0) + " m", none(), true, dark, 1.4});
  }
  double ly = legend(f, 842 - 200, py, 190, es, 595.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 200, ly, 190);
  footer(f, c);
  return f;
}

namespace {
Color viridis_like(double t);  // defined with the map figures below
Figure burial_plot(const FigureContext& c, const FigureContext::Burial& b, const FigureContext::Burial::Field* fd) {
  Figure f(842, 595);
  const auto& inf = c.model->info();
  const bool ages = !b.ages.empty();
  if (fd)
    title(f, (c.title.empty() ? c.model_name : c.title) + " - " + fd->label + " through burial at " + b.location,
          std::string("Burial curves of the horizons coloured by ") + fd->label + " (" + fd->unit + "), linear between horizons and steps; " +
              (ages ? "time in Ma" : "deposition steps (no horizon ages given)") + "; blank where it is not computed");
  else
    title(f, (c.title.empty() ? c.model_name : c.title) + " - burial history at " + b.location,
          std::string("Burial depth of each horizon below the depositional surface, from compaction of every unit by its law; ") +
              (ages ? "time in Ma" : "deposition steps (no horizon ages given)"));
  const double px = 70, py = 70, pw = 842 - 300, ph = 595 - 140;
  double tmin = 0, tmax = static_cast<double>(b.depth.size() - 1), dmax = 1;
  if (ages) { tmin = 1e300; tmax = -1e300; for (double a : b.ages) { tmin = std::min(tmin, a); tmax = std::max(tmax, a); } }
  for (auto& row : b.depth)
    for (int h : b.horizons)
      if (!std::isnan(row[static_cast<std::size_t>(h)])) dmax = std::max(dmax, row[static_cast<std::size_t>(h)]);
  auto ticks_d = nice_ticks(0, dmax, 6);
  dmax = std::max(dmax, ticks_d.back());
  auto X = [&](std::size_t k) {
    double t = ages ? b.ages[k] : static_cast<double>(k);
    double u = ages ? (tmax - t) / std::max(1e-9, tmax - tmin) : (t - tmin) / std::max(1e-9, tmax - tmin);  // old on the left
    return px + u * pw;
  };
  auto Y = [&](double dd) { return py + dd / dmax * ph; };
  f.rect(px, py, pw, ph, none(), kInk, 0.6);
  for (double t : ticks_d) {
    f.line({px, Y(t)}, {px + pw, Y(t)}, kGrid, 0.4);
    f.text({px - 5, Y(t) + 2.5}, fmt_fixed(t, 0), 6.5, Anchor::End, false, kInk);
  }
  for (std::size_t k = 0; k < b.depth.size(); ++k) {
    f.line({X(k), py + ph}, {X(k), py + ph + 3}, kInk, 0.5);
    f.text({X(k), py + ph + 12}, ages ? fmt_shortest(b.ages[k]) : b.step_labels[k], 6.5, Anchor::Middle, false, kInk);
  }
  f.text({px + pw / 2, py + ph + 28}, ages ? "Age (Ma)" : "Horizon at the surface (deposition step)", 7.5, Anchor::Middle, false, kInk);
  f.text({px - 45, py + ph / 2}, "Burial depth (m)", 7.5, Anchor::Middle, false, kInk, -90);
  double fmin = 1e300, fmax = -1e300;
  if (fd) {
    for (auto& row : fd->v)
      for (double v : row)
        if (std::isfinite(v)) { fmin = std::min(fmin, v); fmax = std::max(fmax, v); }
    if (!(fmax > fmin)) fmax = fmin + 1;
    // cells between consecutive horizons (in depth order at the step) and consecutive steps, split
    // into strips coloured by the value interpolated along the cell
    std::vector<int> hs = b.horizons;
    for (std::size_t k = 0; k + 1 < b.depth.size() && k + 1 < fd->v.size(); ++k)
      for (std::size_t i = 0; i + 1 < hs.size(); ++i) {
        const auto h0 = static_cast<std::size_t>(hs[i]), h1 = static_cast<std::size_t>(hs[i + 1]);
        const double d00 = b.depth[k][h0], d01 = b.depth[k][h1], d10 = b.depth[k + 1][h0], d11 = b.depth[k + 1][h1];
        const double v00 = fd->v[k][h0], v01 = fd->v[k][h1], v10 = fd->v[k + 1][h0], v11 = fd->v[k + 1][h1];
        if (!std::isfinite(d10) || !std::isfinite(d11) || !std::isfinite(v10) || !std::isfinite(v11)) continue;
        // a horizon not yet deposited at step k starts at the surface with the younger one's value
        const double a00 = std::isfinite(d00) ? d00 : 0.0, a01 = std::isfinite(d01) ? d01 : 0.0;
        const double w00 = std::isfinite(v00) && std::isfinite(d00) ? v00 : v10, w01 = std::isfinite(v01) && std::isfinite(d01) ? v01 : v11;
        const int ns = 6;
        for (int sI = 0; sI < ns; ++sI) {
          const double t0 = static_cast<double>(sI) / ns, t1 = static_cast<double>(sI + 1) / ns;
          auto D = [&](double tt, double u) { return (1 - tt) * ((1 - u) * a00 + u * a01) + tt * ((1 - u) * d10 + u * d11); };
          auto V = [&](double tt, double u) { return (1 - tt) * ((1 - u) * w00 + u * w01) + tt * ((1 - u) * v10 + u * v11); };
          auto XT = [&](double tt) { return X(k) + tt * (X(k + 1) - X(k)); };
          for (int uI = 0; uI < ns; ++uI) {
            const double u0 = static_cast<double>(uI) / ns, u1 = static_cast<double>(uI + 1) / ns;
            const double vm = V(0.5 * (t0 + t1), 0.5 * (u0 + u1));
            const Color col = viridis_like(std::clamp((vm - fmin) / (fmax - fmin), 0.0, 1.0));
            f.path({{XT(t0), Y(D(t0, u0))}, {XT(t1), Y(D(t1, u0))}, {XT(t1), Y(D(t1, u1))}, {XT(t0), Y(D(t0, u1))}}, true, col, col, 0.2);
          }
        }
      }
  }
  std::vector<LegendEntry> es;
  for (int h : b.horizons) {
    std::vector<Pt> pl;
    for (std::size_t k = 0; k < b.depth.size(); ++k) {
      double v = b.depth[k][static_cast<std::size_t>(h)];
      if (!std::isnan(v)) pl.push_back({X(k), Y(v)});
    }
    Color col = unit_color(c, static_cast<std::size_t>(h) + 1);
    Color dark = rgb(static_cast<std::uint8_t>(col.r * 0.6), static_cast<std::uint8_t>(col.g * 0.6), static_cast<std::uint8_t>(col.b * 0.6));
    if (pl.size() >= 2) f.path(pl, false, none(), dark, 1.3);
    for (auto& p : pl) f.rect(p.x - 1.5, p.y - 1.5, 3, 3, dark, none());
    es.push_back({"Horizons", inf.horizons[static_cast<std::size_t>(h)].name, none(), true, dark, 1.3});
  }
  // stage 4B: tectonic subsidence from backstripping, on the same depth axis
  for (int kind = 0; kind < 2; ++kind) {
    const auto& v = kind == 0 ? b.tectonic_airy : b.tectonic_flexural;
    if (v.size() != b.depth.size()) continue;
    std::vector<Pt> pl;
    for (std::size_t k = 0; k < v.size(); ++k)
      if (std::isfinite(v[k]) && v[k] <= dmax) pl.push_back({X(k), Y(std::max(0.0, v[k]))});
    if (pl.size() < 2) continue;
    f.path(pl, false, none(), kInk, 1.6, kind == 0 ? std::vector<double>{5, 2} : std::vector<double>{1.5, 1.5});
    es.push_back({"Backstripping", kind == 0 ? "Tectonic subsidence, Airy (dashed)" : "Tectonic subsidence, flexural (dotted)", none(), true, kInk, 1.6});
  }
  double ly = py;
  if (fd) {
    const double rx = 842 - 200;
    f.text({rx, ly}, fd->label + " (" + fd->unit + ")", 7.5, Anchor::Start, true, kInk);
    for (int k = 0; k < 100; ++k) f.rect(rx, ly + 6 + k * 1.2, 14, 1.3, viridis_like(1 - k / 99.0), none());
    f.rect(rx, ly + 6, 14, 120, none(), kInk, 0.4);
    const int dec = std::clamp(2 - static_cast<int>(std::floor(std::log10(std::max(1e-12, std::fabs(fmax - fmin))))), 0, 6);
    f.text({rx + 18, ly + 12}, fmt_fixed(fmax, dec), 6.5, Anchor::Start, false, kInk);
    f.text({rx + 18, ly + 68}, fmt_fixed(0.5 * (fmin + fmax), dec), 6.5, Anchor::Start, false, kInk);
    f.text({rx + 18, ly + 126}, fmt_fixed(fmin, dec), 6.5, Anchor::Start, false, kInk);
    ly += 140;
  }
  ly = legend(f, 842 - 200, ly, 190, es, 595.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 200, ly, 190);
  footer(f, c);
  return f;
}
}  // namespace

Figure burial_figure(const FigureContext& c, const FigureContext::Burial& b) { return burial_plot(c, b, nullptr); }
Figure burial_field_figure(const FigureContext& c, const FigureContext::Burial& b, const FigureContext::Burial::Field& fd) { return burial_plot(c, b, &fd); }

Figure isopach_figure(const FigureContext& c, int unit, const Grid2D& g, const std::vector<repr::TerminationLine>& terms) {
  Figure f(842, 595);
  const auto& inf = c.model->info();
  const auto& d = inf.domain;
  double side = std::min(595.0 - 110, 842.0 - 280);
  double sc = side / std::max(d.sx, d.sy);
  double mx0 = 60, my0 = 60, mw = d.sx * sc, mh = d.sy * sc;
  auto P = [&](double x, double y) { return Pt{mx0 + (x - d.x0) * sc, my0 + mh - (y - d.y0) * sc}; };
  double tmax = 0;
  for (double v : g.v)
    if (!std::isnan(v)) tmax = std::max(tmax, v);
  const auto& un = inf.units[static_cast<std::size_t>(unit)];
  title(f, (c.title.empty() ? c.model_name : c.title) + " - " + un.name + " isopach",
        "Vertical thickness of " + un.name + " (" + inf.horizons[static_cast<std::size_t>(unit - 1)].name + " to " +
            inf.horizons[static_cast<std::size_t>(unit)].name + ") in each fault block; fault gaps grey" +
            (un.termination.empty() ? "" : "; terminates by " + un.termination));
  Image im;
  im.w = std::max(1, static_cast<int>(std::lround(mw * 3)));
  im.h = std::max(1, static_cast<int>(std::lround(mh * 3)));
  im.rgba.resize(static_cast<std::size_t>(im.w) * im.h * 4);
  const auto& gs = g.spec;
  for (int yy = 0; yy < im.h; ++yy)
    for (int xx = 0; xx < im.w; ++xx) {
      double x = d.x0 + (xx + 0.5) / im.w * d.sx, y = d.y0 + d.sy - (yy + 0.5) / im.h * d.sy;
      long i = std::lround((x - gs.x0) / gs.dx), j = std::lround((y - gs.y0) / gs.dy);
      i = std::clamp(i, 0L, static_cast<long>(gs.nx) - 1);
      j = std::clamp(j, 0L, static_cast<long>(gs.ny) - 1);
      double v = g.v[gs.index(static_cast<std::size_t>(i), static_cast<std::size_t>(j))];
      std::uint8_t* p = &im.rgba[(static_cast<std::size_t>(yy) * im.w + xx) * 4];
      Color cc = std::isnan(v) ? rgb(235, 235, 235) : v < 0.01 ? rgb(250, 250, 250) : ramp(1 - v / std::max(1.0, tmax));
      p[0] = cc.r; p[1] = cc.g; p[2] = cc.b; p[3] = 255;
    }
  const Image* img = f.own(std::move(im));
  f.clip(mx0, my0, mw, mh);
  f.image(mx0, my0, mw, mh, img);
  auto ticks = nice_ticks(0, std::max(1.0, tmax), 10);
  for (double t : ticks) {
    if (t <= 0) continue;
    for (auto& line : repr::contour_lines(g, t)) {
      if (line.size() < 3) continue;
      std::vector<Pt> pl;
      for (auto& p : line) pl.push_back(P(p.x, p.y));
      f.path(pl, false, none(), rgb(90, 90, 90), 0.4);
    }
  }
  std::vector<std::string> kinds;
  for (auto& t : terms) {
    std::vector<Pt> pl;
    for (auto& p : t.line) pl.push_back(P(p.x, p.y));
    f.path(pl, false, none(), termination_color(t.kind), 1.4, {4, 2});
    if (std::find(kinds.begin(), kinds.end(), t.kind) == kinds.end()) kinds.push_back(t.kind);
  }
  for (auto& w : c.wells) {
    Pt p = P(w.x(), w.y());
    f.path({{p.x - 2.5, p.y}, {p.x, p.y - 2.5}, {p.x + 2.5, p.y}, {p.x, p.y + 2.5}}, true, rgb(255, 255, 255), kWell, 0.9);
    f.text({p.x + 4, p.y - 3}, w.name(), 7, Anchor::Start, true, kWell);
  }
  f.unclip();
  f.rect(mx0, my0, mw, mh, none(), kInk, 0.6);
  for (double t : nice_ticks(d.x0, d.x0 + d.sx, 6)) {
    Pt p = P(t, d.y0);
    f.line(p, {p.x, p.y + 3}, kInk, 0.5);
    f.text({p.x, p.y + 11}, km_or_m(t), 6.5, Anchor::Middle, false, kInk);
  }
  for (double t : nice_ticks(d.y0, d.y0 + d.sy, 6)) {
    Pt p = P(d.x0, t);
    f.line({p.x - 3, p.y}, p, kInk, 0.5);
    f.text({p.x - 5, p.y + 2.5}, km_or_m(t), 6.5, Anchor::End, false, kInk);
  }
  f.text({mx0 + mw / 2, my0 + mh + 24}, "Easting (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({mx0 - 40, my0 + mh / 2}, "Northing (m)", 7.5, Anchor::Middle, false, kInk, -90);
  double rx = 842 - 200, ry = 60;
  f.text({rx, ry}, "Thickness (m)", 7.5, Anchor::Start, true, kInk);
  for (int k = 0; k < 100; ++k) f.rect(rx, ry + 6 + k * 1.6, 14, 1.7, ramp(k / 99.0), none());
  f.rect(rx, ry + 6, 14, 160, none(), kInk, 0.4);
  for (double t : nice_ticks(0, std::max(1.0, tmax), 6)) {
    double yy = ry + 6 + (1 - t / std::max(1.0, tmax)) * 160;
    if (yy < ry + 6 - 0.1 || yy > ry + 166.1) continue;
    f.line({rx + 14, yy}, {rx + 17, yy}, kInk, 0.4);
    f.text({rx + 20, yy + 2.5}, fmt_fixed(t, 0), 6.5, Anchor::Start, false, kInk);
  }
  f.text({rx, ry + 178}, "Contours every " + fmt_shortest(ticks.size() > 1 ? ticks[1] - ticks[0] : tmax) + " m; white: absent",
         6.5, Anchor::Start, false, kInk);
  std::vector<LegendEntry> es;
  es.push_back({"Symbols", "Fault gap", rgb(235, 235, 235)});
  for (auto& k : kinds) es.push_back({"Symbols", termination_label(k), none(), true, termination_color(k), 1.4});
  if (!c.wells.empty()) es.push_back({"Symbols", "Well", rgb(255, 255, 255), false, none(), 1, 1});
  double ly = legend(f, rx, ry + 190, 190, es, 595.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, rx, ly, 190);
  footer(f, c);
  return f;
}

// =========================================================================================
// 3D view
// =========================================================================================
Figure view3d_figure(const FigureContext& c, double az, double el, const std::string& preset) {
  Figure f(842, 560);
  const auto& inf = c.model->info();
  const Domain& d = inf.domain;
  double sp = std::max(d.sx, d.sy) / 150;
  GridSpec g = GridSpec::covering(d, sp);
  std::vector<Mesh3> meshes;
  // Faults and the deepest horizon opaque; shallower horizons see-through, deepest first.
  for (std::size_t ff = 0; ff < inf.faults.size(); ++ff) {
    auto fs = repr::fault_surface(*c.model, static_cast<int>(ff), sp);
    meshes.push_back({fs.mesh, rgb(200, 50, 40), 1.0});
  }
  // With many horizons the see-through layers are made fainter so that at least 20 % of the light
  // from the faults and the deepest horizon passes through all of them.
  const double n_see = std::max<double>(1.0, static_cast<double>(inf.horizons.size()) - 1.0);
  const double alpha = std::min(0.4, 1.0 - std::pow(0.2, 1.0 / n_see));
  for (std::size_t k = inf.horizons.size(); k-- > 0;) {
    repr::FaultedSurface scratch;
    const auto& fs = repr::faulted_surface(*c.model, static_cast<int>(k), g, c.surfaces, scratch);
    meshes.push_back({repr::grid_to_mesh(fs.z), unit_color(c, k + 1), k + 1 == inf.horizons.size() ? 1.0 : alpha});
  }
  BBox box = d.box();
  Image im = render_meshes(meshes, az, el, c.vertical_exaggeration, 1800, 1000, box);
  title(f, (c.title.empty() ? c.model_name : c.title) + " - 3D view (" + preset + ")",
        "Horizons coloured by the unit below them, faults in red; viewed from azimuth " + fmt_fixed(az, 0) + " deg, elevation " +
            fmt_fixed(el, 0) + " deg; vertical exaggeration " + fmt_fixed(c.vertical_exaggeration, 1) + "x");
  const Image* img = f.own(std::move(im));
  double iw = 842 - 250, ih = iw * 1000.0 / 1800.0;
  f.image(20, 70, iw, ih, img);
  f.text({20 + iw / 2, 70 + ih + 14}, "Shallower horizons drawn see-through so faults and deeper horizons remain visible.", 7, Anchor::Middle, false,
         rgb(70, 70, 70));
  std::vector<LegendEntry> es;
  for (std::size_t h = 0; h < inf.horizons.size(); ++h) es.push_back({"Horizons (top of unit)", inf.horizons[h].name, unit_color(c, h + 1)});
  es.push_back({"Faults", "Fault surfaces (clipped to tips)", rgb(200, 50, 40)});
  double ly = legend(f, 842 - 214, 56, 196, es, 560.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 214, ly, 196);
  footer(f, c);
  return f;
}

// =========================================================================================
// Fence diagram
// =========================================================================================
Figure fence_figure(const FigureContext& c, const std::vector<repr::Section>& legs) {
  Figure f(842, 560);
  const auto& inf = c.model->info();
  const Domain& d = inf.domain;
  const double az = deg2rad(210), el = deg2rad(28), ve = c.vertical_exaggeration;
  Vec3 ctr{d.x0 + d.sx / 2, d.y0 + d.sy / 2, 0.5 * (d.z_top + d.z_base)};
  Vec3 fwd{-std::sin(az) * std::cos(el), -std::cos(az) * std::cos(el), -std::sin(el)};
  Vec3 right = normalized(cross(fwd, {0, 0, 1}));
  Vec3 up = cross(right, fwd);
  auto proj = [&](double x, double y, double z, double& u, double& v, double& dep) {
    Vec3 p{x - ctr.x, y - ctr.y, (z - ctr.z) * ve};
    u = dot(p, right);
    v = dot(p, up);
    dep = dot(p, fwd);
  };
  // extents of all legs
  double umin = 1e30, umax = -1e30, vmin = 1e30, vmax = -1e30;
  double zt = -1e30, zb = 1e30;
  for (auto& s : legs)
    for (auto& hz : s.horizon_z)
      for (double z : hz)
        if (!std::isnan(z)) { zt = std::max(zt, z); zb = std::min(zb, z); }
  double span = std::max(100.0, zt - zb);
  zt = std::min(d.z_top, zt + 0.3 * span);
  zb = std::max(d.z_base, zb - 0.3 * span);
  for (auto& s : legs)
    for (auto& p : s.xy)
      for (double z : {zt, zb}) {
        double u, v, dp;
        proj(p.x, p.y, z, u, v, dp);
        umin = std::min(umin, u); umax = std::max(umax, u); vmin = std::min(vmin, v); vmax = std::max(vmax, v);
      }
  double x0 = 60, y0 = 78, W = 842 - 290, H = 560 - 130;
  double sc = std::min(W / (umax - umin), H / (vmax - vmin));
  auto P = [&](double x, double y, double z) {
    double u, v, dp;
    proj(x, y, z, u, v, dp);
    return Pt{x0 + (u - umin) * sc + (W - (umax - umin) * sc) / 2, y0 + (vmax - v) * sc};
  };
  title(f, (c.title.empty() ? c.model_name : c.title) + " - fence diagram through wells",
        "Panels between wells, viewed from azimuth 210 deg, elevation 28 deg; vertical exaggeration " + fmt_fixed(ve, 1) + "x");
  // Painter's order: farthest legs first.
  std::vector<std::pair<double, std::size_t>> order;
  for (std::size_t k = 0; k < legs.size(); ++k) {
    auto& s = legs[k];
    double u, v, dp;
    proj(s.xy[s.ns / 2].x, s.xy[s.ns / 2].y, 0.5 * (zt + zb), u, v, dp);
    order.push_back({-dp, k});
  }
  std::sort(order.begin(), order.end());
  const std::size_t nh = inf.horizons.size();
  for (auto& [dd, k] : order) {
    auto& s = legs[k];
    for (std::size_t col = 0; col + 1 < s.ns; ++col) {
      auto a = s.xy[col], b = s.xy[col + 1];
      std::size_t r = 0;
      while (r < s.nz) {
        std::int32_t u = s.unit[r * s.ns + col];
        std::size_t r2 = r;
        while (r2 + 1 < s.nz && s.unit[(r2 + 1) * s.ns + col] == u) ++r2;
        double ztop = std::min(zt, s.z_top - (static_cast<double>(r) - 0.5) * s.dz), zbot = std::max(zb, s.z_top - (static_cast<double>(r2) + 0.5) * s.dz);
        if (ztop > zbot) {
          Color cc = unit_color(c, static_cast<std::size_t>(u));
          f.path({P(a.x, a.y, ztop), P(b.x, b.y, ztop), P(b.x, b.y, zbot), P(a.x, a.y, zbot)}, true, cc, cc, 0.25);
        }
        r = r2 + 1;
      }
    }
    for (std::size_t h = 0; h < nh; ++h) {
      std::vector<Pt> cur;
      for (std::size_t col = 0; col < s.ns; ++col) {
        double z = s.horizon_z[h][col];
        if (std::isnan(z) || z > zt || z < zb) { if (cur.size() > 1) f.path(cur, false, none(), kInk, 0.5); cur.clear(); continue; }
        cur.push_back(P(s.xy[col].x, s.xy[col].y, z));
      }
      if (cur.size() > 1) f.path(cur, false, none(), kInk, 0.5);
    }
    for (std::size_t ff = 0; ff < s.fault_z.size(); ++ff) {
      std::vector<Pt> cur;
      for (std::size_t col = 0; col < s.ns; ++col) {
        double z = s.fault_z[ff][col];
        if (std::isnan(z) || z > zt || z < zb) { if (cur.size() > 1) f.path(cur, false, none(), kFault, 1.2); cur.clear(); continue; }
        cur.push_back(P(s.xy[col].x, s.xy[col].y, z));
      }
      if (cur.size() > 1) f.path(cur, false, none(), kFault, 1.2);
    }
    // panel outline
    f.path({P(s.xy.front().x, s.xy.front().y, zt), P(s.xy.back().x, s.xy.back().y, zt), P(s.xy.back().x, s.xy.back().y, zb),
            P(s.xy.front().x, s.xy.front().y, zb)}, true, none(), rgb(90, 90, 90), 0.5);
  }
  for (auto& w : c.wells) {
    Pt a = P(w.x(), w.y(), zt), b = P(w.x(), w.y(), zb);
    f.line(a, b, kWell, 1.0);
    f.text({a.x, a.y - 4}, w.name(), 7.5, Anchor::Middle, true, kWell);
    for (auto& t : repr::well_intersections(*c.model, w, 2.0)) {
      if (t.p.z > zt || t.p.z < zb) continue;
      Pt q = P(w.x(), w.y(), t.p.z);
      f.line({q.x - 2.5, q.y}, {q.x + 2.5, q.y}, kWell, 1.4);
    }
  }
  // Depth scale along the first well.
  if (!c.wells.empty()) {
    const auto& w = c.wells.front();
    for (double t : nice_ticks(-zt, -zb, 6)) {
      Pt q = P(w.x(), w.y(), -t);
      f.line({q.x - 5, q.y}, q, kInk, 0.5);
      f.text({q.x - 7, q.y + 2.5}, fmt_fixed(t, 0), 6.5, Anchor::End, false, kInk);
    }
    Pt top = P(w.x(), w.y(), zt);
    f.text({top.x - 12, top.y - 18}, "Depth along " + w.name() + " (m)", 6.5, Anchor::End, false, kInk);
  }
  auto es = unit_legend(c);
  es.push_back({"Lines", "Horizon", none(), true, kInk, 0.5});
  es.push_back({"Lines", "Fault (between tips)", none(), true, kFault, 1.2});
  es.push_back({"Lines", "Well and tops", none(), true, kWell, 1.0});
  double ly = legend(f, 842 - 214, 56, 196, es, 560.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 214, ly, 196);
  footer(f, c);
  return f;
}

// =========================================================================================
// Writing
// =========================================================================================
// =========================================================================================
// Fault systems (stage 2A): throw along strike, and the displacement on a fault surface
// =========================================================================================
Figure strike_throw_figure(const FigureContext& c, const std::vector<StrikeThrowSeries>& series, const std::string& system_name) {
  Figure f(842, 520);
  double x0 = 1e300, x1 = -1e300, y1 = 0;
  for (auto& s : series)
    for (auto& [t, v] : s.points) { x0 = std::min(x0, t); x1 = std::max(x1, t); if (!std::isnan(v)) y1 = std::max(y1, v); }
  if (!(x1 > x0)) { x0 = 0; x1 = 1; }
  if (!(y1 > 0)) y1 = 1;
  y1 *= 1.1;
  title(f, (c.title.empty() ? c.model_name : c.title) + " - throw along strike, " + system_name,
        "Cut-off separation per horizon along the first fault's strike; solid: summed over the system (open where throw passes "
        "between faults), dashed: each fault");
  const double px = 70, py = 60, pw = 842 - 70 - 250, ph = 520 - 60 - 70;
  auto X = [&](double t) { return px + (t - x0) / (x1 - x0) * pw; };
  auto Y = [&](double v) { return py + ph - v / y1 * ph; };
  for (double t : nice_ticks(x0, x1, 8)) {
    f.line({X(t), py}, {X(t), py + ph}, kGrid, 0.4);
    f.text({X(t), py + ph + 11}, km_or_m(t), 7, Anchor::Middle, false, kInk);
  }
  for (double v : nice_ticks(0, y1, 6)) {
    f.line({px, Y(v)}, {px + pw, Y(v)}, kGrid, 0.4);
    f.text({px - 4, Y(v) + 2.5}, fmt_fixed(v, 0), 7, Anchor::End, false, kInk);
  }
  f.rect(px, py, pw, ph, none(), kInk, 0.6);
  f.text({px + pw / 2, py + ph + 26}, "Distance along strike (m)", 8, Anchor::Middle, false, kInk);
  f.text({px - 40, py + ph / 2}, "Throw (m)", 8, Anchor::Middle, false, kInk, -90);
  std::vector<LegendEntry> es;
  for (auto& s : series) {
    std::vector<Pt> pts;
    auto flush = [&] {
      if (pts.size() >= 2) f.path(pts, false, none(), s.color, s.summed ? 1.4 : 0.9, s.summed ? std::vector<double>{} : std::vector<double>{3, 2});
      pts.clear();
    };
    for (auto& [t, v] : s.points) {
      if (std::isnan(v)) { flush(); continue; }
      pts.push_back({X(t), Y(v)});
    }
    flush();
    es.push_back({s.summed ? "Summed throw" : "Each fault", s.label, none(), true, s.color, s.summed ? 1.4 : 0.9});
  }
  double ly = legend(f, 842 - 214, 56, 196, es, 520.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 214, ly, 196);
  footer(f, c);
  return f;
}

Figure fault_plane_figure(const FigureContext& c, int ff, const repr::FaultSurface& fs, const std::vector<double>& disp,
                          const std::vector<std::vector<repr::CutoffLine>>& cutoffs) {
  Figure f(842, 560);
  const auto& inf = c.model->info();
  const auto& fi = inf.faults.at(static_cast<std::size_t>(ff));
  double a0 = fs.a.front(), a1 = fs.a.back(), b0 = fs.b.front(), b1 = fs.b.back();
  // only the part inside the tip line and the domain
  double ia0 = 1e300, ia1 = -1e300, ib0 = 1e300, ib1 = -1e300, dmax = 0;
  for (std::size_t j = 0; j < fs.nb; ++j)
    for (std::size_t i = 0; i < fs.na; ++i) {
      std::size_t k = i + fs.na * j;
      if (fs.extent[k] > 1) continue;
      ia0 = std::min(ia0, fs.a[i]); ia1 = std::max(ia1, fs.a[i]);
      ib0 = std::min(ib0, fs.b[j]); ib1 = std::max(ib1, fs.b[j]);
      dmax = std::max(dmax, disp[k]);
    }
  if (ia1 > ia0) { a0 = ia0; a1 = ia1; }
  if (ib1 > ib0) { b0 = ib0; b1 = ib1; }
  if (!(dmax > 0)) dmax = 1;
  title(f, (c.title.empty() ? c.model_name : c.title) + " - displacement on fault " + fi.name,
        "Fault surface seen from the hanging wall, along strike and down dip; colour: displacement magnitude (m); "
        "lines: hanging-wall (solid) and footwall (dashed) cut-offs of each horizon");
  const double px = 70, py = 60, maxw = 842 - 70 - 260, maxh = 560 - 60 - 80;
  double sc = std::min(maxw / (a1 - a0), maxh / (b1 - b0));
  const double pw = (a1 - a0) * sc, ph = (b1 - b0) * sc;
  auto P = [&](double a, double b) { return Pt{px + (a - a0) * sc, py + (b - b0) * sc}; };
  auto colour = [](double t) {  // light (small) to dark red (large)
    t = std::clamp(t, 0.0, 1.0);
    return rgb(static_cast<int>(255 - 75 * t), static_cast<int>(245 - 205 * t), static_cast<int>(220 - 190 * t));
  };
  for (std::size_t j = 0; j + 1 < fs.nb; ++j)
    for (std::size_t i = 0; i + 1 < fs.na; ++i) {
      std::size_t k = i + fs.na * j;
      if (fs.extent[k] > 1 || fs.extent[k + 1] > 1 || fs.extent[k + fs.na] > 1 || fs.extent[k + fs.na + 1] > 1) continue;
      double v = 0.25 * (disp[k] + disp[k + 1] + disp[k + fs.na] + disp[k + fs.na + 1]);
      Pt q = P(fs.a[i], fs.b[j]), r = P(fs.a[i + 1], fs.b[j + 1]);
      f.rect(q.x, q.y, r.x - q.x + 0.3, r.y - q.y + 0.3, colour(v / dmax), none());
    }
  f.rect(px, py, pw, ph, none(), kInk, 0.6);
  for (double t : nice_ticks(a0, a1, 7)) f.text({P(t, b0).x, py + ph + 11}, km_or_m(t), 7, Anchor::Middle, false, kInk);
  for (double t : nice_ticks(b0, b1, 6)) f.text({px - 4, P(a0, t).y + 2.5}, km_or_m(t), 7, Anchor::End, false, kInk);
  f.text({px + pw / 2, py + ph + 26}, "Along strike from the reference point (m)", 8, Anchor::Middle, false, kInk);
  f.text({px - 44, py + ph / 2}, "Down dip (m)", 8, Anchor::Middle, false, kInk, -90);
  std::vector<LegendEntry> es;
  for (std::size_t h = 0; h < cutoffs.size(); ++h)
    for (auto& cl : cutoffs[h]) {
      for (auto& part : cl.parts) {
        std::vector<Pt> pts;
        for (auto& p : part) pts.push_back(P(dot(p - fi.ref_point, fi.strike_dir), dot(p - fi.ref_point, fi.dip_dir)));
        if (pts.size() >= 2)
          f.path(pts, false, none(), rgb(40, 40, 40), 0.8,
                 cl.hanging_wall ? std::vector<double>{} : std::vector<double>{3, 2});
      }
    }
  // colour bar
  double cbx = 842 - 214, cby = 60;
  f.text({cbx, cby + 8}, "Displacement (m)", 7.5, Anchor::Start, true, kInk);
  for (int k = 0; k < 100; ++k) f.rect(cbx + 1.6 * k, cby + 14, 1.7, 10, colour(k / 99.0), none());
  f.rect(cbx, cby + 14, 160, 10, none(), kInk, 0.4);
  f.text({cbx, cby + 34}, "0", 7, Anchor::Start, false, kInk);
  f.text({cbx + 160, cby + 34}, fmt_fixed(dmax, 1), 7, Anchor::End, false, kInk);
  es.push_back({"Cut-offs", "Hanging wall", none(), true, rgb(40, 40, 40), 0.8});
  es.push_back({"Cut-offs", "Footwall (dashed)", none(), true, rgb(40, 40, 40), 0.8});
  for (auto& r : fi.relations) es.push_back({"Relations", r, none(), true, kFault, 0.8});
  double ly = legend(f, cbx, cby + 44, 196, es, 560.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, cbx, ly, 196);
  footer(f, c);
  return f;
}

namespace {
/// Layout of a fault surface seen from the hanging wall (the part inside the tip line).
struct FaultFrame {
  double a0, a1, b0, b1, px = 70, py = 60, sc = 1, pw = 0, ph = 0;
  Pt P(double a, double b) const { return Pt{px + (a - a0) * sc, py + (b - b0) * sc}; }
};
FaultFrame fault_frame(const repr::FaultSurface& fs) {
  FaultFrame fr{fs.a.front(), fs.a.back(), fs.b.front(), fs.b.back()};
  double ia0 = 1e300, ia1 = -1e300, ib0 = 1e300, ib1 = -1e300;
  for (std::size_t j = 0; j < fs.nb; ++j)
    for (std::size_t i = 0; i < fs.na; ++i) {
      if (fs.extent[i + fs.na * j] > 1) continue;
      ia0 = std::min(ia0, fs.a[i]); ia1 = std::max(ia1, fs.a[i]);
      ib0 = std::min(ib0, fs.b[j]); ib1 = std::max(ib1, fs.b[j]);
    }
  if (ia1 > ia0) { fr.a0 = ia0; fr.a1 = ia1; }
  if (ib1 > ib0) { fr.b0 = ib0; fr.b1 = ib1; }
  const double maxw = 842 - 70 - 260, maxh = 560 - 60 - 80;
  fr.sc = std::min(maxw / (fr.a1 - fr.a0), maxh / (fr.b1 - fr.b0));
  fr.pw = (fr.a1 - fr.a0) * fr.sc;
  fr.ph = (fr.b1 - fr.b0) * fr.sc;
  return fr;
}
void fault_frame_axes_and_cutoffs(Figure& f, const FaultFrame& fr, const FaultInfo& fi, const std::vector<std::vector<repr::CutoffLine>>& cutoffs) {
  for (auto& lines : cutoffs)
    for (auto& cl : lines)
      for (auto& part : cl.parts) {
        std::vector<Pt> pts;
        for (auto& p : part) pts.push_back(fr.P(dot(p - fi.ref_point, fi.strike_dir), dot(p - fi.ref_point, fi.dip_dir)));
        if (pts.size() >= 2) f.path(pts, false, none(), rgb(40, 40, 40), 0.8, cl.hanging_wall ? std::vector<double>{} : std::vector<double>{3, 2});
      }
  f.rect(fr.px, fr.py, fr.pw, fr.ph, none(), kInk, 0.6);
  for (double t : nice_ticks(fr.a0, fr.a1, 7)) f.text({fr.P(t, fr.b0).x, fr.py + fr.ph + 11}, km_or_m(t), 7, Anchor::Middle, false, kInk);
  for (double t : nice_ticks(fr.b0, fr.b1, 6)) f.text({fr.px - 4, fr.P(fr.a0, t).y + 2.5}, km_or_m(t), 7, Anchor::End, false, kInk);
  f.text({fr.px + fr.pw / 2, fr.py + fr.ph + 26}, "Along strike from the reference point (m)", 8, Anchor::Middle, false, kInk);
  f.text({fr.px - 44, fr.py + fr.ph / 2}, "Down dip (m)", 8, Anchor::Middle, false, kInk, -90);
}
Color viridis_like(double t) {  // dark blue (0) through green to yellow (1)
  t = std::clamp(t, 0.0, 1.0);
  const double r = std::clamp(-0.1 + 1.2 * t * t, 0.0, 1.0), g = std::clamp(0.1 + 0.85 * t, 0.0, 1.0), b = std::clamp(0.45 + 0.5 * std::sin(kPi * (0.15 + 0.85 * t)) - 0.4 * t, 0.0, 1.0);
  return rgb(static_cast<int>(40 + 215 * r), static_cast<int>(20 + 215 * g), static_cast<int>(60 + 140 * b));
}
Color diverging(double t) {  // -1 blue, 0 white, +1 red
  t = std::clamp(t, -1.0, 1.0);
  if (t >= 0) return rgb(255, static_cast<int>(255 - 190 * t), static_cast<int>(255 - 210 * t));
  return rgb(static_cast<int>(255 + 210 * t), static_cast<int>(255 + 150 * t), 255);
}
}  // namespace

Figure fault_attribute_figure(const FigureContext& c, int ff, const repr::FaultSurface& fs, const std::vector<double>& v,
                              const std::vector<std::vector<repr::CutoffLine>>& cutoffs, const std::string& attr, const std::string& unit,
                              double vmin, double vmax) {
  Figure f(842, 560);
  const auto& fi = c.model->info().faults.at(static_cast<std::size_t>(ff));
  auto fr = fault_frame(fs);
  double lo = 1e300, hi = -1e300;
  for (std::size_t k = 0; k < v.size(); ++k)
    if (fs.extent[k] <= 1 && std::isfinite(v[k])) { lo = std::min(lo, v[k]); hi = std::max(hi, v[k]); }
  if (!(hi >= lo)) { lo = 0; hi = 1; }
  const bool uniform = hi - lo <= 1e-6 * std::max(1.0, std::fabs(hi));
  const double uval = hi;
  if (!std::isnan(vmin)) lo = vmin;
  if (!std::isnan(vmax)) hi = vmax;
  if (!(hi > lo)) hi = lo + 1;
  title(f, (c.title.empty() ? c.model_name : c.title) + " - " + attr + " on fault " + fi.name,
        "Fault surface seen from the hanging wall; colour: " + attr + " (" + unit + ")" +
            (uniform ? ", uniform at " + fmt_fixed(uval, 3) + " (planar fault, stresses proportional to depth)" : std::string()) +
            ", grey where not defined; lines: hanging-wall (solid) and footwall (dashed) cut-offs");
  for (std::size_t j = 0; j + 1 < fs.nb; ++j)
    for (std::size_t i = 0; i + 1 < fs.na; ++i) {
      std::size_t k = i + fs.na * j;
      if (fs.extent[k] > 1 || fs.extent[k + 1] > 1 || fs.extent[k + fs.na] > 1 || fs.extent[k + fs.na + 1] > 1) continue;
      double s = 0;
      int n = 0;
      for (std::size_t q : {k, k + 1, k + fs.na, k + fs.na + 1})
        if (std::isfinite(v[q])) { s += v[q]; ++n; }
      Pt q0 = fr.P(fs.a[i], fs.b[j]), r = fr.P(fs.a[i + 1], fs.b[j + 1]);
      f.rect(q0.x, q0.y, r.x - q0.x + 0.3, r.y - q0.y + 0.3, n ? viridis_like((s / n - lo) / (hi - lo)) : rgb(215, 215, 215), none());
    }
  fault_frame_axes_and_cutoffs(f, fr, fi, cutoffs);
  double cbx = 842 - 214, cby = 60;
  f.text({cbx, cby + 8}, attr + " (" + unit + ")", 7.5, Anchor::Start, true, kInk);
  for (int k = 0; k < 100; ++k) f.rect(cbx + 1.6 * k, cby + 14, 1.7, 10, viridis_like(k / 99.0), none());
  f.rect(cbx, cby + 14, 160, 10, none(), kInk, 0.4);
  auto num = [](double x) { return std::fabs(x) >= 100 || x == 0 ? fmt_fixed(x, 0) : std::fabs(x) >= 1 ? fmt_fixed(x, 2) : fmt_fixed(x, 3); };
  f.text({cbx, cby + 34}, num(lo), 7, Anchor::Start, false, kInk);
  f.text({cbx + 160, cby + 34}, num(hi), 7, Anchor::End, false, kInk);
  std::vector<LegendEntry> es{{"Cut-offs", "Hanging wall", none(), true, rgb(40, 40, 40), 0.8},
                              {"Cut-offs", "Footwall (dashed)", none(), true, rgb(40, 40, 40), 0.8},
                              {"Not defined", "Missing input or no separation", rgb(215, 215, 215)}};
  double ly = legend(f, cbx, cby + 44, 196, es, 560.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, cbx, ly, 196);
  footer(f, c);
  return f;
}

Figure allan_figure(const FigureContext& c, int ff, const repr::FaultSurface& fs, const std::vector<repr::FaultRockPoint>& pts,
                    const std::vector<std::vector<repr::CutoffLine>>& cutoffs) {
  Figure f(842, 560);
  const auto& fi = c.model->info().faults.at(static_cast<std::size_t>(ff));
  auto fr = fault_frame(fs);
  title(f, (c.title.empty() ? c.model_name : c.title) + " - Allan diagram of fault " + fi.name,
        "Fault surface seen from the hanging wall; each cell: footwall unit (upper left) against hanging-wall unit (lower right); "
        "lines: hanging-wall (solid) and footwall (dashed) cut-offs");
  auto col = [&](int u) { return u >= 0 ? unit_color(c, static_cast<std::size_t>(u)) : rgb(215, 215, 215); };
  for (std::size_t j = 0; j + 1 < fs.nb; ++j)
    for (std::size_t i = 0; i + 1 < fs.na; ++i) {
      std::size_t k = i + fs.na * j;
      if (fs.extent[k] > 1 || fs.extent[k + 1] > 1 || fs.extent[k + fs.na] > 1 || fs.extent[k + fs.na + 1] > 1) continue;
      const auto& p = pts[k];
      if (!p.valid) continue;
      Pt a = fr.P(fs.a[i], fs.b[j]), b = fr.P(fs.a[i + 1], fs.b[j + 1]);
      f.path({{a.x, a.y}, {b.x + 0.3, a.y}, {a.x, b.y + 0.3}}, true, col(p.fw_unit), none());
      f.path({{b.x + 0.3, a.y}, {b.x + 0.3, b.y + 0.3}, {a.x, b.y + 0.3}}, true, col(p.hw_unit), none());
    }
  fault_frame_axes_and_cutoffs(f, fr, fi, cutoffs);
  auto es = unit_legend(c);
  es.push_back({"Cut-offs", "Hanging wall", none(), true, rgb(40, 40, 40), 0.8});
  es.push_back({"Cut-offs", "Footwall (dashed)", none(), true, rgb(40, 40, 40), 0.8});
  double cbx = 842 - 214;
  double ly = legend(f, cbx, 60, 196, es, 560.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, cbx, ly, 196);
  footer(f, c);
  return f;
}

Figure curvature_figure(const FigureContext& c, int h, const repr::CurvatureGrids& cg) {
  Figure f(842, 595);
  const auto& inf = c.model->info();
  const auto& d = inf.domain;
  double side = std::min(595.0 - 110, 842.0 - 280);
  double sc = side / std::max(d.sx, d.sy);
  double mx0 = 60, my0 = 60, mw = d.sx * sc, mh = d.sy * sc;
  auto P = [&](double x, double y) { return Pt{mx0 + (x - d.x0) * sc, my0 + mh - (y - d.y0) * sc}; };
  const auto& g = cg.mean;
  std::vector<double> mags;
  for (double v : g.v)
    if (std::isfinite(v)) mags.push_back(std::fabs(v));
  double kmax = 1e-6;
  if (!mags.empty()) {
    std::sort(mags.begin(), mags.end());
    kmax = std::max(1e-9, mags[static_cast<std::size_t>(0.98 * static_cast<double>(mags.size() - 1))]);
  }
  title(f, (c.title.empty() ? c.model_name : c.title) + " - " + inf.horizons[static_cast<std::size_t>(h)].name + " curvature",
        "Mean curvature of " + inf.horizons[static_cast<std::size_t>(h)].name +
            " (convex up red, concave up blue; fault gaps grey) and fracture trends (ticks) where the curvature passes the threshold");
  Image im;
  im.w = std::max(1, static_cast<int>(std::lround(mw * 3)));
  im.h = std::max(1, static_cast<int>(std::lround(mh * 3)));
  im.rgba.resize(static_cast<std::size_t>(im.w) * im.h * 4);
  const auto& gs = g.spec;
  for (int yy = 0; yy < im.h; ++yy)
    for (int xx = 0; xx < im.w; ++xx) {
      double x = d.x0 + (xx + 0.5) / im.w * d.sx, y = d.y0 + d.sy - (yy + 0.5) / im.h * d.sy;
      long i = std::clamp(std::lround((x - gs.x0) / gs.dx), 0L, static_cast<long>(gs.nx) - 1);
      long j = std::clamp(std::lround((y - gs.y0) / gs.dy), 0L, static_cast<long>(gs.ny) - 1);
      double v = g.v[gs.index(static_cast<std::size_t>(i), static_cast<std::size_t>(j))];
      Color cc = std::isnan(v) ? rgb(225, 225, 225) : diverging(v / kmax);
      std::uint8_t* p = &im.rgba[(static_cast<std::size_t>(yy) * im.w + xx) * 4];
      p[0] = cc.r; p[1] = cc.g; p[2] = cc.b; p[3] = 255;
    }
  const Image* img = f.own(std::move(im));
  f.clip(mx0, my0, mw, mh);
  f.image(mx0, my0, mw, mh, img);
  // fracture-trend ticks on a coarse lattice (about 25 across)
  const std::size_t stride = std::max<std::size_t>(1, std::max(gs.nx, gs.ny) / 25);
  const double len = 0.4 * static_cast<double>(stride) * std::min(gs.dx, gs.dy) * sc;
  std::size_t nticks = 0;
  for (std::size_t j = stride / 2; j < gs.ny; j += stride)
    for (std::size_t i = stride / 2; i < gs.nx; i += stride) {
      const double az = cg.trend.v[gs.index(i, j)];
      if (std::isnan(az)) continue;
      Pt p = P(gs.x0 + i * gs.dx, gs.y0 + j * gs.dy);
      const double ux = std::sin(az * kPi / 180) * len, uy = -std::cos(az * kPi / 180) * len;
      f.line({p.x - ux, p.y - uy}, {p.x + ux, p.y + uy}, rgb(20, 20, 20), 0.9);
      ++nticks;
    }
  for (auto& w : c.wells) {
    Pt p = P(w.x(), w.y());
    f.path({{p.x - 2.5, p.y}, {p.x, p.y - 2.5}, {p.x + 2.5, p.y}, {p.x, p.y + 2.5}}, true, rgb(255, 255, 255), kWell, 0.9);
    f.text({p.x + 4, p.y - 3}, w.name(), 7, Anchor::Start, true, kWell);
  }
  f.unclip();
  f.rect(mx0, my0, mw, mh, none(), kInk, 0.6);
  for (double t : nice_ticks(d.x0, d.x0 + d.sx, 6)) {
    Pt p = P(t, d.y0);
    f.line(p, {p.x, p.y + 3}, kInk, 0.5);
    f.text({p.x, p.y + 11}, km_or_m(t), 6.5, Anchor::Middle, false, kInk);
  }
  for (double t : nice_ticks(d.y0, d.y0 + d.sy, 6)) {
    Pt p = P(d.x0, t);
    f.line({p.x - 3, p.y}, p, kInk, 0.5);
    f.text({p.x - 5, p.y + 2.5}, km_or_m(t), 6.5, Anchor::End, false, kInk);
  }
  f.text({mx0 + mw / 2, my0 + mh + 24}, "Easting (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({mx0 - 40, my0 + mh / 2}, "Northing (m)", 7.5, Anchor::Middle, false, kInk, -90);
  double rx = 842 - 200, ry = 60;
  f.text({rx, ry}, "Mean curvature (1/km)", 7.5, Anchor::Start, true, kInk);
  for (int k = 0; k < 100; ++k) f.rect(rx, ry + 6 + k * 1.6, 14, 1.7, diverging(1 - 2 * k / 99.0), none());
  f.rect(rx, ry + 6, 14, 160, none(), kInk, 0.4);
  f.text({rx + 18, ry + 12}, "+" + fmt_fixed(kmax * 1000, 3), 6.5, Anchor::Start, false, kInk);
  f.text({rx + 18, ry + 88}, "0", 6.5, Anchor::Start, false, kInk);
  f.text({rx + 18, ry + 166}, "-" + fmt_fixed(kmax * 1000, 3), 6.5, Anchor::Start, false, kInk);
  std::vector<LegendEntry> es{{"Fractures", "Trend (" + std::to_string(nticks) + " ticks)", none(), true, rgb(20, 20, 20), 0.9},
                              {"Wells", "Well", rgb(255, 255, 255), false, kWell, 0.9, 1}};
  double ly = legend(f, rx, ry + 186, 180, es, 595.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, rx, ly, 180);
  footer(f, c);
  return f;
}

Figure conditioning_figure(const FigureContext& c, const std::vector<FittedRow>& rows, double chi2, std::size_t ndata) {
  const double H = std::max(240.0, 120.0 + 46.0 * static_cast<double>(rows.size()));
  Figure f(842, H);
  title(f, (c.title.empty() ? c.model_name : c.title) + " - fitted parameters",
        "Each bar is the admissible range; open diamond: start value; filled diamond: fitted value with +-2 sigma. Reduced chi2 " + fmt_fixed(chi2, 3) + " over " +
            std::to_string(ndata) + " data");
  const double x0 = 300, x1 = 760, y0 = 70;
  for (std::size_t k = 0; k < rows.size(); ++k) {
    const auto& r = rows[k];
    const double y = y0 + 46.0 * static_cast<double>(k) + 20;
    auto X = [&](double v) { return x0 + (v - r.min) / (r.max - r.min) * (x1 - x0); };
    f.text({x0 - 10, y - 4}, r.path, 7.5, Anchor::End, true, kInk);
    f.text({x0 - 10, y + 8}, fmt_fixed(r.value, 3) + " +- " + fmt_fixed(r.sigma, 3) + " " + r.unit + " (start " + fmt_shortest(r.start) + ")", 7,
           Anchor::End, false, kInk);
    f.rect(x0, y - 3, x1 - x0, 6, rgb(225, 225, 225), kInk, 0.4);
    f.text({x0, y + 16}, fmt_shortest(r.min), 6.5, Anchor::Start, false, kInk);
    f.text({x1, y + 16}, fmt_shortest(r.max) + " " + r.unit, 6.5, Anchor::End, false, kInk);
    double a = std::clamp(X(r.value - 2 * r.sigma), x0, x1), b = std::clamp(X(r.value + 2 * r.sigma), x0, x1);
    f.line({a, y}, {b, y}, rgb(20, 60, 140), 1.6);
    f.line({a, y - 4}, {a, y + 4}, rgb(20, 60, 140), 1.0);
    f.line({b, y - 4}, {b, y + 4}, rgb(20, 60, 140), 1.0);
    Pt ps{X(r.start), y};
    f.path({{ps.x - 3, ps.y}, {ps.x, ps.y - 3}, {ps.x + 3, ps.y}, {ps.x, ps.y + 3}}, true, rgb(255, 255, 255), kInk, 0.8);
    Pt pv{X(r.value), y};
    f.path({{pv.x - 3.5, pv.y}, {pv.x, pv.y - 3.5}, {pv.x + 3.5, pv.y}, {pv.x, pv.y + 3.5}}, true, rgb(20, 60, 140), rgb(20, 60, 140), 0.8);
  }
  footer(f, c);
  return f;
}

namespace {
Color walker_color(std::size_t w) {
  static const Color k[] = {rgb(20, 60, 140), rgb(200, 90, 20), rgb(30, 130, 70), rgb(150, 40, 120), rgb(90, 90, 90), rgb(0, 140, 160), rgb(170, 140, 0), rgb(120, 70, 30)};
  Color c = k[w % 8];
  c.a = 150;
  return c;
}
std::string diag_text(const PosteriorSamples& s, std::size_t p) {
  return "R-hat " + fmt_fixed(s.rhat[p], 4) + ", ESS " + fmt_fixed(s.ess[p], 0);
}
}  // namespace

Figure anomaly_map_figure(const FigureContext& c, const std::string& kind, const std::string& unit, const std::vector<AnomalyStation>& st, double chi2) {
  Figure f(842, 595);
  const auto& d = c.model->info().domain;
  title(f, (c.title.empty() ? c.model_name : c.title) + " - " + kind + " anomaly",
        "Observed and modelled anomaly (" + unit + ", contrast against the background) at " + std::to_string(st.size()) + " stations; residual = modelled - observed" +
            (std::isfinite(chi2) ? "; reduced chi2 " + fmt_fixed(chi2, 3) : std::string()));
  double vmin = 1e300, vmax = -1e300, rmax = 1e-12;
  for (auto& a : st) {
    vmin = std::min({vmin, a.observed, a.modelled});
    vmax = std::max({vmax, a.observed, a.modelled});
    rmax = std::max(rmax, std::fabs(a.modelled - a.observed));
  }
  if (!(vmax > vmin)) { vmin -= 1; vmax += 1; }
  rmax = std::max(rmax, 1e-3 * (vmax - vmin));  // residuals below 0.1 % of the range show white
  auto sig = [](double v) {
    if (v == 0) return std::string("0");
    const int dec = std::max(0, 2 - static_cast<int>(std::floor(std::log10(std::fabs(v)))));
    return fmt_fixed(v, std::min(dec, 8));
  };
  const double pw = 250, gap = 20, x00 = 30, y00 = 70;
  const double sc = std::min(pw / d.sx, 380.0 / d.sy), mw = d.sx * sc, mh = d.sy * sc;
  const char* names[] = {"Observed", "Modelled", "Residual"};
  for (int panel = 0; panel < 3; ++panel) {
    const double ox = x00 + panel * (pw + gap);
    f.rect(ox, y00, mw, mh, rgb(248, 248, 248), kInk, 0.6);
    f.text({ox + mw / 2, y00 - 6}, names[panel], 8, Anchor::Middle, true, kInk);
    const double ss = std::max(2.5, std::min(8.0, 0.6 * std::sqrt(mw * mh / std::max<std::size_t>(1, st.size()))));
    for (auto& a : st) {
      const double v = panel == 0 ? a.observed : panel == 1 ? a.modelled : a.modelled - a.observed;
      const Color col = panel < 2 ? viridis_like((v - vmin) / (vmax - vmin)) : diverging(v / rmax);
      const double px = ox + (a.x - d.x0) * sc, py = y00 + mh - (a.y - d.y0) * sc;
      f.rect(px - ss / 2, py - ss / 2, ss, ss, col, none());
    }
    // colour bar
    const double by = y00 + mh + 18;
    for (int k = 0; k < 100; ++k) {
      const double t = k / 99.0;
      f.rect(ox + t * (mw - 2), by, mw / 99.0 + 0.3, 8, panel < 2 ? viridis_like(t) : diverging(2 * t - 1), none());
    }
    const double lo = panel < 2 ? vmin : -rmax, hi = panel < 2 ? vmax : rmax;
    f.text({ox, by + 18}, sig(lo), 6.5, Anchor::Start, false, kInk);
    f.text({ox + mw, by + 18}, sig(hi) + " " + unit, 6.5, Anchor::End, false, kInk);
  }
  footer(f, c);
  return f;
}

Figure anomaly_profile_figure(const FigureContext& c, const std::string& kind, const std::string& unit, const AnomalyProfile& p) {
  Figure f(842, 420);
  title(f, (c.title.empty() ? c.model_name : c.title) + " - " + kind + " profile " + p.section,
        "Modelled anomaly along the section (line) and observed stations within " + fmt_fixed(p.band, 0) + " m of it (dots), " + unit);
  const double x0 = 80, x1 = 800, y0 = 60, y1 = 340;
  if (p.s.empty()) { footer(f, c); return f; }
  double lo = 1e300, hi = -1e300;
  for (double v : p.modelled) { lo = std::min(lo, v); hi = std::max(hi, v); }
  for (auto& o : p.observed) { lo = std::min(lo, o.second); hi = std::max(hi, o.second); }
  if (!(hi > lo)) { lo -= 1; hi += 1; }
  const double pad = 0.08 * (hi - lo), L = std::max(1e-9, p.s.back());
  lo -= pad;
  hi += pad;
  auto X = [&](double s) { return x0 + s / L * (x1 - x0); };
  auto Y = [&](double v) { return y1 - (v - lo) / (hi - lo) * (y1 - y0); };
  f.rect(x0, y0, x1 - x0, y1 - y0, rgb(252, 252, 252), kInk, 0.6);
  for (double t : nice_ticks(lo, hi, 6)) {
    f.line({x0, Y(t)}, {x1, Y(t)}, kGrid, 0.3);
    f.text({x0 - 4, Y(t) + 2.5}, fmt_shortest(t), 6.5, Anchor::End, false, kInk);
  }
  for (double t : nice_ticks(0, L, 8)) f.text({X(t), y1 + 11}, fmt_fixed(t, 0), 6.5, Anchor::Middle, false, kInk);
  std::vector<Pt> line;
  for (std::size_t k = 0; k < p.s.size(); ++k) line.push_back({X(p.s[k]), Y(p.modelled[k])});
  f.path(line, false, none(), rgb(20, 60, 140), 1.4);
  for (auto& o : p.observed) f.path({{X(o.first) - 2.5, Y(o.second)}, {X(o.first), Y(o.second) - 2.5}, {X(o.first) + 2.5, Y(o.second)}, {X(o.first), Y(o.second) + 2.5}},
                                    true, rgb(200, 90, 20), kInk, 0.4);
  f.text({(x0 + x1) / 2, y1 + 24}, "Distance along " + p.section + " (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({x0 - 45, (y0 + y1) / 2}, kind + " anomaly (" + unit + ")", 7.5, Anchor::Middle, false, kInk, -90);
  footer(f, c);
  return f;
}

Figure wheeler_figure(const FigureContext& c, const WheelerDiagram& w) {
  Figure f(842, 595);
  const auto& inf = c.model->info();
  title(f, (c.title.empty() ? c.model_name : c.title) + " - Wheeler diagram, section " + w.name,
        std::string("Distance along the section against geological time") + (w.ordinal ? " (horizon ordinals: no ages given)" : " (Ma)") +
            "; blank: hiatus (no record); grey: fault cut-out; dashed: horizon ages");
  const double x0 = 80, x1 = 842 - 190, y0 = 60, y1 = 595 - 70;
  if (w.columns.empty() || w.ages.empty()) {
    footer(f, c);
    return f;
  }
  const double half = w.ages.size() > 1 ? 0.5 * (w.ages.back() - w.ages.front()) / static_cast<double>(w.ages.size() - 1) : 0.5;
  const double tmin = w.ages.front() - half, tmax = w.ages.back() + half, L = std::max(1e-9, w.columns.back().station);
  auto X = [&](double s) { return x0 + s / L * (x1 - x0); };
  auto Y = [&](double t) { return y0 + (t - tmin) / (tmax - tmin) * (y1 - y0); };  // young at the top
  f.rect(x0, y0, x1 - x0, y1 - y0, rgb(255, 255, 255), kInk, 0.6);
  const double cw = (x1 - x0) / static_cast<double>(std::max<std::size_t>(1, w.columns.size() - 1));
  std::set<int> units;
  for (std::size_t k = 0; k < w.columns.size(); ++k)
    for (auto& iv : w.columns[k].intervals) {
      const double xa = std::max(x0, X(w.columns[k].station) - 0.5 * cw), xb = std::min(x1, X(w.columns[k].station) + 0.5 * cw);
      const Color col = iv.unit == kWheelerFault ? rgb(150, 150, 150) : unit_color(c, static_cast<std::size_t>(std::max(0, iv.unit)));
      f.rect(xa, Y(iv.t_young), xb - xa + 0.3, Y(iv.t_old) - Y(iv.t_young), col, none());
      units.insert(iv.unit);
    }
  for (std::size_t h = 0; h < w.ages.size(); ++h) {
    f.line({x0, Y(w.ages[h])}, {x1, Y(w.ages[h])}, rgb(90, 90, 90), 0.5, {3, 2});
    if (h < inf.horizons.size()) f.text({x1 + 4, Y(w.ages[h]) + 2.5}, inf.horizons[h].name, 6.5, Anchor::Start, false, kInk);
  }
  f.rect(x0, y0, x1 - x0, y1 - y0, none(), kInk, 0.6);
  for (double t : nice_ticks(tmin, tmax, 8)) f.text({x0 - 5, Y(t) + 2.5}, fmt_shortest(t), 6.5, Anchor::End, false, kInk);
  for (double s : nice_ticks(0, L, 8)) f.text({X(s), y1 + 11}, fmt_fixed(s, 0), 6.5, Anchor::Middle, false, kInk);
  f.text({(x0 + x1) / 2, y1 + 24}, "Distance along " + w.name + " (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({x0 - 40, (y0 + y1) / 2}, w.ordinal ? "Geological time (horizon ordinal)" : "Geological time (Ma)", 7.5, Anchor::Middle, false, kInk, -90);
  double ly = y0 + 40;
  f.text({x1 + 70, y0 + 26}, "Units", 7.5, Anchor::Start, true, kInk);
  if (units.count(kWheelerFault)) {
    f.rect(x1 + 70, ly - 7, 14, 8, rgb(150, 150, 150), rgb(90, 90, 90), 0.3);
    f.text({x1 + 88, ly}, "fault cut-out (not a hiatus)", 6.5, Anchor::Start, false, kInk);
    ly += 12;
  }
  for (int u : units) {
    if (u < 0 || static_cast<std::size_t>(u) >= inf.units.size()) continue;
    f.rect(x1 + 70, ly - 7, 14, 8, unit_color(c, static_cast<std::size_t>(u)), rgb(90, 90, 90), 0.3);
    f.text({x1 + 88, ly}, inf.units[static_cast<std::size_t>(u)].name, 6.5, Anchor::Start, false, kInk);
    ly += 12;
  }
  if (w.undefined > 0) f.text({x0, y1 + 36}, std::to_string(w.undefined) + " samples outside the defined region (no column)", 6.5, Anchor::Start, false, kFault);
  footer(f, c);
  return f;
}

Figure posterior_trace_figure(const FigureContext& c, const PosteriorSamples& s) {
  const std::size_t np = s.names.size();
  const double rowh = 110, H = 110 + rowh * static_cast<double>(np);
  Figure f(842, H);
  title(f, (c.title.empty() ? c.model_name : c.title) + " - posterior traces",
        std::to_string(s.walkers) + " walkers, " + std::to_string(s.steps) + " steps after " + std::to_string(s.burn_in) + " burn-in, seed " +
            std::to_string(s.seed) + ", acceptance " + fmt_fixed(s.acceptance, 2) + (s.converged ? "; converged (R-hat < 1.01, ESS >= 100)" : "; NOT converged"));
  const double x0 = 90, x1 = 650;
  const std::size_t nw = std::min<std::size_t>(16, static_cast<std::size_t>(s.walkers));
  const std::size_t stride = std::max<std::size_t>(1, static_cast<std::size_t>(s.steps) / 400);
  for (std::size_t p = 0; p < np; ++p) {
    const double ytop = 60 + rowh * static_cast<double>(p), ybot = ytop + rowh - 30;
    double lo = 1e300, hi = -1e300;
    for (auto& v : s.chain) { lo = std::min(lo, v[p]); hi = std::max(hi, v[p]); }
    if (!(hi > lo)) { lo -= 1; hi += 1; }
    const double pad = 0.05 * (hi - lo);
    lo -= pad;
    hi += pad;
    auto X = [&](double st) { return x0 + st / std::max(1.0, static_cast<double>(s.steps - 1)) * (x1 - x0); };
    auto Y = [&](double v) { return ybot - (v - lo) / (hi - lo) * (ybot - ytop); };
    f.rect(x0, ytop, x1 - x0, ybot - ytop, rgb(250, 250, 250), kInk, 0.5);
    for (double t : nice_ticks(lo, hi, 4)) {
      f.line({x0, Y(t)}, {x1, Y(t)}, kGrid, 0.3);
      f.text({x0 - 4, Y(t) + 2.5}, fmt_shortest(t), 6.5, Anchor::End, false, kInk);
    }
    f.clip(x0, ytop, x1 - x0, ybot - ytop);
    for (std::size_t w = 0; w < nw; ++w) {
      std::vector<Pt> pts;
      for (std::size_t st = 0; st < static_cast<std::size_t>(s.steps); st += stride)
        pts.push_back({X(static_cast<double>(st)), Y(s.chain[w * static_cast<std::size_t>(s.steps) + st][p])});
      f.path(pts, false, none(), walker_color(w), 0.5);
    }
    f.line({x0, Y(s.mean[p])}, {x1, Y(s.mean[p])}, kInk, 0.9, {3, 2});
    f.unclip();
    f.text({x0, ytop - 4}, s.names[p] + " (" + s.units[p] + ")", 7.5, Anchor::Start, true, kInk);
    double tx = x1 + 12;
    f.text({tx, ytop + 12}, "mean " + fmt_fixed(s.mean[p], 3) + " " + s.units[p], 7, Anchor::Start, false, kInk);
    f.text({tx, ytop + 24}, "95 %: " + fmt_fixed(s.q025[p], 3) + " to " + fmt_fixed(s.q975[p], 3), 7, Anchor::Start, false, kInk);
    const bool okp = s.rhat[p] < 1.01 && s.ess[p] >= 100;
    f.text({tx, ytop + 36}, diag_text(s, p), 7, Anchor::Start, !okp, okp ? kInk : kFault);
    if (p + 1 == np)
      for (double t : nice_ticks(0, static_cast<double>(s.steps - 1), 6)) f.text({X(t), ybot + 11}, fmt_fixed(t, 0), 6.5, Anchor::Middle, false, kInk);
  }
  f.text({(x0 + x1) / 2, 60 + rowh * static_cast<double>(np) - 4}, "Step after burn-in (" + std::to_string(nw) + " of " + std::to_string(s.walkers) +
                                                                     " walkers drawn; dashed: posterior mean)",
         7.5, Anchor::Middle, false, kInk);
  footer(f, c);
  return f;
}

Figure posterior_corner_figure(const FigureContext& c, const PosteriorSamples& s) {
  const std::size_t np = std::min<std::size_t>(6, s.names.size());
  Figure f(842, 595);
  title(f, (c.title.empty() ? c.model_name : c.title) + " - posterior",
        "Marginal histograms (shaded: 95 % credible interval; dashed: least-squares value) and pairwise samples, every " + std::to_string(s.thin) +
            "th step of every walker" + (s.names.size() > np ? " (first 6 parameters)" : ""));
  const double gx = 90, gy = 60, side = std::min(842.0 - 140, 595.0 - 120), cell = side / static_cast<double>(np), gap = 8;
  std::vector<double> lo(np), hi(np);
  for (std::size_t p = 0; p < np; ++p) {
    lo[p] = 1e300;
    hi[p] = -1e300;
    for (auto& v : s.chain) { lo[p] = std::min(lo[p], v[p]); hi[p] = std::max(hi[p], v[p]); }
    if (!(hi[p] > lo[p])) { lo[p] -= 1; hi[p] += 1; }
  }
  std::vector<std::size_t> pick;
  for (std::size_t w = 0; w < static_cast<std::size_t>(s.walkers); ++w)
    for (std::size_t st = static_cast<std::size_t>(s.thin) - 1; st < static_cast<std::size_t>(s.steps); st += static_cast<std::size_t>(s.thin))
      pick.push_back(w * static_cast<std::size_t>(s.steps) + st);
  const std::size_t every = std::max<std::size_t>(1, pick.size() / 3000);
  for (std::size_t i = 0; i < np; ++i)
    for (std::size_t j = 0; j <= i; ++j) {
      const double cx = gx + cell * static_cast<double>(j), cy = gy + cell * static_cast<double>(i), w = cell - gap, h = cell - gap;
      auto X = [&](double v) { return cx + (v - lo[j]) / (hi[j] - lo[j]) * w; };
      f.rect(cx, cy, w, h, rgb(250, 250, 250), kInk, 0.5);
      if (i == j) {
        const int nb = 30;
        std::vector<double> cnt(nb, 0.0);
        for (std::size_t k : pick) {
          int b = static_cast<int>((s.chain[k][j] - lo[j]) / (hi[j] - lo[j]) * nb);
          cnt[static_cast<std::size_t>(std::clamp(b, 0, nb - 1))] += 1;
        }
        double cm = *std::max_element(cnt.begin(), cnt.end());
        f.rect(X(s.q025[j]), cy, X(s.q975[j]) - X(s.q025[j]), h, rgb(220, 230, 245), none());
        for (int b = 0; b < nb; ++b) {
          double bh = cm > 0 ? cnt[static_cast<std::size_t>(b)] / cm * (h - 6) : 0;
          f.rect(cx + w * b / nb, cy + h - bh, w / nb, bh, rgb(20, 60, 140), none());
        }
        if (j < s.lsq.size() && std::isfinite(s.lsq[j])) f.line({X(s.lsq[j]), cy}, {X(s.lsq[j]), cy + h}, kFault, 0.9, {3, 2});
        f.text({cx + w / 2, cy - 3}, s.names[j], 6.5, Anchor::Middle, true, kInk);
      } else {
        auto Y = [&](double v) { return cy + h - (v - lo[i]) / (hi[i] - lo[i]) * h; };
        f.clip(cx, cy, w, h);
        for (std::size_t q = 0; q < pick.size(); q += every) {
          const auto& v = s.chain[pick[q]];
          f.rect(X(v[j]) - 0.6, Y(v[i]) - 0.6, 1.2, 1.2, Color{20, 60, 140, 90}, none());
        }
        f.unclip();
      }
      if (i + 1 == np)
        for (double t : nice_ticks(lo[j], hi[j], 3)) f.text({X(t), cy + h + 10}, fmt_shortest(t), 6, Anchor::Middle, false, kInk);
      if (j == 0 && i > 0) {
        auto Y = [&](double v) { return cy + h - (v - lo[i]) / (hi[i] - lo[i]) * h; };
        for (double t : nice_ticks(lo[i], hi[i], 3)) f.text({cx - 3, Y(t) + 2}, fmt_shortest(t), 6, Anchor::End, false, kInk);
      }
    }
  const double tx = gx + side + 14;
  double ty = gy + 10;
  for (std::size_t p = 0; p < s.names.size(); ++p) {
    f.text({tx, ty}, s.names[p] + ": " + fmt_fixed(s.mean[p], 3) + " " + s.units[p], 7, Anchor::Start, true, kInk);
    f.text({tx, ty + 10}, "95 % " + fmt_fixed(s.q025[p], 3) + " to " + fmt_fixed(s.q975[p], 3) + "; " + diag_text(s, p), 6.5, Anchor::Start, false, kInk);
    ty += 26;
  }
  footer(f, c);
  return f;
}

Figure misfit_map_figure(const FigureContext& c, const std::string& horizon, const std::vector<MisfitPoint>& pts) {
  Figure f(842, 595);
  const auto& d = c.model->info().domain;
  double side = std::min(595.0 - 110, 842.0 - 280), sc = side / std::max(d.sx, d.sy);
  double mx0 = 60, my0 = 60, mw = d.sx * sc, mh = d.sy * sc;
  auto P = [&](double x, double y) { return Pt{mx0 + (x - d.x0) * sc, my0 + mh - (y - d.y0) * sc}; };
  const bool norm = std::any_of(pts.begin(), pts.end(), [](const MisfitPoint& p) { return !std::isnan(p.normalised); });
  double vmax = 3.0;
  if (!norm) {
    vmax = 1e-3;
    for (auto& p : pts) vmax = std::max(vmax, std::fabs(p.misfit));
  }
  std::size_t beyond = 0;
  for (auto& p : pts) beyond += norm && std::fabs(p.normalised) > 3;
  title(f, (c.title.empty() ? c.model_name : c.title) + " - data misfit on " + horizon,
        std::string(norm ? "Normalised misfit (misfit / stated uncertainty), -3 to 3" : "Misfit (m)") + " at every datum of " + horizon +
            "; positive: the model lies above the datum");
  f.rect(mx0, my0, mw, mh, rgb(248, 248, 248), kInk, 0.6);
  for (auto& p : pts) {
    double v = norm ? p.normalised : p.misfit;
    Color cc = std::isnan(v) ? rgb(150, 150, 150) : diverging(v / vmax);
    Pt q = P(p.x, p.y);
    f.rect(q.x - 1.6, q.y - 1.6, 3.2, 3.2, cc, std::fabs(v) > vmax ? kInk : none(), 0.4);
  }
  for (double t : nice_ticks(d.x0, d.x0 + d.sx, 6)) f.text({P(t, d.y0).x, my0 + mh + 11}, km_or_m(t), 6.5, Anchor::Middle, false, kInk);
  for (double t : nice_ticks(d.y0, d.y0 + d.sy, 6)) f.text({mx0 - 5, P(d.x0, t).y + 2.5}, km_or_m(t), 6.5, Anchor::End, false, kInk);
  f.text({mx0 + mw / 2, my0 + mh + 24}, "Easting (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({mx0 - 40, my0 + mh / 2}, "Northing (m)", 7.5, Anchor::Middle, false, kInk, -90);
  double rx = 842 - 200, ry = 60;
  f.text({rx, ry}, norm ? "Misfit / sigma" : "Misfit (m)", 7.5, Anchor::Start, true, kInk);
  for (int k = 0; k < 100; ++k) f.rect(rx, ry + 6 + k * 1.6, 14, 1.7, diverging(1 - 2 * k / 99.0), none());
  f.rect(rx, ry + 6, 14, 160, none(), kInk, 0.4);
  f.text({rx + 18, ry + 12}, "+" + fmt_shortest(vmax), 6.5, Anchor::Start, false, kInk);
  f.text({rx + 18, ry + 88}, "0", 6.5, Anchor::Start, false, kInk);
  f.text({rx + 18, ry + 166}, "-" + fmt_shortest(vmax), 6.5, Anchor::Start, false, kInk);
  f.text({rx, ry + 186}, std::to_string(pts.size()) + " data" + (norm ? ", " + std::to_string(beyond) + " beyond 3 sigma (outlined)" : std::string()), 7,
         Anchor::Start, false, kInk);
  footer(f, c);
  return f;
}

WrittenFigure write_figure(const Figure& f, const std::filesystem::path& folder, const std::string& base, const std::string& t,
                           const std::string& caption, const std::string& notes) {
  WrittenFigure w;
  w.basename = base;
  w.title = t;
  w.caption = caption;
  std::filesystem::create_directories(folder);
  write_file(folder / (base + ".svg"), to_svg(f, t));
  w.files.push_back(base + ".svg");
  if (raster_available()) {
    Raster r = rasterise(f, 300);
    std::string png = encode_png(r, 300);
    write_file(folder / (base + ".png"), png);
    write_file(folder / (base + ".tif"), encode_tiff_lzw(r, 300, "Stratonix " + build_info().version));
    Slide s{t, caption, notes, png, r.w, r.h};
    write_file(folder / (base + ".pptx"), make_pptx({s}, t));
    w.files.push_back(base + ".png");
    w.files.push_back(base + ".tif");
    w.files.push_back(base + ".pptx");
  } else {
    log(LogLevel::Warn, "Cairo not available: PNG, TIFF and PPTX not written for " + base);
  }
  return w;
}

std::vector<WrittenFigure> write_figure_set(const FigureContext& c, const std::vector<std::pair<std::string, Polyline>>& sections,
                                            const std::vector<std::vector<double>>& detail_windows,
                                            const std::filesystem::path& folder) {
  std::vector<WrittenFigure> out;
  const auto& inf = c.model->info();
  const Domain& d = inf.domain;
  std::string basis = basis_label(inf.basis);
  auto notes_common = "Evidence basis: " + basis + ". Coordinate system: " + inf.crs.label() + "; vertical datum: " + inf.crs.vertical_datum +
                      ". Built by Stratonix " + build_info().version + " from model " + c.model_name + " v" + c.model_version + ".";
  // Sections
  for (std::size_t k = 0; k < sections.size(); ++k) {
    auto s = repr::make_section(*c.model, sections[k].first, sections[k].second, std::max(5.0, c.grid_spacing / 4), 2.0);
    auto fgr = section_figure(c, s, k < detail_windows.size() && !detail_windows[k].empty() ? &detail_windows[k] : nullptr);
    out.push_back(write_figure(fgr, folder, figure_basename(c, "section", sections[k].first), c.model_name + " - section " + sections[k].first,
                               "Section " + sections[k].first + ": units per fault block, fault traces between tips, wells within 300 m.",
                               notes_common + " Section length " + fmt_fixed(s.length, 0) + " m."));
  }
  // Maps per horizon
  GridSpec g = GridSpec::covering(d, std::max(c.grid_spacing, std::max(d.sx, d.sy) / 600));
  // Stratal terminations per unit (stage 2C); empty for models without them.
  std::vector<std::vector<repr::TerminationLine>> terms(inf.units.size());
  for (std::size_t u = 1; u + 1 < inf.units.size(); ++u)
    if (!inf.units[u].termination.empty()) terms[u] = repr::termination_lines(*c.model, static_cast<int>(u), g);
  for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
    std::vector<repr::TerminationLine> adj = terms[h];
    adj.insert(adj.end(), terms[h + 1].begin(), terms[h + 1].end());
    repr::FaultedSurface scratch;
    const auto& fs = repr::faulted_surface(*c.model, static_cast<int>(h), g, c.surfaces, scratch);
    std::vector<std::vector<repr::CutoffLine>> cut;
    for (std::size_t ff = 0; ff < inf.faults.size(); ++ff) {
      auto fsurf = repr::fault_surface(*c.model, static_cast<int>(ff), std::max(10.0, g.dx / 2));
      cut.push_back(repr::cutoff_lines(*c.model, fsurf, static_cast<int>(h)));
    }
    repr::Closure cl;
    bool has = static_cast<int>(h) == c.closure_horizon;
    if (has) cl = repr::closure_analysis(fs.z, c.faults_seal);
    auto fgr = map_figure(c, static_cast<int>(h), fs, cut, has ? &cl : nullptr, adj.empty() ? nullptr : &adj);
    out.push_back(write_figure(fgr, folder, figure_basename(c, "map", inf.horizons[h].name), c.model_name + " - " + inf.horizons[h].name + " depth map",
                               "Depth map of " + inf.horizons[h].name + " with fault polygons and wells" + (has ? ", closure outline and spill point." : ".") +
                                   (adj.empty() ? "" : " Stratal termination lines of the units above and below it."),
                               notes_common));
  }
  // Isopach maps (stage 2C): units that terminate, and units asked for.
  for (std::size_t u = 1; u + 1 < inf.units.size(); ++u) {
    bool asked = std::find(c.isopach_units.begin(), c.isopach_units.end(), static_cast<int>(u)) != c.isopach_units.end();
    if (inf.units[u].termination.empty() && !asked) continue;
    auto t = repr::unit_thickness_grid(*c.model, static_cast<int>(u), g);
    auto fgr = isopach_figure(c, static_cast<int>(u), t, terms[u]);
    out.push_back(write_figure(fgr, folder, figure_basename(c, "isopach", inf.units[u].name), c.model_name + " - " + inf.units[u].name + " isopach",
                               "Isopach (vertical thickness) of " + inf.units[u].name + " per fault block, with its stratal termination lines.",
                               notes_common));
  }
  // Burial history (stage 2D)
  for (auto& b : c.burial) {
    auto fgr = burial_figure(c, b);
    out.push_back(write_figure(fgr, folder, figure_basename(c, "burial", b.location), c.model_name + " - burial history at " + b.location,
                               "Burial-history curves at " + b.location + ": burial depth of each horizon through the deposition steps.",
                               notes_common));
    for (auto& fd : b.fields)
      out.push_back(write_figure(burial_field_figure(c, b, fd), folder, figure_basename(c, "burial_" + fd.kind, b.location),
                                 c.model_name + " - " + fd.label + " through burial at " + b.location,
                                 "Burial-history curves at " + b.location + " coloured by " + fd.label + " (" + fd.unit + ").", notes_common));
  }
  // Curvature maps (stage 2F)
  for (int h : c.curvature_horizons) {
    auto cg = repr::curvature(*c.model, h, g, std::isnan(c.fracture_threshold) ? 1e300 : c.fracture_threshold);
    auto fgr = curvature_figure(c, h, cg);
    const auto& hn = inf.horizons[static_cast<std::size_t>(h)].name;
    out.push_back(write_figure(fgr, folder, figure_basename(c, "curvature", hn), c.model_name + " - " + hn + " curvature",
                               "Mean curvature of " + hn + " (Roberts 2001) with fracture trends where the curvature passes the threshold.",
                               notes_common));
  }
  // 3D views
  for (auto [az, el, name] : {std::tuple<double, double, std::string>{135, 22, "SE-oblique"}, {315, 22, "NW-oblique"}}) {
    auto fgr = view3d_figure(c, az, el, name);
    out.push_back(write_figure(fgr, folder, figure_basename(c, "3d", name), c.model_name + " - 3D view " + name,
                               "Horizons and fault planes, " + name + " view, vertical exaggeration " + fmt_fixed(c.vertical_exaggeration, 1) + "x.",
                               notes_common));
  }
  // Fence through wells (needs two or more wells)
  if (c.wells.size() >= 2) {
    std::vector<repr::Section> legs;
    for (std::size_t k = 0; k + 1 < c.wells.size(); ++k) {
      Polyline pl{{c.wells[k].x(), c.wells[k].y(), 0}, {c.wells[k + 1].x(), c.wells[k + 1].y(), 0}};
      legs.push_back(repr::make_section(*c.model, c.wells[k].name() + "-" + c.wells[k + 1].name(), pl, std::max(10.0, c.grid_spacing / 2), 4.0));
    }
    auto fgr = fence_figure(c, legs);
    out.push_back(write_figure(fgr, folder, figure_basename(c, "fence", "wells"), c.model_name + " - fence diagram through wells",
                               "Fence diagram through the wells in description order.", notes_common));
  }
  // Fault-plane displacement maps and throw along strike of linked fault systems (stage 2A).
  if (!inf.faults.empty()) {
    const double fsp = std::max(10.0, std::max(d.sx, d.sy) / 300);
    std::vector<repr::FaultSurface> surfs;
    std::vector<std::vector<std::vector<repr::CutoffLine>>> cuts(inf.faults.size());
    for (std::size_t ff = 0; ff < inf.faults.size(); ++ff) {
      surfs.push_back(repr::fault_surface(*c.model, static_cast<int>(ff), fsp));
      for (std::size_t h = 0; h < inf.horizons.size(); ++h) cuts[ff].push_back(repr::cutoff_lines(*c.model, surfs[ff], static_cast<int>(h)));
      FaultEval fe;
      c.model->evaluate_fault(static_cast<int>(ff), surfs[ff].nodes, fe);
      auto fgr = fault_plane_figure(c, static_cast<int>(ff), surfs[ff], fe.displacement, cuts[ff]);
      out.push_back(write_figure(fgr, folder, figure_basename(c, "fault-plane", inf.faults[ff].name),
                                 c.model_name + " - displacement on fault " + inf.faults[ff].name,
                                 "Displacement on fault " + inf.faults[ff].name + " with the cut-offs of every horizon.", notes_common));
      // fault rock and stress (stage 2F)
      if (c.fault_rock || c.stress.present) {
        const auto& fname = inf.faults[ff].name;
        auto pts = repr::fault_rock_surface(*c.model, surfs[ff], c.rock, c.stress);
        auto attr = [&](const char* kind, const std::string& label, const std::string& unit, auto get, double scale, double lo, double hi,
                        const std::string& caption) {
          std::vector<double> v(pts.size(), kNaN);
          bool any = false;
          for (std::size_t k = 0; k < pts.size(); ++k)
            if (pts[k].valid) {
              v[k] = get(pts[k]) * scale;
              any = any || std::isfinite(v[k]);
            }
          if (!any) return;
          auto fg = fault_attribute_figure(c, static_cast<int>(ff), surfs[ff], v, cuts[ff], label, unit, lo, hi);
          out.push_back(write_figure(fg, folder, figure_basename(c, kind, fname), c.model_name + " - " + label + " on fault " + fname,
                                     caption, notes_common));
        };
        if (c.fault_rock) {
          auto fg = allan_figure(c, static_cast<int>(ff), surfs[ff], pts, cuts[ff]);
          out.push_back(write_figure(fg, folder, figure_basename(c, "allan", fname), c.model_name + " - Allan diagram of fault " + fname,
                                     "Allan diagram of fault " + fname + ": units juxtaposed across the fault, with the cut-offs of every horizon.",
                                     notes_common));
          attr("fault-sgr", "Shale gouge ratio", "fraction", [](const repr::FaultRockPoint& p) { return p.sgr; }, 1.0, 0.0, 1.0,
               "Shale gouge ratio on fault " + fname + " (Yielding et al. 1997) from the clay fractions of the units slipped past each point.");
          attr("fault-csp", "Clay smear potential", "m", [](const repr::FaultRockPoint& p) { return p.csp; }, 1.0, kNaN, kNaN,
               "Clay smear potential on fault " + fname + " (Bouvier et al. 1989) from the shale beds of both walls.");
          attr("fault-seal", "Seal capacity", "bar", [](const repr::FaultRockPoint& p) { return p.seal_bar; }, 1.0, kNaN, kNaN,
               "Seal capacity (across-fault pressure difference) on fault " + fname + " from the SGR calibration in the description.");
        }
        if (c.stress.present) {
          attr("fault-ts", "Slip tendency", "ratio", [](const repr::FaultRockPoint& p) { return p.ts; }, 1.0, 0.0, kNaN,
               "Slip tendency (Morris et al. 1996) on fault " + fname + " under the stress in the description.");
          attr("fault-td", "Dilation tendency", "ratio", [](const repr::FaultRockPoint& p) { return p.td; }, 1.0, 0.0, 1.0,
               "Dilation tendency (Ferrill et al. 1999) on fault " + fname + " under the stress in the description.");
          attr("fault-dp", "Pressure to reactivate", "MPa", [](const repr::FaultRockPoint& p) { return p.dp; }, 1e-6, kNaN, kNaN,
               "Pore-pressure increase to reactivation (Mohr-Coulomb) on fault " + fname + ".");
        }
      }
    }
    // fault systems: faults joined by recorded relations
    std::vector<int> group(inf.faults.size());
    for (std::size_t ff = 0; ff < group.size(); ++ff) group[ff] = static_cast<int>(ff);
    auto find = [&](int x) { while (group[static_cast<std::size_t>(x)] != x) x = group[static_cast<std::size_t>(x)]; return x; };
    for (std::size_t ff = 0; ff < inf.faults.size(); ++ff)
      for (auto& r : inf.faults[ff].relations) {
        std::string other = r.substr(r.rfind(' ') + 1);
        for (std::size_t gg = 0; gg < inf.faults.size(); ++gg)
          if (inf.faults[gg].name == other) group[static_cast<std::size_t>(find(static_cast<int>(ff)))] = find(static_cast<int>(gg));
      }
    for (std::size_t root = 0; root < inf.faults.size(); ++root) {
      std::vector<int> members;
      for (std::size_t ff = 0; ff < inf.faults.size(); ++ff)
        if (find(static_cast<int>(ff)) == static_cast<int>(root)) members.push_back(static_cast<int>(ff));
      if (members.size() < 2) continue;
      const auto& f0 = inf.faults[static_cast<std::size_t>(members[0])];
      std::vector<StrikeThrowSeries> series;
      std::string name;
      for (int m : members) name += (name.empty() ? "" : "-") + inf.faults[static_cast<std::size_t>(m)].name;
      for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
        // each fault's throw against the strike coordinate of the first fault
        std::vector<std::vector<std::pair<double, double>>> per;
        double lo = 1e300, hi = -1e300;
        for (int m : members) {
          std::vector<std::pair<double, double>> pts;
          const repr::CutoffLine* hw = nullptr;
          const repr::CutoffLine* fw = nullptr;
          for (auto& l : cuts[static_cast<std::size_t>(m)][h]) (l.hanging_wall ? hw : fw) = &l;
          if (hw && fw)
            for (std::size_t i = 0; i < hw->by_column.size() && i < fw->by_column.size(); ++i) {
              const Vec3 &a = hw->by_column[i], &b = fw->by_column[i];
              if (std::isnan(a.x) || std::isnan(b.x)) continue;
              // the footwall cut-off: those of linked faults meet at the branch line
              pts.push_back({dot(b - f0.ref_point, f0.strike_dir), b.z - a.z});
            }
          std::sort(pts.begin(), pts.end());
          if (!pts.empty()) { lo = std::min(lo, pts.front().first); hi = std::max(hi, pts.back().first); }
          per.push_back(pts);
        }
        if (!(hi > lo)) continue;
        Color col = unit_color(c, h + 1);
        col = rgb(col.r * 6 / 10, col.g * 6 / 10, col.b * 6 / 10);
        StrikeThrowSeries sum{inf.horizons[h].name, col, true, {}};
        // Where a fault's throw jumps between neighbouring samples (the branch line of a linked
        // splay hands throw over between faults) the sum of separately sampled faults is not
        // defined: the summed line is left open there instead of drawing a spike.
        std::vector<std::pair<double, double>> jumps;  // (position, half width)
        for (auto& pts : per)
          for (std::size_t q = 1; q < pts.size(); ++q)
            if (std::fabs(pts[q].second - pts[q - 1].second) > 5.0)
              jumps.push_back({0.5 * (pts[q].first + pts[q - 1].first), 1.5 * (pts[q].first - pts[q - 1].first) + 1.0});
        for (int k = 0; k <= 400; ++k) {
          double t = lo + (hi - lo) * k / 400.0, v = 0;
          bool open = false;
          for (auto& [jt, jw] : jumps) open = open || std::fabs(t - jt) <= jw;
          if (open) { sum.points.push_back({t, std::nan("")}); continue; }
          for (auto& pts : per) {
            if (pts.size() < 2 || t < pts.front().first || t > pts.back().first) continue;
            auto it = std::lower_bound(pts.begin(), pts.end(), std::pair<double, double>{t, -1e300});
            if (it == pts.begin()) { v += it->second; continue; }
            auto pv = it - 1;
            double w = it->first > pv->first ? (t - pv->first) / (it->first - pv->first) : 0;
            v += pv->second + w * (it->second - pv->second);
          }
          sum.points.push_back({t, v});
        }
        series.push_back(sum);
        if (h == 0 || per.empty()) {
          for (std::size_t q = 0; q < per.size(); ++q)
            series.push_back({inf.faults[static_cast<std::size_t>(members[q])].name + " at " + inf.horizons[h].name,
                              q % 4 == 0 ? kFault : q % 4 == 1 ? rgb(60, 90, 160) : q % 4 == 2 ? rgb(40, 130, 70) : rgb(150, 90, 20),
                              false, per[q]});
        }
      }
      if (series.empty()) continue;
      auto fgr = strike_throw_figure(c, series, name);
      out.push_back(write_figure(fgr, folder, figure_basename(c, "strike-throw", name), c.model_name + " - throw along strike, " + name,
                                 "Throw along strike of the linked faults " + name + ", summed and per fault.", notes_common));
    }
  }
  // Combined deck of all figures.
  if (raster_available()) {
    std::vector<Slide> slides;
    for (auto& w : out) {
      std::string png = read_file(folder / (w.basename + ".png"));
      auto u32 = [&](std::size_t o) { return (static_cast<unsigned char>(png[o]) << 24) | (static_cast<unsigned char>(png[o + 1]) << 16) |
                                             (static_cast<unsigned char>(png[o + 2]) << 8) | static_cast<unsigned char>(png[o + 3]); };
      slides.push_back({w.title, w.caption, notes_common, png, static_cast<int>(u32(16)), static_cast<int>(u32(20))});
    }
    write_file(folder / (figure_basename(c, "figure-set", "deck") + ".pptx"), make_pptx(slides, c.model_name + " figure set"));
  }
  return out;
}

// =========================================================================================
// Grids (stage 3B, STX-TN-12)
// =========================================================================================
namespace {
struct PropScale {
  double lo = 0, hi = 1;
  bool zone = true;
};
PropScale prop_scale(const repr::CellProperty* p, const repr::CornerPointGrid& g) {
  PropScale s;
  if (!p) return s;
  s.zone = false;
  std::vector<double> v;
  for (std::size_t c = 0; c < p->v.size(); ++c)
    if (g.actnum[c] && std::isfinite(p->v[c])) v.push_back(p->v[c]);
  if (v.empty()) return s;
  std::sort(v.begin(), v.end());
  s.lo = v[static_cast<std::size_t>(0.02 * static_cast<double>(v.size() - 1))];
  s.hi = v[static_cast<std::size_t>(0.98 * static_cast<double>(v.size() - 1))];
  if (!(s.hi > s.lo)) { s.lo = v.front(); s.hi = v.back() > v.front() ? v.back() : v.front() + 1; }
  return s;
}
Color cell_color(const FigureContext& c, const repr::CornerPointGrid& g, const repr::CellProperty* p, const PropScale& s, std::size_t id) {
  if (!p) return unit_color(c, static_cast<std::size_t>(std::max(0, g.zone[id])));
  double v = p->v[id];
  if (!std::isfinite(v)) return rgb(200, 200, 200);
  return viridis_like((v - s.lo) / (s.hi - s.lo));
}
void prop_legend(Figure& f, const FigureContext& c, const repr::CornerPointGrid& g, const repr::CellProperty* p, const PropScale& s, double rx,
                 double ry, double maxy) {
  if (!p) {
    std::vector<LegendEntry> es;
    std::vector<int> seen;
    for (int u : g.layer_unit)
      if (std::find(seen.begin(), seen.end(), u) == seen.end()) {
        seen.push_back(u);
        es.push_back({"Zones", c.model->info().units[static_cast<std::size_t>(u)].name, unit_color(c, static_cast<std::size_t>(u))});
      }
    es.push_back({"Faults", "Fault pillars", none(), true, kFault, 1.2});
    double ly = legend(f, rx, ry, 180, es, maxy) + 8;
    key_box(f, c, rx, ly, 180);
    return;
  }
  f.text({rx, ry}, p->name + " (" + p->unit + ")", 7.5, Anchor::Start, true, kInk);
  for (int k = 0; k < 100; ++k) f.rect(rx, ry + 6 + k * 1.6, 14, 1.7, viridis_like(1 - k / 99.0), none());
  f.rect(rx, ry + 6, 14, 160, none(), kInk, 0.4);
  f.text({rx + 18, ry + 12}, fmt_shortest(std::round(s.hi * 1000) / 1000), 6.5, Anchor::Start, false, kInk);
  f.text({rx + 18, ry + 166}, fmt_shortest(std::round(s.lo * 1000) / 1000), 6.5, Anchor::Start, false, kInk);
  f.rect(rx, ry + 176, 14, 8, rgb(200, 200, 200), kInk, 0.4);
  f.text({rx + 18, ry + 183}, "missing (" + std::to_string(p->missing) + " cells)", 6.5, Anchor::Start, false, kInk);
  // the evidence class, wrapped at about 40 characters
  double ly = ry + 198;
  std::string rest = "Evidence: " + p->evidence;
  while (!rest.empty()) {
    std::size_t cut = rest.size() <= 40 ? rest.size() : rest.rfind(' ', 40);
    if (cut == std::string::npos || cut == 0) cut = std::min<std::size_t>(40, rest.size());
    f.text({rx, ly}, rest.substr(0, cut), 6.5, Anchor::Start, false, kInk);
    rest = cut < rest.size() ? rest.substr(cut + (rest[cut] == ' ' ? 1 : 0)) : "";
    ly += 9;
  }
  ly += 7;
  key_box(f, c, rx, ly, 180);
}
}  // namespace

Figure grid_section_figure(const FigureContext& c, const repr::CornerPointGrid& g, const repr::CellProperty* prop, bool along_x, int line) {
  Figure f(842, 595);
  const auto& inf = c.model->info();
  const auto sc = prop_scale(prop, g);
  // cell polygons in (distance along the line, z), averaged over the two sides of the row
  struct Poly { std::vector<std::pair<double, double>> p; Color col; };
  std::vector<Poly> polys;
  double zmin = 1e300, zmax = -1e300;
  const int n = along_x ? g.nx : g.ny;
  const double a0 = along_x ? g.x0 : g.y0;
  for (int k = 0; k < g.nz; ++k)
    for (int s = 0; s < n; ++s) {
      int i = along_x ? s : line, j = along_x ? line : s;
      auto id = g.cell(i, j, k);
      if (!g.actnum[id]) continue;
      auto cc = g.corners(i, j, k);
      auto mid = [&](int a, int b) { Vec3 m = (cc[a] + cc[b]) / 2.0; return std::make_pair((along_x ? m.x : m.y) - a0, m.z); };
      // corners: lower along-line side first; pairs across the row
      Poly pl;
      if (along_x) pl.p = {mid(0, 2), mid(1, 3), mid(5, 7), mid(4, 6)};
      else pl.p = {mid(0, 1), mid(2, 3), mid(6, 7), mid(4, 5)};
      for (auto& q : pl.p) { zmin = std::min(zmin, q.second); zmax = std::max(zmax, q.second); }
      pl.col = cell_color(c, g, prop, sc, id);
      polys.push_back(std::move(pl));
    }
  const double len = along_x ? g.nx * g.dx : g.ny * g.dy;
  if (zmax < zmin) { zmin = inf.domain.z_base; zmax = inf.domain.z_top; }
  const double px = 60, py = 70, pw = 842 - 280, ph = 595 - 150;
  const double ve = c.vertical_exaggeration;
  // the description's exaggeration, raised where needed so the section fills at least 40 % of the height
  double sx = pw / len, sz = std::min(ph / (zmax - zmin), sx * std::max(ve, 0.4 * ph / ((zmax - zmin) * sx)));
  auto P = [&](double a, double z) { return Pt{px + a * sx, py + (zmax - z) * sz}; };
  const std::string where = along_x ? "row j = " + std::to_string(line + 1) + " (northing " + fmt_fixed(g.y0 + (line + 0.5) * g.dy, 0) + " m)"
                                    : "column i = " + std::to_string(line + 1) + " (easting " + fmt_fixed(g.x0 + (line + 0.5) * g.dx, 0) + " m)";
  title(f, (c.title.empty() ? c.model_name : c.title) + " - grid section, " + where,
        std::string("Active cells of the corner-point grid coloured by ") + (prop ? prop->name : std::string("zone")) +
            "; fault pillars in red; vertical exaggeration " + fmt_fixed(sz / sx, 1) + "x");
  f.clip(px, py, len * sx, (zmax - zmin) * sz);
  for (auto& pl : polys) {
    std::vector<Pt> q;
    for (auto& [a, z] : pl.p) q.push_back(P(a, z));
    f.path(q, true, pl.col, rgb(90, 90, 90), 0.15);
  }
  for (int s = 0; s <= n; ++s) {
    int i = along_x ? s : line, j = along_x ? line : s;
    for (int side = 0; side < 2; ++side) {
      int ii = along_x ? i : i + side, jj = along_x ? j + side : j;
      auto p = g.pillar(ii, jj);
      if (g.pillar_fault[p] < 0) continue;
      const double* cd = &g.coord[6 * p];
      f.line(P((along_x ? cd[0] : cd[1]) - a0, cd[2]), P((along_x ? cd[3] : cd[4]) - a0, cd[5]), kFault, 1.0);
    }
  }
  f.unclip();
  f.rect(px, py, len * sx, (zmax - zmin) * sz, none(), kInk, 0.6);
  for (double t : nice_ticks(0, len, 7)) f.text({P(t, zmin).x, py + (zmax - zmin) * sz + 11}, km_or_m(t), 6.5, Anchor::Middle, false, kInk);
  for (double t : nice_ticks(-zmax, -zmin, 6)) f.text({px - 4, P(0, -t).y + 2.5}, km_or_m(t), 6.5, Anchor::End, false, kInk);
  f.text({px + len * sx / 2, py + (zmax - zmin) * sz + 26}, std::string("Distance along the ") + (along_x ? "row (m, easting from the grid origin)" : "column (m, northing from the grid origin)"),
         7.5, Anchor::Middle, false, kInk);
  f.text({px - 40, py + (zmax - zmin) * sz / 2}, "Depth below datum (m)", 7.5, Anchor::Middle, false, kInk, -90);
  prop_legend(f, c, g, prop, sc, 842 - 200, 60, 595.0 - 30 - key_box_height(c) - 8);
  footer(f, c);
  return f;
}

Figure grid_layer_figure(const FigureContext& c, const repr::CornerPointGrid& g, const repr::CellProperty* prop, int k) {
  Figure f(842, 595);
  const auto& inf = c.model->info();
  const auto& d = inf.domain;
  const auto sc = prop_scale(prop, g);
  double side = std::min(595.0 - 110, 842.0 - 280);
  double s = side / std::max(d.sx, d.sy);
  double mx0 = 60, my0 = 60, mw = d.sx * s, mh = d.sy * s;
  auto P = [&](double x, double y) { return Pt{mx0 + (x - d.x0) * s, my0 + mh - (y - d.y0) * s}; };
  const std::string zn = inf.units[static_cast<std::size_t>(g.layer_unit[static_cast<std::size_t>(k)])].name;
  title(f, (c.title.empty() ? c.model_name : c.title) + " - grid layer " + std::to_string(k + 1) + " (" + zn + ")",
        std::string("Cells of layer ") + std::to_string(k + 1) + " in map view coloured by " + (prop ? prop->name : std::string("zone")) +
            "; inactive cells blank; fault pillars in red");
  f.clip(mx0, my0, mw, mh);
  for (int j = 0; j < g.ny; ++j)
    for (int i = 0; i < g.nx; ++i) {
      auto id = g.cell(i, j, k);
      if (!g.actnum[id]) continue;
      auto cc = g.corners(i, j, k);
      std::vector<Pt> q;
      for (int a : {0, 1, 3, 2}) { Vec3 m = (cc[a] + cc[a + 4]) / 2.0; q.push_back(P(m.x, m.y)); }
      f.path(q, true, cell_color(c, g, prop, sc, id), cell_color(c, g, prop, sc, id), 0.2);
    }
  for (int j = 0; j <= g.ny; ++j)
    for (int i = 0; i <= g.nx; ++i) {
      auto p = g.pillar(i, j);
      if (g.pillar_fault[p] < 0) continue;
      for (int dir = 0; dir < 2; ++dir) {
        int i2 = dir == 0 ? i : i + 1, j2 = dir == 0 ? j + 1 : j;
        if (i2 > g.nx || j2 > g.ny) continue;
        auto p2 = g.pillar(i2, j2);
        if (g.pillar_fault[p2] != g.pillar_fault[p]) continue;
        const double *a = &g.coord[6 * p], *b = &g.coord[6 * p2];
        f.line(P(0.5 * (a[0] + a[3]), 0.5 * (a[1] + a[4])), P(0.5 * (b[0] + b[3]), 0.5 * (b[1] + b[4])), kFault, 1.2);
      }
    }
  f.unclip();
  f.rect(mx0, my0, mw, mh, none(), kInk, 0.6);
  for (double t : nice_ticks(d.x0, d.x0 + d.sx, 6)) f.text({P(t, d.y0).x, my0 + mh + 11}, km_or_m(t), 6.5, Anchor::Middle, false, kInk);
  for (double t : nice_ticks(d.y0, d.y0 + d.sy, 6)) f.text({mx0 - 5, P(d.x0, t).y + 2.5}, km_or_m(t), 6.5, Anchor::End, false, kInk);
  f.text({mx0 + mw / 2, my0 + mh + 24}, "Easting (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({mx0 - 40, my0 + mh / 2}, "Northing (m)", 7.5, Anchor::Middle, false, kInk, -90);
  prop_legend(f, c, g, prop, sc, 842 - 200, 60, 595.0 - 30 - key_box_height(c) - 8);
  footer(f, c);
  return f;
}

Figure grid_view3d_figure(const FigureContext& c, const repr::CornerPointGrid& g, double az, double el) {
  Figure f(842, 560);
  const auto& inf = c.model->info();
  std::map<int, TriMesh> by_zone;
  auto active = [&](int i, int j, int k) { return i >= 0 && j >= 0 && k >= 0 && i < g.nx && j < g.ny && k < g.nz && g.actnum[g.cell(i, j, k)]; };
  auto quad = [&](TriMesh& m, const std::array<Vec3, 8>& cc, int a, int b, int d2, int e) {
    auto base = static_cast<std::uint32_t>(m.v.size());
    for (int q : {a, b, d2, e}) m.v.push_back(cc[q]);
    m.t.push_back({base, base + 1, base + 2});
    m.t.push_back({base, base + 2, base + 3});
  };
  for (int k = 0; k < g.nz; ++k)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i) {
        if (!active(i, j, k)) continue;
        auto cc = g.corners(i, j, k);
        auto& m = by_zone[g.zone[g.cell(i, j, k)]];
        if (!active(i, j, k - 1)) quad(m, cc, 0, 1, 3, 2);
        if (!active(i, j, k + 1)) quad(m, cc, 4, 6, 7, 5);
        // sides at the grid's edge, next to inactive cells and across faults (block changes)
        const BlockCode b = g.block[g.cell(i, j, k)];
        auto open = [&](int a, int bb) { return !active(a, bb, k) || g.block[g.cell(a, bb, k)] != b; };
        if (open(i - 1, j)) quad(m, cc, 0, 2, 6, 4);
        if (open(i + 1, j)) quad(m, cc, 1, 5, 7, 3);
        if (open(i, j - 1)) quad(m, cc, 0, 4, 5, 1);
        if (open(i, j + 1)) quad(m, cc, 2, 3, 7, 6);
      }
  std::vector<Mesh3> meshes;
  for (auto& [u, m] : by_zone) meshes.push_back({m, unit_color(c, static_cast<std::size_t>(std::max(0, u))), 1.0});
  BBox box = inf.domain.box();
  Image im = render_meshes(meshes, az, el, c.vertical_exaggeration, 1800, 1000, box);
  title(f, (c.title.empty() ? c.model_name : c.title) + " - grid, 3D view",
        "Outer faces of the active cells and their faces across faults, coloured by zone; viewed from azimuth " + fmt_fixed(az, 0) + " deg, elevation " + fmt_fixed(el, 0) +
            " deg; vertical exaggeration " + fmt_fixed(c.vertical_exaggeration, 1) + "x");
  const Image* img = f.own(std::move(im));
  double iw = 842 - 250, ih = iw * 1000.0 / 1800.0;
  f.image(20, 70, iw, ih, img);
  std::vector<LegendEntry> es;
  for (auto& [u, m] : by_zone) es.push_back({"Zones", inf.units[static_cast<std::size_t>(std::max(0, u))].name, unit_color(c, static_cast<std::size_t>(std::max(0, u)))});
  double ly = legend(f, 842 - 214, 56, 196, es, 560.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 214, ly, 196);
  footer(f, c);
  return f;
}

// =========================================================================================
// Ensembles (stage 3C, STX-TN-13)
// =========================================================================================
namespace {
std::string ensemble_note(const EnsembleLabel& l) {
  return std::to_string(l.members_ok) + " members ok" + (l.members_failed ? ", " + std::to_string(l.members_failed) + " failed (listed in members.csv)" : std::string()) +
         "; seed " + std::to_string(l.seed);
}
std::string short_num(double v) {
  const double a = std::fabs(v);
  if (a < 1e-12) return "0";
  return fmt_fixed(v, a >= 100 ? 1 : a >= 1 ? 2 : 4);
}
}  // namespace

Figure tornado_figure(const FigureContext& c, const std::string& key, double base, const std::vector<TornadoBar>& bars_in, const EnsembleLabel& l) {
  Figure f(842, 520);
  auto bars = bars_in;
  std::stable_sort(bars.begin(), bars.end(), [](auto& a, auto& b) { return std::fabs(a.high - a.low) > std::fabs(b.high - b.low); });
  double lo = base, hi = base;
  for (auto& b : bars) { lo = std::min({lo, b.low, b.high}); hi = std::max({hi, b.low, b.high}); }
  if (!(hi > lo)) { lo -= 1; hi += 1; }
  const double pad = 0.08 * (hi - lo);
  lo -= pad;
  hi += pad;
  title(f, (c.title.empty() ? c.model_name : c.title) + " - sensitivity of " + key,
        "One at a time: each parameter at its P10 (light) and P90 (dark), the others at their base values; base case line; " + ensemble_note(l));
  const double px = 230, py = 70, pw = 842 - 230 - 230, rowh = std::min(34.0, (520.0 - 150) / std::max<std::size_t>(1, bars.size()));
  const double ph = rowh * static_cast<double>(bars.size());
  auto X = [&](double v) { return px + (v - lo) / (hi - lo) * pw; };
  for (double t : nice_ticks(lo, hi, 6)) {
    f.line({X(t), py}, {X(t), py + ph}, kGrid, 0.4);
    f.text({X(t), py + ph + 11}, short_num(t), 7, Anchor::Middle, false, kInk);
  }
  for (std::size_t k = 0; k < bars.size(); ++k) {
    const double y = py + rowh * static_cast<double>(k) + 0.15 * rowh, h = 0.7 * rowh;
    const auto& b = bars[k];
    auto seg = [&](double v, Color col) {
      double a = X(std::min(base, v)), e = X(std::max(base, v));
      if (e - a < 0.8) { a -= 0.4; e += 0.4; }
      f.rect(a, y, e - a, h, col, kInk, 0.3);
    };
    seg(b.low, rgb(170, 200, 230));
    seg(b.high, rgb(40, 90, 160));
    f.text({px - 6, y + h / 2 + 2.5}, b.parameter, 7, Anchor::End, false, kInk);
    f.text({X(std::max({b.low, b.high, base})) + 4, y + h / 2 + 2.5}, "swing " + short_num(std::fabs(b.high - b.low)), 6.5, Anchor::Start, false, kInk);
  }
  f.line({X(base), py - 4}, {X(base), py + ph + 4}, kFault, 1.0);
  f.text({X(base), py - 7}, "base " + short_num(base), 7, Anchor::Middle, true, kFault);
  f.text({px + pw / 2, py + ph + 26}, key, 8, Anchor::Middle, false, kInk);
  std::vector<LegendEntry> es{{"Parameter at", "P10", rgb(170, 200, 230)}, {"Parameter at", "P90", rgb(40, 90, 160)}};
  double ly = legend(f, 842 - 214, 56, 196, es, 520.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 214, ly, 196);
  footer(f, c);
  return f;
}

Figure histogram_figure(const FigureContext& c, const std::string& key, const std::vector<double>& values, double p10, double p50, double p90,
                        const EnsembleLabel& l) {
  Figure f(842, 520);
  std::vector<double> v;
  for (double x : values)
    if (std::isfinite(x)) v.push_back(x);
  double lo = v.empty() ? 0 : *std::min_element(v.begin(), v.end()), hi = v.empty() ? 1 : *std::max_element(v.begin(), v.end());
  if (!(hi > lo)) { lo -= 0.5; hi += 0.5; }
  const int nb = std::clamp(static_cast<int>(std::lround(std::sqrt(static_cast<double>(v.size())))), 5, 40);
  std::vector<int> cnt(static_cast<std::size_t>(nb), 0);
  for (double x : v) cnt[static_cast<std::size_t>(std::clamp(static_cast<int>((x - lo) / (hi - lo) * nb), 0, nb - 1))]++;
  const int cmax = std::max(1, *std::max_element(cnt.begin(), cnt.end()));
  title(f, (c.title.empty() ? c.model_name : c.title) + " - distribution of " + key,
        "Histogram over the ok members; P10, P50 and P90 marked; " + ensemble_note(l));
  const double px = 70, py = 60, pw = 842 - 70 - 250, ph = 520 - 60 - 80;
  auto X = [&](double t) { return px + (t - lo) / (hi - lo) * pw; };
  auto Y = [&](double n) { return py + ph - n / (cmax * 1.1) * ph; };
  for (double t : nice_ticks(lo, hi, 7)) {
    f.line({X(t), py + ph}, {X(t), py + ph + 3}, kInk, 0.5);
    f.text({X(t), py + ph + 12}, short_num(t), 7, Anchor::Middle, false, kInk);
  }
  for (double n : nice_ticks(0, cmax * 1.1, 5)) {
    f.line({px, Y(n)}, {px + pw, Y(n)}, kGrid, 0.4);
    f.text({px - 4, Y(n) + 2.5}, fmt_fixed(n, 0), 7, Anchor::End, false, kInk);
  }
  for (int b = 0; b < nb; ++b) {
    const double a = lo + (hi - lo) * b / nb, e = lo + (hi - lo) * (b + 1) / nb;
    f.rect(X(a), Y(cnt[static_cast<std::size_t>(b)]), X(e) - X(a), py + ph - Y(cnt[static_cast<std::size_t>(b)]), rgb(120, 160, 200), rgb(60, 90, 130), 0.4);
  }
  for (auto [val, lab] : {std::pair{p10, "P10"}, {p50, "P50"}, {p90, "P90"}}) {
    if (!std::isfinite(val)) continue;
    f.line({X(val), py}, {X(val), py + ph}, kFault, 1.0, {4, 2});
    f.text({X(val), py - 4}, std::string(lab) + " " + short_num(val), 7, Anchor::Middle, true, kFault);
  }
  f.rect(px, py, pw, ph, none(), kInk, 0.6);
  f.text({px + pw / 2, py + ph + 26}, key, 8, Anchor::Middle, false, kInk);
  f.text({px - 40, py + ph / 2}, "Members", 8, Anchor::Middle, false, kInk, -90);
  std::vector<LegendEntry> es{{"Ensemble", "Members per bin", rgb(120, 160, 200)}, {"Percentiles", "P10, P50, P90", none(), true, kFault, 1.0}};
  double ly = legend(f, 842 - 214, 56, 196, es, 520.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 214, ly, 196);
  footer(f, c);
  return f;
}

Figure sobol_figure(const FigureContext& c, const std::string& key, const std::vector<std::pair<std::string, double>>& idx_in, double threshold,
                    const EnsembleLabel& l) {
  Figure f(842, 520);
  auto idx = idx_in;
  std::stable_sort(idx.begin(), idx.end(), [](auto& a, auto& b) { return a.second > b.second; });
  title(f, (c.title.empty() ? c.model_name : c.title) + " - first-order indices of " + key,
        "Share of the variance of " + key + " explained by each parameter alone (given-data estimate); below the dashed line not distinguishable from 0; " +
            ensemble_note(l));
  const double px = 230, py = 70, pw = 842 - 230 - 230, rowh = std::min(30.0, (520.0 - 150) / std::max<std::size_t>(1, idx.size()));
  const double ph = rowh * static_cast<double>(idx.size());
  auto X = [&](double v) { return px + std::clamp(v, 0.0, 1.0) * pw; };
  for (double t : {0.0, 0.25, 0.5, 0.75, 1.0}) {
    f.line({X(t), py}, {X(t), py + ph}, kGrid, 0.4);
    f.text({X(t), py + ph + 11}, fmt_fixed(t, 2), 7, Anchor::Middle, false, kInk);
  }
  for (std::size_t k = 0; k < idx.size(); ++k) {
    const double y = py + rowh * static_cast<double>(k) + 0.15 * rowh, h = 0.7 * rowh;
    f.rect(px, y, std::max(0.5, X(idx[k].second) - px), h, idx[k].second >= threshold ? rgb(40, 90, 160) : rgb(190, 190, 190), kInk, 0.3);
    f.text({px - 6, y + h / 2 + 2.5}, idx[k].first, 7, Anchor::End, false, kInk);
    f.text({X(idx[k].second) + 4, y + h / 2 + 2.5}, fmt_fixed(idx[k].second, 3), 6.5, Anchor::Start, false, kInk);
  }
  f.line({X(threshold), py - 4}, {X(threshold), py + ph + 4}, kInk, 0.8, {3, 2});
  f.text({px + pw / 2, py + ph + 26}, "First-order index S_i", 8, Anchor::Middle, false, kInk);
  std::vector<LegendEntry> es{{"Index", "Distinguishable from 0", rgb(40, 90, 160)}, {"Index", "Not distinguishable (< 2/sqrt(N))", rgb(190, 190, 190)}};
  double ly = legend(f, 842 - 214, 56, 196, es, 520.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 214, ly, 196);
  footer(f, c);
  return f;
}

Figure probability_map_figure(const FigureContext& c, const std::string& what, const Grid2D& p, const EnsembleLabel& l) {
  Figure f(842, 595);
  const auto& d = c.model->info().domain;
  double side = std::min(595.0 - 110, 842.0 - 280);
  double sc = side / std::max(d.sx, d.sy);
  double mx0 = 60, my0 = 60, mw = d.sx * sc, mh = d.sy * sc;
  auto P = [&](double x, double y) { return Pt{mx0 + (x - d.x0) * sc, my0 + mh - (y - d.y0) * sc}; };
  title(f, (c.title.empty() ? c.model_name : c.title) + " - probability: " + what, "Fraction of the ok members; " + ensemble_note(l));
  const auto& g = p.spec;
  f.clip(mx0, my0, mw, mh);
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i) {
      const double v = p.v[g.index(i, j)];
      Pt a = P(g.x(i) - g.dx / 2, g.y(j) + g.dy / 2), b = P(g.x(i) + g.dx / 2, g.y(j) - g.dy / 2);
      Color col = std::isfinite(v) ? viridis_like(v) : rgb(220, 220, 220);
      f.rect(a.x, a.y, b.x - a.x, b.y - a.y, col, col, 0.2);
    }
  for (auto& w : c.wells) {
    Pt q = P(w.x(), w.y());
    f.path({{q.x - 2.5, q.y}, {q.x, q.y - 2.5}, {q.x + 2.5, q.y}, {q.x, q.y + 2.5}}, true, rgb(255, 255, 255), kWell, 0.9);
    f.text({q.x + 4, q.y - 3}, w.name(), 7, Anchor::Start, true, kWell);
  }
  f.unclip();
  f.rect(mx0, my0, mw, mh, none(), kInk, 0.6);
  for (double t : nice_ticks(d.x0, d.x0 + d.sx, 6)) f.text({P(t, d.y0).x, my0 + mh + 11}, km_or_m(t), 6.5, Anchor::Middle, false, kInk);
  for (double t : nice_ticks(d.y0, d.y0 + d.sy, 6)) f.text({mx0 - 5, P(d.x0, t).y + 2.5}, km_or_m(t), 6.5, Anchor::End, false, kInk);
  f.text({mx0 + mw / 2, my0 + mh + 24}, "Easting (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({mx0 - 40, my0 + mh / 2}, "Northing (m)", 7.5, Anchor::Middle, false, kInk, -90);
  double rx = 842 - 200, ry = 60;
  f.text({rx, ry}, "Probability", 7.5, Anchor::Start, true, kInk);
  for (int k = 0; k < 100; ++k) f.rect(rx, ry + 6 + k * 1.6, 14, 1.7, viridis_like(1 - k / 99.0), none());
  f.rect(rx, ry + 6, 14, 160, none(), kInk, 0.4);
  f.text({rx + 18, ry + 12}, "1", 6.5, Anchor::Start, false, kInk);
  f.text({rx + 18, ry + 88}, "0.5", 6.5, Anchor::Start, false, kInk);
  f.text({rx + 18, ry + 166}, "0", 6.5, Anchor::Start, false, kInk);
  key_box(f, c, rx, ry + 186, 180);
  footer(f, c);
  return f;
}

Figure scalar_map_figure(const FigureContext& c, const std::string& ttl, const std::string& subtitle, const std::string& bar_label, const Grid2D& p,
                         double vmin, double vmax) {
  Figure f(842, 595);
  const auto& d = c.model->info().domain;
  double side = std::min(595.0 - 110, 842.0 - 280);
  double sc = side / std::max(d.sx, d.sy);
  double mx0 = 60, my0 = 60, mw = d.sx * sc, mh = d.sy * sc;
  auto P = [&](double x, double y) { return Pt{mx0 + (x - d.x0) * sc, my0 + mh - (y - d.y0) * sc}; };
  title(f, (c.title.empty() ? c.model_name : c.title) + " - " + ttl, subtitle);
  if (!(vmax > vmin)) vmax = vmin + 1;
  const auto& g = p.spec;
  f.clip(mx0, my0, mw, mh);
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i) {
      const double v = p.v[g.index(i, j)];
      Pt a = P(g.x(i) - g.dx / 2, g.y(j) + g.dy / 2), b = P(g.x(i) + g.dx / 2, g.y(j) - g.dy / 2);
      Color col = std::isfinite(v) ? viridis_like(std::clamp((v - vmin) / (vmax - vmin), 0.0, 1.0)) : rgb(220, 220, 220);
      f.rect(a.x, a.y, b.x - a.x, b.y - a.y, col, col, 0.2);
    }
  for (auto& w : c.wells) {
    Pt q = P(w.x(), w.y());
    f.path({{q.x - 2.5, q.y}, {q.x, q.y - 2.5}, {q.x + 2.5, q.y}, {q.x, q.y + 2.5}}, true, rgb(255, 255, 255), kWell, 0.9);
    f.text({q.x + 4, q.y - 3}, w.name(), 7, Anchor::Start, true, kWell);
  }
  f.unclip();
  f.rect(mx0, my0, mw, mh, none(), kInk, 0.6);
  for (double t : nice_ticks(d.x0, d.x0 + d.sx, 6)) f.text({P(t, d.y0).x, my0 + mh + 11}, km_or_m(t), 6.5, Anchor::Middle, false, kInk);
  for (double t : nice_ticks(d.y0, d.y0 + d.sy, 6)) f.text({mx0 - 5, P(d.x0, t).y + 2.5}, km_or_m(t), 6.5, Anchor::End, false, kInk);
  f.text({mx0 + mw / 2, my0 + mh + 24}, "Easting (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({mx0 - 40, my0 + mh / 2}, "Northing (m)", 7.5, Anchor::Middle, false, kInk, -90);
  double rx = 842 - 200, ry = 60;
  f.text({rx, ry}, bar_label, 7.5, Anchor::Start, true, kInk);
  for (int k = 0; k < 100; ++k) f.rect(rx, ry + 6 + k * 1.6, 14, 1.7, viridis_like(1 - k / 99.0), none());
  f.rect(rx, ry + 6, 14, 160, none(), kInk, 0.4);
  // three significant figures on the colour bar
  auto sig3 = [](double v) {
    if (v == 0 || !std::isfinite(v)) return fmt_fixed(v, 0);
    const int dec = std::clamp(2 - static_cast<int>(std::floor(std::log10(std::fabs(v)))), 0, 8);
    return fmt_fixed(v, dec);
  };
  f.text({rx + 18, ry + 12}, sig3(vmax), 6.5, Anchor::Start, false, kInk);
  f.text({rx + 18, ry + 88}, sig3(0.5 * (vmin + vmax)), 6.5, Anchor::Start, false, kInk);
  f.text({rx + 18, ry + 166}, sig3(vmin), 6.5, Anchor::Start, false, kInk);
  key_box(f, c, rx, ry + 186, 180);
  footer(f, c);
  return f;
}

}  // namespace stratonix::fig

// ---- stage 4C: well log and correlation panels (STX-TN-16 section 1.6) ---------------------
namespace stratonix::fig {

namespace {
void log_track(Figure& f, double x, double y, double w, double h, const WellPanel::Track& t, const std::vector<double>& idx, double d0, double d1,
               Color col) {
  f.rect(x, y, w, h, none(), kInk, 0.5);
  double vmin = 1e300, vmax = -1e300;
  for (double v : t.v)
    if (std::isfinite(v)) { vmin = std::min(vmin, v); vmax = std::max(vmax, v); }
  f.text({x + w / 2, y - 14}, t.name, 7, Anchor::Middle, true, kInk);
  f.text({x + w / 2, y - 5}, t.unit.empty() ? "(no unit)" : t.unit, 6, Anchor::Middle, false, kInk);
  if (!(vmax > vmin)) { vmax = vmin + 1; }
  f.text({x, y + h + 9}, fmt_shortest(vmin), 5.5, Anchor::Start, false, kInk);
  f.text({x + w, y + h + 9}, fmt_shortest(vmax), 5.5, Anchor::End, false, kInk);
  std::vector<Pt> seg;
  auto flush = [&]() { if (seg.size() >= 2) f.path(seg, false, none(), col, 0.7); seg.clear(); };
  for (std::size_t i = 0; i < idx.size() && i < t.v.size(); ++i) {
    if (!std::isfinite(t.v[i]) || !std::isfinite(idx[i]) || idx[i] < d0 || idx[i] > d1) { flush(); continue; }
    seg.push_back({x + (t.v[i] - vmin) / (vmax - vmin) * w, y + (idx[i] - d0) / (d1 - d0) * h});
  }
  flush();
}
Color track_color(std::size_t k) {
  static const Color cs[] = {rgb(31, 119, 180), rgb(214, 39, 40), rgb(44, 160, 44), rgb(148, 103, 189), rgb(140, 86, 75), rgb(23, 190, 207)};
  return cs[k % 6];
}
}  // namespace

Figure well_log_figure(const FigureContext& c, const WellPanel& p) {
  Figure f(842, 595);
  title(f, (c.title.empty() ? c.model_name : c.title) + " - logs of " + p.well,
        "Curves as read from the well's log files (units and nulls from the files), with tops, age data and the model's units along the well");
  double d0 = 1e300, d1 = -1e300;
  for (double v : p.index)
    if (std::isfinite(v)) { d0 = std::min(d0, v); d1 = std::max(d1, v); }
  for (auto& t : p.tops) { d0 = std::min(d0, t.second); d1 = std::max(d1, t.second); }
  if (!(d1 > d0)) { d0 = 0; d1 = 1; }
  const double pad = 0.02 * (d1 - d0);
  d0 -= pad;
  d1 += pad;
  const double py = 95, ph = 595 - 95 - 60, px = 70;
  const std::size_t nt = std::min<std::size_t>(p.tracks.size(), 6);
  const double tw = nt ? std::min(110.0, (842 - 230 - px - 40) / static_cast<double>(nt) - 8) : 0;
  auto Y = [&](double d) { return py + (d - d0) / (d1 - d0) * ph; };
  // depth axis and the unit column
  f.rect(px, py, 18, ph, none(), kInk, 0.5);
  for (auto& u : p.units) {
    const double a = std::max(d0, u.top), b = std::min(d1, u.base);
    if (b <= a) continue;
    f.rect(px, Y(a), 18, Y(b) - Y(a), unit_color(c, static_cast<std::size_t>(std::max(0, u.unit))), none());
  }
  for (double t : nice_ticks(d0, d1, 8)) {
    f.line({px - 3, Y(t)}, {px, Y(t)}, kInk, 0.5);
    f.text({px - 5, Y(t) + 2.5}, fmt_fixed(t, 0), 6.5, Anchor::End, false, kInk);
  }
  f.text({px - 45, py + ph / 2}, p.index_label, 7.5, Anchor::Middle, false, kInk, -90);
  f.text({px + 9, py - 5}, "Units", 6, Anchor::Middle, false, kInk);
  double x = px + 30;
  for (std::size_t k = 0; k < nt; ++k) {
    log_track(f, x, py, tw, ph, p.tracks[k], p.index, d0, d1, track_color(k));
    x += tw + 8;
  }
  const double xr = x - 8;
  for (auto& t : p.tops) {
    f.line({px, Y(t.second)}, {xr, Y(t.second)}, rgb(200, 30, 30), 0.8, {4, 2});
    f.text({xr + 3, Y(t.second) + 2.5}, t.first, 6.5, Anchor::Start, true, rgb(200, 30, 30));
  }
  for (auto& a : p.ages) {
    f.path({{px + 18, Y(a.depth) - 3}, {px + 24, Y(a.depth)}, {px + 18, Y(a.depth) + 3}}, true, rgb(60, 60, 60), none());
    f.text({px + 25, Y(a.depth) - 3}, a.name + " " + fmt_shortest(a.age) + " Ma", 5.5, Anchor::Start, false, rgb(60, 60, 60));
  }
  std::vector<LegendEntry> es;
  es.push_back({"Markers", "Top (dashed)", none(), true, rgb(200, 30, 30), 0.8});
  if (!p.ages.empty()) es.push_back({"Markers", "Age datum (triangle)", rgb(60, 60, 60), false, none(), 1, 2});
  std::vector<int> seen;
  for (auto& u : p.units)
    if (u.unit >= 0 && std::find(seen.begin(), seen.end(), u.unit) == seen.end()) {
      seen.push_back(u.unit);
      const auto& inf = c.model->info();
      es.push_back({"Model units", static_cast<std::size_t>(u.unit) < inf.units.size() ? inf.units[static_cast<std::size_t>(u.unit)].name : "unit",
                    unit_color(c, static_cast<std::size_t>(u.unit))});
    }
  double ly = legend(f, 842 - 200, py, 190, es, 595.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 200, ly, 190);
  footer(f, c);
  return f;
}

Figure correlation_figure(const FigureContext& c, const CorrelationPanel& p) {
  Figure f(842, 595);
  title(f, (c.title.empty() ? c.model_name : c.title) + " - well correlation",
        "Wells by TVD below the model datum; each line labelled with the ages and rule, the tops or the model it rests on" +
            (p.note.empty() ? std::string() : " (" + p.note + ")"));
  double d0 = 1e300, d1 = -1e300;
  for (auto& w : p.wells)
    for (double v : w.index)
      if (std::isfinite(v)) { d0 = std::min(d0, v); d1 = std::max(d1, v); }
  for (auto& row : p.depth)
    for (double v : row)
      if (std::isfinite(v)) { d0 = std::min(d0, v); d1 = std::max(d1, v); }
  if (!(d1 > d0)) { d0 = 0; d1 = 1; }
  const double pad = 0.03 * (d1 - d0);
  d0 -= pad;
  d1 += pad;
  const double py = 110, ph = 595 - 110 - 60, px = 75, pw = 842 - 230 - px;
  auto Y = [&](double d) { return py + (d - d0) / (d1 - d0) * ph; };
  const std::size_t nw = p.wells.size();
  const double slot = nw ? pw / static_cast<double>(nw) : pw, tw = std::min(60.0, slot * 0.45);
  for (double t : nice_ticks(d0, d1, 8)) {
    f.line({px - 3, Y(t)}, {px + pw, Y(t)}, kGrid, 0.3);
    f.text({px - 5, Y(t) + 2.5}, fmt_fixed(t, 0), 6.5, Anchor::End, false, kInk);
  }
  f.text({px - 48, py + ph / 2}, "TVD below the model datum (m)", 7.5, Anchor::Middle, false, kInk, -90);
  std::vector<double> xs;
  for (std::size_t k = 0; k < nw; ++k) {
    const double x = px + slot * (static_cast<double>(k) + 0.5) - tw / 2;
    xs.push_back(x);
    const auto& w = p.wells[k];
    f.rect(x, py, tw, ph, rgb(250, 250, 250), kInk, 0.5);
    f.text({x + tw / 2, py - 30}, w.well, 8.5, Anchor::Middle, true, kInk);
    for (auto& u : w.units) {
      const double a = std::max(d0, u.top), b = std::min(d1, u.base);
      if (b > a) f.rect(x, Y(a), 6, Y(b) - Y(a), unit_color(c, static_cast<std::size_t>(std::max(0, u.unit))), none());
    }
    if (!w.tracks.empty()) log_track(f, x + 6, py, tw - 6, ph, w.tracks[0], w.index, d0, d1, track_color(0));
    for (auto& a : w.ages) {
      f.path({{x - 6, Y(a.depth) - 3}, {x, Y(a.depth)}, {x - 6, Y(a.depth) + 3}}, true, rgb(60, 60, 60), none());
      f.text({x - 7, Y(a.depth) + 2}, fmt_shortest(a.age), 5, Anchor::End, false, rgb(60, 60, 60));
    }
  }
  for (std::size_t h = 0; h < p.horizons.size(); ++h) {
    Color col = unit_color(c, h + 1);
    Color dark = rgb(static_cast<std::uint8_t>(col.r * 0.55), static_cast<std::uint8_t>(col.g * 0.55), static_cast<std::uint8_t>(col.b * 0.55));
    std::vector<Pt> seg;
    for (std::size_t k = 0; k < nw && k < p.depth[h].size(); ++k) {
      const double d = p.depth[h][k];
      if (!std::isfinite(d)) {
        if (seg.size() >= 2) f.path(seg, false, none(), dark, 1.1);
        seg.clear();
        continue;
      }
      if (!seg.empty()) seg.push_back({xs[k], Y(d)});
      seg.push_back({xs[k], Y(d)});
      seg.push_back({xs[k] + tw, Y(d)});
    }
    if (seg.size() >= 2) f.path(seg, false, none(), dark, 1.1);
    for (std::size_t k = nw; k-- > 0;)
      if (k < p.depth[h].size() && std::isfinite(p.depth[h][k])) {
        f.text({xs[k] + tw + 3, Y(p.depth[h][k]) + 2.5}, p.horizons[h], 6.5, Anchor::Start, true, dark);
        break;
      }
  }
  std::vector<LegendEntry> es;
  for (std::size_t h = 0; h < p.horizons.size(); ++h) {
    Color col = unit_color(c, h + 1);
    Color dark = rgb(static_cast<std::uint8_t>(col.r * 0.55), static_cast<std::uint8_t>(col.g * 0.55), static_cast<std::uint8_t>(col.b * 0.55));
    es.push_back({"Correlation lines", p.horizons[h] + ": " + (h < p.basis.size() ? p.basis[h] : std::string()), none(), true, dark, 1.1});
  }
  if (!p.wells.empty() && !p.wells[0].tracks.empty()) es.push_back({"Track", p.wells[0].tracks[0].name + " (" + p.wells[0].tracks[0].unit + ")", none(), true, track_color(0), 0.7});
  double ly = legend(f, 842 - 200, py, 190, es, 595.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 200, ly, 190);
  footer(f, c);
  return f;
}

}  // namespace stratonix::fig

namespace stratonix::fig {

Figure seismic_section_figure(const FigureContext& c, const repr::Section& s, const std::vector<float>& amp, const std::string& how) {
  Figure f(842, 540);
  const auto& inf = c.model->info();
  const double z0 = s.z_top, z1 = s.z_base;
  const double x0 = 70, x1 = 842 - 232, y0 = 58;
  const double h_ve = (x1 - x0) * (z0 - z1) / s.length * c.vertical_exaggeration;
  const double y1 = y0 + std::clamp(h_ve, 120.0, 540.0 - 58 - 92);
  const Frame2 fr{x0, y0, x1, y1, 0, s.length, z0, z1};
  title(f, (c.title.empty() ? c.model_name : c.title) + " - seismic along " + s.name,
        "Amplitude sampled from the SEG-Y volume " + how + "; the model's horizons and faults drawn over it");
  // amplitude image, symmetric colour scale at the 98th percentile of |amplitude|
  std::vector<float> mags;
  for (float a : amp) if (std::isfinite(a)) mags.push_back(std::fabs(a));
  float clip = 1;
  if (!mags.empty()) {
    std::nth_element(mags.begin(), mags.begin() + static_cast<long>(0.98 * static_cast<double>(mags.size() - 1)), mags.end());
    clip = std::max(1e-12f, mags[static_cast<std::size_t>(0.98 * static_cast<double>(mags.size() - 1))]);
  }
  Image im;
  im.w = static_cast<int>(s.ns);
  im.h = static_cast<int>(s.nz);
  im.rgba.assign(static_cast<std::size_t>(im.w) * static_cast<std::size_t>(im.h) * 4, 0);
  for (std::size_t r = 0; r < s.nz; ++r)
    for (std::size_t q = 0; q < s.ns; ++q) {
      const float a = r * s.ns + q < amp.size() ? amp[r * s.ns + q] : NAN;
      std::uint8_t* px = &im.rgba[4 * (r * s.ns + q)];
      if (!std::isfinite(a)) { px[0] = px[1] = px[2] = 235; px[3] = 255; continue; }
      const double t = std::clamp(static_cast<double>(a / clip), -1.0, 1.0);
      // blue (negative) - white - red (positive)
      px[0] = static_cast<std::uint8_t>(t < 0 ? 255 * (1 + t) : 255);
      px[1] = static_cast<std::uint8_t>(255 * (1 - std::fabs(t)));
      px[2] = static_cast<std::uint8_t>(t > 0 ? 255 * (1 - t) : 255);
      px[3] = 255;
    }
  f.image(x0, y0, x1 - x0, y1 - y0, f.own(std::move(im)));
  f.clip(fr.x0, fr.y0, fr.x1 - fr.x0, fr.y1 - fr.y0);
  std::vector<LegendEntry> es;
  for (std::size_t h = 0; h < s.horizon_z.size() && h < inf.horizons.size(); ++h) {
    Color col = unit_color(c, h + 1);
    Color dark = rgb(static_cast<std::uint8_t>(col.r * 0.45), static_cast<std::uint8_t>(col.g * 0.45), static_cast<std::uint8_t>(col.b * 0.45));
    std::vector<Pt> seg;
    for (std::size_t q = 0; q < s.ns; ++q) {
      const double z = s.horizon_z[h][q];
      if (!std::isfinite(z)) { if (seg.size() >= 2) f.path(seg, false, none(), dark, 1.2); seg.clear(); continue; }
      seg.push_back(fr.p(s.station[q], z));
    }
    if (seg.size() >= 2) f.path(seg, false, none(), dark, 1.2);
    es.push_back({"Model horizons", inf.horizons[h].name, none(), true, dark, 1.2});
  }
  for (std::size_t k = 0; k < s.fault_z.size(); ++k) {
    std::vector<Pt> seg;
    for (std::size_t q = 0; q < s.ns; ++q) {
      const double z = s.fault_z[k][q];
      if (!std::isfinite(z)) { if (seg.size() >= 2) f.path(seg, false, none(), kInk, 1.4, {5, 2}); seg.clear(); continue; }
      seg.push_back(fr.p(s.station[q], z));
    }
    if (seg.size() >= 2) f.path(seg, false, none(), kInk, 1.4, {5, 2});
  }
  if (!s.fault_z.empty()) es.push_back({"Model faults", "Fault (dashed)", none(), true, kInk, 1.4});
  f.unclip();
  f.rect(x0, y0, x1 - x0, y1 - y0, none(), kInk, 0.6);
  for (double t : nice_ticks(0, s.length, 7)) f.text({fr.p(t, z1).x, y1 + 11}, km_or_m(t), 6.5, Anchor::Middle, false, kInk);
  for (double t : nice_ticks(z1, z0, 6)) f.text({x0 - 4, fr.p(0, t).y + 2.5}, fmt_fixed(-t, 0), 6.5, Anchor::End, false, kInk);
  f.text({(x0 + x1) / 2, y1 + 24}, "Distance along the section (m)", 7.5, Anchor::Middle, false, kInk);
  f.text({x0 - 42, (y0 + y1) / 2}, "Depth below the model datum (m)", 7.5, Anchor::Middle, false, kInk, -90);
  es.push_back({"Amplitude", "Red positive, blue negative (98 % clip)", none(), false, rgb(230, 80, 80)});
  double ly = legend(f, 842 - 200, y0, 190, es, 540.0 - 30 - key_box_height(c) - 8) + 8;
  key_box(f, c, 842 - 200, ly, 190);
  footer(f, c);
  return f;
}

}  // namespace stratonix::fig
