// Figure set: sections, maps, 3D views, fence diagram; writing in four formats.
#include <algorithm>
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
  if (u + 1 == inf.units.size()) return Color::hex(kPalette[11]);
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

Figure map_figure(const FigureContext& c, int h, const repr::FaultedSurface& surf,
                  const std::vector<std::vector<repr::CutoffLine>>& cutoffs, const repr::Closure* cl) {
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
    auto fs = repr::horizon_grid_faulted(*c.model, static_cast<int>(k), g);
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
  for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
    auto fs = repr::horizon_grid_faulted(*c.model, static_cast<int>(h), g);
    std::vector<std::vector<repr::CutoffLine>> cut;
    for (std::size_t ff = 0; ff < inf.faults.size(); ++ff) {
      auto fsurf = repr::fault_surface(*c.model, static_cast<int>(ff), std::max(10.0, g.dx / 2));
      cut.push_back(repr::cutoff_lines(*c.model, fsurf, static_cast<int>(h)));
    }
    repr::Closure cl;
    bool has = static_cast<int>(h) == c.closure_horizon;
    if (has) cl = repr::closure_analysis(fs.z, c.faults_seal);
    auto fgr = map_figure(c, static_cast<int>(h), fs, cut, has ? &cl : nullptr);
    out.push_back(write_figure(fgr, folder, figure_basename(c, "map", inf.horizons[h].name), c.model_name + " - " + inf.horizons[h].name + " depth map",
                               "Depth map of " + inf.horizons[h].name + " with fault polygons and wells" + (has ? ", closure outline and spill point." : "."),
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

}  // namespace stratonix::fig
