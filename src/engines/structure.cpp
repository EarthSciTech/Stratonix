// Unfaulted structure: fold primitives, stratigraphy, imported surfaces, isochores, corrections.
#include "stratonix/engines/structure.hpp"

#include <algorithm>
#include <cmath>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/wells.hpp"

namespace stratonix {

// ---- primitives ----------------------------------------------------------------------------
double FoldPrimitive::bump(double t, double w) {
  if (w <= 0) return 0;
  double a = std::fabs(t);
  return a < w ? 0.5 * (1 + std::cos(kPi * a / w)) : 0.0;
}
void FoldPrimitive::axis_coords(double x, double y, double& u, double& v) const {
  double ax = std::sin(azimuth), ay = std::cos(azimuth);  // along axis
  double rx = ay, ry = -ax;                               // to the right looking along the axis
  double dx = x - cx, dy = y - cy;
  u = dx * ax + dy * ay;
  v = dx * rx + dy * ry;
}
double FoldPrimitive::shape(double x, double y) const {
  double u, v;
  axis_coords(x, y, u, v);
  if (type == Type::Anticline || type == Type::Syncline) {
    double gu = half_length > 0 ? bump(u, half_length) : 1.0;
    double gv = bump(v, v >= 0 ? w_right : w_left);
    return gu * gv;
  }
  double rho = std::sqrt((u / ra) * (u / ra) + (v / rb) * (v / rb));
  return bump(rho, 1.0);
}
double FoldPrimitive::rim(double x, double y) const {
  if (rim_depth == 0 || rim_width <= 0) return 0;
  double u, v;
  axis_coords(x, y, u, v);
  if (type == Type::Anticline || type == Type::Syncline) {
    double w = v >= 0 ? w_right : w_left;
    double gu = half_length > 0 ? bump(u, half_length + 2 * rim_width) : 1.0;
    return gu * bump(std::fabs(v) - w - rim_width, rim_width);
  }
  double rho = std::sqrt((u / ra) * (u / ra) + (v / rb) * (v / rb));
  return bump((rho - 1.0) * ra - rim_width, rim_width);
}
double FoldPrimitive::tilt(double x, double y) const {
  if (plunge == 0) return 0;
  double u, v;
  axis_coords(x, y, u, v);
  return -std::tan(plunge) * u;
}
double FoldPrimitive::Offset::eval(const FoldPrimitive& p, double x, double y) const {
  if (cyl) {
    double u, v;
    p.axis_coords(x, y, u, v);
    double f = (v - v0) / dv;
    if (f < 0 || f >= static_cast<double>(val.size() - 1)) return outside;
    std::size_t i = static_cast<std::size_t>(f);
    double w = f - static_cast<double>(i);
    return val[i] * (1 - w) + val[i + 1] * w;
  }
  double fx = (x - x0) / dx, fy = (y - y0) / dx;
  if (fx < 0 || fy < 0 || fx >= static_cast<double>(nx - 1) || fy >= static_cast<double>(ny - 1))
    return p.tilt(x, y) - t * std::sqrt(1 + std::tan(p.plunge) * std::tan(p.plunge)) + t;
  std::size_t i = static_cast<std::size_t>(fx), j = static_cast<std::size_t>(fy);
  double a = fx - static_cast<double>(i), b = fy - static_cast<double>(j);
  auto at = [&](std::size_t ii, std::size_t jj) { return val[jj * nx + ii]; };
  return (at(i, j) * (1 - a) + at(i + 1, j) * a) * (1 - b) + (at(i, j + 1) * (1 - a) + at(i + 1, j + 1) * a) * b;
}

double FoldPrimitive::contribution(std::size_t h, double x, double y) const {
  if (h < offset.size() && offset[h]) return offset[h]->eval(*this, x, y);
  double r = relief[h];
  double c = r * shape(x, y) + tilt(x, y);
  if (rim_depth != 0 && relief_max > 0) c -= rim_depth * (std::fabs(r) / relief_max) * rim(x, y) * (r >= 0 ? 1 : -1);
  return c;
}

// ---- small dense solver (Gaussian elimination with partial pivoting) -----------------------
static std::vector<double> solve_dense(std::vector<double> a, std::vector<double> b, std::size_t n) {
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; ++i)
      if (std::fabs(a[i * n + k]) > std::fabs(a[p * n + k])) p = i;
    if (std::fabs(a[p * n + k]) < 1e-14) throw Error(ErrorCode::Build, "well-correction system is singular (wells too close?)");
    if (p != k) {
      for (std::size_t j = 0; j < n; ++j) std::swap(a[k * n + j], a[p * n + j]);
      std::swap(b[k], b[p]);
    }
    for (std::size_t i = k + 1; i < n; ++i) {
      double f = a[i * n + k] / a[k * n + k];
      for (std::size_t j = k; j < n; ++j) a[i * n + j] -= f * a[k * n + j];
      b[i] -= f * b[k];
    }
  }
  std::vector<double> x(n);
  for (std::size_t i = n; i-- > 0;) {
    double s = b[i];
    for (std::size_t j = i + 1; j < n; ++j) s -= a[i * n + j] * x[j];
    x[i] = s / a[i * n + i];
  }
  return x;
}

// ---- structure -----------------------------------------------------------------------------
Structure::Structure(const desc::Description& d, const BuildInputs& in, BuildLog& log) {
  const std::size_t nh = d.horizons.size();
  for (auto& [k, g] : in.grids) owned_[k] = g;
  kind_.resize(nh);
  zref_.assign(nh, 0.0);
  ref_depth_.assign(nh, kNaN);
  surf_.assign(nh, nullptr);
  surf_sign_.assign(nh, -1.0);
  iso_from_.assign(nh, -1);
  iso_const_.assign(nh, 0.0);
  iso_grid_.assign(nh, nullptr);
  corr_.resize(nh);
  max_corr_.assign(nh, 0.0);
  hseries_.assign(nh, 0);
  tilt_gx_.assign(nh, 0.0);
  tilt_gy_.assign(nh, 0.0);
  trend_.assign(nh, 0);
  trend_a_.assign(nh, 0.0);
  trend_b_.assign(nh, 0.0);
  trend_cx_.assign(nh, 0.0);
  trend_cy_.assign(nh, 0.0);
  term_.assign(nh + 1, std::string());

  // Reference location.
  const desc::WellDesc* rw = d.ref_well.empty() ? nullptr : d.well(d.ref_well);
  double rx = rw ? rw->head.si(0) : d.ref_point.si(0);
  double ry = rw ? rw->head.si(1) : d.ref_point.si(1);
  reg_x0_ = rx;
  reg_y0_ = ry;
  if (d.regional_type == "planar") {
    double dip = d.regional_dip.si(), az = d.regional_dip_azimuth.si();
    reg_gx_ = -std::tan(dip) * std::sin(az);
    reg_gy_ = -std::tan(dip) * std::cos(az);
  }

  // Horizon kinds and grids.
  auto grid_of = [&](const std::string& id, const std::string& what) -> const Grid2D* {
    auto it = owned_.find(id);
    if (it == owned_.end()) throw Error(ErrorCode::Build, "input grid '" + id + "' for " + what + " was not loaded");
    return &it->second;
  };
  for (std::size_t h = 0; h < nh; ++h) {
    auto& hd = d.horizons[h];
    if (!hd.surface_input.empty()) {
      kind_[h] = Kind::Surface;
      surf_[h] = grid_of(hd.surface_input, "horizon " + hd.name);
      surf_sign_[h] = hd.surface_z_is == "elevation" ? 1.0 : -1.0;
    } else if (!hd.iso_from.empty()) {
      kind_[h] = Kind::Isochore;
      iso_from_[h] = d.horizon_index(hd.iso_from);
      if (hd.iso_thickness.present) iso_const_[h] = hd.iso_thickness.si();
      else iso_grid_[h] = grid_of(hd.iso_input, "isochore of " + hd.name);
    } else {
      kind_[h] = Kind::Parametric;
    }
    if (!hd.iso_trend.empty()) has_trend_ = true;
    if (hd.iso_trend == "linear") {
      trend_[h] = 1;
      trend_a_[h] = hd.iso_gradient.si(0);
      trend_b_[h] = hd.iso_gradient.si(1);
    } else if (hd.iso_trend == "radial") {
      trend_[h] = 2;
      trend_a_[h] = hd.iso_gradient.si();
      trend_cx_[h] = hd.iso_centre.si(0);
      trend_cy_[h] = hd.iso_centre.si(1);
    }
  }
  // Stratigraphic series (STX-TN-07 section 1), youngest first.
  if (d.series.empty()) {
    SeriesInfo si;
    si.first = 0;
    si.last = nh - 1;
    series_.push_back(si);
  }
  for (std::size_t k = 0; k < d.series.size(); ++k) {
    const auto& sd = d.series[k];
    SeriesInfo si;
    si.first = static_cast<std::size_t>(d.horizon_index(sd.horizons.front()));
    si.last = static_cast<std::size_t>(d.horizon_index(sd.horizons.back()));
    si.base = sd.base == "erosional" ? 1 : sd.base == "onlap" ? 2 : 0;
    if (!sd.clino_profile.empty()) {
      si.clino = true;
      si.sigmoid = sd.clino_profile == "sigmoidal";
      si.height = sd.clino_height.si();
      si.width = sd.clino_width.si();
      si.ax = std::sin(sd.clino_direction.si());
      si.ay = std::cos(sd.clino_direction.si());
      si.bx = sd.clino_first_break.si(0);
      si.by = sd.clino_first_break.si(1);
      si.rate = sd.clino_progradation.si();
      for (std::size_t h = si.first; h <= si.last; ++h) kind_[h] = Kind::Clinoform;
      log.value(sd.name + ".clinoform_foreset_dip_max", fmt_fixed(rad2deg(std::atan(si.height / si.width)), 3));
    }
    if (si.base != 0 || si.clino) has_rel_ = true;
    for (std::size_t h = si.first; h <= si.last; ++h) hseries_[h] = static_cast<int>(k);
    if (sd.tilt_dip.present && sd.tilt_dip.si() > 0) {
      // tilted after deposition: this series and every older one, about the reference location
      double gx = -std::tan(sd.tilt_dip.si()) * std::sin(sd.tilt_azimuth.si()), gy = -std::tan(sd.tilt_dip.si()) * std::cos(sd.tilt_azimuth.si());
      for (std::size_t h = si.first; h < nh; ++h) {
        tilt_gx_[h] += gx;
        tilt_gy_[h] += gy;
      }
      has_tilt_ = true;
    }
    series_.push_back(si);
  }

  // Reference depths (repeated until settled: clinoforms hang from the horizon below them).
  for (int pass = 0; pass < static_cast<int>(nh) + 1; ++pass)
  for (std::size_t h = 0; h < nh; ++h) {
    if (!std::isnan(ref_depth_[h])) continue;
    auto& hd = d.horizons[h];
    if (kind_[h] == Kind::Clinoform) {
      std::size_t b = series_[static_cast<std::size_t>(hseries_[h])].last + 1;
      if (b < nh && !std::isnan(ref_depth_[b])) ref_depth_[h] = ref_depth_[b] - clino_offset(h, rx, ry);
      continue;
    }
    if (rw)
      for (auto& [n, p] : rw->tops)
        if (n == hd.name) ref_depth_[h] = p.si();
    if (std::isnan(ref_depth_[h]) && hd.depth.present) ref_depth_[h] = hd.depth.si();
    if (std::isnan(ref_depth_[h]) && hd.thickness_above.present && h > 0) ref_depth_[h] = ref_depth_[h - 1] + hd.thickness_above.si();
    if (std::isnan(ref_depth_[h]) && kind_[h] == Kind::Surface) ref_depth_[h] = -surf_sign_[h] * surf_[h]->sample_clamped(rx, ry);
    if (std::isnan(ref_depth_[h]) && kind_[h] == Kind::Isochore && !std::isnan(ref_depth_[static_cast<std::size_t>(iso_from_[h])]))
      ref_depth_[h] = ref_depth_[static_cast<std::size_t>(iso_from_[h])] + iso_thickness(h, rx, ry, true);
  }
  for (std::size_t h = 0; h < nh; ++h)
    if (std::isnan(ref_depth_[h]))
      throw Error(ErrorCode::Build, "cannot determine the reference depth of horizon " + d.horizons[h].name, "stratigraphy.horizons");

  // Primitives with relief interpolated from anchors in reference depth.
  for (auto& pd : d.primitives) {
    FoldPrimitive p;
    p.name = pd.name;
    p.type = pd.type == "anticline" ? FoldPrimitive::Type::Anticline
             : pd.type == "syncline" ? FoldPrimitive::Type::Syncline
             : pd.type == "dome"     ? FoldPrimitive::Type::Dome
                                     : FoldPrimitive::Type::Basin;
    const bool det = pd.type == "detachment_fold";
    if (det) p.type = FoldPrimitive::Type::Anticline;
    double sign = (p.type == FoldPrimitive::Type::Syncline || p.type == FoldPrimitive::Type::Basin) ? -1.0 : 1.0;
    p.cx = pd.crest.si(0);
    p.cy = pd.crest.si(1);
    p.azimuth = pd.axis_azimuth.present ? pd.axis_azimuth.si() : 0.0;
    p.plunge = pd.plunge.present ? pd.plunge.si() : 0.0;
    p.half_length = pd.half_length.present ? pd.half_length.si() : 0.0;
    p.w_left = pd.half_width_left.present ? pd.half_width_left.si() : 0.0;
    p.w_right = pd.half_width_right.present ? pd.half_width_right.si() : 0.0;
    p.ra = pd.radius_major.present ? pd.radius_major.si() : 0.0;
    p.rb = pd.radius_minor.present ? pd.radius_minor.si() : p.ra;
    p.rim_depth = pd.rim_depth.present ? pd.rim_depth.si() : 0.0;
    p.rim_width = pd.rim_width.present ? pd.rim_width.si() : 0.0;
    p.relief.assign(nh, 0.0);
    if (det) {
      // STX-TN-06 section 3: area balance; relief S * D_h / w falls to zero at the detachment.
      p.w_left = p.w_right = pd.half_width.si();
      const double S = pd.shortening.si(), w = p.w_left, zd = ref_depth_[static_cast<std::size_t>(d.horizon_index(pd.detachment))];
      for (std::size_t h = 0; h < nh; ++h) {
        double r = ref_depth_[h] < zd ? S * (zd - ref_depth_[h]) / w : 0.0;
        p.relief[h] = r;
        p.relief_max = std::max(p.relief_max, r);
        if (r > 0) log.value("primitive." + p.name + "." + d.horizons[h].name + ".relief", fmt_fixed(r, 3));
      }
      add_copies(pd, p, d, log);
      continue;
    }
    if (!pd.mechanics.empty()) {
      // STX-TN-06 section 4: the anchor is the shallowest horizon given relief.
      std::size_t ha = nh;
      double ra = 0;
      for (auto& [hn, q] : pd.relief) {
        std::size_t h = static_cast<std::size_t>(d.horizon_index(hn));
        if (ha == nh || ref_depth_[h] < ref_depth_[ha]) { ha = h; ra = q.si(); }
      }
      for (std::size_t h = 0; h < nh; ++h) p.relief[h] = sign * ra;
      p.relief_max = std::fabs(ra);
      if (pd.mechanics == "parallel") {
        p.offset.assign(nh, nullptr);
        for (std::size_t h = ha + 1; h < nh; ++h)
          if (kind_[h] == Kind::Parametric) p.offset[h] = parallel_offset(p, ha, ref_depth_[h] - ref_depth_[ha], d);
      }
      add_copies(pd, p, d, log);
      continue;
    }
    std::vector<std::pair<double, double>> anchors;  // (depth, relief)
    for (auto& [hn, q] : pd.relief) {
      int h = d.horizon_index(hn);
      anchors.emplace_back(ref_depth_[static_cast<std::size_t>(h)], q.si());
    }
    if (d.base_type == "flat") {
      bool has = false;
      for (auto& [hn, q] : pd.relief) has = has || d.horizon_index(hn) == static_cast<int>(nh) - 1;
      if (!has) anchors.emplace_back(ref_depth_[nh - 1], 0.0);
    }
    std::sort(anchors.begin(), anchors.end());
    for (std::size_t h = 0; h < nh; ++h) {
      if (kind_[h] != Kind::Parametric) continue;
      double z = ref_depth_[h];
      double r;
      if (z <= anchors.front().first) r = anchors.front().second;
      else if (z >= anchors.back().first) r = anchors.back().second;
      else {
        std::size_t k = 1;
        while (anchors[k].first < z) ++k;
        auto [z0, r0] = anchors[k - 1];
        auto [z1, r1] = anchors[k];
        r = z1 > z0 ? r0 + (r1 - r0) * (z - z0) / (z1 - z0) : r1;
      }
      if (d.base_type == "flat" && h == nh - 1) r = 0;
      p.relief[h] = sign * r;
      p.relief_max = std::max(p.relief_max, std::fabs(r));
    }
    add_copies(pd, p, d, log);
  }

  // Constant terms: honour the reference location by construction.
  for (std::size_t h = 0; h < nh; ++h) {
    if (kind_[h] != Kind::Parametric) continue;
    double s = regional(rx, ry);
    for (auto& p : prims_) s += p.contribution(h, rx, ry);
    zref_[h] = -ref_depth_[h] - s;
  }

  // Correction field for control wells (the reference well is included with its residual).
  corr_range_ = d.correction_range.present ? d.correction_range.si() : 0.0;
  for (std::size_t h = 0; h < nh; ++h) {
    std::vector<double> xs, ys, rs;
    for (auto& w : d.wells) {
      if (w.role != "control") continue;
      for (auto& [n, p] : w.tops) {
        if (n != d.horizons[h].name) continue;
        WellPath wp = make_well_path(w, d);
        Vec3 at = wp.at_elevation(-p.si());
        xs.push_back(at.x);
        ys.push_back(at.y);
        rs.push_back(-p.si() - z_uncorrected(h, at.x, at.y));
      }
    }
    double maxr = 0;
    for (double r : rs) maxr = std::max(maxr, std::fabs(r));
    if (xs.empty() || maxr < 1e-9) continue;
    if (corr_range_ <= 0)
      throw Error(ErrorCode::Build, "horizon " + d.horizons[h].name + " misses a control well by " + fmt_fixed(maxr, 3) +
                                        " m; give structure.correction_range so it can be corrected", "structure.correction_range");
    std::size_t n = xs.size();
    std::vector<double> a(n * n);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j) {
        double dd = std::hypot(xs[i] - xs[j], ys[i] - ys[j]) / corr_range_;
        a[i * n + j] = std::exp(-dd * dd);
      }
    corr_[h] = {xs, ys, solve_dense(a, rs, n)};
    max_corr_[h] = maxr;
    log.note("horizon " + d.horizons[h].name + ": control wells honoured by a correction field (max " +
             fmt_fixed(maxr, 2) + " m at the wells, range " + fmt_fixed(corr_range_, 0) + " m)");
    log.value("correction_max_m." + d.horizons[h].name, fmt_fixed(maxr, 4));
  }
  if (has_rel_ || std::any_of(trend_.begin(), trend_.end(), [](int t) { return t != 0; })) report_terminations(d, log);
}

void Structure::add_copies(const desc::PrimitiveDesc& pd, FoldPrimitive p, const desc::Description& d, BuildLog& log) {
  // STX-TN-07: a fold of an older series leaves younger series unfolded; growth scales the relief.
  if (!pd.series.empty()) {
    std::size_t first = 0;
    for (std::size_t k = 0; k < d.series.size(); ++k)
      if (d.series[k].name == pd.series) first = static_cast<std::size_t>(d.horizon_index(d.series[k].horizons.front()));
    for (std::size_t h = 0; h < first && h < p.relief.size(); ++h) {
      p.relief[h] = 0;
      if (h < p.offset.size()) p.offset[h] = nullptr;
    }
  }
  if (pd.growth.present()) {
    const int on = d.horizon_index(pd.growth.onset), en = d.horizon_index(pd.growth.end);
    for (std::size_t h = 0; h < p.relief.size(); ++h) {
      p.relief[h] *= growth_fraction(h, on, en);
      log.value("primitive." + p.name + "." + d.horizons[h].name + ".growth_fraction", fmt_fixed(growth_fraction(h, on, en), 4));
    }
  }
  // STX-TN-06 section 5: en echelon copies, each moved along its axis and to its right.
  const int n = pd.array_count;
  if (n <= 1) { prims_.push_back(std::move(p)); return; }
  const double sa = pd.array_step_along.si(), sc = pd.array_step_across.si();
  const double ax = std::sin(p.azimuth), ay = std::cos(p.azimuth), rx = ay, ry = -ax;
  const double cx = p.cx, cy = p.cy;
  for (int k = 0; k < n; ++k) {
    FoldPrimitive q = p;
    q.cx = cx + k * (sa * ax + sc * rx);
    q.cy = cy + k * (sa * ay + sc * ry);
    q.name = p.name + "#" + std::to_string(k + 1);
    for (std::size_t h = 0; h < q.offset.size(); ++h)
      if (q.offset[h] && !q.offset[h]->cyl) q.offset[h] = parallel_offset(q, offset_anchor_, q.offset[h]->t, d);
    log.value("primitive." + q.name + ".crest", fmt_fixed(q.cx, 2) + ", " + fmt_fixed(q.cy, 2));
    prims_.push_back(std::move(q));
  }
}

std::shared_ptr<const FoldPrimitive::Offset> Structure::parallel_offset(const FoldPrimitive& p, std::size_t ha, double t,
                                                                        const desc::Description& d) {
  // Lower offset of the anchor surface c at normal distance t (STX-TN-06 section 4):
  // z(x) = min over |x' - x| <= t of c(x') - sqrt(t^2 - |x' - x|^2); cusps are kept.
  offset_anchor_ = ha;
  auto o = std::make_shared<FoldPrimitive::Offset>();
  o->t = t;
  const double ra = p.relief[ha];
  const bool cylindrical = (p.type == FoldPrimitive::Type::Anticline || p.type == FoldPrimitive::Type::Syncline) &&
                           p.half_length <= 0 && p.plunge == 0;
  if (t <= 0) {
    o->cyl = true;
    o->v0 = 0;
    o->dv = 1;
    o->val = {0, 0};
    return o;
  }
  if (cylindrical) {
    auto c = [&](double v) { return ra * FoldPrimitive::bump(v, v >= 0 ? p.w_right : p.w_left); };
    o->cyl = true;
    o->dv = 2.0;
    o->v0 = -p.w_left - t - 20;
    const double v1 = p.w_right + t + 20, step = 0.5;
    const std::size_t n = static_cast<std::size_t>(std::ceil((v1 - o->v0) / o->dv)) + 1;
    o->val.resize(n);
    const int m = static_cast<int>(std::ceil(t / step));
    for (std::size_t i = 0; i < n; ++i) {
      double v = o->v0 + o->dv * static_cast<double>(i), best = 1e300;
      for (int k = -m; k <= m; ++k) {
        double dd = std::min(std::fabs(k * step), t);
        best = std::min(best, c(v + k * step) - std::sqrt(t * t - dd * dd));
      }
      o->val[i] = best + t;
    }
    o->outside = 0;
    return o;
  }
  const Domain dom = d.domain();
  const double dx = 0.5 * dom.spacing, margin = t + 4 * dx;
  o->cyl = false;
  o->dx = dx;
  o->x0 = dom.x0 - margin;
  o->y0 = dom.y0 - margin;
  o->nx = static_cast<std::size_t>(std::ceil((dom.sx + 2 * margin) / dx)) + 1;
  o->ny = static_cast<std::size_t>(std::ceil((dom.sy + 2 * margin) / dx)) + 1;
  o->val.resize(o->nx * o->ny);
  const double kt = std::sqrt(1 + std::tan(p.plunge) * std::tan(p.plunge));
  auto c = [&](double x, double y) { return ra * p.shape(x, y) + p.tilt(x, y); };
  const int m = static_cast<int>(std::ceil(t / dx));
  // support of the shape, expanded by t: nodes beyond it see a plane
  auto near = [&](double x, double y) {
    double u, v;
    p.axis_coords(x, y, u, v);
    if (p.type == FoldPrimitive::Type::Anticline || p.type == FoldPrimitive::Type::Syncline)
      return v > -p.w_left - t - dx && v < p.w_right + t + dx && (p.half_length <= 0 || std::fabs(u) < p.half_length + t + dx);
    return std::hypot(u, v) < std::max(p.ra, p.rb) + t + dx;
  };
  for (std::size_t j = 0; j < o->ny; ++j)
    for (std::size_t i = 0; i < o->nx; ++i) {
      double x = o->x0 + dx * static_cast<double>(i), y = o->y0 + dx * static_cast<double>(j);
      if (!near(x, y)) { o->val[j * o->nx + i] = p.tilt(x, y) - t * kt + t; continue; }
      double best = 1e300;
      for (int b = -m; b <= m; ++b)
        for (int a = -m; a <= m; ++a) {
          double dd2 = (a * a + b * b) * dx * dx;
          if (dd2 > t * t) continue;
          best = std::min(best, c(x + a * dx, y + b * dx) - std::sqrt(t * t - dd2));
        }
      o->val[j * o->nx + i] = best + t;
    }
  return o;
}

double Structure::missing_section(int u, double x, double y) const {
  const std::size_t nh = kind_.size();
  if (u < 1 || static_cast<std::size_t>(u) >= nh || !has_rel_) return 0.0;
  const std::size_t a = static_cast<std::size_t>(u) - 1, b = static_cast<std::size_t>(u);
  if (hseries_[a] != hseries_[b]) return 0.0;
  std::vector<double> raw(nh), ero(nh);
  z_all_raw(x, y, raw.data());
  ero = raw;
  apply_relations(x, y, ero.data(), nullptr, true);
  return std::max(0.0, (raw[a] - raw[b]) - (ero[a] - ero[b]));
}

std::string Structure::termination_at(int u, double x0, double y0, double radius) const {
  const std::size_t nh = kind_.size();
  if (u < 1 || static_cast<std::size_t>(u) >= nh) return "";
  static const char* names[] = {"", "clinoform", "pinch-out", "onlap", "truncation"};
  std::vector<double> raw(nh), eff(nh);
  std::vector<unsigned char> fl(nh);
  const std::size_t a = static_cast<std::size_t>(u) - 1, b = static_cast<std::size_t>(u);
  int best = 0;
  const double off[][2] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}, {0.7, 0.7}, {-0.7, 0.7}, {0.7, -0.7}, {-0.7, -0.7}};
  for (auto& o : off) {
    const double x = x0 + o[0] * radius, y = y0 + o[1] * radius;
    z_detail(x, y, raw.data(), eff.data(), fl.data());
    if (eff[a] - eff[b] >= 0.01) continue;
    int k = 0;
    if ((fl[a] | fl[b]) & 1) k = 4;
    else if ((fl[a] | fl[b]) & 2) k = 3;
    else if (trend_clamped(b, x, y)) k = 2;
    else if (kind_[a] == Kind::Clinoform || kind_[b] == Kind::Clinoform) k = 1;
    best = std::max(best, k);
  }
  return names[best];
}

void Structure::report_terminations(const desc::Description& d, BuildLog& log) {
  // Sample the domain: which units terminate, how, and the largest missing section.
  const std::size_t nh = kind_.size();
  const Domain dom = d.domain();
  const int n = 121;
  std::vector<double> raw(nh), eff(nh);
  std::vector<unsigned char> fl(nh);
  std::vector<unsigned> kind(nh + 1, 0);  // bits: 8 truncation, 4 onlap, 2 pinch-out, 1 clinoform
  std::vector<double> miss(nh + 1, 0.0), mx(nh + 1, 0.0), my(nh + 1, 0.0);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
      const double x = dom.x0 + dom.sx * i / (n - 1), y = dom.y0 + dom.sy * j / (n - 1);
      z_detail(x, y, raw.data(), eff.data(), fl.data());
      for (std::size_t u = 1; u < nh; ++u) {
        const std::size_t a = u - 1, b = u;
        if (eff[a] - eff[b] < 0.01) {
          unsigned k = 0;
          if ((fl[a] | fl[b]) & 1) k = 8;
          else if ((fl[a] | fl[b]) & 2) k = 4;
          else if (trend_clamped(b, x, y)) k = 2;
          else if (kind_[a] == Kind::Clinoform || kind_[b] == Kind::Clinoform) k = 1;
          kind[u] |= k;
        }
        if (has_rel_) {
          double m = missing_section(static_cast<int>(u), x, y);
          if (m > miss[u]) { miss[u] = m; mx[u] = x; my[u] = y; }
        }
      }
    }
  for (std::size_t u = 1; u < nh; ++u) {
    term_[u].clear();
    for (auto [bit, nm] : {std::pair<unsigned, const char*>{8, "truncation"}, {4, "onlap"}, {2, "pinch-out"}, {1, "clinoform"}})
      if (kind[u] & bit) term_[u] += (term_[u].empty() ? "" : ", ") + std::string(nm);
    if (kind[u]) log.value("unit." + d.horizons[u - 1].name + "_" + d.horizons[u].name + ".termination", term_[u]);
    if (miss[u] > 0.005)
      log.value("unit." + d.horizons[u - 1].name + "_" + d.horizons[u].name + ".missing_section_max",
                fmt_fixed(miss[u], 2) + " m at (" + fmt_fixed(mx[u], 0) + ", " + fmt_fixed(my[u], 0) + ")");
  }
  // Missing section per well column (at the well head; stage 2C, STX-TN-07 section 1).
  if (has_rel_)
    for (auto& w : d.wells)
      for (std::size_t u = 1; u < nh; ++u) {
        double m = missing_section(static_cast<int>(u), w.head.si(0), w.head.si(1));
        if (m > 0.005)
          log.value("well." + w.name + "." + d.horizons[u - 1].name + "_" + d.horizons[u].name + ".missing_section", fmt_fixed(m, 2));
      }
  for (std::size_t h = 0; h < nh; ++h)
    if (trend_[h]) {
      double lo = 1e300, hi = -1e300;
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
          double t = iso_thickness(h, dom.x0 + dom.sx * i / (n - 1), dom.y0 + dom.sy * j / (n - 1), true);
          lo = std::min(lo, t);
          hi = std::max(hi, t);
        }
      log.value(d.horizons[h].name + ".isochore_min", fmt_fixed(lo, 3));
      log.value(d.horizons[h].name + ".isochore_max", fmt_fixed(hi, 3));
      if (lo <= 0) log.note("horizon " + d.horizons[h].name + ": its isochore trend pinches out inside the model");
    }
}

double Structure::regional(double x, double y) const { return reg_gx_ * (x - reg_x0_) + reg_gy_ * (y - reg_y0_); }

double Structure::iso_thickness(std::size_t h, double x, double y, bool clamp) const {
  double t = iso_grid_[h] ? iso_grid_[h]->sample_clamped(x, y) : iso_const_[h];
  if (trend_[h] == 1) t += trend_a_[h] * (x - reg_x0_) + trend_b_[h] * (y - reg_y0_);
  else if (trend_[h] == 2) t += trend_a_[h] * std::hypot(x - trend_cx_[h], y - trend_cy_[h]);
  else return t;
  return clamp ? std::max(0.0, t) : t;
}

bool Structure::trend_clamped(std::size_t h, double x, double y) const {
  return trend_[h] != 0 && iso_thickness(h, x, y, false) <= 0;
}

double Structure::clino_offset(std::size_t h, double x, double y) const {
  const SeriesInfo& S = series_[static_cast<std::size_t>(hseries_[h])];
  const double i = static_cast<double>(S.last - h);  // 0 for the oldest surface
  const double u = ((x - S.bx) * S.ax + (y - S.by) * S.ay - i * S.rate) / S.width;
  const double p = S.sigmoid ? 0.5 * (1 - std::tanh(2 * u)) : std::min(1.0, std::max(0.0, 0.5 - u));
  return S.height * p;
}

double Structure::growth_fraction(std::size_t h, int onset, int end) const {
  if (static_cast<int>(h) >= onset) return 1.0;
  if (static_cast<int>(h) <= end) return 0.0;
  const double don = ref_depth_[static_cast<std::size_t>(onset)], den = ref_depth_[static_cast<std::size_t>(end)];
  return std::clamp((ref_depth_[h] - den) / (don - den), 0.0, 1.0);
}

double Structure::z_uncorrected(std::size_t h, double x, double y) const {
  switch (kind_[h]) {
    case Kind::Parametric: {
      double s = zref_[h] + regional(x, y);
      if (has_tilt_) s += tilt_gx_[h] * (x - reg_x0_) + tilt_gy_[h] * (y - reg_y0_);
      for (auto& p : prims_) s += p.contribution(h, x, y);
      return s;
    }
    case Kind::Surface: return surf_sign_[h] * surf_[h]->sample_clamped(x, y);
    case Kind::Isochore: return z_raw(static_cast<std::size_t>(iso_from_[h]), x, y) - iso_thickness(h, x, y, true);
    case Kind::Clinoform: return z_raw(series_[static_cast<std::size_t>(hseries_[h])].last + 1, x, y) + clino_offset(h, x, y);
  }
  return kNaN;
}

void Structure::apply_relations(double x, double y, double* z, unsigned char* flags, bool erosion_only) const {
  const std::size_t nh = kind_.size();
  for (std::size_t k = series_.size(); k-- > 0;) {  // oldest series first
    const SeriesInfo& S = series_[k];
    if (S.clino && !erosion_only && S.last + 1 < nh)
      for (std::size_t h = S.last + 1; h-- > S.first;) z[h] = z[S.last + 1] + clino_offset(h, x, y);
    if (S.base == 2 && !erosion_only && S.last + 1 < nh) {
      const double b = z[S.last + 1];
      for (std::size_t h = S.first; h <= S.last; ++h)
        if (z[h] < b) { z[h] = b; if (flags) flags[h] |= 2; }
    }
    if (S.base == 1) {
      const double u = z[S.last];
      for (std::size_t h = S.last + 1; h < nh; ++h)
        if (z[h] > u) { z[h] = u; if (flags) flags[h] |= 1; }
    }
  }
}

void Structure::z_detail(double x, double y, double* raw, double* eff, unsigned char* flags) const {
  const std::size_t nh = kind_.size();
  z_all_raw(x, y, raw);
  for (std::size_t h = 0; h < nh; ++h) { eff[h] = raw[h]; flags[h] = 0; }
  if (has_rel_) apply_relations(x, y, eff, flags, false);
}

bool Structure::absent(std::size_t h, double x, double y) const {
  if (!has_rel_) return false;
  const std::size_t nh = kind_.size();
  std::vector<double> raw(nh), eff(nh);
  std::vector<unsigned char> fl(nh);
  z_detail(x, y, raw.data(), eff.data(), fl.data());
  return fl[h] != 0;
}

double Structure::correction(std::size_t h, double x, double y) const {
  const auto& c = corr_[h];
  double s = 0;
  for (std::size_t i = 0; i < c.w.size(); ++i) {
    double dd = std::hypot(x - c.x[i], y - c.y[i]) / corr_range_;
    s += c.w[i] * std::exp(-dd * dd);
  }
  return s;
}

double Structure::z_raw(std::size_t h, double x, double y) const {
  double v = z_uncorrected(h, x, y);
  return corr_[h].w.empty() ? v : v + correction(h, x, y);
}

double Structure::z(std::size_t h, double x, double y) const {
  if (!has_rel_) return z_raw(h, x, y);
  double zs[64];
  std::vector<double> big;
  double* zz = zs;
  if (kind_.size() > 64) { big.resize(kind_.size()); zz = big.data(); }
  z_all(x, y, zz);
  return zz[h];
}

void Structure::z_all(double x, double y, double* out) const {
  z_all_raw(x, y, out);
  if (has_rel_) apply_relations(x, y, out, nullptr, false);
}

void Structure::z_all_raw(double x, double y, double* out) const {
  if (!optimised_kernels()) {  // reference: each horizon evaluated independently
    for (std::size_t h = 0; h < kind_.size(); ++h) out[h] = z_raw(h, x, y);
    return;
  }
  // Optimised: shape, tilt and rim of each primitive once per (x, y); the per-horizon arithmetic
  // is the same expression as FoldPrimitive::contribution, so results are bitwise identical.
  const std::size_t np = prims_.size();
  double shp[16], tlt[16], rm[16];
  bool small = np <= 16;
  if (small)
    for (std::size_t k = 0; k < np; ++k) {
      shp[k] = prims_[k].shape(x, y);
      tlt[k] = prims_[k].tilt(x, y);
      rm[k] = (prims_[k].rim_depth != 0 && prims_[k].relief_max > 0) ? prims_[k].rim(x, y) : 0.0;
    }
  for (std::size_t h = 0; h < kind_.size(); ++h) {
    if (kind_[h] != Kind::Parametric || !small) { out[h] = z_raw(h, x, y); continue; }
    double sacc = zref_[h] + regional(x, y);
    if (has_tilt_) sacc += tilt_gx_[h] * (x - reg_x0_) + tilt_gy_[h] * (y - reg_y0_);
    for (std::size_t k = 0; k < np; ++k) {
      const auto& p = prims_[k];
      if (h < p.offset.size() && p.offset[h]) { sacc += p.offset[h]->eval(p, x, y); continue; }
      double r = p.relief[h];
      double c = r * shp[k] + tlt[k];
      if (p.rim_depth != 0 && p.relief_max > 0) c -= p.rim_depth * (std::fabs(r) / p.relief_max) * rm[k] * (r >= 0 ? 1 : -1);
      sacc += c;
    }
    out[h] = corr_[h].w.empty() ? sacc : sacc + correction(h, x, y);
  }
}

}  // namespace stratonix
