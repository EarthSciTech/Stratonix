// Stage 4D (STX-TN-17 section 1, decision 24): geostatistical properties of the description on the
// corner-point grid. Each unit's cells are mapped to geological-time coordinates (u, v, w) - w the
// relative stratigraphic position scaled by the unit's mean thickness - so ranges are depositional
// distances; conditioning data are the upscaled logs of stage 4C.
#include "stratonix/package/geostat_props.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>

#include "stratonix/core/error.hpp"
#include "stratonix/core/units.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/repr/geotime.hpp"

namespace stratonix {
namespace {

constexpr double kPiP = 3.14159265358979323846;

geo::Variogram to_variogram(const desc::GridDesc::VariogramDesc& v, double scale = 1.0) {
  geo::Variogram g;
  g.nugget = (v.nugget.present ? v.nugget.si() : 0.0) * scale;
  for (auto& st : v.structures) {
    geo::VariogramStructure s;
    s.type = st.type;
    s.contribution = st.contribution.si() * scale;
    s.range_major = st.ranges.si(0);
    s.range_minor = st.ranges.values.size() > 1 ? st.ranges.si(1) : s.range_major;
    s.range_vertical = st.ranges.values.size() > 2 ? st.ranges.si(2) : s.range_major;
    s.azimuth = st.azimuth.present ? st.azimuth.values.at(0) * (st.azimuth.to_si / (kPiP / 180.0)) : 0.0;  // degrees
    g.structures.push_back(s);
  }
  return g;
}
std::string g4(double x) {
  char b[32];
  std::snprintf(b, sizeof b, "%.4g", x);
  return b;
}
std::string describe(const geo::Variogram& v) {
  std::string t = "nugget " + g4(v.nugget);
  for (auto& s : v.structures)
    t += " + " + g4(s.contribution) + " " + s.type + " (" + fmt_fixed(s.range_major, 0) + ", " + fmt_fixed(s.range_minor, 0) + ", " +
         fmt_fixed(s.range_vertical, 1) + " m, azimuth " + fmt_fixed(s.azimuth, 0) + " deg)";
  return t;
}
double norm_q(double p) {
  // Acklam with one Halley step (as stage 3C); only used to make a stated target distribution
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
  double x;
  if (p < 0.02425) {
    const double q = std::sqrt(-2 * std::log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else if (p > 1 - 0.02425) {
    const double q = std::sqrt(-2 * std::log(1 - p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else {
    const double q = p - 0.5, r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }
  return x;
}
/// Kolmogorov-Smirnov distance between two samples.
double ks_distance(std::vector<double> a, std::vector<double> b) {
  if (a.empty() || b.empty()) return 0;
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());
  std::size_t i = 0, j = 0;
  double dmax = 0;
  while (i < a.size() && j < b.size()) {
    const double x = std::min(a[i], b[j]);
    while (i < a.size() && a[i] <= x) ++i;
    while (j < b.size() && b[j] <= x) ++j;
    dmax = std::max(dmax, std::fabs(static_cast<double>(i) / a.size() - static_cast<double>(j) / b.size()));
  }
  return dmax;
}

}  // namespace

GeostatBuild geostat_properties(const desc::Description& d, const Model& m, const repr::CornerPointGrid& g, const repr::PropertySet& existing) {
  GeostatBuild out;
  if (d.grid.geostatistics.empty()) return out;
  out.present = true;
  const auto& info = m.info();
  // geological time: horizon ages when every horizon states one, else ordinals
  std::vector<double> ages;
  bool all_ages = !d.horizons.empty();
  for (auto& h : d.horizons) all_ages = all_ages && h.age.present;
  if (all_ages)
    for (auto& h : d.horizons) ages.push_back(h.age.si());
  GeologicalTime gt(m, all_ages ? ages : std::vector<double>{});
  const auto& tages = gt.ages();
  std::string rep = "Stratonix geostatistical properties (stage 4D, STX-TN-17 section 1; decision 24)\n";
  rep += std::string("Cells are simulated in geological-time coordinates (u, v, w): w is the relative stratigraphic position scaled by the unit's mean thickness; ") +
         (all_ages ? "horizon ages" : "horizon ordinals (no ages stated)") + " give the time axis.\n";
  std::string csv = "property,unit,direction,lag_m,gamma,pairs,model_gamma\n";

  for (const auto& q : d.grid.geostatistics) {
    int ui = -1;
    for (std::size_t u = 0; u < info.units.size(); ++u)
      if (info.units[u].name == q.unit) ui = static_cast<int>(u);
    if (ui < 0) throw Error(ErrorCode::Build, "grid.geostatistics " + q.name + ": unit '" + q.unit + "' is not in the model");
    const int ht = info.units[static_cast<std::size_t>(ui)].top_horizon, hb = info.units[static_cast<std::size_t>(ui)].base_horizon;
    // the unit's cells and its mean thickness per column
    std::vector<std::size_t> cells;
    double hsum = 0;
    std::size_t hn = 0, nlay = 1;
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i) {
        double zt = -1e300, zb = 1e300;
        std::size_t lay = 0;
        for (int k = 0; k < g.nz; ++k) {
          const auto c = g.cell(i, j, k);
          if (!g.actnum[c] || g.zone[c] != ui) continue;
          cells.push_back(c);
          ++lay;
          const auto cn = g.corners(i, j, k);
          zt = std::max(zt, 0.25 * (cn[0].z + cn[1].z + cn[2].z + cn[3].z));
          zb = std::min(zb, 0.25 * (cn[4].z + cn[5].z + cn[6].z + cn[7].z));
        }
        if (zt > zb) { hsum += zt - zb; ++hn; }
        nlay = std::max(nlay, lay);
      }
    std::sort(cells.begin(), cells.end());
    const double H = hn ? hsum / static_cast<double>(hn) : 1.0;
    const double tt = ht >= 0 && static_cast<std::size_t>(ht) < tages.size() ? tages[static_cast<std::size_t>(ht)] : kNaN;
    const double tb = hb >= 0 && static_cast<std::size_t>(hb) < tages.size() ? tages[static_cast<std::size_t>(hb)] : kNaN;
    std::vector<Vec3> nodes(cells.size());
    std::size_t fallback = 0;
    for (std::size_t n = 0; n < cells.size(); ++n) {
      const std::size_t c = cells[n];
      const int i = static_cast<int>(c % static_cast<std::size_t>(g.nx)), j = static_cast<int>((c / static_cast<std::size_t>(g.nx)) % static_cast<std::size_t>(g.ny)),
                k = static_cast<int>(c / (static_cast<std::size_t>(g.nx) * static_cast<std::size_t>(g.ny)));
      const Vec3 x = repr::cell_centre(g.corners(i, j, k));
      const GeoTimePoint p = gt.forward(x);
      if (p.defined && std::isfinite(p.u) && std::isfinite(p.v) && std::isfinite(p.t) && std::isfinite(tt) && std::isfinite(tb) && tb != tt) {
        nodes[n] = {p.u, p.v, -H * (p.t - tt) / (tb - tt)};
      } else {
        nodes[n] = x;  // present-day coordinates (counted and reported)
        ++fallback;
      }
    }
    out.cells += cells.size();
    out.fallback += fallback;
    // conditioning data: the upscaled log's values in the unit's cells
    const repr::CellProperty* cond = q.condition_to.empty() ? nullptr : existing.find(q.condition_to);
    if (!q.condition_to.empty() && !cond)
      throw Error(ErrorCode::Build, "grid.geostatistics " + q.name + ": the conditioning property '" + q.condition_to + "' is not on the grid");
    std::vector<std::pair<std::size_t, double>> data;  // node, value
    if (cond)
      for (std::size_t n = 0; n < cells.size(); ++n)
        if (cells[n] < cond->v.size() && std::isfinite(cond->v[cells[n]])) data.push_back({n, cond->v[cells[n]]});
    geo::SearchOptions so;
    if (q.search_radii.present) {
      so.radius_major = q.search_radii.si(0);
      so.radius_minor = q.search_radii.values.size() > 1 ? q.search_radii.si(1) : so.radius_major;
      so.radius_vertical = q.search_radii.values.size() > 2 ? q.search_radii.si(2) : so.radius_major;
      so.azimuth = q.search_azimuth.present ? q.search_azimuth.values.at(0) * (q.search_azimuth.to_si / (kPiP / 180.0)) : 0.0;
      so.max_data = static_cast<std::size_t>(std::llround(q.max_data.si()));
    }
    const std::uint64_t seed = static_cast<std::uint64_t>(std::llround(q.seed.si()));
    repr::CellProperty prop;
    prop.name = q.name;
    prop.unit = q.property_unit;
    prop.v.assign(g.ncell(), kNaN);
    rep += "\n" + q.name + " (" + q.kind + ") in " + q.unit + ": " + std::to_string(cells.size()) + " cells, mean thickness " + fmt_fixed(H, 1) + " m, " +
           std::to_string(fallback) + " cells in present-day coordinates; seed " + std::to_string(seed) + "\n";
    // nodes for the experimental variograms: at most 2500, every k-th (deterministic)
    std::vector<std::size_t> sub;
    const std::size_t stride = std::max<std::size_t>(1, cells.size() / 2500);
    for (std::size_t n = 0; n < cells.size(); n += stride) sub.push_back(n);
    auto add_variograms = [&](const std::vector<double>& vals, const geo::Variogram& model, double az) {
      std::vector<Vec3> pts;
      std::vector<double> vv;
      for (auto n : sub) { pts.push_back(nodes[n]); vv.push_back(vals[n]); }
      double rmax = 0, rv = 0;
      for (auto& s : model.structures) { rmax = std::max(rmax, s.range_major); rv = std::max(rv, s.range_vertical); }
      for (int dir = 0; dir < 2; ++dir) {
        geo::VariogramDirection vd;
        vd.azimuth = az;
        vd.dip = dir == 1 ? 90 : 0;
        vd.angle_tolerance = dir == 1 ? 10 : 22.5;
        vd.bandwidth = dir == 1 ? 0.5 * std::max(g.dx, g.dy) : 0.25 * rmax;
        vd.lags = 12;
        vd.lag = (dir == 1 ? 1.5 * rv : 1.5 * rmax) / 12.0;
        vd.lag_tolerance = 0.5 * vd.lag;
        auto e = geo::experimental_variogram(pts, vv, vd);
        GeostatBuild::VariogramPlot vp{q.name, q.unit, dir == 1 ? "vertical" : "horizontal (major)", e.distance, e.gamma, e.pairs, model, az, dir == 1};
        for (std::size_t l = 0; l < e.distance.size(); ++l) {
          const double h = e.distance[l], a = az * kPiP / 180.0;
          const double mg = model.gamma(dir == 1 ? Vec3{0, 0, h} : Vec3{h * std::sin(a), h * std::cos(a), 0});
          csv += q.name + "," + q.unit + "," + vp.direction + "," + fmt_fixed(h, 2) + "," + (std::isfinite(e.gamma[l]) ? fmt_fixed(e.gamma[l], 5) : std::string()) +
                 "," + std::to_string(e.pairs[l]) + "," + fmt_fixed(mg, 5) + "\n";
        }
        out.variograms.push_back(std::move(vp));
      }
    };

    if (q.kind == "sgs") {
      // the distribution's parameters in the property's own unit (the conditioning log's unit,
      // else property_unit): a log in g/cm3 is never compared with limits converted to kg/m3
      const std::string vunit = cond && !cond->unit.empty() ? cond->unit : q.property_unit;
      const auto ud = find_unit(vunit);
      auto in_unit = [&](const desc::Param& p) {
        if (ud && (p.dim == ud->dim || p.dim == Dimension::Any)) return p.si() / ud->to_si;
        return p.values.at(0);  // no known unit: as written
      };
      const double zmin = in_unit(q.min), zmax = in_unit(q.max);
      // target distribution: the conditioning data (declustered when a cell is stated), else stated
      std::vector<double> tv, tw;
      if (!data.empty()) {
        std::vector<Vec3> dp;
        for (auto& [n, v] : data) { tv.push_back(v); dp.push_back(nodes[n]); }
        tw = q.decluster_cell.present ? geo::decluster_weights(dp, q.decluster_cell.si()) : std::vector<double>(tv.size(), 1.0);
      } else {
        for (int i = 0; i < 1000; ++i) tv.push_back(std::clamp(in_unit(q.mean) + in_unit(q.sd) * norm_q((i + 0.5) / 1000.0), zmin, zmax));
        tw.assign(tv.size(), 1.0);
      }
      geo::NormalScore ns(tv, tw, zmin, zmax);
      geo::Variogram v = to_variogram(q.variogram);
      // vertical range from the wells (the usual practice: vertical structure from logs)
      std::vector<geo::ConditioningPoint> cp;
      for (auto& [n, val] : data) cp.push_back({n, ns.forward(val)});
      if (q.variogram.fit_vertical && cp.size() >= 10) {
        std::vector<Vec3> pts;
        std::vector<double> vals;
        for (auto& c : cp) { pts.push_back(nodes[c.node]); vals.push_back(c.value); }
        double rv = 0;
        for (auto& s : v.structures) rv = std::max(rv, s.range_vertical);
        geo::VariogramDirection vd;
        vd.dip = 90;
        vd.angle_tolerance = 10;
        vd.bandwidth = 0.5 * std::max(g.dx, g.dy);
        vd.lag = H / static_cast<double>(nlay);  // the layer thickness: the wells' vertical sampling on the grid
        vd.lags = std::max<std::size_t>(2, nlay - 1);
        vd.lag_tolerance = 0.5 * vd.lag;
        auto e = geo::experimental_variogram(pts, vals, vd);
        std::size_t usable = 0;
        for (std::size_t l = 0; l < e.gamma.size(); ++l) usable += e.pairs[l] >= 5 && std::isfinite(e.gamma[l]);
        // one factor on every structure's vertical range, by weighted least squares over a fixed grid
        double best = 1, bc = 1e300;
        for (int s = -40; s <= 40; ++s) {
          const double f = std::pow(2.0, s / 20.0);
          double cst = 0;
          for (std::size_t l = 0; l < e.gamma.size(); ++l) {
            if (e.pairs[l] < 5 || !std::isfinite(e.gamma[l])) continue;
            geo::Variogram t = v;
            for (auto& st : t.structures) st.range_vertical *= f;
            const double gm = std::max(1e-6, t.gamma({0, 0, e.distance[l]}));
            cst += e.pairs[l] * (e.gamma[l] - gm) * (e.gamma[l] - gm) / (gm * gm);
          }
          if (cst < bc) { bc = cst; best = f; }
        }
        const bool bound = best <= std::pow(2.0, -2.0) + 1e-12 || best >= std::pow(2.0, 2.0) - 1e-12;
        if (usable < 3 || bound) {
          rep += "  vertical range not fitted: " + std::string(usable < 3 ? "fewer than 3 lags with 5 pairs in the wells" : "the fit reached the search bound (x0.25 to x4)") +
                 "; the stated range is used\n";
        } else {
          for (auto& st : v.structures) st.range_vertical *= best;
          rep += "  vertical range fitted to the wells (calibrated): x" + fmt_fixed(best, 3) + " of the stated value\n";
        }
        out.variograms.push_back({q.name, q.unit, "wells (vertical)", e.distance, e.gamma, e.pairs, v, 0, true});
      }
      std::vector<double> lm;
      if (!q.trend.empty()) {
        const repr::CellProperty* tr = existing.find(q.trend);
        if (!tr) throw Error(ErrorCode::Build, "grid.geostatistics " + q.name + ": the trend property '" + q.trend + "' is not on the grid");
        lm.resize(cells.size());
        for (std::size_t n = 0; n < cells.size(); ++n) lm[n] = std::isfinite(tr->v[cells[n]]) ? ns.forward(tr->v[cells[n]]) : 0.0;
        rep += "  trend: " + q.trend + " as the local mean (simple kriging with a varying mean)\n";
      }
      const auto y = geo::sgs(nodes, cp, v, so, seed, lm);
      std::vector<double> z(cells.size());
      for (std::size_t n = 0; n < cells.size(); ++n) z[n] = ns.back(y[n]);
      for (auto& [n, val] : data) z[n] = val;  // conditioning values exactly at their cells
      for (std::size_t n = 0; n < cells.size(); ++n) prop.v[cells[n]] = z[n];
      prop.evidence = std::string("simulated (sequential Gaussian") + (data.empty() ? ", unconditioned" : ", conditioned to " + q.condition_to) + ")";
      prop.source = "SGS seed " + std::to_string(seed) + ", normal-score variogram " + describe(v) + ", range [" + fmt_shortest(zmin) + ", " + fmt_shortest(zmax) + "]";
      rep += "  variogram (normal scores): " + describe(v) + "\n  conditioning: " + std::to_string(data.size()) + " cells of " +
             (q.condition_to.empty() ? std::string("(none)") : q.condition_to) + "\n";
      // conditioning and statistics
      GeostatBuild::Conditioning c{q.name, data.size(), 0, 0, true};
      for (auto& [n, val] : data) {
        const double e = std::fabs(prop.v[cells[n]] - val);
        c.largest = std::max(c.largest, e);
        c.mismatched += e > 1e-9 * (1 + std::fabs(val));
      }
      out.conditioning.push_back(c);
      const double ks = ks_distance(z, tv);
      out.statistics.push_back({q.name, "Kolmogorov-Smirnov distance to the target distribution", 0.0, ks, 0.1});
      rep += "  histogram: KS distance " + fmt_fixed(ks, 3) + " to the " + (data.empty() ? "stated" : "conditioning") + " distribution\n";
      add_variograms(y, v, v.structures.empty() ? 0.0 : v.structures[0].azimuth);
      GeostatBuild::Histogram hg{q.name, q.unit, q.property_unit, z, {}, {}};
      for (auto& [n, val] : data) hg.wells.push_back(val);
      out.histograms.push_back(std::move(hg));
    } else if (q.kind == "sis") {
      std::vector<double> props;
      std::vector<geo::Variogram> vs;
      std::map<int, int> code_index;
      std::vector<std::string> names;
      for (std::size_t k = 0; k < q.codes.size(); ++k) {
        const double p = q.codes[k].proportion.si();
        props.push_back(p);
        vs.push_back(to_variogram(q.codes[k].variogram, p * (1 - p)));  // indicator variogram sill p (1 - p)
        code_index[q.codes[k].code] = static_cast<int>(k);
        names.push_back(std::to_string(q.codes[k].code) + " " + q.codes[k].name);
      }
      std::vector<std::pair<std::size_t, int>> cd;
      std::size_t unknown = 0;
      for (auto& [n, val] : data) {
        auto it = code_index.find(static_cast<int>(std::lround(val)));
        if (it == code_index.end()) { ++unknown; continue; }
        cd.push_back({n, it->second});
      }
      const auto f = geo::sis(nodes, cd, props, vs, so, seed);
      std::vector<double> got(q.codes.size(), 0.0);
      for (std::size_t n = 0; n < cells.size(); ++n) {
        prop.v[cells[n]] = q.codes[static_cast<std::size_t>(f[n])].code;
        got[static_cast<std::size_t>(f[n])] += 1.0 / static_cast<double>(std::max<std::size_t>(1, cells.size()));
      }
      prop.integer = true;
      prop.evidence = std::string("simulated (sequential indicator") + (cd.empty() ? ", unconditioned" : ", conditioned to " + q.condition_to) + ")";
      prop.source = "SIS seed " + std::to_string(seed) + " with the stated proportions and indicator variograms";
      rep += "  conditioning: " + std::to_string(cd.size()) + " cells" + (unknown ? " (" + std::to_string(unknown) + " with codes not in the list, left out)" : "") + "\n";
      GeostatBuild::Conditioning c{q.name, cd.size(), 0, 0, true};
      for (auto& [n, k] : cd) c.mismatched += f[n] != k;
      out.conditioning.push_back(c);
      for (std::size_t k = 0; k < q.codes.size(); ++k) {
        out.statistics.push_back({q.name, "proportion of " + names[k], props[k], got[k], 0.05});
        rep += "  " + names[k] + ": proportion " + fmt_fixed(got[k], 3) + " (target " + fmt_fixed(props[k], 3) + "), variogram " + describe(vs[k]) + "\n";
      }
      std::vector<double> ind(cells.size());
      for (std::size_t n = 0; n < cells.size(); ++n) ind[n] = f[n] == 0 ? 1.0 : 0.0;
      add_variograms(ind, vs[0], vs[0].structures.empty() ? 0.0 : vs[0].structures[0].azimuth);
      GeostatBuild::Histogram hg{q.name, q.unit, "code", {}, {}, names};
      for (std::size_t n = 0; n < cells.size(); ++n) hg.simulated.push_back(prop.v[cells[n]]);
      for (auto& [n, val] : data) hg.wells.push_back(std::lround(val));
      out.histograms.push_back(std::move(hg));
    } else if (q.kind == "objects") {
      geo::ObjectSpec o;
      o.kind = q.object_kind;
      o.code = q.object_code;
      const desc::GeobodyDesc* gb = nullptr;
      for (auto& b : d.geobodies)
        if (b.name == q.geobody) gb = &b;
      auto pick = [&](const desc::Param& own, const desc::Param* fromgb, double dflt_none) {
        if (own.present) return own.si();
        if (fromgb && fromgb->present) return fromgb->si();
        return dflt_none;
      };
      o.width = pick(q.width, gb ? &gb->width : nullptr, 0.0);
      o.thickness = pick(q.thickness, gb ? &gb->thickness : nullptr, 0.0);
      o.wavelength = pick(q.wavelength, gb ? &gb->wavelength : nullptr, 0.0);
      o.radius = pick(q.radius, gb ? &gb->radius : nullptr, 0.0);
      o.sinuosity = q.sinuosity.present ? q.sinuosity.si() : 1.0;
      o.azimuth = q.azimuth.present ? q.azimuth.values.at(0) * (q.azimuth.to_si / (kPiP / 180.0))
                                    : (gb && gb->orientation.present ? gb->orientation.values.at(0) * (gb->orientation.to_si / (kPiP / 180.0)) : 0.0);
      o.proportion = q.proportion.si();
      Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
      for (auto& p : nodes) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
      }
      std::vector<std::pair<std::size_t, int>> wells;
      for (auto& [n, val] : data) wells.push_back({n, static_cast<int>(std::lround(val))});
      const std::size_t iters = q.max_iterations.present ? static_cast<std::size_t>(std::llround(q.max_iterations.si())) : 20000;
      const auto r = geo::object_facies(nodes, lo, hi, o, wells, seed, iters);
      for (std::size_t n = 0; n < cells.size(); ++n) prop.v[cells[n]] = r.facies[n];
      prop.integer = true;
      prop.evidence = std::string("simulated (object-based ") + o.kind + "s" + (wells.empty() ? "" : ", conditioned to " + q.condition_to) + ")";
      prop.source = "objects seed " + std::to_string(seed) + ": " + o.kind + " thickness " + fmt_shortest(o.thickness) + " m" +
                    (o.kind == "channel" ? ", width " + fmt_shortest(o.width) + " m, sinuosity " + fmt_shortest(o.sinuosity) : ", radius " + fmt_shortest(o.radius) + " m") +
                    (gb ? " (from geobody " + gb->name + " where not stated)" : "");
      rep += "  " + std::to_string(r.objects) + " objects, proportion " + fmt_fixed(r.proportion, 3) + " (target " + fmt_fixed(o.proportion, 3) + "); wells: " +
             std::to_string(r.well_cells) + " cells, " + std::to_string(r.well_misfit) + " not honoured\n";
      out.conditioning.push_back({q.name, r.well_cells, r.well_misfit, 0, false});
      out.statistics.push_back({q.name, "proportion of code " + std::to_string(o.code), o.proportion, r.proportion, 0.05});
      GeostatBuild::Histogram hg{q.name, q.unit, "code", {}, {}, {"0 background", std::to_string(o.code) + " " + o.kind}};
      for (std::size_t n = 0; n < cells.size(); ++n) hg.simulated.push_back(r.facies[n]);
      for (auto& [n, val] : data) hg.wells.push_back(std::lround(val));
      out.histograms.push_back(std::move(hg));
    } else {
      throw Error(ErrorCode::Build, "grid.geostatistics " + q.name + ": unknown kind '" + q.kind + "'");
    }
    for (double v : prop.v) prop.missing += std::isfinite(v) ? 0 : 1;
    if (prop.missing) prop.missing_reason = "inactive or outside unit " + q.unit;
    out.props.push_back(std::move(prop));
  }
  if (out.fallback) rep += "\n" + std::to_string(out.fallback) + " of " + std::to_string(out.cells) + " cells were simulated in present-day coordinates (geological time undefined there)\n";
  out.report = rep;
  out.csv = csv;
  return out;
}

}  // namespace stratonix
