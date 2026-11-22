// Stage 4C (STX-TN-16 section 1): well logs, surveys, tops and ages in a build - logs read as
// delivered, placed on the well path, upscaled to the grid, correlated by age, checked (C35, C38)
// and drawn (well log and correlation panels).
#include "stratonix/package/welldata.hpp"

#include <algorithm>
#include <cmath>
#include <map>

#include "stratonix/core/error.hpp"
#include "stratonix/core/units.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/ages.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/io/welldata.hpp"
#include "stratonix/repr/upscale.hpp"

namespace stratonix {

namespace {
std::string lower(std::string s) {
  for (auto& ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  return s;
}
/// Length unit of a log index (m, ft, ...) to metres; refused when the file states none.
double index_scale(const io::WellLogSet& s, const std::string& file) {
  std::string u = lower(s.index.unit);
  if (u == "f" || u == "feet") u = "ft";
  if (u == "meter" || u == "meters" || u == "metres") u = "m";
  const auto def = find_unit(u);
  if (!def || def->dim != Dimension::Length)
    throw Error(ErrorCode::Validation, file + ": the index unit '" + s.index.unit + "' is not a length unit (time-indexed logs are read with Prompt 56)",
                "wells.logs");
  return def->to_si;
}
AgeRule age_rule(const std::string& r) { return r == "constant_rate_per_unit" ? AgeRule::ConstantRatePerUnit : AgeRule::LinearInDepth; }
std::string rule_text(const std::string& r) { return r == "constant_rate_per_unit" ? "constant rate per unit" : "linear in depth between dated samples"; }
repr::UpscaleMethod upscale_method(const std::string& m) {
  using M = repr::UpscaleMethod;
  return m == "harmonic" ? M::Harmonic : m == "geometric" ? M::Geometric : m == "rms" ? M::Rms : m == "minimum" ? M::Minimum
       : m == "maximum" ? M::Maximum : m == "most_frequent" ? M::MostFrequent : M::Arithmetic;
}
}  // namespace

WellDataBuild well_data(const desc::Description& d, const Model& m, const repr::CornerPointGrid* grid) {
  WellDataBuild out;
  bool any = !d.correlation_rule.empty();
  for (auto& w : d.wells) any = any || !w.logs.empty() || !w.ages.empty() || !w.survey_method.empty();
  if (!any) return out;
  out.present = true;
  std::string rep = "Stratonix well data (stage 4C, STX-TN-16 section 1)\n";
  struct PerWell { WellPath path; std::vector<io::WellLogSet> logs; std::vector<double> md, tvd; std::vector<Vec3> p; };
  std::vector<PerWell> pw;
  for (auto& w : d.wells) {
    PerWell P{make_well_path(w, d), {}, {}, {}, {}};
    CheckContext::WellDataResults::Well R;
    R.well = w.name;
    R.survey = !w.survey_method.empty() ? "survey by " + w.survey_method + " (" + w.azimuth_reference + ")"
             : w.stations.present ? "trajectory as written in v3 (minimum curvature, grid azimuths)" : "vertical well";
    rep += "\nwell " + w.name + ": " + R.survey + "\n";
    const double td_md = P.path.total_depth_md();
    for (auto& ld : w.logs) {
      const std::string ext = lower(ld.resolved.extension().string());
      std::vector<io::WellLogSet> sets;
      if (ext == ".dlis" || ext == ".dls") sets = io::read_dlis(ld.resolved);
      else sets.push_back(io::read_las(ld.resolved));
      const io::WellLogSet* s = nullptr;
      for (auto& x : sets)
        if (ld.frame.empty() || x.source.find("frame " + ld.frame + ")") != std::string::npos) { s = &x; break; }
      if (!s) throw Error(ErrorCode::Validation, ld.file + ": no frame '" + ld.frame + "'", "wells.logs");
      std::string kind = ld.index.empty() ? s->index_type : (ld.index == "tvd" ? "TVD" : "MD");
      if (kind != "MD" && kind != "TVD")
        throw Error(ErrorCode::Validation, ld.file + ": the index type of '" + s->index.mnemonic + "' is not stated; write index: md or index: tvd",
                    "wells.logs");
      const double sc = index_scale(*s, ld.file);
      io::WellLogSet sel = *s;
      sel.curves.clear();
      if (ld.curves.empty()) sel.curves = s->curves;
      for (auto& cn : ld.curves) {
        const auto* c = s->find(cn);
        if (!c || c == &s->index) throw Error(ErrorCode::Validation, ld.file + ": no curve '" + cn + "'", "wells.logs");
        sel.curves.push_back(*c);
      }
      for (auto& v : sel.index.values) v *= sc;
      sel.index_type = kind;
      rep += "  " + ld.file + " (" + s->version + ", index " + s->index.mnemonic + " [" + s->index.unit + "] as " + kind + ")\n";
      for (auto& line : s->report) rep += "    " + line + "\n";
      double last = 0;
      for (double v : sel.index.values)
        if (std::isfinite(v)) last = std::max(last, kind == "MD" ? v : P.path.md_at_elevation(P.path.datum_z() - v));
      R.beyond = std::max(R.beyond, last - td_md);
      for (auto& c : sel.curves) R.notes.push_back(c.mnemonic + " [" + (c.unit.empty() ? "no unit" : c.unit) + "] " + std::to_string(c.nulls) + " nulls");
      R.curves += sel.curves.size();
      ++R.logs;
      P.logs.push_back(std::move(sel));
    }
    for (auto& [h, dp] : w.tops) {
      const double bottom = -P.path.points().back().z;
      R.tops_off_path = std::max(R.tops_off_path, dp.si() - bottom);
    }
    out.results.wells.push_back(R);
    pw.push_back(std::move(P));
  }
  // logs along the paths: csv, panels
  for (std::size_t k = 0; k < d.wells.size(); ++k) {
    const auto& w = d.wells[k];
    auto& P = pw[k];
    // units of the model along the well, by MD
    std::vector<fig::WellPanel::Interval> units;
    {
      const double td = P.path.total_depth_md(), step = std::max(1.0, td / 2000);
      int cur = -2;
      for (double md = 0; md <= td + 1e-9; md += step) {
        const int u = m.unit_at(P.path.at_md(md));
        if (u != cur) {
          if (!units.empty()) units.back().base = md;
          units.push_back({u, md, md});
          cur = u;
        }
      }
      if (!units.empty()) units.back().base = td;
    }
    for (std::size_t q = 0; q < P.logs.size(); ++q) {
      const auto& s = P.logs[q];
      fig::WellPanel pan;
      pan.well = w.name + (P.logs.size() > 1 ? " (" + w.logs[q].file + ")" : "");
      pan.index_label = "MD (m)";
      std::string csv = "md_m,tvd_m,x,y,z";
      for (auto& c : s.curves) csv += "," + c.mnemonic + (c.unit.empty() ? "" : " [" + c.unit + "]");
      csv += "\n";
      for (std::size_t i = 0; i < s.index.values.size(); ++i) {
        const double v = s.index.values[i];
        const double md = s.index_type == "MD" ? v : P.path.md_at_elevation(P.path.datum_z() - v);
        const Vec3 p = s.index_type == "MD" ? P.path.at_md(v) : P.path.at_elevation(P.path.datum_z() - v);
        pan.index.push_back(md);
        csv += fmt_fixed(md, 3) + "," + fmt_fixed(-p.z, 3) + "," + fmt_fixed(p.x, 2) + "," + fmt_fixed(p.y, 2) + "," + fmt_fixed(p.z, 3);
        for (auto& c : s.curves) csv += "," + (i < c.values.size() && std::isfinite(c.values[i]) ? fmt_shortest(c.values[i]) : std::string());
        csv += "\n";
        if (q == 0) { P.md.push_back(md); P.tvd.push_back(-p.z); P.p.push_back(p); }
      }
      for (auto& c : s.curves) pan.tracks.push_back({c.mnemonic, c.unit, c.values});
      for (auto& [h, dp] : w.tops) pan.tops.push_back({h, P.path.md_at_elevation(-dp.si())});
      for (auto& a : w.ages) pan.ages.push_back({a.name, P.path.md_at_elevation(-a.depth.si()), a.age.si()});
      pan.units = units;
      out.panels.push_back(pan);
      out.csv["well_logs_" + w.name + (P.logs.size() > 1 ? "_" + std::to_string(q + 1) : "") + ".csv"] = csv;
    }
  }
  // ages and correlation (C38)
  std::vector<WellAges> wa;
  for (auto& w : d.wells) {
    WellAges a{w.name, {}};
    for (auto& x : w.ages)
      a.data.push_back({x.kind, x.name, x.depth.si(), x.depth_base.present ? x.depth_base.si() : -1.0, x.age.si(), x.sigma.present ? x.sigma.si() : 0.0,
                        x.timescale, x.source});
    out.results.dated += a.data.size();
    wa.push_back(a);
  }
  std::vector<std::pair<std::string, double>> hages;
  for (auto& h : d.horizons)
    if (h.age.present) hages.push_back({h.name, h.age.si()});
  Correlation cor;
  const bool dated = out.results.dated > 0;
  if (dated) {
    out.results.ages = true;
    cor = correlate_by_age(wa, hages, age_rule(d.correlation_rule));
    for (auto& cf : cor.conflicts)
      if (cf.what.find("inversion") != std::string::npos) out.results.conflicts.push_back({cf.well, cf.what});
    // tops against the age data at their depths (2 sigma of the bracketing data)
    for (std::size_t k = 0; k < d.wells.size(); ++k)
      for (auto& [h, dp] : d.wells[k].tops) {
        const int hi = d.horizon_index(h);
        if (hi < 0 || !d.horizons[static_cast<std::size_t>(hi)].age.present || wa[k].data.size() < 2) continue;
        const double a = age_at_depth(wa[k].data, dp.si(), age_rule(d.correlation_rule));
        if (!std::isfinite(a)) continue;
        double sg = 0;
        for (auto& x : wa[k].data) sg = std::max(sg, x.sigma);
        const double ah = d.horizons[static_cast<std::size_t>(hi)].age.si();
        if (std::fabs(a - ah) > 2 * sg + 1e-9)
          out.results.conflicts.push_back({d.wells[k].name, "top " + h + " at " + fmt_fixed(dp.si(), 1) + " m dates " + fmt_fixed(a, 2) + " Ma; the horizon is " +
                                                                fmt_fixed(ah, 2) + " Ma (beyond 2 sigma, " + fmt_fixed(2 * sg, 2) + " Ma)"});
      }
  }
  // correlation panel: ages where the horizon is dated, else tops, else the model
  std::size_t with = 0;
  for (auto& w : d.wells) with += (!w.tops.empty() || !w.ages.empty() || !w.logs.empty()) ? 1 : 0;
  if (with >= 2) {
    out.correlation = true;
    std::string cc = "horizon,well,tvd_m,basis\n";
    for (std::size_t h = 0; h < d.horizons.size(); ++h) {
      const auto& H = d.horizons[h];
      std::vector<double> row(d.wells.size(), kNaN);
      std::vector<std::string> by(d.wells.size());  // what each well's point rests on
      std::size_t hp = hages.size();
      for (std::size_t q = 0; q < hages.size(); ++q)
        if (hages[q].first == H.name) hp = q;
      if (dated && hp < hages.size()) {
        row = cor.depth[hp];
        for (std::size_t k = 0; k < d.wells.size(); ++k)
          if (std::isfinite(row[k])) by[k] = "age " + fmt_shortest(H.age.si()) + " Ma, " + rule_text(d.correlation_rule);
      }
      for (std::size_t k = 0; k < d.wells.size(); ++k)
        if (!std::isfinite(row[k]))
          for (auto& [th, dp] : d.wells[k].tops)
            if (th == H.name) { row[k] = dp.si(); by[k] = "top"; }
      for (std::size_t k = 0; k < d.wells.size(); ++k)
        if (!std::isfinite(row[k])) {
          // the model horizon along the path (first crossing of the horizon's base unit)
          const auto& P = pw[k].path;
          const double td = P.total_depth_md(), step = std::max(1.0, td / 4000);
          int prev = m.unit_at(P.at_md(0));
          for (double md = step; md <= td; md += step) {
            const int u = m.unit_at(P.at_md(md));
            if (prev <= static_cast<int>(h) && u > static_cast<int>(h)) { row[k] = -P.at_md(md - step / 2).z; by[k] = "model"; break; }
            prev = u;
          }
        }
      // label: the wells on each basis ("age 23 Ma in W2, W3; top in W1")
      std::string basis;
      for (const char* kind : {"age", "top", "model"}) {
        std::string ws;
        for (std::size_t k = 0; k < d.wells.size(); ++k)
          if (by[k].rfind(kind, 0) == 0) ws += std::string(ws.empty() ? "" : ", ") + d.wells[k].name;
        if (ws.empty()) continue;
        const std::string what = std::string(kind) == "age" ? "age " + fmt_shortest(H.age.si()) + " Ma" : std::string(kind) == "top" ? "tops" : "model";
        basis += std::string(basis.empty() ? "" : "; ") + what + " in " + ws;
      }
      out.corr.horizons.push_back(H.name);
      out.corr.basis.push_back(basis.empty() ? "not placed" : basis);
      out.corr.depth.push_back(row);
      for (std::size_t k = 0; k < d.wells.size(); ++k)
        if (std::isfinite(row[k])) cc += H.name + "," + d.wells[k].name + "," + fmt_fixed(row[k], 2) + ",\"" + by[k] + "\"\n";
    }
    for (std::size_t k = 0; k < d.wells.size(); ++k) {
      fig::WellPanel pan;
      pan.well = d.wells[k].name;
      pan.index_label = "TVD (m)";
      pan.index = pw[k].tvd;
      if (!pw[k].logs.empty() && !pw[k].logs[0].curves.empty()) {
        const auto& c = pw[k].logs[0].curves[0];
        pan.tracks.push_back({c.mnemonic, c.unit, c.values});
      }
      for (auto& a : d.wells[k].ages) pan.ages.push_back({a.name, a.depth.si(), a.age.si()});
      const auto& P = pw[k].path;
      const double td = P.total_depth_md(), step = std::max(1.0, td / 2000);
      int cur = -2;
      for (double md = 0; md <= td + 1e-9; md += step) {
        const Vec3 p = P.at_md(md);
        const int u = m.unit_at(p);
        if (u != cur) {
          if (!pan.units.empty()) pan.units.back().base = -p.z;
          pan.units.push_back({u, -p.z, -p.z});
          cur = u;
        }
      }
      if (!pan.units.empty()) pan.units.back().base = -P.points().back().z;
      out.corr.wells.push_back(pan);
    }
    out.csv["well_correlation.csv"] = cc;
    if (dated) out.corr.note = "ages correlated " + rule_text(d.correlation_rule);
  }
  // logs upscaled to the grid
  if (grid && !d.grid.upscaled_logs.empty()) {
    for (auto& up : d.grid.upscaled_logs) {
      repr::CellProperty prop;
      prop.name = up.name.empty() ? up.curve : up.name;
      prop.evidence = "derived";
      prop.v.assign(grid->ncell(), kNaN);
      std::vector<int> count(grid->ncell(), 0);
      std::string wells;
      std::size_t filled = 0, few = 0;
      const auto meth = upscale_method(up.method);
      for (std::size_t k = 0; k < d.wells.size(); ++k)
        for (auto& s : pw[k].logs) {
          const auto* c = s.find(up.curve);
          if (!c || c == &s.index) continue;
          if (prop.unit.empty()) prop.unit = c->unit;
          repr::PathLog pl{d.wells[k].name, up.curve, {}, {}, c->values};
          for (std::size_t i = 0; i < s.index.values.size(); ++i) {
            const double v = s.index.values[i];
            pl.p.push_back(s.index_type == "MD" ? pw[k].path.at_md(v) : pw[k].path.at_elevation(pw[k].path.datum_z() - v));
            pl.md.push_back(s.index_type == "MD" ? v : pw[k].path.md_at_elevation(pw[k].path.datum_z() - v));
          }
          auto u = repr::upscale_log(*grid, pl, meth, up.min_valid.present ? up.min_valid.si() : 0.5);
          for (std::size_t q = 0; q < u.cells.size(); ++q) {
            const auto cell = u.cells[q];
            prop.v[cell] = count[cell] ? (prop.v[cell] * count[cell] + u.value[q]) / (count[cell] + 1) : u.value[q];
            ++count[cell];
          }
          few += u.missing.size();
          if (!u.cells.empty()) wells += std::string(wells.empty() ? "" : ", ") + d.wells[k].name;
        }
      for (auto n : count) filled += n ? 1 : 0;
      prop.source = "log " + up.curve + " upscaled by " + repr::upscale_method_name(meth) + " (weighted by measured-depth length, valid fraction >= " +
                    fmt_shortest(up.min_valid.present ? up.min_valid.si() : 0.5) + ") from " + (wells.empty() ? std::string("no well") : wells);
      for (double v : prop.v) prop.missing += std::isfinite(v) ? 0 : 1;
      prop.missing_reason = "not crossed by a well with the curve" + std::string(few ? " or too few valid samples" : "") + " (stage 4D fills cells by geostatistics)";
      rep += "\nupscaled " + prop.name + ": " + std::to_string(filled) + " cells filled; " + prop.source + "\n";
      out.upscaled.push_back(std::move(prop));
    }
  }
  if (!out.results.conflicts.empty()) {
    rep += "\nage conflicts:\n";
    for (auto& cf : out.results.conflicts) rep += "  " + cf.well + ": " + cf.what + "\n";
  }
  out.report = rep;
  return out;
}

}  // namespace stratonix
