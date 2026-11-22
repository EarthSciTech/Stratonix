// Pipeline shared by the command line and the C API.
#include "stratonix/package/pipeline.hpp"
#include "stratonix/package/welldata.hpp"
#include "stratonix/package/geostat_props.hpp"
#include "stratonix/io/maps.hpp"
#include "stratonix/io/seismic.hpp"
#include "stratonix/io/velocity.hpp"

#include <algorithm>
#include <sstream>

#include "stratonix/checks/checks.hpp"
#include "stratonix/checks/manifest.hpp"
#include "stratonix/engines/conditioning.hpp"
#include "stratonix/engines/geophysics.hpp"
#include "stratonix/engines/posterior.hpp"
#include "stratonix/engines/restoration3d.hpp"
#include "stratonix/engines/fractures.hpp"
#include "stratonix/engines/thermal.hpp"
#include "stratonix/repr/faultrock.hpp"
#include "stratonix/repr/geotime.hpp"
#include "stratonix/repr/topology.hpp"
#include "stratonix/engines/history.hpp"
#include "stratonix/engines/structure.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/figures/figures.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/io/export.hpp"
#include "stratonix/io/io.hpp"

namespace stratonix {

void register_all_plugins() {
  static bool done = false;
  if (done) return;
  done = true;
  register_engine_plugins();
  io::register_io_plugins();
  register_check_plugins();
  repr::register_grid_plugins();
}

void ensure_new_output_dir(const std::filesystem::path& out_dir) {
  std::error_code ec;
  if (std::filesystem::exists(out_dir, ec) && !std::filesystem::is_empty(out_dir, ec))
    throw Error(ErrorCode::Io, "output folder exists and is not empty: " + out_dir.generic_string() +
                                   " (issued outputs are never overwritten; choose a new folder or version)");
  std::filesystem::create_directories(out_dir);
}

BuildOutcome build_only(const desc::Description& d, const BuildOptions& opt) {
  register_all_plugins();
  if (opt.threads) set_thread_count(opt.threads);
  set_optimised_kernels(!opt.reference_kernels);
  BuildOutcome o;
  BuildInputs in;
  io::load_build_inputs(d, in);
  o.model = build_model(d, in, o.log, opt.builder);
  return o;
}

static std::string summary_text(const desc::Description& d, const Model& m, const BuildLog& log) {
  const auto& inf = m.info();
  std::ostringstream s;
  s << "Stratonix model summary\n";
  s << "model: " << inf.name << " " << inf.version << " (" << basis_label(inf.basis) << ")\n";
  s << "engine: " << inf.engine << "\n";
  s << "coordinate system: " << inf.crs.label() << "; vertical datum: " << inf.crs.vertical_datum << "\n";
  s << "domain: x " << fmt_fixed(inf.domain.x0, 1) << " to " << fmt_fixed(inf.domain.x0 + inf.domain.sx, 1)
    << ", y " << fmt_fixed(inf.domain.y0, 1) << " to " << fmt_fixed(inf.domain.y0 + inf.domain.sy, 1)
    << ", depth " << fmt_fixed(-inf.domain.z_top, 1) << " to " << fmt_fixed(-inf.domain.z_base, 1) << " m\n";
  s << "horizons:";
  for (auto& h : inf.horizons) s << " " << h.name;
  s << "\nunits:";
  for (auto& u : inf.units) s << " " << u.name;
  s << "\nfaults:";
  for (auto& f : inf.faults) s << " " << f.name;
  s << "\nparameters: " << d.all_params().size() << " (all with value, unit, evidence class and source)\n";
  for (auto& n : log.notes) s << "note: " << n << "\n";
  return s.str();
}

std::vector<std::pair<std::string, std::string>> key_numbers(const desc::Description& d, const Model& m, const CheckReport& rep,
                                                             CheckContext& ctx) {
  std::vector<std::pair<std::string, std::string>> k;
  const auto& inf = m.info();
  k.push_back({"Evidence basis", basis_label(inf.basis)});
  k.push_back({"Horizons / faults", std::to_string(inf.horizons.size()) + " / " + std::to_string(inf.faults.size())});
  for (auto& r : rep.results)
    if (r.id == "C9" && r.info) k.push_back({std::string(1, static_cast<char>(std::toupper(static_cast<unsigned char>(r.name[0])))) + r.name.substr(1) + " (" + r.subject + ")", fmt_fixed(r.value, r.unit == "km2" ? 2 : 1) + " " + r.unit});
  // per fault; a network of many faults (stage 4B fault systems) lists its six largest throws and
  // sums up the rest in one line, so the panel stays on the page
  struct Row { std::string name, value; double best; };
  std::vector<Row> rows;
  for (std::size_t f = 0; f < inf.faults.size(); ++f) {
    double best = 0;
    std::string hn;
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      double t = std::fabs(ctx.throw_row(static_cast<int>(f), static_cast<int>(h)).throw_max);
      if (t > best) { best = t; hn = inf.horizons[h].name; }
    }
    rows.push_back({inf.faults[f].name, hn.empty() ? "0 m" : fmt_fixed(best, 1) + " m at " + hn, best});
  }
  constexpr std::size_t kListed = 6;
  if (rows.size() > kListed + 2) {
    std::stable_sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.best > b.best; });
    double rest = 0;
    for (std::size_t i = kListed; i < rows.size(); ++i) rest = std::max(rest, rows[i].best);
    const std::size_t n_rest = rows.size() - kListed;
    rows.resize(kListed);
    rows.push_back({"other " + std::to_string(n_rest) + " faults", "at most " + fmt_fixed(rest, 1) + " m", rest});
  }
  for (auto& r : rows) k.push_back({"Max throw " + r.name, r.value});
  k.push_back({"Checks", std::to_string(rep.results.size()) + " run, " + std::to_string(rep.failures()) + " failed"});
  (void)d;
  return k;
}

std::string data_misfit_csv(const Model& m) {
  auto* mis = static_cast<const std::vector<DataMisfit>*>(m.capability("data_misfit"));
  if (!mis) return {};
  // stage 3A columns (STX-TN-11 section 5) only when some datum has an uncertainty or a leave-one-out error
  bool extra = std::any_of(mis->begin(), mis->end(), [](const DataMisfit& r) { return !std::isnan(r.sigma) || !std::isnan(r.loo); });
  auto num = [](double v, int dec) { return std::isnan(v) ? std::string() : fmt_fixed(v, dec); };
  std::string s = std::string("input,evidence,horizon,kind,x,y,z,misfit,misfit_unit") + (extra ? ",sigma,normalised_misfit,leave_one_out_m" : "") + "\n";
  for (auto& r : *mis) {
    s += r.input + "," + r.evidence + "," + r.label + "," + r.kind + "," + fmt_fixed(r.p.x, 2) + "," + fmt_fixed(r.p.y, 2) + "," + fmt_fixed(r.p.z, 2) +
         "," + fmt_fixed(r.misfit, 4) + "," + (r.kind == "orientation" ? "deg" : "m");
    if (extra) s += "," + num(r.sigma, 4) + "," + num(r.normalised, 4) + "," + num(r.loo, 4);
    s += "\n";
  }
  return s;
}

std::string summary_csv(const std::vector<std::pair<std::string, std::string>>& k) {
  std::string s = "key,value\n";
  for (auto& [a, b] : k) s += "\"" + a + "\",\"" + b + "\"\n";
  return s;
}

// ---- stage 4B (STX-TN-15 section 4): restoration, backstripping and eroded-section helpers -------
namespace {
BackstripSpec backstrip_spec(const desc::Description& d) {
  const auto& B = d.backstripping;
  BackstripSpec s;
  s.basement = d.horizon_index(B.basement);
  s.rho_m = B.mantle_density.si();
  s.rho_w = B.water_density.si();
  s.rho_g = B.grain_density.si();
  s.te = B.te.si();
  if (B.young.present) s.young = B.young.si();
  if (B.poisson.present) s.poisson = B.poisson.si();
  s.water_depth.assign(d.horizons.size(), kNaN);
  s.sea_level.assign(d.horizons.size(), kNaN);
  bool present_given = false;
  for (auto& x : B.steps) {
    if (x.horizon == "present") {
      s.water_depth_present = x.water_depth.si();
      s.sea_level_present = x.sea_level.present ? x.sea_level.si() : 0.0;
      present_given = true;
      continue;
    }
    const int h = d.horizon_index(x.horizon);
    if (h < 0) continue;
    s.water_depth[static_cast<std::size_t>(h)] = x.water_depth.si();
    s.sea_level[static_cast<std::size_t>(h)] = x.sea_level.present ? x.sea_level.si() : 0.0;
  }
  if (!present_given && !d.horizons.empty()) {  // the present step takes the top horizon's water depth
    s.water_depth_present = std::isfinite(s.water_depth[0]) ? s.water_depth[0] : 0.0;
    s.sea_level_present = std::isfinite(s.sea_level[0]) ? s.sea_level[0] : 0.0;
  }
  return s;
}
/// Elevation of horizon h on a vertical at (x, y): where the unit changes from h to h + 1 downward.
double horizon_elevation(const Model& m, int h, double x, double y) {
  const auto& dom = m.info().domain;
  const int n = 600;
  std::vector<Vec3> pts(n);
  for (int k = 0; k < n; ++k) pts[static_cast<std::size_t>(k)] = {x, y, dom.z_top - (dom.z_top - dom.z_base) * (k + 0.5) / n};
  EvalResult er;
  m.evaluate(pts, kEvalUnit, er);
  for (int k = 1; k < n; ++k)
    if (er.unit[static_cast<std::size_t>(k - 1)] <= h && er.unit[static_cast<std::size_t>(k)] > h) {
      double a = pts[static_cast<std::size_t>(k - 1)].z, b = pts[static_cast<std::size_t>(k)].z;
      for (int it = 0; it < 30; ++it) {
        const double zm = 0.5 * (a + b);
        (m.unit_at({x, y, zm}) <= h ? a : b) = zm;
      }
      return 0.5 * (a + b);
    }
  return kNaN;
}
}  // namespace

fig::FigureContext figure_context(const desc::Description& d, const Model& m, const std::vector<std::pair<std::string, std::string>>& keys) {
  fig::FigureContext c;
  c.model = &m;
  c.model_name = d.model.name;
  c.model_version = d.model.version;
  c.title = d.model.title;
  c.wells = make_well_paths(d);
  c.key_numbers = keys;
  c.vertical_exaggeration = d.vertical_exaggeration.present ? d.vertical_exaggeration.si() : 1.0;
  c.grid_spacing = d.grid_spacing.si();
  c.closure_horizon = d.closure_horizon.empty() ? -1 : d.horizon_index(d.closure_horizon);
  c.faults_seal = d.closure_faults_seal.present && d.closure_faults_seal.flag;
  for (auto& u : m.info().units) {
    std::string role, color, label;
    for (auto& du : d.units)
      if (du.name == u.name) { role = du.role; color = du.color; label = du.label; }
    for (auto& b : d.bodies)  // stage 4B: bodies carry their own colour
      if (b.name == u.name) { role = "Bodies"; color = b.color; }
    for (auto& g : d.geobodies)  // stage 4B: geobody sets too, listed as facies
      if (g.name == u.name) { role = "Geobodies"; color = g.color; }
    c.unit_roles.push_back(role);
    c.unit_colors.push_back(color);
    c.unit_labels.push_back(label);
  }
  for (auto& n : d.isopachs)
    for (std::size_t u = 0; u < m.info().units.size(); ++u)
      if (m.info().units[u].name == n) c.isopach_units.push_back(static_cast<int>(u));
  // Burial history at every well head (stage 2D), when the description has compaction.
  if (auto* st = static_cast<const Structure*>(m.capability("structure")); st && st->compaction_mode() != 0) {
    BackstripMaps maps;  // stage 4B: flexural maps, computed once for all wells when needed
    for (auto& w : c.wells) {
      auto hist = st->burial_history(w.x(), w.y());
      if (hist.empty()) continue;
      fig::FigureContext::Burial b;
      b.location = w.name();
      bool all_ages = true;
      for (auto& stp : hist) {
        b.step_labels.push_back(stp.label);
        b.depth.push_back(stp.depth);
        all_ages = all_ages && !std::isnan(stp.age);
      }
      if (all_ages)
        for (auto& stp : hist) b.ages.push_back(stp.age);
      for (std::size_t h = 0; h < hist.back().depth.size(); ++h)
        if (!std::isnan(hist.back().depth[h])) b.horizons.push_back(static_cast<int>(h));
      if (d.backstripping.present) {
        // stage 4B: tectonic subsidence on the same steps (steps without the basement stay blank)
        const BackstripSpec bs = backstrip_spec(d);
        auto ts = backstrip_column(*st, bs, w.x(), w.y());
        b.tectonic_airy.assign(hist.size(), kNaN);
        b.tectonic_flexural.assign(hist.size(), kNaN);
        for (auto& t : ts)
          for (std::size_t k = 0; k < hist.size(); ++k)
            if (hist[k].label == t.label) b.tectonic_airy[k] = t.airy;
        if (bs.te > 0) {
          if (maps.labels.empty()) maps = backstrip_maps(*st, bs, m.info().domain, 32);
          for (std::size_t q = 0; q < maps.labels.size(); ++q)
            for (std::size_t k = 0; k < hist.size(); ++k)
              if (hist[k].label == maps.labels[q]) b.tectonic_flexural[k] = backstrip_map_value(maps, maps.flexural[q], m.info().domain, w.x(), w.y());
        } else {
          b.tectonic_flexural.clear();
        }
      }
      c.burial.push_back(std::move(b));
    }
  }
  c.fault_rock = d.fault_rock.present;
  if (d.fault_rock.present || d.stress.present) {
    c.rock = fault_rock_input(d, m);
    c.stress = stress_input(d, m);
  }
  for (auto& hn : d.curvature_horizons)
    for (std::size_t h = 0; h < m.info().horizons.size(); ++h)
      if (m.info().horizons[h].name == hn) c.curvature_horizons.push_back(static_cast<int>(h));
  if (d.fracture_threshold.present) c.fracture_threshold = d.fracture_threshold.si();
  return c;
}

fig::SheetSpec sheet_spec(const desc::Description& d, const Model& m, const desc::SheetDesc& sh) {
  fig::SheetSpec sp;
  const auto& inf = m.info();
  sp.name = sh.name;
  sp.title = sh.title;
  sp.along_label = sh.along_label.empty() ? "Distance along section (m)" : sh.along_label;
  sp.depth_label = sh.depth_label.empty() ? "Depth below datum (m)" : sh.depth_label;
  sp.end_left = sh.end_left;
  sp.end_right = sh.end_right;
  for (auto& s : d.sections)
    if (s.name == sh.section)
      for (std::size_t r = 0; r < s.polyline.rows(); ++r) sp.line.push_back({s.polyline.si(2 * r), s.polyline.si(2 * r + 1), 0});
  auto unit_of = [&](const std::string& n) {
    for (std::size_t u = 0; u < inf.units.size(); ++u)
      if (inf.units[u].name == n) return static_cast<int>(u);
    return -1;
  };
  auto panel = [&](const desc::SheetPanelDesc& pd, fig::SheetPanel& p) {
    p.present = pd.present;
    p.title = pd.title;
    p.ve = pd.vertical_exaggeration.present ? pd.vertical_exaggeration.si() : 1.0;
    if (pd.window.present && pd.window.values.size() == 4) {
      p.s0 = pd.window.si(0); p.s1 = pd.window.si(1); p.d0 = pd.window.si(2); p.d1 = pd.window.si(3);
    }
    p.label_at = pd.label_at.present ? pd.label_at.si() : -1;
    for (auto& u : pd.label_units)
      if (int k = unit_of(u); k >= 0) p.label_units.push_back(k);
    for (auto& a : pd.annotations) {
      fig::SheetAnnotation x;
      x.type = a.type;
      x.text = a.text;
      x.horizon = a.horizon.empty() ? -1 : d.horizon_index(a.horizon);
      x.fault = a.fault.empty() ? -1 : d.fault_index(a.fault);
      x.footwall = a.side != "hanging_wall";
      x.well = a.well;
      x.completion = a.completion;
      x.at = a.at.present ? a.at.si() : 0;
      if (a.text_at.present && a.text_at.values.size() == 2) { x.text_s = a.text_at.si(0); x.text_depth = a.text_at.si(1); }
      p.annotations.push_back(std::move(x));
    }
  };
  panel(sh.main, sp.main);
  panel(sh.detail, sp.detail);
  sp.well_legend = sh.well_legend;
  sp.completion_legend = sh.completion_legend;
  sp.fault_legend = sh.fault_legend;
  sp.notes = sh.notes;
  for (auto& k : sh.keys) {
    fig::SheetKey x;
    x.kind = k.kind;
    x.label = k.label;
    x.horizon = k.horizon.empty() ? -1 : d.horizon_index(k.horizon);
    x.unit = k.unit.empty() ? -1 : unit_of(k.unit);
    x.fault = k.fault.empty() ? -1 : d.fault_index(k.fault);
    for (auto& f : k.faults) x.faults.push_back(d.fault_index(f));
    x.well = k.well;
    x.end = k.end;
    x.fraction = k.fraction.present ? k.fraction.si() : 0.5;
    x.decimals = k.decimals >= 0                                                                    ? k.decimals
                 : (k.kind == "fault_tips" || k.kind == "well_fault_distance" || k.kind == "fault_cutoff_position" ||
                    k.kind == "fault_gap" || k.kind == "fault_dip" || k.kind == "completion")          ? 0
                 : (k.kind == "fault_overlap" || k.kind == "throw_ratio")                              ? 2
                                                                                                      : 1;
    x.basis = k.basis;
    x.on_sheet = k.on_sheet;
    x.value = k.value;
    x.completion = k.completion;
    x.side = k.side;
    x.coordinate = k.coordinate;
    x.from = k.from.empty() ? -1 : d.horizon_index(k.from);
    x.to = k.to.empty() ? -1 : d.horizon_index(k.to);
    sp.keys.push_back(std::move(x));
  }
  sp.summary_name = d.model.name + "_v" + d.model.version + "_" + sh.name + "_summary.csv";
  sp.key_title = sh.key_title;
  for (auto& w : d.wells)
    for (auto& cp : w.completions) sp.completions.push_back({w.name, cp.name, cp.top.si(), cp.base.si()});
  if (auto* st = static_cast<const Structure*>(m.capability("structure")))
    sp.unfaulted = [st](int h, double x, double y) { return st->z(static_cast<std::size_t>(h), x, y); };
  return sp;
}

void write_sheets(const desc::Description& d, const Model& m, const fig::FigureContext& c, const std::filesystem::path& folder) {
  for (auto& sh : d.sheets) {
    auto sp = sheet_spec(d, m, sh);
    auto res = fig::section_sheet(c, sp);
    auto q = [](const std::string& t) {
      std::string o = "\"";
      for (char ch : t) o += ch == '"' ? std::string("\"\"") : std::string(1, ch);
      return o + "\"";
    };
    std::string csv = "quantity,value,basis\n";
    for (auto& k : res.keys) csv += q(k.label) + "," + q(k.text) + "," + q(k.basis) + "\n";
    std::filesystem::create_directories(folder);
    write_file(folder / sp.summary_name, csv);
    std::string notes = std::string("Evidence basis: ") + basis_label(m.info().basis) + ". Key numbers computed by Stratonix from the model (" +
                        sp.summary_name + ").";
    fig::write_figure(res.figure, folder, fig::figure_basename(c, "sheet", sh.name), c.model_name + " - section sheet " + sh.name,
                      "Section sheet " + sh.name + ": main panel, detail panel, legend by role and key numbers.", notes);
  }
}

/// Horizon ages (Ma) in model order when the description gives one to every horizon, else empty.
std::vector<double> horizon_ages(const desc::Description& d, const Model& m) {
  std::vector<double> a;
  for (auto& h : m.info().horizons) {
    int k = d.horizon_index(h.name);
    if (k < 0 || !d.horizons[static_cast<std::size_t>(k)].age.present) return {};
    a.push_back(d.horizons[static_cast<std::size_t>(k)].age.si());
  }
  for (std::size_t k = 1; k < a.size(); ++k)
    if (!(a[k] > a[k - 1])) return {};
  return a;
}

void write_figures(const desc::Description& d, const Model& m, const std::vector<std::pair<std::string, std::string>>& keys,
                   const std::filesystem::path& folder, repr::FaultedSurfaceCache* surfaces, const BurialFields* fields) {
  auto c = figure_context(d, m, keys);
  c.surfaces = surfaces;
  if (fields)
    for (auto& b : c.burial) {
      auto it = fields->find(b.location);
      if (it == fields->end()) continue;
      // the field's steps are matched to the burial steps by label
      for (auto& [labels, fd] : it->second) {
        fig::FigureContext::Burial::Field g{fd.kind, fd.label, fd.unit, {}};
        for (auto& lb : b.step_labels) {
          auto p = std::find(labels.begin(), labels.end(), lb);
          g.v.push_back(p == labels.end() ? std::vector<double>(b.depth.front().size(), kNaN) : fd.v[static_cast<std::size_t>(p - labels.begin())]);
        }
        b.fields.push_back(std::move(g));
      }
    }
  std::vector<std::pair<std::string, Polyline>> secs;
  std::vector<std::vector<double>> det;
  for (auto& s : d.sections) {
    Polyline pl;
    for (std::size_t r = 0; r < s.polyline.rows(); ++r) pl.push_back({s.polyline.si(2 * r), s.polyline.si(2 * r + 1), 0});
    secs.push_back({s.name, pl});
    det.push_back(s.detail_window.present ? s.detail_window.si_list() : std::vector<double>{});
  }
  if (secs.empty()) {  // default: a section through the domain centre along x
    const Domain& dm = m.info().domain;
    secs.push_back({"X_centre", {{dm.x0, dm.y0 + dm.sy / 2, 0}, {dm.x0 + dm.sx, dm.y0 + dm.sy / 2, 0}}});
    det.push_back({});
  }
  fig::write_figure_set(c, secs, det, folder);
  write_sheets(d, m, c, folder);
  // stage 4A (STX-TN-14 section 5): Wheeler diagrams along the sections, when every horizon has an age
  if (auto ages = horizon_ages(d, m); !ages.empty()) {
    GeologicalTime gt(m, ages);
    for (auto& [name, pl] : secs) {
      auto w = wheeler_diagram(gt, name, pl);
      fig::write_figure(fig::wheeler_figure(c, w), folder, fig::figure_basename(c, "wheeler", name), c.model_name + " - Wheeler diagram " + name,
                        "Distance along the section against geological time, coloured by unit; hiatuses blank.",
                        "Geological time from the model's restoration and the stated horizon ages.");
    }
  }
}

std::vector<std::string> resolve_export_ids(const desc::Description& d, const std::vector<std::string>& requested) {
  std::vector<std::string> want = requested.empty() ? d.exports : requested;
  std::vector<std::string> ids;
  auto& reg = registry::Registry::global();
  for (auto& w : want) {
    if (w == "petrel") {
      for (auto& p : io::petrel_set_ids())
        if (std::find(ids.begin(), ids.end(), p) == ids.end()) ids.push_back(p);
      continue;
    }
    if (!reg.find(registry::Kind::Exporter, w)) {
      log(LogLevel::Warn, "export format '" + w + "' is not available in this version; skipped");
      continue;
    }
    if (std::find(ids.begin(), ids.end(), w) == ids.end()) ids.push_back(w);
  }
  return ids;
}

std::vector<io::ExportedFile> write_exports(const desc::Description& d, const Model& m, const std::vector<std::string>& ids,
                                            const std::filesystem::path& folder, double grid_spacing, const GridProducts* gp = nullptr,
                                            repr::FaultedSurfaceCache* surfaces = nullptr) {
  double sp = grid_spacing > 0 ? grid_spacing : d.grid_spacing.si();
  double sticks = d.fault_stick_spacing.present ? d.fault_stick_spacing.si() : 250.0;  // Stratonix default when not given
  io::ExportContext c(m, make_well_paths(d), GridSpec::covering(m.info().domain, sp), std::max(10.0, sp / 2), sticks);
  c.share_surfaces(surfaces);
  c.vertical_exaggeration = d.vertical_exaggeration.present ? d.vertical_exaggeration.si() : 1.0;
  c.fault_rock = d.fault_rock.present;
  if (d.fault_rock.present || d.stress.present) {
    c.rock = fault_rock_input(d, m);
    c.stress = stress_input(d, m);
  }
  for (auto& sr : d.series) c.series.push_back({sr.name, sr.base, sr.horizons});
  if (d.implicit.present) c.isochore_horizons = d.implicit.isochore_horizons;
  for (auto& hn : d.curvature_horizons) {
    const int h = d.horizon_index(hn);
    if (h >= 0) c.curvature_horizons.push_back(h);
  }
  if (d.fracture_threshold.present) c.fracture_threshold = d.fracture_threshold.si();
  if (d.voxel_spacing.present) {
    auto v = d.voxel_spacing.si_list();
    if (v.size() == 3) { c.voxel_dx = v[0]; c.voxel_dy = v[1]; c.voxel_dz = v[2]; }
  }
  c.sample_spacing = std::max(250.0, 5 * sp);  // Stratonix default: sparse enough for implicit modellers
  std::vector<repr::FaultFace> faces;
  if (gp) {
    c.cells.grid = &gp->grid;
    c.cells.props = &gp->props;
    if (gp->has_trans) {
      c.cells.trans = &gp->trans;
      c.cells.faces = &gp->trans.faces;
    } else {
      faces = repr::grid_fault_faces(m, gp->grid);
      c.cells.faces = &faces;
    }
    if (gp->has_tets) c.cells.tets = &gp->tets;
  }
  std::string note = "Every parameter of the description carries a value, unit, evidence class and source (see manifest.json).";
  return io::export_set(c, ids, folder, note);
}

/// Summary rows for the events of a timeline (stage 2E): each event's displacement.
std::vector<std::pair<std::string, std::string>> event_rows(const desc::Description& d) {
  std::vector<std::pair<std::string, std::string>> k;
  for (auto& e : d.events) {
    const std::string when = e.from.present ? fmt_shortest(e.from.si()) + " to " + fmt_shortest(e.to.si()) + " Ma" : fmt_shortest(e.age.si()) + " Ma";
    std::string what;
    if (e.type == "deposition") what = "deposition of " + e.ref;
    else if (e.type == "fault") {
      const auto& fd = d.faults[static_cast<std::size_t>(d.fault_index(e.ref))];
      const double t = e.throw_.present ? e.throw_.si() : fd.throw_max.si();
      const std::string sense = !e.sense.empty() ? e.sense : (fd.sense.empty() ? "normal" : fd.sense);
      what = "fault " + e.ref + ": throw " + fmt_fixed(t, 1) + " m " + sense;
    } else if (e.type == "fold") {
      double r = 0;
      for (auto& p : d.primitives)
        if (p.name == e.ref)
          for (auto& [hn, q] : p.relief) r = std::max(r, std::fabs(q.si()));
      what = "fold " + e.ref + ": relief " + fmt_fixed(r, 1) + " m";
    } else {
      for (auto& b : d.salt.bodies)
        if (b.name == e.ref) what = "salt body " + e.ref + ": height " + fmt_fixed(b.height.si(), 1) + " m";
    }
    k.push_back({"Event " + e.name + " (" + when + ")", what});
  }
  return k;
}

/// Summary rows for fault rock and stress (stage 2F): the inputs, labelled, then per fault the
/// range of each attribute over its surface.
std::vector<std::pair<std::string, std::string>> fault_rock_rows(const desc::Description& d, const Model& m) {
  std::vector<std::pair<std::string, std::string>> k;
  if (!d.fault_rock.present && !d.stress.present) return k;
  auto lab = [](const desc::Param& p) { return p.present ? std::string(p.text) + " (" + desc::evidence_name(p.evidence) + (p.source.empty() ? std::string() : ", " + p.source) + ")" : std::string("not given"); };
  if (d.fault_rock.present) {
    std::string clay;
    for (auto& u : d.units)
      if (u.clay_fraction.present) clay += (clay.empty() ? "" : "; ") + u.name + " " + u.clay_fraction.text;
    k.push_back({"Fault rock: clay fractions", clay.empty() ? "not given" : clay});
    k.push_back({"Fault rock: shale threshold", lab(d.fault_rock.shale_threshold)});
    k.push_back({"Fault rock: thickness / displacement", lab(d.fault_rock.thickness_ratio)});
    k.push_back({"Fault rock: seal calibration", d.fault_rock.seal_type.empty() ? "none (no seal capacity reported)"
                                                                                  : d.fault_rock.seal_type + ", a " + lab(d.fault_rock.seal_a) + ", c " + lab(d.fault_rock.seal_c)});
  }
  if (d.stress.present) {
    const auto& st = d.stress;
    if (st.regime())
      k.push_back({"Stress: regime", "Sv " + lab(st.sv_gradient) + "; SHmax/Sv " + lab(st.shmax_ratio) + "; Shmin/Sv " + lab(st.shmin_ratio) +
                                         "; SHmax azimuth " + lab(st.shmax_azimuth)});
    else
      k.push_back({"Stress: tensor gradient (xx, yy, zz, xy, xz, yz; east, north, up)", lab(st.tensor_gradient)});
    k.push_back({"Stress: pore pressure gradient", lab(st.pore_pressure_gradient)});
    k.push_back({"Stress: friction", lab(st.friction)});
    k.push_back({"Stress: cohesion", st.cohesion.present ? lab(st.cohesion) : "0 (not given)"});
  }
  const auto in = fault_rock_input(d, m);
  const auto sin = stress_input(d, m);
  const double spacing = std::max(m.info().domain.spacing * 2, 50.0);
  for (std::size_t f = 0; f < m.info().faults.size(); ++f) {
    auto pts = repr::fault_rock_surface(m, repr::fault_surface(m, static_cast<int>(f), spacing), in, sin);
    auto range = [&](auto get, double scale, int dec, const std::string& unit) {
      double lo = 1e300, hi = -1e300;
      for (auto& p : pts)
        if (p.valid && std::isfinite(get(p))) { lo = std::min(lo, get(p) * scale); hi = std::max(hi, get(p) * scale); }
      return hi >= lo ? fmt_fixed(lo, dec) + " to " + fmt_fixed(hi, dec) + (unit.empty() ? "" : " " + unit) : std::string("missing");
    };
    const std::string fn = "Fault " + m.info().faults[f].name + ": ";
    if (d.fault_rock.present) {
      std::size_t n = 0, self = 0, sand = 0;
      for (auto& p : pts) {
        if (!p.valid || p.fw_unit < 0 || p.hw_unit < 0) continue;
        ++n;
        if (p.fw_unit == p.hw_unit) ++self;
        auto clay = [&](int u) { return static_cast<std::size_t>(u) < in.unit_clay.size() ? in.unit_clay[static_cast<std::size_t>(u)] : kNaN; };
        const double cf = clay(p.fw_unit), ch = clay(p.hw_unit);
        if (!std::isnan(in.shale_threshold) && cf < in.shale_threshold && ch < in.shale_threshold) ++sand;
      }
      auto pct = [&](std::size_t x) { return n ? fmt_fixed(100.0 * static_cast<double>(x) / static_cast<double>(n), 1) + " %" : std::string("-"); };
      k.push_back({fn + "juxtaposition", "self-juxtaposed " + pct(self) + "; sand on sand (both below the shale threshold) " +
                                             (std::isnan(in.shale_threshold) ? std::string("missing (no shale threshold)") : pct(sand)) +
                                             "; of the fault surface area"});
      k.push_back({fn + "SGR", range([](const repr::FaultRockPoint& p) { return p.sgr; }, 1.0, 3, "")});
      k.push_back({fn + "CSP", range([](const repr::FaultRockPoint& p) { return p.csp; }, 1.0, 1, "m")});
      k.push_back({fn + "fault-zone thickness", range([](const repr::FaultRockPoint& p) { return p.thickness; }, 1.0, 2, "m")});
      if (in.seal) k.push_back({fn + "seal capacity", range([](const repr::FaultRockPoint& p) { return p.seal_bar; }, 1.0, 2, "bar")});
    }
    if (sin.present) {
      k.push_back({fn + "slip tendency", range([](const repr::FaultRockPoint& p) { return p.ts; }, 1.0, 3, "")});
      k.push_back({fn + "dilation tendency", range([](const repr::FaultRockPoint& p) { return p.td; }, 1.0, 3, "")});
      k.push_back({fn + "pressure to reactivate", range([](const repr::FaultRockPoint& p) { return p.dp; }, 1e-6, 2, "MPa")});
    }
  }
  return k;
}

void write_restoration(const desc::Description& d, const Model& m, const std::filesystem::path& folder) {
  if (d.events.empty()) return;
  std::filesystem::create_directories(folder);
  auto ages = event_ages(d);
  ages.push_back(0.0);
  Polyline line;
  std::string sname = "X_centre";
  if (!d.sections.empty()) {
    sname = d.sections.front().name;
    for (std::size_t r = 0; r < d.sections.front().polyline.rows(); ++r)
      line.push_back({d.sections.front().polyline.si(2 * r), d.sections.front().polyline.si(2 * r + 1), 0});
  } else {
    const Domain& dm = m.info().domain;
    line = {{dm.x0, dm.y0 + dm.sy / 2, 0}, {dm.x0 + dm.sx, dm.y0 + dm.sy / 2, 0}};
  }
  std::string csv = "age_ma,events,horizons,faults,checks_run,checks_failed\n";
  std::vector<fig::Slide> slides;
  for (std::size_t k = 0; k < ages.size(); ++k) {
    const double a = ages[k];
    auto sd = state_description(d, m, a);
    if (sd.horizons.empty()) continue;
    BuildLog lg;
    auto sm = build_model(sd, BuildInputs{}, lg, "parametric");
    CheckContext ctx(*sm, sd, make_well_paths(sd), sd.grid_spacing.si());
    ctx.timeline_checks = false;
    auto rep = run_checks(ctx);
    std::string names;
    for (auto& e : d.events)
      if (std::fabs(e.end() - a) < 1e-9 || std::fabs(e.start() - a) < 1e-9) names += (names.empty() ? "" : " ") + e.name;
    csv += fmt_shortest(a) + ",\"" + names + "\"," + std::to_string(sd.horizons.size()) + "," + std::to_string(sd.faults.size()) + "," +
           std::to_string(rep.results.size()) + "," + std::to_string(rep.failures()) + "\n";
    std::vector<std::pair<std::string, std::string>> keys{{"State", fmt_shortest(a) + " Ma"}, {"Events", names.empty() ? "-" : names},
                                                         {"Checks", std::to_string(rep.results.size()) + " run, " + std::to_string(rep.failures()) + " failed"}};
    auto c = figure_context(sd, *sm, keys);
    std::string tag = (k + 1 < 10 ? "0" : "") + std::to_string(k + 1) + "_" + fmt_shortest(a) + "Ma";
    std::replace(tag.begin(), tag.end(), '.', 'p');
    auto s = repr::make_section(*sm, sname, line, std::max(5.0, c.grid_spacing / 4), 2.0);
    auto sf = fig::section_figure(c, s, nullptr);
    const std::string title = c.model_name + " at " + fmt_shortest(a) + " Ma - section " + sname;
    auto w = fig::write_figure(sf, folder, fig::figure_basename(c, "evolution", tag + "_section"), title,
                               "State at " + fmt_shortest(a) + " Ma" + (names.empty() ? "" : " (after " + names + ")") + ": section " + sname + ".",
                               "Restored state built from the description; horizons deposited later and events after it removed.");
    if (fig::raster_available()) {
      std::string png = read_file(folder / (w.basename + ".png"));
      auto u32 = [&](std::size_t o) { return (static_cast<unsigned char>(png[o]) << 24) | (static_cast<unsigned char>(png[o + 1]) << 16) |
                                             (static_cast<unsigned char>(png[o + 2]) << 8) | static_cast<unsigned char>(png[o + 3]); };
      slides.push_back({title, w.caption, "State " + std::to_string(k + 1) + " of " + std::to_string(ages.size()), png,
                        static_cast<int>(u32(16)), static_cast<int>(u32(20))});
    }
    const Domain& dm = sm->info().domain;
    GridSpec g = GridSpec::covering(dm, std::max(c.grid_spacing, std::max(dm.sx, dm.sy) / 300));
    auto fs = repr::horizon_grid_faulted(*sm, 0, g);
    std::vector<std::vector<repr::CutoffLine>> cut;
    for (std::size_t ff = 0; ff < sm->info().faults.size(); ++ff)
      cut.push_back(repr::cutoff_lines(*sm, repr::fault_surface(*sm, static_cast<int>(ff), std::max(10.0, g.dx / 2)), 0));
    auto mf = fig::map_figure(c, 0, fs, cut, nullptr);
    fig::write_figure(mf, folder, fig::figure_basename(c, "evolution", tag + "_map"),
                      c.model_name + " at " + fmt_shortest(a) + " Ma - " + sd.horizons.front().name + " depth map",
                      "State at " + fmt_shortest(a) + " Ma: depth map of the youngest horizon, " + sd.horizons.front().name + ".",
                      "Restored state built from the description.");
  }
  write_file(folder / "states_checks.csv", csv);
  if (!slides.empty()) write_file(folder / (m.info().name + "_v" + m.info().version + "_evolution-sequence.pptx"),
                                  fig::make_pptx(slides, m.info().name + " evolution sequence"));
  // Balance and the flexural-slip unfolded section of the present model.
  std::string b = "section,subject,method,present,restored,mismatch,tested\n";
  for (auto& r : balance_rows(m, d))
    b += r.section + "," + r.subject + ",\"" + r.method + "\"," + fmt_fixed(r.present, 3) + "," + fmt_fixed(r.restored, 3) + "," +
         fmt_fixed(r.mismatch, 6) + "," + (r.tested ? "yes" : "no") + "\n";
  write_file(folder / "balance.csv", b);
  auto c = figure_context(d, m, {{"Method", "flexural slip"}});
  auto s = repr::make_section(m, sname, line, std::max(5.0, c.grid_spacing / 4), 2.0);
  auto uf = fig::unfolded_figure(c, sname, repr::unfold_flexural(m, s));
  fig::write_figure(uf, folder, fig::figure_basename(c, "unfolded", sname), c.model_name + " - section " + sname + " unfolded",
                    "Section " + sname + " unfolded by flexural slip.", "Bed lengths conserved by construction; compare them for line-length balance.");
}

/// Stage 4A posterior samples in the parameters' units, for the figures (STX-TN-14 section 4).
fig::PosteriorSamples posterior_samples(const desc::Description& d, const PosteriorResult& pr) {
  fig::PosteriorSamples s;
  std::vector<double> to;
  for (auto& p : pr.paths) {
    const desc::Param* q = d.param_by_name(p);
    s.names.push_back(p);
    s.units.push_back(q ? q->unit : "");
    to.push_back(q ? q->to_si : 1.0);
  }
  s.walkers = pr.walkers;
  s.steps = pr.steps;
  s.burn_in = pr.burn_in;
  s.thin = pr.thin;
  s.seed = pr.seed;
  s.acceptance = pr.acceptance;
  s.converged = pr.converged;
  s.rhat = pr.rhat;
  s.ess = pr.ess;
  for (std::size_t k = 0; k < to.size(); ++k) {
    s.mean.push_back(pr.mean[k] / to[k]);
    s.q025.push_back(pr.q025[k] / to[k]);
    s.q975.push_back(pr.q975[k] / to[k]);
    s.lsq.push_back(k < pr.lsq.size() ? pr.lsq[k] / to[k] : kNaN);
  }
  s.chain.reserve(pr.chain.size());
  for (auto& v : pr.chain) {
    std::vector<double> u(v.size());
    for (std::size_t k = 0; k < v.size(); ++k) u[k] = v[k] / to[k];
    s.chain.push_back(std::move(u));
  }
  return s;
}

/// posterior_samples.csv (every thin-th step of every walker) and posterior_summary.csv.
void write_posterior_tables(const desc::Description& d, const PosteriorResult& pr, const std::filesystem::path& dir) {
  auto s = posterior_samples(d, pr);
  std::string t = "walker,step";
  for (std::size_t k = 0; k < s.names.size(); ++k) t += "," + s.names[k] + " (" + s.units[k] + ")";
  t += ",log_posterior\n";
  for (std::size_t w = 0; w < static_cast<std::size_t>(s.walkers); ++w)
    for (std::size_t st = static_cast<std::size_t>(s.thin) - 1; st < static_cast<std::size_t>(s.steps); st += static_cast<std::size_t>(s.thin)) {
      const std::size_t i = w * static_cast<std::size_t>(s.steps) + st;
      t += std::to_string(w + 1) + "," + std::to_string(st + 1);
      for (double v : s.chain[i]) t += "," + fmt_fixed(v, 6);
      t += "," + fmt_fixed(pr.log_posterior[i], 6) + "\n";
    }
  write_file(dir / "posterior_samples.csv", t);
  std::string u = "path,unit,prior,least_squares,laplace_sigma,mean,sd,q025,q975,rhat,ess,walkers,steps,burn_in,thin,seed,acceptance,converged\n";
  for (std::size_t k = 0; k < s.names.size(); ++k) {
    const double to = d.param_by_name(s.names[k]) ? d.param_by_name(s.names[k])->to_si : 1.0;
    const auto& P = pr.priors.at(k);
    std::string prior = P.distribution + "(" + fmt_shortest(P.a / to) + " " + fmt_shortest(P.b / to) +
                        (P.distribution == "triangular" ? " " + fmt_shortest(P.c / to) : std::string()) + ") on [" + fmt_shortest(P.min / to) + " " +
                        fmt_shortest(P.max / to) + "]";
    u += s.names[k] + "," + s.units[k] + "," + prior + "," + (k < pr.lsq.size() ? fmt_fixed(pr.lsq[k] / to, 6) : std::string()) + "," +
         (k < pr.lsq_sigma.size() ? fmt_fixed(pr.lsq_sigma[k] / to, 6) : std::string()) + "," + fmt_fixed(s.mean[k], 6) + "," + fmt_fixed(pr.sd[k] / to, 6) +
         "," + fmt_fixed(s.q025[k], 6) + "," + fmt_fixed(s.q975[k], 6) + "," + fmt_fixed(pr.rhat[k], 5) + "," + fmt_fixed(pr.ess[k], 1) + "," +
         std::to_string(pr.walkers) + "," + std::to_string(pr.steps) + "," + std::to_string(pr.burn_in) + "," + std::to_string(pr.thin) + "," +
         std::to_string(pr.seed) + "," + fmt_fixed(pr.acceptance, 4) + "," + (pr.converged ? "yes" : "no") + "\n";
  }
  write_file(dir / "posterior_summary.csv", u);
}

/// Stage 4A figures: anomaly maps per kind, and profiles along the sections (STX-TN-14 section 7).
void write_geophysics_figures(const desc::Description& d, const Model& m, const std::vector<std::pair<std::string, std::string>>& keys,
                              const PotentialFieldResult& pf, const std::filesystem::path& folder) {
  auto c = figure_context(d, m, keys);
  const auto spec = potential_field_spec(d, m.info());
  for (int kind = 0; kind < 2; ++kind) {
    const std::string k = kind == 0 ? "gravity" : "magnetics", unit = kind == 0 ? "mGal" : "nT";
    std::vector<fig::AnomalyStation> st;
    double zmean = 0;
    for (auto& s : pf.stations)
      if (s.kind == k) {
        st.push_back({s.p.x, s.p.y, s.observed, s.modelled, s.sigma});
        zmean += s.p.z;
      }
    if (st.empty()) continue;
    zmean /= static_cast<double>(st.size());
    const std::string ev = "Densities and susceptibilities as stated in the description; station data as delivered (no corrections applied by Stratonix).";
    fig::write_figure(fig::anomaly_map_figure(c, k, unit, st, kind == 0 ? pf.chi2_gravity : pf.chi2_magnetics), folder, fig::figure_basename(c, k + "_map", ""),
                      c.model_name + " - " + k + " anomaly", "Observed and modelled " + k + " anomaly at the stations, and the residual.", ev);
    // profiles along the description's sections (modelled at the mean station elevation)
    const auto& dom = m.info().domain;
    const double band = std::max(spec.cell, 0.05 * std::max(dom.sx, dom.sy));
    for (auto& sd : d.sections) {
      Polyline pl;
      for (std::size_t r = 0; r < sd.polyline.rows(); ++r) pl.push_back({sd.polyline.si(2 * r), sd.polyline.si(2 * r + 1), 0});
      if (pl.size() < 2) continue;
      fig::AnomalyProfile prof;
      prof.section = sd.name;
      prof.band = band;
      std::vector<double> cum{0};
      for (std::size_t i = 1; i < pl.size(); ++i) cum.push_back(cum.back() + std::hypot(pl[i].x - pl[i - 1].x, pl[i].y - pl[i - 1].y));
      std::vector<Vec3> pts;
      const int np = 120;
      for (int q = 0; q <= np; ++q) {
        const double dd = cum.back() * q / np;
        std::size_t sg = 1;
        while (sg + 1 < pl.size() && cum[sg] < dd) ++sg;
        const double t = cum[sg] > cum[sg - 1] ? (dd - cum[sg - 1]) / (cum[sg] - cum[sg - 1]) : 0;
        pts.push_back({pl[sg - 1].x + t * (pl[sg].x - pl[sg - 1].x), pl[sg - 1].y + t * (pl[sg].y - pl[sg - 1].y), zmean});
        prof.s.push_back(dd);
      }
      prof.modelled = kind == 0 ? gravity_anomaly(m, spec, pts) : magnetic_anomaly(m, spec, pts);
      for (auto& a : st) {
        // nearest point on the line
        double best = 1e300, bs = 0;
        for (std::size_t i = 1; i < pl.size(); ++i) {
          const double ex = pl[i].x - pl[i - 1].x, ey = pl[i].y - pl[i - 1].y, l2 = ex * ex + ey * ey;
          const double t = l2 > 0 ? std::clamp(((a.x - pl[i - 1].x) * ex + (a.y - pl[i - 1].y) * ey) / l2, 0.0, 1.0) : 0;
          const double dx = pl[i - 1].x + t * ex - a.x, dy = pl[i - 1].y + t * ey - a.y, dist = std::hypot(dx, dy);
          if (dist < best) { best = dist; bs = cum[i - 1] + t * std::sqrt(l2); }
        }
        if (best <= band) prof.observed.push_back({bs, a.observed});
      }
      fig::write_figure(fig::anomaly_profile_figure(c, k, unit, prof), folder, fig::figure_basename(c, k + "_profile", sd.name),
                        c.model_name + " - " + k + " profile " + sd.name, "Modelled " + k + " anomaly along the section, with the observed stations near it.", ev);
    }
  }
}

/// Stage 3A figures: fitted parameters, and the data misfit per horizon where uncertainties are
/// stated or the model was conditioned (STX-TN-11 section 6).
void write_data_figures(const desc::Description& d, const Model& m, const std::vector<std::pair<std::string, std::string>>& keys,
                        const std::filesystem::path& folder) {
  auto* cr = static_cast<const ConditioningResult*>(m.capability("conditioning"));
  auto* mis = static_cast<const std::vector<DataMisfit>*>(m.capability("data_misfit"));
  const bool with_sigma = mis && std::any_of(mis->begin(), mis->end(), [](const DataMisfit& r) { return !std::isnan(r.sigma); });
  if (!cr && !with_sigma) return;
  auto c = figure_context(d, m, keys);
  if (cr) {
    std::vector<fig::FittedRow> rows;
    for (auto& p : cr->parameters) rows.push_back({p.path, p.unit, p.start, p.value, p.sigma, p.min, p.max});
    fig::write_figure(fig::conditioning_figure(c, rows, cr->chi2_reduced, cr->residuals.size()), folder, fig::figure_basename(c, "conditioning", ""),
                      c.model_name + " - fitted parameters", "Parameters fitted to the data, with their admissible ranges and uncertainties.",
                      "Fitted values are labelled calibrated; the data controlling each are listed in influence.csv.");
  }
  if (auto* pr = static_cast<const PosteriorResult*>(m.capability("posterior"))) {
    auto ps = posterior_samples(d, *pr);
    const std::string ev = "Posterior of conditioning.method: bayesian; seed " + std::to_string(pr->seed) + "; samples in posterior_samples.csv.";
    fig::write_figure(fig::posterior_trace_figure(c, ps), folder, fig::figure_basename(c, "posterior_trace", ""), c.model_name + " - posterior traces",
                      "Walker paths of every fitted parameter after burn-in, with split R-hat and effective sample size.", ev);
    fig::write_figure(fig::posterior_corner_figure(c, ps), folder, fig::figure_basename(c, "posterior_corner", ""), c.model_name + " - posterior",
                      "Marginal and pairwise posterior of the fitted parameters.", ev);
  }
  std::map<std::string, std::vector<fig::MisfitPoint>> per;
  for (auto& r : *mis)
    if (r.kind != "orientation" && r.kind != "gravity" && r.kind != "magnetics") per[r.label].push_back({r.p.x, r.p.y, r.misfit, r.normalised});
  for (auto& [h, pts] : per)
    fig::write_figure(fig::misfit_map_figure(c, h, pts), folder, fig::figure_basename(c, "data-misfit", h), c.model_name + " - data misfit on " + h,
                      "Misfit of the model at every datum of " + h + ".", "Data are never altered; their misfit is reported, normalised by the stated uncertainty.");
}

std::vector<std::pair<std::string, std::string>> grid_rows(const GridProducts& gp) {
  std::vector<std::pair<std::string, std::string>> k;
  const auto& g = gp.grid;
  std::size_t act = 0;
  double vol = 0;
  for (int kk = 0; kk < g.nz; ++kk)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i)
        if (g.actnum[g.cell(i, j, kk)]) { ++act; vol += repr::cell_volume(g.corners(i, j, kk)); }
  std::size_t fp = 0;
  for (auto f : g.pillar_fault) fp += f >= 0;
  k.push_back({"Grid", std::to_string(g.nx) + " x " + std::to_string(g.ny) + " x " + std::to_string(g.nz) + " cells (" + fmt_fixed(g.dx, 1) + " x " +
                           fmt_fixed(g.dy, 1) + " m), " + std::to_string(act) + " active, " + std::to_string(fp) + " fault pillars"});
  k.push_back({"Grid bulk volume (active cells)", fmt_fixed(vol / 1e6, 3) + " million m3"});
  for (auto& p : gp.props.props)
    if (p.missing) k.push_back({"Grid property " + p.name + " missing", std::to_string(p.missing) + " cells: " + p.missing_reason});
  if (gp.has_trans) {
    std::size_t n = 0, missing = 0, limit = 0;
    double lo = 1e300;
    for (auto& f : gp.trans.faces) {
      ++n;
      if (!f.missing.empty()) ++missing;
      else if (!f.note.empty()) ++limit;
      else lo = std::min(lo, f.mult);
    }
    k.push_back({"Fault transmissibility (Manzocchi et al. 1999)", std::to_string(n) + " fault faces: " + std::to_string(n - missing - limit) +
                                                                     " computed (lowest multiplier " + (lo < 1e299 ? fmt_shortest(std::round(lo * 1e4) / 1e4) : std::string("-")) + "), " +
                                                                     std::to_string(limit) + " at 1 (no displacement), " + std::to_string(missing) + " missing"});
  }
  return k;
}

/// Grid tables and figures (stage 3B): grid/cells.csv, grid/fault_faces.csv and the grid figures.
void write_grid_outputs(const desc::Description& d, const Model& m, const GridProducts& gp, const std::vector<std::pair<std::string, std::string>>& keys,
                        const std::filesystem::path& out_dir, bool figures) {
  const auto& g = gp.grid;
  const auto& inf = m.info();
  std::string t = "i,j,k,active,zone,block,volume_m3,x,y,z";
  for (auto& p : gp.props.props)
    if (p.name != "ZONE" && p.name != "BLOCK") t += "," + p.name;
  t += "\n";
  for (int k = 0; k < g.nz; ++k)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i) {
        auto id = g.cell(i, j, k);
        auto cc = g.corners(i, j, k);
        Vec3 ctr = repr::cell_centre(cc);
        t += std::to_string(i + 1) + "," + std::to_string(j + 1) + "," + std::to_string(k + 1) + "," + std::to_string(g.actnum[id]) + "," +
             (g.zone[id] >= 0 ? inf.units[static_cast<std::size_t>(g.zone[id])].name : std::string()) + "," + std::to_string(g.block[id]) + "," +
             fmt_fixed(repr::cell_volume(cc), 3) + "," + fmt_fixed(ctr.x, 3) + "," + fmt_fixed(ctr.y, 3) + "," + fmt_fixed(ctr.z, 3);
        for (auto& p : gp.props.props)
          if (p.name != "ZONE" && p.name != "BLOCK") t += "," + (std::isnan(p.v[id]) ? std::string() : (p.integer ? std::to_string(std::llround(p.v[id])) : fmt_shortest(p.v[id])));
        t += "\n";
      }
  std::filesystem::create_directories(out_dir / "grid");
  write_file(out_dir / "grid" / "cells.csv", t);
  if (gp.has_trans) {
    std::string f = "fault,i,j,k,face,x,y,z,sgr,thickness_m,displacement_m,k_fault_mD,k_i_mD,k_j_mD,L_i_m,L_j_m,multiplier,missing,note\n";
    for (auto& ff : gp.trans.faces)
      f += gp.fault_names[static_cast<std::size_t>(ff.fault)] + "," + std::to_string(ff.i + 1) + "," + std::to_string(ff.j + 1) + "," + std::to_string(ff.k + 1) +
           "," + ff.dir + "+," + fmt_fixed(ff.centre.x, 2) + "," + fmt_fixed(ff.centre.y, 2) + "," + fmt_fixed(ff.centre.z, 2) + "," + fmt_shortest(ff.sgr) + "," +
           fmt_shortest(ff.thickness) + "," + fmt_shortest(ff.displacement) + "," + fmt_shortest(ff.k_fault) + "," + fmt_shortest(ff.k_i) + "," +
           fmt_shortest(ff.k_j) + "," + fmt_shortest(ff.L_i) + "," + fmt_shortest(ff.L_j) + "," + fmt_shortest(ff.mult) + ",\"" + ff.missing + "\",\"" + ff.note + "\"\n";
    write_file(out_dir / "grid" / "fault_faces.csv", f);
  }
  if (!figures) return;
  auto c = figure_context(d, m, keys);
  const std::string ev = std::string("Evidence basis: ") + basis_label(inf.basis) + ".";
  auto folder = out_dir / "figures";
  auto put = [&](const fig::Figure& f, const std::string& kind, const std::string& name, const std::string& t2, const std::string& cap) {
    fig::write_figure(f, folder, fig::figure_basename(c, kind, name), c.model_name + " - " + t2, cap, ev);
  };
  put(fig::grid_section_figure(c, g, nullptr, true, g.ny / 2), "grid-section", "row_zone", "grid section along a row",
      "Corner-point grid along the middle row, cells by zone, fault pillars in red.");
  put(fig::grid_section_figure(c, g, nullptr, false, g.nx / 2), "grid-section", "column_zone", "grid section along a column",
      "Corner-point grid along the middle column, cells by zone, fault pillars in red.");
  if (const auto* b = gp.props.find("BLOCK"))
    put(fig::grid_layer_figure(c, g, b, 0), "grid-layer", "1_block", "fault blocks on layer 1", "Fault-block code of the cells of layer 1.");
  for (auto& p : gp.props.props) {
    if (p.integer) continue;
    const std::string nm = io::petrel_name(p.name);
    put(fig::grid_section_figure(c, g, &p, true, g.ny / 2), "grid-section", "row_" + nm, "grid section by " + p.name,
        p.name + " (" + p.unit + ", " + p.evidence + ") along the middle row; grey cells are missing.");
    put(fig::grid_layer_figure(c, g, &p, g.nz / 2), "grid-layer", std::to_string(g.nz / 2 + 1) + "_" + nm, p.name + " on a layer",
        p.name + " (" + p.unit + ") on layer " + std::to_string(g.nz / 2 + 1) + "; grey cells are missing.");
  }
  if (gp.has_trans)
    put(fig::grid_layer_figure(c, g, &gp.trans.multx, g.nz / 2), "grid-layer", std::to_string(g.nz / 2 + 1) + "_MULTX", "fault transmissibility multipliers",
        "MULTX (x+ faces) on layer " + std::to_string(g.nz / 2 + 1) + ": 1 away from faults; grey where a multiplier is missing.");
  put(fig::grid_view3d_figure(c, g, 210, 30), "grid-3d", "SW", "grid 3D view from the south-west", "Outer faces of the active cells and their faces across faults, by zone.");
  put(fig::grid_view3d_figure(c, g, 135, 30), "grid-3d", "SE", "grid 3D view from the south-east", "Outer faces of the active cells and their faces across faults, by zone.");
}

BuildOutcome run_build(const std::filesystem::path& description, const std::filesystem::path& out_dir,
                       const BuildOptions& opt, const std::string& command) {
  auto d = desc::load_description(description);
  ensure_new_output_dir(out_dir);
  BuildOutcome o = build_only(d, opt);
  write_file(out_dir / "model_summary.txt", summary_text(d, *o.model, o.log));
  if (auto mc = data_misfit_csv(*o.model); !mc.empty()) write_file(out_dir / "data_misfit.csv", mc);
  // stage 3A: fitted parameters and the data that control them (STX-TN-11 section 6)
  if (auto* cr = static_cast<const ConditioningResult*>(o.model->capability("conditioning"))) {
    std::string t = "path,unit,start,fitted,sigma,min,max,evidence,source\n", inf = "path,rank,input,horizon,kind,x,y,z,share_of_variance,leverage\n";
    for (auto& p : cr->parameters) {
      const desc::Param* fp = cr->fitted.param_by_name(p.path);
      t += p.path + "," + p.unit + "," + fmt_shortest(p.start) + "," + fmt_shortest(p.value) + "," + fmt_shortest(p.sigma) + "," + fmt_shortest(p.min) +
           "," + fmt_shortest(p.max) + ",calibrated,\"" + (fp ? fp->source : std::string()) + "\"\n";
      for (std::size_t q = 0; q < p.controlling.size(); ++q) {
        auto& r = cr->residuals[p.controlling[q].first];
        inf += p.path + "," + std::to_string(q + 1) + "," + r.input + "," + r.label + "," + r.kind + "," + fmt_fixed(r.p.x, 2) + "," + fmt_fixed(r.p.y, 2) +
               "," + fmt_fixed(r.p.z, 2) + "," + fmt_fixed(p.controlling[q].second, 4) + "," + fmt_fixed(cr->leverage[p.controlling[q].first], 4) + "\n";
      }
    }
    write_file(out_dir / "conditioning.csv", t);
    write_file(out_dir / "influence.csv", inf);
  }
  if (auto* pr = static_cast<const PosteriorResult*>(o.model->capability("posterior"))) write_posterior_tables(d, *pr, out_dir);
  CheckContext ctx(*o.model, d, make_well_paths(d), d.grid_spacing.si());
  // stage 4A (STX-TN-14 section 7): gravity and magnetic stations, observed against modelled
  PotentialFieldResult pf;
  if (d.geophysics.present && (!d.geophysics.gravity.input.empty() || !d.geophysics.magnetics.input.empty())) {
    BuildInputs in;
    io::load_build_inputs(d, in);
    pf = potential_field_misfit(*o.model, d, in);
    std::string t = "kind,input,x,y,z,observed,modelled,residual,sigma,unit\n";
    for (auto& st : pf.stations)
      t += st.kind + "," + st.input + "," + fmt_fixed(st.p.x, 2) + "," + fmt_fixed(st.p.y, 2) + "," + fmt_fixed(st.p.z, 2) + "," + fmt_fixed(st.observed, 5) + "," +
           fmt_fixed(st.modelled, 5) + "," + fmt_fixed(st.modelled - st.observed, 5) + "," + (std::isfinite(st.sigma) ? fmt_fixed(st.sigma, 5) : std::string()) + "," +
           (st.kind == "gravity" ? "mGal" : "nT") + "\n";
    write_file(out_dir / "geophysics_stations.csv", t);
    ctx.potential_fields = &pf;
  }
  // stage 4B (STX-TN-15 section 4.2): volumetric restoration of the model to a horizon (C31)
  VolumetricResult vr;
  if (d.restoration.volumetric) {
    const auto& R = d.restoration;
    VolumetricSpec vs;
    vs.horizon = R.horizon;
    vs.datum_z = -R.datum.si();
    vs.young = R.young.si();
    vs.poisson = R.poisson.si();
    vs.pin_given = true;
    vs.pin_x = R.pin.si(0);
    vs.pin_y = R.pin.si(1);
    if (R.mesh_size.present) vs.mesh_size = R.mesh_size.si();
    if (R.max_strain.present) vs.max_strain = R.max_strain.si();
    const auto& units = o.model->info().units;
    vs.unit_props.assign(units.size(), {vs.young, vs.poisson});
    for (auto& u : R.units)
      for (std::size_t k = 0; k < units.size(); ++k)
        if (units[k].name == u.unit) vs.unit_props[k] = {u.young.si(), u.poisson.si()};
    vr = volumetric_restore(*o.model, vs);
    ctx.volumetric = &vr;
    std::string t = "present_x,present_y,present_z,restored_x,restored_y,restored_z,dilatation,max_shear,max_principal,bed_stretch,unit\n";
    for (auto& sm : vr.samples)
      if (sm.at_horizon)
        t += fmt_fixed(sm.present.x, 2) + "," + fmt_fixed(sm.present.y, 2) + "," + fmt_fixed(sm.present.z, 2) + "," + fmt_fixed(sm.restored.x, 2) + "," +
             fmt_fixed(sm.restored.y, 2) + "," + fmt_fixed(sm.restored.z, 2) + "," + fmt_fixed(sm.dilatation, 6) + "," + fmt_fixed(sm.max_shear, 6) + "," +
             fmt_fixed(sm.max_abs, 6) + "," + fmt_fixed(sm.bed_stretch, 6) + "," + (sm.unit >= 0 && static_cast<std::size_t>(sm.unit) < units.size() ? units[static_cast<std::size_t>(sm.unit)].name : "") + "\n";
    write_file(out_dir / "restoration_strain.csv", t);
    o.log.value("restoration.volumetric", vr.ok ? "max strain " + fmt_fixed(vr.max_strain, 5) + ", flatness " + fmt_fixed(vr.flatness_rms, 3) + " m rms, " +
                                                      std::to_string(vr.tetrahedra) + " tetrahedra, " + std::to_string(vr.cg_iterations) + " CG iterations"
                                                : "not restored: " + vr.reason);
  }
  // stage 4B (section 4.4): tectonic subsidence at the wells (and maps with flexure)
  if (d.backstripping.present)
    if (auto* st = static_cast<const Structure*>(o.model->capability("structure")); st && st->compaction_mode() != 0) {
      const BackstripSpec bs = backstrip_spec(d);
      std::string t = "location,step,age_Ma,sediment_m,mean_density_kg_m3,water_depth_m,sea_level_m,tectonic_airy_m,tectonic_flexural_m\n";
      BackstripMaps maps;
      if (bs.te > 0) maps = backstrip_maps(*st, bs, o.model->info().domain, 32);
      for (auto& w : make_well_paths(d))
        for (auto& ts : backstrip_column(*st, bs, w.x(), w.y())) {
          double fl = kNaN;
          for (std::size_t q = 0; q < maps.labels.size(); ++q)
            if (maps.labels[q] == ts.label) fl = backstrip_map_value(maps, maps.flexural[q], o.model->info().domain, w.x(), w.y());
          t += w.name() + "," + ts.label + "," + (std::isfinite(ts.age) ? fmt_shortest(ts.age) : std::string()) + "," + fmt_fixed(ts.sediment, 2) + "," +
               fmt_fixed(ts.density, 1) + "," + fmt_fixed(ts.water_depth, 2) + "," + fmt_fixed(ts.sea_level, 2) + "," + fmt_fixed(ts.airy, 2) + "," +
               (std::isfinite(fl) ? fmt_fixed(fl, 2) : std::string()) + "\n";
        }
      write_file(out_dir / "tectonic_subsidence.csv", t);
      if (!maps.labels.empty()) {
        std::string mt = "step,x,y,tectonic_airy_m,tectonic_flexural_m\n";
        const auto& dom = o.model->info().domain;
        for (std::size_t q = 0; q < maps.labels.size(); ++q)
          for (std::size_t j = 0; j < maps.n; ++j)
            for (std::size_t i = 0; i < maps.n; ++i)
              mt += maps.labels[q] + "," + fmt_fixed(dom.x0 + (i + 0.5) * dom.sx / maps.n, 1) + "," + fmt_fixed(dom.y0 + (j + 0.5) * dom.sy / maps.n, 1) + "," +
                    fmt_fixed(maps.airy[q][j * maps.n + i], 2) + "," + fmt_fixed(maps.flexural[q][j * maps.n + i], 2) + "\n";
        write_file(out_dir / "tectonic_subsidence_maps.csv", mt);
      }
    }
  // stage 4B (section 4.5): eroded section at unconformities from porosity logs
  if (!d.eroded_sections.empty()) {
    std::string t = "name,well,unconformity,unconformity_depth_m,samples,eroded_m,residual_rms_m,vitrinite_estimate_m,note\n";
    for (auto& e : d.eroded_sections) {
      const desc::InputDesc* in = nullptr;
      for (auto& x : d.inputs)
        if (x.id == e.log) in = &x;
      const desc::WellDesc* w = d.well(e.well);
      if (!in || !w) continue;
      const auto tab = io::read_csv(read_file(in->resolved));
      const double zu = horizon_elevation(*o.model, d.horizon_index(e.unconformity), w->head.si(0), w->head.si(1));
      const double du = -zu;  // depth of the unconformity below the datum (log depths are below the datum)
      std::vector<std::pair<double, double>> log;
      for (auto& row : tab.rows)
        if (row.size() >= 2) {
          const double z = std::stod(row[0]) - du, phi = std::stod(row[1]);
          if (z > 0 && phi > 0) log.push_back({z, phi});
        }
      std::string note;
      double est = kNaN, rms = kNaN;
      if (log.empty() || !std::isfinite(du)) {
        note = "no log samples below the unconformity";
      } else {
        const double phi0 = e.surface_porosity.si(), c = e.decay.si();
        double sum = 0, ss = 0;
        for (auto& [z, phi] : log) sum += -std::log(phi / phi0) / c - z;
        const double raw = sum / static_cast<double>(log.size());
        for (auto& [z, phi] : log) {
          const double r = -std::log(phi / phi0) / c - z - raw;
          ss += r * r;
        }
        rms = std::sqrt(ss / static_cast<double>(log.size()));
        est = std::max(0.0, raw);
        note = raw < 0 ? "underburied against the normal trend (" + fmt_fixed(raw, 1) + " m): reported as zero erosion" : "Magara (1976) shift of the normal trend";
      }
      t += e.name + "," + e.well + "," + e.unconformity + "," + fmt_fixed(du, 2) + "," + std::to_string(log.size()) + "," +
           (std::isfinite(est) ? fmt_fixed(est, 1) : std::string()) + "," + (std::isfinite(rms) ? fmt_fixed(rms, 1) : std::string()) + "," +
           (e.vitrinite_estimate.present ? fmt_fixed(e.vitrinite_estimate.si(), 1) : std::string()) + "," + note + "\n";
      if (std::isfinite(est)) o.log.value("eroded_section." + e.name, fmt_fixed(est, 1) + " m (residual " + fmt_fixed(rms, 1) + " m rms)");
    }
    write_file(out_dir / "eroded_section.csv", t);
  }
  // stage 4B (sections 5 and 6): mechanics, fractures, heat and pressure
  CheckContext::BasinResults basin;
  std::vector<std::pair<std::string, Grid2D>> basin_maps;  // figures written after the checks
  BurialFields burial_fields;  // stage 4B: temperature and overpressure per well and step, drawn over the burial curves
  {
    const Model& m = *o.model;
    const auto& inf = m.info();
    const auto wells = make_well_paths(d);
    auto* st = static_cast<const Structure*>(m.capability("structure"));
    const bool columns = st && st->compaction_mode() != 0;
    auto unit_props = [&](ColumnProps& cp) {
      const std::size_t nu = inf.units.size();
      cp.conductivity.assign(nu, kNaN);
      cp.heat_production.assign(nu, kNaN);
      cp.perm_law.assign(nu, 0);
      cp.k0.assign(nu, kNaN);
      cp.b.assign(nu, kNaN);
      cp.phi0.assign(nu, kNaN);
      for (auto& u : d.heat.units) {
        const int k = model_unit_of(d, m, u.unit);
        if (k >= 0) { cp.conductivity[static_cast<std::size_t>(k)] = u.conductivity.si(); cp.heat_production[static_cast<std::size_t>(k)] = u.heat_production.si(); }
      }
      for (auto& u : d.pressure.units) {
        const int k = model_unit_of(d, m, u.unit);
        if (k < 0) continue;
        cp.perm_law[static_cast<std::size_t>(k)] = u.law == "kozeny_carman" ? 2 : 1;
        cp.k0[static_cast<std::size_t>(k)] = u.k0.si();
        cp.b[static_cast<std::size_t>(k)] = u.b.present ? u.b.si() : kNaN;
        cp.phi0[static_cast<std::size_t>(k)] = u.phi0.si();
      }
    };
    const Domain& dom = inf.domain;
    const std::size_t nm = 16;  // map columns per side
    auto map_grid = [&]() { return Grid2D(GridSpec{dom.x0 + dom.sx / nm / 2, dom.y0 + dom.sy / nm / 2, dom.sx / nm, dom.sy / nm, nm, nm}); };
    if (d.heat.present && columns) {
      ThermalColumnSpec ts;
      ts.t_surface = d.heat.surface_temperature.si();
      ts.q_base = d.heat.basal_heat_flow.si();
      for (auto& [a, q] : d.heat.history) ts.q_history.push_back({a.si(), q.si()});
      ts.k_water = d.heat.water_conductivity.si();
      if (d.heat.column_depth.present) ts.column_depth = d.heat.column_depth.si();
      ts.vitrinite = d.heat.vitrinite;
      unit_props(ts.props);
      std::string t = "location,step,age_Ma,horizon,temperature_degC" + std::string(ts.vitrinite ? ",ro_percent" : "") + "\n";
      std::string prof = "location,depth_m,temperature_degC\n";
      for (auto& w : wells) {
        auto th = thermal_history(*st, ts, w.x(), w.y());
        basin.heat_balance.push_back({w.name(), th.balance_error});
        burial_fields[w.name()].push_back({th.labels, {"temperature", "temperature", "degC", th.temperature}});
        for (std::size_t k = 0; k < th.labels.size(); ++k)
          for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
            if (!std::isfinite(th.temperature[k][h])) continue;
            t += w.name() + "," + th.labels[k] + "," + (std::isfinite(th.ages[k]) ? fmt_shortest(th.ages[k]) : std::string()) + "," + inf.horizons[h].name + "," +
                 fmt_fixed(th.temperature[k][h], 2) + (ts.vitrinite ? "," + (std::isfinite(th.ro[k][h]) ? fmt_fixed(th.ro[k][h], 3) : std::string()) : "") + "\n";
          }
        for (std::size_t i = 0; i < th.present_depth.size(); i += 4) prof += w.name() + "," + fmt_fixed(th.present_depth[i], 1) + "," + fmt_fixed(th.present_temperature[i], 2) + "\n";
        o.log.value("heat." + w.name() + ".surface_heat_flow_mW_m2", fmt_fixed(1e3 * th.surface_flow_present, 2));
      }
      write_file(out_dir / "temperature_history.csv", t);
      write_file(out_dir / "temperature_profile.csv", prof);
      // present temperature (and %Ro) at each horizon on a coarse map
      std::vector<Grid2D> gt(inf.horizons.size(), map_grid()), gr(inf.horizons.size(), map_grid());
      std::vector<ThermalHistory> cols(nm * nm);
      parallel_for(nm * nm, 4, [&](std::size_t a0, std::size_t a1) {
        for (std::size_t q = a0; q < a1; ++q) cols[q] = thermal_history(*st, ts, gt[0].spec.x(q % nm), gt[0].spec.y(q / nm));
      });
      for (std::size_t q = 0; q < nm * nm; ++q)
        for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
          const auto& c = cols[q];
          gt[h].v[q] = c.temperature.empty() ? kNaN : c.temperature.back()[h];
          gr[h].v[q] = c.ro.empty() ? kNaN : c.ro.back()[h];
        }
      for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
        basin_maps.push_back({"temperature:" + inf.horizons[h].name, gt[h]});
        if (ts.vitrinite) basin_maps.push_back({"ro:" + inf.horizons[h].name, gr[h]});
      }
    }
    if (d.pressure.present && columns) {
      PressureColumnSpec ps;
      ps.viscosity = d.pressure.viscosity.si();
      ps.rho_w = d.pressure.rho_water.si();
      ps.rho_g = d.pressure.rho_grain.si();
      unit_props(ps.props);
      std::string t = "location,step,age_Ma,horizon,overpressure_MPa\n", prof = "location,depth_m,pore_MPa,hydrostatic_MPa,lithostatic_MPa\n";
      for (auto& w : wells) {
        auto ph = pressure_history(*st, ps, w.x(), w.y());
        {
          auto mpa = ph.overpressure;
          for (auto& row : mpa)
            for (double& v : row) v /= 1e6;
          burial_fields[w.name()].push_back({ph.labels, {"overpressure", "overpressure", "MPa", mpa}});
        }
        for (std::size_t k = 0; k < ph.labels.size(); ++k)
          for (std::size_t h = 0; h < inf.horizons.size(); ++h)
            if (std::isfinite(ph.overpressure[k][h]))
              t += w.name() + "," + ph.labels[k] + "," + (std::isfinite(ph.ages[k]) ? fmt_shortest(ph.ages[k]) : std::string()) + "," + inf.horizons[h].name + "," +
                   fmt_fixed(ph.overpressure[k][h] / 1e6, 4) + "\n";
        CheckContext::BasinResults::Bounds b;
        b.location = w.name();
        for (std::size_t i = 0; i < ph.present_depth.size(); ++i) {
          const double p = ph.present_pore[i], hy = ph.present_hydrostatic[i], li = ph.present_lithostatic[i];
          if (p < hy - 1e3) { ++b.below; b.worst = std::max(b.worst, hy - p); }
          if (i > 0 && p > li + 1e3) { ++b.above; b.worst = std::max(b.worst, p - li); }
          if (i % 4 == 0)
            prof += w.name() + "," + fmt_fixed(ph.present_depth[i], 1) + "," + fmt_fixed(p / 1e6, 4) + "," + fmt_fixed(hy / 1e6, 4) + "," + fmt_fixed(li / 1e6, 4) + "\n";
        }
        basin.pressure.push_back(b);
      }
      write_file(out_dir / "pressure_history.csv", t);
      write_file(out_dir / "pressure_profile.csv", prof);
      std::vector<Grid2D> gp(inf.horizons.size(), map_grid());
      std::vector<PressureHistory> cols(nm * nm);
      parallel_for(nm * nm, 4, [&](std::size_t a0, std::size_t a1) {
        for (std::size_t q = a0; q < a1; ++q) cols[q] = pressure_history(*st, ps, gp[0].spec.x(q % nm), gp[0].spec.y(q / nm));
      });
      for (std::size_t q = 0; q < nm * nm; ++q)
        for (std::size_t h = 0; h < inf.horizons.size(); ++h) gp[h].v[q] = cols[q].overpressure.empty() ? kNaN : cols[q].overpressure.back()[h] / 1e6;
      for (std::size_t h = 0; h < inf.horizons.size(); ++h) basin_maps.push_back({"overpressure:" + inf.horizons[h].name, gp[h]});
    }
    // stress along the wells with the mechanical stratigraphy (stage 4B, I-033)
    if (d.mechanics.present && d.stress.present) {
      const auto sin = stress_input(d, m);
      std::string t = "location,depth_m,unit,sigma_v_MPa,pore_MPa,sigma_h_MPa,sigma_H_MPa\n";
      for (auto& w : wells)
        for (double z = 50; z < -dom.z_base; z += 50) {
          const int u = m.unit_at({w.x(), w.y(), -z});
          double S6[6], pp;
          repr::stress_at(sin, u, z, S6, pp);
          const double a = 0.5 * (S6[0] + S6[1]), r = std::sqrt(0.25 * (S6[0] - S6[1]) * (S6[0] - S6[1]) + S6[3] * S6[3]);
          t += w.name() + "," + fmt_fixed(z, 1) + "," + (u >= 0 && static_cast<std::size_t>(u) < inf.units.size() ? inf.units[static_cast<std::size_t>(u)].name : "") + "," +
               fmt_fixed(S6[2] / 1e6, 3) + "," + fmt_fixed(pp / 1e6, 3) + "," + fmt_fixed((a - r) / 1e6, 3) + "," + fmt_fixed((a + r) / 1e6, 3) + "\n";
        }
      write_file(out_dir / "stress_profile.csv", t);
    }
    // fracture intensity per unit from the curvature of its top horizon, and the network
    if (d.fractures.present) {
      std::vector<Disk> dfn_disks;
      const std::size_t nf = 60;
      GridSpec gs{dom.x0 + dom.sx / nf / 2, dom.y0 + dom.sy / nf / 2, dom.sx / nf, dom.sy / nf, nf, nf};
      std::string t = "unit,x,y,kmax_abs_1_per_m,p32_1_per_m\n";
      for (auto& fu : d.fractures.units) {
        const int ui = d.unit_index(fu.unit);
        if (ui < 0) continue;
        const int h = d.horizon_index(d.units[static_cast<std::size_t>(ui)].top);
        if (h < 0) continue;
        const auto cg = repr::curvature(m, h, gs, 0.0);
        FractureLaw law{fu.p32_0.si(), fu.a_k.si(), fu.a_e.present ? fu.a_e.si() : 0.0, fu.p32_max.si()};
        Grid2D g(gs);
        double sum = 0;
        std::size_t cnt = 0;
        for (std::size_t q = 0; q < g.v.size(); ++q) {
          const double k1 = cg.kmax.v[q], k2 = cg.kmin.v[q];
          if (!std::isfinite(k1) || !std::isfinite(k2)) { g.v[q] = kNaN; continue; }  // missing: no curvature here
          const double k = std::max(std::fabs(k1), std::fabs(k2));
          g.v[q] = fracture_intensity(k, 0.0, law);
          t += fu.unit + "," + fmt_fixed(gs.x(q % nf), 1) + "," + fmt_fixed(gs.y(q / nf), 1) + "," + fmt_fixed(k, 8) + "," + fmt_fixed(g.v[q], 5) + "\n";
          if (d.fractures.dfn.present && fu.unit == d.fractures.dfn.unit) {
            const auto& A = d.fractures.dfn.area;
            const double x = gs.x(q % nf), y = gs.y(q / nf);
            if (x >= A.si(0) && x <= A.si(2) && y >= A.si(1) && y <= A.si(3)) { sum += g.v[q]; ++cnt; }
          }
        }
        basin_maps.push_back({"p32:" + fu.unit, g});
        if (d.fractures.dfn.present && fu.unit == d.fractures.dfn.unit) {
          const auto& N = d.fractures.dfn;
          basin.dfn = true;
          basin.dfn_unit = fu.unit;
          basin.dfn_target = cnt ? sum / static_cast<double>(cnt) : 0.0;
          const double cx = 0.5 * (N.area.si(0) + N.area.si(2)), cy = 0.5 * (N.area.si(1) + N.area.si(3));
          const double zt = horizon_elevation(m, h, cx, cy), zb = horizon_elevation(m, d.horizon_index(d.units[static_cast<std::size_t>(ui)].base), cx, cy);
          BBox box;
          box.add({N.area.si(0), N.area.si(1), zb});
          box.add({N.area.si(2), N.area.si(3), zt});
          std::vector<FractureSet> sets;
          for (auto& s : N.sets) {
            const double tr = s.pole.si(0), pl = s.pole.si(1);
            FractureSet fs;
            fs.mean_pole = {std::sin(tr) * std::cos(pl), std::cos(tr) * std::cos(pl), -std::sin(pl)};
            fs.kappa = s.kappa.si();
            fs.exponent = s.exponent.si();
            fs.rmin = s.rmin.si();
            fs.rmax = s.rmax.si();
            fs.weight = s.weight.si();
            sets.push_back(fs);
          }
          if (!std::isfinite(zt) || !std::isfinite(zb) || !(zt > zb)) {
            basin.dfn_reason = "the unit has no thickness at the area's centre";
          } else {
            const auto lim = N.limit.present ? static_cast<std::size_t>(N.limit.si()) : std::size_t{1000000};
            auto r = sample_dfn(box, basin.dfn_target, sets, static_cast<std::uint64_t>(N.seed.si()), lim);
            basin.dfn_ok = r.ok;
            basin.dfn_reason = r.reason;
            if (r.ok) {
              basin.dfn_count = r.disks.size();
              basin.dfn_realised = realised_p32(r.disks, box);
              dfn_disks = r.disks;
              std::string dt = "x,y,z,normal_x,normal_y,normal_z,radius_m,set\n";
              for (auto& dk : r.disks)
                dt += fmt_fixed(dk.c.x, 2) + "," + fmt_fixed(dk.c.y, 2) + "," + fmt_fixed(dk.c.z, 2) + "," + fmt_fixed(dk.n.x, 5) + "," + fmt_fixed(dk.n.y, 5) + "," +
                      fmt_fixed(dk.n.z, 5) + "," + fmt_fixed(dk.r, 3) + "," + N.sets[static_cast<std::size_t>(dk.set)].name + "\n";
              write_file(out_dir / "dfn_disks.csv", dt);
            }
          }
        }
      }
      write_file(out_dir / "fracture_intensity.csv", t);
      // as grid properties when the description has a grid (stage 4B definition of done)
      if (const GridProducts* gp = ctx.grid_products()) {
        const auto& G = gp->grid;
        const GridSpec cs{G.x0 + 0.5 * G.dx, G.y0 + 0.5 * G.dy, G.dx, G.dy, static_cast<std::size_t>(G.nx), static_cast<std::size_t>(G.ny)};
        repr::CellProperty P{"P32", "1/m", "derived",
                             "fracture intensity law of the unit (fractures block) at the maximum absolute principal curvature of the unit's top horizon at the column centre"};
        P.v.assign(G.ncell(), kNaN);
        for (auto& fu : d.fractures.units) {
          const int ui = d.unit_index(fu.unit);
          if (ui < 0) continue;
          const int h = d.horizon_index(d.units[static_cast<std::size_t>(ui)].top);
          const int mu = model_unit_of(d, m, fu.unit);
          if (h < 0 || mu < 0) continue;
          const auto cg = repr::curvature(m, h, cs, 0.0);
          FractureLaw law{fu.p32_0.si(), fu.a_k.si(), fu.a_e.present ? fu.a_e.si() : 0.0, fu.p32_max.si()};
          for (int k = 0; k < G.nz; ++k)
            for (int j = 0; j < G.ny; ++j)
              for (int i = 0; i < G.nx; ++i) {
                const auto c = G.cell(i, j, k);
                if (!G.actnum[c] || G.zone[c] != mu) continue;
                const auto q = cs.index(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
                const double k1 = cg.kmax.v[q], k2 = cg.kmin.v[q];
                if (std::isfinite(k1) && std::isfinite(k2)) P.v[c] = fracture_intensity(std::max(std::fabs(k1), std::fabs(k2)), 0.0, law);
              }
        }
        for (double v : P.v) P.missing += std::isfinite(v) ? 0 : 1;
        if (P.missing) P.missing_reason = "inactive, outside the units with a fracture law, or no curvature at the column";
        ctx.add_grid_property(std::move(P));
        if (!dfn_disks.empty()) {
          repr::CellProperty R{"P32_DFN", "1/m", "derived", "realised P32 of the sampled network: disk area inside the cell's bounding box over its volume"};
          R.v.assign(G.ncell(), kNaN);
          const int mu = model_unit_of(d, m, basin.dfn_unit);
          const auto& A = d.fractures.dfn.area;
          // disks binned by column so each cell visits only the disks that can reach it
          std::vector<std::vector<std::size_t>> bins(G.ncell() / std::max(1, G.nz));
          for (std::size_t n = 0; n < dfn_disks.size(); ++n) {
            const auto& dk = dfn_disks[n];
            const int i0 = std::max(0, static_cast<int>(std::floor((dk.c.x - dk.r - G.x0) / G.dx))), i1 = std::min(G.nx - 1, static_cast<int>(std::floor((dk.c.x + dk.r - G.x0) / G.dx)));
            const int j0 = std::max(0, static_cast<int>(std::floor((dk.c.y - dk.r - G.y0) / G.dy))), j1 = std::min(G.ny - 1, static_cast<int>(std::floor((dk.c.y + dk.r - G.y0) / G.dy)));
            for (int j = j0; j <= j1; ++j)
              for (int i = i0; i <= i1; ++i) bins[static_cast<std::size_t>(i) + static_cast<std::size_t>(G.nx) * static_cast<std::size_t>(j)].push_back(n);
          }
          std::vector<Disk> local;
          for (int k = 0; k < G.nz; ++k)
            for (int j = 0; j < G.ny; ++j)
              for (int i = 0; i < G.nx; ++i) {
                const auto c = G.cell(i, j, k);
                if (!G.actnum[c] || G.zone[c] != mu) continue;
                BBox b;
                for (auto& v : G.corners(i, j, k)) b.add(v);
                const double xc = 0.5 * (b.lo.x + b.hi.x), yc = 0.5 * (b.lo.y + b.hi.y);
                if (xc < A.si(0) || xc > A.si(2) || yc < A.si(1) || yc > A.si(3)) continue;
                local.clear();
                for (auto n : bins[static_cast<std::size_t>(i) + static_cast<std::size_t>(G.nx) * static_cast<std::size_t>(j)]) local.push_back(dfn_disks[n]);
                R.v[c] = realised_p32(local, b);
              }
          for (double v : R.v) R.missing += std::isfinite(v) ? 0 : 1;
          if (R.missing) R.missing_reason = "inactive, outside the network's unit, or outside its area";
          ctx.add_grid_property(std::move(R));
        }
      }
    }
    if (!basin.heat_balance.empty() || !basin.pressure.empty() || basin.dfn) ctx.basin = &basin;
  }
  // stage 4C (STX-TN-16 section 1): logs, surveys, ages; upscaled logs as grid properties
  const GridProducts* gp_wd = ctx.grid_products();
  WellDataBuild wdb = well_data(d, *o.model, gp_wd ? &gp_wd->grid : nullptr);
  if (wdb.present) {
    for (auto& prop : wdb.upscaled) ctx.add_grid_property(prop);
    ctx.well_data = &wdb.results;
    write_file(out_dir / "well_data_report.txt", wdb.report);
    for (auto& [name, text] : wdb.csv) write_file(out_dir / name, text);
  }
  // stage 4D (STX-TN-17 section 1): geostatistical properties, conditioned to the upscaled logs
  GeostatBuild gsb;
  CheckContext::GeostatResults gsr;
  if (!d.grid.geostatistics.empty() && ctx.grid_products()) {
    const GridProducts* gp = ctx.grid_products();
    gsb = geostat_properties(d, *o.model, gp->grid, gp->props);
    for (auto& prop : gsb.props) ctx.add_grid_property(prop);
    for (auto& k : gsb.conditioning) gsr.conditioning.push_back({k.property, k.cells, k.mismatched, k.largest, k.exact});
    for (auto& st : gsb.statistics) gsr.statistics.push_back({st.property, st.what, st.target, st.got, st.tolerance});
    gsr.cells = gsb.cells;
    gsr.fallback = gsb.fallback;
    ctx.geostat = &gsr;
    write_file(out_dir / "geostatistics_report.txt", gsb.report);
    write_file(out_dir / "variograms.csv", gsb.csv);
  }
  // stage 4C (STX-TN-16 section 2.5): what the geological map gave the builder
  if (d.geological_map.present) {
    const io::ModelData md = io::geological_map_data(d);
    std::map<std::string, std::size_t> per;
    for (auto& [h, p] : md.contact_points) ++per[h];
    std::string t = "Stratonix geological map (stage 4C, STX-TN-16 section 2.5)\n";
    for (auto& line : md.report) t += line + "\n";
    t += "contact points per horizon (base of the younger unit):\n";
    for (auto& [h, n] : per) t += "  " + h + ": " + std::to_string(n) + "\n";
    std::string pts = "horizon,x,y,z\n";
    for (auto& [h, p] : md.contact_points) pts += h + "," + fmt_fixed(p.x, 1) + "," + fmt_fixed(p.y, 1) + "," + fmt_fixed(p.z, 1) + "\n";
    std::string ori = "horizon,x,y,z,dip_deg,dip_direction_deg,source\n";
    for (auto& m : md.orientations)
      ori += m.horizon + "," + fmt_fixed(m.p.x, 1) + "," + fmt_fixed(m.p.y, 1) + "," + fmt_fixed(m.p.z, 1) + "," + fmt_fixed(rad2deg(m.dip), 1) + "," +
             fmt_fixed(rad2deg(m.azimuth), 1) + "," + m.source + "\n";
    std::string ft = "fault,x,y,z\n";
    for (auto& [n, tr] : md.fault_traces)
      for (auto& q : tr) ft += n + "," + fmt_fixed(q.x, 1) + "," + fmt_fixed(q.y, 1) + "," + fmt_fixed(q.z, 1) + "\n";
    write_file(out_dir / "geological_map_report.txt", t);
    write_file(out_dir / "map_contact_points.csv", pts);
    write_file(out_dir / "map_orientations.csv", ori);
    write_file(out_dir / "map_fault_traces.csv", ft);
  }
  // stage 4C (STX-TN-16 section 2.3): the velocity model and its ties (C36)
  CheckContext::TimeDepthResults tdr;
  io::VelocityBuild vb = io::velocity_build(d);
  if (vb.present) {
    for (auto& t : vb.ties) tdr.ties.push_back({t.well, t.layer, t.residual, t.loo, t.excluded});
    if (d.velocity.tolerance.present) {
      tdr.tolerance = d.velocity.tolerance.si();
      tdr.tolerance_source = std::string("velocity_model.tolerance (") + desc::evidence_name(d.velocity.tolerance.evidence) + ")";
    } else {
      tdr.tolerance_source = "Stratonix default: 5 m";
    }
    ctx.time_depth = &tdr;
    write_file(out_dir / "velocity_report.txt", vb.report);
    write_file(out_dir / "velocity_calibration.csv", vb.csv);
  }
  CheckReport rep = run_checks(ctx);
  o.checks_passed = rep.passed();
  write_file(out_dir / "topology.json", ctx.contact_graph().json());  // stage 4A (STX-TN-14 section 6)
  write_file(out_dir / "checks_report.txt", rep.text());
  write_file(out_dir / "checks_report.csv", rep.csv());
  auto keys = key_numbers(d, *o.model, rep, ctx);
  auto rows = keys;
  for (auto& r : event_rows(d)) rows.push_back(r);
  for (auto& r : fault_rock_rows(d, *o.model)) rows.push_back(r);
  if (auto* cr = static_cast<const ConditioningResult*>(o.model->capability("conditioning"))) {
    for (auto& p : cr->parameters)
      rows.push_back({"Fitted " + p.path, fmt_fixed(p.value, 3) + " +- " + fmt_fixed(p.sigma, 3) + " " + p.unit + " (calibrated; start " + fmt_shortest(p.start) +
                                              ", range " + fmt_shortest(p.min) + " to " + fmt_shortest(p.max) + ")"});
    rows.push_back({"Conditioning", std::to_string(cr->residuals.size()) + " data, reduced chi2 " + fmt_fixed(cr->chi2_reduced, 3) + ", " +
                                        std::to_string(cr->iterations) + " iterations"});
  }
  if (auto* pr = static_cast<const PosteriorResult*>(o.model->capability("posterior"))) {
    auto ps = posterior_samples(d, *pr);
    for (std::size_t k = 0; k < ps.names.size(); ++k)
      rows.push_back({"Posterior " + ps.names[k], fmt_fixed(ps.mean[k], 3) + " " + ps.units[k] + " (95 % " + fmt_fixed(ps.q025[k], 3) + " to " +
                                                     fmt_fixed(ps.q975[k], 3) + "; R-hat " + fmt_fixed(pr->rhat[k], 4) + ", ESS " + fmt_fixed(pr->ess[k], 0) + ")"});
    rows.push_back({"Posterior sampler", std::to_string(pr->walkers) + " walkers x " + std::to_string(pr->steps) + " steps, seed " + std::to_string(pr->seed) +
                                             (pr->converged ? ", converged" : ", not converged")});
  }
  const GridProducts* gp = ctx.grid_products();
  if (gp)
    for (auto& r : grid_rows(*gp)) rows.push_back(r);
  write_file(out_dir / "summary_table.csv", summary_csv(rows));
  if (!o.checks_passed) {
    log(LogLevel::Error, "checks failed: exports and figures are not written (see checks_report.txt)");
  } else {
    if (opt.figures) write_figures(d, *o.model, keys, out_dir / "figures", ctx.surfaces(), &burial_fields);
    if (opt.figures && d.seismic.present) {
      // stage 4C (STX-TN-16 section 2.2): the SEG-Y volume along each section, in depth
      const desc::InputDesc* in = nullptr;
      for (auto& x : d.inputs)
        if (x.id == d.seismic.input) in = &x;
      io::SegyLayout lay{d.seismic.il, d.seismic.xl, d.seismic.x, d.seismic.y};
      const io::SegyVolume vol = io::read_segy(in->resolved, lay);
      std::string rep = "Stratonix seismic volume (stage 4C, STX-TN-16 section 2.2)\n" + in->resolved.generic_string() + "\n";
      for (auto& line : vol.report) rep += line + "\n";
      write_file(out_dir / "seismic_report.txt", rep);
      auto c = figure_context(d, *o.model, keys);
      const std::string how = d.seismic.domain == "time" ? "in time, drawn in depth through the velocity model (velocity_report.txt)" : "in depth";
      for (auto& sd : d.sections) {
        Polyline pl;
        for (std::size_t r = 0; r < sd.polyline.rows(); ++r) pl.push_back({sd.polyline.si(2 * r), sd.polyline.si(2 * r + 1), 0});
        auto s = repr::make_section(*o.model, sd.name, pl, std::max(10.0, c.grid_spacing / 2), 4.0);
        std::vector<float> amp(s.nz * s.ns, NAN);
        for (std::size_t q = 0; q < s.ns; ++q)
          for (std::size_t r = 0; r < s.nz; ++r) {
            const double z = s.z_top - (static_cast<double>(r) + 0.5) * s.dz, depth = -z;
            const double t = d.seismic.domain == "time" ? (vb.present ? io::depth_to_time_safe(vb.model, s.xy[q].x, s.xy[q].y, depth) : NAN) : depth;
            amp[r * s.ns + q] = static_cast<float>(vol.sample(s.xy[q].x, s.xy[q].y, t));
          }
        fig::write_figure(fig::seismic_section_figure(c, s, amp, how), out_dir / "figures", fig::figure_basename(c, "seismic_section", sd.name),
                          c.model_name + " - seismic along " + sd.name, "Seismic amplitude along " + sd.name + " with the model drawn over it.",
                          "Amplitude from the SEG-Y volume; see seismic_report.txt and velocity_report.txt.");
      }
    }
    if (opt.figures && wdb.present) {
      // stage 4C: well log panels and the correlation panel
      auto c = figure_context(d, *o.model, keys);
      for (auto& pan : wdb.panels)
        fig::write_figure(fig::well_log_figure(c, pan), out_dir / "figures", fig::figure_basename(c, "well_logs", pan.well),
                          c.model_name + " - logs of " + pan.well, "Logs of " + pan.well + " as read from its files, with tops, age data and the model's units.",
                          "Units and null values from the log files; see well_data_report.txt.");
      if (wdb.correlation)
        fig::write_figure(fig::correlation_figure(c, wdb.corr), out_dir / "figures", fig::figure_basename(c, "well_correlation", "wells"),
                          c.model_name + " - well correlation", "Wells side by side with the correlation lines and what each rests on.",
                          "Correlation by age where the horizon is dated and the wells have age data, else by tops, else by the model.");
    }
    if (opt.figures && gsb.present) {
      // stage 4D: variograms per property and direction, histograms against the wells
      auto c = figure_context(d, *o.model, keys);
      std::map<std::string, std::vector<fig::VariogramSeries>> by;
      std::map<std::string, std::pair<std::string, double>> meta;
      for (auto& vp : gsb.variograms) {
        fig::VariogramSeries s{vp.direction, vp.lag, vp.gamma, vp.pairs, {}, {}};
        double hmax = 0;
        for (double h : vp.lag) hmax = std::max(hmax, h);
        const double a = vp.azimuth * 3.14159265358979323846 / 180.0;
        for (int k = 0; k <= 60; ++k) {
          const double h = hmax * k / 60.0;
          s.model_lag.push_back(h);
          s.model_gamma.push_back(vp.model.gamma(vp.vertical ? Vec3{0, 0, h} : Vec3{h * std::sin(a), h * std::cos(a), 0}));
        }
        by[vp.property].push_back(std::move(s));
        meta[vp.property] = {vp.unit, vp.model.sill()};
      }
      for (auto& [prop, series] : by)
        fig::write_figure(fig::variogram_figure(c, prop, meta[prop].first, meta[prop].second, series), out_dir / "figures", fig::figure_basename(c, "variogram", prop),
                          c.model_name + " - variograms of " + prop, "Experimental and model variograms of " + prop + " in geological-time coordinates.",
                          "Model parameters labelled in the description; see geostatistics_report.txt and variograms.csv.");
      for (auto& hg : gsb.histograms)
        fig::write_figure(fig::property_histogram_figure(c, hg.property, hg.unit, hg.value_unit, hg.simulated, hg.wells, hg.codes), out_dir / "figures",
                          fig::figure_basename(c, "property_histogram", hg.property), c.model_name + " - " + hg.property + " against the wells",
                          "Distribution of " + hg.property + " in the grid cells against the wells' upscaled cells.",
                          "One realisation (seed in geostatistics_report.txt).");
    }
    if (opt.figures && !d.events.empty()) write_restoration(d, *o.model, out_dir / "restoration");
    if (opt.figures) write_data_figures(d, *o.model, keys, out_dir / "figures");
    if (opt.figures && !pf.stations.empty()) write_geophysics_figures(d, *o.model, keys, pf, out_dir / "figures");
    if (opt.figures && !basin_maps.empty()) {
      // stage 4B: present temperature, %Ro, overpressure per horizon and fracture intensity per unit
      auto c = figure_context(d, *o.model, keys);
      for (auto& [key, g] : basin_maps) {
        const auto colon = key.find(':');
        const std::string kind = key.substr(0, colon), name = key.substr(colon + 1);
        double vmin = 1e300, vmax = -1e300;
        for (double v : g.v)
          if (std::isfinite(v)) { vmin = std::min(vmin, v); vmax = std::max(vmax, v); }
        if (vmin > vmax) continue;
        const std::string what = kind == "temperature" ? "present temperature at " + name : kind == "ro" ? "vitrinite reflectance (EASY%Ro) at " + name
                                 : kind == "overpressure" ? "present overpressure at " + name : "fracture intensity P32 of " + name;
        const std::string bar = kind == "temperature" ? "Temperature (degC)" : kind == "ro" ? "%Ro" : kind == "overpressure" ? "Overpressure (MPa)" : "P32 (1/m)";
        const std::string sub = kind == "p32" ? "From the maximum absolute principal curvature of the unit's top horizon by the stated fracture law; grey where it is missing"
                                              : "1D columns through the burial history on a " + std::to_string(g.spec.nx) + " x " + std::to_string(g.spec.ny) + " map";
        fig::write_figure(fig::scalar_map_figure(c, what, sub, bar, g, vmin, vmax), out_dir / "figures", fig::figure_basename(c, kind + "_map", name),
                          c.model_name + " - " + what, what + ".", "Properties as stated in the description.");
      }
    }
    if (opt.figures && vr.ok) {
      // stage 4B: strain maps of the restored horizon (by present map position)
      auto c = figure_context(d, *o.model, keys);
      const auto& dom = o.model->info().domain;
      // one map cell per mesh cell, so every cell holds the tetrahedra the horizon crosses there
      const std::size_t n = std::clamp<std::size_t>(static_cast<std::size_t>(std::llround(std::max(dom.sx, dom.sy) / std::max(1.0, vr.mesh_spacing))), 8, 200);
      for (int kind = 0; kind < 2; ++kind) {
        Grid2D g(GridSpec{dom.x0 + dom.sx / n / 2, dom.y0 + dom.sy / n / 2, dom.sx / n, dom.sy / n, n, n});
        std::vector<double> sum(n * n, 0.0), cnt(n * n, 0.0);
        double vmin = 1e300, vmax = -1e300;
        for (auto& sm : vr.samples) {
          if (!sm.at_horizon) continue;
          const auto i = static_cast<std::size_t>(std::clamp((sm.present.x - dom.x0) / dom.sx * n, 0.0, n - 1.0));
          const auto j = static_cast<std::size_t>(std::clamp((sm.present.y - dom.y0) / dom.sy * n, 0.0, n - 1.0));
          const double v = kind == 0 ? sm.max_shear : sm.dilatation;
          sum[j * n + i] += v;
          cnt[j * n + i] += 1;
        }
        for (std::size_t q = 0; q < n * n; ++q) {
          g.v[q] = cnt[q] > 0 ? sum[q] / cnt[q] : kNaN;
          if (cnt[q] > 0) { vmin = std::min(vmin, g.v[q]); vmax = std::max(vmax, g.v[q]); }
        }
        const std::string what = kind == 0 ? "maximum shear" : "dilatation";
        fig::write_figure(fig::scalar_map_figure(c, "restoration strain of " + d.restoration.horizon + ": " + what,
                                                 "Volumetric restoration to " + fmt_fixed(d.restoration.datum.si(), 0) + " m depth (linear elasticity, faults sliding); "
                                                 "Green-Lagrange strain of the tetrahedra the horizon crosses, mean per cell",
                                                 kind == 0 ? "Maximum shear" : "Dilatation", g, vmin, vmax),
                          out_dir / "figures", fig::figure_basename(c, "restoration_strain", d.restoration.horizon + (kind == 0 ? "_shear" : "_dilatation")),
                          c.model_name + " - restoration strain (" + what + ")", "Strain of the volumetric restoration at the restored horizon.",
                          "Elastic properties as stated in the description; strain is reported, never hidden.");
      }
    }
    if (gp) write_grid_outputs(d, *o.model, *gp, keys, out_dir, opt.figures);
    auto ids = resolve_export_ids(d, opt.exports);
    if (!ids.empty()) write_exports(d, *o.model, ids, out_dir / "exports", 0, gp, ctx.surfaces());
  }
  Manifest man = Manifest::from_description(d, command);
  for (auto& [k, v] : o.log.values) man.extra.emplace_back(k, v);
  man.add_outputs(out_dir);
  man.write(out_dir / "manifest.json");
  return o;
}

BuildOutcome run_section(const std::filesystem::path& description, const std::filesystem::path& out_dir, const std::string& name,
                         const Polyline& line, const BuildOptions& opt) {
  auto d = desc::load_description(description);
  ensure_new_output_dir(out_dir);
  BuildOutcome o = build_only(d, opt);
  CheckContext ctx(*o.model, d, make_well_paths(d), d.grid_spacing.si());
  CheckReport rep = run_checks(ctx);
  o.checks_passed = rep.passed();
  write_file(out_dir / "checks_report.txt", rep.text());
  if (!o.checks_passed) {
    log(LogLevel::Error, "checks failed: the section is not drawn (see checks_report.txt)");
  } else {
    auto keys = key_numbers(d, *o.model, rep, ctx);
    auto c = figure_context(d, *o.model, keys);
    Polyline pl = line;
    std::vector<double> det;
    std::string nm = name.empty() ? "line" : name;
    if (pl.empty()) {
      bool found = false;
      for (auto& s : d.sections)
        if (s.name == name) {
          for (std::size_t r = 0; r < s.polyline.rows(); ++r) pl.push_back({s.polyline.si(2 * r), s.polyline.si(2 * r + 1), 0});
          if (s.detail_window.present) det = s.detail_window.si_list();
          found = true;
        }
      if (!found) throw Error(ErrorCode::InvalidArgument, "no section named '" + name + "' in the description", "outputs.sections");
    }
    if (pl.size() < 2) throw Error(ErrorCode::InvalidArgument, "a section needs at least two vertices", "line");
    auto s = repr::make_section(*o.model, nm, pl, std::max(5.0, c.grid_spacing / 4), 2.0);
    auto f = fig::section_figure(c, s, det.empty() ? nullptr : &det);
    fig::write_figure(f, out_dir, fig::figure_basename(c, "section", nm), c.model_name + " - section " + nm,
                      "Section " + nm + ": units per fault block, fault traces between tips, wells within 300 m.",
                      std::string("Evidence basis: ") + basis_label(o.model->info().basis) + ". Section length " + fmt_fixed(s.length, 0) + " m.");
  }
  Manifest man = Manifest::from_description(d, "section");
  man.add_outputs(out_dir);
  man.write(out_dir / "manifest.json");
  return o;
}

BuildOutcome run_export(const std::filesystem::path& description, const std::filesystem::path& out_dir, const std::vector<std::string>& formats,
                        double grid_spacing, const BuildOptions& opt) {
  auto d = desc::load_description(description);
  ensure_new_output_dir(out_dir);
  BuildOutcome o = build_only(d, opt);
  CheckContext ctx(*o.model, d, make_well_paths(d), d.grid_spacing.si());
  CheckReport rep = run_checks(ctx);
  o.checks_passed = rep.passed();
  write_file(out_dir / "checks_report.txt", rep.text());
  if (!o.checks_passed) {
    log(LogLevel::Error, "checks failed: nothing exported (see checks_report.txt)");
  } else {
    auto ids = resolve_export_ids(d, formats);
    if (ids.empty()) throw Error(ErrorCode::InvalidArgument, "no available export format requested");
    write_exports(d, *o.model, ids, out_dir / "exports", grid_spacing, ctx.grid_products(), ctx.surfaces());
  }
  Manifest man = Manifest::from_description(d, "export");
  man.add_outputs(out_dir);
  man.write(out_dir / "manifest.json");
  return o;
}

BuildOutcome run_state(const std::filesystem::path& description, const std::filesystem::path& out_dir, double age,
                       const BuildOptions& opt) {
  auto d = desc::load_description(description);
  ensure_new_output_dir(out_dir);
  BuildOutcome full = build_only(d, opt);
  auto sd = state_description(d, *full.model, age);
  BuildOutcome o;
  o.model = build_model(sd, BuildInputs{}, o.log, "parametric");
  write_file(out_dir / "model_summary.txt", summary_text(sd, *o.model, o.log));
  CheckContext ctx(*o.model, sd, make_well_paths(sd), sd.grid_spacing.si());
  ctx.timeline_checks = false;
  CheckReport rep = run_checks(ctx);
  o.checks_passed = rep.passed();
  write_file(out_dir / "checks_report.txt", rep.text());
  write_file(out_dir / "checks_report.csv", rep.csv());
  auto keys = key_numbers(sd, *o.model, rep, ctx);
  keys.insert(keys.begin(), {"State", fmt_shortest(age) + " Ma"});
  write_file(out_dir / "summary_table.csv", summary_csv(keys));
  if (o.checks_passed) {
    if (opt.figures) write_figures(sd, *o.model, keys, out_dir / "figures");
    auto ids = resolve_export_ids(sd, opt.exports);
    if (!ids.empty()) write_exports(sd, *o.model, ids, out_dir / "exports", 0);
  }
  Manifest man = Manifest::from_description(d, "state --age " + fmt_shortest(age));
  for (auto& [k, v] : o.log.values) man.extra.emplace_back(k, v);
  man.add_outputs(out_dir);
  man.write(out_dir / "manifest.json");
  return o;
}

// Stage 2F (STX-TN-10): implemented in Prompts 32 and 33.
}  // namespace stratonix

// ---- speed and scale (stage 4D, STX-TN-17 section 4) ----------------------------------------------
namespace stratonix {
/// Whole builds memoised by the hash of the resolved description and the options (a build of the
/// same description is copied, not repeated), and the input cache switched on while it lives.
struct BuildCache {
  std::mutex mx;
  std::map<std::string, std::filesystem::path> done;  ///< key -> folder of the first build
  bool previous_input_cache = false;
  ~BuildCache() { io::set_input_cache(previous_input_cache); }
};
std::shared_ptr<BuildCache> make_build_cache() {
  auto c = std::make_shared<BuildCache>();
  c->previous_input_cache = io::set_input_cache(true);
  return c;
}
BuildOutcome run_build_cached(const std::filesystem::path& description, const std::filesystem::path& out_dir, const BuildOptions& opt, BuildCache& cache) {
  const auto d = desc::load_description(description);
  std::string key = d.file_sha256 + "|" + sha256_hex(desc::to_json(d.resolved)) + "|" + (opt.figures ? "f" : "-") + "|" + opt.builder + "|" +
                    (opt.reference_kernels ? "r" : "o");
  for (auto& e : opt.exports) key += "|" + e;
  std::filesystem::path from;
  {
    std::lock_guard<std::mutex> lk(cache.mx);
    auto it = cache.done.find(key);
    if (it != cache.done.end()) from = it->second;
  }
  if (from.empty()) {
    auto out = run_build(description, out_dir, opt);
    std::lock_guard<std::mutex> lk(cache.mx);
    cache.done.emplace(key, out_dir);
    return out;
  }
  // served from the cache: the first build's files, the manifest stating the copy
  if (std::filesystem::exists(out_dir) && !std::filesystem::is_empty(out_dir))
    throw Error(ErrorCode::Io, "output folder " + out_dir.string() + " is not empty (an issued output is never overwritten)");
  std::filesystem::create_directories(out_dir);
  BuildOutcome res;
  for (auto& e : std::filesystem::recursive_directory_iterator(from)) {
    const auto rel = std::filesystem::relative(e.path(), from);
    if (e.is_directory()) {
      std::filesystem::create_directories(out_dir / rel);
      continue;
    }
    std::filesystem::copy_file(e.path(), out_dir / rel);
    res.written.push_back(rel.generic_string());
  }
  res.checks_passed = true;
  res.model = build_only(d, opt).model;
  return res;
}

std::size_t stream_grid(const std::filesystem::path& description, const std::filesystem::path& out_dir, const GridStreamOptions& o) {
  // the grid is built in memory (compact arrays) and its files are written as they are made
  // (GRDECL never held whole); a memory limit refuses a grid whose arrays would not fit
  const auto d = desc::load_description(description);
  if (!d.grid.present) throw Error(ErrorCode::InvalidArgument, "stream_grid: the description has no grid block");
  if (o.memory_limit) {
    const auto& dom = d.size;
    const double sp0 = d.grid.spacing.si(0), sp1 = d.grid.spacing.values.size() > 1 ? d.grid.spacing.si(1) : sp0;
    double layers = 0;
    for (auto& z : d.grid.zones) layers += z.layers.present ? z.layers.si() : 1;
    const double cells = std::ceil(dom.si(0) / sp0) * std::ceil(dom.si(1) / sp1) * std::max(1.0, layers);
    const double need = cells * 130.0;  // corners 64 B, codes 16 B, reason 32 B, ZONE and BLOCK 16 B
    if (need > static_cast<double>(o.memory_limit))
      throw Error(ErrorCode::InvalidArgument, "stream_grid: about " + fmt_fixed(cells / 1e6, 1) + " million cells need about " + fmt_fixed(need / 1e9, 2) +
                                                  " GB, above the stated limit of " + fmt_fixed(static_cast<double>(o.memory_limit) / 1e9, 2) + " GB");
  }
  BuildOptions bo;
  bo.figures = false;
  bo.exports = o.formats.empty() ? std::vector<std::string>{"grdecl"} : o.formats;
  run_build(description, out_dir, bo, "grid");
  return peak_rss_bytes();
}

std::vector<std::pair<std::string, std::string>> parallel_build_digest(const std::filesystem::path& description, unsigned threads) {
  // a build at the given thread count into a scratch folder; every file but the manifest (which
  // records the thread count and the time) by its SHA-256, in path order
  const auto tmp = std::filesystem::temp_directory_path() / ("stratonix_digest_" + std::to_string(threads) + "_" +
                                                              sha256_hex(std::filesystem::absolute(description).generic_string()).substr(0, 12));
  std::filesystem::remove_all(tmp);
  BuildOptions bo;
  bo.threads = threads;
  bo.figures = false;
  run_build(description, tmp, bo);
  std::vector<std::pair<std::string, std::string>> out;
  for (auto& e : std::filesystem::recursive_directory_iterator(tmp))
    if (e.is_regular_file() && e.path().filename() != "manifest.json") out.push_back({std::filesystem::relative(e.path(), tmp).generic_string(), sha256_file(e.path())});
  std::sort(out.begin(), out.end());
  std::filesystem::remove_all(tmp);
  return out;
}
}  // namespace stratonix
