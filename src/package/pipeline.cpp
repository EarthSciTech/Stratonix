// Pipeline shared by the command line and the C API.
#include "stratonix/package/pipeline.hpp"

#include <algorithm>
#include <sstream>

#include "stratonix/checks/checks.hpp"
#include "stratonix/checks/manifest.hpp"
#include "stratonix/engines/conditioning.hpp"
#include "stratonix/engines/geophysics.hpp"
#include "stratonix/engines/posterior.hpp"
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
  for (std::size_t f = 0; f < inf.faults.size(); ++f) {
    double best = 0;
    std::string hn;
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      double t = std::fabs(ctx.throw_row(static_cast<int>(f), static_cast<int>(h)).throw_max);
      if (t > best) { best = t; hn = inf.horizons[h].name; }
    }
    k.push_back({"Max throw " + inf.faults[f].name, hn.empty() ? "0 m" : fmt_fixed(best, 1) + " m at " + hn});
  }
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
    c.unit_roles.push_back(role);
    c.unit_colors.push_back(color);
    c.unit_labels.push_back(label);
  }
  for (auto& n : d.isopachs)
    for (std::size_t u = 0; u < m.info().units.size(); ++u)
      if (m.info().units[u].name == n) c.isopach_units.push_back(static_cast<int>(u));
  // Burial history at every well head (stage 2D), when the description has compaction.
  if (auto* st = static_cast<const Structure*>(m.capability("structure")); st && st->compaction_mode() != 0) {
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
      c.burial.push_back(std::move(b));
    }
  }
  c.fault_rock = d.fault_rock.present;
  if (d.fault_rock.present || d.stress.present) {
    c.rock = fault_rock_input(d, m);
    c.stress = stress_input(d);
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
                   const std::filesystem::path& folder) {
  auto c = figure_context(d, m, keys);
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
                                            const std::filesystem::path& folder, double grid_spacing, const GridProducts* gp = nullptr) {
  double sp = grid_spacing > 0 ? grid_spacing : d.grid_spacing.si();
  double sticks = d.fault_stick_spacing.present ? d.fault_stick_spacing.si() : 250.0;  // Stratonix default when not given
  io::ExportContext c(m, make_well_paths(d), GridSpec::covering(m.info().domain, sp), std::max(10.0, sp / 2), sticks);
  c.vertical_exaggeration = d.vertical_exaggeration.present ? d.vertical_exaggeration.si() : 1.0;
  c.fault_rock = d.fault_rock.present;
  if (d.fault_rock.present || d.stress.present) {
    c.rock = fault_rock_input(d, m);
    c.stress = stress_input(d);
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
  const auto sin = stress_input(d);
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
    if (opt.figures) write_figures(d, *o.model, keys, out_dir / "figures");
    if (opt.figures && !d.events.empty()) write_restoration(d, *o.model, out_dir / "restoration");
    if (opt.figures) write_data_figures(d, *o.model, keys, out_dir / "figures");
    if (opt.figures && !pf.stations.empty()) write_geophysics_figures(d, *o.model, keys, pf, out_dir / "figures");
    if (gp) write_grid_outputs(d, *o.model, *gp, keys, out_dir, opt.figures);
    auto ids = resolve_export_ids(d, opt.exports);
    if (!ids.empty()) write_exports(d, *o.model, ids, out_dir / "exports", 0, gp);
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
    write_exports(d, *o.model, ids, out_dir / "exports", grid_spacing, ctx.grid_products());
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
