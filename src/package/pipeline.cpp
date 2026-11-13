// Pipeline shared by the command line and the C API.
#include "stratonix/package/pipeline.hpp"

#include <sstream>

#include "stratonix/checks/checks.hpp"
#include "stratonix/checks/manifest.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/figures/figures.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/io/io.hpp"

namespace stratonix {

void register_all_plugins() {
  static bool done = false;
  if (done) return;
  done = true;
  register_engine_plugins();
  io::register_io_plugins();
  register_check_plugins();
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
    std::string role, color;
    for (auto& du : d.units)
      if (du.name == u.name) { role = du.role; color = du.color; }
    c.unit_roles.push_back(role);
    c.unit_colors.push_back(color);
  }
  return c;
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
}

BuildOutcome run_build(const std::filesystem::path& description, const std::filesystem::path& out_dir,
                       const BuildOptions& opt, const std::string& command) {
  auto d = desc::load_description(description);
  ensure_new_output_dir(out_dir);
  BuildOutcome o = build_only(d, opt);
  write_file(out_dir / "model_summary.txt", summary_text(d, *o.model, o.log));
  CheckContext ctx(*o.model, d, make_well_paths(d), d.grid_spacing.si());
  CheckReport rep = run_checks(ctx);
  o.checks_passed = rep.passed();
  write_file(out_dir / "checks_report.txt", rep.text());
  write_file(out_dir / "checks_report.csv", rep.csv());
  auto keys = key_numbers(d, *o.model, rep, ctx);
  write_file(out_dir / "summary_table.csv", summary_csv(keys));
  if (!o.checks_passed) {
    log(LogLevel::Error, "checks failed: exports and figures are not written (see checks_report.txt)");
  } else {
    if (opt.figures) write_figures(d, *o.model, keys, out_dir / "figures");
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

}  // namespace stratonix
