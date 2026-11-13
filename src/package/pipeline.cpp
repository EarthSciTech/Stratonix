// Pipeline shared by the command line and the C API.
#include "stratonix/package/pipeline.hpp"

#include <sstream>

#include "stratonix/checks/manifest.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"

namespace stratonix {

void register_all_plugins() {
  static bool done = false;
  if (done) return;
  done = true;
  register_engine_plugins();
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
  o.model = build_model(d, o.log, opt.builder);
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

BuildOutcome run_build(const std::filesystem::path& description, const std::filesystem::path& out_dir,
                       const BuildOptions& opt, const std::string& command) {
  auto d = desc::load_description(description);
  ensure_new_output_dir(out_dir);
  BuildOutcome o = build_only(d, opt);
  write_file(out_dir / "model_summary.txt", summary_text(d, *o.model, o.log));
  Manifest man = Manifest::from_description(d, command);
  for (auto& [k, v] : o.log.values) man.extra.emplace_back(k, v);
  man.add_outputs(out_dir);
  man.write(out_dir / "manifest.json");
  return o;
}

}  // namespace stratonix
