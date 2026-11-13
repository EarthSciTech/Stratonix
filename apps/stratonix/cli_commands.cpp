#include "cli_commands.hpp"

#include <cctype>
#include <iostream>

#include "stratonix/checks/manifest.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/model/model_suite.hpp"
#include "stratonix/package/bench.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/registry/registry.hpp"

namespace stratonix::cli {

int cmd_validate(Args& a) {
  auto file = a.pos(0, "description file");
  desc::Description d;
  auto rep = desc::validate_description(file, &d);
  std::cout << rep.text();
  if (rep.ok()) {
    std::cout << "valid: " << file << " (model " << d.model.name << " " << d.model.version << ", "
              << d.model.status << ", " << d.horizons.size() << " horizons, " << d.faults.size() << " faults, "
              << d.all_params().size() << " labelled parameters)\n";
    return 0;
  }
  std::size_t n = 0;
  for (auto& i : rep.issues) n += i.error ? 1 : 0;
  std::cout << "invalid: " << file << " (" << n << " error(s))\n";
  return 3;
}

static std::vector<std::string> csv(const std::string& s) {
  std::vector<std::string> r;
  for (auto& x : split(s, ','))
    if (!trim(x).empty()) r.push_back(trim(x));
  return r;
}

int cmd_build(Args& a) {
  BuildOptions o;
  o.reference_kernels = a.flag("reference-kernels");
  o.figures = !a.flag("no-figures");
  o.exports = csv(a.get("formats"));
  o.builder = a.get("builder");
  auto r = run_build(a.pos(0, "description file"), a.out_dir(), o);
  std::cout << "built " << r.model->info().name << " -> " << a.out_dir() << "\n";
  for (auto& n : r.log.notes) std::cout << "  " << n << "\n";
  if (!r.checks_passed) {
    std::cout << "checks FAILED (see checks_report.txt)\n";
    return 4;
  }
  return 0;
}

int cmd_section(Args&) { throw Error(ErrorCode::NotImplemented, "sections arrive with the figures stage"); }
int cmd_export(Args&) { throw Error(ErrorCode::NotImplemented, "exporters arrive with the export stages"); }
int cmd_package(Args&) { throw Error(ErrorCode::NotImplemented, "packaging arrives with the packaging stage"); }

int cmd_plugins(Args&) {
  for (auto* e : registry::Registry::global().all())
    std::cout << registry::kind_name(e->kind) << "  " << e->name << "  " << e->version << "  " << e->description << "\n";
  return 0;
}

int cmd_suite(Args& a) {
  auto d = desc::load_description(a.pos(0, "description file"));
  BuildOptions o;
  auto r = build_only(d, o);
  auto rep = run_model_suite(*r.model);
  std::cout << rep.text();
  return rep.passed() ? 0 : 4;
}

int cmd_verify(Args& a) {
  auto probs = verify_manifest(a.pos(0, "manifest.json"));
  for (auto& p : probs) std::cout << "problem: " << p << "\n";
  std::cout << (probs.empty() ? "manifest verified\n" : "manifest NOT verified\n");
  return probs.empty() ? 0 : 4;
}

int cmd_schema(Args&) {
  std::cout << desc::schema_v1_text();
  return 0;
}

int cmd_version(Args&) {
  auto& b = build_info();
  std::cout << "stratonix " << b.version << "\ncommit " << b.git_commit << "\ncompiler " << b.compiler
            << "\nplatform " << b.platform << "\nbuild " << b.build_type << "\nschema " << desc::kSchemaVersion
            << "\nreprojection " << (crs_reprojection_available() ? "PROJ" : "not available (built without PROJ)") << "\n";
  return 0;
}

}  // namespace stratonix::cli

namespace stratonix::cli {
int cmd_bench(Args& a) {
  if (a.has("write-model")) {
    write_file(a.get("write-model"), production_model_yaml(std::stod(a.get("spacing", "50"))));
    std::cout << "wrote " << a.get("write-model") << "\n";
    return 0;
  }
  BenchOptions o;
  o.spacings.clear();
  for (auto& x : csv(a.get("spacing", "50"))) o.spacings.push_back(parse_double(x));
  std::string mode = a.get("mode", "both");
  o.reference = mode != "optimised";
  o.optimised = mode != "reference";
  if (a.has("voxel")) {
    auto v = csv(a.get("voxel"));
    o.voxel_xy = parse_double(v.at(0));
    o.voxel_z = parse_double(v.at(1));
  }
  if (!a.positional().empty()) o.description = a.positional()[0];
  auto rows = run_benchmark(o);
  std::string t = bench_table(rows);
  std::cout << t;
  if (a.has("o")) {
    std::filesystem::create_directories(a.get("o"));
    write_file(std::filesystem::path(a.get("o")) / "bench_report.txt", t);
    write_file(std::filesystem::path(a.get("o")) / "bench_report.csv", bench_csv(rows));
  }
  return 0;
}
}  // namespace stratonix::cli
