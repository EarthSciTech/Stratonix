#include "cli_commands.hpp"
#include "stratonix/package/uncertainty.hpp"

#include <fstream>

#include <cctype>
#include <cmath>
#include <iostream>

#include "stratonix/checks/manifest.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/foldframe.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/model/model_suite.hpp"
#include "stratonix/package/bench.hpp"
#include "stratonix/package/delivery.hpp"
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

int cmd_state(Args& a) {
  if (a.positional().empty() || !a.has("o") || !a.has("age"))
    throw Error(ErrorCode::InvalidArgument, "usage: stratonix state <description.yaml> --age <Ma> -o <new folder>");
  BuildOptions o;
  o.figures = !a.flag("no-figures");
  o.exports = csv(a.get("formats"));
  const double age = std::stod(a.get("age"));
  auto r = run_state(a.positional()[0], a.get("o"), age, o);
  std::cout << "state at " << a.get("age") << " Ma of " << r.model->info().name << " -> " << a.get("o") << "\n";
  if (!r.checks_passed) {
    std::cout << "checks FAILED (see checks_report.txt)\n";
    return 4;
  }
  return 0;
}

int cmd_section(Args& a) {
  if (a.positional().empty() || !a.has("o")) throw Error(ErrorCode::InvalidArgument, "usage: stratonix section <description.yaml> -o <new folder> (--name N | --line x1,y1;x2,y2;...)");
  Polyline pl;
  if (a.has("line"))
    for (auto& v : split(a.get("line"), ';')) {
      auto xy = csv(v);
      if (xy.size() != 2) throw Error(ErrorCode::InvalidArgument, "--line expects x,y pairs separated by ';'");
      pl.push_back({parse_double(xy[0]), parse_double(xy[1]), 0});
    }
  if (pl.empty() && !a.has("name")) throw Error(ErrorCode::InvalidArgument, "give --name or --line");
  BuildOptions o;
  o.reference_kernels = a.flag("reference-kernels");
  auto r = run_section(a.positional()[0], a.get("o"), a.get("name", pl.empty() ? "" : "line"), pl, o);
  std::cout << (r.checks_passed ? "section written to " : "checks FAILED; no section drawn: ") << a.get("o") << "\n";
  return r.checks_passed ? 0 : 2;
}
int cmd_export(Args& a) {
  if (a.positional().empty() || !a.has("o"))
    throw Error(ErrorCode::InvalidArgument, "usage: stratonix export <description.yaml> -o <new folder> [--formats petrel|id,id,...] [--grid-spacing m]");
  BuildOptions o;
  o.reference_kernels = a.flag("reference-kernels");
  double sp = a.has("grid-spacing") ? parse_double(a.get("grid-spacing")) : 0.0;
  auto r = run_export(a.positional()[0], a.get("o"), csv(a.get("formats")), sp, o);
  std::cout << (r.checks_passed ? "exported to " : "checks FAILED; nothing exported: ") << a.get("o") << "\n";
  return r.checks_passed ? 0 : 2;
}
int cmd_package(Args& a) {
  if (a.positional().empty() || !a.has("o"))
    throw Error(ErrorCode::InvalidArgument, "usage: stratonix package <description.yaml> -o <new folder> [--config package.yaml] [--zip] [--formats ...]");
  PackageConfig cfg = a.has("config") ? PackageConfig::load(a.get("config")) : PackageConfig{};
  if (a.flag("zip")) cfg.zip = true;
  if (a.flag("no-figures")) cfg.figures = false;
  if (a.has("formats")) cfg.formats = csv(a.get("formats"));
  BuildOptions o;
  o.reference_kernels = a.flag("reference-kernels");
  auto d = desc::load_description(a.positional()[0]);
  std::cout << "inputs and their sharing status:\n";
  if (d.inputs.empty()) std::cout << "  (none recorded in the description)\n";
  for (auto& in : d.inputs) std::cout << "  " << in.id << ": " << in.sharing << "  (" << in.path << ")\n";
  auto r = run_package(a.positional()[0], a.get("o"), cfg, o);
  if (!r.scan_passed) {
    std::cout << "identifier scan FAILED (" << r.hits.size() << " hits); package rejected: " << r.folder.string() << "\n";
    for (std::size_t i = 0; i < r.hits.size() && i < 10; ++i)
      std::cout << "  " << r.hits[i].file << ":" << r.hits[i].line << "  " << r.hits[i].pattern << "\n";
    return 6;
  }
  std::cout << "package written: " << r.folder.string() << (r.zip.empty() ? "" : " and " + r.zip.string()) << "\n";
  return 0;
}

int cmd_ensemble(Args& a) {
  if (a.positional().empty() || !a.has("o"))
    throw Error(ErrorCode::InvalidArgument, "usage: stratonix ensemble <description.yaml> -o <new folder> [--members N] [--seed S] [--threads N] [--no-figures]");
  ensure_new_output_dir(a.get("o"));
  unc::EnsembleOptions o;
  if (a.has("members")) o.members = std::stoi(a.get("members"));
  if (a.has("seed")) o.seed = std::stoll(a.get("seed"));
  if (a.has("threads")) o.threads = static_cast<unsigned>(std::stoul(a.get("threads")));
  o.figures = !a.flag("no-figures");
  o.base_build = true;
  auto r = unc::run_ensemble(a.positional()[0], a.get("o"), o);
  std::size_t failed = 0;
  for (auto& m : r.members) failed += m.sample.index > 0 && m.status != "ok";
  std::cout << "ensemble of " << r.spec.members << " members (seed " << r.spec.seed << "): " << r.ok() << " ok, " << failed << " failed -> " << a.get("o") << "\n";
  for (auto& k : r.stats)
    std::cout << "  " << k.key << ": P10 " << fmt_fixed(k.stats.p10, 2) << ", P50 " << fmt_fixed(k.stats.p50, 2) << ", P90 " << fmt_fixed(k.stats.p90, 2) << " (n "
              << k.stats.n << ")\n";
  return 0;
}
int cmd_fold_frame(Args& a) {
  if (a.positional().empty() || !a.has("o"))
    throw Error(ErrorCode::InvalidArgument, "usage: stratonix fold-frame <bedding.csv> [--foliation <axial_foliation.csv>] [--at <points.csv>] -o <frame.json>");
  auto& reg = registry::Registry::global();
  auto read_orient = [&](const std::string& f) {
    auto d = reg.create<io::Importer>(registry::Kind::Importer, "csv_orientations")->read(f);
    std::vector<OrientationDatum> v;
    for (std::size_t i = 0; i < d.points.size(); ++i) v.push_back({d.points[i], d.normals[i]});
    return v;
  };
  auto bed = read_orient(a.positional()[0]);
  std::vector<OrientationDatum> fol;
  if (a.has("foliation")) fol = read_orient(a.get("foliation"));
  FoldFrameOptions o;
  auto F = fit_fold_frame(bed, fol, o);
  std::string j = "{\n  \"ok\": " + std::string(F.ok ? "true" : "false") + ",\n  \"reason\": \"" + F.reason + "\",\n";
  auto num = [](double v) { return std::isfinite(v) ? fmt_shortest(v) : std::string("null"); };
  j += "  \"axis\": [" + num(F.axis.x) + ", " + num(F.axis.y) + ", " + num(F.axis.z) + "],\n";
  j += "  \"trend_deg\": " + num(rad2deg(F.trend)) + ",\n  \"plunge_deg\": " + num(rad2deg(F.plunge)) + ",\n";
  j += "  \"interlimb_deg\": " + num(rad2deg(F.interlimb)) + ",\n  \"wavelength_m\": " + num(F.wavelength) + ",\n";
  j += "  \"order\": " + std::to_string(F.order) + ",\n  \"cylindricity\": " + num(F.cylindricity) + ",\n";
  j += "  \"rms_deg\": " + num(rad2deg(F.rms)) + ",\n  \"data\": " + std::to_string(F.data);
  if (a.has("at") && F.ok) {
    auto d = reg.create<io::Importer>(registry::Kind::Importer, "csv_points")->read(a.get("at"));
    j += ",\n  \"predicted_normals\": [";
    for (std::size_t i = 0; i < d.points.size(); ++i) {
      Vec3 n = F.bedding_normal(d.points[i]);
      j += std::string(i ? ", " : "") + "[" + num(n.x) + ", " + num(n.y) + ", " + num(n.z) + "]";
    }
    j += "]";
  }
  j += "\n}\n";
  std::ofstream(a.get("o"), std::ios::binary) << j;
  std::cout << "fold frame " << (F.ok ? "fitted" : "NOT fitted: " + F.reason) << ": axis trend " << fmt_fixed(rad2deg(F.trend), 1) << " deg, plunge "
            << fmt_fixed(rad2deg(F.plunge), 1) << " deg, interlimb " << fmt_fixed(rad2deg(F.interlimb), 1) << " deg, wavelength "
            << fmt_fixed(F.wavelength, 0) << " m -> " << a.get("o") << "\n";
  return F.ok ? 0 : 4;
}
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
