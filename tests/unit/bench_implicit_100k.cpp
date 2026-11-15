// Stage 3A benchmark (STX-TN-11 section 1): 100 000 interface points and 10 000 orientations,
// no faults, surfaces at 100 m over 8 km x 6 km, build within 300 s on the CI Linux release
// machine. Debug and sanitiser builds skip it (exit code 77).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>

#include "../stx_test.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/package/pipeline.hpp"

namespace fs = std::filesystem;
using namespace stratonix;

namespace {
double surf(int h, double x, double y) {
  double r2 = (x - 4000) * (x - 4000) + (y - 3000) * (y - 3000);
  return -(900.0 + 180.0 * h) + 160.0 * std::exp(-r2 / (2 * 1800.0 * 1800.0)) + 0.01 * (x - 4000);
}
}  // namespace

int main() {
#if !defined(NDEBUG) || defined(STX_SANITIZE)
  std::printf("SKIP bench_implicit_100k: optimised (release) builds only\n");
  return 77;
#else
  register_all_plugins();
  fs::path dir = fs::path(stxtest::tmp_dir()) / "bench_implicit_100k";
  fs::remove_all(dir);
  fs::create_directories(dir);
  std::mt19937_64 rng(100000);
  std::uniform_real_distribution<double> UX(100, 7900), UY(100, 5900);
  {
    std::ofstream p(dir / "points.csv"), o(dir / "orient.csv");
    p << "X,Y,Z,surface\n";
    o << "X,Y,Z,G_x,G_y,G_z,surface\n";
    for (int i = 0; i < 100000; ++i) {
      int h = i % 5;
      double x = UX(rng), y = UY(rng), z = surf(h, x, y);
      p << fmt_fixed(x, 2) << "," << fmt_fixed(y, 2) << "," << fmt_fixed(z, 3) << ",H" << h + 1 << "\n";
      if (i % 10 == 0) {
        const double e = 1;
        double zx = (surf(h, x + e, y) - surf(h, x - e, y)) / (2 * e), zy = (surf(h, x, y + e) - surf(h, x, y - e)) / (2 * e);
        double n = std::sqrt(zx * zx + zy * zy + 1);
        o << fmt_fixed(x, 2) << "," << fmt_fixed(y, 2) << "," << fmt_fixed(z, 3) << "," << fmt_fixed(-zx / n, 6) << "," << fmt_fixed(-zy / n, 6) << ","
          << fmt_fixed(1 / n, 6) << ",H" << h + 1 << "\n";
      }
    }
  }
  const char* u = " | user_instruction | benchmark";
  auto q = [&](const std::string& v) { return "\"" + v + u + "\""; };
  std::string y = "stratonix_schema: 1\nmodel: {name: bench, version: \"1.0.0\", status: DESIGN}\ncrs: {name: local, vertical_datum: datum, length_unit: m}\n"
                  "domain:\n  origin: " + q("[0, 0] m") + "\n  size: " + q("[8000, 6000] m") + "\n  depth_range: " + q("[0, 2500] m") +
                  "\n  grid_spacing: " + q("50 m") + "\n"
                  "inputs:\n  - {id: P, path: points.csv, role: surface points, sharing: shareable, source: benchmark, evidence: user_instruction}\n"
                  "  - {id: O, path: orient.csv, role: orientations, sharing: shareable, source: benchmark, evidence: user_instruction}\n"
                  "wells:\n  - {name: W1, role: check, head: " + q("[4000, 3000] m") + "}\n"
                  "stratigraphy:\n  reference: {well: W1}\n  horizons: [{name: H1}, {name: H2}, {name: H3}, {name: H4}, {name: H5}]\n"
                  "  units: [{name: U1, top: H1, base: H2}, {name: U2, top: H2, base: H3}, {name: U3, top: H3, base: H4}, {name: U4, top: H4, base: H5}]\n"
                  "implicit:\n  data:\n    - {input: P, kind: points, format: csv_points}\n    - {input: O, kind: orientations, format: csv_orientations}\n"
                  "  use_well_tops: false\n  surface_spacing: " + q("100 m") + "\n";
  desc::Description d;
  auto rep = desc::validate_description_text(y, dir, &d);
  if (!rep.ok()) {
    std::printf("%s", rep.text().c_str());
    return 1;
  }
  auto t0 = std::chrono::steady_clock::now();
  BuildInputs in;
  io::load_build_inputs(d, in);
  BuildLog log;
  auto m = build_model(d, in, log);
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::string solver;
  for (auto& [k, v] : log.values)
    if (k.rfind("implicit.", 0) == 0 && k.find("iteration") == std::string::npos) std::printf("  %s = %s\n", k.c_str(), v.c_str());
  std::printf("  100 000 points, 10 000 orientations: built in %.1f s on %u threads (limit 300 s)\n", secs, thread_count());
  if (!(secs <= 300.0)) {
    std::printf("FAIL bench_implicit_100k: %.1f s exceeds 300 s\n", secs);
    return 1;
  }
  std::printf("PASS bench_implicit_100k\n");
  return 0;
#endif
}
