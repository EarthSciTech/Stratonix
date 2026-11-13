// Regression against the prototype's reference tables (skipped until they are supplied).
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "../stx_test.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

using namespace stratonix;
namespace fs = std::filesystem;

int main() {
  const char* dir = std::getenv("STRATONIX_PROTOTYPE_REF");
  if (!dir || !fs::exists(fs::path(dir) / "prototype_description.yaml")) {
    std::cout << "SKIP: prototype reference files not supplied (set STRATONIX_PROTOTYPE_REF; see tests/regression/README.md)\n";
    return 77;
  }
  register_all_plugins();
  auto d = desc::load_description(fs::path(dir) / "prototype_description.yaml");
  BuildOptions o;
  auto m = build_only(d, o).model;
  std::ifstream in(fs::path(dir) / "horizons.csv");
  std::string line;
  std::getline(in, line);
  double worst = 0;
  std::size_t n = 0;
  while (std::getline(in, line)) {
    auto t = split(line, ',');
    if (t.size() < 4) continue;
    int h = d.horizon_index(trim(t[0]));
    double x = parse_double(t[1]), y = parse_double(t[2]), z = parse_double(t[3]);
    auto g = repr::horizon_grid_faulted(*m, h, GridSpec{x, y, 1, 1, 1, 1});
    if (std::isnan(g.z.v[0])) continue;
    worst = std::max(worst, std::fabs(g.z.v[0] - z));
    ++n;
  }
  std::cout << n << " horizon samples, max difference " << fmt_fixed(worst, 4) << " m\n";
  return (n > 0 && worst <= 0.01) ? 0 : 1;
}
