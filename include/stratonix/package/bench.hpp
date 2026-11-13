/// @file bench.hpp
/// @brief Timing and memory benchmarks of a full build and of each kernel, with reference and
/// optimised kernels, at chosen sizes (used by `stratonix bench` and the benchmark executable).
#pragma once
#include <filesystem>
#include <string>
#include <vector>

namespace stratonix {

/// Synthetic production-size description: 40 x 40 km, six horizons, two folds, ten faults.
/// Every parameter is labelled user_instruction, source "synthetic benchmark".
std::string production_model_yaml(double grid_spacing_m);

struct BenchRow {
  std::string size, mode, step;
  double seconds = 0;
  double peak_rss_mb = 0;
  std::string detail;
};

struct BenchOptions {
  std::vector<double> spacings = {50};    ///< grid spacings (m) of the production model
  bool reference = true;                  ///< also run the reference kernels
  bool optimised = true;
  unsigned threads = 0;
  double voxel_xy = 200, voxel_z = 25;    ///< voxel grid resolution (m)
  std::filesystem::path description;      ///< optional: benchmark this description instead
};

/// Run the benchmark and return one row per step and mode.
std::vector<BenchRow> run_benchmark(const BenchOptions& o);
std::string bench_table(const std::vector<BenchRow>& rows);
std::string bench_csv(const std::vector<BenchRow>& rows);

}  // namespace stratonix
