/// @file upscale.hpp
/// @brief Stage 4C (STX-TN-16 section 1.5): logs upscaled to the cells a well passes through.
#pragma once
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/repr/grid.hpp"

namespace stratonix::repr {

enum class UpscaleMethod { Arithmetic, Harmonic, Geometric, Rms, Minimum, Maximum, MostFrequent };
std::string upscale_method_name(UpscaleMethod m);

/// A log sampled along a path: positions (x, y, z), measured depths and values (NaN at nulls).
struct PathLog { std::string well, curve; std::vector<Vec3> p; std::vector<double> md, v; };
struct UpscaledCells {
  std::string method;
  std::vector<std::size_t> cells;          ///< cell indices filled
  std::vector<double> value;               ///< per filled cell
  std::vector<double> md_top, md_base;     ///< interval of the well in each cell
  std::vector<std::size_t> missing;        ///< cells crossed but left missing (too few valid samples)
  std::vector<std::string> notes;          ///< ties of the most frequent value, etc.
};
UpscaledCells upscale_log(const CornerPointGrid& g, const PathLog& log, UpscaleMethod method, double min_valid_fraction = 0.5);

}  // namespace stratonix::repr
