// Shared kernels of the representation layer.
#pragma once
#include <vector>

#include "stratonix/model/model.hpp"

namespace stratonix::repr {

/// Roots of horizon h's scalar along verticals at (xs, ys) in block `block`, searched in
/// [z_lo, z_hi]. NaN where the scalar does not change sign in the range.
/// `guess` (optional, NaN = none) narrows the initial bracket to guess +/- 40 m when the scalar
/// changes sign there; otherwise the full bracket is used.
void vertical_roots(const Model& m, int h, BlockCode block, const std::vector<double>& xs, const std::vector<double>& ys,
                    double z_hi, double z_lo, std::vector<double>& z_out, const std::vector<double>* guess = nullptr);

/// Faulted horizon elevations at points (NaN in gaps or outside the domain's elevation range).
void faulted_roots(const Model& m, int h, const std::vector<double>& xs, const std::vector<double>& ys,
                   std::vector<double>& z, std::vector<BlockCode>& blk, std::size_t* gaps, std::size_t* overlaps,
                   const std::vector<double>* guess = nullptr);

/// Elevation where fault f's distance is zero along verticals (NaN if none in [z_lo, z_hi]).
void fault_vertical_roots(const Model& m, int f, const std::vector<double>& xs, const std::vector<double>& ys,
                          double z_hi, double z_lo, std::vector<double>& z_out, std::vector<double>* extent);

}  // namespace stratonix::repr
