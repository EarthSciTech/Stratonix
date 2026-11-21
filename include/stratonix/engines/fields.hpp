/// @file fields.hpp
/// @brief Stage 4A (STX-TN-14 section 1): discrete implicit field on a regular grid with free
/// horizon values, solved by preconditioned conjugate-gradient least squares, coarse to fine.
#pragma once
#include <cstddef>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/engines/implicit.hpp"

namespace stratonix {

/// An interface point of horizon `horizon` (index in the series' own horizon list).
struct HorizonDatum {
  Vec3 p;
  int horizon = -1;
  double sigma = -1;  ///< stated uncertainty (m); < 0 = the field's default weight
};

/// Options of the discrete field (`implicit.discrete`).
struct DiscreteOptions {
  double spacing = 0;           ///< m, horizontal node spacing; 0 = 1/80 of the largest horizontal extent
  double vertical_spacing = 0;  ///< m; 0 = 1/40 of the vertical extent
  double regularisation = 0.2;  ///< weight of the second-difference rows
  double orientation_weight = 10.0;  ///< weight of the gradient rows relative to interface rows
  int max_iterations = 5000;    ///< CGLS iterations per level
  double tolerance = 1e-6;      ///< relative residual of the (scaled) normal equations
  int levels = 0;               ///< coarse-to-fine levels; 0 = automatic
  int threads = 0;              ///< 0 = hardware concurrency
};

/// Field s(x) = trilinear interpolation of node values on a regular grid covering `box`, fitted by
/// least squares to interface rows s(p) - c_k = 0, orientation rows (two tangent rows and one
/// unit-magnitude row) and second-difference regularisation. The horizon values c_k are unknowns,
/// except the anchor's and those given in `fixed` (NaN = free).
class DiscreteField {
 public:
  DiscreteField(const std::vector<HorizonDatum>& data, const std::vector<RbfField::Gradient>& grads, int anchor,
                double anchor_value, const BBox& box, const DiscreteOptions& o, const std::vector<double>& fixed = {});
  double value(const Vec3& x) const;
  Vec3 gradient(const Vec3& x) const;
  /// Value of every horizon (fitted or fixed), by horizon index; NaN for a horizon without data.
  const std::vector<double>& horizon_values() const { return values_; }
  std::size_t unknowns() const { return u_.size() + free_count_; }
  std::size_t nx() const { return nx_; }
  std::size_t ny() const { return ny_; }
  std::size_t nz() const { return nz_; }
  int iterations() const { return iterations_; }       ///< total over levels
  double relative_residual() const { return residual_; }
  bool converged() const { return converged_; }
  /// Rows written to the log: grid, levels, iterations, residual, values.
  std::string summary() const;

 private:
  std::vector<double> u_;       ///< node values
  std::vector<double> values_;  ///< horizon values
  std::size_t free_count_ = 0;
  std::size_t nx_ = 0, ny_ = 0, nz_ = 0;
  Vec3 o_{0, 0, 0}, d_{1, 1, 1};
  int iterations_ = 0, levels_ = 0;
  double residual_ = 0;
  bool converged_ = false;
};

/// Rule F34: horizon values must decrease from the youngest horizon (index 0) to the oldest.
/// Returns "" when they do, else the first pair out of order, named.
std::string horizon_order_problem(const std::vector<double>& values, const std::vector<std::string>& names);

}  // namespace stratonix
