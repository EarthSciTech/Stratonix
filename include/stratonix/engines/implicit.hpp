/// @file implicit.hpp
/// @brief Data-conditioned (implicit) structure builder: a stratigraphic scalar field interpolated
/// from interface points, orientations and well tops by Hermite radial basis functions (cubic
/// kernel, linear drift), with faults from the v1 fault network. The field is interpolated in the
/// restored frame; data are restored through the faults and the field is re-solved until the
/// fault tips (which refer to the horizons) and the data agree. The result is a model behind the
/// same interface as the parametric engine, so every section, export, check and package command
/// works on it unchanged.
#pragma once
#include <memory>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/engines/parametric.hpp"

namespace stratonix {

/// Hermite RBF interpolant s(x) = sum_i a_i phi(x - x_i) - sum_j beta_j e_j . grad phi(x - y_j) + c0 + c . x
/// with phi(d) = |d|^3, fitted to values at interface points x_i and, at each orientation, to
/// two tangent constraints e . grad s = 0 (the direction of the gradient, not its magnitude).
class RbfField {
 public:
  /// sigma: the point's stated uncertainty (m), its nugget on the diagonal; < 0 = the field's
  /// nugget (stage 3A, STX-TN-11 section 5).
  struct Value { Vec3 p; double v; double sigma = -1; };
  struct Gradient { Vec3 p; Vec3 g; };
  /// Fit. nugget (m) is added to the diagonal of the value rows (0 = exact interpolation).
  /// Coordinates are centred and scaled by `scale` internally for conditioning.
  /// vertical_stretch > 1 stretches z before the kernel is applied (anisotropy: horizons vary
  /// slowly laterally compared with the spacing of the horizons).
  RbfField(std::vector<Value> values, std::vector<Gradient> grads, double nugget, Vec3 centre, double scale, double vertical_stretch = 1.0);
  double value(const Vec3& x) const;
  Vec3 gradient(const Vec3& x) const;
  std::size_t unknowns() const { return a_.size() + beta_.size() + 4; }
  const std::string& solver() const { return solver_; }

 private:
  Vec3 local(const Vec3& x) const {
    Vec3 l = (x - centre_) * (1.0 / scale_);
    l.z *= stretch_;
    return l;
  }
  std::vector<Vec3> xi_, yj_, ej_;  // local coordinates; tangent directions
  std::vector<double> a_, beta_;
  double c0_ = 0;
  Vec3 c_{0, 0, 0};
  Vec3 centre_;
  double scale_ = 1, stretch_ = 1;
  std::string solver_;
};

/// Large-data solver options (stage 3A, STX-TN-11 section 1).
struct SolverOptions {
  std::string type = "auto";           ///< auto | dense | partition_of_unity
  std::size_t points_per_patch = 300;  ///< target interface data per patch
  double overlap = 1.25;               ///< patch radius = overlap * h / sqrt(2)
  std::size_t dense_limit = 3000;      ///< auto: dense up to this many interface data
};

/// Partition of unity of Hermite RBF fields on overlapping vertical columns, blended with
/// normalised Wendland C2 weights (STX-TN-11 section 1). One patch equals RbfField.
class PuField {
 public:
  /// [x0, x1] x [y0, y1]: the area the patches cover.
  PuField(std::vector<RbfField::Value> values, std::vector<RbfField::Gradient> grads, double nugget, Vec3 centre, double scale,
          double vertical_stretch, double x0, double x1, double y0, double y1, const SolverOptions& o);
  double value(const Vec3& x) const;
  Vec3 gradient(const Vec3& x) const;
  std::size_t patches() const { return patches_.size(); }
  std::size_t largest_patch() const;  ///< unknowns of the largest patch system
  /// Leave-one-out error of every interface value (Rippa 1999 per patch, blended with the patch
  /// weights at the datum), in the order the values were given.
  std::vector<double> leave_one_out() const;

 private:
  struct Patch { double cx, cy, r; std::vector<std::size_t> vi, gi; std::shared_ptr<RbfField> f; };
  std::vector<Patch> patches_;
  std::vector<RbfField::Value> values_;
  std::vector<RbfField::Gradient> grads_;
  double nugget_ = 0, scale_ = 1, stretch_ = 1;
  Vec3 centre_;
  double h_ = 0;              ///< lattice spacing
  std::size_t nx_ = 0, ny_ = 0;
  double x0_ = 0, y0_ = 0;
};

/// Leave-one-out errors of the interface values of the dense field (Rippa 1999):
/// e_i = alpha_i / (A^-1)_ii, in the order given.
std::vector<double> rippa_loo(const std::vector<RbfField::Value>& values, const std::vector<RbfField::Gradient>& grads, double nugget,
                              Vec3 centre, double scale, double vertical_stretch);

/// Plane fitted to fault sticks: centroid, unit normal pointing up, strike and dip (radians,
/// right-hand rule: dip direction = strike + 90 degrees) and RMS distance of the points.
struct FittedPlane { Vec3 centroid, normal; double strike = 0, dip = 0, rms = 0; std::size_t points = 0; };
FittedPlane fit_plane(const std::vector<Polyline>& sticks);

/// The model returned by the implicit builder: a parametric model whose horizons are the
/// implicit surfaces, plus provenance and the misfit at every data point.
class ImplicitModel final : public Model {
 public:
  ImplicitModel(std::shared_ptr<ParametricModel> inner, ModelInfo info, std::vector<DataMisfit> misfit);
  const ModelInfo& info() const override { return info_; }
  void evaluate(std::span<const Vec3> pts, unsigned fields, EvalResult& out) const override { inner_->evaluate(pts, fields, out); }
  void scalar_in_block(std::span<const Vec3> pts, BlockCode block, int horizon, std::span<double> out) const override {
    inner_->scalar_in_block(pts, block, horizon, out);
  }
  void evaluate_fault(int fault, std::span<const Vec3> pts, FaultEval& out) const override { inner_->evaluate_fault(fault, pts, out); }
  void evaluate_fault_geometry(int fault, std::span<const Vec3> pts, FaultEval& out) const override { inner_->evaluate_fault_geometry(fault, pts, out); }
  void restore(std::span<const Vec3> pts, std::span<Vec3> out) const override { inner_->restore(pts, out); }
  /// "structure", "faults" (as the parametric engine), "data_misfit" -> const std::vector<DataMisfit>*.
  const void* capability(std::string_view name) const override;
  const std::vector<DataMisfit>& misfit() const { return misfit_; }

 private:
  std::shared_ptr<ParametricModel> inner_;
  ModelInfo info_;
  std::vector<DataMisfit> misfit_;
};

/// Build an implicit model (also available as the "implicit-rbf" structure-builder plugin).
ModelPtr build_implicit_model(const desc::Description& d, const BuildInputs& in, BuildLog& log);

}  // namespace stratonix
