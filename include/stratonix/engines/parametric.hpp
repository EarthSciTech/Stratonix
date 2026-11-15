/// @file parametric.hpp
/// @brief The v1 parametric engine: unfaulted structure + fault network + frame rule, exposed
/// through the Model interface.
#pragma once
#include <algorithm>
#include <memory>

#include "stratonix/engines/faults.hpp"
#include "stratonix/engines/structure.hpp"
#include "stratonix/model/model.hpp"

namespace stratonix {

class ParametricModel final : public Model {
 public:
  ParametricModel(const desc::Description& d, const BuildInputs& in, BuildLog& log);
  const ModelInfo& info() const override { return info_; }
  void evaluate(std::span<const Vec3> pts, unsigned fields, EvalResult& out) const override;
  void scalar_in_block(std::span<const Vec3> pts, BlockCode block, int horizon,
                       std::span<double> out) const override;
  void evaluate_fault(int fault, std::span<const Vec3> pts, FaultEval& out) const override;
  void restore(std::span<const Vec3> pts, std::span<Vec3> out) const override;
  /// "structure" -> const Structure*, "faults" -> const FaultNetwork*.
  const void* capability(std::string_view name) const override;

  const Structure& structure() const { return *structure_; }
  const FaultNetwork& faults() const { return *faults_; }

 private:
  /// Evaluate one point: restored position, block code and scalars (reference kernel).
  void eval_point(const Vec3& p, BlockCode* block, Vec3* restored, double* scalars) const;
  /// capability("porosity") (stage 3B): the compaction law of each unit over restored depths.
  struct Porosity final : PorosityQuery {
    const ParametricModel* m = nullptr;
    bool has_law(int unit) const override;
    double mean_porosity(int unit, const Vec3& top, const Vec3& base) const override;
  };
  Porosity porosity_;
  ModelInfo info_;
  std::shared_ptr<Structure> structure_;
  std::shared_ptr<FaultNetwork> faults_;
  /// Erosional surfaces that no fault displaces (stage 2E): they truncate the older horizons in
  /// the present geometry, after the faulting they postdate.
  std::vector<std::size_t> post_erosion_;
  /// True when p lies above an erosional surface no fault displaces: younger material that never
  /// moved (its restoration is the identity).
  bool above_post_erosion(const Vec3& p) const {
    for (std::size_t u : post_erosion_)
      if (p.z >= structure_->z(u, p.x, p.y) - 1e-6) return true;
    return false;
  }
  void post_truncate(double* scalars) const {
    for (std::size_t u : post_erosion_)
      for (std::size_t k = u + 1; k < structure_->horizons(); ++k) scalars[k] = std::max(scalars[k], scalars[u]);
  }
};

}  // namespace stratonix
