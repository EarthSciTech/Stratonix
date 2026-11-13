/// @file parametric.hpp
/// @brief The v1 parametric engine: unfaulted structure + fault network + frame rule, exposed
/// through the Model interface.
#pragma once
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
  ModelInfo info_;
  std::shared_ptr<Structure> structure_;
  std::shared_ptr<FaultNetwork> faults_;
};

}  // namespace stratonix
