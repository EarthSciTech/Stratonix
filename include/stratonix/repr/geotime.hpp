/// @file geotime.hpp
/// @brief Stage 4A (STX-TN-14 section 5): geological-time coordinates. (u, v) from the model's
/// restoration, t from the restored-frame stratigraphic column and horizon ages; inverse by
/// fixed-point iteration on the restoration; hiatuses where horizons are absent.
#pragma once
#include <optional>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/model/model.hpp"

namespace stratonix {

struct GeoTimePoint {
  double u = kNaN, v = kNaN, t = kNaN;  ///< m, m, Ma (or ordinal when ages are missing)
  bool defined = false;
  int unit = -1;                        ///< unit the point lies in
};

/// A gap in the record at a column: between the ages of two horizons the column lacks.
struct Hiatus { double t_young = kNaN, t_old = kNaN; int horizon_missing = -1; };

class GeologicalTime {
 public:
  /// ages: one per horizon (Ma, young to old); empty = ordinals (the horizon index).
  GeologicalTime(const Model& m, std::vector<double> ages);
  GeoTimePoint forward(const Vec3& x) const;
  /// Present-day point of (u, v, t); nullopt where the iteration does not converge.
  std::optional<Vec3> inverse(double u, double v, double t, double tolerance = 0.01) const;
  std::vector<Hiatus> hiatuses(double u, double v) const;
  bool ordinal() const { return ordinal_; }
  const std::vector<double>& ages() const { return ages_; }

 private:
  const Model& m_;
  std::vector<double> ages_;
  bool ordinal_ = false;
};

}  // namespace stratonix
