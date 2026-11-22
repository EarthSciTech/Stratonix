/// @file well.hpp
/// @brief Well paths: vertical or deviated (stations of measured depth, inclination, azimuth),
/// computed by the minimum-curvature method.
#pragma once
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"

namespace stratonix {

struct WellStation { double md, inc, azi; };  ///< metres, radians, radians

/// A well path from its head (x, y) at elevation `datum_z`. Depths are positive down below the
/// model datum (z = 0); elevations z = -depth.
class WellPath {
 public:
  WellPath() = default;
  /// Vertical well.
  WellPath(std::string name, double x, double y, double datum_z, double total_depth);
  /// Deviated well from survey stations (first station at md 0).
  WellPath(std::string name, double x, double y, double datum_z, std::vector<WellStation> stations);
  /// Stage 4C: a well from positions computed by a stated survey method: measured depth and
  /// (north, east, TVD) offsets from the head, TVD positive down.
  static WellPath from_positions(std::string name, double x, double y, double datum_z, const std::vector<double>& md,
                                 const std::vector<Vec3>& north_east_tvd, std::string method);
  /// Point at a measured depth (interpolated along the path; vertical beyond total depth).
  Vec3 at_md(double md) const;
  const std::string& survey_method() const { return method_; }
  const std::string& name() const { return name_; }
  double x() const { return x_; }
  double y() const { return y_; }
  double datum_z() const { return datum_z_; }
  bool vertical() const { return stations_.empty(); }
  /// Path points (x, y, z) from the head to total depth.
  const Polyline& points() const { return pts_; }
  /// Horizontal position at an elevation (interpolated along the path).
  Vec3 at_elevation(double z) const;
  /// Measured depth at a point index-interpolated along the path (for tops reports).
  double md_at_elevation(double z) const;
  double total_depth_md() const { return md_.empty() ? 0 : md_.back(); }

 private:
  std::string name_;
  double x_ = 0, y_ = 0, datum_z_ = 0;
  std::vector<WellStation> stations_;
  Polyline pts_;
  std::vector<double> md_;
  std::string method_;
};

}  // namespace stratonix
