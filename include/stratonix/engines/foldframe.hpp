/// @file foldframe.hpp
/// @brief Stage 4A (STX-TN-14 section 2): fold frames fitted to bedding (and axial foliation)
/// orientations: fold axis, axial surface and the fold rotation angle as a Fourier series of the
/// distance to the axial surface; predicted bedding normals constrain the series' field.
#pragma once
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"

namespace stratonix {

/// An orientation measurement: position and unit normal (pointing up for bedding).
struct OrientationDatum {
  Vec3 p, n;
};

struct FoldFrameOptions {
  double max_cylindricity = 0.2;  ///< above: axis varies along the fold (reported)
  int max_order = 3;              ///< Fourier order chosen by AIC up to this
  double wavelength = 0;          ///< m start value; 0 = estimated from the hinge spacing
};

struct FoldFrame {
  bool ok = false;
  std::string reason;            ///< why not ok
  Vec3 axis{0, 1, 0};            ///< unit fold axis (plunging down, or horizontal)
  double plunge = 0, trend = 0;  ///< rad
  double cylindricity = 0;       ///< smallest / middle eigenvalue of the orientation tensor
  Vec3 axial_point{0, 0, 0};     ///< a point on the axial surface
  Vec3 axial_normal{1, 0, 0};    ///< unit normal of the (planar) axial surface
  double wavelength = 0;         ///< m, fitted
  double interlimb = 0;          ///< rad, interlimb angle from the limb dips
  std::vector<double> fourier;   ///< a0, a1, b1, a2, b2, ... of the fold rotation angle (rad)
  int order = 0;
  double rms = 0;                ///< rad, misfit of the bedding normals
  std::size_t data = 0;
  /// Signed distance to the axial surface (the fold-frame coordinate).
  double frame_coordinate(const Vec3& x) const;
  /// Fold rotation angle (rad) at fold-frame coordinate c.
  double rotation_angle(double c) const;
  /// Predicted unit bedding normal at x (pointing up).
  Vec3 bedding_normal(const Vec3& x) const;
};

/// Fit a fold frame to bedding orientations, with axial foliation when measured.
FoldFrame fit_fold_frame(const std::vector<OrientationDatum>& bedding, const std::vector<OrientationDatum>& foliation,
                         const FoldFrameOptions& o);

}  // namespace stratonix
