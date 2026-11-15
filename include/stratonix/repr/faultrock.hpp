// Faults as rock and the stress field (stage 2F, STX-TN-10): juxtaposition, shale gouge ratio,
// clay smear potential, fault-zone thickness, seal capacity, stress, slip and dilation tendency,
// and horizon curvature with fracture trends. Inputs are plain values (no description types).
#pragma once
#include <string>
#include <vector>

#include "stratonix/core/grid.hpp"
#include "stratonix/model/model.hpp"

namespace stratonix::repr {

struct FaultRockInput {
  std::vector<double> unit_clay;  ///< clay fraction per model unit (NaN: not given)
  double shale_threshold = kNaN;  ///< clay fraction of a shale bed for clay smear
  double thickness_ratio = kNaN;  ///< fault-zone thickness / displacement
  bool seal = false;              ///< sgr_power calibration given
  double seal_a = kNaN, seal_c = kNaN;
};
struct StressInput {
  bool present = false;
  double g[6] = {0, 0, 0, 0, 0, 0};  ///< total stress per metre of depth (Pa/m): xx, yy, zz, xy, xz, yz (east, north, up)
  double pp = 0;                     ///< pore pressure per metre of depth (Pa/m)
  double friction = kNaN, cohesion = 0;
};
/// Attributes at one point of a fault surface. Missing values are NaN, with the reasons in `missing`.
struct FaultRockPoint {
  int fw_unit = -1, hw_unit = -1;
  double d_fw = kNaN, d_hw = kNaN, separation = kNaN;  ///< stratigraphic positions (m below the top horizon)
  double sgr = kNaN, csp = kNaN, thickness = kNaN, seal_bar = kNaN;
  double sigma_n = kNaN, tau = kNaN, ts = kNaN, td = kNaN, dp = kNaN;  ///< Pa; tendencies dimensionless
  double ts_max = kNaN;
  int fw_unit_col = -1, hw_unit_col = -1;  ///< units at d_fw, d_hw in the restored columns (C19)
  bool near_horizon = false;               ///< within 0.5 m of a horizon on either side (C19 skips)
  bool valid = true;                       ///< false for lattice nodes outside the tip line or the domain
  std::string missing;
};
FaultRockPoint fault_rock_at(const Model& m, int fault, const Vec3& p, const FaultRockInput& in, const StressInput& s);

struct FaultSurface;
/// Attributes at every lattice node of a fault surface (invalid outside `max_extent` or the domain).
std::vector<FaultRockPoint> fault_rock_surface(const Model& m, const FaultSurface& fs, const FaultRockInput& in, const StressInput& s,
                                               double max_extent = 1.0);

/// Curvature attributes of a horizon (Roberts 2001), 1/m, and the fracture-trend azimuth (deg;
/// NaN where the curvature is below `fracture_threshold`).
struct CurvatureGrids {
  Grid2D mean, gaussian, kmax, kmin, most_pos, most_neg, dip, strike, trend;
};
CurvatureGrids curvature(const Model& m, int horizon, const GridSpec& g, double fracture_threshold);

}  // namespace stratonix::repr
