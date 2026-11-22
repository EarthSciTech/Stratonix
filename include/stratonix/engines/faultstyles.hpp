// Stage 4B (STX-TN-15 section 2): polygonal faults, flower structures, duplexes, breached relays,
// trishear on curved faults and growth strata on every hanging wall. The network planners
// (polygonal, flower, duplex) live in core/faultsystems.hpp: descriptions expand fault systems with them.
#pragma once
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "stratonix/core/faultsystems.hpp"
#include "stratonix/core/geometry.hpp"

namespace stratonix {

struct RelaySpec { double length_a = 0, length_b = 0, overlap = 0, separation = 0, throw_a = 0, throw_b = 0; };
/// Summed throw along strike across a breached relay: (along-strike position, summed throw).
std::vector<std::pair<double, double>> breached_relay_throw(const RelaySpec& s, int samples);

struct TrishearParams { double apical_deg = 60, ps_ratio = 1.5, slip = 0, increment = 5; };
/// Points (x, z) displaced by trishear with the tip propagating along the fault polyline (x, z).
std::vector<Vec3> trishear_curved(const std::vector<Vec3>& fault, const TrishearParams& p, const std::vector<Vec3>& points);

struct GrowthSpec {
  std::string hanging_wall = "rigid";  ///< layer_rigid | rigid | shear | fault_bend | fault_propagation | trishear
  double dip_deg = 60, slip_rate = 0;  ///< m per unit time
  std::vector<double> deposition_times;  ///< growth unit boundaries (time)
};
/// Thickness of each growth unit across the fault (hanging wall minus footwall).
std::vector<double> growth_thickness(const GrowthSpec& s);

}  // namespace stratonix
