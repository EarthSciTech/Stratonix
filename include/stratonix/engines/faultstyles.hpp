// Stage 4B (STX-TN-15 section 2): polygonal faults, flower structures, duplexes, breached relays,
// trishear on curved faults and growth strata on every hanging wall.
#pragma once
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "stratonix/core/geometry.hpp"

namespace stratonix {

struct PolygonalSpec {
  double x0 = 0, y0 = 0, x1 = 0, y1 = 0;        ///< area
  double unit_top = 0, unit_base = 0;            ///< elevations of the host unit
  double spacing = 500;                          ///< labelled mean fault spacing (m)
  double dip_mean = 50, dip_sd = 5;              ///< degrees
  double throw_ratio = 0.01;                     ///< throw / length (C10)
};
struct PolygonalFault { Vec3 a, b; double dip_deg = 0, throw_m = 0, length = 0; };
struct PolygonalNetwork { std::vector<PolygonalFault> faults; double mean_spacing = 0; };
PolygonalNetwork polygonal_network(const PolygonalSpec& s, std::uint64_t seed);

struct FlowerSpec {
  double horizontal_slip = 0, vertical_component = 0;  ///< m; vertical > 0 positive flower (uplift)
  double spread_deg = 30;
  int splays = 3;
  double branch_depth = 0;
};
struct FlowerResult {
  std::vector<double> dips_deg, vertical_throw, horizontal_slip;  ///< per splay
};
FlowerResult flower_splays(const FlowerSpec& s);

struct DuplexSpec { int horses = 3; double horse_length = 1000, ramp_angle_deg = 30, slip = 600, thickness = 200; };
struct DuplexResult { double shortening = 0, stack_thickness = 0; std::string style; };
DuplexResult duplex(const DuplexSpec& s);

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
