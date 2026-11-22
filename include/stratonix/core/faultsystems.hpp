// Stage 4B (STX-TN-15 sections 2.1 to 2.3): planners of fault networks. A description's fault
// systems expand into ordinary faults with these (desc layer), so they sit in core.
#pragma once
#include <cstdint>
#include <string>
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

}  // namespace stratonix
