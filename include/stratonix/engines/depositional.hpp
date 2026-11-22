// Stage 4B (STX-TN-15 section 3): stratal terminations, clinoform sets, channel belts and lobes.
#pragma once
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/core/grid.hpp"

namespace stratonix {

/// The termination line where horizon z = h(x, y) meets surface z = s(x, y), sampled on a grid.
std::vector<Vec3> onlap_line(const std::function<double(double, double)>& horizon, const std::function<double(double, double)>& surface, const GridSpec& g);

struct ClinoformSpec {
  int count = 5;
  double spacing = 500;          ///< shelf-break advance per clinoform (m)
  double foreset_dip_deg = 5, trajectory_deg = 0;
  double shelf_z = 0, base_z = -200, x_start = 0;
  std::string shape = "oblique";  ///< oblique | sigmoidal | tangential
};
struct Clinoform { double shelf_break_x = 0, shelf_break_z = 0, foreset_dip_deg = 0; std::vector<Vec3> profile; };
std::vector<Clinoform> clinoform_set(const ClinoformSpec& s);

struct UnitBox { double x0 = 0, y0 = 0, x1 = 0, y1 = 0, z_base = 0, z_top = 0; };
struct ChannelSpec {
  int count = 5;
  double width = 200, thickness = 10, wavelength = 2000, amplitude = 0, orientation_deg = 0;
  std::string shape = "rectangular";  ///< rectangular | parabolic
  bool random_placement = true;
};
struct Geobody {
  std::string kind;               ///< channel | lobe
  std::vector<Vec3> centreline;   ///< for channels
  double width = 0, thickness = 0, z_base = 0;
  bool contains(const Vec3& p) const;
};
std::vector<Geobody> place_channels(const ChannelSpec& s, const UnitBox& u, std::uint64_t seed);
/// Volume fraction of the unit inside any geobody, on a lattice of the given spacing.
double net_to_gross(const std::vector<Geobody>& b, const UnitBox& u, double spacing);

}  // namespace stratonix
