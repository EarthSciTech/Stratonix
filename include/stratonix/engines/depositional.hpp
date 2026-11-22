// Stage 4B (STX-TN-15 section 3): stratal terminations, clinoform sets, channel belts and lobes.
#pragma once
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/core/grid.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/model/model.hpp"

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
  std::string shape = "rectangular";  ///< channels: rectangular | parabolic (thickness at the axis)
  Vec3 centre{};                  ///< lobes: plan centre
  double radius = 0, aspect = 1, orientation = 0;  ///< lobes: semi-major axis, major/minor ratio, azimuth (rad)
  double x0 = 0, y0 = 0, x1 = 0, y1 = 0;  ///< plan bounding box (fast rejection)
  bool contains(const Vec3& p) const;
  double volume() const;          ///< closed form (unclipped)
};
std::vector<Geobody> place_channels(const ChannelSpec& s, const UnitBox& u, std::uint64_t seed);

struct LobeSpec {
  int count = 5;
  double radius = 1000, aspect = 2, thickness = 5, orientation_deg = 0;
};
/// Lobes: elliptical in plan, a flat base and a lens-shaped top (thickness t (1 - r^2)), seeded.
std::vector<Geobody> place_lobes(const LobeSpec& s, const UnitBox& u, std::uint64_t seed);
/// Volume fraction of the unit inside any geobody, on a lattice of the given spacing.
double net_to_gross(const std::vector<Geobody>& b, const UnitBox& u, double spacing);

/// A set of geobodies in a model (capability "geobodies"): one facies unit per set.
struct GeobodyReport {
  std::string name, type, host;
  int unit = -1;              ///< the set's unit index in the model
  int count = 0;
  double ntg = 0;             ///< volume fraction of the host unit in the set (lattice, geological-time box)
  double ntg_closed = 0;      ///< sum of closed-form volumes / box volume (no overlap or clipping)
  double box_thickness = 0;   ///< m: the unit's restored thickness at the domain centre (the box height)
};
/// Wrap a built model with the description's geobodies (STX-TN-15 section 3.3). A point of the host
/// unit is mapped to its restored position and its fraction of the unit's restored thickness; the
/// geobodies are placed in that box, so they follow the unit's deformation. Later sets take a point
/// before earlier ones (erosional).
ModelPtr wrap_geobodies(const desc::Description& d, ModelPtr base, BuildLog& log);

}  // namespace stratonix
