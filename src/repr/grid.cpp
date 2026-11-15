// Grids and properties (stage 3B, STX-TN-12). Interfaces opened by Prompt 37; corner-point grids
// are implemented by Prompt 38, properties, transmissibility and tetrahedral grids by Prompt 39.
#include "stratonix/repr/grid.hpp"

#include <stdexcept>

namespace stratonix::repr {
namespace {
[[noreturn]] void todo(const char* what, int prompt) {
  throw std::runtime_error(std::string("not implemented yet: ") + what + " (STX-TN-12, Prompt " + std::to_string(prompt) + ")");
}
}  // namespace

const char* layering_name(Layering l) {
  switch (l) {
    case Layering::Proportional: return "proportional";
    case Layering::FollowTop: return "follow_top";
    case Layering::FollowBase: return "follow_base";
  }
  return "?";
}

std::size_t CornerPointGrid::zcorn_index(int, int, int, int) const { todo("corner-point ZCORN", 38); }
Vec3 CornerPointGrid::corner(int, int, int, int) const { todo("corner-point corners", 38); }
std::array<Vec3, 8> CornerPointGrid::corners(int, int, int) const { todo("corner-point corners", 38); }
CornerPointGrid corner_point_grid(const Model&, const GridOptions&) { todo("corner-point grid extraction", 38); }
double cell_volume(const std::array<Vec3, 8>&) { todo("cell volume", 38); }
Vec3 cell_centre(const std::array<Vec3, 8>&) { todo("cell centre", 38); }
std::array<std::array<int, 4>, 6> hex_tets_fixed() { todo("hexahedron split", 38); }
double tet_volume(const Vec3&, const Vec3&, const Vec3&, const Vec3&) { todo("tetrahedron volume", 38); }
GridDefects grid_defects(const Model&, const CornerPointGrid&, double) { todo("grid checks C23-C25", 38); }

const CellProperty* PropertySet::find(const std::string& name) const {
  for (auto& p : props)
    if (p.name == name) return &p;
  return nullptr;
}
CellProperty* PropertySet::find(const std::string& name) {
  for (auto& p : props)
    if (p.name == name) return &p;
  return nullptr;
}
PropertySet grid_properties(const Model&, const CornerPointGrid&, const PropertyInputs&) { todo("grid properties", 39); }
CellProperty apply_transform(const CellProperty&, const TransformSpec&) { todo("property transforms", 39); }
double manzocchi_fault_permeability(double, double) { todo("fault permeability (Manzocchi 1999)", 39); }
double manzocchi_multiplier(double, double, double, double, double, double) { todo("transmissibility multiplier", 39); }
TransmissibilityResult fault_transmissibility(const Model&, const CornerPointGrid&, const PropertySet&, const FaultRockInput&) {
  todo("fault transmissibility", 39);
}
CellProperty upscale_volume_weighted_mean(const CornerPointGrid&, const CellProperty&, const CornerPointGrid&) {
  todo("upscaling (volume-weighted mean)", 39);
}
TetGrid tetrahedral_grid(const CornerPointGrid&) { todo("tetrahedral grid", 39); }

}  // namespace stratonix::repr
