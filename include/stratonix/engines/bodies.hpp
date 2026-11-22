// Stage 4B (STX-TN-15 section 1): multi-valued and cross-cutting bodies - salt sheets, sills,
// saucer sills, dykes and sand injectites.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"

namespace stratonix {

/// An allochthonous salt sheet: footprint polygon (x, y), flat base, thickness tapering linearly to
/// zero over `taper` metres inside the footprint edge.
struct SheetBody {
  std::vector<Vec3> footprint;  ///< polygon vertices (z ignored)
  double base_z = 0, thickness = 0, taper = 0;
  /// thickness at (x, y): 0 outside the footprint
  double thickness_at(double x, double y) const;
};
/// Volume of a sheet integrated on a grid of the given spacing (m^3).
double sheet_volume(const SheetBody& s, double spacing);

/// A circular sill with the elliptical crack profile h(r) = h0 sqrt(1 - r^2 / a^2) and the forced
/// fold above it, w(r) = w0 (1 - r^2 / a^2)^2, with w0 set by volume and the accommodation factor.
struct SillBody {
  Vec3 centre;              ///< centre of the sill's mid-plane
  double radius = 0;        ///< a (m)
  double thickness = 0;     ///< h0, central thickness (m)
  double accommodation = 1; ///< f in [0, 1]: share of the sill volume made by roof uplift
  double thickness_at(double r) const;
};
double sill_volume(const SillBody& s);                       ///< closed form (2/3) pi a^2 h0
double forced_fold_uplift(const SillBody& s, double r);      ///< w(r) at the sill roof
double forced_fold_volume(const SillBody& s, double spacing); ///< uplift integrated on a grid

/// Sand injectite complex sourced from a parent unit (sheets seeded and reproducible).
struct InjectiteSpec {
  double parent_top = 0, ceiling = 0;      ///< elevations: parent top and ceiling horizon
  double parent_sand_volume = 0;           ///< m^3
  double max_fraction = 0.3;               ///< injected volume / parent sand, at most
  int sheets = 10;
  double length_mean = 500, length_sd = 100, thickness_mean = 5, thickness_sd = 1;  ///< m
  double climb_deg = 30;                   ///< wing climb angle
  BBox area;                               ///< placement area (x, y)
};
struct InjectiteSheet { Vec3 start, end; double thickness = 0, width = 0; };
struct InjectiteComplex {
  bool ok = false;
  std::string reason;
  std::vector<InjectiteSheet> sheets;
  double injected_volume = 0, parent_volume_after = 0;
};
InjectiteComplex build_injectites(const InjectiteSpec& s, std::uint64_t seed);

}  // namespace stratonix
