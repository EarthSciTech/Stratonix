/// @file wells.hpp
/// @brief Well paths from description entries.
#pragma once
#include <vector>

#include "stratonix/core/well.hpp"
#include "stratonix/desc/description.hpp"

namespace stratonix {

/// Path of a described well. Tops are depths below the model datum (z = 0). Measured depth is
/// counted from the well datum (datum_elevation, or the model datum if not given). Trajectory
/// stations are [md (m), inclination (deg), azimuth (deg)].
WellPath make_well_path(const desc::WellDesc& w, const desc::Description& d);
std::vector<WellPath> make_well_paths(const desc::Description& d);

}  // namespace stratonix
