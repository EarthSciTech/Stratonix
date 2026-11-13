/// @file crs.hpp
/// @brief Coordinate reference system record. Transformations need PROJ (optional build).
#pragma once
#include <string>

namespace stratonix {

/// Coordinate reference system of a model. A local frame has epsg == 0 and no WKT.
struct CrsInfo {
  std::string name;            ///< human-readable name
  int epsg = 0;                ///< EPSG code of the projected CRS, 0 = local frame
  std::string wkt;             ///< optional WKT
  std::string vertical_datum;  ///< description of the vertical datum (z = 0)
  std::string length_unit = "m";
  bool is_local() const { return epsg == 0 && wkt.empty(); }
  /// One-line label used in figure footers and README files.
  std::string label() const {
    std::string s = name;
    if (epsg > 0) s += " (EPSG:" + std::to_string(epsg) + ")";
    return s;
  }
};

/// True when the library was built with PROJ (reprojection available).
bool crs_reprojection_available();

}  // namespace stratonix
