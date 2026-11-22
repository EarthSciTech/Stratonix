/// @file crs.hpp
/// @brief Coordinate reference system record. Transformations need PROJ (optional build).
#pragma once
#include <string>
#include <vector>

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

/// Stage 4C (STX-TN-16 section 2.4): a transformation between two CRS chosen by PROJ, with its
/// name, method and stated accuracy (negative when PROJ states none). Refused without PROJ.
struct CrsOperation { std::string name, method; double accuracy = -1; };
class CrsTransform {
 public:
  CrsTransform(const CrsInfo& from, const CrsInfo& to);
  ~CrsTransform();
  CrsTransform(const CrsTransform&) = delete;
  CrsTransform& operator=(const CrsTransform&) = delete;
  const CrsOperation& operation() const { return op_; }
  /// Transform (x, y, z) in place; x, y in the CRS's axis order as easting, northing (or
  /// longitude, latitude in degrees for a geographic CRS).
  void apply(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z) const;
 private:
  CrsOperation op_;
  void* pj_ = nullptr;
  void* ctx_ = nullptr;
};

}  // namespace stratonix
