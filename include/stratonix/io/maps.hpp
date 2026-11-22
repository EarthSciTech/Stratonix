/// @file maps.hpp
/// @brief Stage 4C (STX-TN-16 section 2.5): geological maps (shapefile, GeoPackage), outcrop
/// measurements and elevation models as model inputs.
#pragma once
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/core/grid.hpp"

namespace stratonix::io {

struct MapFeature {
  int id = 0;
  int type = 0;                                  ///< 1 point, 3 polyline, 5 polygon (shapefile codes)
  std::vector<std::vector<Vec3>> parts;          ///< rings or lines
  std::map<std::string, std::string> attributes;
};
struct MapLayer { std::string source, name; std::vector<MapFeature> features; std::vector<std::string> report; };
MapLayer read_shapefile(const std::filesystem::path& shp);
MapLayer read_geopackage(const std::filesystem::path& gpkg, const std::string& layer);
Grid2D read_esri_ascii_grid(const std::filesystem::path& file);
Grid2D read_geotiff_dem(const std::filesystem::path& file);
/// Write a shapefile (with .shx and .dbf) - used by the tests and the exports.
void write_shapefile(const std::filesystem::path& shp, const MapLayer& layer);

struct Measurement { Vec3 p; double dip = 0, azimuth = 0; std::string kind, unit, source, date; };  ///< azimuth = dip direction (rad)
struct ModelData {
  std::vector<std::pair<std::string, Vec3>> contact_points;    ///< horizon name, point on the elevation model
  std::vector<Measurement> orientations;
  std::vector<std::pair<std::string, std::vector<Vec3>>> fault_traces;
  std::vector<std::string> report;
};
/// Contacts between map units (youngest first in `unit_order`) draped on the elevation model become
/// points on the horizon at the base of the younger unit; faults become traces.
ModelData map_to_model_data(const MapLayer& units, const std::string& unit_field, const std::vector<std::string>& unit_order,
                            const MapLayer* faults, const std::vector<Measurement>& measurements, const Grid2D& dem);

}  // namespace stratonix::io
