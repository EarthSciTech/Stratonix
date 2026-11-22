// Stage 4C stubs (STX-TN-16 section 2.5): implemented by Prompt 56.
#include "stratonix/io/maps.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix::io {

MapLayer read_shapefile(const std::filesystem::path&) { throw Error(ErrorCode::NotImplemented, "read_shapefile: not implemented yet (Prompt 56)"); }
MapLayer read_geopackage(const std::filesystem::path&, const std::string&) {
  throw Error(ErrorCode::NotImplemented, "read_geopackage: not implemented yet (Prompt 56)");
}
Grid2D read_esri_ascii_grid(const std::filesystem::path&) { throw Error(ErrorCode::NotImplemented, "read_esri_ascii_grid: not implemented yet (Prompt 56)"); }
Grid2D read_geotiff_dem(const std::filesystem::path&) { throw Error(ErrorCode::NotImplemented, "read_geotiff_dem: not implemented yet (Prompt 56)"); }
void write_shapefile(const std::filesystem::path&, const MapLayer&) { throw Error(ErrorCode::NotImplemented, "write_shapefile: not implemented yet (Prompt 56)"); }
ModelData map_to_model_data(const MapLayer&, const std::string&, const std::vector<std::string>&, const MapLayer*, const std::vector<Measurement>&,
                            const Grid2D&) {
  throw Error(ErrorCode::NotImplemented, "map_to_model_data: not implemented yet (Prompt 56)");
}

}  // namespace stratonix::io
