// Functions the C API forwards to; they grow as stages add representations, exports and packaging.
#pragma once
#include <cstddef>
#include <cstdint>
#include <string>

#include "stratonix/model/model.hpp"

namespace stratonix::capi_ext {
void horizon_grid(const Model& m, int horizon, double x0, double y0, double dx, double dy, std::size_t nx,
                  std::size_t ny, std::int64_t block, double* z_out);
void export_formats(const std::string& description, const std::string& out_dir, const std::string& formats);
void package(const std::string& description, const std::string& out_dir, const std::string& options_json);
}  // namespace stratonix::capi_ext
