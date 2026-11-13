#include "capi_extensions.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix::capi_ext {
void horizon_grid(const Model&, int, double, double, double, double, std::size_t, std::size_t, std::int64_t, double*) {
  throw Error(ErrorCode::NotImplemented, "surface extraction arrives with the representation stage");
}
void export_formats(const std::string&, const std::string&, const std::string&) {
  throw Error(ErrorCode::NotImplemented, "exporters arrive with the export stages");
}
void package(const std::string&, const std::string&, const std::string&) {
  throw Error(ErrorCode::NotImplemented, "packaging arrives with the packaging stage");
}
}  // namespace stratonix::capi_ext
