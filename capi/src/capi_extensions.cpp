#include "capi_extensions.hpp"

#include <cstring>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/package/delivery.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

namespace stratonix::capi_ext {

void horizon_grid(const Model& m, int horizon, double x0, double y0, double dx, double dy, std::size_t nx, std::size_t ny,
                  std::int64_t block, double* z_out) {
  if (horizon < 0 || horizon >= static_cast<int>(m.info().horizons.size()))
    throw Error(ErrorCode::InvalidArgument, "horizon index out of range");
  if (!(dx > 0) || !(dy > 0) || nx == 0 || ny == 0) throw Error(ErrorCode::InvalidArgument, "grid needs dx, dy > 0 and nx, ny >= 1");
  GridSpec g;
  g.x0 = x0; g.y0 = y0; g.dx = dx; g.dy = dy; g.nx = nx; g.ny = ny;
  if (block < 0) {
    auto s = repr::horizon_grid_faulted(m, horizon, g);
    std::memcpy(z_out, s.z.v.data(), g.size() * sizeof(double));
  } else {
    auto s = repr::horizon_grid_block(m, horizon, static_cast<BlockCode>(block), g);
    std::memcpy(z_out, s.v.data(), g.size() * sizeof(double));
  }
}

void export_formats(const std::string& description, const std::string& out_dir, const std::string& formats) {
  std::vector<std::string> f;
  for (auto& x : split(formats, ','))
    if (!trim(x).empty()) f.push_back(trim(x));
  BuildOptions o;
  auto r = run_export(description, out_dir, f, 0, o);
  if (!r.checks_passed) throw Error(ErrorCode::CheckFailed, "one or more checks failed; nothing exported (see checks_report.txt)");
}

void package(const std::string& description, const std::string& out_dir, const std::string& options_json) {
  PackageConfig cfg = options_json.empty() ? PackageConfig{} : PackageConfig::parse(options_json);
  BuildOptions o;
  auto r = run_package(description, out_dir, cfg, o);
  if (!r.scan_passed)
    throw Error(ErrorCode::Validation, "identifier scan failed (" + std::to_string(r.hits.size()) + " hits); package rejected: " + r.folder.string());
}

}  // namespace stratonix::capi_ext
