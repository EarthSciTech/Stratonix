// C API implementation: converts exceptions to status codes and a thread-local message.
#include "stratonix_c.h"

#include <cmath>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/package/pipeline.hpp"
#include "capi_extensions.hpp"

using namespace stratonix;

struct stx_model {
  ModelPtr model;
  desc::Description description;
};

namespace {
thread_local std::string g_err;

template <class F>
stx_status guard(F&& f) {
  g_err.clear();
  try {
    f();
    return STX_OK;
  } catch (const Error& e) {
    g_err = e.what();
    return static_cast<stx_status>(static_cast<int>(e.code()));
  } catch (const std::exception& e) {
    g_err = e.what();
    return STX_E_INTERNAL;
  } catch (...) {
    g_err = "unknown error";
    return STX_E_INTERNAL;
  }
}
}  // namespace

extern "C" {

int stx_capi_version(void) { return STX_CAPI_VERSION_MAJOR * 100 + STX_CAPI_VERSION_MINOR; }
const char* stx_library_version(void) { return build_info().version.c_str(); }
const char* stx_last_error(void) { return g_err.c_str(); }

stx_status stx_validate(const char* path) {
  return guard([&] {
    if (!path) throw Error(ErrorCode::InvalidArgument, "path is NULL");
    auto r = desc::validate_description(path);
    if (!r.ok()) throw Error(ErrorCode::Validation, r.text());
  });
}

stx_status stx_model_build(const char* path, stx_model** out) {
  return guard([&] {
    if (!path || !out) throw Error(ErrorCode::InvalidArgument, "NULL argument");
    *out = nullptr;
    auto h = std::make_unique<stx_model>();
    h->description = desc::load_description(path);
    BuildOptions o;
    h->model = build_only(h->description, o).model;
    *out = h.release();
  });
}

void stx_model_free(stx_model* m) { delete m; }

stx_status stx_model_counts(const stx_model* m, int* h, int* u, int* f) {
  return guard([&] {
    if (!m) throw Error(ErrorCode::InvalidArgument, "model is NULL");
    auto& i = m->model->info();
    if (h) *h = static_cast<int>(i.horizons.size());
    if (u) *u = static_cast<int>(i.units.size());
    if (f) *f = static_cast<int>(i.faults.size());
  });
}

const char* stx_model_horizon_name(const stx_model* m, int i) {
  if (!m || i < 0 || i >= static_cast<int>(m->model->info().horizons.size())) return nullptr;
  return m->model->info().horizons[static_cast<std::size_t>(i)].name.c_str();
}

stx_status stx_model_evaluate(const stx_model* m, size_t n, const double* xyz, int32_t* unit_out,
                              uint64_t* block_out, double* scalar_out) {
  return guard([&] {
    if (!m || (!xyz && n)) throw Error(ErrorCode::InvalidArgument, "NULL argument");
    std::vector<Vec3> pts(n);
    for (size_t i = 0; i < n; ++i) pts[i] = {xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]};
    unsigned f = (unit_out ? kEvalUnit : 0u) | (block_out ? kEvalBlock : 0u) | (scalar_out ? kEvalScalars : 0u);
    EvalResult r;
    m->model->evaluate(pts, f, r);
    if (unit_out) std::memcpy(unit_out, r.unit.data(), n * sizeof(int32_t));
    if (block_out) std::memcpy(block_out, r.block.data(), n * sizeof(uint64_t));
    if (scalar_out) std::memcpy(scalar_out, r.scalar.data(), r.scalar.size() * sizeof(double));
  });
}

stx_status stx_model_horizon_grid(const stx_model* m, int horizon, double x0, double y0, double dx, double dy,
                                  size_t nx, size_t ny, int64_t block, double* z_out) {
  return guard([&] {
    if (!m || !z_out) throw Error(ErrorCode::InvalidArgument, "NULL argument");
    capi_ext::horizon_grid(*m->model, horizon, x0, y0, dx, dy, nx, ny, block, z_out);
  });
}

stx_status stx_build_to_folder(const char* path, const char* out_dir) {
  return guard([&] {
    if (!path || !out_dir) throw Error(ErrorCode::InvalidArgument, "NULL argument");
    BuildOptions o;
    auto r = run_build(path, out_dir, o, "capi:build");
    if (!r.checks_passed) throw Error(ErrorCode::CheckFailed, "one or more checks failed; see checks_report.txt");
  });
}

stx_status stx_export(const char* path, const char* out_dir, const char* formats) {
  return guard([&] {
    if (!path || !out_dir) throw Error(ErrorCode::InvalidArgument, "NULL argument");
    capi_ext::export_formats(path, out_dir, formats ? formats : "");
  });
}

stx_status stx_package(const char* path, const char* out_dir, const char* options_json) {
  return guard([&] {
    if (!path || !out_dir) throw Error(ErrorCode::InvalidArgument, "NULL argument");
    capi_ext::package(path, out_dir, options_json ? options_json : "");
  });
}

}  // extern "C"
