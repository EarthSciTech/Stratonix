/* stratonix_c.h - stable C interface to the Stratonix library.
 *
 * Opaque handles, plain C types, status codes; nothing throws across this boundary. The
 * message of the last error on the calling thread is available from stx_last_error().
 * Binary interface: functions are only ever added; STX_CAPI_VERSION_MAJOR changes if a
 * signature or struct layout changes.
 */
#ifndef STRATONIX_C_H
#define STRATONIX_C_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define STX_CAPI_VERSION_MAJOR 1
#define STX_CAPI_VERSION_MINOR 1

/** Status codes (mirror stratonix::ErrorCode; values never change).
 *  STX_OK                  success
 *  STX_E_INVALID_ARGUMENT  NULL pointer, index out of range, bad grid or option
 *  STX_E_IO                file or folder could not be read or written (also: output folder not empty)
 *  STX_E_PARSE             a file could not be parsed (description, grid, options)
 *  STX_E_VALIDATION        the description breaks a schema or semantic rule; input sharing status
 *                          unknown; identifier scan of a package failed
 *  STX_E_BUILD             the model could not be built (for example a tip rule cannot be met)
 *  STX_E_CHECK_FAILED      the model was built but a check failed; nothing was exported
 *  STX_E_NOT_IMPLEMENTED   the function is not available in this version
 *  STX_E_UNSUPPORTED       a format feature that v1 does not support (for example rotated grids)
 *  STX_E_INTERNAL          unexpected failure (please report with stx_last_error()) */
typedef enum stx_status {
  STX_OK = 0,
  STX_E_INVALID_ARGUMENT = 1,
  STX_E_IO = 2,
  STX_E_PARSE = 3,
  STX_E_VALIDATION = 4,
  STX_E_BUILD = 5,
  STX_E_CHECK_FAILED = 6,
  STX_E_NOT_IMPLEMENTED = 7,
  STX_E_UNSUPPORTED = 8,
  STX_E_INTERNAL = 9
} stx_status;

/** Opaque handle to a built model. */
typedef struct stx_model stx_model;

/** Version of the C API (major * 100 + minor) and of the library ("1.0.0"). */
int stx_capi_version(void);
const char* stx_library_version(void);

/** Message of the last error on this thread ("" if none). Valid until the next call. */
const char* stx_last_error(void);
/** Name of a status code, e.g. "STX_E_VALIDATION" (static string). Since C API 1.1. */
const char* stx_status_name(stx_status s);

/** Validate a description file. Returns STX_OK or STX_E_VALIDATION (details in stx_last_error). */
stx_status stx_validate(const char* description_path);

/** Build a model from a description file. On success *out holds a handle to free with stx_model_free. */
stx_status stx_model_build(const char* description_path, stx_model** out);
void stx_model_free(stx_model* m);

/** Counts of horizons, units (horizons + 1) and faults. */
stx_status stx_model_counts(const stx_model* m, int* horizons, int* units, int* faults);
/** Name of horizon i (pointer valid while the model lives). */
const char* stx_model_horizon_name(const stx_model* m, int i);
/** Name of fault i, oldest first (pointer valid while the model lives). Since C API 1.1. */
const char* stx_model_fault_name(const stx_model* m, int i);

/** Evaluate n points (xyz interleaved, metres, z elevation). Any output may be NULL.
 *  unit_out[n], block_out[n], scalar_out[n * horizons] (horizon-major). */
stx_status stx_model_evaluate(const stx_model* m, size_t n, const double* xyz, int32_t* unit_out,
                              uint64_t* block_out, double* scalar_out);

/** Horizon elevations on a regular grid (nx * ny, x fastest); NaN in fault gaps.
 *  block: -1 for the faulted surface with gaps, otherwise a fault-block code. */
stx_status stx_model_horizon_grid(const stx_model* m, int horizon, double x0, double y0, double dx, double dy,
                                  size_t nx, size_t ny, int64_t block, double* z_out);

/** Build, check and write all outputs of a description into a new folder (as `stratonix build`). */
stx_status stx_build_to_folder(const char* description_path, const char* out_dir);
/** Export a description's model into a new folder; formats is a comma-separated list of exporter
 *  ids or sets ("petrel", "vtk,csv", ...); NULL or "" takes the description's list. */
stx_status stx_export(const char* description_path, const char* out_dir, const char* formats);
/** Write a delivery package (as `stratonix package`). options_json may be NULL or a JSON object
 *  with the keys of the package configuration (formats, figures, zip, terms_of_use, contact,
 *  assumptions, limitations, banned_patterns, translations, default_patterns). */
stx_status stx_package(const char* description_path, const char* out_dir, const char* options_json);

#ifdef __cplusplus
}
#endif
#endif /* STRATONIX_C_H */
