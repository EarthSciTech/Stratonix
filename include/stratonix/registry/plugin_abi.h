/* plugin_abi.h - stable C interface for Stratonix plugins loaded at run time (from v2).
 *
 * A plugin shared library (.dll or .so) exports one function:
 *
 *   STX_PLUGIN_EXPORT int stratonix_plugin_entry(const stx_host_api* host, stx_plugin_info* info);
 *
 * It fills `info` and returns 0 on success. The host refuses a plugin whose abi_major differs
 * from STX_PLUGIN_ABI_MAJOR. No C++ types cross this boundary; strings are UTF-8 and owned by
 * the plugin for the life of the process. In Stratonix v1 the loader is not enabled; the header
 * fixes the contract so that plugins written later do not need a host rebuild.
 */
#ifndef STRATONIX_PLUGIN_ABI_H
#define STRATONIX_PLUGIN_ABI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define STX_PLUGIN_ABI_MAJOR 1
#define STX_PLUGIN_ABI_MINOR 0

#if defined(_WIN32)
#define STX_PLUGIN_EXPORT __declspec(dllexport)
#else
#define STX_PLUGIN_EXPORT __attribute__((visibility("default")))
#endif

typedef enum stx_plugin_kind {
  STX_PLUGIN_STRUCTURE_BUILDER = 0,
  STX_PLUGIN_KINEMATICS = 1,
  STX_PLUGIN_HANGING_WALL = 2,
  STX_PLUGIN_EXPORTER = 3,
  STX_PLUGIN_IMPORTER = 4,
  STX_PLUGIN_EXTRACTOR = 5,
  STX_PLUGIN_CHECK = 6
} stx_plugin_kind;

/* Services the host offers to plugins. */
typedef struct stx_host_api {
  uint32_t abi_major, abi_minor;
  void (*log)(int level, const char* message);
  /* Evaluate the current model: n points (xyz interleaved), returns unit index per point. */
  int (*evaluate_units)(void* model, size_t n, const double* xyz, int32_t* unit_out);
} stx_host_api;

/* Exporter plugins: write one model to a folder. */
typedef struct stx_exporter_vtable {
  int (*write)(void* model, const char* folder, const char* options_json, char* err, size_t err_len);
} stx_exporter_vtable;

typedef struct stx_plugin_info {
  uint32_t abi_major, abi_minor;
  stx_plugin_kind kind;
  const char* name;
  const char* version;      /* semantic version of the plugin */
  const char* description;
  const void* vtable;       /* kind-specific table, e.g. stx_exporter_vtable */
} stx_plugin_info;

typedef int (*stx_plugin_entry_fn)(const stx_host_api* host, stx_plugin_info* info);

#ifdef __cplusplus
}
#endif
#endif /* STRATONIX_PLUGIN_ABI_H */
