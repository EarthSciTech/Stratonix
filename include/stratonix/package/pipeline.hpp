/// @file pipeline.hpp
/// @brief High-level operations shared by the command line and the C API.
#pragma once
#include <filesystem>
#include <string>
#include <vector>

#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"

namespace stratonix {

/// Register every built-in plugin of every layer (idempotent). Call once at start-up.
void register_all_plugins();

struct BuildOptions {
  unsigned threads = 0;           ///< 0 = all cores; results never depend on it
  bool reference_kernels = false; ///< use plain reference kernels instead of optimised ones
  bool figures = true;            ///< write the figure set (when the figures layer is available)
  std::vector<std::string> exports;  ///< export formats; empty = description's list
  std::string builder;            ///< preferred structure builder ("" = automatic)
};

struct BuildOutcome {
  ModelPtr model;
  BuildLog log;
  bool checks_passed = true;
  std::vector<std::string> written;  ///< files written, relative to the output folder
};

/// Validate, build and check a model, write summary, check report, exports and figures into
/// `out_dir`, and finish with manifest.json. An existing non-empty output folder is refused, so
/// an issued output is never overwritten.
BuildOutcome run_build(const std::filesystem::path& description, const std::filesystem::path& out_dir,
                       const BuildOptions& opt, const std::string& command = "build");

/// Build only (no files written); used by the C API and tests.
BuildOutcome build_only(const desc::Description& d, const BuildOptions& opt);

/// Refuse to write into an existing non-empty folder.
void ensure_new_output_dir(const std::filesystem::path& out_dir);

}  // namespace stratonix
