/// @file pipeline.hpp
/// @brief High-level operations shared by the command line and the C API.
#pragma once
#include <filesystem>
#include <string>
#include <vector>

#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/repr/faultrock.hpp"

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

/// Key numbers of a built model (summary table; also read by the figures).
std::vector<std::pair<std::string, std::string>> key_numbers(const desc::Description& d, const Model& m, const struct CheckReport& rep,
                                                             class CheckContext& ctx);
std::string summary_csv(const std::vector<std::pair<std::string, std::string>>& keys);
/// Misfit at every data point of a data-conditioned model (empty if the model has no data).
std::string data_misfit_csv(const Model& m);
/// Write the figure set (sections from the description, maps, 3D views, fence) into a folder.
void write_figures(const desc::Description& d, const Model& m, const std::vector<std::pair<std::string, std::string>>& keys,
                   const std::filesystem::path& folder);

/// Draw one section (named in the description, or along `line` when it is not empty) in the
/// four figure formats; the model is built and checked first and nothing is drawn if a check fails.
BuildOutcome run_section(const std::filesystem::path& description, const std::filesystem::path& out_dir, const std::string& name,
                         const Polyline& line, const BuildOptions& opt);
/// Build the model at time `age` (Ma) of a description with events, and write its package
/// (stage 2E, STX-TN-09 section 2).
BuildOutcome run_state(const std::filesystem::path& description, const std::filesystem::path& out_dir, double age,
                       const BuildOptions& opt);
/// Restoration report of a timeline model into `folder`: a section and a map at every event,
/// the sequence as one deck, the flexural-slip unfolded section, balance.csv and the checks of
/// every state (STX-TN-09 sections 2 and 3). Nothing is written for a model without events.
void write_restoration(const desc::Description& d, const Model& m, const std::filesystem::path& folder);

struct PosteriorResult;
/// Stage 4A: posterior_samples.csv and posterior_summary.csv of conditioning.method: bayesian.
void write_posterior_tables(const desc::Description& d, const PosteriorResult& pr, const std::filesystem::path& dir);

/// Export ids from a request ("petrel" expands to the Petrel set; unknown ids are skipped with a
/// warning); an empty request takes the description's outputs.exports.
std::vector<std::string> resolve_export_ids(const desc::Description& d, const std::vector<std::string>& requested);
/// Build, check and export (into out_dir/exports, with README); nothing is exported if a check fails.
/// grid_spacing <= 0 takes the description's grid spacing.
BuildOutcome run_export(const std::filesystem::path& description, const std::filesystem::path& out_dir, const std::vector<std::string>& formats,
                        double grid_spacing, const BuildOptions& opt);

/// Refuse to write into an existing non-empty folder.
void ensure_new_output_dir(const std::filesystem::path& out_dir);

}  // namespace stratonix
