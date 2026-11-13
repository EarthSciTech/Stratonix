/// @file manifest.hpp
/// @brief Run manifest: inputs with SHA-256, resolved description, build information, outputs.
/// Timestamps appear only in the manifest, so every other output is byte-identical between runs.
#pragma once
#include <filesystem>
#include <string>
#include <vector>

#include "stratonix/desc/description.hpp"

namespace stratonix {

struct ManifestFile {
  std::string path;     ///< relative to the output folder, forward slashes
  std::string sha256;
  std::uintmax_t bytes = 0;
  std::string role;     ///< e.g. "description", "input", "export", "figure", "report"
  std::string sharing;  ///< inputs only
};

struct Manifest {
  std::string tool = "stratonix";
  std::string command;
  std::string model_name, model_version, evidence_basis;
  std::string description_path, description_sha256;
  desc::Node resolved_description;
  std::vector<ManifestFile> inputs;
  std::vector<ManifestFile> outputs;
  std::vector<std::string> plugins;  ///< "kind/name@version"
  std::vector<std::pair<std::string, std::string>> extra;  ///< key/value notes (e.g. corrections)
  std::string timestamp_utc;         ///< ISO 8601; the only time-dependent field

  /// Fill description, inputs and build information from a description.
  static Manifest from_description(const desc::Description& d, const std::string& command);
  /// Add every regular file under `root` (except the manifest itself) as an output, sorted.
  void add_outputs(const std::filesystem::path& root, const std::string& role_default = "output");
  std::string to_json() const;
  void write(const std::filesystem::path& file) const;
};

/// Current UTC time as ISO 8601 (seconds).
std::string utc_timestamp();

/// Verify a manifest.json against the files in its folder. Returns a list of problems (empty = ok).
std::vector<std::string> verify_manifest(const std::filesystem::path& manifest_file);

}  // namespace stratonix
