#pragma once
// Delivery package: a self-contained, checksummed folder (and optional zip) for an external
// recipient, with an identifier scan that refuses internal labels.
#include <filesystem>
#include <string>
#include <vector>

#include "stratonix/desc/description.hpp"
#include "stratonix/package/pipeline.hpp"

namespace stratonix {

/// Package settings, read from a YAML file (all keys optional; see docs/package_config.md).
struct PackageConfig {
  std::vector<std::string> formats = {"petrel"};  ///< export ids or sets ("petrel")
  bool figures = true;
  bool zip = false;
  std::string terms_of_use;                        ///< text approved by the model owner; empty = not supplied
  std::string contact;                             ///< overrides model.contact
  std::vector<std::string> assumptions, limitations;
  /// Regular expressions (ECMAScript) that must not appear in any text file or file name.
  std::vector<std::string> banned_patterns;
  /// Internal label -> external wording, applied to every text file before the scan.
  std::vector<std::pair<std::string, std::string>> translations;
  bool default_patterns = true;                    ///< also apply the built-in patterns below

  static PackageConfig load(const std::filesystem::path& file);
  /// Built-in banned patterns (Stratonix default): register identifiers (D-001, I-004, O-010,
  /// R-0.9.0, S-002), WBS codes, internal note identifiers (STX-TN-02, STX-RI-01), absolute
  /// local paths (C:\..., /home/...).
  static std::vector<std::string> builtin_patterns();
  /// Built-in translations of internal labels used by Stratonix's own texts.
  static std::vector<std::pair<std::string, std::string>> builtin_translations();
};

struct ScanHit { std::string file, pattern, excerpt; std::size_t line = 0; };

/// Scan every file name and every text file under root for the patterns.
std::vector<ScanHit> scan_folder(const std::filesystem::path& root, const std::vector<std::string>& patterns);

struct PackageOutcome {
  bool checks_passed = false, scan_passed = false;
  std::vector<std::string> input_status;  ///< "id: sharing (path)" lines
  std::vector<ScanHit> hits;
  std::filesystem::path folder, zip;
};

/// Build, check, export and package a model. Refuses (throws) when an input's sharing status is
/// unknown or a check fails. When the scan finds a banned pattern the folder is renamed to
/// "<folder>_REJECTED_SCAN", no zip is written and scan_passed is false.
PackageOutcome run_package(const std::filesystem::path& description, const std::filesystem::path& out_dir, const PackageConfig& cfg,
                           const BuildOptions& opt);

/// Plain-words explanation of an evidence class.
std::string evidence_in_plain_words(desc::Evidence e);

}  // namespace stratonix
