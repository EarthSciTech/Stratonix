/// @file welldata.hpp
/// @brief Stage 4C (STX-TN-16 section 1): well logs read from LAS 2.0, LAS 3.0 and DLIS (RP66 V1)
/// files as delivered. Units and null values come from the file; every read reports them.
#pragma once
#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace stratonix::io {

struct LogCurve {
  std::string mnemonic, unit, description;
  double null_value = -999.25;       ///< the file's null; samples equal to it are NaN in `values`
  std::size_t nulls = 0;             ///< samples at the null value
  std::vector<double> values;        ///< one per index sample (NaN at nulls)
  std::vector<std::string> text;     ///< string columns (LAS 3.0), empty otherwise
};

/// One LAS 3.0 section group or the parameters of a LAS 2.0 file: mnemonic -> (unit, value, description).
struct LasParameter { std::string unit, value, description; };

struct WellLogSet {
  std::string source;                ///< the file, and for DLIS the logical file and frame
  std::string version;               ///< "2.0", "3.0" or "RP66 V1"
  std::string well;                  ///< well name from the file (WELL, or ORIGIN for DLIS)
  std::string index_type;            ///< "MD", "TVD", "TIME" or "" when the file does not say
  LogCurve index;                    ///< depth or time index
  std::vector<LogCurve> curves;
  std::map<std::string, LasParameter> well_info, parameters;
  /// LAS 3.0 tops and inclinometry sections when present (curves of their data sections)
  std::vector<LogCurve> tops, inclinometry;
  std::vector<std::string> report;   ///< what was read and what was not, with line or byte offsets
  const LogCurve* find(const std::string& mnemonic) const;
};

/// LAS 2.0 (wrapped or not) and LAS 3.0 (log, tops and inclinometry groups).
WellLogSet read_las(const std::filesystem::path& file);
WellLogSet parse_las(const std::string& text, const std::string& source);
/// DLIS (RP66 V1, decision 21): one set per logical file and frame.
std::vector<WellLogSet> read_dlis(const std::filesystem::path& file);

}  // namespace stratonix::io
