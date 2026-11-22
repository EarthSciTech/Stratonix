/// @file seismic.hpp
/// @brief Stage 4C (STX-TN-16 section 2): SEG-Y volumes read for sampling and display, horizon
/// picks and fault sticks from common interpretation exports, each in its stated domain.
#pragma once
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"

namespace stratonix::io {

/// Trace-header byte positions (1-based, as in the standard); 0 = from the file's revision.
struct SegyLayout { int inline_byte = 0, crossline_byte = 0, x_byte = 0, y_byte = 0; };
struct SegyTrace { int il = 0, xl = 0; double x = 0, y = 0; std::vector<float> s; };
struct SegyVolume {
  std::string source, textual_header;
  int revision = 0;                 ///< 0, 1 or 2 (major) from the binary header
  int format = 0;                   ///< sample format code (1 IBM, 5 IEEE, 2/3/8 integers ...)
  bool little_endian = false;
  double dt = 0;                    ///< sample interval (s or m, per `domain`)
  double t0 = 0;
  std::size_t samples = 0;
  std::string domain;               ///< "time" or "depth" from the file or the description
  std::vector<SegyTrace> traces;
  std::vector<std::string> report;
  /// Amplitude at a point by bilinear interpolation between traces and linear in time; NaN outside.
  double sample(double x, double y, double t) const;
};
SegyVolume read_segy(const std::filesystem::path& file, const SegyLayout& layout = {});

/// Picks and sticks; `domain` is "time" (two-way, ms unless the unit says otherwise) or "depth".
struct Pick { int il = 0, xl = 0; Vec3 p; };  ///< p.z holds the time or depth value as read
struct PickSet { std::string name, source, domain, unit; std::vector<Pick> picks; };
struct FaultStick { std::string fault; int id = 0; std::vector<Vec3> p; };
struct StickSet { std::string source, domain, unit; std::vector<FaultStick> sticks; };
/// Formats: "charisma", "iesx", "petrel_ascii", "csv". Refused when `domain` is empty.
std::vector<PickSet> read_picks(const std::filesystem::path& file, const std::string& format, const std::string& domain,
                                const std::string& unit);
StickSet read_fault_sticks(const std::filesystem::path& file, const std::string& format, const std::string& domain, const std::string& unit);

}  // namespace stratonix::io
