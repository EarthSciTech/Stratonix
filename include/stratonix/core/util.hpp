/// @file util.hpp
/// @brief SHA-256, locale-independent number formatting, files, threads, version, logging.
#pragma once
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace stratonix {

// ---- SHA-256 -------------------------------------------------------------------------------
/// Hex SHA-256 of a byte string.
std::string sha256_hex(std::string_view data);
/// Hex SHA-256 of a file's contents.
std::string sha256_file(const std::filesystem::path& p);

// ---- Number formatting (independent of locale and platform) --------------------------------
/// Fixed notation with `decimals` digits after the point ("-0.000" is normalised to "0.000").
std::string fmt_fixed(double v, int decimals);
/// Shortest round-trip representation.
std::string fmt_shortest(double v);
/// Parse a double with '.' as decimal separator; throws Error on failure.
double parse_double(std::string_view s);

// ---- Files ---------------------------------------------------------------------------------
std::string read_file(const std::filesystem::path& p);
/// Write bytes exactly (binary mode, so '\n' stays '\n' on Windows). Creates parent folders.
void write_file(const std::filesystem::path& p, std::string_view data);
/// Path with forward slashes, relative to base if possible (used in manifests).
std::string portable_path(const std::filesystem::path& p, const std::filesystem::path& base = {});

// ---- Threads -------------------------------------------------------------------------------
/// Set the number of worker threads (0 = hardware concurrency). Results never depend on it.
void set_thread_count(unsigned n);
unsigned thread_count();
/// Run fn(begin, end) over [0, n) in fixed chunks of `chunk` items. The chunking does not depend
/// on the thread count, and each index is processed by exactly one call, so results are identical
/// for any number of threads. Exceptions from workers are rethrown on the caller.
void parallel_for(std::size_t n, std::size_t chunk,
                  const std::function<void(std::size_t, std::size_t)>& fn);
/// Global switch used by tests: when false, optimised kernels fall back to reference versions.
void set_optimised_kernels(bool on);
bool optimised_kernels();

/// Peak resident memory of this process in bytes (0 if unknown on this platform).
std::size_t peak_rss_bytes();

// ---- Version and build information ---------------------------------------------------------
struct BuildInfo {
  std::string version, git_commit, compiler, platform, build_type;
};
const BuildInfo& build_info();

// ---- Logging -------------------------------------------------------------------------------
enum class LogLevel { Debug, Info, Warn, Error };
void set_log_level(LogLevel l);
LogLevel log_level();
void log(LogLevel l, std::string_view msg);

// ---- Strings -------------------------------------------------------------------------------
std::string trim(std::string_view s);
std::vector<std::string> split(std::string_view s, char sep);
bool starts_with(std::string_view s, std::string_view p);
std::string to_lower(std::string_view s);
/// Escape for XML/SVG text and attributes.
std::string xml_escape(std::string_view s);
/// Escape for JSON strings (without quotes).
std::string json_escape(std::string_view s);

}  // namespace stratonix
