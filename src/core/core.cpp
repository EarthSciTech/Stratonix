// Core utilities: mesh quality, units, SHA-256, formatting, files, threads, build info, logging.
#include <atomic>
#include <cctype>
#include <charconv>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <thread>
#include <utility>

#include "stratonix/core/crs.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/mesh.hpp"
#include "stratonix/core/units.hpp"
#include "stratonix/core/util.hpp"

#ifndef STRATONIX_VERSION
#define STRATONIX_VERSION "0.0.0"
#endif
#ifndef STRATONIX_GIT_COMMIT
#define STRATONIX_GIT_COMMIT "unknown"
#endif
#ifndef STRATONIX_BUILD_TYPE
#define STRATONIX_BUILD_TYPE "unknown"
#endif

namespace stratonix {

// ---- mesh ----------------------------------------------------------------------------------
double TriMesh::area() const {
  double a = 0;
  for (auto& tr : t) a += 0.5 * norm(cross(v[tr[1]] - v[tr[0]], v[tr[2]] - v[tr[0]]));
  return a;
}

MeshQuality mesh_quality(const TriMesh& m, double area_eps) {
  MeshQuality q;
  q.triangles = m.t.size();
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> edges;
  bool first = true;
  for (auto& tr : m.t) {
    double a = 0.5 * norm(cross(m.v[tr[1]] - m.v[tr[0]], m.v[tr[2]] - m.v[tr[0]]));
    if (tr[0] == tr[1] || tr[1] == tr[2] || tr[0] == tr[2] || a < area_eps) ++q.degenerate;
    if (first) { q.min_area = q.max_area = a; first = false; }
    q.min_area = std::min(q.min_area, a);
    q.max_area = std::max(q.max_area, a);
    for (int k = 0; k < 3; ++k) {
      std::uint32_t a0 = tr[k], b0 = tr[(k + 1) % 3];
      edges[{std::min(a0, b0), std::max(a0, b0)}]++;
    }
  }
  for (auto& [e, c] : edges) {
    if (c == 1) ++q.boundary_edges;
    if (c > 2) ++q.nonmanifold_edges;
  }
  return q;
}

// ---- units ---------------------------------------------------------------------------------
std::optional<UnitDef> find_unit(std::string_view name) {
  static const UnitDef table[] = {
      {"m", Dimension::Length, 1.0},          {"km", Dimension::Length, 1000.0},
      {"ft", Dimension::Length, 0.3048},      {"deg", Dimension::Angle, kPi / 180.0},
      {"rad", Dimension::Angle, 1.0},         {"fraction", Dimension::Dimensionless, 1.0},
      {"ratio", Dimension::Dimensionless, 1.0}, {"percent", Dimension::Dimensionless, 0.01},
      {"ordinal", Dimension::Ordinal, 1.0},   {"flag", Dimension::Flag, 1.0},
  };
  for (auto& u : table)
    if (u.name == name) return u;
  return std::nullopt;
}
const char* dimension_name(Dimension d) {
  switch (d) {
    case Dimension::Length: return "length";
    case Dimension::Angle: return "angle";
    case Dimension::Dimensionless: return "dimensionless";
    case Dimension::Ordinal: return "ordinal";
    case Dimension::Flag: return "flag";
    default: return "any";
  }
}

bool crs_reprojection_available() {
#ifdef STRATONIX_WITH_PROJ
  return true;
#else
  return false;
#endif
}

// ---- SHA-256 (FIPS 180-4) -------------------------------------------------------------------
namespace {
struct Sha256 {
  std::uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
  unsigned char buf[64];
  std::size_t blen = 0;
  std::uint64_t total = 0;
  static std::uint32_t rotr(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
  void block(const unsigned char* p) {
    static const std::uint32_t k[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2};
    std::uint32_t w[64];
    for (int i = 0; i < 16; ++i)
      w[i] = (std::uint32_t(p[4 * i]) << 24) | (std::uint32_t(p[4 * i + 1]) << 16) |
             (std::uint32_t(p[4 * i + 2]) << 8) | std::uint32_t(p[4 * i + 3]);
    for (int i = 16; i < 64; ++i) {
      std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
      std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    std::uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
    for (int i = 0; i < 64; ++i) {
      std::uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
      std::uint32_t ch = (e & f) ^ (~e & g);
      std::uint32_t t1 = hh + S1 + ch + k[i] + w[i];
      std::uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
      std::uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
      std::uint32_t t2 = S0 + mj;
      hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
  }
  void update(const unsigned char* p, std::size_t n) {
    total += n;
    while (n > 0) {
      std::size_t take = std::min<std::size_t>(64 - blen, n);
      std::copy(p, p + take, buf + blen);
      blen += take; p += take; n -= take;
      if (blen == 64) { block(buf); blen = 0; }
    }
  }
  std::string hex() {
    std::uint64_t bits = total * 8;
    unsigned char pad = 0x80;
    update(&pad, 1);
    unsigned char z = 0;
    while (blen != 56) update(&z, 1);
    unsigned char len[8];
    for (int i = 0; i < 8; ++i) len[i] = static_cast<unsigned char>(bits >> (56 - 8 * i));
    update(len, 8);
    static const char* hx = "0123456789abcdef";
    std::string out;
    for (auto v : h)
      for (int s = 28; s >= 0; s -= 4) out += hx[(v >> s) & 0xF];
    return out;
  }
};
}  // namespace

std::string sha256_hex(std::string_view data) {
  Sha256 s;
  s.update(reinterpret_cast<const unsigned char*>(data.data()), data.size());
  return s.hex();
}
std::string sha256_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  if (!f) throw Error(ErrorCode::Io, "cannot open " + p.string());
  Sha256 s;
  std::vector<char> b(1 << 16);
  while (f) {
    f.read(b.data(), static_cast<std::streamsize>(b.size()));
    s.update(reinterpret_cast<const unsigned char*>(b.data()), static_cast<std::size_t>(f.gcount()));
  }
  return s.hex();
}

// ---- formatting ----------------------------------------------------------------------------
std::string fmt_fixed(double v, int decimals) {
  if (std::isnan(v)) return "nan";
  char buf[64];
  auto r = std::to_chars(buf, buf + sizeof buf, v, std::chars_format::fixed, decimals);
  std::string s(buf, r.ptr);
  if (!s.empty() && s[0] == '-') {
    bool zero = true;
    for (char c : s.substr(1))
      if (c != '0' && c != '.') zero = false;
    if (zero) s.erase(0, 1);
  }
  return s;
}
std::string fmt_shortest(double v) {
  if (std::isnan(v)) return "nan";
  char buf[64];
  auto r = std::to_chars(buf, buf + sizeof buf, v);
  return std::string(buf, r.ptr);
}
double parse_double(std::string_view s) {
  std::string t = trim(s);
  if (!t.empty() && t[0] == '+') t.erase(0, 1);
  double v = 0;
  auto r = std::from_chars(t.data(), t.data() + t.size(), v);
  if (r.ec != std::errc() || r.ptr != t.data() + t.size())
    throw Error(ErrorCode::Parse, "not a number: '" + std::string(s) + "'");
  return v;
}

// ---- files ---------------------------------------------------------------------------------
std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  if (!f) throw Error(ErrorCode::Io, "cannot open " + p.string());
  std::ostringstream ss;
  ss << f.rdbuf();
  return ss.str();
}
void write_file(const std::filesystem::path& p, std::string_view data) {
  if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());
  std::ofstream f(p, std::ios::binary | std::ios::trunc);
  if (!f) throw Error(ErrorCode::Io, "cannot write " + p.string());
  f.write(data.data(), static_cast<std::streamsize>(data.size()));
  if (!f) throw Error(ErrorCode::Io, "write failed " + p.string());
}
std::string portable_path(const std::filesystem::path& p, const std::filesystem::path& base) {
  std::filesystem::path q = p;
  if (!base.empty()) {
    std::error_code ec;
    auto r = std::filesystem::relative(p, base, ec);
    if (!ec && !r.empty()) q = r;
  }
  return q.generic_string();
}

// ---- threads -------------------------------------------------------------------------------
namespace {
std::atomic<unsigned> g_threads{0};
std::atomic<bool> g_optimised{true};
}  // namespace
void set_thread_count(unsigned n) { g_threads = n; }
unsigned thread_count() {
  unsigned n = g_threads.load();
  if (n == 0) n = std::max(1u, std::thread::hardware_concurrency());
  return n;
}
void set_optimised_kernels(bool on) { g_optimised = on; }
bool optimised_kernels() { return g_optimised.load(); }

void parallel_for(std::size_t n, std::size_t chunk,
                  const std::function<void(std::size_t, std::size_t)>& fn) {
  if (n == 0) return;
  if (chunk == 0) chunk = 1;
  std::size_t nchunks = (n + chunk - 1) / chunk;
  unsigned nt = static_cast<unsigned>(std::min<std::size_t>(thread_count(), nchunks));
  if (nt <= 1) {
    for (std::size_t c = 0; c < nchunks; ++c) fn(c * chunk, std::min(n, (c + 1) * chunk));
    return;
  }
  std::atomic<std::size_t> next{0};
  std::exception_ptr err;
  std::mutex err_m;
  auto worker = [&] {
    for (;;) {
      std::size_t c = next.fetch_add(1);
      if (c >= nchunks) return;
      try {
        fn(c * chunk, std::min(n, (c + 1) * chunk));
      } catch (...) {
        std::lock_guard<std::mutex> lk(err_m);
        if (!err) err = std::current_exception();
      }
    }
  };
  std::vector<std::thread> ts;
  ts.reserve(nt - 1);
  for (unsigned i = 1; i < nt; ++i) ts.emplace_back(worker);
  worker();
  for (auto& t : ts) t.join();
  if (err) std::rethrow_exception(err);
}

// ---- build info ----------------------------------------------------------------------------
const BuildInfo& build_info() {
  static const BuildInfo bi = [] {
    BuildInfo b;
    b.version = STRATONIX_VERSION;
    b.git_commit = STRATONIX_GIT_COMMIT;
    b.build_type = STRATONIX_BUILD_TYPE;
#if defined(_MSC_VER)
    b.compiler = "MSVC " + std::to_string(_MSC_VER);
#elif defined(__clang__)
    b.compiler = std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
    b.compiler = "GCC " + std::to_string(__GNUC__) + "." + std::to_string(__GNUC_MINOR__) + "." +
                 std::to_string(__GNUC_PATCHLEVEL__);
#else
    b.compiler = "unknown";
#endif
#if defined(_WIN32)
    b.platform = "windows-x64";
#elif defined(__linux__)
    b.platform = "linux-x64";
#elif defined(__APPLE__)
    b.platform = "macos";
#else
    b.platform = "unknown";
#endif
    return b;
  }();
  return bi;
}

// ---- logging -------------------------------------------------------------------------------
namespace {
std::atomic<int> g_level{static_cast<int>(LogLevel::Info)};
std::mutex g_log_m;
}  // namespace
void set_log_level(LogLevel l) { g_level = static_cast<int>(l); }
void log(LogLevel l, std::string_view msg) {
  if (static_cast<int>(l) < g_level.load()) return;
  static const char* tag[] = {"debug", "info", "warning", "error"};
  std::lock_guard<std::mutex> lk(g_log_m);
  std::cerr << "[stratonix " << tag[static_cast<int>(l)] << "] " << msg << "\n";
}

// ---- strings -------------------------------------------------------------------------------
std::string trim(std::string_view s) {
  std::size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return std::string(s.substr(a, b - a));
}
std::vector<std::string> split(std::string_view s, char sep) {
  std::vector<std::string> out;
  std::size_t a = 0;
  for (std::size_t i = 0; i <= s.size(); ++i)
    if (i == s.size() || s[i] == sep) { out.emplace_back(s.substr(a, i - a)); a = i + 1; }
  return out;
}
bool starts_with(std::string_view s, std::string_view p) { return s.substr(0, p.size()) == p; }
std::string to_lower(std::string_view s) {
  std::string r(s);
  for (auto& c : r) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return r;
}
std::string xml_escape(std::string_view s) {
  std::string r;
  for (char c : s) {
    switch (c) {
      case '&': r += "&amp;"; break;
      case '<': r += "&lt;"; break;
      case '>': r += "&gt;"; break;
      case '"': r += "&quot;"; break;
      case '\'': r += "&apos;"; break;
      default: r += c;
    }
  }
  return r;
}
std::string json_escape(std::string_view s) {
  std::string r;
  for (char c : s) {
    switch (c) {
      case '"': r += "\\\""; break;
      case '\\': r += "\\\\"; break;
      case '\n': r += "\\n"; break;
      case '\r': r += "\\r"; break;
      case '\t': r += "\\t"; break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char b[8];
          std::snprintf(b, sizeof b, "\\u%04x", c);
          r += b;
        } else {
          r += c;
        }
    }
  }
  return r;
}

}  // namespace stratonix
