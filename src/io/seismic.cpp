// Stage 4C (STX-TN-16 section 2, decision 22): SEG-Y revisions 0 to 2.1 read for sampling and
// display; horizon picks and fault sticks from common interpretation exports, in their stated domain.
#include "stratonix/io/seismic.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/geometry.hpp"

namespace stratonix::io {

namespace {
// EBCDIC (code page 037) to ASCII for the printable range of textual headers
char ebcdic(unsigned char c) {
  static const char t[] =
      "                                                                "
      "           .<(+|&         !$*); -/         ,%_>?         `:#@'=\""
      " abcdefghi       jklmnopqr       ~stuvwxyz      ^         []    "
      "{ABCDEFGHI      }JKLMNOPQR      \\ STUVWXYZ      0123456789      ";
  return t[c];
}
double ibm_to_double(std::uint32_t x) {
  if ((x & 0x7FFFFFFF) == 0) return 0.0;
  const double s = (x >> 31) ? -1.0 : 1.0;
  const int e = static_cast<int>((x >> 24) & 0x7F) - 64;
  return s * static_cast<double>(x & 0xFFFFFF) / 16777216.0 * std::pow(16.0, e);
}
struct Bin {
  const std::vector<unsigned char>& b;
  bool le;
  std::uint64_t u(std::size_t off, int n) const {
    std::uint64_t v = 0;
    for (int i = 0; i < n; ++i) v |= std::uint64_t(b[off + static_cast<std::size_t>(le ? i : n - 1 - i)]) << (8 * i);
    return v;
  }
  std::int64_t s(std::size_t off, int n) const {
    const std::uint64_t v = u(off, n);
    const std::uint64_t sign = std::uint64_t(1) << (8 * n - 1);
    return static_cast<std::int64_t>((v ^ sign) - sign);
  }
};
int sample_bytes(int fmt) {
  switch (fmt) {
    case 1: case 2: case 4: case 5: case 10: return 4;
    case 3: case 11: return 2;
    case 6: case 9: case 12: return 8;
    case 7: case 15: return 3;
    case 8: case 16: return 1;
    default: return 0;
  }
}
}  // namespace

SegyVolume read_segy(const std::filesystem::path& file, const SegyLayout& layout) {
  std::ifstream in(file, std::ios::binary);
  if (!in) throw Error(ErrorCode::Io, "cannot open " + file.string());
  std::vector<unsigned char> b((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  if (b.size() < 3600) throw Error(ErrorCode::Parse, file.string() + ": shorter than the 3600-byte SEG-Y headers");
  SegyVolume v;
  v.source = file.string();
  // textual header: EBCDIC unless it reads as ASCII ('C' is 0x43 in ASCII, 0xC3 in EBCDIC)
  const bool ascii = b[0] == 'C' || (b[0] >= 0x20 && b[0] < 0x7F && b[1] >= 0x20 && b[1] < 0x7F);
  for (std::size_t i = 0; i < 3200; ++i) {
    const char c = ascii ? static_cast<char>(b[i]) : ebcdic(b[i]);
    v.textual_header += (i % 80 == 79) ? '\n' : (c >= 0x20 && c < 0x7F ? c : ' ');
  }
  v.report.push_back(std::string("textual header in ") + (ascii ? "ASCII" : "EBCDIC"));
  // byte order: revision 2 marker at bytes 3297-3300; otherwise big-endian
  {
    const std::uint32_t be = (std::uint32_t(b[3296]) << 24) | (std::uint32_t(b[3297]) << 16) | (std::uint32_t(b[3298]) << 8) | b[3299];
    v.little_endian = be == 0x04030201u;
    if (be == 0x02010403u || be == 0x03040102u) throw Error(ErrorCode::Unsupported, file.string() + ": pair-wise byte-swapped SEG-Y is not read");
    // without the marker (revisions 0 and 1, and writers that omit it): the sample format code
    // (bytes 3225-3226, 1 to 16) read both ways tells the byte order
    if (be != 0x01020304u && be != 0x04030201u) {
      const unsigned fb = (unsigned(b[3224]) << 8) | b[3225], fl = (unsigned(b[3225]) << 8) | b[3224];
      if ((fb < 1 || fb > 16) && fl >= 1 && fl <= 16) {
        v.little_endian = true;
        v.report.push_back("byte order: little-endian, from the sample format code (no revision 2 marker)");
      }
    }
  }
  Bin B{b, v.little_endian};
  const int rev = static_cast<int>(B.u(3200 + 300, 2));
  v.revision = rev == 0 ? 0 : (rev >> 8) != 0 ? (rev >> 8) : rev;  // 0x0100 -> 1, 0x0200 -> 2 (some writers put 1 or 2)
  const double dt_us = static_cast<double>(B.u(3200 + 16, 2));
  std::size_t ns = static_cast<std::size_t>(B.u(3200 + 20, 2));
  v.format = static_cast<int>(B.u(3200 + 24, 2));
  if (v.revision >= 2) {
    const std::uint64_t ns_ext = B.u(3200 + 268, 4);  // extended number of samples
    if (ns_ext) ns = static_cast<std::size_t>(ns_ext);
  }
  const int ext_text = static_cast<int>(B.s(3200 + 304, 2));
  const int nb = sample_bytes(v.format);
  if (!nb || v.format == 4 || v.format == 7 || v.format == 9 || v.format == 12 || v.format == 15)
    throw Error(ErrorCode::Unsupported, file.string() + ": sample format code " + std::to_string(v.format) + " is not read");
  if (!ns || dt_us <= 0) throw Error(ErrorCode::Parse, file.string() + ": the binary header gives no samples or no sample interval");
  v.samples = ns;
  v.dt = dt_us * 1e-6;
  v.domain = "time";
  std::size_t pos = 3600 + 3200 * static_cast<std::size_t>(std::max(0, ext_text));
  if (ext_text < 0) v.report.push_back("a variable number of extended textual headers is not read; traces assumed after the binary header");
  if (ext_text > 0) v.report.push_back(std::to_string(ext_text) + " extended textual headers skipped");
  const std::size_t il_off = layout.inline_byte ? static_cast<std::size_t>(layout.inline_byte - 1) : 188;
  const std::size_t xl_off = layout.crossline_byte ? static_cast<std::size_t>(layout.crossline_byte - 1) : 192;
  std::size_t nfallback = 0;
  while (pos + 240 <= b.size()) {
    SegyTrace t;
    t.il = static_cast<int>(B.s(pos + il_off, 4));
    t.xl = static_cast<int>(B.s(pos + xl_off, 4));
    double sc = static_cast<double>(B.s(pos + 70, 2));
    sc = sc == 0 ? 1.0 : sc < 0 ? -1.0 / sc : sc;
    double x, y;
    if (layout.x_byte && layout.y_byte) {
      x = static_cast<double>(B.s(pos + static_cast<std::size_t>(layout.x_byte - 1), 4));
      y = static_cast<double>(B.s(pos + static_cast<std::size_t>(layout.y_byte - 1), 4));
    } else {
      x = static_cast<double>(B.s(pos + 180, 4));  // CDP X and Y (bytes 181-188)
      y = static_cast<double>(B.s(pos + 184, 4));
      if (x == 0 && y == 0) {                      // source X and Y (bytes 73-80) when they are given
        const double sx = static_cast<double>(B.s(pos + 72, 4)), sy = static_cast<double>(B.s(pos + 76, 4));
        if (sx != 0 || sy != 0) { x = sx; y = sy; ++nfallback; }
      }
    }
    t.x = x * sc;
    t.y = y * sc;
    if (v.traces.empty()) v.t0 = static_cast<double>(B.s(pos + 108, 2)) * 1e-3;  // delay recording time (ms)
    const std::size_t tns = B.u(pos + 114, 2) ? static_cast<std::size_t>(B.u(pos + 114, 2)) : ns;
    if (tns != ns && v.report.size() < 20) v.report.push_back("trace " + std::to_string(v.traces.size() + 1) + " has " + std::to_string(tns) + " samples; the file's " + std::to_string(ns) + " used");
    pos += 240;
    if (pos + ns * static_cast<std::size_t>(nb) > b.size()) {
      v.report.push_back("the last trace is truncated: not read");
      break;
    }
    t.s.resize(ns);
    for (std::size_t k = 0; k < ns; ++k) {
      const std::size_t o = pos + k * static_cast<std::size_t>(nb);
      double a = 0;
      switch (v.format) {
        case 1: a = ibm_to_double(static_cast<std::uint32_t>(B.u(o, 4))); break;
        case 2: a = static_cast<double>(B.s(o, 4)); break;
        case 3: a = static_cast<double>(B.s(o, 2)); break;
        case 5: { std::uint32_t w = static_cast<std::uint32_t>(B.u(o, 4)); float f; std::memcpy(&f, &w, 4); a = f; break; }
        case 6: { std::uint64_t w = B.u(o, 8); double d; std::memcpy(&d, &w, 8); a = d; break; }
        case 8: a = static_cast<double>(B.s(o, 1)); break;
        case 10: a = static_cast<double>(B.u(o, 4)); break;
        case 11: a = static_cast<double>(B.u(o, 2)); break;
        case 16: a = static_cast<double>(B.u(o, 1)); break;
      }
      t.s[k] = static_cast<float>(a);
    }
    pos += ns * static_cast<std::size_t>(nb);
    v.traces.push_back(std::move(t));
  }
  if (nfallback) v.report.push_back(std::to_string(nfallback) + " traces have no CDP coordinates; source coordinates (bytes 73-80) used");
  v.report.push_back("SEG-Y revision " + std::to_string(v.revision) + ", format " + std::to_string(v.format) + (v.little_endian ? ", little-endian" : ", big-endian") +
                     ", " + std::to_string(v.traces.size()) + " traces of " + std::to_string(ns) + " samples at " + std::to_string(dt_us) + " us");
  return v;
}

double SegyVolume::sample(double x, double y, double t) const {
  if (traces.empty() || dt <= 0) return kNaN;
  // the survey grid: affine map from (inline, crossline) to (x, y) by least squares, fitted once
  if (affine_.empty()) {
    double S[3][3] = {{0}}, Rx[3] = {0}, Ry[3] = {0};
    for (std::size_t i = 0; i < traces.size(); ++i) {
      at_[{traces[i].il, traces[i].xl}] = i;
      const double r[3] = {1, static_cast<double>(traces[i].il), static_cast<double>(traces[i].xl)};
      for (int a = 0; a < 3; ++a) {
        for (int c = 0; c < 3; ++c) S[a][c] += r[a] * r[c];
        Rx[a] += r[a] * traces[i].x;
        Ry[a] += r[a] * traces[i].y;
      }
    }
    auto solve = [&](double M[3][3], double r[3], double out[3]) {
      double A[3][4];
      for (int a = 0; a < 3; ++a) { for (int c = 0; c < 3; ++c) A[a][c] = M[a][c]; A[a][3] = r[a]; }
      for (int p = 0; p < 3; ++p) {
        int piv = p;
        for (int q = p + 1; q < 3; ++q) if (std::fabs(A[q][p]) > std::fabs(A[piv][p])) piv = q;
        std::swap(A[p], A[piv]);
        if (std::fabs(A[p][p]) < 1e-12) return false;
        for (int q = 0; q < 3; ++q) {
          if (q == p) continue;
          const double f = A[q][p] / A[p][p];
          for (int c = p; c < 4; ++c) A[q][c] -= f * A[p][c];
        }
      }
      for (int a = 0; a < 3; ++a) out[a] = A[a][3] / A[a][a];
      return true;
    };
    double ax[3], ay[3];
    if (!solve(S, Rx, ax) || !solve(S, Ry, ay)) { affine_ = {kNaN}; return kNaN; }  // a line or a single trace
    affine_ = {ax[0], ax[1], ax[2], ay[0], ay[1], ay[2]};
  }
  if (affine_.size() < 6) return kNaN;
  const double ax[3] = {affine_[0], affine_[1], affine_[2]}, ay[3] = {affine_[3], affine_[4], affine_[5]};
  const auto& at = at_;
  const double det = ax[1] * ay[2] - ax[2] * ay[1];
  if (std::fabs(det) < 1e-12) return kNaN;
  const double dx = x - ax[0], dy = y - ay[0];
  const double fi = (dx * ay[2] - dy * ax[2]) / det, fj = (ax[1] * dy - ay[1] * dx) / det;
  const int i0 = static_cast<int>(std::floor(fi + 1e-9)), j0 = static_cast<int>(std::floor(fj + 1e-9));
  const double u = std::clamp(fi - i0, 0.0, 1.0), w = std::clamp(fj - j0, 0.0, 1.0);
  const double kf = (t - t0) / dt;
  if (kf < 0 || kf > static_cast<double>(samples - 1)) return kNaN;
  const std::size_t k0 = std::min(static_cast<std::size_t>(kf), samples - 1), k1 = std::min(k0 + 1, samples - 1);
  const double tk = kf - static_cast<double>(k0);
  auto amp = [&](int i, int j, double& out) {
    auto it = at.find({i, j});
    if (it == at.end()) return false;
    const auto& s = traces[it->second].s;
    out = (1 - tk) * s[k0] + tk * s[k1];
    return true;
  };
  double a00, a10, a01, a11;
  const bool e00 = amp(i0, j0, a00);
  const bool e10 = u < 1e-9 ? (a10 = 0, true) : amp(i0 + 1, j0, a10);
  const bool e01 = w < 1e-9 ? (a01 = 0, true) : amp(i0, j0 + 1, a01);
  const bool e11 = (u < 1e-9 || w < 1e-9) ? (a11 = 0, true) : amp(i0 + 1, j0 + 1, a11);
  if (!e00 || !e10 || !e01 || !e11) return kNaN;
  return (1 - u) * (1 - w) * a00 + u * (1 - w) * a10 + (1 - u) * w * a01 + u * w * a11;
}

// ---- picks and sticks -------------------------------------------------------------------------
namespace {
std::vector<std::string> read_lines(const std::filesystem::path& f) {
  std::ifstream in(f, std::ios::binary);
  if (!in) throw Error(ErrorCode::Io, "cannot open " + f.string());
  std::vector<std::string> out;
  std::string l;
  while (std::getline(in, l)) {
    if (!l.empty() && l.back() == '\r') l.pop_back();
    out.push_back(l);
  }
  return out;
}
std::vector<std::string> tokens(const std::string& l, bool comma) {
  std::vector<std::string> t;
  std::string cur;
  for (char c : l) {
    if ((comma && c == ',') || (!comma && std::isspace(static_cast<unsigned char>(c)))) {
      if (comma || !cur.empty()) t.push_back(cur);
      cur.clear();
    } else if (!(comma && std::isspace(static_cast<unsigned char>(c)) && cur.empty())) {
      cur += c;
    }
  }
  if (!cur.empty() || comma) t.push_back(cur);
  for (auto& x : t)
    while (!x.empty() && std::isspace(static_cast<unsigned char>(x.back()))) x.pop_back();
  return t;
}
bool num(const std::string& s, double& v) {
  char* e = nullptr;
  v = std::strtod(s.c_str(), &e);
  return !s.empty() && e && *e == '\0';
}
void need_domain(const std::string& domain, const std::filesystem::path& f) {
  if (domain != "time" && domain != "depth")
    throw Error(ErrorCode::Validation, f.string() + ": the domain (time or depth) is not stated; it is never assumed (rule F45)");
}
}  // namespace

std::vector<PickSet> read_picks(const std::filesystem::path& file, const std::string& format, const std::string& domain, const std::string& unit) {
  need_domain(domain, file);
  std::map<std::string, PickSet> sets;
  std::vector<std::string> order;
  auto set_of = [&](const std::string& name) -> PickSet& {
    if (!sets.count(name)) { sets[name] = PickSet{name, file.string(), domain, unit, {}}; order.push_back(name); }
    return sets[name];
  };
  const auto lines = read_lines(file);
  const std::string stem = file.stem().string();
  if (format == "charisma") {
    // INLINE : il XLINE : xl x y z
    for (auto& l : lines) {
      auto t = tokens(l, false);
      std::vector<double> v;
      for (auto& s : t) { double x; if (num(s, x)) v.push_back(x); }
      if (l.find("INLINE") == std::string::npos || v.size() < 5) continue;
      Pick p;
      p.il = static_cast<int>(v[0]);
      p.xl = static_cast<int>(v[1]);
      p.p = {v[2], v[3], v[4]};
      set_of(stem).picks.push_back(p);
    }
  } else if (format == "csv") {
    if (lines.empty()) throw Error(ErrorCode::Parse, file.string() + ": empty");
    auto h = tokens(lines[0], true);
    int cx = -1, cy = -1, cz = -1, ci = -1, cj = -1, cn = -1;
    for (std::size_t i = 0; i < h.size(); ++i) {
      std::string n = h[i];
      for (auto& c : n) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      if (n == "x") cx = static_cast<int>(i);
      else if (n == "y") cy = static_cast<int>(i);
      else if (n == "z" || n == "time" || n == "twt" || n == "depth") cz = static_cast<int>(i);
      else if (n == "inline" || n == "il") ci = static_cast<int>(i);
      else if (n == "crossline" || n == "xl") cj = static_cast<int>(i);
      else if (n == "horizon" || n == "surface" || n == "name") cn = static_cast<int>(i);
    }
    if (cx < 0 || cy < 0 || cz < 0) throw Error(ErrorCode::Parse, file.string() + ": a csv of picks needs x, y and z (or time, twt, depth) columns");
    for (std::size_t r = 1; r < lines.size(); ++r) {
      auto t = tokens(lines[r], true);
      double x, y, z;
      if (t.size() <= static_cast<std::size_t>(std::max({cx, cy, cz})) || !num(t[static_cast<std::size_t>(cx)], x) || !num(t[static_cast<std::size_t>(cy)], y) ||
          !num(t[static_cast<std::size_t>(cz)], z))
        continue;
      Pick p;
      p.p = {x, y, z};
      double a;
      if (ci >= 0 && static_cast<std::size_t>(ci) < t.size() && num(t[static_cast<std::size_t>(ci)], a)) p.il = static_cast<int>(a);
      if (cj >= 0 && static_cast<std::size_t>(cj) < t.size() && num(t[static_cast<std::size_t>(cj)], a)) p.xl = static_cast<int>(a);
      set_of(cn >= 0 && static_cast<std::size_t>(cn) < t.size() ? t[static_cast<std::size_t>(cn)] : stem).picks.push_back(p);
    }
  } else if (format == "iesx" || format == "petrel_ascii" || format == "xyz") {
    for (auto& l : lines) {
      auto t = tokens(l, false);
      double x, y, z;
      if (t.size() < 3 || !num(t[0], x) || !num(t[1], y) || !num(t[2], z)) continue;
      Pick p;
      p.p = {x, y, z};
      set_of(stem).picks.push_back(p);
    }
  } else {
    throw Error(ErrorCode::Unsupported, "pick format '" + format + "' (charisma, csv, iesx, petrel_ascii, xyz)");
  }
  std::vector<PickSet> out;
  for (auto& n : order) out.push_back(sets[n]);
  return out;
}

StickSet read_fault_sticks(const std::filesystem::path& file, const std::string& format, const std::string& domain, const std::string& unit) {
  need_domain(domain, file);
  StickSet s{file.string(), domain, unit, {}};
  std::map<std::pair<std::string, int>, std::size_t> idx;
  auto add = [&](const std::string& fault, int id, const Vec3& p) {
    auto k = std::make_pair(fault, id);
    if (!idx.count(k)) { idx[k] = s.sticks.size(); s.sticks.push_back({fault, id, {}}); }
    s.sticks[idx[k]].p.push_back(p);
  };
  const auto lines = read_lines(file);
  const std::string stem = file.stem().string();
  if (format == "csv") {
    if (lines.empty()) throw Error(ErrorCode::Parse, file.string() + ": empty");
    auto h = tokens(lines[0], true);
    int c[5] = {-1, -1, -1, -1, -1};
    for (std::size_t i = 0; i < h.size(); ++i) {
      std::string n = h[i];
      for (auto& ch : n) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
      if (n == "x") c[0] = static_cast<int>(i);
      else if (n == "y") c[1] = static_cast<int>(i);
      else if (n == "z" || n == "time" || n == "twt" || n == "depth") c[2] = static_cast<int>(i);
      else if (n == "fault" || n == "name") c[3] = static_cast<int>(i);
      else if (n == "stick" || n == "stick_id" || n == "id") c[4] = static_cast<int>(i);
    }
    if (c[0] < 0 || c[1] < 0 || c[2] < 0 || c[4] < 0) throw Error(ErrorCode::Parse, file.string() + ": fault sticks need x, y, z and stick columns");
    for (std::size_t r = 1; r < lines.size(); ++r) {
      auto t = tokens(lines[r], true);
      double v[3], id;
      bool ok = true;
      for (int q = 0; q < 3; ++q) ok = ok && static_cast<std::size_t>(c[q]) < t.size() && num(t[static_cast<std::size_t>(c[q])], v[q]);
      ok = ok && static_cast<std::size_t>(c[4]) < t.size() && num(t[static_cast<std::size_t>(c[4])], id);
      if (!ok) continue;
      add(c[3] >= 0 && static_cast<std::size_t>(c[3]) < t.size() ? t[static_cast<std::size_t>(c[3])] : stem, static_cast<int>(id), {v[0], v[1], v[2]});
    }
  } else if (format == "charisma") {
    // ... x y z fault_name stick_id  (the last two fields)
    for (auto& l : lines) {
      auto t = tokens(l, false);
      if (t.size() < 5) continue;
      double x, y, z, id;
      const std::size_t n = t.size();
      if (!num(t[n - 1], id) || !num(t[n - 3], z) || !num(t[n - 4], y) || !num(t[n - 5], x)) continue;
      add(t[n - 2], static_cast<int>(id), {x, y, z});
    }
  } else if (format == "petrel_ascii") {
    // x y z stick_id, one fault per file (the file name)
    for (auto& l : lines) {
      auto t = tokens(l, false);
      double x, y, z, id;
      if (t.size() < 4 || !num(t[0], x) || !num(t[1], y) || !num(t[2], z) || !num(t[3], id)) continue;
      add(stem, static_cast<int>(id), {x, y, z});
    }
  } else {
    throw Error(ErrorCode::Unsupported, "fault stick format '" + format + "' (csv, charisma, petrel_ascii)");
  }
  return s;
}

}  // namespace stratonix::io
