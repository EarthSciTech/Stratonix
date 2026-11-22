// Stage 4C (STX-TN-16 sections 1.1 and 1.2): well logs from LAS 2.0, LAS 3.0 and DLIS (RP66 V1).
// Files are read as delivered and never rewritten; what cannot be read is reported with its line or
// byte offset, never guessed.
#include "stratonix/io/welldata.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/geometry.hpp"

namespace stratonix::io {

namespace {

std::string trim(const std::string& s) {
  std::size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}
std::string upper(std::string s) {
  for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}
bool parse_number(const std::string& t, double& v) {
  const std::string s = trim(t);
  if (s.empty()) return false;
  char* end = nullptr;
  v = std::strtod(s.c_str(), &end);
  return end && *end == '\0';
}
std::string index_type_of(const std::string& mnem, const std::string& unit) {
  const std::string m = upper(mnem), u = upper(unit);
  if (u == "S" || u == "MS" || u == "SEC" || u == "US" || m == "TIME" || m == "ETIM") return "TIME";
  if (m.find("TVD") != std::string::npos) return "TVD";
  if (m == "DEPT" || m == "DEPTH" || m == "MD" || m == "TDEP") return "MD";
  return "";
}
std::string fmt(double v) {
  std::ostringstream o;
  o.precision(10);
  o << v;
  return o.str();
}

// ---- LAS ------------------------------------------------------------------------------------
struct HeaderLine { std::string mnem, unit, value, desc, format, assoc; };
/// MNEM.UNIT  VALUE : DESCRIPTION {format} | association (CWLS LAS 2.0 and 3.0). The mnemonic ends
/// at the first period, the unit at the first space after it, the value at the last colon.
bool parse_header_line(const std::string& line, bool las3, HeaderLine& h) {
  const std::size_t dot = line.find('.');
  if (dot == std::string::npos) return false;
  std::string rest = line.substr(dot + 1);
  std::size_t colon = rest.rfind(':');
  if (las3) {
    // in LAS 3.0 the description may be followed by {format} and | association, which may not
    // contain colons; the value ends at the last colon before them
    const std::size_t brace = rest.find('{'), bar = rest.find('|');
    const std::size_t lim = std::min(brace, bar);
    if (lim != std::string::npos) colon = rest.substr(0, lim).rfind(':');
    if (brace != std::string::npos) {
      const std::size_t e = rest.find('}', brace);
      h.format = trim(rest.substr(brace + 1, e == std::string::npos ? std::string::npos : e - brace - 1));
    }
    if (bar != std::string::npos) h.assoc = trim(rest.substr(bar + 1));
  }
  if (colon == std::string::npos) colon = rest.size();  // "MNEM.UNIT" without a colon: no value, no description
  h.mnem = trim(line.substr(0, dot));
  std::size_t sp = 0;
  while (sp < colon && !std::isspace(static_cast<unsigned char>(rest[sp]))) ++sp;
  h.unit = rest.substr(0, sp);
  h.value = trim(rest.substr(sp, colon - sp));
  std::string d = colon < rest.size() ? rest.substr(colon + 1) : std::string();
  if (las3) {
    const std::size_t lim = std::min(d.find('{'), d.find('|'));
    if (lim != std::string::npos) d = d.substr(0, lim);
  }
  h.desc = colon < rest.size() ? trim(d) : std::string();
  if (h.mnem.empty()) h.mnem = "UNKNOWN";  // a definition without a mnemonic (reported by the caller)
  return true;
}
std::vector<std::string> split_data(const std::string& line, char dlm) {
  std::vector<std::string> out;
  if (dlm == ' ') {
    std::string cur;
    bool q = false;
    for (char c : line) {
      if (c == '"') { q = !q; continue; }
      if (!q && std::isspace(static_cast<unsigned char>(c))) {
        if (!cur.empty()) out.push_back(cur);
        cur.clear();
      } else {
        cur += c;
      }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
  }
  std::string cur;
  bool q = false;
  for (char c : line) {
    if (c == '"') { q = !q; continue; }
    if (!q && c == dlm) { out.push_back(trim(cur)); cur.clear(); }
    else cur += c;
  }
  out.push_back(trim(cur));
  return out;
}

struct Group {
  std::vector<HeaderLine> defs;
  std::vector<std::vector<std::string>> rows;
};
}  // namespace

const LogCurve* WellLogSet::find(const std::string& m) const {
  for (auto& c : curves)
    if (c.mnemonic == m) return &c;
  if (index.mnemonic == m) return &index;
  return nullptr;
}

WellLogSet parse_las(const std::string& text, const std::string& source) {
  WellLogSet w;
  w.source = source;
  std::vector<std::string> lines;
  {
    std::string cur;
    for (char c : text) {
      if (c == '\n') { if (!cur.empty() && cur.back() == '\r') cur.pop_back(); lines.push_back(cur); cur.clear(); }
      else cur += c;
    }
    if (!cur.empty()) { if (cur.back() == '\r') cur.pop_back(); lines.push_back(cur); }
  }
  std::string section;       // canonical: V, W, C, P, O, A (LAS 2.0) or the LAS 3.0 title
  std::string data_of;       // the definition section a LAS 3.0 data section belongs to
  bool wrap = false, las3 = false;
  char dlm = ' ';
  std::map<std::string, Group> groups;  // LAS 3.0 groups and the LAS 2.0 log ("Log")
  std::vector<std::string> data_lines;  // LAS 2.0 ~A lines with their line numbers
  std::vector<std::size_t> data_line_no;
  for (std::size_t n = 0; n < lines.size(); ++n) {
    const std::string& raw = lines[n];
    const std::string line = trim(raw);
    const std::size_t lno = n + 1;
    if (line.empty() || line[0] == '#') continue;
    if (line[0] == '~') {
      std::string title = line.substr(1);
      std::size_t e = 0;
      while (e < title.size() && !std::isspace(static_cast<unsigned char>(title[e])) && title[e] != '|') ++e;
      std::string t = title.substr(0, e);
      data_of.clear();
      const std::size_t bar = title.find('|');
      if (bar != std::string::npos) data_of = trim(title.substr(bar + 1));
      const std::string T = upper(t);
      // LAS 2.0 sections by their first letter; LAS 3.0 keeps the 2.0 names for the log (~Curve,
      // ~Parameter, ~Other, ~A or a bare ~) and names its other groups in full
      if (T.empty()) {
        section = "A";
        w.report.push_back("line " + std::to_string(lno) + ": a section title without a name ('~" + title + "'), read as the log data (~A)");
      }
      else if (!las3 || T.size() == 1) section = std::string(1, T[0]);
      else if (T == "VERSION" || T == "WELL" || T == "CURVE" || T == "PARAMETER" || T == "OTHER" || T == "ASCII") section = std::string(1, T[0]);
      else section = T;
      continue;
    }
    if (section == "A" && !las3) {
      data_lines.push_back(line);
      data_line_no.push_back(lno);
      continue;
    }
    // LAS 3.0 data sections: ~Name_Data | Name_Definition
    if (las3 && section.size() > 5 && section.compare(section.size() - 5, 5, "_DATA") == 0) {
      const std::string g = section.substr(0, section.size() - 5);
      groups[g].rows.push_back(split_data(line, dlm));
      continue;
    }
    if (las3 && section == "A") {  // ~A in a LAS 3.0 file is the log data
      groups["LOG"].rows.push_back(split_data(line, dlm));
      continue;
    }
    if (section == "O") continue;  // ~Other: free text, kept out of the data
    HeaderLine h;
    if (!parse_header_line(line, las3, h)) {
      w.report.push_back("line " + std::to_string(lno) + ": not a header line (no period or colon) in ~" + section + ", not read: " + line);
      continue;
    }
    const std::string M = upper(h.mnem);
    if (h.mnem == "UNKNOWN" && trim(line.substr(0, line.find('.'))).empty())
      w.report.push_back("line " + std::to_string(lno) + ": a definition without a mnemonic in ~" + section + ", named UNKNOWN");
    if (section == "V") {
      if (M == "VERS") {
        w.version = h.value;
        double v = 0;
        if (parse_number(h.value, v)) {
          if (v >= 2.95) { las3 = true; w.version = "3.0"; }
          else if (v >= 1.95) w.version = "2.0";
          else w.version = h.value;
        }
      } else if (M == "WRAP") {
        wrap = upper(h.value) == "YES";
      } else if (M == "DLM") {
        const std::string d = upper(h.value);
        dlm = d == "COMMA" ? ',' : d == "TAB" ? '\t' : ' ';
      }
    } else if (section == "W") {
      w.well_info[M] = {h.unit, h.value, h.desc};
      if (M == "WELL") w.well = h.value;
    } else if (section == "C") {
      groups["LOG"].defs.push_back(h);
    } else if (section == "P") {
      w.parameters[M] = {h.unit, h.value, h.desc};
    } else if (section == "O") {
      // ~Other: free text, kept out of the data
    } else if (las3 && section.size() > 11 && section.compare(section.size() - 11, 11, "_DEFINITION") == 0) {
      groups[section.substr(0, section.size() - 11)].defs.push_back(h);
    } else if (las3 && section.size() > 10 && section.compare(section.size() - 10, 10, "_PARAMETER") == 0) {
      w.parameters[M] = {h.unit, h.value, h.desc};
    } else {
      w.report.push_back("line " + std::to_string(lno) + ": section ~" + section + " is not read");
    }
  }
  // null value
  double null_value = -999.25;
  if (auto it = w.well_info.find("NULL"); it != w.well_info.end()) {
    if (!parse_number(it->second.value, null_value)) {
      w.report.push_back("NULL '" + it->second.value + "' is not a number; -999.25 not assumed: no null value");
      null_value = kNaN;
    }
  } else {
    w.report.push_back("no NULL in the well section: no null value");
    null_value = kNaN;
  }
  // LAS 2.0 data
  if (!las3) {
    auto& g = groups["LOG"];
    const std::size_t nc = g.defs.size();
    if (wrap) {
      std::vector<std::string> tok;
      for (auto& l : data_lines)
        for (auto& t : split_data(l, ' ')) tok.push_back(t);
      if (nc && tok.size() % nc)
        w.report.push_back("wrapped data: " + std::to_string(tok.size()) + " values is not a multiple of " + std::to_string(nc) + " curves; the last " +
                           std::to_string(tok.size() % nc) + " not read");
      for (std::size_t a = 0; nc && a + nc <= tok.size(); a += nc) g.rows.emplace_back(tok.begin() + static_cast<long>(a), tok.begin() + static_cast<long>(a + nc));
    } else {
      // values beyond the defined curves become curves named UNKNOWN (no unit), reported
      std::size_t widest = 0;
      for (auto& l : data_lines) widest = std::max(widest, split_data(l, ' ').size());
      if (widest > nc) {
        w.report.push_back(std::to_string(widest - nc) + " data column" + (widest - nc == 1 ? "" : "s") + " beyond the " + std::to_string(nc) +
                           " defined curves: read as curves named UNKNOWN without a unit");
        for (std::size_t k = nc; k < widest; ++k) g.defs.push_back({"UNKNOWN", "", "", "", "", ""});
      }
      const std::size_t ncols = g.defs.size();
      std::size_t mismatched = 0;
      for (std::size_t r = 0; r < data_lines.size(); ++r) {
        auto t = split_data(data_lines[r], ' ');
        if (t.size() != ncols) {
          // read as delivered: the values go to the first curves; the others have none on this
          // line (line by line for the first five, then counted)
          if (++mismatched <= 5)
            w.report.push_back("line " + std::to_string(data_line_no[r]) + ": " + std::to_string(t.size()) + " values for " + std::to_string(ncols) +
                               " curves: the last " + std::to_string(ncols - t.size()) + " curves have no value on this line");
          t.resize(ncols);
        }
        g.rows.push_back(t);
      }
      if (mismatched > 5) w.report.push_back(std::to_string(mismatched) + " data lines in all have a number of values other than the curves'");
    }
  }
  // curves of each group
  auto make_curves = [&](Group& g, std::vector<LogCurve>& out, const std::string& gname) {
    for (std::size_t c = 0; c < g.defs.size(); ++c) {
      LogCurve lc;
      lc.mnemonic = g.defs[c].mnem;
      lc.unit = g.defs[c].unit;
      lc.description = g.defs[c].desc;
      lc.null_value = null_value;
      const bool is_text = !g.defs[c].format.empty() && std::toupper(static_cast<unsigned char>(g.defs[c].format[0])) == 'S';
      const bool is_index = gname == "LOG" && c == 0;  // an index value is a depth or time, never a null
      std::size_t bad = 0, commas = 0;
      for (auto& row : g.rows) {
        const std::string cell = c < row.size() ? row[c] : std::string();
        double v = kNaN;
        if (is_text) {
          lc.text.push_back(cell);
          lc.values.push_back(kNaN);
          continue;
        }
        bool ok = parse_number(cell, v);
        if (!ok && dlm != ',') {
          // a decimal comma ("123,42") where the delimiter is not a comma
          const std::size_t cm = cell.find(',');
          if (cm != std::string::npos && cell.find(',', cm + 1) == std::string::npos) {
            std::string t = cell;
            t[cm] = '.';
            ok = parse_number(t, v);
            commas += ok;
          }
        }
        if (!ok) {
          ++bad;
          lc.text.push_back(cell);
          lc.values.push_back(kNaN);
          continue;
        }
        if (!lc.text.empty()) lc.text.push_back("");
        if (!is_index && std::isfinite(null_value) && std::fabs(v - null_value) <= 1e-9 * std::max(1.0, std::fabs(null_value))) {
          ++lc.nulls;
          v = kNaN;
        }
        lc.values.push_back(v);
      }
      if (bad) w.report.push_back(gname + " " + lc.mnemonic + ": " + std::to_string(bad) + " values are not numbers, kept as text");
      if (commas) w.report.push_back(gname + " " + lc.mnemonic + ": " + std::to_string(commas) + " values with a decimal comma, read as decimal points");
      out.push_back(std::move(lc));
    }
  };
  for (auto& [name, g] : groups) {
    std::vector<LogCurve> cs;
    make_curves(g, cs, name);
    if (name == "LOG") {
      if (cs.empty()) { w.report.push_back("no curves"); continue; }
      w.index = cs.front();
      w.curves.assign(cs.begin() + 1, cs.end());
      w.index_type = index_type_of(w.index.mnemonic, w.index.unit);
    } else if (name == "TOPS") {
      w.tops = cs;
    } else if (name == "INCLINOMETRY") {
      w.inclinometry = cs;
    } else if (!g.defs.empty()) {
      w.report.push_back("LAS 3.0 group " + name + " (" + std::to_string(g.rows.size()) + " rows) is not used");
    }
  }
  if (w.version.empty()) w.report.push_back("no VERS in a version section");
  // report: units, nulls, header against data
  w.report.push_back("LAS " + w.version + (wrap ? ", wrapped" : "") + ": index " + w.index.mnemonic + " (" + (w.index.unit.empty() ? "no unit" : w.index.unit) +
                     ", " + (w.index_type.empty() ? "index type not stated" : w.index_type) + "), " + std::to_string(w.index.values.size()) + " samples, " +
                     std::to_string(w.curves.size()) + " curves; null value " + (std::isfinite(null_value) ? fmt(null_value) : std::string("none")));
  for (auto& c : w.curves)
    w.report.push_back("curve " + c.mnemonic + " [" + (c.unit.empty() ? "no unit" : c.unit) + "]: " + std::to_string(c.nulls) + " null" + (c.nulls == 1 ? "" : "s"));
  if (!w.index.values.empty()) {
    double strt = 0, stop = 0, step = 0;
    auto get = [&](const char* k, double& v) { auto it = w.well_info.find(k); return it != w.well_info.end() && parse_number(it->second.value, v); };
    if (get("STRT", strt) && std::fabs(strt - w.index.values.front()) > 1e-6) w.report.push_back("STRT " + fmt(strt) + " differs from the first index " + fmt(w.index.values.front()));
    if (get("STOP", stop) && std::fabs(stop - w.index.values.back()) > 1e-6) w.report.push_back("STOP " + fmt(stop) + " differs from the last index " + fmt(w.index.values.back()));
    if (get("STEP", step) && step != 0 && w.index.values.size() > 1) {
      std::size_t off = 0;
      for (std::size_t i = 1; i < w.index.values.size(); ++i)
        if (std::fabs(w.index.values[i] - w.index.values[i - 1] - step) > 1e-6 * std::max(1.0, std::fabs(step))) ++off;
      if (off) w.report.push_back(std::to_string(off) + " index steps differ from STEP " + fmt(step));
    }
    std::size_t nonmono = 0;
    for (std::size_t i = 1; i < w.index.values.size(); ++i) {
      const double d = w.index.values[i] - w.index.values[i - 1];
      if (!(d * (w.index.values.back() - w.index.values.front()) > 0)) ++nonmono;
    }
    if (nonmono) w.report.push_back(std::to_string(nonmono) + " index samples repeat or go backwards");
  }
  return w;
}

WellLogSet read_las(const std::filesystem::path& file) {
  std::ifstream in(file, std::ios::binary);
  if (!in) throw Error(ErrorCode::Io, "cannot open " + file.string());
  std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  std::vector<std::string> enc;
  // UTF-16 (with a byte-order mark, or recognised by its zero bytes) is converted to UTF-8
  {
    const auto b = [&](std::size_t i) { return static_cast<unsigned char>(s[i]); };
    int order = 0;  // 1 little-endian, 2 big-endian
    std::size_t start = 0;
    if (s.size() >= 2 && b(0) == 0xFF && b(1) == 0xFE) { order = 1; start = 2; }
    else if (s.size() >= 2 && b(0) == 0xFE && b(1) == 0xFF) { order = 2; start = 2; }
    else if (s.size() >= 4) {
      std::size_t z_even = 0, z_odd = 0;
      for (std::size_t i = 0; i + 1 < std::min<std::size_t>(s.size(), 400); i += 2) { z_even += b(i) == 0; z_odd += b(i + 1) == 0; }
      if (z_odd > 50 && z_even == 0) order = 1;
      if (z_even > 50 && z_odd == 0) order = 2;
    }
    if (order) {
      std::string u;
      for (std::size_t i = start; i + 1 < s.size(); i += 2) {
        std::uint32_t cp = order == 1 ? (b(i) | (b(i + 1) << 8)) : ((b(i) << 8) | b(i + 1));
        if (cp >= 0xD800 && cp < 0xDC00 && i + 3 < s.size()) {
          const std::uint32_t lo = order == 1 ? (b(i + 2) | (b(i + 3) << 8)) : ((b(i + 2) << 8) | b(i + 3));
          cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
          i += 2;
        }
        if (cp < 0x80) u += static_cast<char>(cp);
        else if (cp < 0x800) { u += static_cast<char>(0xC0 | (cp >> 6)); u += static_cast<char>(0x80 | (cp & 0x3F)); }
        else if (cp < 0x10000) { u += static_cast<char>(0xE0 | (cp >> 12)); u += static_cast<char>(0x80 | ((cp >> 6) & 0x3F)); u += static_cast<char>(0x80 | (cp & 0x3F)); }
        else { u += static_cast<char>(0xF0 | (cp >> 18)); u += static_cast<char>(0x80 | ((cp >> 12) & 0x3F)); u += static_cast<char>(0x80 | ((cp >> 6) & 0x3F)); u += static_cast<char>(0x80 | (cp & 0x3F)); }
      }
      s = u;
      enc.push_back(std::string("encoding: UTF-16 ") + (order == 1 ? "little" : "big") + "-endian" + (start ? " with a byte-order mark" : "") + ", converted to UTF-8");
    } else if (s.size() >= 3 && b(0) == 0xEF && b(1) == 0xBB && b(2) == 0xBF) {
      s.erase(0, 3);
      enc.push_back("encoding: UTF-8 with a byte-order mark");
    }
  }
  // encoding: ASCII or UTF-8; another byte is reported with its offset
  for (std::size_t i = 0; i < s.size();) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if (c < 0x80) { ++i; continue; }
    int n = (c >> 5) == 6 ? 1 : (c >> 4) == 14 ? 2 : (c >> 3) == 30 ? 3 : -1;
    bool ok = n > 0 && i + static_cast<std::size_t>(n) < s.size() + 0;
    for (int k = 1; ok && k <= n; ++k) ok = (static_cast<unsigned char>(s[i + static_cast<std::size_t>(k)]) >> 6) == 2;
    if (!ok) {
      if (enc.size() < 5) enc.push_back("byte " + std::to_string(i) + " (0x" + [&] { std::ostringstream o; o << std::hex << int(c); return o.str(); }() +
                                        ") is neither ASCII nor UTF-8");
      ++i;
    } else {
      i += static_cast<std::size_t>(n) + 1;
    }
  }
  auto w = parse_las(s, file.string());
  w.report.insert(w.report.end(), enc.begin(), enc.end());
  return w;
}

// ---- DLIS (RP66 V1, decision 21) ---------------------------------------------------------------
namespace {
struct Reader {
  const std::uint8_t* p;
  std::size_t n, i = 0;
  bool ok = true;
  bool need(std::size_t k) { if (i + k > n) { ok = false; return false; } return true; }
  unsigned u8() { return need(1) ? p[i++] : 0; }
  unsigned u16() { if (!need(2)) return 0; unsigned v = (unsigned(p[i]) << 8) | p[i + 1]; i += 2; return v; }
  std::uint32_t u32() { if (!need(4)) return 0; std::uint32_t v = (std::uint32_t(p[i]) << 24) | (std::uint32_t(p[i + 1]) << 16) | (std::uint32_t(p[i + 2]) << 8) | p[i + 3]; i += 4; return v; }
  std::uint32_t uvari() {
    if (!need(1)) return 0;
    const unsigned b = p[i];
    if ((b & 0x80) == 0) { ++i; return b; }
    if ((b & 0xC0) == 0x80) return u16() & 0x3FFF;
    return u32() & 0x3FFFFFFF;
  }
  std::string ident() { const unsigned k = u8(); if (!need(k)) return {}; std::string s(reinterpret_cast<const char*>(p + i), k); i += k; return s; }
  std::string ascii() { const auto k = uvari(); if (!need(k)) return {}; std::string s(reinterpret_cast<const char*>(p + i), k); i += k; return s; }
};
float f32(std::uint32_t x) { float f; std::memcpy(&f, &x, 4); return f; }
double ibm32(std::uint32_t x) {
  if ((x & 0x7FFFFFFF) == 0) return 0.0;
  const double s = (x >> 31) ? -1.0 : 1.0;
  const int e = static_cast<int>((x >> 24) & 0x7F) - 64;
  return s * static_cast<double>(x & 0xFFFFFF) / 16777216.0 * std::pow(16.0, e);
}
double vax32(std::uint32_t raw) {  // VAX F-floating, stored as two 16-bit words, low word first
  const std::uint32_t x = ((raw & 0xFFFF) << 16) | (raw >> 16);
  const int e = static_cast<int>((x >> 23) & 0xFF);
  if (e == 0) return 0.0;
  const double s = (x >> 31) ? -1.0 : 1.0;
  return s * (0.5 + static_cast<double>(x & 0x7FFFFF) / 16777216.0) * std::pow(2.0, e - 128);
}
struct Value { std::vector<double> num; std::vector<std::string> str; };
/// One value of a representation code (API RP66 V1, Appendix B). Strings for identifiers, names
/// and references; numbers otherwise (complex values: real parts; validated: the value only).
bool read_rep(Reader& r, unsigned code, Value& out) {
  switch (code) {
    case 1: { const unsigned v = r.u16(); const int frac = static_cast<int>(static_cast<std::int16_t>(v & 0xFFF0)) >> 4; out.num.push_back(frac / 2048.0 * std::pow(2.0, v & 0xF)); return r.ok; }
    case 2: out.num.push_back(f32(r.u32())); return r.ok;
    case 3: { out.num.push_back(f32(r.u32())); r.u32(); return r.ok; }
    case 4: { out.num.push_back(f32(r.u32())); r.u32(); r.u32(); return r.ok; }
    case 5: out.num.push_back(ibm32(r.u32())); return r.ok;
    case 6: { const std::uint32_t b = r.u32(); out.num.push_back(vax32(((b & 0xFF) << 24) | ((b & 0xFF00) << 8) | ((b >> 8) & 0xFF00) | (b >> 24))); return r.ok; }
    case 7: { std::uint64_t x = (std::uint64_t(r.u32()) << 32) | r.u32(); double d; std::memcpy(&d, &x, 8); out.num.push_back(d); return r.ok; }
    case 8: { std::uint64_t x = (std::uint64_t(r.u32()) << 32) | r.u32(); double d; std::memcpy(&d, &x, 8); out.num.push_back(d); r.u32(); r.u32(); return r.ok; }
    case 9: { std::uint64_t x = (std::uint64_t(r.u32()) << 32) | r.u32(); double d; std::memcpy(&d, &x, 8); out.num.push_back(d); for (int k = 0; k < 4; ++k) r.u32(); return r.ok; }
    case 10: { out.num.push_back(f32(r.u32())); r.u32(); return r.ok; }
    case 11: { std::uint64_t x = (std::uint64_t(r.u32()) << 32) | r.u32(); double d; std::memcpy(&d, &x, 8); out.num.push_back(d); r.u32(); r.u32(); return r.ok; }
    case 12: out.num.push_back(static_cast<std::int8_t>(r.u8())); return r.ok;
    case 13: out.num.push_back(static_cast<std::int16_t>(r.u16())); return r.ok;
    case 14: out.num.push_back(static_cast<std::int32_t>(r.u32())); return r.ok;
    case 15: out.num.push_back(r.u8()); return r.ok;
    case 16: out.num.push_back(r.u16()); return r.ok;
    case 17: out.num.push_back(r.u32()); return r.ok;
    case 18: case 22: out.num.push_back(r.uvari()); return r.ok;
    case 19: case 27: out.str.push_back(r.ident()); return r.ok;
    case 20: out.str.push_back(r.ascii()); return r.ok;
    case 21: { for (int k = 0; k < 8; ++k) r.u8(); out.str.push_back("dtime"); return r.ok; }
    case 23: { r.uvari(); r.u8(); out.str.push_back(r.ident()); return r.ok; }
    case 24: { r.ident(); r.uvari(); r.u8(); out.str.push_back(r.ident()); return r.ok; }
    case 25: { r.ident(); r.uvari(); r.u8(); r.ident(); out.str.push_back(r.ident()); return r.ok; }
    case 26: out.num.push_back(r.u8()); return r.ok;
    default: return false;
  }
}
std::size_t rep_size(unsigned code) {  // fixed sizes of frame-data codes; 0 for variable
  static const std::size_t s[28] = {0, 2, 4, 8, 12, 4, 4, 8, 16, 24, 8, 16, 1, 2, 4, 1, 2, 4, 0, 0, 0, 8, 0, 0, 0, 0, 1, 0};
  return code < 28 ? s[code] : 0;
}
struct Attr { std::string label; std::uint32_t count = 1; unsigned rep = 19; std::string units; Value value; bool present = false; };
struct Object { std::string name; std::map<std::string, Attr> attrs; };
struct Set { std::string type, name; std::vector<Object> objects; };
/// An explicitly formatted logical record: set, template, objects.
bool parse_eflr(const std::vector<std::uint8_t>& lr, Set& set, std::string& err) {
  Reader r{lr.data(), lr.size()};
  const unsigned d = r.u8();
  const unsigned role = d >> 5;
  if (role < 5) { err = "EFLR does not start with a set component"; return false; }
  if (d & 0x10) set.type = r.ident();
  if (d & 0x08) set.name = r.ident();
  std::vector<Attr> tmpl;
  // template
  while (r.ok && r.i < r.n) {
    const unsigned c = r.p[r.i];
    if ((c >> 5) == 3) break;  // first object
    r.u8();
    if ((c >> 5) != 1 && (c >> 5) != 2) { err = "unexpected component in the template"; return false; }
    Attr a;
    if (c & 0x10) a.label = r.ident();
    if (c & 0x08) a.count = r.uvari();
    if (c & 0x04) a.rep = r.u8();
    if (c & 0x02) a.units = r.ident();
    if (c & 0x01) { for (std::uint32_t k = 0; k < a.count; ++k) if (!read_rep(r, a.rep, a.value)) { err = "template value of representation code " + std::to_string(a.rep); return false; } a.present = true; }
    tmpl.push_back(a);
  }
  while (r.ok && r.i < r.n) {
    const unsigned c = r.u8();
    if ((c >> 5) != 3) { err = "expected an object component"; return false; }
    Object o;
    if (c & 0x10) { r.uvari(); r.u8(); o.name = r.ident(); }
    for (auto& t : tmpl) o.attrs[t.label] = t;
    for (std::size_t k = 0; k < tmpl.size() && r.ok && r.i < r.n; ++k) {
      const unsigned a = r.p[r.i];
      if ((a >> 5) == 3) break;  // the next object: remaining attributes take the template defaults
      r.u8();
      if ((a >> 5) == 0) { o.attrs[tmpl[k].label].present = false; o.attrs[tmpl[k].label].value = {}; continue; }  // absent
      Attr at = tmpl[k];
      if (a & 0x10) at.label = r.ident();
      if (a & 0x08) at.count = r.uvari();
      if (a & 0x04) at.rep = r.u8();
      if (a & 0x02) at.units = r.ident();
      if (a & 0x01) {
        at.value = {};
        for (std::uint32_t q = 0; q < at.count; ++q)
          if (!read_rep(r, at.rep, at.value)) { err = "value of " + tmpl[k].label + " (representation code " + std::to_string(at.rep) + ")"; return false; }
        at.present = true;
      }
      o.attrs[tmpl[k].label] = at;
    }
    set.objects.push_back(std::move(o));
  }
  if (!r.ok) { err = "record ends inside a component"; return false; }
  return true;
}
}  // namespace

std::vector<WellLogSet> read_dlis(const std::filesystem::path& file) {
  std::ifstream in(file, std::ios::binary);
  if (!in) throw Error(ErrorCode::Io, "cannot open " + file.string());
  std::vector<std::uint8_t> raw((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  std::vector<std::string> report;
  // tape-image (TIF) wrapper: 12-byte markers of type, previous and next offsets (little-endian)
  auto le32 = [&](std::size_t o) { return std::uint32_t(raw[o]) | (std::uint32_t(raw[o + 1]) << 8) | (std::uint32_t(raw[o + 2]) << 16) | (std::uint32_t(raw[o + 3]) << 24); };
  std::vector<std::uint8_t> b;
  if (raw.size() > 21 && std::memcmp(raw.data() + 16, "V1.00", 5) == 0 && le32(0) == 0) {
    std::size_t pos = 0;
    while (pos + 12 <= raw.size()) {
      const std::uint32_t type = le32(pos), next = le32(pos + 8);
      if (type != 0 || next <= pos + 12 || next > raw.size()) break;
      b.insert(b.end(), raw.begin() + static_cast<long>(pos + 12), raw.begin() + static_cast<long>(next));
      pos = next;
    }
    report.push_back("tape-image (TIF) markers removed");
  } else {
    b = std::move(raw);
  }
  if (b.size() < 80 || std::memcmp(b.data() + 4, "V1.00", 5) != 0)
    throw Error(ErrorCode::Parse, file.string() + ": no RP66 V1 storage unit label (bytes 5 to 9 are not 'V1.00')");
  report.push_back("storage unit label: " + std::string(reinterpret_cast<const char*>(b.data() + 15), 5) + " maximum record length, set '" +
                   trim(std::string(reinterpret_cast<const char*>(b.data() + 20), 60)) + "'");
  // visible records -> logical record segments -> logical records
  struct LR { bool eflr; unsigned type; std::vector<std::uint8_t> body; std::size_t offset; bool encrypted; };
  std::vector<LR> lrs;
  std::size_t pos = 80;
  LR cur{};
  bool open = false;
  while (pos + 4 <= b.size()) {
    const std::size_t vlen = (std::size_t(b[pos]) << 8) | b[pos + 1];
    if (vlen < 4 || pos + vlen > b.size() || b[pos + 2] != 0xFF) {
      report.push_back("byte " + std::to_string(pos) + ": invalid visible record header; the rest of the file is not read");
      break;
    }
    std::size_t s = pos + 4;
    const std::size_t vend = pos + vlen;
    while (s + 4 <= vend) {
      const std::size_t slen = (std::size_t(b[s]) << 8) | b[s + 1];
      const unsigned attr = b[s + 2], type = b[s + 3];
      if (slen < 4 || s + slen > vend) {
        report.push_back("byte " + std::to_string(s) + ": invalid segment length " + std::to_string(slen));
        break;
      }
      std::size_t e = s + slen;
      if (attr & 0x02) e -= 2;  // trailing length
      if (attr & 0x04) e -= 2;  // checksum
      if (attr & 0x01) {        // padding: the last byte counts the pad bytes
        const std::size_t pad = b[e - 1];
        e = pad <= e - s - 4 ? e - pad : s + 4;
      }
      if (!(attr & 0x40)) {     // no predecessor: a new logical record
        if (open) report.push_back("byte " + std::to_string(s) + ": a logical record ended without its last segment");
        cur = LR{(attr & 0x80) != 0, type, {}, s, (attr & 0x10) != 0};
        open = true;
      }
      if (open) cur.body.insert(cur.body.end(), b.begin() + static_cast<long>(s + 4 + ((attr & 0x08) ? 0 : 0)), b.begin() + static_cast<long>(e));
      if (!(attr & 0x20) && open) { lrs.push_back(std::move(cur)); open = false; }
      s += slen;
    }
    pos = vend;
  }
  // logical files: a FILE-HEADER starts each
  std::vector<WellLogSet> out;
  struct Chan { std::string unit; unsigned rep = 2; std::uint32_t dim = 1; };
  std::map<std::string, Chan> chans;
  struct Frame { std::string name, index_type; std::vector<std::string> channels; std::vector<std::vector<double>> data; std::size_t skipped = 0; };
  std::vector<Frame> frames;
  std::string well, fileid;
  std::size_t nfile = 0, encrypted = 0;
  auto flush = [&]() {
    for (auto& f : frames) {
      if (f.channels.empty()) continue;
      WellLogSet w;
      w.version = "RP66 V1";
      w.well = well;
      w.source = file.string() + " (logical file " + std::to_string(nfile) + (fileid.empty() ? "" : " '" + fileid + "'") + ", frame " + f.name + ")";
      const std::string it = upper(f.index_type);
      w.index_type = it == "BOREHOLE-DEPTH" ? "MD" : it == "VERTICAL-DEPTH" ? "TVD" : it == "TIME" ? "TIME" : "";
      for (std::size_t c = 0; c < f.channels.size(); ++c) {
        LogCurve lc;
        lc.mnemonic = f.channels[c];
        lc.unit = chans[lc.mnemonic].unit;
        lc.null_value = kNaN;  // RP66 V1 has no file-wide null value
        for (auto& row : f.data) lc.values.push_back(c < row.size() ? row[c] : kNaN);
        if (c == 0) w.index = std::move(lc);
        else w.curves.push_back(std::move(lc));
      }
      w.report = report;
      w.report.push_back("frame " + f.name + ": index " + w.index.mnemonic + " (" + w.index.unit + ", " + (w.index_type.empty() ? f.index_type : w.index_type) +
                         "), " + std::to_string(f.data.size()) + " frames, " + std::to_string(w.curves.size()) + " channels; no null value (RP66 V1 has none)");
      if (f.skipped) w.report.push_back(std::to_string(f.skipped) + " frames not read (channel arrays or unknown representation codes)");
      out.push_back(std::move(w));
    }
    frames.clear();
    chans.clear();
  };
  for (auto& lr : lrs) {
    if (lr.encrypted) { ++encrypted; continue; }
    if (lr.eflr) {
      Set s;
      std::string err;
      if (!parse_eflr(lr.body, s, err)) {
        report.push_back("byte " + std::to_string(lr.offset) + ": explicitly formatted record not read: " + err);
        continue;
      }
      if (s.type == "FILE-HEADER") {
        if (nfile) flush();
        ++nfile;
        fileid.clear();
        for (auto& o : s.objects)
          if (o.attrs.count("ID") && !o.attrs["ID"].value.str.empty()) fileid = trim(o.attrs["ID"].value.str[0]);
      } else if (s.type == "ORIGIN") {
        for (auto& o : s.objects)
          if (well.empty() && o.attrs.count("WELL-NAME") && !o.attrs["WELL-NAME"].value.str.empty()) well = trim(o.attrs["WELL-NAME"].value.str[0]);
      } else if (s.type == "CHANNEL") {
        for (auto& o : s.objects) {
          Chan c;
          if (o.attrs.count("UNITS") && !o.attrs["UNITS"].value.str.empty()) c.unit = o.attrs["UNITS"].value.str[0];
          if (o.attrs.count("REPRESENTATION-CODE") && !o.attrs["REPRESENTATION-CODE"].value.num.empty())
            c.rep = static_cast<unsigned>(o.attrs["REPRESENTATION-CODE"].value.num[0]);
          if (o.attrs.count("DIMENSION")) {
            c.dim = 1;
            for (double v : o.attrs["DIMENSION"].value.num) c.dim *= static_cast<std::uint32_t>(std::max(1.0, v));
          }
          chans[o.name] = c;
        }
      } else if (s.type == "FRAME") {
        for (auto& o : s.objects) {
          Frame f;
          f.name = o.name;
          if (o.attrs.count("CHANNELS")) f.channels = o.attrs["CHANNELS"].value.str;
          if (o.attrs.count("INDEX-TYPE") && !o.attrs["INDEX-TYPE"].value.str.empty()) f.index_type = o.attrs["INDEX-TYPE"].value.str[0];
          frames.push_back(f);
        }
      }
    } else if (lr.type == 0) {  // FDATA
      Reader r{lr.body.data(), lr.body.size()};
      r.uvari(); r.u8();
      const std::string fname = r.ident();
      r.uvari();  // frame number
      Frame* f = nullptr;
      for (auto& x : frames)
        if (x.name == fname) f = &x;
      if (!f) { report.push_back("frame data for an undeclared frame " + fname); continue; }
      std::vector<double> row;
      bool ok = true;
      for (auto& cn : f->channels) {
        auto it = chans.find(cn);
        if (it == chans.end() || it->second.dim != 1) { ok = false; break; }
        Value v;
        if (!read_rep(r, it->second.rep, v) || v.num.empty()) { ok = false; break; }
        row.push_back(v.num[0]);
      }
      if (ok) f->data.push_back(row);
      else ++f->skipped;
    }
  }
  if (encrypted) report.push_back(std::to_string(encrypted) + " encrypted logical records not read");
  if (!nfile) nfile = 1;
  flush();
  (void)rep_size;
  return out;
}

}  // namespace stratonix::io
