// Line, surface, point and well file formats: writers and readers.
#include <algorithm>
#include <cmath>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/io/export.hpp"

namespace stratonix::io {

namespace {
std::vector<std::string> lines_of(const std::string& text) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : text) {
    if (c == '\n') { if (!cur.empty() && cur.back() == '\r') cur.pop_back(); out.push_back(cur); cur.clear(); }
    else cur += c;
  }
  if (!cur.empty()) { if (cur.back() == '\r') cur.pop_back(); out.push_back(cur); }
  return out;
}
/// Whitespace tokens; double-quoted strings are one token (quotes removed).
std::vector<std::string> fields(const std::string& line) {
  std::vector<std::string> t;
  std::size_t i = 0;
  while (i < line.size()) {
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    if (i >= line.size()) break;
    if (line[i] == '"') {
      std::size_t j = line.find('"', i + 1);
      if (j == std::string::npos) throw Error(ErrorCode::Parse, "unterminated quoted string: " + line);
      t.push_back(line.substr(i + 1, j - i - 1));
      i = j + 1;
    } else {
      std::size_t j = i;
      while (j < line.size() && !std::isspace(static_cast<unsigned char>(line[j]))) ++j;
      t.push_back(line.substr(i, j - i));
      i = j;
    }
  }
  return t;
}
std::string quote(const std::string& s) {
  if (s.find('"') != std::string::npos) throw Error(ErrorCode::InvalidArgument, "names may not contain '\"': " + s);
  return "\"" + s + "\"";
}
std::string xyz(const Vec3& p, int d) { return fmt_fixed(p.x, d) + " " + fmt_fixed(p.y, d) + " " + fmt_fixed(p.z, d); }
}  // namespace

std::string petrel_name(const std::string& s) {
  std::string r;
  for (char c : s) r += (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-') ? c : '_';
  if (r.empty()) r = "_";
  return r;
}

// ---- IRAP classic lines -------------------------------------------------------------------
std::string write_irap_lines(const std::vector<Polyline>& lines, int d) {
  std::string o;
  const std::string sep = fmt_fixed(kIrapLineSeparator, 6) + " " + fmt_fixed(kIrapLineSeparator, 6) + " " + fmt_fixed(kIrapLineSeparator, 6) + "\n";
  for (auto& l : lines) {
    if (l.size() < 2) continue;
    for (auto& p : l) o += xyz(p, d) + "\n";
    o += sep;
  }
  return o;
}

std::vector<Polyline> read_irap_lines(const std::string& text) {
  std::vector<Polyline> out;
  Polyline cur;
  for (auto& ln : lines_of(text)) {
    auto t = fields(ln);
    if (t.empty() || t[0][0] == '#' || t[0][0] == '!') continue;
    if (t.size() < 3) throw Error(ErrorCode::Parse, "IRAP lines: expected X Y Z, got: " + ln);
    Vec3 p{parse_double(t[0]), parse_double(t[1]), parse_double(t[2])};
    if (p.x == kIrapLineSeparator && p.y == kIrapLineSeparator && p.z == kIrapLineSeparator) {
      if (!cur.empty()) out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(p);
    }
  }
  if (!cur.empty()) out.push_back(cur);
  return out;
}

// ---- GOCAD TSurf ---------------------------------------------------------------------------
std::string write_tsurf(const TriMesh& m, const TSurfInfo& in, int d) {
  const bool prop = !in.property.empty() && in.values.size() == m.v.size();
  std::string o = "GOCAD TSurf 1\nHEADER {\nname:" + in.name + "\n*solid*color:" + fmt_fixed(in.r, 4) + " " + fmt_fixed(in.g, 4) + " " +
                  fmt_fixed(in.b, 4) + " 1\n}\n";
  o += "GOCAD_ORIGINAL_COORDINATE_SYSTEM\nNAME " + in.crs_name + "\nAXIS_NAME \"X\" \"Y\" \"Z\"\nAXIS_UNIT \"m\" \"m\" \"m\"\nZPOSITIVE " +
       in.z_positive + "\nEND_ORIGINAL_COORDINATE_SYSTEM\n";
  if (prop) o += "PROPERTIES " + in.property + "\nPROPERTY_CLASSES " + in.property + "\nPROPERTY_KINDS \"unknown\"\nESIZES 1\nUNITS m\n";
  o += "TFACE\n";
  // only vertices used by triangles, renumbered from 1
  std::vector<long> id(m.v.size(), -1);
  long n = 0;
  for (auto& t : m.t)
    for (auto k : t)
      if (id[k] < 0) id[k] = 0;
  for (std::size_t i = 0; i < m.v.size(); ++i)
    if (id[i] == 0) {
      id[i] = ++n;
      o += (prop ? "PVRTX " : "VRTX ") + std::to_string(n) + " " + xyz(m.v[i], d) + (prop ? " " + fmt_fixed(in.values[i], d) : "") + "\n";
    }
  for (auto& t : m.t) o += "TRGL " + std::to_string(id[t[0]]) + " " + std::to_string(id[t[1]]) + " " + std::to_string(id[t[2]]) + "\n";
  o += "END\n";
  return o;
}

TriMesh read_tsurf(const std::string& text, TSurfInfo* info) {
  TriMesh m;
  std::map<long, std::size_t> idx;
  bool header = false, gocad = false;
  for (auto& ln : lines_of(text)) {
    std::string tl = trim(ln);
    if (tl.empty()) continue;
    if (tl.rfind("GOCAD TSurf", 0) == 0) { gocad = true; continue; }
    if (tl.rfind("HEADER", 0) == 0) { header = true; continue; }
    if (header) {
      if (tl[0] == '}') { header = false; continue; }
      if (info && tl.rfind("name:", 0) == 0) info->name = tl.substr(5);
      continue;
    }
    auto t = fields(tl);
    if (t[0] == "ZPOSITIVE" && info && t.size() > 1) info->z_positive = t[1];
    if (t[0] == "NAME" && info && t.size() > 1) info->crs_name = t[1];
    if (t[0] == "PROPERTIES" && info && t.size() > 1) info->property = t[1];
    if (t[0] == "VRTX" || t[0] == "PVRTX") {
      if (t.size() < 5) throw Error(ErrorCode::Parse, "TSurf: short VRTX line: " + tl);
      idx[std::stol(t[1])] = m.v.size();
      m.v.push_back({parse_double(t[2]), parse_double(t[3]), parse_double(t[4])});
      if (t[0] == "PVRTX" && info && t.size() > 5) info->values.push_back(parse_double(t[5]));
    } else if (t[0] == "TRGL") {
      if (t.size() < 4) throw Error(ErrorCode::Parse, "TSurf: short TRGL line: " + tl);
      std::array<std::uint32_t, 3> tri{};
      for (int k = 0; k < 3; ++k) {
        auto it = idx.find(std::stol(t[1 + k]));
        if (it == idx.end()) throw Error(ErrorCode::Parse, "TSurf: triangle refers to an unknown vertex: " + tl);
        tri[k] = static_cast<std::uint32_t>(it->second);
      }
      m.t.push_back(tri);
    }
  }
  if (!gocad) throw Error(ErrorCode::Parse, "not a GOCAD TSurf file");
  return m;
}

// ---- Petrel points with attributes ---------------------------------------------------------
std::string write_petrel_points(const PointsTable& t, const std::vector<std::string>& comments, int d) {
  std::string o = "# Petrel Points with attributes\n";
  for (auto& c : comments) o += "# " + c + "\n";
  o += "VERSION 1\nBEGIN HEADER\nX\nY\nZ\n";
  for (auto& [name, type] : t.attr) o += (type == "Float" ? name : type + "," + name) + "\n";
  o += "END HEADER\n";
  for (std::size_t i = 0; i < t.pts.size(); ++i) {
    o += xyz(t.pts[i], d);
    for (std::size_t a = 0; a < t.attr.size(); ++a) o += " " + (t.attr[a].second == "String" ? quote(t.values[i][a]) : t.values[i][a]);
    o += "\n";
  }
  return o;
}

PointsTable read_petrel_points(const std::string& text) {
  PointsTable t;
  bool in_header = false, done_header = false;
  int coords = 0;
  for (auto& ln : lines_of(text)) {
    std::string tl = trim(ln);
    if (tl.empty() || tl[0] == '#') continue;
    if (tl == "VERSION 1") continue;
    if (tl == "BEGIN HEADER") { in_header = true; continue; }
    if (tl == "END HEADER") { in_header = false; done_header = true; continue; }
    if (in_header) {
      if (coords < 3) { ++coords; continue; }  // X, Y, Z
      auto c = tl.find(',');
      if (c == std::string::npos) t.attr.push_back({tl, "Float"});
      else t.attr.push_back({tl.substr(c + 1), tl.substr(0, c)});
      continue;
    }
    if (!done_header) throw Error(ErrorCode::Parse, "Petrel points: data before END HEADER");
    auto f = fields(tl);
    if (f.size() < 3 + t.attr.size()) throw Error(ErrorCode::Parse, "Petrel points: short row: " + tl);
    t.pts.push_back({parse_double(f[0]), parse_double(f[1]), parse_double(f[2])});
    t.values.emplace_back(f.begin() + 3, f.begin() + 3 + static_cast<long>(t.attr.size()));
  }
  return t;
}

// ---- Wells ---------------------------------------------------------------------------------
std::string write_well_heads(const std::vector<WellHead>& w, const std::vector<std::string>& comments, int d) {
  std::string o = "# Petrel well heads\n";
  for (auto& c : comments) o += "# " + c + "\n";
  o += "VERSION 1\nBEGIN HEADER\nWell name\nSurface X\nSurface Y\nKB\nTD (MD)\nEND HEADER\n";
  for (auto& h : w)
    o += quote(h.name) + " " + fmt_fixed(h.x, d) + " " + fmt_fixed(h.y, d) + " " + fmt_fixed(h.kb, d) + " " + fmt_fixed(h.td_md, d) + "\n";
  return o;
}

std::vector<WellHead> read_well_heads(const std::string& text) {
  std::vector<WellHead> out;
  bool data = false;
  for (auto& ln : lines_of(text)) {
    std::string tl = trim(ln);
    if (tl.empty() || tl[0] == '#') continue;
    if (tl == "END HEADER") { data = true; continue; }
    if (!data) continue;
    auto f = fields(tl);
    if (f.size() < 5) throw Error(ErrorCode::Parse, "well heads: short row: " + tl);
    out.push_back({f[0], parse_double(f[1]), parse_double(f[2]), parse_double(f[3]), parse_double(f[4])});
  }
  return out;
}

std::string write_well_tops(const std::vector<WellTopRow>& t, const std::vector<std::string>& comments, int d) {
  std::string o = "# Petrel well tops\n";
  for (auto& c : comments) o += "# " + c + "\n";
  o += "VERSION 1\nBEGIN HEADER\nX\nY\nZ\nMD\nWell\nSurface\nEND HEADER\n";
  for (auto& r : t)
    o += fmt_fixed(r.x, d) + " " + fmt_fixed(r.y, d) + " " + fmt_fixed(r.z, d) + " " + fmt_fixed(r.md, d) + " " + quote(r.well) + " " +
         quote(r.surface) + "\n";
  return o;
}

std::vector<WellTopRow> read_well_tops(const std::string& text) {
  std::vector<WellTopRow> out;
  bool data = false;
  for (auto& ln : lines_of(text)) {
    std::string tl = trim(ln);
    if (tl.empty() || tl[0] == '#') continue;
    if (tl == "END HEADER") { data = true; continue; }
    if (!data) continue;
    auto f = fields(tl);
    if (f.size() < 6) throw Error(ErrorCode::Parse, "well tops: short row: " + tl);
    out.push_back({f[4], f[5], parse_double(f[0]), parse_double(f[1]), parse_double(f[2]), parse_double(f[3])});
  }
  return out;
}

std::string write_well_path(const WellPath& w, int d) {
  std::string o = "# Well path (deviation) of " + w.name() + "\n# Columns: MD X Y Z; Z is elevation (m, negative below the datum)\n";
  o += "# Well head: " + fmt_fixed(w.x(), d) + " " + fmt_fixed(w.y(), d) + "; datum (KB) elevation " + fmt_fixed(w.datum_z(), d) + "\n";
  o += "MD X Y Z\n";
  const auto& p = w.points();
  double md = 0;
  for (std::size_t i = 0; i < p.size(); ++i) {
    if (i) md += norm(p[i] - p[i - 1]);
    o += fmt_fixed(md, d) + " " + xyz(p[i], d) + "\n";
  }
  return o;
}

}  // namespace stratonix::io
