// Grid formats (IRAP classic ASCII, ZMAP+, regular XYZ) and the grid importer plugins.
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/io/io.hpp"

// Plugin anchors called below (C linkage, defined by STRATONIX_REGISTER in each plugin file).
STRATONIX_DECLARE_PLUGIN_ANCHOR(imp_irap_ascii);
STRATONIX_DECLARE_PLUGIN_ANCHOR(imp_zmap);
STRATONIX_DECLARE_PLUGIN_ANCHOR(imp_xyz);

namespace stratonix::io {

static std::vector<std::string> tokens(const std::string& text, bool comma = false) {
  std::vector<std::string> t;
  std::string cur;
  for (char c : text) {
    bool sep = std::isspace(static_cast<unsigned char>(c)) || (comma && c == ',');
    if (sep) { if (!cur.empty()) { t.push_back(cur); cur.clear(); } }
    else cur += c;
  }
  if (!cur.empty()) t.push_back(cur);
  return t;
}

// ---- IRAP classic ASCII --------------------------------------------------------------------
std::string write_irap_ascii(const Grid2D& g, int decimals) {
  const auto& s = g.spec;
  std::string o;
  o += "-996 " + std::to_string(s.ny) + " " + fmt_fixed(s.dx, 6) + " " + fmt_fixed(s.dy, 6) + "\n";
  o += fmt_fixed(s.x0, 6) + " " + fmt_fixed(s.x(s.nx - 1), 6) + " " + fmt_fixed(s.y0, 6) + " " + fmt_fixed(s.y(s.ny - 1), 6) + "\n";
  o += std::to_string(s.nx) + " 0.000000 " + fmt_fixed(s.x0, 6) + " " + fmt_fixed(s.y0, 6) + "\n";
  o += "0 0 0 0 0 0 0\n";
  std::size_t k = 0;
  for (std::size_t j = 0; j < s.ny; ++j)
    for (std::size_t i = 0; i < s.nx; ++i) {
      double v = g.at(i, j);
      o += std::isnan(v) ? fmt_fixed(kIrapUndef, decimals) : fmt_fixed(v, decimals);
      o += (++k % 6 == 0) ? "\n" : " ";
    }
  if (k % 6 != 0) o.back() = '\n';
  return o;
}

Grid2D read_irap_ascii(const std::string& text) {
  auto t = tokens(text);
  if (t.size() < 19 || t[0] != "-996") throw Error(ErrorCode::Parse, "not an IRAP classic ASCII grid (expected -996 header)");
  GridSpec s;
  s.ny = static_cast<std::size_t>(parse_double(t[1]));
  s.dx = parse_double(t[2]);
  s.dy = parse_double(t[3]);
  s.x0 = parse_double(t[4]);
  s.y0 = parse_double(t[6]);
  s.nx = static_cast<std::size_t>(parse_double(t[8]));
  double rot = parse_double(t[9]);
  if (std::fabs(rot) > 1e-9) throw Error(ErrorCode::Unsupported, "rotated IRAP grids are not supported in v1");
  if (t.size() < 19 + s.nx * s.ny) throw Error(ErrorCode::Parse, "IRAP grid has fewer values than ncol x nrow");
  Grid2D g(s);
  for (std::size_t k = 0; k < s.size(); ++k) {
    double v = parse_double(t[19 + k]);
    g.v[k] = v >= kIrapUndef - 1 ? kNaN : v;
  }
  return g;
}

// ---- ZMAP+ ---------------------------------------------------------------------------------
std::string write_zmap(const Grid2D& g, const std::string& name, const std::vector<std::string>& comments, int decimals) {
  const auto& s = g.spec;
  const int per_line = 5, width = 20;
  std::string o;
  for (auto& c : comments) o += "! " + c + "\n";
  o += "@" + name + ", GRID, " + std::to_string(per_line) + "\n";
  o += std::to_string(width) + ", " + fmt_fixed(kZmapNull, 1) + ", , " + std::to_string(decimals) + ", 1\n";
  o += std::to_string(s.ny) + ", " + std::to_string(s.nx) + ", " + fmt_fixed(s.x0, 4) + ", " + fmt_fixed(s.x(s.nx - 1), 4) + ", " +
       fmt_fixed(s.y0, 4) + ", " + fmt_fixed(s.y(s.ny - 1), 4) + "\n";
  o += "0.0, 0.0, 0.0\n@\n";
  for (std::size_t i = 0; i < s.nx; ++i) {
    int k = 0;
    for (std::size_t jj = 0; jj < s.ny; ++jj) {
      std::size_t j = s.ny - 1 - jj;  // from the top (ymax) down
      double v = g.at(i, j);
      std::string f = std::isnan(v) ? fmt_fixed(kZmapNull, decimals) : fmt_fixed(v, decimals);
      if (static_cast<int>(f.size()) < width) f = std::string(static_cast<std::size_t>(width) - f.size(), ' ') + f;
      o += f;
      if (++k % per_line == 0) o += "\n";
    }
    if (k % per_line != 0) o += "\n";
  }
  return o;
}

Grid2D read_zmap(const std::string& text) {
  std::istringstream in(text);
  std::string line;
  std::vector<std::string> header;
  int ats = 0;
  std::string data;
  while (std::getline(in, line)) {
    std::string tl = trim(line);
    if (tl.empty() || tl[0] == '!') continue;
    if (ats < 2) {
      if (tl[0] == '@') { ++ats; if (ats == 1) header.push_back(tl.substr(1)); continue; }
      header.push_back(tl);
    } else {
      data += line + "\n";
    }
  }
  if (ats < 2 || header.size() < 4) throw Error(ErrorCode::Parse, "not a ZMAP+ grid (header between '@' lines expected)");
  auto h1 = split(header[1], ',');
  double nullv = parse_double(h1.at(1));
  auto h2 = split(header[2], ',');
  GridSpec s;
  s.ny = static_cast<std::size_t>(parse_double(h2.at(0)));
  s.nx = static_cast<std::size_t>(parse_double(h2.at(1)));
  s.x0 = parse_double(h2.at(2));
  double xmax = parse_double(h2.at(3));
  s.y0 = parse_double(h2.at(4));
  double ymax = parse_double(h2.at(5));
  s.dx = s.nx > 1 ? (xmax - s.x0) / (s.nx - 1) : 1;
  s.dy = s.ny > 1 ? (ymax - s.y0) / (s.ny - 1) : 1;
  auto t = tokens(data);
  if (t.size() < s.size()) throw Error(ErrorCode::Parse, "ZMAP+ grid has fewer values than rows x columns");
  Grid2D g(s);
  std::size_t k = 0;
  for (std::size_t i = 0; i < s.nx; ++i)
    for (std::size_t jj = 0; jj < s.ny; ++jj) {
      double v = parse_double(t[k++]);
      g.at(i, s.ny - 1 - jj) = std::fabs(v - nullv) < 1e-6 ? kNaN : v;
    }
  return g;
}

// ---- XYZ -----------------------------------------------------------------------------------
Grid2D read_xyz_grid(const std::string& text) {
  std::istringstream in(text);
  std::string line;
  std::vector<Vec3> p;
  while (std::getline(in, line)) {
    auto t = tokens(line, true);
    if (t.size() < 3) continue;
    try {
      p.push_back({parse_double(t[0]), parse_double(t[1]), parse_double(t[2])});
    } catch (const Error&) {
      continue;  // header or comment line
    }
  }
  if (p.size() < 4) throw Error(ErrorCode::Parse, "XYZ file has fewer than four points");
  std::set<double> xs, ys;
  for (auto& q : p) { xs.insert(q.x); ys.insert(q.y); }
  std::vector<double> X(xs.begin(), xs.end()), Y(ys.begin(), ys.end());
  auto regular = [](const std::vector<double>& v) {
    if (v.size() < 2) return false;
    double d = v[1] - v[0];
    for (std::size_t i = 2; i < v.size(); ++i)
      if (std::fabs(v[i] - v[i - 1] - d) > 1e-6 * std::max(1.0, std::fabs(d))) return false;
    return true;
  };
  if (!regular(X) || !regular(Y)) throw Error(ErrorCode::Unsupported, "XYZ input must be a regular grid in v1");
  GridSpec s{X[0], Y[0], X[1] - X[0], Y[1] - Y[0], X.size(), Y.size()};
  Grid2D g(s);
  for (auto& q : p) {
    auto i = static_cast<std::size_t>(std::llround((q.x - s.x0) / s.dx));
    auto j = static_cast<std::size_t>(std::llround((q.y - s.y0) / s.dy));
    g.at(i, j) = q.z;
  }
  return g;
}

// ---- importer plugins ----------------------------------------------------------------------
namespace {
template <Grid2D (*Fn)(const std::string&)>
class GridImporter final : public Importer {
 public:
  std::string data_kind() const override { return "grid"; }
  ImportedData read(const std::filesystem::path& f) const override {
    ImportedData d;
    d.kind = "grid";
    d.grid = Fn(read_file(f));
    return d;
  }
};
}  // namespace

STRATONIX_REGISTER(imp_irap_ascii, registry::Kind::Importer, "irap_ascii", "1.0.0",
                   "IRAP classic ASCII grid (unrotated)", [] { return std::make_unique<GridImporter<read_irap_ascii>>(); })
STRATONIX_REGISTER(imp_zmap, registry::Kind::Importer, "zmap", "1.0.0", "ZMAP+ ASCII grid",
                   [] { return std::make_unique<GridImporter<read_zmap>>(); })
STRATONIX_REGISTER(imp_xyz, registry::Kind::Importer, "xyz", "1.0.0", "XYZ points on a regular grid",
                   [] { return std::make_unique<GridImporter<read_xyz_grid>>(); })

void load_build_inputs(const desc::Description& d, BuildInputs& in) {
  auto& reg = registry::Registry::global();
  auto load = [&](const std::string& id, const std::string& fmt) {
    if (id.empty() || in.grids.count(id)) return;
    const desc::InputDesc* inp = nullptr;
    for (auto& x : d.inputs)
      if (x.id == id) inp = &x;
    if (!inp) throw Error(ErrorCode::Build, "input '" + id + "' is not listed under inputs");
    auto imp = reg.create<Importer>(registry::Kind::Importer, fmt);
    auto data = imp->read(inp->resolved);
    if (!data.grid) throw Error(ErrorCode::Build, "input '" + id + "' did not contain a grid");
    in.grids[id] = std::move(*data.grid);
  };
  for (auto& h : d.horizons) {
    load(h.surface_input, h.surface_format);
    load(h.iso_input, h.iso_format);
  }
  auto input_of = [&](const std::string& id) {
    for (auto& x : d.inputs)
      if (x.id == id) return &x;
    throw Error(ErrorCode::Build, "input '" + id + "' is not listed under inputs");
  };
  // Data for the implicit builder and for hybrid conditioning: points and orientations.
  std::vector<desc::ImplicitData> all_data = d.implicit.data;
  for (auto& x : d.conditioning.data)
    if (std::none_of(all_data.begin(), all_data.end(), [&](const desc::ImplicitData& y) { return y.input == x.input; })) all_data.push_back(x);
  // stage 4A: gravity and magnetic stations
  for (auto* g : {&d.geophysics.gravity, &d.geophysics.magnetics})
    if (!g->input.empty() && std::none_of(all_data.begin(), all_data.end(), [&](const desc::ImplicitData& y) { return y.input == g->input; }))
      all_data.push_back({g->input, "stations", g->format.empty() ? std::string("csv_stations") : g->format, "", desc::Param{}});
  for (auto& x : all_data) {
    auto data = reg.create<Importer>(registry::Kind::Importer, x.format)->read(input_of(x.input)->resolved);
    if (data.z_is == "depth")
      for (auto& p : data.points) p.z = -p.z;
    auto& out = in.data[x.input];
    for (std::size_t i = 0; i < data.points.size(); ++i) {
      DataPoint dp;
      dp.p = data.points[i];
      if (i < data.normals.size()) { dp.n = data.normals[i]; dp.oriented = true; }
      dp.label = !x.horizon.empty() ? x.horizon : (i < data.labels.size() ? data.labels[i] : std::string());
      if (i < data.sigmas.size()) dp.sigma = data.sigmas[i];
      if (i < data.values.size()) dp.value = data.values[i];
      out.push_back(dp);
    }
    if (out.empty()) throw Error(ErrorCode::Build, "input '" + x.input + "' holds no data points");
  }
  // Fault sticks for faults given by sticks.
  for (auto& f : d.faults) {
    if (f.geom_type != "sticks" || in.lines.count(f.geom_input)) continue;
    auto data = reg.create<Importer>(registry::Kind::Importer, f.geom_input_format)->read(input_of(f.geom_input)->resolved);
    std::vector<Polyline> ls;
    for (auto& l : data.lines)
      if (l.name == f.name || data.lines.size() == 1 || f.geom_input_format == "irap_lines") ls.push_back(l.pts);
    if (ls.empty()) throw Error(ErrorCode::Build, "input '" + f.geom_input + "' holds no sticks for fault " + f.name);
    in.lines[f.geom_input] = std::move(ls);
  }
}

void register_grid_importers() {
  STRATONIX_PLUGIN_ANCHOR(imp_irap_ascii);
  STRATONIX_PLUGIN_ANCHOR(imp_zmap);
  STRATONIX_PLUGIN_ANCHOR(imp_xyz);
}

}  // namespace stratonix::io
