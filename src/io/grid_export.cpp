// Grid exporters (stage 3B, STX-TN-12 section 6): GRDECL, VTK unstructured grids of the cells
// and of the tetrahedra, and RESQML 2.0.1 through FESAPI where this build has it.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <ostream>
#include <tuple>

#include "stratonix/core/crs.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/io/export.hpp"

STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_grdecl);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_grid_vtk);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_resqml);

namespace stratonix::io {

namespace fs = std::filesystem;

#ifdef STRATONIX_HAVE_FESAPI
std::vector<ExportedFile> write_resqml_grid(ExportContext& c, const fs::path& folder);  // resqml_fesapi.cpp
#endif

namespace {

/// Output buffered in memory and flushed to a stream every megabyte (stage 4D, STX-TN-17
/// section 4.2: a grid's text is never held whole).
struct Out {
  std::ostream* os = nullptr;
  std::string buf;
  Out& operator+=(const std::string& t) {
    buf += t;
    if (os && buf.size() > (1u << 20)) flush();
    return *this;
  }
  Out& operator+=(const char* t) { return *this += std::string(t); }
  void flush() {
    if (os) { os->write(buf.data(), static_cast<std::streamsize>(buf.size())); buf.clear(); }
  }
};
template <class F>
void put_values(Out& o, std::size_t count, F value, int dec, std::size_t per_line) {
  for (std::size_t n = 0; n < count; ++n) {
    o += fmt_fixed(value(n), dec);
    o += (n + 1) % per_line == 0 || n + 1 == count ? "\n" : " ";
  }
}
/// Integers with run-length repeats (n*v), as Eclipse reads them.
template <class F>
void put_ints(Out& o, std::size_t count, F value) {
  std::size_t n = 0, on_line = 0;
  while (n < count) {
    const long long vn = value(n);
    std::size_t r = n;
    while (r < count && value(r) == vn) ++r;
    std::size_t len = r - n;
    o += len > 1 ? std::to_string(len) + "*" + std::to_string(vn) : std::to_string(vn);
    n = r;
    o += ++on_line % 12 == 0 || n == count ? "\n" : " ";
  }
}

void grdecl_to(Out& o, const repr::CornerPointGrid& g, const repr::PropertySet* props, const repr::TransmissibilityResult* trans, const GrdeclInfo& info,
               int dec);

}  // namespace

namespace {
void grdecl_to(Out& o, const repr::CornerPointGrid& g, const repr::PropertySet* props, const repr::TransmissibilityResult* trans, const GrdeclInfo& info,
               int dec) {
  for (auto& h : info.header) o += "-- " + h + "\n";
  o += "-- Stratonix corner-point grid: depth positive down (m) = -elevation; i fastest, then j, then k (k = 1 at the top).\n";
  o += "-- Real properties hold " + fmt_fixed(kGrdeclMissing, 2) + " where a value could not be computed (never defaulted); see the reasons below.\n\n";
  o += "MAPUNITS\n  METRES /\n\nGRIDUNIT\n  METRES  /\n\n";
  o += "SPECGRID\n  " + std::to_string(g.nx) + " " + std::to_string(g.ny) + " " + std::to_string(g.nz) + " 1 F /\n\n";
  o += "COORD\n";
  put_values(o, g.coord.size(), [&](std::size_t n) { return n % 3 == 2 ? -g.coord[n] : g.coord[n]; }, dec, 6);
  o += "/\n\nZCORN\n";
  put_values(o, g.zcorn.size(), [&](std::size_t n) { return -g.zcorn[n]; }, dec, 8);
  o += "/\n\nACTNUM\n";
  put_ints(o, g.actnum.size(), [&](std::size_t n) { return static_cast<long long>(g.actnum[n]); });
  o += "/\n";
  auto write_prop = [&](const repr::CellProperty& p, const std::string& kw) {
    o += "\n-- " + kw + ": " + p.unit + "; evidence " + p.evidence + "; source " + p.source + "\n";
    if (p.missing) o += "-- " + std::to_string(p.missing) + " cells missing: " + p.missing_reason + "\n";
    o += kw + "\n";
    if (p.integer) {
      put_ints(o, p.v.size(), [&](std::size_t n) { return std::isnan(p.v[n]) ? 0LL : static_cast<long long>(std::llround(p.v[n])); });
    } else {
      put_values(o, p.v.size(), [&](std::size_t n) { return std::isnan(p.v[n]) ? kGrdeclMissing : p.v[n]; }, 6, 8);
    }
    o += "/\n";
  };
  if (props)
    for (auto& p : props->props) write_prop(p, p.name);
  if (trans) {
    write_prop(trans->multx, "MULTX");
    write_prop(trans->multy, "MULTY");
  }
  if (info.faces && !info.faces->empty()) {
    // FAULTS: per fault and face, consecutive layers in one line
    o += "\n-- FAULTS: faces between cells of different fault blocks inside the fault's tip line\nFAULTS\n";
    std::map<std::tuple<int, char, int, int>, std::vector<int>> runs;
    for (auto& f : *info.faces) runs[{f.fault, f.dir, f.j, f.i}].push_back(f.k);
    for (auto& [key, ks] : runs) {
      auto [fault, dir, j, i] = key;
      std::sort(ks.begin(), ks.end());
      std::string name = fault < static_cast<int>(info.fault_names.size()) ? petrel_name(info.fault_names[static_cast<std::size_t>(fault)])
                                                                          : "F" + std::to_string(fault + 1);
      for (std::size_t a = 0; a < ks.size();) {
        std::size_t b = a;
        while (b + 1 < ks.size() && ks[b + 1] == ks[b] + 1) ++b;
        o += "  '" + name + "' " + std::to_string(i + 1) + " " + std::to_string(i + 1) + " " + std::to_string(j + 1) + " " + std::to_string(j + 1) + " " +
             std::to_string(ks[a] + 1) + " " + std::to_string(ks[b] + 1) + " " + dir + " /\n";
        a = b + 1;
      }
    }
    o += "/\n";
  }
}
}  // namespace

std::string write_grdecl(const repr::CornerPointGrid& g, const repr::PropertySet* props, const repr::TransmissibilityResult* trans,
                         const GrdeclInfo& info, int dec) {
  Out o;
  grdecl_to(o, g, props, trans, info, dec);
  return std::move(o.buf);
}
void write_grdecl(std::ostream& os, const repr::CornerPointGrid& g, const repr::PropertySet* props, const repr::TransmissibilityResult* trans,
                  const GrdeclInfo& info, int dec) {
  Out o;
  o.os = &os;
  grdecl_to(o, g, props, trans, info, dec);
  o.flush();
}

namespace {
/// Cell data from a property set: integers as Int32, the rest as Float64, for the cells `ids`.
void add_cell_data(VtuData& d, const repr::CornerPointGrid& g, const repr::PropertySet* props, const std::vector<std::size_t>& ids) {
  std::vector<std::int32_t> zone, block;
  for (auto c : ids) { zone.push_back(g.zone[c]); block.push_back(static_cast<std::int32_t>(g.block[c])); }
  d.cell_data.push_back({"ZONE", zone});
  d.cell_data.push_back({"BLOCK", block});
  if (!props) return;
  for (auto& p : props->props) {
    if (p.name == "ZONE" || p.name == "BLOCK") continue;
    if (p.integer) {
      std::vector<std::int32_t> v;
      for (auto c : ids) v.push_back(std::isnan(p.v[c]) ? -1 : static_cast<std::int32_t>(std::llround(p.v[c])));
      d.cell_data.push_back({p.name, v});
    } else {
      std::vector<double> v;
      for (auto c : ids) v.push_back(p.v[c]);
      d.cell_data_real.push_back({p.name, v});
    }
  }
}
}  // namespace

VtuData grid_vtu(const repr::CornerPointGrid& g, const repr::PropertySet* props) {
  VtuData d;
  std::map<std::pair<std::size_t, long long>, std::int64_t> node;
  auto id = [&](int i, int j, int k, int c) {
    std::size_t p = g.pillar(i + (c & 1), j + ((c >> 1) & 1));
    double z = g.zcorn[g.zcorn_index(i, j, k, c)];
    auto key = std::make_pair(p, std::llround(z * 1000.0));
    auto it = node.find(key);
    if (it != node.end()) return it->second;
    d.pts.push_back(g.corner(i, j, k, c));
    return node[key] = static_cast<std::int64_t>(d.pts.size() - 1);
  };
  std::vector<std::size_t> ids;
  static const int order[8] = {4, 5, 7, 6, 0, 1, 3, 2};  // VTK hexahedron: bottom face, then top face
  for (int k = 0; k < g.nz; ++k)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i) {
        auto c = g.cell(i, j, k);
        if (!g.actnum[c]) continue;
        std::vector<std::int64_t> cell;
        for (int q : order) cell.push_back(id(i, j, k, q));
        d.cells.push_back(cell);
        d.types.push_back(12);
        ids.push_back(c);
      }
  add_cell_data(d, g, props, ids);
  return d;
}

VtuData tet_vtu(const repr::TetGrid& t, const repr::CornerPointGrid& g, const repr::PropertySet* props) {
  VtuData d;
  d.pts = t.nodes;
  std::vector<std::size_t> ids;
  for (std::size_t n = 0; n < t.tets.size(); ++n) {
    d.cells.push_back({t.tets[n][0], t.tets[n][1], t.tets[n][2], t.tets[n][3]});
    d.types.push_back(10);
    ids.push_back(static_cast<std::size_t>(t.cell[n]));
  }
  add_cell_data(d, g, props, ids);
  return d;
}

namespace {

void put(const fs::path& folder, const fs::path& rel, const std::string& text) {
  fs::create_directories((folder / rel).parent_path());
  write_file(folder / rel, text);
}
const repr::CornerPointGrid& need_grid(const ExportContext& c, const char* fmt) {
  if (!c.cells.grid)
    throw Error(ErrorCode::InvalidArgument, std::string("the ") + fmt + " export needs a grid block in the description (STX-TN-12 section 10)", "grid");
  return *c.cells.grid;
}

class GrdeclExporter final : public Exporter {
 public:
  std::string id() const override { return "grdecl"; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    const auto& g = need_grid(c, "grdecl");
    GrdeclInfo info;
    info.header = c.comments();
    for (auto& f : c.info().faults) info.fault_names.push_back(f.name);
    info.faces = c.cells.faces;
    std::string name = petrel_name(c.info().name) + "_grid.grdecl";
    {  // written as it is made (STX-TN-17 section 4.2)
      fs::create_directories(folder / "grid");
      std::ofstream os(folder / "grid" / name, std::ios::binary);
      if (!os) throw Error(ErrorCode::Io, "cannot write " + (folder / "grid" / name).string());
      write_grdecl(os, g, c.cells.props, c.cells.trans, info, std::max(2, c.decimals - 1));
    }
    std::size_t act = 0;
    for (auto a : g.actnum) act += a;
    return {{"grid/" + name, "Eclipse GRDECL (corner-point grid)", "Eclipse grid (*.GRDECL)",
             std::to_string(g.nx) + " x " + std::to_string(g.ny) + " x " + std::to_string(g.nz) + " cells, " + std::to_string(act) +
                 " active; ZONE, BLOCK and the grid properties; MULTX/MULTY and FAULTS when computed"}};
  }
};

class GridVtkExporter final : public Exporter {
 public:
  std::string id() const override { return "grid_vtk"; }
  bool petrel() const override { return false; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    const auto& g = need_grid(c, "grid_vtk");
    std::vector<ExportedFile> out;
    put(folder, "grid/grid_cells.vtu", write_vtu(grid_vtu(g, c.cells.props), c.decimals));
    out.push_back({"grid/grid_cells.vtu", "VTK XML unstructured grid", "-", "active hexahedral cells (type 12) with ZONE, BLOCK and the properties"});
    if (c.cells.tets) {
      put(folder, "grid/grid_tetrahedra.vtu", write_vtu(tet_vtu(*c.cells.tets, g, c.cells.props), c.decimals));
      out.push_back({"grid/grid_tetrahedra.vtu", "VTK XML unstructured grid", "-",
                     std::to_string(c.cells.tets->tets.size()) + " tetrahedra (type 10) with the cell's ZONE, BLOCK and properties"});
    }
    return out;
  }
};

class ResqmlExporter final : public Exporter {
 public:
  std::string id() const override { return "resqml"; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    need_grid(c, "resqml");
#ifdef STRATONIX_HAVE_FESAPI
    return write_resqml_grid(c, folder);
#else
    (void)folder;
    throw Error(ErrorCode::InvalidArgument,
                "RESQML export needs a build with FESAPI (STRATONIX_WITH_FESAPI); this build has none. Use grdecl, or the release built by the CI job 'resqml'",
                "outputs.exports");
#endif
  }
};

}  // namespace

STRATONIX_REGISTER(exp_grdecl, registry::Kind::Exporter, "grdecl", "1.0.0", "Eclipse GRDECL corner-point grid with properties and faults",
                   [] { return std::make_unique<GrdeclExporter>(); })
STRATONIX_REGISTER(exp_grid_vtk, registry::Kind::Exporter, "grid_vtk", "1.0.0", "VTK unstructured grids of the grid cells and tetrahedra",
                   [] { return std::make_unique<GridVtkExporter>(); })
STRATONIX_REGISTER(exp_resqml, registry::Kind::Exporter, "resqml", "1.0.0", "RESQML 2.0.1 IJK grid with properties (needs FESAPI)",
                   [] { return std::make_unique<ResqmlExporter>(); })

void register_grid_exporters() {
  STRATONIX_PLUGIN_ANCHOR(exp_grdecl);
  STRATONIX_PLUGIN_ANCHOR(exp_grid_vtk);
  STRATONIX_PLUGIN_ANCHOR(exp_resqml);
}

}  // namespace stratonix::io

namespace stratonix::io {
std::vector<std::string> optional_components() {
  std::vector<std::string> c;
#ifdef STRATONIX_HAVE_FESAPI
  c.push_back("RESQML 2.0.1 export (FESAPI): included");
#else
  c.push_back("RESQML 2.0.1 export (FESAPI): not included");
#endif
  c.push_back(std::string("coordinate transformations (PROJ): ") + (crs_reprojection_available() ? "included" : "not included"));
#ifdef STRATONIX_HAVE_SQLITE
  c.push_back("GeoPackage maps (SQLite): included");
#else
  c.push_back("GeoPackage maps (SQLite): not included");
#endif
  c.push_back("LAS 2.0/3.0 and DLIS (RP66 V1) well data: included");
  c.push_back("SEG-Y revisions 0 to 2.1: included");
  return c;
}
}  // namespace stratonix::io
