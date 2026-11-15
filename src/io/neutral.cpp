// Neutral formats (VTK XML, OBJ, STL, CSV/XYZ, GeoJSON) and the GemPy / LoopStructural input sets.
#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/io/export.hpp"
#include "stratonix/io/io.hpp"

// Plugin anchors called below (C linkage, defined by STRATONIX_REGISTER in each plugin file).
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_vtk);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_obj);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_stl);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_csv);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_geojson);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_gempy);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_loop);

namespace stratonix::io {

namespace fs = std::filesystem;

namespace {
std::string num(double v, int d) { return fmt_fixed(v, d); }

/// Text between <tag ...> and </tag> for the first DataArray whose attributes contain `key`.
std::string data_array(const std::string& xml, const std::string& key, std::size_t from = 0, std::size_t to = std::string::npos) {
  std::size_t p = from;
  for (;;) {
    p = xml.find("<DataArray", p);
    if (p == std::string::npos || p >= to) return {};
    std::size_t e = xml.find('>', p);
    std::string head = xml.substr(p, e - p);
    if (head.find(key) != std::string::npos) {
      std::size_t c = xml.find("</DataArray>", e);
      return xml.substr(e + 1, c - e - 1);
    }
    p = e;
  }
}
std::string attr(const std::string& xml, const std::string& name, std::size_t from = 0) {
  std::size_t p = xml.find(name + "=\"", from);
  if (p == std::string::npos) return {};
  p += name.size() + 2;
  return xml.substr(p, xml.find('"', p) - p);
}
std::vector<double> numbers(const std::string& s) {
  std::vector<double> v;
  std::istringstream in(s);
  std::string t;
  while (in >> t) v.push_back(parse_double(t));
  return v;
}
std::vector<std::string> lines_of(const std::string& text) {
  std::vector<std::string> out;
  std::istringstream in(text);
  std::string l;
  while (std::getline(in, l)) {
    if (!l.empty() && l.back() == '\r') l.pop_back();
    out.push_back(l);
  }
  return out;
}
template <class T>
std::string join_nums(const std::vector<T>& v, int d, std::size_t per_line = 12) {
  std::string o;
  for (std::size_t i = 0; i < v.size(); ++i) {
    if constexpr (std::is_floating_point_v<T>) o += num(v[i], d);
    else o += std::to_string(v[i]);
    o += ((i + 1) % per_line == 0 || i + 1 == v.size()) ? "\n" : " ";
  }
  return o;
}
}  // namespace

// ---- VTK XML -------------------------------------------------------------------------------
std::string write_vtu(const VtuData& d, int dec) {
  std::string o = "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  o += "<UnstructuredGrid>\n<Piece NumberOfPoints=\"" + std::to_string(d.pts.size()) + "\" NumberOfCells=\"" + std::to_string(d.cells.size()) + "\">\n";
  o += "<PointData>\n";
  for (auto& [n, v] : d.point_data) o += "<DataArray type=\"Float64\" Name=\"" + n + "\" format=\"ascii\">\n" + join_nums(v, dec) + "</DataArray>\n";
  o += "</PointData>\n<CellData>\n";
  for (auto& [n, v] : d.cell_data) o += "<DataArray type=\"Int32\" Name=\"" + n + "\" format=\"ascii\">\n" + join_nums(v, 0) + "</DataArray>\n";
  o += "</CellData>\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (auto& p : d.pts) o += num(p.x, dec) + " " + num(p.y, dec) + " " + num(p.z, dec) + "\n";
  o += "</DataArray>\n</Points>\n<Cells>\n";
  std::vector<std::int64_t> conn, off;
  std::vector<int> types;
  for (std::size_t c = 0; c < d.cells.size(); ++c) {
    conn.insert(conn.end(), d.cells[c].begin(), d.cells[c].end());
    off.push_back(static_cast<std::int64_t>(conn.size()));
    types.push_back(d.types[c]);
  }
  o += "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n" + join_nums(conn, 0) + "</DataArray>\n";
  o += "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n" + join_nums(off, 0) + "</DataArray>\n";
  o += "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n" + join_nums(types, 0) + "</DataArray>\n";
  o += "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  return o;
}

VtuData read_vtu(const std::string& x) {
  if (x.find("type=\"UnstructuredGrid\"") == std::string::npos) throw Error(ErrorCode::Parse, "not a VTK XML unstructured grid");
  VtuData d;
  std::size_t np = static_cast<std::size_t>(parse_double(attr(x, "NumberOfPoints")));
  std::size_t pts_at = x.find("<Points>");
  auto p = numbers(data_array(x, "NumberOfComponents=\"3\"", pts_at));
  if (p.size() != 3 * np) throw Error(ErrorCode::Parse, "VTU: point count mismatch");
  for (std::size_t i = 0; i < np; ++i) d.pts.push_back({p[3 * i], p[3 * i + 1], p[3 * i + 2]});
  auto conn = numbers(data_array(x, "Name=\"connectivity\""));
  auto off = numbers(data_array(x, "Name=\"offsets\""));
  auto ty = numbers(data_array(x, "Name=\"types\""));
  std::size_t s = 0;
  for (std::size_t c = 0; c < off.size(); ++c) {
    std::vector<std::int64_t> cell;
    for (auto k = s; k < static_cast<std::size_t>(off[c]); ++k) cell.push_back(static_cast<std::int64_t>(conn[k]));
    d.cells.push_back(cell);
    d.types.push_back(static_cast<std::uint8_t>(ty.at(c)));
    s = static_cast<std::size_t>(off[c]);
  }
  std::size_t pd = x.find("<PointData>"), pde = x.find("</PointData>");
  for (std::size_t q = x.find("<DataArray", pd); q != std::string::npos && q < pde; q = x.find("<DataArray", q + 1)) {
    std::string name = attr(x, "Name", q);
    d.point_data.push_back({name, numbers(data_array(x, "Name=\"" + name + "\"", q))});
  }
  std::size_t cd = x.find("<CellData>"), cde = x.find("</CellData>");
  for (std::size_t q = x.find("<DataArray", cd); q != std::string::npos && q < cde; q = x.find("<DataArray", q + 1)) {
    std::string name = attr(x, "Name", q);
    std::vector<std::int32_t> v;
    for (double t : numbers(data_array(x, "Name=\"" + name + "\"", q))) v.push_back(static_cast<std::int32_t>(t));
    d.cell_data.push_back({name, v});
  }
  return d;
}

std::string write_vti(const repr::VoxelGrid& v, int dec) {
  std::string ext = "0 " + std::to_string(v.nx) + " 0 " + std::to_string(v.ny) + " 0 " + std::to_string(v.nz);
  std::string o = "<?xml version=\"1.0\"?>\n<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  o += "<ImageData WholeExtent=\"" + ext + "\" Origin=\"" + num(v.x0 - v.dx / 2, dec) + " " + num(v.y0 - v.dy / 2, dec) + " " + num(v.z0 - v.dz / 2, dec) +
       "\" Spacing=\"" + num(v.dx, dec) + " " + num(v.dy, dec) + " " + num(v.dz, dec) + "\">\n";
  o += "<Piece Extent=\"" + ext + "\">\n<CellData Scalars=\"unit\">\n";
  o += "<DataArray type=\"Int32\" Name=\"unit\" format=\"ascii\">\n" + join_nums(v.unit, 0, 40) + "</DataArray>\n";
  std::vector<std::int64_t> b(v.block.begin(), v.block.end());
  if (!b.empty()) o += "<DataArray type=\"Int64\" Name=\"block\" format=\"ascii\">\n" + join_nums(b, 0, 40) + "</DataArray>\n";
  o += "</CellData>\n</Piece>\n</ImageData>\n</VTKFile>\n";
  return o;
}

repr::VoxelGrid read_vti(const std::string& x) {
  if (x.find("type=\"ImageData\"") == std::string::npos) throw Error(ErrorCode::Parse, "not a VTK XML image data file");
  repr::VoxelGrid v;
  auto e = numbers(attr(x, "WholeExtent"));
  auto o = numbers(attr(x, "Origin"));
  auto s = numbers(attr(x, "Spacing"));
  v.nx = static_cast<std::size_t>(e.at(1));
  v.ny = static_cast<std::size_t>(e.at(3));
  v.nz = static_cast<std::size_t>(e.at(5));
  v.dx = s.at(0); v.dy = s.at(1); v.dz = s.at(2);
  v.x0 = o.at(0) + v.dx / 2; v.y0 = o.at(1) + v.dy / 2; v.z0 = o.at(2) + v.dz / 2;
  for (double t : numbers(data_array(x, "Name=\"unit\""))) v.unit.push_back(static_cast<std::int32_t>(t));
  for (double t : numbers(data_array(x, "Name=\"block\""))) v.block.push_back(static_cast<BlockCode>(t));
  if (v.unit.size() != v.nx * v.ny * v.nz) throw Error(ErrorCode::Parse, "VTI: cell count mismatch");
  return v;
}

// ---- OBJ / STL -----------------------------------------------------------------------------
std::string write_obj(const TriMesh& m, const std::string& name, const std::vector<std::string>& comments, int d) {
  std::string o;
  for (auto& c : comments) o += "# " + c + "\n";
  o += "o " + name + "\n";
  for (auto& p : m.v) o += "v " + num(p.x, d) + " " + num(p.y, d) + " " + num(p.z, d) + "\n";
  for (auto& t : m.t) o += "f " + std::to_string(t[0] + 1) + " " + std::to_string(t[1] + 1) + " " + std::to_string(t[2] + 1) + "\n";
  return o;
}
TriMesh read_obj(const std::string& text) {
  TriMesh m;
  for (auto& l : lines_of(text)) {
    std::istringstream in(l);
    std::string k;
    in >> k;
    if (k == "v") {
      std::string a, b, c;
      in >> a >> b >> c;
      m.v.push_back({parse_double(a), parse_double(b), parse_double(c)});
    } else if (k == "f") {
      std::array<std::uint32_t, 3> t{};
      for (auto& x : t) {
        std::string s;
        in >> s;
        x = static_cast<std::uint32_t>(std::stoul(s.substr(0, s.find('/'))) - 1);
      }
      m.t.push_back(t);
    }
  }
  return m;
}
std::string write_stl(const TriMesh& m, const std::string& name, int d) {
  std::string o = "solid " + name + "\n";
  for (auto& t : m.t) {
    Vec3 a = m.v[t[0]], b = m.v[t[1]], c = m.v[t[2]];
    Vec3 n = normalized(cross(b - a, c - a));
    o += " facet normal " + num(n.x, 6) + " " + num(n.y, 6) + " " + num(n.z, 6) + "\n  outer loop\n";
    for (auto& p : {a, b, c}) o += "   vertex " + num(p.x, d) + " " + num(p.y, d) + " " + num(p.z, d) + "\n";
    o += "  endloop\n endfacet\n";
  }
  o += "endsolid " + name + "\n";
  return o;
}
TriMesh read_stl(const std::string& text) {
  if (text.rfind("solid", 0) != 0) throw Error(ErrorCode::Parse, "not an ASCII STL file");
  TriMesh m;
  std::vector<Vec3> cur;
  for (auto& l : lines_of(text)) {
    std::istringstream in(l);
    std::string k;
    in >> k;
    if (k == "vertex") {
      std::string a, b, c;
      in >> a >> b >> c;
      cur.push_back({parse_double(a), parse_double(b), parse_double(c)});
      if (cur.size() == 3) {
        auto i = static_cast<std::uint32_t>(m.v.size());
        m.v.insert(m.v.end(), cur.begin(), cur.end());
        m.t.push_back({i, i + 1, i + 2});
        cur.clear();
      }
    }
  }
  return m;
}

// ---- CSV / XYZ -----------------------------------------------------------------------------
std::string write_csv(const CsvTable& t, const std::vector<std::string>& comments) {
  std::string o;
  for (auto& c : comments) o += "# " + c + "\n";
  for (std::size_t i = 0; i < t.header.size(); ++i) o += t.header[i] + (i + 1 < t.header.size() ? "," : "\n");
  for (auto& r : t.rows)
    for (std::size_t i = 0; i < r.size(); ++i) {
      if (r[i].find_first_of(",\"\n") != std::string::npos) throw Error(ErrorCode::InvalidArgument, "CSV value with a comma or quote: " + r[i]);
      o += r[i] + (i + 1 < r.size() ? "," : "\n");
    }
  return o;
}
CsvTable read_csv(const std::string& text) {
  CsvTable t;
  for (auto& l : lines_of(text)) {
    if (l.empty() || l[0] == '#') continue;
    auto f = split(l, ',');
    for (auto& x : f) x = trim(x);
    if (t.header.empty()) t.header = f;
    else t.rows.push_back(f);
  }
  return t;
}
std::vector<Vec3> read_xyz_points(const std::string& text) {
  std::vector<Vec3> p;
  for (auto& l : lines_of(text)) {
    std::string tl = trim(l);
    if (tl.empty() || tl[0] == '#') continue;
    std::istringstream in(tl);
    std::string a, b, c;
    if (!(in >> a >> b >> c)) throw Error(ErrorCode::Parse, "XYZ: expected three numbers: " + tl);
    p.push_back({parse_double(a), parse_double(b), parse_double(c)});
  }
  return p;
}

// ---- exporters -----------------------------------------------------------------------------
namespace {
void put(const fs::path& folder, const fs::path& rel, const std::string& text) {
  fs::create_directories((folder / rel).parent_path());
  write_file(folder / rel, text);
}
std::string g(const fs::path& p) { return p.generic_string(); }

VtuData mesh_vtu(const TriMesh& m, const std::string& prop, const std::vector<double>& vals, int id_value, const std::string& id_name) {
  VtuData d;
  d.pts = m.v;
  for (auto& t : m.t) {
    d.cells.push_back({t[0], t[1], t[2]});
    d.types.push_back(5);
  }
  if (!prop.empty()) d.point_data.push_back({prop, vals});
  d.cell_data.push_back({id_name, std::vector<std::int32_t>(m.t.size(), id_value)});
  return d;
}

class VtkExporter final : public Exporter {
 public:
  std::string id() const override { return "vtk"; }
  bool petrel() const override { return false; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    std::vector<ExportedFile> out;
    const auto& inf = c.info();
    std::vector<std::string> datasets;
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      TriMesh m = repr::grid_to_mesh(c.faulted(static_cast<int>(h)).z);
      std::vector<double> depth;
      for (auto& p : m.v) depth.push_back(-p.z);
      fs::path p = fs::path("vtk") / (petrel_name(inf.horizons[h].name) + ".vtu");
      put(folder, p, write_vtu(mesh_vtu(m, "depth_m", depth, static_cast<int>(h), "horizon"), c.decimals));
      out.push_back({g(p), "VTK XML unstructured grid", "-", inf.horizons[h].name + " faulted surface; point data depth_m (positive down)"});
      datasets.push_back(p.filename().string());
    }
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      const auto& m = c.fault_surface(static_cast<int>(f)).mesh;
      FaultEval fe;
      c.model().evaluate_fault(static_cast<int>(f), m.v, fe);
      fs::path p = fs::path("vtk") / (petrel_name(inf.faults[f].name) + ".vtu");
      auto vd = mesh_vtu(m, "displacement_m", fe.displacement, static_cast<int>(f), "fault");
      std::string extra;
      if (c.fault_rock || c.stress.present) {  // stage 2F attributes as further point data
        std::vector<repr::FaultRockPoint> pr;
        for (auto& q : m.v) pr.push_back(repr::fault_rock_at(c.model(), static_cast<int>(f), q, c.rock, c.stress));
        auto add = [&](const char* name, auto get, double scale) {
          std::vector<double> v;
          for (auto& r : pr) v.push_back(get(r) * scale);
          vd.point_data.push_back({name, v});
          extra += std::string(", ") + name;
        };
        add("footwall_unit", [](const repr::FaultRockPoint& r) { return static_cast<double>(r.fw_unit); }, 1.0);
        add("hanging_wall_unit", [](const repr::FaultRockPoint& r) { return static_cast<double>(r.hw_unit); }, 1.0);
        if (c.fault_rock) {
          add("sgr", [](const repr::FaultRockPoint& r) { return r.sgr; }, 1.0);
          add("csp", [](const repr::FaultRockPoint& r) { return r.csp; }, 1.0);
          add("fault_thickness_m", [](const repr::FaultRockPoint& r) { return r.thickness; }, 1.0);
          add("seal_capacity_bar", [](const repr::FaultRockPoint& r) { return r.seal_bar; }, 1.0);
        }
        if (c.stress.present) {
          add("slip_tendency", [](const repr::FaultRockPoint& r) { return r.ts; }, 1.0);
          add("dilation_tendency", [](const repr::FaultRockPoint& r) { return r.td; }, 1.0);
          add("reactivation_dp_MPa", [](const repr::FaultRockPoint& r) { return r.dp; }, 1e-6);
        }
      }
      put(folder, p, write_vtu(vd, c.decimals));
      out.push_back({g(p), "VTK XML unstructured grid", "-", inf.faults[f].name + " plane; point data displacement_m" + extra});
      datasets.push_back(p.filename().string());
    }
    if (!c.wells().empty()) {
      VtuData d;
      std::vector<std::int32_t> ids;
      for (std::size_t w = 0; w < c.wells().size(); ++w) {
        std::vector<std::int64_t> cell;
        for (auto& q : c.wells()[w].points()) {
          cell.push_back(static_cast<std::int64_t>(d.pts.size()));
          d.pts.push_back(q);
        }
        d.cells.push_back(cell);
        d.types.push_back(4);
        ids.push_back(static_cast<std::int32_t>(w));
      }
      d.cell_data.push_back({"well", ids});
      fs::path p = fs::path("vtk") / "wells.vtu";
      put(folder, p, write_vtu(d, c.decimals));
      out.push_back({g(p), "VTK XML unstructured grid", "-", "well paths as polylines (cell data well = index in README order)"});
      datasets.push_back(p.filename().string());
    }
    std::string vox;
    if (c.voxel_dx > 0) {
      auto v = repr::voxel_grid(c.model(), c.voxel_dx, c.voxel_dy, c.voxel_dz, true);
      fs::path p = fs::path("vtk") / "voxels.vti";
      put(folder, p, write_vti(v, c.decimals));
      out.push_back({g(p), "VTK XML image data", "-",
                     "voxel grid " + std::to_string(v.nx) + " x " + std::to_string(v.ny) + " x " + std::to_string(v.nz) + "; cell data unit, block"});
      vox = p.filename().string();
    }
    // Collection and a ParaView Python state.
    std::string pvd = "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n<Collection>\n";
    int part = 0;
    for (auto& n : datasets) pvd += "<DataSet timestep=\"0\" part=\"" + std::to_string(part++) + "\" file=\"" + n + "\"/>\n";
    if (!vox.empty()) pvd += "<DataSet timestep=\"0\" part=\"" + std::to_string(part++) + "\" file=\"" + vox + "\"/>\n";
    pvd += "</Collection>\n</VTKFile>\n";
    put(folder, "vtk/model.pvd", pvd);
    out.push_back({"vtk/model.pvd", "ParaView data collection", "-", "all VTK files above in one collection"});
    std::string py = "# ParaView state (Python). ParaView 5.9 or later: File > Load State, choose this file;\n"
                     "# or run: pvpython paraview_state.py. Written by Stratonix " + build_info().version + ".\n"
                     "# Z is elevation (m); vertical exaggeration " + fmt_shortest(c.vertical_exaggeration) + " is applied by a Transform.\n"
                     "import os\nfrom paraview.simple import *\nhere = os.path.dirname(os.path.abspath(__file__))\nVE = " +
                     fmt_shortest(c.vertical_exaggeration) + "\nview = GetActiveViewOrCreate('RenderView')\n";
    py += "surfaces = [";
    for (auto& n : datasets) py += "'" + n + "', ";
    py += "]\nfor name in surfaces:\n    r = XMLUnstructuredGridReader(registrationName=name, FileName=[os.path.join(here, name)])\n"
          "    t = Transform(registrationName=name + ' (VE)', Input=r)\n    t.Transform.Scale = [1.0, 1.0, VE]\n"
          "    d = Show(t, view)\n    if name.startswith('wells'):\n        d.LineWidth = 3.0\n        ColorBy(d, None)\n"
          "        d.AmbientColor = d.DiffuseColor = [0.0, 0.0, 0.0]\n    elif 'depth_m' in r.PointData.keys():\n        ColorBy(d, ('POINTS', 'depth_m'))\n"
          "    else:\n        ColorBy(d, None)\n        d.DiffuseColor = [0.8, 0.2, 0.2]\n";
    if (!vox.empty())
      py += "v = XMLImageDataReader(registrationName='" + vox + "', FileName=[os.path.join(here, '" + vox +
            "')])\nvt = Transform(registrationName='voxels (VE)', Input=v)\nvt.Transform.Scale = [1.0, 1.0, VE]\n"
            "vd = Show(vt, view)\nColorBy(vd, ('CELLS', 'unit'))\nHide(vt, view)  # switch on in the pipeline browser\n";
    py += "view.ResetCamera()\nRender()\n";
    put(folder, "vtk/paraview_state.py", py);
    out.push_back({"vtk/paraview_state.py", "ParaView Python state", "-", "loads every file with vertical exaggeration and colours"});
    return out;
  }
};

class ObjExporter final : public Exporter {
 public:
  std::string id() const override { return "obj"; }
  bool petrel() const override { return false; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    std::vector<ExportedFile> out;
    const auto& inf = c.info();
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      fs::path p = fs::path("obj") / (petrel_name(inf.horizons[h].name) + ".obj");
      put(folder, p, write_obj(repr::grid_to_mesh(c.faulted(static_cast<int>(h)).z), petrel_name(inf.horizons[h].name), c.comments(), c.decimals));
      out.push_back({g(p), "Wavefront OBJ", "-", inf.horizons[h].name + " faulted surface"});
    }
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      fs::path p = fs::path("obj") / (petrel_name(inf.faults[f].name) + ".obj");
      put(folder, p, write_obj(c.fault_surface(static_cast<int>(f)).mesh, petrel_name(inf.faults[f].name), c.comments(), c.decimals));
      out.push_back({g(p), "Wavefront OBJ", "-", inf.faults[f].name + " plane"});
    }
    return out;
  }
};

class StlExporter final : public Exporter {
 public:
  std::string id() const override { return "stl"; }
  bool petrel() const override { return false; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    std::vector<ExportedFile> out;
    const auto& inf = c.info();
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      fs::path p = fs::path("stl") / (petrel_name(inf.horizons[h].name) + ".stl");
      put(folder, p, write_stl(repr::grid_to_mesh(c.faulted(static_cast<int>(h)).z), petrel_name(inf.horizons[h].name), c.decimals));
      out.push_back({g(p), "STL (ASCII)", "-", inf.horizons[h].name + " faulted surface"});
    }
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      fs::path p = fs::path("stl") / (petrel_name(inf.faults[f].name) + ".stl");
      put(folder, p, write_stl(c.fault_surface(static_cast<int>(f)).mesh, petrel_name(inf.faults[f].name), c.decimals));
      out.push_back({g(p), "STL (ASCII)", "-", inf.faults[f].name + " plane"});
    }
    return out;
  }
};

class CsvExporter final : public Exporter {
 public:
  std::string id() const override { return "csv"; }
  bool petrel() const override { return false; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    std::vector<ExportedFile> out;
    const auto& inf = c.info();
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      const auto& fsu = c.faulted(static_cast<int>(h));
      CsvTable t{{"x", "y", "elevation_m", "depth_m", "block"}, {}};
      std::string xyz;
      for (std::size_t j = 0; j < c.grid().ny; ++j)
        for (std::size_t i = 0; i < c.grid().nx; ++i) {
          std::size_t k = c.grid().index(i, j);
          double z = fsu.z.v[k];
          if (std::isnan(z)) continue;
          t.rows.push_back({num(c.grid().x(i), 2), num(c.grid().y(j), 2), num(z, c.decimals), num(-z, c.decimals), std::to_string(fsu.block[k])});
          xyz += num(c.grid().x(i), 2) + " " + num(c.grid().y(j), 2) + " " + num(z, c.decimals) + "\n";
        }
      std::string hn = petrel_name(inf.horizons[h].name);
      put(folder, fs::path("csv") / (hn + ".csv"), write_csv(t, c.comments()));
      put(folder, fs::path("xyz") / (hn + ".xyz"), xyz);
      out.push_back({"csv/" + hn + ".csv", "CSV", "-", inf.horizons[h].name + " defined nodes: x, y, elevation, depth (positive down), block"});
      out.push_back({"xyz/" + hn + ".xyz", "XYZ points", "General lines/points (ASCII)", inf.horizons[h].name + " defined nodes, z elevation"});
    }
    if (!c.wells().empty()) {
      CsvTable t{{"well", "horizon", "x", "y", "elevation_m", "depth_m", "md_m", "block"}, {}};
      for (std::size_t w = 0; w < c.wells().size(); ++w)
        for (auto& tp : c.tops(w))
          t.rows.push_back({petrel_name(c.wells()[w].name()), petrel_name(tp.horizon), num(tp.p.x, 2), num(tp.p.y, 2), num(tp.p.z, 2), num(-tp.p.z, 2),
                            num(tp.md, 2), std::to_string(tp.block)});
      put(folder, "csv/well_tops.csv", write_csv(t, c.comments()));
      out.push_back({"csv/well_tops.csv", "CSV", "-", "model tops on all wells"});
    }
    // Stage 2F: fault-surface attributes at the lattice nodes inside the tip line (STX-TN-10).
    if (c.fault_rock || c.stress.present)
      for (std::size_t f = 0; f < inf.faults.size(); ++f) {
        const auto& fsu = c.fault_surface(static_cast<int>(f));
        auto pr = repr::fault_rock_surface(c.model(), fsu, c.rock, c.stress);
        CsvTable t{{"x", "y", "elevation_m", "footwall_unit", "hanging_wall_unit", "d_fw_m", "d_hw_m", "sgr", "csp", "thickness_m",
                    "seal_capacity_bar", "sigma_n_MPa", "tau_MPa", "slip_tendency", "dilation_tendency", "reactivation_dp_MPa", "missing"},
                   {}};
        auto nm = [&](double v, int dec) { return std::isfinite(v) ? num(v, dec) : std::string(); };
        auto un = [&](int u) { return u >= 0 && static_cast<std::size_t>(u) < inf.units.size() ? petrel_name(inf.units[static_cast<std::size_t>(u)].name) : std::string(); };
        for (std::size_t k = 0; k < pr.size(); ++k) {
          const auto& r = pr[k];
          if (!r.valid) continue;
          const auto& q = fsu.nodes[k];
          t.rows.push_back({num(q.x, 2), num(q.y, 2), num(q.z, 2), un(r.fw_unit), un(r.hw_unit), nm(r.d_fw, 2), nm(r.d_hw, 2), nm(r.sgr, 4),
                            nm(r.csp, 3), nm(r.thickness, 4), nm(r.seal_bar, 3), nm(r.sigma_n * 1e-6, 4), nm(r.tau * 1e-6, 4), nm(r.ts, 4),
                            nm(r.td, 4), nm(r.dp * 1e-6, 4), r.missing});
        }
        const std::string fn = petrel_name(inf.faults[f].name);
        put(folder, fs::path("csv") / (fn + "_fault_attributes.csv"), write_csv(t, c.comments()));
        out.push_back({"csv/" + fn + "_fault_attributes.csv", "CSV", "-",
                       inf.faults[f].name + " surface nodes: juxtaposed units, SGR, CSP, thickness, seal capacity, stresses, tendencies"});
      }
    return out;
  }
};

class GeoJsonExporter final : public Exporter {
 public:
  std::string id() const override { return "geojson"; }
  bool petrel() const override { return false; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    const auto& inf = c.info();
    const auto& d = inf.domain;
    auto ring = [&](const Polyline& p) {
      std::string s = "[";
      for (std::size_t i = 0; i < p.size(); ++i) s += (i ? "," : "") + std::string("[") + num(p[i].x, 2) + "," + num(p[i].y, 2) + "]";
      return s + "]";
    };
    std::vector<std::string> feats;
    Polyline dom{{d.x0, d.y0, 0}, {d.x0 + d.sx, d.y0, 0}, {d.x0 + d.sx, d.y0 + d.sy, 0}, {d.x0, d.y0 + d.sy, 0}, {d.x0, d.y0, 0}};
    feats.push_back("{\"type\":\"Feature\",\"properties\":{\"kind\":\"domain\",\"name\":\"" + petrel_name(inf.name) +
                    "\"},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[" + ring(dom) + "]}}");
    for (std::size_t h = 0; h < inf.horizons.size(); ++h)
      for (std::size_t f = 0; f < inf.faults.size(); ++f)
        for (auto& poly : c.fault_polygons(static_cast<int>(f), static_cast<int>(h)))
          feats.push_back("{\"type\":\"Feature\",\"properties\":{\"kind\":\"fault_polygon\",\"horizon\":\"" + petrel_name(inf.horizons[h].name) +
                          "\",\"fault\":\"" + petrel_name(inf.faults[f].name) + "\"},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[" + ring(poly) + "]}}");
    for (auto& w : c.wells())
      feats.push_back("{\"type\":\"Feature\",\"properties\":{\"kind\":\"well\",\"name\":\"" + petrel_name(w.name()) +
                      "\"},\"geometry\":{\"type\":\"Point\",\"coordinates\":[" + num(w.x(), 2) + "," + num(w.y(), 2) + "]}}");
    std::string o = "{\"type\":\"FeatureCollection\",\n\"stratonix\":{\"model\":\"" + petrel_name(inf.name) + "\",\"version\":\"" + inf.version +
                    "\",\"evidence_basis\":\"" + basis_label(inf.basis) + "\",\"coordinate_system\":\"" + inf.crs.label() +
                    "\",\"note\":\"coordinates are in the model coordinate system, not WGS 84 as RFC 7946 expects; reproject before use in web maps\"},\n"
                    "\"features\":[\n";
    for (std::size_t i = 0; i < feats.size(); ++i) o += feats[i] + (i + 1 < feats.size() ? ",\n" : "\n");
    o += "]}\n";
    put(folder, "geojson/outlines.geojson", o);
    return {{"geojson/outlines.geojson", "GeoJSON", "-", "domain outline, fault polygons per horizon, well locations (model coordinates)"}};
  }
};

/// Sample points and orientations of horizons from the faulted grids, away from fault gaps.
struct Samples {
  CsvTable points{{"X", "Y", "Z", "surface", "stx_block"}, {}};
  CsvTable orient{{"X", "Y", "Z", "G_x", "G_y", "G_z", "azimuth", "dip", "polarity", "surface"}, {}};
};
Samples sample_horizons(ExportContext& c, int d) {
  Samples s;
  const auto& inf = c.info();
  const auto& gs = c.grid();
  std::size_t step = std::max<std::size_t>(1, static_cast<std::size_t>(std::lround(c.sample_spacing / gs.dx)));
  for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
    const auto& fsu = c.faulted(static_cast<int>(h));
    auto ok = [&](std::size_t i, std::size_t j) {  // node and its neighbours defined and in one block
      if (i == 0 || j == 0 || i + 1 >= gs.nx || j + 1 >= gs.ny) return false;
      std::size_t k = gs.index(i, j);
      for (std::size_t kk : {k, gs.index(i - 1, j), gs.index(i + 1, j), gs.index(i, j - 1), gs.index(i, j + 1)})
        if (std::isnan(fsu.z.v[kk]) || fsu.block[kk] != fsu.block[k]) return false;
      return true;
    };
    // Where a series relation removes the horizon (truncation, onlap: stage 2C) there is no pick.
    auto* sq = static_cast<const StratigraphyQuery*>(c.model().capability("stratigraphy"));
    auto present = [&](std::size_t i, std::size_t j) {
      if (!sq) return true;
      std::size_t k = gs.index(i, j);
      Vec3 p{gs.x(i), gs.y(j), fsu.z.v[k]}, r;
      c.model().restore(std::span<const Vec3>(&p, 1), std::span<Vec3>(&r, 1));
      return !sq->horizon_absent(static_cast<int>(h), r.x, r.y);
    };
    // Where the unit above has zero thickness the pick would coincide with a shallower horizon's
    // pick at another value: a contradiction for an interpolator. Only the shallower is written.
    auto coincides = [&](std::size_t i, std::size_t j) {
      std::size_t k = gs.index(i, j);
      for (std::size_t a = 0; a < h; ++a) {
        double za = c.faulted(static_cast<int>(a)).z.v[k];
        if (!std::isnan(za) && std::fabs(za - fsu.z.v[k]) < 1e-6) return true;
      }
      return false;
    };
    std::string hn = petrel_name(inf.horizons[h].name);
    for (std::size_t j = step / 2; j < gs.ny; j += step)
      for (std::size_t i = step / 2; i < gs.nx; i += step) {
        if (!ok(i, j) || !present(i, j) || coincides(i, j)) continue;
        std::size_t k = gs.index(i, j);
        s.points.rows.push_back({num(gs.x(i), 2), num(gs.y(j), 2), num(fsu.z.v[k], d), hn, std::to_string(fsu.block[k])});
        if ((i / step + j / step) % 2 != 0) continue;  // orientations at every other sample
        double gx = -(fsu.z.v[gs.index(i + 1, j)] - fsu.z.v[gs.index(i - 1, j)]) / (2 * gs.dx);
        double gy = -(fsu.z.v[gs.index(i, j + 1)] - fsu.z.v[gs.index(i, j - 1)]) / (2 * gs.dy);
        Vec3 n = normalized(Vec3{gx, gy, 1});
        double dip = rad2deg(std::acos(std::clamp(n.z, -1.0, 1.0)));
        double az = std::fmod(rad2deg(std::atan2(n.x, n.y)) + 360.0, 360.0);
        s.orient.rows.push_back({num(gs.x(i), 2), num(gs.y(j), 2), num(fsu.z.v[k], d), num(n.x, 6), num(n.y, 6), num(n.z, 6), num(az, 3), num(dip, 3),
                                 "1", hn});
      }
  }
  return s;
}

class GemPyExporter final : public Exporter {
 public:
  std::string id() const override { return "gempy"; }
  bool petrel() const override { return false; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    const auto& inf = c.info();
    Samples s = sample_horizons(c, c.decimals);
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      auto sticks = repr::fault_sticks(c.fault_surface(static_cast<int>(f)), std::max(c.sample_spacing, c.stick_spacing()));
      std::string fn = petrel_name(inf.faults[f].name);
      for (auto& st : sticks)
        for (std::size_t i = 0; i < st.size(); i += std::max<std::size_t>(1, st.size() / 4))
          s.points.rows.push_back({num(st[i].x, 2), num(st[i].y, 2), num(st[i].z, c.decimals), fn, "-1"});
      // one orientation at the middle of the middle stick: the plane normal, pointing up
      if (!sticks.empty()) {
        auto& st = sticks[sticks.size() / 2];
        Vec3 p = st[st.size() / 2];
        Vec3 n = normalized(cross(st.back() - st.front(), sticks.back().front() - sticks.front().front()));
        if (n.z < 0) n = n * -1.0;
        double dip = rad2deg(std::acos(std::clamp(n.z, -1.0, 1.0)));
        double az = std::fmod(rad2deg(std::atan2(n.x, n.y)) + 360.0, 360.0);
        s.orient.rows.push_back({num(p.x, 2), num(p.y, 2), num(p.z, c.decimals), num(n.x, 6), num(n.y, 6), num(n.z, 6), num(az, 3), num(dip, 3), "1", fn});
      }
    }
    // No comment lines: GemPy and pandas read these files directly. Provenance is in model_setup.json.
    put(folder, "gempy/surface_points.csv", write_csv(s.points, {}));
    put(folder, "gempy/orientations.csv", write_csv(s.orient, {}));
    const auto& d = inf.domain;
    std::string js = "{\n  \"model\": \"" + petrel_name(inf.name) + "\",\n  \"extent\": [" + num(d.x0, 2) + ", " + num(d.x0 + d.sx, 2) + ", " + num(d.y0, 2) +
                     ", " + num(d.y0 + d.sy, 2) + ", " + num(d.z_base, 2) + ", " + num(d.z_top, 2) + "],\n  \"surfaces_top_to_base\": [";
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) js += (h ? ", " : "") + std::string("\"") + petrel_name(inf.horizons[h].name) + "\"";
    js += "],\n  \"faults_oldest_first\": [";
    for (std::size_t f = 0; f < inf.faults.size(); ++f) js += (f ? ", " : "") + std::string("\"") + petrel_name(inf.faults[f].name) + "\"";
    js += "],\n  \"z\": \"elevation (m, negative below the datum)\",\n  \"sampling_m\": " + num(c.sample_spacing, 1) +
          ",\n  \"provenance\": \"" + c.comments()[0] + "; " + c.comments()[1] + "; " + c.comments()[3] +
          "\",\n  \"note\": \"points sampled from the model away from fault gaps; stx_block = Stratonix fault block (-1 for fault points); column names follow GemPy 2 (X, Y, Z, surface; G_x, G_y, G_z, azimuth, dip, polarity); check against your GemPy version\"\n}\n";
    put(folder, "gempy/model_setup.json", js);
    return {{"gempy/surface_points.csv", "GemPy surface points", "-", std::to_string(s.points.rows.size()) + " points (horizons and faults)"},
            {"gempy/orientations.csv", "GemPy orientations", "-", std::to_string(s.orient.rows.size()) + " orientations"},
            {"gempy/model_setup.json", "JSON", "-", "extent, surface order, fault order"}};
  }
};

class LoopExporter final : public Exporter {
 public:
  std::string id() const override { return "loopstructural"; }
  bool petrel() const override { return false; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    const auto& inf = c.info();
    Samples s = sample_horizons(c, c.decimals);
    // Stratigraphic value: mean elevation of each horizon (constant on a horizon, larger above).
    std::map<std::string, double> val;
    std::map<std::string, int> cnt;
    for (auto& r : s.points.rows) { val[r[3]] += parse_double(r[2]); ++cnt[r[3]]; }
    for (auto& [k, v] : val) v /= cnt[k];
    CsvTable t{{"X", "Y", "Z", "feature_name", "val", "gx", "gy", "gz", "coord", "stx_horizon", "stx_block"}, {}};
    for (auto& r : s.points.rows) t.rows.push_back({r[0], r[1], r[2], "strati", num(val[r[3]], 3), "", "", "", "0", r[3], r[4]});
    for (auto& r : s.orient.rows) t.rows.push_back({r[0], r[1], r[2], "strati", "", r[3], r[4], r[5], "0", r[9], ""});
    std::string faults_json = "[";
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      std::string fn = petrel_name(inf.faults[f].name);
      auto sticks = repr::fault_sticks(c.fault_surface(static_cast<int>(f)), std::max(c.sample_spacing, c.stick_spacing()));
      Vec3 n{0, 0, 1};
      if (sticks.size() > 1) n = normalized(cross(sticks.front().back() - sticks.front().front(), sticks.back().front() - sticks.front().front()));
      for (auto& st : sticks)
        for (std::size_t i = 0; i < st.size(); i += std::max<std::size_t>(1, st.size() / 4))
          t.rows.push_back({num(st[i].x, 2), num(st[i].y, 2), num(st[i].z, c.decimals), fn, "0", num(n.x, 6), num(n.y, 6), num(n.z, 6), "0", "", "-1"});
      FaultEval fe;
      const auto& m = c.fault_surface(static_cast<int>(f)).mesh;
      c.model().evaluate_fault(static_cast<int>(f), m.v, fe);
      double dmax = 0;
      for (double v : fe.displacement) dmax = std::max(dmax, v);
      faults_json += std::string(f ? ", " : "") + "{\"name\": \"" + fn + "\", \"displacement_max_m\": " + num(dmax, 3) + "}";
    }
    faults_json += "]";
    put(folder, "loopstructural/data.csv", write_csv(t, {}));  // no comments: read directly by pandas
    const auto& d = inf.domain;
    std::string js = "{\n  \"origin\": [" + num(d.x0, 2) + ", " + num(d.y0, 2) + ", " + num(d.z_base, 2) + "],\n  \"maximum\": [" + num(d.x0 + d.sx, 2) + ", " +
                     num(d.y0 + d.sy, 2) + ", " + num(d.z_top, 2) + "],\n  \"strati_values\": {";
    bool first = true;
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      std::string hn = petrel_name(inf.horizons[h].name);
      if (!val.count(hn)) continue;
      js += std::string(first ? "" : ", ") + "\"" + hn + "\": " + num(val[hn], 3);
      first = false;
    }
    js += "},\n  \"faults_oldest_first\": " + faults_json + ",\n  \"provenance\": \"" + c.comments()[0] + "; " + c.comments()[1] + "; " +
          c.comments()[3] + "\",\n  \"note\": \"columns X, Y, Z, feature_name, val, gx, gy, gz, coord as in the LoopStructural input-data guide; stx_horizon and stx_block are Stratonix extras for checks\"\n}\n";
    put(folder, "loopstructural/model_setup.json", js);
    return {{"loopstructural/data.csv", "LoopStructural data", "-", std::to_string(t.rows.size()) + " constraints (strati interfaces and normals, fault points)"},
            {"loopstructural/model_setup.json", "JSON", "-", "bounding box, strati value per horizon, faults with maximum displacement"}};
  }
};
}  // namespace

STRATONIX_REGISTER(exp_vtk, registry::Kind::Exporter, "vtk", "1.0.0", "VTK XML (VTU surfaces and wells, VTI voxels), PVD and ParaView state",
                   [] { return std::make_unique<VtkExporter>(); })
STRATONIX_REGISTER(exp_obj, registry::Kind::Exporter, "obj", "1.0.0", "Wavefront OBJ surfaces", [] { return std::make_unique<ObjExporter>(); })
STRATONIX_REGISTER(exp_stl, registry::Kind::Exporter, "stl", "1.0.0", "ASCII STL surfaces", [] { return std::make_unique<StlExporter>(); })
STRATONIX_REGISTER(exp_csv, registry::Kind::Exporter, "csv", "1.0.0", "CSV and XYZ horizon nodes; CSV well tops", [] { return std::make_unique<CsvExporter>(); })
STRATONIX_REGISTER(exp_geojson, registry::Kind::Exporter, "geojson", "1.0.0", "GeoJSON outlines: domain, fault polygons, wells",
                   [] { return std::make_unique<GeoJsonExporter>(); })
STRATONIX_REGISTER(exp_gempy, registry::Kind::Exporter, "gempy", "1.0.0", "GemPy input set sampled from the model", [] { return std::make_unique<GemPyExporter>(); })
STRATONIX_REGISTER(exp_loop, registry::Kind::Exporter, "loopstructural", "1.0.0", "LoopStructural input set sampled from the model",
                   [] { return std::make_unique<LoopExporter>(); })

void register_neutral_exporters() {
  STRATONIX_PLUGIN_ANCHOR(exp_vtk);
  STRATONIX_PLUGIN_ANCHOR(exp_obj);
  STRATONIX_PLUGIN_ANCHOR(exp_stl);
  STRATONIX_PLUGIN_ANCHOR(exp_csv);
  STRATONIX_PLUGIN_ANCHOR(exp_geojson);
  STRATONIX_PLUGIN_ANCHOR(exp_gempy);
  STRATONIX_PLUGIN_ANCHOR(exp_loop);
}

}  // namespace stratonix::io
