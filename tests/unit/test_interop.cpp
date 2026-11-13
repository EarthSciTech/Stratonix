// Neutral exporters and importers: round trips against the model.
#include <cmath>
#include <filesystem>

#include "../stx_test.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/desc/document.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/io/export.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/registry/registry.hpp"
#include "stratonix/repr/repr.hpp"

using namespace stratonix;
namespace fs = std::filesystem;

namespace {
struct Fixture {
  desc::Description d;
  ModelPtr m;
  std::unique_ptr<io::ExportContext> c;
  fs::path out;
};
Fixture& fx() {
  static Fixture f = [] {
    register_all_plugins();
    Fixture x;
    x.d = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples/synthetic_fold_two_faults.yaml");
    BuildLog log;
    x.m = build_model(x.d, {}, log);
    x.c = std::make_unique<io::ExportContext>(*x.m, make_well_paths(x.d), GridSpec::covering(x.m->info().domain, 50), 25, 250);
    x.c->voxel_dx = 200; x.c->voxel_dy = 200; x.c->voxel_dz = 50;
    x.c->vertical_exaggeration = 3;
    x.out = fs::path(STX_TMP_DIR) / "export_interop";
    fs::remove_all(x.out);
    auto ids = io::petrel_set_ids();
    for (const char* e : {"vtk", "obj", "stl", "csv", "geojson", "gempy", "loopstructural"}) ids.push_back(e);
    io::export_set(*x.c, ids, x.out, "test");
    return x;
  }();
  return f;
}
double max_dist(const TriMesh& a, const TriMesh& b) {
  double w = 0;
  for (std::size_t i = 0; i < a.v.size(); ++i) w = std::max(w, norm(a.v[i] - b.v[i]));
  return w;
}
io::ImportedData import(const std::string& id, const fs::path& f) {
  return registry::Registry::global().create<io::Importer>(registry::Kind::Importer, id)->read(f);
}
}  // namespace

TEST_CASE("I1 VTU surfaces, wells and VTI voxels read back equal to the model") {
  auto& f = fx();
  auto h3 = io::read_vtu(read_file(f.out / "vtk/H3.vtu"));
  auto ref = repr::grid_to_mesh(f.c->faulted(2).z);
  REQUIRE(h3.pts.size() == ref.v.size());
  REQUIRE(h3.cells.size() == ref.t.size());
  double w = 0;
  for (std::size_t i = 0; i < ref.v.size(); ++i) w = std::max(w, norm(h3.pts[i] - ref.v[i]));
  CHECK(w < 1e-4);
  REQUIRE(h3.point_data.size() == 1);
  CHECK(h3.point_data[0].first == "depth_m");
  CHECK_NEAR(h3.point_data[0].second[10], -ref.v[10].z, 1e-4);
  CHECK(h3.types[0] == 5);
  CHECK(h3.cells[7][2] == static_cast<std::int64_t>(ref.t[7][2]));
  auto wells = io::read_vtu(read_file(f.out / "vtk/wells.vtu"));
  CHECK(wells.cells.size() == f.c->wells().size());
  CHECK(wells.types[0] == 4);
  auto v = io::read_vti(read_file(f.out / "vtk/voxels.vti"));
  auto vr = repr::voxel_grid(*f.m, 200, 200, 50, true);
  CHECK(v.nx == vr.nx && v.ny == vr.ny && v.nz == vr.nz);
  CHECK(v.unit == vr.unit);
  CHECK(v.block == vr.block);
  CHECK_NEAR(v.z0, vr.z0, 1e-4);
  auto pvd = read_file(f.out / "vtk/model.pvd");
  CHECK(pvd.find("file=\"voxels.vti\"") != std::string::npos);
  CHECK(read_file(f.out / "vtk/paraview_state.py").find("VE = 3") != std::string::npos);
}

TEST_CASE("I2 OBJ and STL read back equal to the meshes") {
  auto& f = fx();
  const auto& ref = f.c->fault_surface(0).mesh;
  auto o = io::read_obj(read_file(f.out / "obj/F1.obj"));
  CHECK(o.t.size() == ref.t.size());
  CHECK(o.v.size() == ref.v.size());
  CHECK(max_dist(o, ref) < 1e-4);
  CHECK(o.t[5] == ref.t[5]);
  auto s = io::read_stl(read_file(f.out / "stl/F1.stl"));
  CHECK(s.t.size() == ref.t.size());
  double w = 0;
  for (std::size_t t = 0; t < ref.t.size(); ++t)
    for (int k = 0; k < 3; ++k) w = std::max(w, norm(s.v[s.t[t][static_cast<std::size_t>(k)]] - ref.v[ref.t[t][static_cast<std::size_t>(k)]]));
  CHECK(w < 1e-4);
}

TEST_CASE("I3 CSV, XYZ and GeoJSON read back") {
  auto& f = fx();
  const auto& fsu = f.c->faulted(1);
  std::size_t defined = 0;
  for (double z : fsu.z.v) defined += !std::isnan(z);
  auto t = io::read_csv(read_file(f.out / "csv/H2.csv"));
  CHECK(t.header.size() == 5 && t.header[3] == "depth_m");
  CHECK(t.rows.size() == defined);
  CHECK_NEAR(parse_double(t.rows[0][2]), -parse_double(t.rows[0][3]), 1e-9);
  auto x = io::read_xyz_points(read_file(f.out / "xyz/H2.xyz"));
  CHECK(x.size() == defined);
  CHECK_NEAR(x[3].z, parse_double(t.rows[3][2]), 1e-9);
  auto gj = desc::parse_document(read_file(f.out / "geojson/outlines.geojson"));
  REQUIRE(gj.is_map());
  const auto* feats = gj.get("features");
  REQUIRE(feats != nullptr);
  std::size_t polys = 0;
  for (int h = 0; h < 5; ++h)
    for (int ff = 0; ff < 2; ++ff) polys += f.c->fault_polygons(ff, h).size();
  CHECK(feats->seq.size() == 1 + polys + f.c->wells().size());
}

TEST_CASE("I4 GemPy and LoopStructural points lie on the model's horizons and faults") {
  auto& f = fx();
  const auto& inf = f.m->info();
  auto check_points = [&](const io::CsvTable& t, int name_col) {
    std::vector<Vec3> hp, fp;
    std::vector<int> hi, fi;
    for (auto& r : t.rows) {
      Vec3 p{parse_double(r[0]), parse_double(r[1]), parse_double(r[2])};
      const std::string& n = r[static_cast<std::size_t>(name_col)];
      for (std::size_t h = 0; h < inf.horizons.size(); ++h)
        if (inf.horizons[h].name == n) { hp.push_back(p); hi.push_back(static_cast<int>(h)); }
      for (std::size_t ff = 0; ff < inf.faults.size(); ++ff)
        if (inf.faults[ff].name == n) { fp.push_back(p); fi.push_back(static_cast<int>(ff)); }
    }
    CHECK(hp.size() > 100);
    CHECK(!fp.empty());
    EvalResult er;
    f.m->evaluate(hp, kEvalScalars, er);
    double w = 0;
    for (std::size_t i = 0; i < hp.size(); ++i) w = std::max(w, std::fabs(er.s(static_cast<std::size_t>(hi[i]), i)));
    CHECK(w < 1e-3);
    double wf = 0;
    for (std::size_t i = 0; i < fp.size(); ++i) {
      FaultEval fe;
      std::vector<Vec3> one{fp[i]};
      f.m->evaluate_fault(fi[i], one, fe);
      wf = std::max(wf, std::fabs(fe.distance[0]));
    }
    CHECK(wf < 0.01);
  };
  auto g = io::read_csv(read_file(f.out / "gempy/surface_points.csv"));
  CHECK(g.header == std::vector<std::string>({"X", "Y", "Z", "surface", "stx_block"}));
  check_points(g, 3);
  auto o = io::read_csv(read_file(f.out / "gempy/orientations.csv"));
  CHECK(o.rows.size() > 20);
  for (auto& r : o.rows) CHECK_NEAR(std::hypot(parse_double(r[3]), parse_double(r[4]), parse_double(r[5])), 1.0, 1e-5);
  auto l = io::read_csv(read_file(f.out / "loopstructural/data.csv"));
  CHECK(l.header[3] == "feature_name" && l.header[4] == "val");
  // LoopStructural rows: horizon name in stx_horizon (col 9) for strati rows, feature_name for faults
  io::CsvTable lp{l.header, {}};
  for (auto& r : l.rows) {
    if (r[4].empty()) continue;  // normals
    auto rr = r;
    if (rr[3] == "strati") rr[3] = rr[9];
    lp.rows.push_back(rr);
  }
  check_points(lp, 3);
  auto js = desc::parse_document(read_file(f.out / "loopstructural/model_setup.json"));
  CHECK(js.get("strati_values") != nullptr);
}

TEST_CASE("I5 importers read every exported format back") {
  auto& f = fx();
  auto pts = import("xyz_points", f.out / "xyz/H3.xyz");
  CHECK(pts.kind == "points");
  CHECK(!pts.points.empty());
  auto ref = repr::fault_sticks(f.c->fault_surface(1), 250);
  auto st = import("irap_lines", f.out / "faults/F2_sticks.dat");
  REQUIRE(st.lines.size() == ref.size());
  CHECK(norm(st.lines[3].pts[1] - ref[3][1]) < 1e-4);
  auto pp = import("petrel_points", f.out / "faults/fault_sticks_points.txt");
  CHECK(pp.kind == "lines");
  std::size_t n2 = 0;
  for (auto& l : pp.lines)
    if (l.name == "F2") {
      ++n2;
      CHECK(l.pts.size() == ref[static_cast<std::size_t>(l.id - 1)].size());
    }
  CHECK(n2 == ref.size());
  auto ts = import("gocad_tsurf", f.out / "faults/F1.ts");
  REQUIRE(ts.meshes.size() == 1);
  CHECK(ts.meshes[0].name == "F1");
  CHECK(ts.z_is == "elevation");
  {
    const auto& a = ts.meshes[0].mesh;
    const auto& b = f.c->fault_surface(0).mesh;  // TSurf keeps only vertices used by triangles
    REQUIRE(a.t.size() == b.t.size());
    double w = 0;
    for (std::size_t t = 0; t < a.t.size(); ++t)
      for (std::size_t k = 0; k < 3; ++k) w = std::max(w, norm(a.v[a.t[t][k]] - b.v[b.t[t][k]]));
    CHECK(w < 1e-4);
  }
  io::TSurfInfo ti;
  io::read_tsurf(read_file(f.out / "faults/F1.ts"), &ti);
  CHECK(ti.property == "displacement");
  CHECK(ti.values.size() == ts.meshes[0].mesh.v.size());
  CHECK(ti.crs_name == "local_engineering_frame");
  auto tp = import("well_tops", f.out / "wells/well_tops.txt");
  CHECK(tp.tops.size() == 15);
  CHECK(tp.tops[0].well == "W1");
  auto gi = import("irap_ascii", f.out / "horizons/H4_faulted.irap");
  REQUIRE(gi.grid.has_value());
  CHECK(gi.grid->spec.nx == 161);
  auto gz = import("zmap", f.out / "horizons/H4_faulted.zmap");
  REQUIRE(gz.grid.has_value());
  std::size_t diff = 0;
  for (std::size_t k = 0; k < gi.grid->v.size(); ++k)
    if (std::isnan(gi.grid->v[k]) != std::isnan(gz.grid->v[k]) || (!std::isnan(gi.grid->v[k]) && std::fabs(gi.grid->v[k] - gz.grid->v[k]) > 1e-9)) ++diff;
  CHECK(diff == 0);
}

STX_TEST_MAIN()
