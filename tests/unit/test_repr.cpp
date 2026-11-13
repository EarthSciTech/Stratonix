// Representation layer (Prompt 5): extracted surfaces, meshes, voxels, wells and sections agree
// with point evaluation on the synthetic models.
#include <filesystem>
#include <random>

#include "../stx_test.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

using namespace stratonix;

static ModelPtr synthetic(desc::Description* dd = nullptr) {
  register_all_plugins();
  auto d = desc::load_description(std::string(stxtest::source_dir()) + "/examples/synthetic_fold_two_faults.yaml");
  BuildOptions o;
  auto m = build_only(d, o).model;
  if (dd) *dd = d;
  return m;
}

TEST_CASE("P5 faulted surfaces lie on the horizon and in the block they report") {
  auto m = synthetic();
  GridSpec g = GridSpec::covering(m->info().domain, 100);
  for (int h = 0; h < 5; ++h) {
    auto fs = repr::horizon_grid_faulted(*m, h, g);
    std::vector<Vec3> p;
    std::vector<BlockCode> want;
    for (std::size_t k = 0; k < g.size(); ++k)
      if (!std::isnan(fs.z.v[k])) { p.push_back({g.x(k % g.nx), g.y(k / g.nx), fs.z.v[k]}); want.push_back(fs.block[k]); }
    EvalResult r;
    m->evaluate(p, kEvalScalars | kEvalBlock, r);
    double worst = 0;
    std::size_t wrong_block = 0;
    for (std::size_t i = 0; i < p.size(); ++i) {
      worst = std::max(worst, std::fabs(r.s(static_cast<std::size_t>(h), i)));
      wrong_block += r.block[i] != want[i];
    }
    CHECK(worst < 1e-6);
    CHECK(wrong_block == 0);
    if (h == 2) CHECK(fs.gap_nodes > 0);  // H3 is cut by both faults: fault gaps exist
  }
}

TEST_CASE("P5 per-block surfaces have zero scalar in their block") {
  auto m = synthetic();
  auto blocks = repr::list_blocks(*m);
  CHECK(blocks.codes.size() >= 3);
  GridSpec g = GridSpec::covering(m->info().domain, 200);
  for (auto b : blocks.codes) {
    auto z = repr::horizon_grid_block(*m, 2, b, g);
    std::vector<Vec3> p;
    for (std::size_t k = 0; k < g.size(); ++k)
      if (!std::isnan(z.v[k])) p.push_back({g.x(k % g.nx), g.y(k / g.nx), z.v[k]});
    std::vector<double> s(p.size());
    m->scalar_in_block(p, b, 2, s);
    double worst = 0;
    for (double v : s) worst = std::max(worst, std::fabs(v));
    CHECK(worst < 1e-6);
  }
}

TEST_CASE("P5 voxel grid agrees with per-block horizon surfaces") {
  auto m = synthetic();
  auto v = repr::voxel_grid(*m, 200, 200, 50, true);
  CHECK(v.nx == 40);
  CHECK(v.ny == 30);
  CHECK(v.nz == 50);
  std::mt19937 rng(7);
  std::size_t mismatch = 0, n = 0, near = 0;
  for (int s = 0; s < 3000; ++s) {
    std::size_t i = rng() % v.nx, j = rng() % v.ny, k = rng() % v.nz;
    std::size_t c = i + v.nx * (j + v.ny * k);
    Vec3 p = v.centre(i, j, k);
    int unit = 0;
    double dmin = 1e9;
    for (int h = 0; h < 5; ++h) {
      auto z = repr::horizon_grid_block(*m, h, v.block[c], GridSpec{p.x, p.y, 1, 1, 1, 1});
      if (p.z < z.v[0]) unit = h + 1;
      dmin = std::min(dmin, std::fabs(p.z - z.v[0]));
    }
    if (dmin < 1e-6) { ++near; continue; }
    ++n;
    mismatch += unit != v.unit[c];
  }
  CHECK(n > 2900);
  CHECK(n + near == 3000);
  CHECK(mismatch == 0);
}

TEST_CASE("P5 horizon meshes and fault meshes are well formed") {
  auto m = synthetic();
  GridSpec g = GridSpec::covering(m->info().domain, 100);
  auto fs = repr::horizon_grid_faulted(*m, 2, g);
  auto mesh = repr::grid_to_mesh(fs.z);
  auto q = mesh_quality(mesh, 1e-3);
  CHECK(q.degenerate == 0);
  CHECK(q.nonmanifold_edges == 0);
  CHECK(!q.watertight());  // open surfaces by design (gaps and domain edges)
  for (int f = 0; f < 2; ++f) {
    auto fm = repr::fault_surface(*m, f, 50);
    auto fq = mesh_quality(fm.mesh, 1e-6);
    CHECK(fq.degenerate == 0);
    CHECK(fq.nonmanifold_edges == 0);
    CHECK(fq.triangles > 100);
    // every mesh vertex lies on the fault surface
    FaultEval e;
    m->evaluate_fault(f, fm.mesh.v, e);
    double worst = 0;
    for (double dd : e.distance) worst = std::max(worst, std::fabs(dd));
    CHECK(worst < 1e-3);
    auto sticks = repr::fault_sticks(fm, 250);
    CHECK(sticks.size() >= 5);
  }
}

TEST_CASE("P5 well intersections along vertical and deviated wells") {
  desc::Description d;
  auto m = synthetic(&d);
  auto w1 = make_well_path(d.wells[0], d);
  auto tops = repr::well_intersections(*m, w1);
  REQUIRE(tops.size() == 5);
  CHECK_NEAR(tops[0].p.z, -800, 1e-6);
  CHECK_NEAR(tops[4].p.z, -1600, 1e-6);
  // A deviated well crossing F2: each top lies on its horizon (scalar zero) in its own block.
  std::vector<WellStation> st{{0, 0, 0}, {800, 0, deg2rad(120)}, {1400, deg2rad(40), deg2rad(120)}, {2400, deg2rad(55), deg2rad(120)}};
  WellPath dev("DEV", 4300, 3600, 0, st);
  auto t2 = repr::well_intersections(*m, dev, 2.0);
  CHECK(t2.size() >= 4);
  std::vector<Vec3> p;
  for (auto& t : t2) p.push_back(t.p);
  EvalResult r;
  m->evaluate(p, kEvalScalars | kEvalBlock, r);
  for (std::size_t i = 0; i < t2.size(); ++i) {
    CHECK(std::fabs(r.s(static_cast<std::size_t>(t2[i].horizon_index), i)) < 1e-4);
    CHECK(r.block[i] == t2[i].block);
    CHECK(t2[i].md > 0);
  }
}

TEST_CASE("P5 sections are slices of the 3D model") {
  desc::Description d;
  auto m = synthetic(&d);
  Polyline line{{1500, 1500, 0}, {6500, 4500, 0}};
  auto s = repr::make_section(*m, "XS", line, 25, 5);
  CHECK(s.ns > 200);
  // horizon line points agree with point evaluation
  std::vector<Vec3> p;
  for (std::size_t c = 0; c < s.ns; ++c)
    if (!std::isnan(s.horizon_z[2][c])) p.push_back({s.xy[c].x, s.xy[c].y, s.horizon_z[2][c]});
  EvalResult r;
  m->evaluate(p, kEvalScalars, r);
  double worst = 0;
  for (std::size_t i = 0; i < p.size(); ++i) worst = std::max(worst, std::fabs(r.s(2, i)));
  CHECK(worst < 1e-6);
  // F2 trace present on the section
  std::size_t nf = 0;
  for (double z : s.fault_z[1]) nf += !std::isnan(z);
  CHECK(nf > 0);
}

TEST_CASE("P5 contour lines and closure mask") {
  auto m = synthetic();
  GridSpec g = GridSpec::covering(m->info().domain, 50);
  auto fs = repr::horizon_grid_faulted(*m, 1, g);
  auto lines = repr::contour_lines(fs.z, -950);
  CHECK(!lines.empty());
  for (auto& l : lines)
    for (auto& p : l) {
      double v = fs.z.sample(p.x, p.y);
      if (!std::isnan(v)) CHECK_NEAR(v, -950, 2.0);  // points on the edge of a fault gap have no bilinear value
    }
}

STX_TEST_MAIN()
