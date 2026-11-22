// Stage 4D tests (STX-TN-17 section 9): cases S1 to S13 of properties, uncertainty and scale.
// Written first (Prompt 58): they fail with "not implemented" until Prompts 59 to 61 build the
// features. Every random field is seeded; analytical values are stated in each case.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../stx_test.hpp"
#include "stratonix/checks/checks.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/geostat.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/package/uncertainty.hpp"
#include "stratonix/repr/grid.hpp"

using namespace stratonix;
namespace fs = std::filesystem;

namespace {
const double kPiD = 3.14159265358979323846;
fs::path tmp(const std::string& name) {
  fs::path p = fs::path(STX_TMP_DIR) / "v4_scale" / name;
  fs::remove_all(p);
  fs::create_directories(p.parent_path());
  return p;
}
fs::path example(const std::string& name) { return fs::path(STX_SOURCE_DIR) / "examples" / (name + ".yaml"); }
geo::Variogram sph(double sill, double range, double nugget = 0) {
  geo::Variogram v;
  v.nugget = nugget;
  v.structures.push_back({"spherical", sill - nugget, range, range, range, 0});
  return v;
}
std::vector<Vec3> lattice(std::size_t nx, std::size_t ny, std::size_t nz, double h, double hz) {
  std::vector<Vec3> p;
  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i < nx; ++i) p.push_back({i * h, j * h, k * hz});
  return p;
}
std::string file_bytes(const fs::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}
}  // namespace

// ---- 1. Geostatistics (Prompt 59) ----------------------------------------------------------------

TEST_CASE("S1 variogram models: values, sill and practical ranges") {
  auto v = sph(1.0, 100.0, 0.2);
  CHECK(std::fabs(v.sill() - 1.0) < 1e-12);
  CHECK(std::fabs(v.gamma({0, 0, 0})) < 1e-12);                 // no nugget at h = 0 exactly
  CHECK(std::fabs(v.gamma({50, 0, 0}) - (0.2 + 0.8 * (1.5 * 0.5 - 0.5 * 0.125))) < 1e-12);
  CHECK(std::fabs(v.gamma({150, 0, 0}) - 1.0) < 1e-12);
  geo::Variogram e;
  e.structures.push_back({"exponential", 1.0, 300, 300, 30, 0});
  CHECK(std::fabs(e.gamma({300, 0, 0}) - 0.95) < 1e-3);           // practical range: 95 % of the sill
  CHECK(std::fabs(e.gamma({0, 0, 30}) - 0.95) < 1e-3);            // vertical range
  geo::Variogram g;
  g.structures.push_back({"gaussian", 2.0, 200, 100, 200, 90});   // major range east-west
  CHECK(std::fabs(g.gamma({200, 0, 0}) / 2.0 - 0.95) < 1e-3);
  CHECK(std::fabs(g.gamma({0, 100, 0}) / 2.0 - 0.95) < 1e-3);
  CHECK(std::fabs(g.covariance({0, 0, 0}) - 2.0) < 1e-12);
}

TEST_CASE("S2 experimental variogram of a known series") {
  // z = 0, 1, 0, 1, ... at unit spacing: gamma(1) = 0.5, gamma(2) = 0
  std::vector<Vec3> p;
  std::vector<double> z;
  for (int i = 0; i < 11; ++i) { p.push_back({double(i), 0, 0}); z.push_back(i % 2); }
  geo::VariogramDirection d;
  d.lag = 1;
  d.lag_tolerance = 0.25;
  d.lags = 3;
  auto e = geo::experimental_variogram(p, z, d);
  REQUIRE(e.gamma.size() == 3);
  CHECK(e.pairs[0] == 10 && e.pairs[1] == 9 && e.pairs[2] == 8);
  CHECK(std::fabs(e.gamma[0] - 0.5) < 1e-12);
  CHECK(std::fabs(e.gamma[1]) < 1e-12);
  CHECK(std::fabs(e.distance[2] - 3.0) < 1e-12);
}

TEST_CASE("S3 normal-score transform and declustering") {
  std::vector<double> z;
  for (int i = 1; i <= 99; ++i) z.push_back(std::exp(0.03 * i));
  geo::NormalScore ns(z, std::vector<double>(z.size(), 1.0), 0.5, 30.0);
  for (double x : {1.1, 2.0, 5.0, 15.0}) CHECK(std::fabs(ns.back(ns.forward(x)) - x) / x < 1e-6);
  CHECK(std::fabs(ns.forward(std::exp(0.03 * 50))) < 0.02);        // the median maps to 0
  // a cluster of 9 points and one isolated point: the isolated one weighs as much as the cluster
  std::vector<Vec3> p;
  for (int i = 0; i < 9; ++i) p.push_back({1.0 + 0.01 * i, 1.0, 0});
  p.push_back({50, 50, 0});
  auto w = geo::decluster_weights(p, 10.0);
  REQUIRE(w.size() == 10);
  double wc = 0;
  for (int i = 0; i < 9; ++i) wc += w[i];
  CHECK(std::fabs(wc - w[9]) / w[9] < 1e-9);
}

TEST_CASE("S4 sequential Gaussian simulation: conditioning reproduced, model variogram reproduced") {
  const auto nodes = lattice(60, 60, 1, 10, 1);
  const auto v = sph(1.0, 150.0);
  geo::SearchOptions s;
  s.radius_major = s.radius_minor = 300;
  s.radius_vertical = 10;
  s.max_data = 16;
  std::vector<geo::ConditioningPoint> data = {{0, 1.5}, {1830, -0.7}, {3599, 0.3}};
  // average experimental variogram along x over 20 realisations
  geo::VariogramDirection d;
  d.azimuth = 90;
  d.angle_tolerance = 10;
  d.lag = 10;
  d.lag_tolerance = 5;
  d.lags = 15;
  std::vector<double> mean_g(d.lags, 0.0);
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    auto y = geo::sgs(nodes, data, v, s, seed);
    REQUIRE(y.size() == nodes.size());
    for (auto& c : data) CHECK(y[c.node] == c.value);
    auto e = geo::experimental_variogram(nodes, y, d);
    for (std::size_t l = 0; l < d.lags; ++l) mean_g[l] += e.gamma[l] / 20.0;
  }
  for (std::size_t l = 0; l < d.lags; ++l) {
    const double h = (l + 1) * 10.0;
    std::printf("  lag %4.0f m: gamma %.3f, model %.3f\n", h, mean_g[l], v.gamma({h, 0, 0}));
    CHECK(std::fabs(mean_g[l] - v.gamma({h, 0, 0})) < 0.1 * v.sill());  // STX-TN-17 section 1.7
  }
}

TEST_CASE("S5 sequential indicator simulation reproduces the proportions") {
  const auto nodes = lattice(50, 50, 4, 20, 2);
  const std::vector<double> prop = {0.6, 0.3, 0.1};
  std::vector<geo::Variogram> v = {sph(0.24, 200), sph(0.21, 150), sph(0.09, 100)};
  geo::SearchOptions s;
  s.radius_major = s.radius_minor = 400;
  s.radius_vertical = 8;
  std::vector<double> got(3, 0.0);
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    auto f = geo::sis(nodes, {{0, 2}, {5000, 1}}, prop, v, s, seed);
    REQUIRE(f.size() == nodes.size());
    CHECK(f[0] == 2 && f[5000] == 1);
    for (int c : f) got[static_cast<std::size_t>(c)] += 1.0 / (20.0 * f.size());
  }
  for (std::size_t c = 0; c < 3; ++c) {
    std::printf("  facies %zu: proportion %.3f, target %.3f\n", c, got[c], prop[c]);
    CHECK(std::fabs(got[c] - prop[c]) < 0.03);
  }
}

TEST_CASE("S6 object-based facies: target proportion and wells honoured") {
  const auto nodes = lattice(80, 80, 10, 25, 2);
  geo::ObjectSpec o;
  o.kind = "channel";
  o.width = 150;
  o.thickness = 6;
  o.sinuosity = 1.3;
  o.wavelength = 800;
  o.azimuth = 30;
  o.proportion = 0.25;
  std::vector<std::pair<std::size_t, int>> wells;
  for (std::size_t k = 0; k < 10; ++k) wells.push_back({k * 6400 + 40 * 80 + 40, k < 3 ? 1 : 0});
  auto r = geo::object_facies(nodes, {0, 0, 0}, {2000, 2000, 20}, o, wells, 7, 10000);
  CHECK(std::fabs(r.proportion - 0.25) < 0.05);
  CHECK(r.well_cells == 10);
  CHECK(r.well_misfit == 0);
  CHECK(r.objects > 0);
}

TEST_CASE("S7 realisations: same seed, same bytes; another seed, another field") {
  const auto nodes = lattice(30, 30, 3, 10, 2);
  const auto v = sph(1.0, 80.0);
  geo::SearchOptions s;
  s.radius_major = s.radius_minor = 160;
  s.radius_vertical = 6;
  auto a = geo::sgs(nodes, {}, v, s, 42), b = geo::sgs(nodes, {}, v, s, 42), c = geo::sgs(nodes, {}, v, s, 43);
  CHECK(std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
  CHECK(std::memcmp(a.data(), c.data(), a.size() * sizeof(double)) != 0);
  double sum = 0;
  for (double x : a) sum += x;
  std::printf("  checksum of seed 42: %.17g\n", sum);
  // the same bytes on Windows and Linux (STX-TN-17 section 1.3): set from the Linux build
  CHECK(sum == -404.38799494961626);
}

TEST_CASE("S8 simulation in geological time: correlation follows a folded layer") {
  // nodes on a folded layer: in present-day space neighbours along the layer are far apart in z;
  // in (u, v, w) they are adjacent. The variogram along the layer must reproduce the model.
  std::vector<Vec3> uvw, xyz;
  for (int i = 0; i < 100; ++i)
    for (int k = 0; k < 5; ++k) {
      uvw.push_back({i * 10.0, 0, k * 2.0});
      xyz.push_back({i * 10.0, 0, 300 * std::sin(2 * kPiD * i / 50.0) + k * 2.0});
    }
  const auto v = sph(1.0, 100.0);
  geo::SearchOptions s;
  s.radius_major = s.radius_minor = 200;
  s.radius_vertical = 10;
  geo::VariogramDirection d;
  d.azimuth = 90;
  d.lag = 10;
  d.lag_tolerance = 5;
  d.lags = 8;
  // the variogram along the layer, averaged over 20 realisations, reproduces the model; measured
  // along x in present-day space the same values decorrelate across the fold's limbs
  double g_layer = 0, g_space = 0;
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    auto y = geo::sgs(uvw, {}, v, s, seed);
    g_layer += geo::experimental_variogram(uvw, y, d).gamma[4] / 20.0;
    g_space += geo::experimental_variogram(xyz, y, d).gamma[4] / 20.0;
  }
  std::printf("  gamma at 50 m: along the layer %.3f (model %.3f); along x in present-day space %.3f\n", g_layer, v.gamma({50, 0, 0}), g_space);
  CHECK(std::fabs(g_layer - v.gamma({50, 0, 0})) < 0.1);
  (void)xyz;
}

// ---- 2. Grids and sensitivity (Prompt 60) --------------------------------------------------------

TEST_CASE("S9 fault-following grid: faces on the oblique fault within 1 % of the cell size") {
  auto d = desc::load_description(example("fault_system_relay_splay"));
  auto m = build_only(d, BuildOptions{}).model;
  auto o = grid_options(d, *m);
  o.follow_fault_traces = true;
  auto g = repr::corner_point_grid(*m, o);
  auto def = repr::grid_defects(*m, g, 0.01 * std::min(o.dx, o.dy));
  std::printf("  %zu fault faces, %zu off, %zu stair-step fall-backs, worst %.3f m\n", def.fault_faces, def.fault_face_off, def.fallback_faces,
              def.worst_fault_distance);
  CHECK(def.fault_faces > 0);
  CHECK(def.fault_face_off == 0);
  CHECK(def.inverted == 0);
  CHECK(def.nonpositive_volume == 0);
}

TEST_CASE("S10 total-effect indices of the Ishigami function within their intervals") {
  const std::size_t N = 4096;
  auto des = unc::saltelli_design(3, N, 2026);
  REQUIRE(des.rows.size() == N * 5);
  std::vector<double> y;
  for (auto& r : des.rows) {
    const double x1 = -kPiD + 2 * kPiD * r[0], x2 = -kPiD + 2 * kPiD * r[1], x3 = -kPiD + 2 * kPiD * r[2];
    y.push_back(std::sin(x1) + 7 * std::sin(x2) * std::sin(x2) + 0.1 * std::pow(x3, 4) * std::sin(x1));
  }
  auto rows = unc::total_effect_indices(des, y, {"x1", "x2", "x3"}, "f", 500, 0.95, 11);
  REQUIRE(rows.size() == 3);
  const double s1[3] = {0.3139, 0.4424, 0.0}, st[3] = {0.5576, 0.4424, 0.2437};
  for (std::size_t i = 0; i < 3; ++i) {
    std::printf("  %s: S1 %.3f [%.3f, %.3f], ST %.3f [%.3f, %.3f]\n", rows[i].parameter.c_str(), rows[i].s1, rows[i].s1_lo, rows[i].s1_hi, rows[i].st,
                rows[i].st_lo, rows[i].st_hi);
    CHECK(rows[i].runs == N * 5);
    CHECK(rows[i].st_lo <= st[i] && st[i] <= rows[i].st_hi);
    CHECK(std::fabs(rows[i].st - st[i]) < 0.03);
    CHECK(std::fabs(rows[i].s1 - s1[i]) < 0.05);
  }
}

// ---- 3. Speed and scale (Prompt 61) ----------------------------------------------------------------

TEST_CASE("S11 cached and uncached builds give the same bytes") {
  auto cache = make_build_cache();
  BuildOptions o;
  o.figures = false;
  const auto a = tmp("uncached"), b = tmp("cached1"), c = tmp("cached2");
  run_build(example("fault_system_relay_splay"), a, o);
  run_build_cached(example("fault_system_relay_splay"), b, o, *cache);
  run_build_cached(example("fault_system_relay_splay"), c, o, *cache);  // served from the cache
  std::size_t files = 0;
  for (auto& e : fs::recursive_directory_iterator(a)) {
    if (!e.is_regular_file() || e.path().filename() == "manifest.json") continue;
    const auto rel = fs::relative(e.path(), a);
    CHECK(file_bytes(e.path()) == file_bytes(b / rel));
    CHECK(file_bytes(e.path()) == file_bytes(c / rel));
    ++files;
  }
  CHECK(files > 3);
}

TEST_CASE("S12 parallel builds at 1, 2 and 8 threads give the same bytes") {
  auto d1 = parallel_build_digest(example("fault_system_relay_splay"), 1);
  auto d2 = parallel_build_digest(example("fault_system_relay_splay"), 2);
  auto d8 = parallel_build_digest(example("fault_system_relay_splay"), 8);
  CHECK(!d1.empty());
  CHECK(d1 == d2);
  CHECK(d1 == d8);
}

TEST_CASE("S13 a grid written as it is made equals the grid written from memory") {
  BuildOptions o;
  o.figures = false;
  o.exports = {"grdecl"};
  const auto a = tmp("in_memory"), b = tmp("streamed");
  run_build(example("fault_system_relay_splay"), a, o);
  repr::GridDefects none;
  (void)none;
  const auto peak = stream_grid(example("fault_system_relay_splay"), b, {64ull << 20, 2, {"grdecl"}});
  CHECK(peak > 0);
  std::size_t files = 0;
  for (auto& e : fs::recursive_directory_iterator(b))
    if (e.is_regular_file() && e.path().extension() == ".grdecl") {
      CHECK(file_bytes(e.path()) == file_bytes(a / "exports" / "grid" / e.path().filename()));
      ++files;
    }
  CHECK(files > 0);
}

STX_TEST_MAIN()
