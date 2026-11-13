// Implicit (data-conditioned) builder: solvers, RBF field, plane fit, and the synthetic example
// rebuilt from sampled data compared with the parametric model it was sampled from.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <random>

#include "../stx_test.hpp"
#include "stratonix/core/linalg.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/engines/implicit.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/model/model_suite.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

using namespace stratonix;
namespace fs = std::filesystem;

TEST_CASE("J1 blocked LU equals the reference solver; independent of the thread count") {
  for (std::size_t n : {7u, 130u, 300u}) {
    std::mt19937 rng(static_cast<unsigned>(n));
    std::uniform_real_distribution<double> u(-1, 1);
    std::vector<double> A(n * n), b(n);
    for (auto& x : A) x = u(rng);
    for (std::size_t i = 0; i < n; ++i) A[i * n + i] += (i % 3 == 0) ? 0.0 : 2.0;  // some weak diagonals: pivoting matters
    for (auto& x : b) x = u(rng);
    auto A1 = A, A2 = A, A3 = A;
    auto b1 = b, b2 = b, b3 = b;
    solve_dense_reference(A1, b1, n);
    set_thread_count(1);
    solve_dense_blocked(A2, b2, n, 32);
    set_thread_count(3);
    solve_dense_blocked(A3, b3, n, 32);
    set_thread_count(0);
    double worst = 0, res = 0;
    for (std::size_t i = 0; i < n; ++i) worst = std::max(worst, std::fabs(b1[i] - b2[i]) / (1 + std::fabs(b1[i])));
    for (std::size_t i = 0; i < n; ++i) {
      double r = -b[i];
      for (std::size_t j = 0; j < n; ++j) r += A[i * n + j] * b2[j];
      res = std::max(res, std::fabs(r));
    }
    CHECK(worst < 1e-9);
    CHECK(res < 1e-9);
    CHECK(b2 == b3);  // bitwise identical for any thread count
  }
}

TEST_CASE("J2 RBF field interpolates values and gradients and reproduces a plane") {
  std::vector<RbfField::Value> vs;
  std::vector<RbfField::Gradient> gs;
  // plane z = -900 + 0.05 x - 0.02 y  ->  s = z - (-900 + 0.05 x - 0.02 y) + const, unit normal of the plane
  Vec3 n = normalized(Vec3{-0.05, 0.02, 1});
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> ux(0, 5000), uy(0, 4000);
  for (int i = 0; i < 40; ++i) {  // two parallel planes 100 m apart, values -900 and -1000
    double x = ux(rng), y = uy(rng), off = (i % 2) ? -100.0 : 0.0;
    vs.push_back({{x, y, -900 + off + 0.05 * x - 0.02 * y}, -900 + off});
  }
  for (int i = 0; i < 10; ++i) {
    double x = ux(rng), y = uy(rng);
    gs.push_back({{x, y, -900 + 0.05 * x - 0.02 * y}, n * (1.0 / n.z)});  // gradient of s = z - plane(x, y)
  }
  RbfField f(vs, gs, 0, {2500, 2000, -900}, 2500);
  double rv = 0, rg = 0;
  for (auto& v : vs) rv = std::max(rv, std::fabs(f.value(v.p) - v.v));
  for (auto& g : gs) rg = std::max(rg, norm(cross(normalized(f.gradient(g.p)), normalized(g.g))));  // direction honoured
  CHECK(rv < 1e-6);
  CHECK(rg < 1e-6);
  // away from the data the field is the plane (linear drift reproduces it exactly)
  Vec3 q{4500, 300, -950};
  CHECK_NEAR(f.value(q), q.z - (0.05 * 4500 - 0.02 * 300), 1e-5);
}

TEST_CASE("J3 a plane fitted to fault sticks recovers strike and dip") {
  // plane through (1000, 2000, -1200), strike 30 deg, dip 60 deg (dip direction 120 deg)
  double st = deg2rad(30), dp = deg2rad(60);
  Vec3 s{std::sin(st), std::cos(st), 0}, dd{std::sin(st + kPi / 2) * std::cos(dp), std::cos(st + kPi / 2) * std::cos(dp), -std::sin(dp)};
  std::vector<Polyline> sticks;
  for (int k = -3; k <= 3; ++k) {
    Polyline l;
    for (int j = 0; j < 5; ++j) l.push_back(Vec3{1000, 2000, -1200} + s * (200.0 * k) + dd * (80.0 * j));
    sticks.push_back(l);
  }
  auto fp = fit_plane(sticks);
  CHECK_NEAR(rad2deg(fp.strike), 30, 1e-6);
  CHECK_NEAR(rad2deg(fp.dip), 60, 1e-6);
  CHECK(fp.rms < 1e-9);
}

TEST_CASE("J4 the synthetic example rebuilt from sampled data reproduces the parametric model") {
  register_all_plugins();
  auto dp = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples/synthetic_fold_two_faults.yaml");
  auto di = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples/implicit_from_samples.yaml");
  BuildLog l1, l2;
  BuildInputs in;
  io::load_build_inputs(di, in);
  auto mp = build_model(dp, {}, l1);
  auto mi = build_model(di, in, l2);
  CHECK(mi->info().engine == "implicit-rbf");
  CHECK(!mi->info().provenance.empty());
  CHECK(run_model_suite(*mi).passed());
  // misfit at every data point
  auto* mis = static_cast<const std::vector<DataMisfit>*>(mi->capability("data_misfit"));
  REQUIRE(mis != nullptr);
  double mp_max = 0, mp_sum = 0, mo_max = 0;
  std::size_t np = 0;
  for (auto& r : *mis) {
    CHECK(!r.evidence.empty());
    if (r.kind == "orientation") mo_max = std::max(mo_max, r.misfit);
    else { mp_max = std::max(mp_max, std::fabs(r.misfit)); mp_sum += std::fabs(r.misfit); ++np; }
  }
  std::printf("  data misfit: points mean %.4f m, max %.4f m; orientations max %.3f deg\n", mp_sum / np, mp_max, mo_max);
  CHECK(mp_sum / np < 0.05);
  CHECK(mp_max < 1.0);
  CHECK(mo_max < 2.0);
  // horizon depths against the parametric model, away from and near the faults
  GridSpec g = GridSpec::covering(mp->info().domain, 50);
  for (int h = 0; h < 5; ++h) {
    auto a = repr::horizon_grid_faulted(*mp, h, g), b = repr::horizon_grid_faulted(*mi, h, g);
    // distance (in nodes) to the nearest gap of the parametric grid
    std::vector<int> dist(g.size(), 1 << 20);
    for (std::size_t k = 0; k < g.size(); ++k)
      if (std::isnan(a.z.v[k])) dist[k] = 0;
    for (int pass = 0; pass < 2; ++pass)  // two-pass chamfer distance (city block)
      for (std::size_t jj = 0; jj < g.ny; ++jj)
        for (std::size_t ii = 0; ii < g.nx; ++ii) {
          std::size_t i = pass ? g.nx - 1 - ii : ii, j = pass ? g.ny - 1 - jj : jj;
          int& d = dist[g.index(i, j)];
          int di1 = pass ? 1 : -1;
          long ni = static_cast<long>(i) + di1, nj = static_cast<long>(j) + di1;
          if (ni >= 0 && ni < static_cast<long>(g.nx)) d = std::min(d, dist[g.index(static_cast<std::size_t>(ni), j)] + 1);
          if (nj >= 0 && nj < static_cast<long>(g.ny)) d = std::min(d, dist[g.index(i, static_cast<std::size_t>(nj))] + 1);
        }
    double away_sum = 0, away_max = 0, near_sum = 0, near_max = 0, edge_max = 0;
    std::size_t na = 0, nn = 0, mismatch = 0;
    const Domain& dm = mp->info().domain;
    for (std::size_t k = 0; k < g.size(); ++k) {
      if (std::isnan(a.z.v[k]) || std::isnan(b.z.v[k])) { mismatch += std::isnan(a.z.v[k]) != std::isnan(b.z.v[k]); continue; }
      double x = g.x(k % g.nx), y = g.y(k / g.nx);
      // outside the data footprint (data start 250 m inside the domain): extrapolation, reported apart
      if (x < dm.x0 + 250 || x > dm.x0 + dm.sx - 250 || y < dm.y0 + 250 || y > dm.y0 + dm.sy - 250) {
        edge_max = std::max(edge_max, std::fabs(a.z.v[k] - b.z.v[k]));
        continue;
      }
      double e = std::fabs(a.z.v[k] - b.z.v[k]);
      if (dist[k] > 10) { away_sum += e; away_max = std::max(away_max, e); ++na; }  // more than 500 m from a gap
      else { near_sum += e; near_max = std::max(near_max, e); ++nn; }
    }
    std::printf("  H%d: away mean %.3f max %.3f m | near mean %.3f max %.3f m | edge band max %.3f m | gap mismatches %zu\n", h + 1,
                away_sum / std::max<std::size_t>(na, 1), away_max, nn ? near_sum / nn : 0.0, near_max, edge_max, mismatch);
    CHECK(away_sum / std::max<std::size_t>(na, 1) < 0.2);
    CHECK(away_max < 10.0);
    if (nn) CHECK(near_sum / nn < 0.5);
  }
}

TEST_CASE("J5 implicit inputs are validated") {
  register_all_plugins();
  std::string y = read_file(fs::path(STX_SOURCE_DIR) / "examples/implicit_from_samples.yaml");
  auto p = y.find("{input: ORIENT, kind: orientations, format: csv_orientations}");
  REQUIRE(p != std::string::npos);
  y.replace(p, 61, "{input: MISSING, kind: orientations, format: csv_orientations}");
  auto rep = desc::validate_description_text(y, fs::path(STX_SOURCE_DIR) / "examples");
  bool i1 = false;
  for (auto& is : rep.issues) i1 |= is.rule == "I1";
  CHECK(i1);
}

STX_TEST_MAIN()
