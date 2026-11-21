// Stage 4A analytical tests (STX-TN-14 section 12), written before the features (Prompt 42) and
// made to pass by Prompts 43 to 47: free horizon values on a discrete field, fold frames, fault
// displacement fields and relations, the posterior sampler and Broyden least squares, geological
// time, topology and watertightness, gravity and magnetics, and joint inversion.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <random>

#include "v2_fixture.hpp"
#include "stratonix/engines/conditioning.hpp"
#include "stratonix/engines/fields.hpp"
#include "stratonix/engines/foldframe.hpp"
#include "stratonix/engines/geophysics.hpp"
#include "stratonix/engines/implicit.hpp"
#include "stratonix/engines/posterior.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/repr/geotime.hpp"
#include "stratonix/repr/topology.hpp"

namespace fs = std::filesystem;
using namespace stratonix;
using namespace v2fx;

namespace {
fs::path tmpdir(const std::string& name) {
  fs::path p = fs::path(STX_TMP_DIR) / "v4_engine" / name;
  fs::remove_all(p);
  fs::create_directories(p);
  return p;
}
void write_text(const fs::path& p, const std::string& s) { std::ofstream(p, std::ios::binary) << s; }
const char* kU = " | user_instruction | test";
std::string qq(const std::string& v) { return "\"" + v + kU + "\""; }

double percentile95(std::vector<double> v) {
  if (v.empty()) return kNaN;
  std::sort(v.begin(), v.end());
  return v[static_cast<std::size_t>(0.95 * static_cast<double>(v.size() - 1))];
}

// ---- the thinning series of decision 13 (STX-TN-14 section 9) --------------------------------------
constexpr double kL = 10000.0;
double dome(double x, double y) {
  double r2 = (x - 5000) * (x - 5000) + (y - 5000) * (y - 5000);
  return -1000.0 + 150.0 * std::exp(-r2 / (2 * 2500.0 * 2500.0));
}
double thin_truth(int k, double x, double y, double thin) {
  double z = dome(x, y);
  if (k >= 1) z -= 300.0 - thin * x / kL;
  if (k >= 2) z -= 250.0;
  return z;
}
Vec3 thin_normal(int k, double x, double y, double thin) {
  const double e = 1.0;
  double zx = (thin_truth(k, x + e, y, thin) - thin_truth(k, x - e, y, thin)) / (2 * e);
  double zy = (thin_truth(k, x, y + e, thin) - thin_truth(k, x, y - e, thin)) / (2 * e);
  return normalized(Vec3{-zx, -zy, 1.0});
}
struct ThinData {
  std::vector<HorizonDatum> pts;
  std::vector<RbfField::Gradient> grads;
};
ThinData thin_data(std::size_t n, double thin, unsigned seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> U(200.0, kL - 200.0);
  ThinData d;
  for (std::size_t i = 0; i < n; ++i) {
    int k = static_cast<int>(i % 3);
    double x = U(rng), y = U(rng);
    d.pts.push_back({{x, y, thin_truth(k, x, y, thin)}, k});
    if (i % 10 == 0) d.grads.push_back({{x, y, thin_truth(k, x, y, thin)}, thin_normal(k, x, y, thin)});
  }
  return d;
}
/// Depth error of horizon k of a field with level value v, on a 20 x 20 lattice of columns.
template <class F>
double depth_p95(const F& f, double v, int k, double thin) {
  std::vector<double> err;
  for (int j = 0; j < 20; ++j)
    for (int i = 0; i < 20; ++i) {
      double x = 500 + i * 9000.0 / 19, y = 500 + j * 9000.0 / 19;
      double za = -700, fa = f.value({x, y, za}) - v;
      bool found = false;
      for (int s = 1; s <= 440 && !found; ++s) {
        double zb = -700 - 2.5 * s, fb = f.value({x, y, zb}) - v;
        if ((fa >= 0) != (fb >= 0)) {
          double z = za + (zb - za) * fa / (fa - fb);
          err.push_back(std::fabs(z - thin_truth(k, x, y, thin)));
          found = true;
        }
        za = zb;
        fa = fb;
      }
      if (!found) err.push_back(1e9);
    }
  return percentile95(err);
}
}  // namespace

// ---- E1, E2: free horizon values ------------------------------------------------------------------
TEST_CASE("E1 discrete field: a layer-cake field is reproduced; free values follow a thinning unit where fixed values fail") {
  // (a) dipping layer cake: values recovered within 2 m
  {
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> U(200.0, 9800.0);
    std::vector<HorizonDatum> pts;
    std::vector<RbfField::Gradient> grads;
    auto z = [](int k, double x, double y) { return -1000.0 - 200.0 * k + 0.02 * x - 0.01 * y; };
    Vec3 n = normalized(Vec3{-0.02, 0.01, 1.0});
    for (int i = 0; i < 300; ++i) {
      int k = i % 3;
      double x = U(rng), y = U(rng);
      pts.push_back({{x, y, z(k, x, y)}, k});
      if (i % 10 == 0) grads.push_back({{x, y, z(k, x, y)}, n});
    }
    BBox box;
    box.add({0, 0, -1900});
    box.add({10000, 10000, -700});
    DiscreteOptions o;
    DiscreteField f(pts, grads, 0, -1000.0, box, o);
    auto& v = f.horizon_values();
    REQUIRE(v.size() == 3);
    std::printf("  layer cake: values %.2f %.2f %.2f, %d iterations, residual %.2e\n", v[0], v[1], v[2], f.iterations(), f.relative_residual());
    CHECK_NEAR(v[1] - v[0], -200.0, 2.0);
    CHECK_NEAR(v[2] - v[0], -400.0, 2.0);
    CHECK(f.converged());
    CHECK_NEAR(f.value({5000, 5000, z(1, 5000, 5000)}), v[1], 2.0);
  }
  // (b) the middle unit thins from 300 m to 10 m: decision 13's case
  {
    const double thin = 290;
    auto d = thin_data(600, thin, 20261121);
    BBox box;
    box.add({0, 0, -1800});
    box.add({kL, kL, -700});
    DiscreteOptions o;
    o.spacing = 250;
    o.vertical_spacing = 25;
    DiscreteField f(d.pts, d.grads, 0, thin_truth(0, 500, 5000, thin), box, o);
    auto& v = f.horizon_values();
    double e2 = depth_p95(f, v[2], 2, thin), e1 = depth_p95(f, v[1], 1, thin);
    std::printf("  thinning unit, free values: p95 H1 %.1f m, H2 %.1f m\n", e1, e2);
    CHECK(e1 <= 15.0);
    CHECK(e2 <= 15.0);
    // the 3A formulation with values fixed at the reference well distorts H2 (I-037)
    std::vector<RbfField::Value> vals;
    for (auto& p : d.pts) vals.push_back({p.p, thin_truth(p.horizon, 500, 5000, thin)});
    RbfField fixed(vals, d.grads, 0.0, {5000, 5000, -1250}, 5000, 5);
    double ef = depth_p95(fixed, thin_truth(2, 500, 5000, thin), 2, thin);
    std::printf("  thinning unit, fixed values: p95 H2 %.1f m\n", ef);
    CHECK(ef > 100.0);
  }
}

TEST_CASE("E2 free horizon values out of stratigraphic order are refused (F34); free values need the discrete solver (F35)") {
  CHECK(horizon_order_problem({-1000, -1200, -1400}, {"H1", "H2", "H3"}).empty());
  auto p = horizon_order_problem({-1000, -1450, -1400}, {"H1", "H2", "H3"});
  CHECK(p.find("H2") != std::string::npos);
  CHECK(p.find("H3") != std::string::npos);
  // F35: a description asking for free values with another solver is refused
  auto y = std::string("stratonix_schema: 1\nmodel: {name: t, version: \"1.0.0\", status: DESIGN}\n"
                       "crs: {name: local, vertical_datum: datum, length_unit: m}\n"
                       "domain:\n  origin: ") + qq("[0, 0] m") + "\n  size: " + qq("[10000, 10000] m") + "\n  depth_range: " +
           qq("[0, 3000] m") + "\n  grid_spacing: " + qq("100 m") +
           "\ninputs:\n  - {id: P, path: points.csv, role: surface points, sharing: shareable, source: test, evidence: user_instruction}\n"
           "stratigraphy:\n  horizons: [{name: H1}, {name: H2}]\n  units: [{name: U1, top: H1, base: H2}]\n"
           "structure:\n  builder: implicit-rbf\n"
           "implicit:\n  data:\n    - {input: P, kind: points, format: csv_points}\n  horizon_values: free\n  solver: {type: partition_of_unity}\n";
  auto dir = tmpdir("f35");
  write_text(dir / "points.csv", "X,Y,Z,surface\n1000,1000,-1000,H1\n2000,1000,-1000,H1\n1000,2000,-1200,H2\n3000,3000,-1200,H2\n");
  auto rep = desc::validate_description_text(y, dir, nullptr);
  bool f35 = false;
  for (auto& i : rep.issues) f35 = f35 || i.rule == "F35";
  if (!f35) std::cerr << rep.text();
  CHECK(f35);
}

// ---- E3: fold frames ----------------------------------------------------------------------------
TEST_CASE("E3 fold frame: axis within 2 degrees, interlimb angle within 5 degrees, wavelength within 10 %") {
  // a sine fold z' = A sin(2 pi x' / W) in a frame whose axis trends 030 and plunges 10 degrees
  const double A = 400, W = 4000, trend = deg2rad(30), plunge = deg2rad(10);
  Vec3 axis{std::sin(trend) * std::cos(plunge), std::cos(trend) * std::cos(plunge), -std::sin(plunge)};
  Vec3 across = normalized(cross(axis, Vec3{0, 0, 1}));  // horizontal, normal to the axis
  Vec3 up = normalized(cross(across, axis));
  std::mt19937_64 rng(11);
  std::uniform_real_distribution<double> U(-6000, 6000), V(-3000, 3000);
  std::vector<OrientationDatum> bed;
  for (int i = 0; i < 400; ++i) {
    double xp = U(rng), along = V(rng);
    double slope = A * 2 * kPi / W * std::cos(2 * kPi * xp / W);
    Vec3 nloc = normalized(up - across * slope);  // normal of z' = f(x') in the (across, up) plane
    Vec3 p = across * xp + axis * along + up * (A * std::sin(2 * kPi * xp / W));
    if (dot(nloc, Vec3{0, 0, 1}) < 0) nloc = nloc * -1.0;
    bed.push_back({p + Vec3{20000, 20000, -2000}, nloc});
  }
  FoldFrameOptions o;
  auto f = fit_fold_frame(bed, {}, o);
  REQUIRE(f.ok);
  double ang = rad2deg(std::acos(std::min(1.0, std::fabs(dot(f.axis, axis)))));
  const double interlimb_true = 180.0 - 2 * rad2deg(std::atan(A * 2 * kPi / W));
  std::printf("  axis error %.2f deg, interlimb %.1f deg (true %.1f), wavelength %.0f m (true %.0f), order %d, cylindricity %.3f\n",
              ang, rad2deg(f.interlimb), interlimb_true, f.wavelength, W, f.order, f.cylindricity);
  CHECK(ang <= 2.0);
  CHECK(std::fabs(rad2deg(f.interlimb) - interlimb_true) <= 5.0);
  CHECK(std::fabs(f.wavelength - W) <= 0.1 * W);
  // the frame predicts the bedding it was fitted to
  double worst = 0;
  for (int i = 0; i < 400; i += 7) worst = std::max(worst, rad2deg(std::acos(std::clamp(dot(f.bedding_normal(bed[i].p), bed[i].n), -1.0, 1.0))));
  std::printf("  worst predicted bedding normal error %.1f deg\n", worst);
  CHECK(worst <= 5.0);
}

// ---- E4, E5: faults from data ----------------------------------------------------------------------
TEST_CASE("E4 displacement field: maximum within 10 % or 5 m, centre within 300 m") {
  // T(s, d) = 100 m, elliptical along strike (length 6 km, centre s = 400 m) and down dip
  // (height 2 km, centre d = 300 m), observed at five horizons
  auto truth = [](double s, double d) {
    double a = 2 * (s - 400) / 6000, b = 2 * (d - 300) / 2000;
    double r = a * a + b * b;
    return r >= 1 ? 0.0 : 100.0 * std::sqrt(1 - a * a) * std::sqrt(1 - b * b);
  };
  std::vector<DisplacementObservation> obs;
  std::mt19937_64 rng(4);
  std::normal_distribution<double> N(0, 2.0);
  for (int h = 0; h < 5; ++h)
    for (double s = -2400; s <= 3200; s += 250) {
      double d = -300 + 200.0 * h;
      double t = truth(s, d);
      if (t > 0) obs.push_back({s, d, t + N(rng), h});
    }
  DisplacementOptions o;
  auto f = fit_displacement_field(obs, o);
  REQUIRE(f.ok);
  std::printf("  max %.1f m (100), s0 %.0f m (400), d0 %.0f m (300), length %.0f m (6000), height %.0f m (2000), rms %.2f m\n", f.max_throw,
              f.s0, f.d0, f.length, f.height, f.rms);
  CHECK(std::fabs(f.max_throw - 100) <= std::max(10.0, 5.0));
  CHECK(std::fabs(f.s0 - 400) <= 300);
  CHECK(std::fabs(f.at(400, 300) - 100) <= 10);
  CHECK(f.at(5000, 300) == 0.0);
}

TEST_CASE("E5 fault relations: abutting and splay are recovered from data and rules") {
  auto plane = [](Vec3 c, double strike_deg, double dip_deg) {
    FittedPlane p;
    p.centroid = c;
    p.strike = deg2rad(strike_deg);
    p.dip = deg2rad(dip_deg);
    Vec3 dd = azimuth_vector(p.strike + kPi / 2);  // dip direction
    p.normal = normalized(dd * std::sin(p.dip) + Vec3{0, 0, std::cos(p.dip)});
    return p;
  };
  // B: strike 000, dips east 60; A: strike 020, dips east 60, data stop 40 m short of B on its
  // west side (a splay); C: strike 090, dips south 60, data on both sides of B (crosses)
  FaultDataSet B{"B", plane({5000, 5000, -1500}, 0, 60), {}, 4000};
  FaultDataSet A{"A", plane({3500, 3000, -1500}, 20, 60), {}, 2500};
  FaultDataSet C{"C", plane({5000, 7000, -1500}, 90, 60), {}, 3000};
  auto on_plane = [](const FittedPlane& p, double s, double d) {
    Vec3 st = azimuth_vector(p.strike), dd = azimuth_vector(p.strike + kPi / 2);
    Vec3 down = normalized(dd * std::cos(p.dip) - Vec3{0, 0, std::sin(p.dip)});
    return p.centroid + st * s + down * d;
  };
  for (double d = -500; d <= 500; d += 100) {
    for (double s = -2000; s <= 2000; s += 200) B.points.push_back(on_plane(B.plane, s, d));
    for (double s = -1500; s <= 1500; s += 200) C.points.push_back(on_plane(C.plane, s, d));
  }
  // A's data: along strike towards B, stopping 40 m from B's surface
  for (double d = -500; d <= 500; d += 100)
    for (double s = -1500; s <= 3000; s += 100) {
      Vec3 p = on_plane(A.plane, s, d);
      double dist = dot(p - B.plane.centroid, B.plane.normal);
      if (dist < -40) A.points.push_back(p);
    }
  auto rel = infer_fault_relations({A, B, C}, 150);
  bool a_splay_b = false, c_crosses_b = false;
  for (auto& r : rel) {
    std::printf("  %s %s %s (%s)\n", r.a.c_str(), r.kind.c_str(), r.b.c_str(), r.rule.c_str());
    if (r.a == "A" && r.b == "B" && r.kind == "splay") a_splay_b = true;
    if (((r.a == "C" && r.b == "B") || (r.a == "B" && r.b == "C")) && r.kind == "crosses") c_crosses_b = true;
    CHECK(!r.rule.empty());
  }
  CHECK(a_splay_b);
  CHECK(c_crosses_b);
}

// ---- E6, E7, E8: posterior and least squares -------------------------------------------------------
namespace {
/// A correlated bivariate normal: mean (3, -2), sd (2, 0.5), correlation 0.8.
double gauss2(const std::vector<double>& x) {
  const double m0 = 3, m1 = -2, s0 = 2, s1 = 0.5, rho = 0.8;
  double a = (x[0] - m0) / s0, b = (x[1] - m1) / s1;
  return -0.5 * (a * a - 2 * rho * a * b + b * b) / (1 - rho * rho);
}
std::vector<std::vector<double>> starts(int n, unsigned seed) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> N(0, 0.1);
  std::vector<std::vector<double>> s;
  for (int i = 0; i < n; ++i) s.push_back({3 + N(rng), -2 + N(rng)});
  return s;
}
}  // namespace

TEST_CASE("E6 posterior: the ensemble sampler matches an analytical Gaussian posterior") {
  PosteriorOptions o;
  o.walkers = 32;
  o.steps = 4000;
  o.burn_in = 500;
  o.thin = 1;
  o.seed = 20261121;
  auto r = sample_ensemble(gauss2, starts(32, 1), o);
  REQUIRE(r.ok);
  std::printf("  mean %.3f %.3f (3, -2), sd %.3f %.3f (2, 0.5), rhat %.4f %.4f, ess %.0f %.0f, acceptance %.2f\n", r.mean[0], r.mean[1],
              r.sd[0], r.sd[1], r.rhat[0], r.rhat[1], r.ess[0], r.ess[1], r.acceptance);
  CHECK(std::fabs(r.mean[0] - 3) <= 0.1 * 2);
  CHECK(std::fabs(r.mean[1] + 2) <= 0.1 * 0.5);
  CHECK(std::fabs(r.sd[0] - 2) <= 0.1 * 2);
  CHECK(std::fabs(r.sd[1] - 0.5) <= 0.1 * 0.5);
  CHECK(r.converged);
  CHECK(r.acceptance > 0.2);
  CHECK(r.acceptance < 0.9);
  // R-hat flags chains that disagree
  std::vector<std::vector<double>> bad{{0, 0.1, -0.1, 0.05}, {5, 5.1, 4.9, 5.05}};
  CHECK(split_rhat(bad) > 1.5);
}

TEST_CASE("E7 posterior: the same seed gives the same samples with 1 and 4 threads") {
  PosteriorOptions o;
  o.walkers = 16;
  o.steps = 300;
  o.burn_in = 50;
  o.thin = 1;
  o.seed = 7;
  o.threads = 1;
  auto a = sample_ensemble(gauss2, starts(16, 2), o);
  o.threads = 4;
  auto b = sample_ensemble(gauss2, starts(16, 2), o);
  REQUIRE(a.ok);
  REQUIRE(b.ok);
  REQUIRE(a.chain.size() == b.chain.size());
  bool same = true;
  for (std::size_t i = 0; i < a.chain.size(); ++i) same = same && a.chain[i] == b.chain[i];
  CHECK(same);
  CHECK(a.log_posterior == b.log_posterior);
  CHECK(a.rhat.size() == 2);
  CHECK(a.ess.size() == 2);
}

namespace {
/// The stage 3A hybrid model: anticline A1 (relief on H1 and H2) and planar fault F1, with H1
/// picks sampled from the truth (relief 250 m, throw 120 m) with N(0, 3 m) noise.
std::string hybrid_yaml(double relief, double throw_m, const std::string& tail) {
  std::string structure = "structure:\n  regional: {type: flat}\n  primitives:\n    - type: anticline\n      name: A1\n      crest: " +
                          qq("[5000, 5000] m") + "\n      axis_azimuth: " + qq("0 deg") + "\n      plunge: " + qq("0 deg") +
                          "\n      half_length: " + qq("4000 m") + "\n      half_width_left: " + qq("2500 m") + "\n      half_width_right: " +
                          qq("2500 m") + "\n      relief:\n        H1: " + qq(fmt_fixed(relief, 3) + " m") + "\n        H2: " +
                          qq(fmt_fixed(relief, 3) + " m") + "\n";
  std::string faults = planar_fault("F1", "1", "[6000, 5000] m", "1200 m", "0 deg", "60 deg", "20000 m",
                                    "max: " + q(fmt_fixed(throw_m, 3) + " m"), "dip_slip");
  return model_yaml(faults, structure, {1000, 1200, 1400}, tail);
}
void write_hybrid_picks(const fs::path& dir) {
  auto truth = build_text(hybrid_yaml(250, 120, ""));
  std::mt19937_64 rng(77);
  std::uniform_real_distribution<double> U(1500, 8500);
  std::normal_distribution<double> N(0, 3.0);
  std::string s = "X,Y,Z,surface,sigma\n";
  int n = 0;
  while (n < 240) {
    double x = U(rng), y = U(rng);
    double z = -depth_faulted(*truth, 0, x, y);
    if (std::isnan(z)) continue;
    s += fmt_fixed(x, 2) + "," + fmt_fixed(y, 2) + "," + fmt_fixed(z + N(rng), 3) + ",H1,3.0\n";
    ++n;
  }
  write_text(dir / "data.csv", s);
}
std::string hybrid_conditioning(const std::string& method) {
  return "inputs:\n  - {id: D, path: data.csv, role: surface points, sharing: shareable, source: test, evidence: interpreted}\n"
         "conditioning:\n  data:\n    - {input: D, kind: points, format: csv_points}\n  use_well_tops: false\n  uncertainty: " +
         qq("3 m") + "\n  parameters:\n    - {path: faults.F1.throw.max, min: " + qq("20 m") + ", max: " + qq("300 m") +
         "}\n    - {path: structure.primitives.A1.relief.H1, min: " + qq("50 m") + ", max: " + qq("500 m") + "}\n  max_iterations: " +
         qq("50 ordinal") + "\n" + method;
}
}  // namespace

TEST_CASE("E8 least squares with Broyden updates: one model build per iteration after the first Jacobian") {
  auto dir = tmpdir("broyden");
  write_hybrid_picks(dir);
  desc::Description d;
  auto rep = desc::validate_description_text(hybrid_yaml(175, 84, hybrid_conditioning("")), dir, &d);
  REQUIRE(rep.ok());
  BuildInputs in;
  io::load_build_inputs(d, in);
  BuildLog log;
  auto r = condition_to_data(d, in, log);
  std::printf("  converged %d, %d iterations, %ld builds, %d Jacobians\n", r.converged, r.iterations, r.builds, r.jacobians);
  REQUIRE(r.converged);
  const long np = 2;
  CHECK(r.jacobians >= 1);
  CHECK(r.builds <= np * r.jacobians + r.iterations + 3);
  CHECK(std::fabs(r.parameters[0].value - 120) <= 2 * r.parameters[0].sigma + 0.5);
}

#if defined(NDEBUG) && !defined(__SANITIZE_ADDRESS__)
TEST_CASE("E6b posterior of the hybrid model: truth within the 95 % credible interval, converged (release builds)") {
  auto dir = tmpdir("posterior");
  write_hybrid_picks(dir);
  desc::Description d;
  auto rep = desc::validate_description_text(
      hybrid_yaml(175, 84, hybrid_conditioning("  method: bayesian\n  sampler: {walkers: 16, steps: 1500, burn_in: 300, thin: 10, seed: 20261121}\n")),
      dir, &d);
  if (!rep.ok()) std::cerr << rep.text();
  REQUIRE(rep.ok());
  BuildInputs in;
  io::load_build_inputs(d, in);
  BuildLog log;
  PosteriorOptions o;
  o.walkers = 16;
  o.steps = 1500;
  o.burn_in = 300;
  o.thin = 10;
  o.seed = 20261121;
  auto r = sample_posterior(d, in, o, log);
  REQUIRE(r.ok);
  const double truth[] = {120, 250};
  for (std::size_t k = 0; k < 2; ++k) {
    std::printf("  %s: mean %.2f sd %.2f [%.2f, %.2f] rhat %.4f ess %.0f\n", r.paths[k].c_str(), r.mean[k], r.sd[k], r.q025[k], r.q975[k],
                r.rhat[k], r.ess[k]);
    CHECK(truth[k] >= r.q025[k]);
    CHECK(truth[k] <= r.q975[k]);
  }
  CHECK(r.converged);
  CHECK(r.builds == static_cast<long>(o.walkers) * (o.steps + o.burn_in) + o.walkers);
}
#endif

// ---- E9, E10: geological time and topology ---------------------------------------------------------
namespace {
std::string faulted_layercake() {
  std::string faults = planar_fault("F1", "1", "[5000, 5000] m", "1200 m", "0 deg", "60 deg", "40000 m", "max: " + q("100 m"), "dip_slip");
  return model_yaml(faults, "", {1000, 1200, 1400, 1600});
}
}  // namespace

TEST_CASE("E9 geological time: horizons map to one time each; the fault offset vanishes; the round trip returns within 0.5 m") {
  auto m = build_text(faulted_layercake());
  GeologicalTime gt(*m, {10, 20, 30, 40});
  for (int h = 0; h < 4; ++h)
    for (double x : {2000.0, 4000.0, 7000.0, 9000.0}) {
      double z = -depth_faulted(*m, h, x, 3000);
      if (std::isnan(z)) continue;
      auto p = gt.forward({x, 3000, z - 0.01});
      REQUIRE(p.defined);
      CHECK_NEAR(p.t, 10.0 * (h + 1), 0.05);
    }
  // a point between two horizons: interpolated; the round trip
  for (double x : {2500.0, 8000.0}) {
    double za = -depth_faulted(*m, 1, x, 6000), zb = -depth_faulted(*m, 2, x, 6000);
    Vec3 p{x, 6000, 0.5 * (za + zb)};
    auto g = gt.forward(p);
    REQUIRE(g.defined);
    CHECK_NEAR(g.t, 25.0, 0.5);
    auto back = gt.inverse(g.u, g.v, g.t);
    REQUIRE(back.has_value());
    CHECK(norm(*back - p) <= 0.5);
  }
  CHECK(!gt.ordinal());
}

TEST_CASE("E10 topology: the contact graph of a faulted layer cake; C26 passes; watertight within tolerance") {
  auto m = build_text(faulted_layercake());
  auto g = contact_graph(*m, 100);
  std::printf("  %zu nodes, %zu edges\n", g.nodes.size(), g.edges.size());
  // five units in two fault blocks
  CHECK(g.nodes.size() == 10);
  int horizon_edges = 0, fault_edges = 0;
  for (auto& e : g.edges) {
    if (e.surface.rfind("horizon:", 0) == 0) ++horizon_edges;
    if (e.surface.rfind("fault:", 0) == 0) ++fault_edges;
  }
  CHECK(horizon_edges == 8);  // four horizons in each block
  CHECK(fault_edges >= 5);    // each unit meets at least its own counterpart across F1
  CHECK(same_topology(g, contact_graph(*m, 100)));
  auto w = watertightness(*m, 100, 50);
  std::printf("  watertight: max gap %.2f m (tolerance %.0f m)\n", w.max_gap, w.tolerance);
  CHECK(w.ok());
  desc::Description d;
  validate(faulted_layercake(), &d);
  auto rep = checks_for(*m, d);
  auto* c26 = find_check(rep, "C26");
  REQUIRE(c26 != nullptr);
  CHECK(c26->passed);
  auto* c27 = find_check(rep, "C27");
  REQUIRE(c27 != nullptr);
  CHECK(c27->passed);
}

// ---- E11, E12, E13: gravity and magnetics ----------------------------------------------------------
TEST_CASE("E11 gravity: the prism closed form, and a sphere discretised into prisms against the point mass within 2 %") {
  // closed form against a numerical sum of small sub-prisms (point masses)
  Vec3 obs{30, -20, 10};
  double g = prism_gz(obs, -100, 100, -50, 150, -400, -200, 500);
  double sum = 0;
  const double G = 6.674e-11;
  const int n = 40;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k) {
        double x = -100 + (i + 0.5) * 200.0 / n, y = -50 + (j + 0.5) * 200.0 / n, z = -400 + (k + 0.5) * 200.0 / n;
        double dm = 500 * (200.0 / n) * (200.0 / n) * (200.0 / n);
        Vec3 r{x - obs.x, y - obs.y, z - obs.z};
        double d = norm(r);
        sum += G * dm * (-r.z) / (d * d * d) * 1e5;  // mGal, positive down
      }
  std::printf("  prism %.6f mGal, numerical %.6f mGal\n", g, sum);
  CHECK(std::fabs(g - sum) <= 0.002 * std::fabs(sum));
  // a sphere of radius 500 m, centre 1500 m deep, density contrast 300 kg/m^3, in 50 m prisms
  const double R = 500, zc = -1500, rho = 300, c = 50;
  double gs = 0;
  for (double x = -R; x < R; x += c)
    for (double y = -R; y < R; y += c)
      for (double z = zc - R; z < zc + R; z += c) {
        Vec3 m{x + c / 2, y + c / 2, z + c / 2};
        if (norm(m - Vec3{0, 0, zc}) <= R) gs += prism_gz({400, 0, 0}, x, x + c, y, y + c, z, z + c, rho);
      }
  double M = rho * 4.0 / 3.0 * kPi * R * R * R;
  double dx = 400, dz = 1500, rr = std::sqrt(dx * dx + dz * dz);
  double gp = G * M * dz / (rr * rr * rr) * 1e5;
  std::printf("  sphere %.5f mGal, point mass %.5f mGal\n", gs, gp);
  CHECK(std::fabs(gs - gp) <= 0.02 * gp);
}

TEST_CASE("E12 magnetics: a sphere discretised into prisms against the dipole within 3 %") {
  const double R = 400, zc = -1200, c = 40, k = 0.01, B0 = 50000;  // nT
  Vec3 F = field_direction(deg2rad(60), deg2rad(10));
  Vec3 Mv = F * (k * B0 / (4e-7 * kPi) * 1e-9);  // induced magnetisation (A/m)
  Vec3 obs{300, 200, 0};
  double dt = 0;
  for (double x = -R; x < R; x += c)
    for (double y = -R; y < R; y += c)
      for (double z = zc - R; z < zc + R; z += c) {
        Vec3 m{x + c / 2, y + c / 2, z + c / 2};
        if (norm(m - Vec3{0, 0, zc}) <= R) dt += prism_dt(obs, x, x + c, y, y + c, z, z + c, Mv, F);
      }
  // dipole: B = mu0/(4 pi) (3 (m . r) r / r^5 - m / r^3), m = M V
  Vec3 mom = Mv * (4.0 / 3.0 * kPi * R * R * R);
  Vec3 r = obs - Vec3{0, 0, zc};
  double rn = norm(r);
  Vec3 B = (r * (3 * dot(mom, r) / std::pow(rn, 5)) - mom / std::pow(rn, 3)) * (1e-7 * 1e9);
  double dip = dot(B, F);
  std::printf("  sphere %.4f nT, dipole %.4f nT\n", dt, dip);
  CHECK(std::fabs(dt - dip) <= 0.03 * std::fabs(dip));
}

#if defined(NDEBUG) && !defined(__SANITIZE_ADDRESS__)
TEST_CASE("E13 joint inversion: a fold relief is recovered from its own gravity anomaly within 5 % (release builds)") {
  auto dir = tmpdir("gravity");
  const std::string geo = "geophysics:\n  cell: " + qq("200 m") + "\n  background_density: " + qq("2650 kg/m3") +
                          "\n  units:\n    - {unit: U1, density: " + qq("2300 kg/m3") + "}\n    - {unit: U2, density: " + qq("2500 kg/m3") +
                          "}\n  gravity: {input: G, format: csv_stations}\n";
  // truth: relief 250 m; stations on a 1 km lattice at the surface; units not listed take the
  // background density
  auto truth = build_text(hybrid_yaml(250, 120, ""));
  PotentialFieldSpec s;
  s.cell = 200;
  s.background_density = 2650;
  for (auto& u : truth->info().units) s.density.push_back(u.name == "U1" ? 2300.0 : (u.name == "U2" ? 2500.0 : kNaN));
  std::vector<Vec3> st;
  for (double x = 1000; x <= 9000; x += 1000)
    for (double y = 1000; y <= 9000; y += 1000) st.push_back({x, y, 0});
  auto gz = gravity_anomaly(*truth, s, st);
  std::string csv = "X,Y,Z,value,sigma\n";
  for (std::size_t i = 0; i < st.size(); ++i)
    csv += fmt_fixed(st[i].x, 1) + "," + fmt_fixed(st[i].y, 1) + ",0," + fmt_fixed(gz[i], 5) + ",0.01\n";
  write_text(dir / "gravity.csv", csv);
  std::string tail = "inputs:\n  - {id: G, path: gravity.csv, role: gravity stations, sharing: shareable, source: test, evidence: interpreted}\n" + geo +
                     "conditioning:\n  data:\n    - {input: G, kind: gravity, format: csv_stations}\n  use_well_tops: false\n"
                     "  parameters:\n    - {path: structure.primitives.A1.relief.H1, min: " + qq("50 m") + ", max: " + qq("500 m") +
                     "}\n    - {path: structure.primitives.A1.relief.H2, min: " + qq("50 m") + ", max: " + qq("500 m") + "}\n";
  desc::Description d;
  auto rep = desc::validate_description_text(hybrid_yaml(175, 120, tail), dir, &d);
  if (!rep.ok()) std::cerr << rep.text();
  REQUIRE(rep.ok());
  BuildInputs in;
  io::load_build_inputs(d, in);
  BuildLog log;
  auto r = condition_to_data(d, in, log);
  std::printf("  converged %d: relief H1 %.1f m, H2 %.1f m (truth 250)\n", r.converged, r.parameters.empty() ? kNaN : r.parameters[0].value,
              r.parameters.size() < 2 ? kNaN : r.parameters[1].value);
  REQUIRE(r.converged);
  CHECK(std::fabs(r.parameters[0].value - 250) <= 0.05 * 250);
}
#endif

STX_TEST_MAIN()
