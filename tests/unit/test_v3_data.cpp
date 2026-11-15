// Stage 3A analytical tests (STX-TN-11 section 10), written before the features (Prompt 34) and
// made to pass by Prompts 35 and 36: the partition-of-unity solver, leave-one-out errors,
// per-point uncertainty, series from data, throws and fault profiles from data, and hybrid
// conditioning with its refusals.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <random>

#include "v2_fixture.hpp"
#include "stratonix/engines/conditioning.hpp"
#include "stratonix/engines/implicit.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/model/model_suite.hpp"

namespace fs = std::filesystem;
using namespace stratonix;
using namespace v2fx;

namespace {
// ---- analytic fold: z_h(x, y) = -(1000 + 200 h) + 150 exp(-r^2 / (2 * 2000^2)), r from (5000, 5000)
double fold_z(int h, double x, double y) {
  double r2 = (x - 5000) * (x - 5000) + (y - 5000) * (y - 5000);
  return -(1000.0 + 200.0 * h) + 150.0 * std::exp(-r2 / (2 * 2000.0 * 2000.0));
}
Vec3 fold_normal(int h, double x, double y) {
  const double e = 1.0;
  double zx = (fold_z(h, x + e, y) - fold_z(h, x - e, y)) / (2 * e), zy = (fold_z(h, x, y + e) - fold_z(h, x, y - e)) / (2 * e);
  return normalized(Vec3{-zx, -zy, 1.0});
}
struct FoldData {
  std::vector<RbfField::Value> values;
  std::vector<RbfField::Gradient> grads;
};
/// n interface points spread over 3 horizons (field value = the horizon's elevation at the
/// centre, as the builder uses the reference well), orientations at every tenth.
FoldData fold_data(std::size_t n, unsigned seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> U(200.0, 9800.0);
  FoldData d;
  for (std::size_t i = 0; i < n; ++i) {
    int h = static_cast<int>(i % 3);
    double x = U(rng), y = U(rng);
    d.values.push_back({{x, y, fold_z(h, x, y)}, fold_z(h, 5000, 5000)});
    if (i % 10 == 0) d.grads.push_back({{x, y, fold_z(h, x, y)}, fold_normal(h, x, y)});
  }
  return d;
}
const Vec3 kCentre{5000, 5000, -1250};
const double kScale = 5000, kStretch = 5;

/// Elevation where field f crosses value v on the column at (x, y), scanning down from 0 m.
template <class F>
double crossing(const F& f, double v, double x, double y) {
  double za = 0, fa = f.value({x, y, za}) - v;
  for (int k = 1; k <= 120; ++k) {
    double zb = -k * 25.0, fb = f.value({x, y, zb}) - v;
    if ((fa >= 0) != (fb >= 0)) {
      for (int it = 0; it < 60; ++it) {
        double zm = 0.5 * (za + zb), fm = f.value({x, y, zm}) - v;
        if ((fm >= 0) == (fa >= 0)) { za = zm; fa = fm; } else zb = zm;
      }
      return 0.5 * (za + zb);
    }
    za = zb;
    fa = fb;
  }
  return std::nan("");
}

double percentile95(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return v.empty() ? std::nan("") : v[static_cast<std::size_t>(0.95 * static_cast<double>(v.size() - 1))];
}

fs::path tmpdir(const std::string& name) {
  fs::path p = fs::path(stxtest::tmp_dir()) / ("v3_" + name);
  fs::remove_all(p);
  fs::create_directories(p);
  return p;
}
void write_text(const fs::path& p, const std::string& s) { std::ofstream(p, std::ios::binary) << s; }

ModelPtr build_with_inputs(const std::string& text, const fs::path& dir, desc::Description* dout = nullptr, BuildLog* log = nullptr) {
  register_all_plugins();
  desc::Description d;
  auto rep = desc::validate_description_text(text, dir, &d);
  if (!rep.ok()) {
    std::cerr << "  refused by validation (expected until stage 3A is built):\n" << rep.text();
    throw std::runtime_error("description refused");
  }
  BuildInputs in;
  io::load_build_inputs(d, in);
  BuildLog l;
  auto m = build_model(d, in, l);
  if (dout) *dout = d;
  if (log) *log = l;
  return m;
}
const char* kU = " | user_instruction | test";
std::string qq(const std::string& v) { return "\"" + v + kU + "\""; }
}  // namespace

TEST_CASE("V1 partition of unity: one patch is the dense field; many patches interpolate the data and agree with it") {
  auto d = fold_data(2000, 11);
  RbfField dense(d.values, d.grads, 0, kCentre, kScale, kStretch);
  SolverOptions one;
  one.type = "partition_of_unity";
  one.points_per_patch = 1000000;
  PuField p1(d.values, d.grads, 0, kCentre, kScale, kStretch, 0, 10000, 0, 10000, one);
  CHECK(p1.patches() == 1);
  std::mt19937_64 rng(5);
  std::uniform_real_distribution<double> U(0, 10000), Z(-2000, -500);
  double worst = 0;
  for (int k = 0; k < 200; ++k) {
    Vec3 x{U(rng), U(rng), Z(rng)};
    double a = dense.value(x), b = p1.value(x);
    worst = std::max(worst, std::fabs(a - b) / std::max(1.0, std::fabs(a)));
  }
  std::printf("  one patch against dense: worst relative difference %.3e\n", worst);
  CHECK(worst < 1e-9);
  SolverOptions many;
  many.type = "partition_of_unity";
  many.points_per_patch = 300;
  PuField pm(d.values, d.grads, 0, kCentre, kScale, kStretch, 0, 10000, 0, 10000, many);
  std::printf("  %zu patches, largest system %zu unknowns\n", pm.patches(), pm.largest_patch());
  CHECK(pm.patches() > 4);
  double miss = 0;
  for (auto& v : d.values) miss = std::max(miss, std::fabs(pm.value(v.p) - v.v));
  std::printf("  multi-patch: largest miss at a datum %.3e m\n", miss);
  CHECK(miss < 1e-6);
  std::vector<double> diff;
  for (int j = 0; j <= 20; ++j)
    for (int i = 0; i <= 20; ++i)
      for (int h = 0; h < 3; ++h) {
        double x = 500 + 450.0 * i, y = 500 + 450.0 * j, v = fold_z(h, 5000, 5000);
        diff.push_back(std::fabs(crossing(pm, v, x, y) - crossing(dense, v, x, y)));
      }
  std::printf("  surfaces, multi-patch against dense: p95 %.3f m\n", percentile95(diff));
  CHECK(percentile95(diff) <= 0.5);
}

TEST_CASE("V2 leave-one-out errors by Rippa's formula equal refitting without each datum") {
  auto d = fold_data(60, 21);
  auto loo = rippa_loo(d.values, d.grads, 0, kCentre, kScale, kStretch);
  REQUIRE(loo.size() == d.values.size());
  double worst = 0;
  for (std::size_t i = 0; i < d.values.size(); ++i) {
    auto vs = d.values;
    vs.erase(vs.begin() + static_cast<long>(i));
    RbfField f(vs, d.grads, 0, kCentre, kScale, kStretch);
    double brute = d.values[i].v - f.value(d.values[i].p);
    worst = std::max(worst, std::fabs(loo[i] - brute) / (1.0 + std::fabs(brute)));
  }
  std::printf("  worst difference to refitting: %.3e\n", worst);
  CHECK(worst < 1e-6);
  // the partition of unity with one patch gives the same errors
  SolverOptions one;
  one.type = "partition_of_unity";
  one.points_per_patch = 1000000;
  PuField p1(d.values, d.grads, 0, kCentre, kScale, kStretch, 0, 10000, 0, 10000, one);
  auto l1 = p1.leave_one_out();
  REQUIRE(l1.size() == loo.size());
  double w1 = 0;
  for (std::size_t i = 0; i < loo.size(); ++i) w1 = std::max(w1, std::fabs(l1[i] - loo[i]) / (1.0 + std::fabs(loo[i])));
  CHECK(w1 < 1e-6);
}

TEST_CASE("V3 per-point uncertainty: exact where sigma is 0, smoothed within sigma where it is not") {
  auto d = fold_data(300, 31);
  // displace one point by 15 m and give it sigma 20 m; pin another with sigma 0
  d.values[7].p.z += 15.0;
  d.values[7].sigma = 20.0;
  d.values[8].sigma = 0.0;
  RbfField f(d.values, d.grads, 0, kCentre, kScale, kStretch);
  double r7 = d.values[7].v - f.value(d.values[7].p), r8 = d.values[8].v - f.value(d.values[8].p);
  std::printf("  displaced point (sigma 20 m): field misses it by %.3f (field units, m); pinned point: %.3e\n", r7, r8);
  CHECK(std::fabs(r7) > 1.0);
  CHECK(std::fabs(r7) < 20.0);
  CHECK(std::fabs(r8) < 1e-6);
}

TEST_CASE("V4 series from data: one field per series, an erosional base truncates the older series") {
  // S2 (younger): H1 flat at 800 m, H2 flat at 1200 m (erosional base); S1 (older): H3, H4
  // dipping east at 0.1 (H3 at 1400 m and H4 at 1700 m at x = 5000), truncated by H2 in the west.
  auto zt = [](int h, double x) {
    switch (h) {
      case 0: return -800.0;
      case 1: return -1200.0;
      case 2: return -1400.0 - 0.1 * (x - 5000);
      default: return -1700.0 - 0.1 * (x - 5000);
    }
  };
  auto dir = tmpdir("series");
  std::string pts = "X,Y,Z,surface\n", ori = "X,Y,Z,G_x,G_y,G_z,surface\n";
  const char* hn[] = {"H1", "H2", "H3", "H4"};
  for (int h = 0; h < 4; ++h)
    for (double y = 500; y < 10000; y += 1000)
      for (double x = 500; x < 10000; x += 1000) {
        double z = zt(h, x);
        if (h >= 2 && z > -1200.0) continue;  // truncated: no data above the unconformity
        pts += fmt_fixed(x, 1) + "," + fmt_fixed(y, 1) + "," + fmt_fixed(z, 3) + "," + hn[h] + "\n";
        if (static_cast<int>(x) % 2000 == 500)
          ori += fmt_fixed(x, 1) + "," + fmt_fixed(y, 1) + "," + fmt_fixed(z, 3) + (h >= 2 ? ",0.0995037,0,0.9950372," : ",0,0,1,") + hn[h] + "\n";
      }
  write_text(dir / "points.csv", pts);
  write_text(dir / "orient.csv", ori);
  std::string y = "stratonix_schema: 1\nmodel: {name: t, version: \"1.0.0\", status: DESIGN}\ncrs: {name: local, vertical_datum: datum, length_unit: m}\n"
                  "domain:\n  origin: " + qq("[0, 0] m") + "\n  size: " + qq("[10000, 10000] m") + "\n  depth_range: " + qq("[0, 3000] m") +
                  "\n  grid_spacing: " + qq("100 m") + "\n"
                  "inputs:\n  - {id: P, path: points.csv, role: surface points, sharing: shareable, source: test, evidence: user_instruction}\n"
                  "  - {id: O, path: orient.csv, role: orientations, sharing: shareable, source: test, evidence: user_instruction}\n"
                  "wells:\n  - {name: W1, role: control, head: " + qq("[8000, 5000] m") + ", tops: {H1: " + qq("800 m") + ", H2: " + qq("1200 m") +
                  ", H3: " + qq("1700 m") + ", H4: " + qq("2000 m") + "}}\n"
                  "stratigraphy:\n  reference: {well: W1}\n  horizons: [{name: H1}, {name: H2}, {name: H3}, {name: H4}]\n"
                  "  units: [{name: U1, top: H1, base: H2}, {name: U2, top: H2, base: H3}, {name: U3, top: H3, base: H4}]\n"
                  "  series:\n    - {name: S2, horizons: [H1, H2], base: erosional}\n    - {name: S1, horizons: [H3, H4], base: conformable}\n"
                  "implicit:\n  data:\n    - {input: P, kind: points, format: csv_points}\n    - {input: O, kind: orientations, format: csv_orientations}\n"
                  "  use_well_tops: true\n  surface_spacing: " + qq("100 m") + "\n";
  desc::Description d;
  BuildLog log;
  auto m = build_with_inputs(y, dir, &d, &log);
  CHECK(log_value(log, "implicit.fields") == 2);
  double worst = 0;
  for (double x = 4000; x <= 9000; x += 500)
    for (int h = 2; h < 4; ++h) worst = std::max(worst, std::fabs(depth_faulted(*m, h, x, 5000) + zt(h, x)));
  std::printf("  older series where present: worst depth error %.3f m\n", worst);
  CHECK(worst <= 2.0);
  auto* sq = static_cast<const StratigraphyQuery*>(m->capability("stratigraphy"));
  REQUIRE(sq != nullptr);
  CHECK(sq->horizon_absent(2, 1000, 5000));   // H3 truncated in the west
  CHECK(!sq->horizon_absent(2, 8000, 5000));
  CHECK(std::fabs(depth_faulted(*m, 1, 1000, 5000) - 1200.0) <= 2.0);
}

TEST_CASE("V5 throws, length and centre estimated from horizon cut-offs across a fault") {
  // Fault striking north, dipping 60 degrees east, through (5000, 5000, -1500). Elliptical throw
  // profile: T_max 80 m, length 4000 m, centre 500 m north of the centroid.
  FittedPlane pl;
  pl.centroid = {5000, 5000, -1500};
  pl.strike = 0;
  pl.dip = deg2rad(60.0);
  pl.normal = normalized(Vec3{std::sin(pl.dip), 0, std::cos(pl.dip)});
  auto T = [](double s) {
    double u = 2 * (s - 500) / 4000;
    return std::fabs(u) < 1 ? 80.0 * std::sqrt(1 - u * u) : 0.0;
  };
  std::vector<HorizonPoint> pts;
  for (int h = 0; h < 2; ++h) {
    double zf = h == 0 ? -1200.0 : -1500.0;
    for (double y = 1500; y <= 8500; y += 100)
      for (double x = 3500; x <= 6500; x += 100) {
        double xf = 5000 + (-1500 - zf) / std::tan(pl.dip);  // fault x at the footwall level
        bool hw = x > xf;
        pts.push_back({{x, y, hw ? zf - T(y - 5000) : zf}, h});
      }
  }
  ThrowOptions o;
  o.band = 800;
  o.station_spacing = 250;
  o.profile = "elliptical";
  auto est = estimate_throws(pl, pts, 2, o);
  REQUIRE(est.size() == 2);
  for (auto& e : est) {
    std::printf("  horizon %d: T_max %.2f m, length %.1f m, centre %.1f m, RMS %.3f m, %zu stations\n", e.horizon, e.max_throw, e.length,
                e.centre_s, e.rms, e.stations.size());
    CHECK(e.ok);
    CHECK(std::fabs(e.max_throw - 80) <= 0.05 * 80);
    CHECK(std::fabs(e.length - 4000) <= 0.10 * 4000);
    CHECK(std::fabs(e.centre_s - 500) <= 200);
  }
}

TEST_CASE("V6 the dip profile of a listric fault fitted to its sticks") {
  // dip 60 degrees at the surface decreasing linearly to 20 degrees at 2000 m; strike 30 degrees
  auto dip = [](double z) { return deg2rad(60.0 - 20.0 * z / 1000.0); };
  const double strike = deg2rad(30.0);
  Vec3 sdir{std::sin(strike), std::cos(strike), 0}, ddir{std::cos(strike), -std::sin(strike), 0};
  std::vector<Polyline> sticks;
  for (int k = 0; k < 6; ++k) {
    Polyline st;
    double u = 0;
    for (double z = 0; z <= 2000.0 + 1e-9; z += 10.0) {
      if (static_cast<int>(z) % 100 == 0 && z >= 100) st.push_back(Vec3{4000, 4000, 0} + sdir * (500.0 * k) + ddir * u + Vec3{0, 0, -z});
      u += 10.0 / std::tan(dip(z + 5.0));
    }
    sticks.push_back(st);
  }
  auto pr = profile_from_sticks(sticks);
  std::printf("  strike %.2f deg, RMS %.3f m, %zu profile nodes\n", rad2deg(pr.strike), pr.rms, pr.depth.size());
  CHECK(std::fabs(rad2deg(pr.strike) - 30.0) < 1.0);
  for (double z : {500.0, 1000.0, 1500.0}) {
    // dip at depth z from the profile (linear interpolation)
    double v = std::nan("");
    for (std::size_t i = 0; i + 1 < pr.depth.size(); ++i)
      if (pr.depth[i] <= z && z <= pr.depth[i + 1])
        v = pr.dip[i] + (pr.dip[i + 1] - pr.dip[i]) * (z - pr.depth[i]) / (pr.depth[i + 1] - pr.depth[i]);
    std::printf("  dip at %.0f m: %.2f deg (truth %.2f)\n", z, rad2deg(v), rad2deg(dip(z)));
    CHECK(std::fabs(rad2deg(v) - rad2deg(dip(z))) <= 2.0);
  }
}

namespace {
/// The model whose parameters V7 and V8 recover: an anticline (A1, relief 250 m on H1 and H2) cut
/// by a planar normal fault F1 (throw 120 m), with conditioning data from `data` and the given
/// parameter lines.
std::string hybrid_yaml(double relief, double throw_m, const std::string& conditioning) {
  std::string structure = "structure:\n  regional: {type: flat}\n  primitives:\n    - type: anticline\n      name: A1\n      crest: " +
                          qq("[5000, 5000] m") + "\n      axis_azimuth: " + qq("0 deg") + "\n      plunge: " + qq("0 deg") +
                          "\n      half_length: " + qq("4000 m") + "\n      half_width_left: " + qq("2500 m") + "\n      half_width_right: " +
                          qq("2500 m") + "\n      relief:\n        H1: " + qq(fmt_fixed(relief, 3) + " m") + "\n        H2: " +
                          qq(fmt_fixed(relief, 3) + " m") + "\n";
  std::string faults = planar_fault("F1", "1", "[6000, 5000] m", "1200 m", "0 deg", "60 deg", "20000 m",
                                    "max: " + q(fmt_fixed(throw_m, 3) + " m"), "dip_slip");
  std::string inputs = conditioning.empty() ? "" : "inputs:\n  - {id: D, path: data.csv, role: surface points, sharing: shareable, source: test, evidence: interpreted}\n";
  return model_yaml(faults, structure, {1000, 1200, 1400}, inputs + conditioning);
}
/// Sample H1 of the true model at random positions, with N(0, 3 m) noise, sigma 3 m.
void write_hybrid_data(const fs::path& dir) {
  auto truth = build_text(hybrid_yaml(250, 120, ""));
  std::mt19937_64 rng(77);
  std::uniform_real_distribution<double> U(1500, 8500);
  std::normal_distribution<double> N(0, 3.0);
  std::string s = "X,Y,Z,surface,sigma\n";
  int n = 0;
  while (n < 240) {
    double x = U(rng), y = U(rng);
    const int h = 0;
    double z = -depth_faulted(*truth, h, x, y);
    if (std::isnan(z)) continue;
    s += fmt_fixed(x, 2) + "," + fmt_fixed(y, 2) + "," + fmt_fixed(z + N(rng), 3) + "," + (h ? "H2" : "H1") + ",3.0\n";
    ++n;
  }
  write_text(dir / "data.csv", s);
}
std::string conditioning_block(const std::string& throw_range, const std::string& relief_range, const std::string& iters = "50") {
  return "conditioning:\n  data:\n    - {input: D, kind: points, format: csv_points}\n  use_well_tops: false\n  uncertainty: " +
         qq("3 m") + "\n  parameters:\n    - {path: faults.F1.throw.max, " + throw_range + "}\n    - {path: structure.primitives.A1.relief.H1, " +
         relief_range + "}\n  max_iterations: " + qq(iters + " ordinal") + "\n";
}
}  // namespace

TEST_CASE("V7 hybrid conditioning recovers a fault throw and a fold relief within their uncertainty") {
  auto dir = tmpdir("hybrid");
  write_hybrid_data(dir);
  std::string y = hybrid_yaml(175, 84, conditioning_block("min: " + qq("20 m") + ", max: " + qq("300 m"),
                                                          "min: " + qq("50 m") + ", max: " + qq("500 m")));
  desc::Description d;
  auto rep = desc::validate_description_text(y, dir, &d);
  if (!rep.ok()) std::cerr << rep.text();
  REQUIRE(rep.ok());
  BuildInputs in;
  io::load_build_inputs(d, in);
  BuildLog log;
  auto r = condition_to_data(d, in, log);
  std::printf("  converged %d after %d iterations, reduced chi2 %.3f (%s)\n", r.converged, r.iterations, r.chi2_reduced, r.reason.c_str());
  REQUIRE(r.converged);
  REQUIRE(r.parameters.size() == 2);
  const double truth[] = {120, 250};
  for (std::size_t k = 0; k < 2; ++k) {
    auto& p = r.parameters[k];
    std::printf("  %s: %.3f +- %.3f %s (start %.1f, truth %.1f)\n", p.path.c_str(), p.value, p.sigma, p.unit.c_str(), p.start, truth[k]);
    CHECK(std::fabs(p.value - truth[k]) <= 2 * p.sigma);
    CHECK(p.sigma > 0);
    CHECK(!p.controlling.empty());
  }
  CHECK(r.chi2_reduced >= 0.5);
  CHECK(r.chi2_reduced <= 2.0);
  CHECK(r.leverage.size() == r.residuals.size());
  // the fitted values are labelled calibrated, with the data set as source
  auto& fp = r.fitted.faults[0].throw_max;
  CHECK(fp.evidence == desc::Evidence::Calibrated);
  CHECK(fp.source.find("D") != std::string::npos);
  CHECK(fp.source.find("Levenberg") != std::string::npos);
  // and the conditioned description builds, with C22 passing
  BuildInputs in2;
  io::load_build_inputs(d, in2);
  BuildLog l2;
  auto m = build_model(d, in2, l2);
  REQUIRE(m->capability("conditioning") != nullptr);
  auto rep2 = checks_for(*m, d);
  auto* c22 = find_check(rep2, "C22");
  REQUIRE(c22 != nullptr);
  CHECK(c22->passed);
  auto* c21 = find_check(rep2, "C21");
  REQUIRE(c21 != nullptr);
  CHECK(c21->passed);
}

TEST_CASE("V8 hybrid conditioning refuses fits that do not converge inside their range or are not controlled by the data") {
  auto dir = tmpdir("hybrid_refuse");
  write_hybrid_data(dir);
  auto run = [&](const std::string& y) {
    desc::Description d;
    auto rep = desc::validate_description_text(y, dir, &d);
    if (!rep.ok()) std::cerr << rep.text();
    REQUIRE(rep.ok());
    BuildInputs in;
    io::load_build_inputs(d, in);
    BuildLog log;
    return condition_to_data(d, in, log);
  };
  // the range excludes the true throw (120 m): converges on the bound
  auto r1 = run(hybrid_yaml(175, 84, conditioning_block("min: " + qq("20 m") + ", max: " + qq("100 m"), "min: " + qq("50 m") + ", max: " + qq("500 m"))));
  std::printf("  range excluding the truth: %s\n", r1.reason.c_str());
  CHECK(!r1.converged);
  CHECK(r1.reason.find("range") != std::string::npos);
  // one iteration only
  auto r2 = run(hybrid_yaml(175, 84, conditioning_block("min: " + qq("20 m") + ", max: " + qq("300 m"), "min: " + qq("50 m") + ", max: " + qq("500 m"), "1")));
  std::printf("  max_iterations 1: %s\n", r2.reason.c_str());
  CHECK(!r2.converged);
  CHECK(r2.reason.find("converge") != std::string::npos);
  // a parameter the data do not control: the relief of H3, which has no data
  std::string y3 = hybrid_yaml(175, 84, conditioning_block("min: " + qq("20 m") + ", max: " + qq("300 m"), "min: " + qq("50 m") + ", max: " + qq("500 m")));
  auto p = y3.find("structure.primitives.A1.relief.H1");
  REQUIRE(p != std::string::npos);
  y3.replace(p, 33, "structure.primitives.A1.relief.H3");
  auto p2 = y3.find("        H2: ");
  REQUIRE(p2 != std::string::npos);
  y3.insert(p2, "        H3: " + qq("100 m") + "\n");
  auto r3 = run(y3);
  std::printf("  uncontrolled parameter: %s\n", r3.reason.c_str());
  CHECK(!r3.converged);
  CHECK(r3.reason.find("relief.H3") != std::string::npos);
  // and a refused fit builds no model
  desc::Description d1;
  REQUIRE(desc::validate_description_text(hybrid_yaml(175, 84, conditioning_block("min: " + qq("20 m") + ", max: " + qq("100 m"),
                                                                                 "min: " + qq("50 m") + ", max: " + qq("500 m"))),
                                          dir, &d1)
              .ok());
  BuildInputs in1;
  io::load_build_inputs(d1, in1);
  BuildLog l1;
  bool threw = false;
  try {
    build_model(d1, in1, l1);
  } catch (const std::exception& e) {
    threw = std::string(e.what()).find("range") != std::string::npos;
  }
  CHECK(threw);
}

TEST_CASE("V9 an implicit model with series and data-derived faults builds and passes its checks") {
  register_all_plugins();
  auto f = fs::path(STX_SOURCE_DIR) / "examples/implicit_faults_series_from_data.yaml";
  REQUIRE(fs::exists(f));
  auto d = desc::load_description(f);
  BuildInputs in;
  io::load_build_inputs(d, in);
  BuildLog log;
  auto m = build_model(d, in, log);
  CHECK(m->info().engine == "implicit-rbf");
  CHECK(run_model_suite(*m).passed());
  auto rep = checks_for(*m, d);
  for (auto& c : rep.results)
    if (!c.passed && !c.warn_only) std::printf("  FAIL %s %s: %s\n", c.id.c_str(), c.subject.c_str(), c.detail.c_str());
  CHECK(rep.passed());
  // the estimated throws are labelled calibrated
  bool estimated = false;
  for (auto& [k, v] : log.values) estimated |= k.find("estimated_throw") != std::string::npos;
  CHECK(estimated);
}

namespace {
/// Depth differences between an implicit model and its parametric source over a 100 m lattice,
/// inside the data footprint (250 m from the edges): away from faults (more than 500 m from a gap
/// of the source) and near them. Returns {p95 away, p95 near}.
std::pair<double, double> rebuild_errors(const Model& src, const Model& imp) {
  GridSpec g = GridSpec::covering(src.info().domain, 100);
  std::vector<double> away, near;
  const Domain& dm = src.info().domain;
  for (std::size_t h = 0; h < src.info().horizons.size(); ++h) {
    auto a = repr::horizon_grid_faulted(src, static_cast<int>(h), g), b = repr::horizon_grid_faulted(imp, static_cast<int>(h), g);
    std::vector<int> dist(g.size(), 1 << 20);  // chamfer distance (nodes) to a gap of the source
    for (std::size_t k = 0; k < g.size(); ++k)
      if (std::isnan(a.z.v[k])) dist[k] = 0;
    for (int pass = 0; pass < 2; ++pass)
      for (std::size_t jj = 0; jj < g.ny; ++jj)
        for (std::size_t ii = 0; ii < g.nx; ++ii) {
          std::size_t i = pass ? g.nx - 1 - ii : ii, j = pass ? g.ny - 1 - jj : jj;
          int& d = dist[g.index(i, j)];
          long ni = static_cast<long>(i) + (pass ? 1 : -1), nj = static_cast<long>(j) + (pass ? 1 : -1);
          if (ni >= 0 && ni < static_cast<long>(g.nx)) d = std::min(d, dist[g.index(static_cast<std::size_t>(ni), j)] + 1);
          if (nj >= 0 && nj < static_cast<long>(g.ny)) d = std::min(d, dist[g.index(i, static_cast<std::size_t>(nj))] + 1);
        }
    for (std::size_t k = 0; k < g.size(); ++k) {
      double x = g.x(k % g.nx), y = g.y(k / g.nx);
      if (x < dm.x0 + 250 || x > dm.x0 + dm.sx - 250 || y < dm.y0 + 250 || y > dm.y0 + dm.sy - 250) continue;
      if (std::isnan(a.z.v[k]) || std::isnan(b.z.v[k])) continue;
      (dist[k] > 5 ? away : near).push_back(std::fabs(a.z.v[k] - b.z.v[k]));
    }
  }
  return {percentile95(away), near.empty() ? 0.0 : percentile95(near)};
}
}  // namespace

TEST_CASE("V10 synthetic v2 examples sampled into data are rebuilt within the stated tolerances; throws within max(15 %, 5 m)") {
  register_all_plugins();
  struct Case { const char* implicit; const char* source; std::vector<std::pair<std::string, double>> throws; };
  const std::vector<Case> cases = {
      {"implicit_unconformity_from_samples", "unconformity_onlap_pinchout", {}},
      {"implicit_relay_splay_from_samples", "fault_system_relay_splay", {{"F1", 120.0}, {"F3", 120.0}}},
      {"implicit_detachment_folds_from_samples", "detachment_folds", {}},
      {"implicit_faults_series_from_data", "data/implicit_faults_series/truth", {{"F1", 80.0}}},
  };
  for (auto& c : cases) {
    auto ds = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples" / (std::string(c.source) + ".yaml"));
    auto di = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples" / (std::string(c.implicit) + ".yaml"));
    BuildInputs ins, ini;
    io::load_build_inputs(ds, ins);
    io::load_build_inputs(di, ini);
    BuildLog ls, li;
    auto ms = build_model(ds, ins, ls);
    auto mi = build_model(di, ini, li);
    auto [away, near] = rebuild_errors(*ms, *mi);
    std::printf("  %s: p95 away from faults %.3f m, near faults %.3f m\n", c.implicit, away, near);
    CHECK(away <= 10.0);
    CHECK(near <= 30.0);
    for (auto& [f, truth] : c.throws) {
      double t = log_value(li, "implicit.fault." + f + ".estimated_throw");
      std::printf("    %s estimated throw %.2f m (source %.1f m)\n", f.c_str(), t, truth);
      CHECK(std::fabs(t - truth) <= std::max(0.15 * truth, 5.0));
    }
  }
}

STX_TEST_MAIN()
