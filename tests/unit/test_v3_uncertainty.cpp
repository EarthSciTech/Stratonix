// Stage 3C analytical tests (STX-TN-13 section 8), written before the features (Prompt 40) and
// made to pass by Prompt 41: the generator, distributions and correlations, the linear key number,
// sensitivity ranking, scenarios, failures reported, probability maps, byte-identical reruns and
// rule F33.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>

#include "v2_fixture.hpp"
#include "stratonix/package/uncertainty.hpp"

namespace fs = std::filesystem;
using namespace stratonix;
using namespace v2fx;

namespace {

fs::path tmpdir(const std::string& name) {
  fs::path p = fs::path(stxtest::tmp_dir()) / ("v3u_" + name);
  fs::remove_all(p);
  fs::create_directories(p);
  return p;
}
fs::path write_doc(const std::string& name, const std::string& text) {
  auto dir = tmpdir(name);
  fs::path f = dir / "model.yaml";
  std::ofstream(f) << text;
  return f;
}
std::string read(const fs::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::stringstream s;
  s << f.rdbuf();
  return s.str();
}
const unc::KeyStats* key_stats(const unc::EnsembleResult& r, const std::string& key) {
  for (auto& k : r.stats)
    if (k.key == key) return &k;
  std::cerr << "  no statistics for key '" << key << "'\n";
  return nullptr;
}

/// One planar normal fault F1 at x = 5 km (constant throw `throw_max`), optionally F2 at x = 8 km.
std::string faults(const std::string& t1, const std::string& t2 = "") {
  std::string f = planar_fault("F1", "1", "[3000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", "max: " + q(t1), "dip_slip");
  if (!t2.empty()) f += planar_fault("F2", "2", "[7500, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", "max: " + q(t2), "dip_slip");
  return f;
}
std::string block(int members, const std::string& body) {
  // the test faults offset the top and base horizons on purpose (C6 reported, not failed)
  return "checks:\n  warn_only: [C6]\nuncertainty:\n  members: " + q(std::to_string(members) + " ordinal") + "\n  seed: " + q("20261115 ordinal") + "\n" + body;
}
std::string normal(const std::string& path, double mean, double sd, const std::string& unit = "m") {
  return "    - {path: " + path + ", distribution: normal, mean: " + q(fmt_fixed(mean, 3) + " " + unit) + ", sd: " + q(fmt_fixed(sd, 3) + " " + unit) + "}\n";
}
std::string uniform(const std::string& path, double lo, double hi, const std::string& unit = "m") {
  return "    - {path: " + path + ", distribution: uniform, min: " + q(fmt_fixed(lo, 3) + " " + unit) + ", max: " + q(fmt_fixed(hi, 3) + " " + unit) + "}\n";
}
/// Ensembles of hundreds of builds run in optimised builds only; debug and sanitiser builds skip
/// them (stating so) and keep the generator, distribution, rerun and description cases.
#if defined(NDEBUG) && !defined(__SANITIZE_ADDRESS__)
#define STX_FULL_ENSEMBLES 1
#else
#define STX_FULL_ENSEMBLES 0
#endif
bool full_ensembles(const char* name) {
  if (STX_FULL_ENSEMBLES) return true;
  std::printf("  skipped in this build (debug or sanitiser): %s runs in the release builds\n", name);
  return false;
}
unc::EnsembleOptions quick(unsigned threads = 0) {
  unc::EnsembleOptions o;
  o.figures = false;
  o.threads = threads;
  return o;
}

}  // namespace

TEST_CASE("U1 generator: reproducible streams, member seeds independent of order") {
  unc::Rng a(unc::member_seed(42, 7)), b(unc::member_seed(42, 7)), c(unc::member_seed(42, 8));
  bool same = true, differ = false;
  for (int i = 0; i < 1000; ++i) {
    auto x = a.next(), y = b.next(), z = c.next();
    same = same && x == y;
    differ = differ || x != z;
  }
  CHECK(same);
  CHECK(differ);
  // golden values: the first draws of SplitMix64-seeded xoshiro256** for seed 0, member 0
  unc::Rng g(unc::member_seed(0, 0));
  const std::uint64_t first = g.next();
  std::printf("  first draw for (0, 0): %llu\n", static_cast<unsigned long long>(first));
  CHECK(first == 7421629122807502682ull);  // golden value, the same on every platform
  unc::Rng g2(unc::member_seed(0, 0));
  CHECK(g2.next() == first);
  double lo = 1, hi = 0;
  for (int i = 0; i < 100000; ++i) {
    double u = a.uniform();
    lo = std::min(lo, u);
    hi = std::max(hi, u);
  }
  CHECK(lo > 0 && hi < 1);
}

TEST_CASE("U2 distributions: inverse CDFs at known quantiles, sample moments within 4 standard errors") {
  CHECK_NEAR(unc::normal_quantile(0.975), 1.959963984540054, 1e-9);
  CHECK_NEAR(unc::normal_quantile(0.5), 0.0, 1e-12);
  CHECK_NEAR(unc::normal_cdf(1.959963984540054), 0.975, 1e-9);
  unc::Distribution uni{"uniform", 2, 5}, nor{"normal", 100, 10}, tri{"triangular", 1, 2, 4}, dis{"discrete"};
  dis.values = {1, 2, 3};
  dis.weights = {0.2, 0.5, 0.3};
  CHECK_NEAR(uni.quantile(0.25), 2.75, 1e-12);
  CHECK_NEAR(nor.quantile(0.975), 100 + 10 * 1.959963984540054, 1e-7);
  CHECK_NEAR(tri.quantile(1.0 / 3.0), 2.0, 1e-9);  // F(mode) = (mode - a) / (c - a) = 1/3
  CHECK_NEAR(dis.quantile(0.6), 2.0, 0);
  CHECK_NEAR(tri.mean(), 7.0 / 3.0, 1e-12);
  CHECK_NEAR(tri.variance(), (1.0 + 4 + 16 - 2 - 4 - 8) / 18.0, 1e-12);
  unc::Rng r(unc::member_seed(1, 1));
  const int n = 20000;
  for (auto* d : {&uni, &nor, &tri, &dis}) {
    double s = 0, s2 = 0;
    for (int i = 0; i < n; ++i) {
      double x = d->quantile(r.uniform());
      s += x;
      s2 += x * x;
    }
    double m = s / n, v = s2 / n - m * m, se = std::sqrt(d->variance() / n);
    std::printf("  %s: mean %.4f (%.4f), variance %.4f (%.4f)\n", d->type.c_str(), m, d->mean(), v, d->variance());
    CHECK(std::fabs(m - d->mean()) <= 4 * se);
    CHECK(std::fabs(v - d->variance()) <= 0.05 * d->variance());
  }
}

TEST_CASE("U3 correlation through a Gaussian copula; a non-positive-definite matrix is refused") {
  unc::EnsembleSpec s;
  s.parameters.resize(2);
  for (auto& p : s.parameters) p.dist = {"normal", 0, 1};
  s.correlations = {{{0, 1}, 0.7}};
  unc::Rng r(unc::member_seed(3, 0));
  const int n = 5000;
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (int i = 0; i < n; ++i) {
    auto u = unc::correlated_uniforms(s, r);
    double x = unc::normal_quantile(u[0]), y = unc::normal_quantile(u[1]);
    sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
  }
  double cov = sxy / n - sx / n * sy / n, rho = cov / std::sqrt((sxx / n - sx * sx / n / n) * (syy / n - sy * sy / n / n));
  std::printf("  sample correlation %.4f\n", rho);
  CHECK_NEAR(rho, 0.7, 0.03);
  CHECK_THROWS(unc::cholesky({1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1}, 3));
}

TEST_CASE("U4 a key number linear in one normal parameter: mean and standard deviation of the ensemble") {
  if (!full_ensembles("U4")) return;
  auto f = write_doc("linear", model_yaml(faults("100 m"), "", {1000, 1200, 1400, 1600},
                                          block(200, "  parameters:\n" + normal("faults.F1.throw.max", 100, 10))));
  auto r = unc::run_ensemble(f, f.parent_path() / "out", quick());
  auto* k = key_stats(r, "Max throw F1");
  REQUIRE(k != nullptr);
  std::printf("  Max throw F1: n %zu, mean %.3f, sd %.3f, P10 %.2f, P50 %.2f, P90 %.2f\n", k->stats.n, k->stats.mean, k->stats.sd, k->stats.p10,
              k->stats.p50, k->stats.p90);
  CHECK(k->stats.n == 200);
  CHECK_NEAR(k->stats.mean, 100, 2.2);  // 3 standard errors (10 / sqrt(200))
  CHECK_NEAR(k->stats.sd, 10, 1.5);
  std::string st = read(f.parent_path() / "out" / "ensemble" / "statistics.csv");
  CHECK(st.find("20261115") != std::string::npos);  // the seed travels with every range
}

TEST_CASE("U5 sensitivity ranks the controlling parameter (first-order index and one at a time)") {
  if (!full_ensembles("U5")) return;
  auto f = write_doc("sens", model_yaml(faults("100 m", "41 m"), "", {1000, 1200, 1400, 1600},
                                        block(200, "  one_at_a_time: true\n  parameters:\n" + uniform("faults.F1.throw.max", 80, 120) +
                                                       uniform("faults.F2.throw.max", 40, 42))));
  auto r = unc::run_ensemble(f, f.parent_path() / "out", quick());
  double s1 = -1, s2 = -1;
  bool d2 = true;
  for (auto& s : r.sobol)
    if (s.key == "Max throw F1") {
      if (s.parameter == "faults.F1.throw.max") s1 = s.index;
      if (s.parameter == "faults.F2.throw.max") { s2 = s.index; d2 = s.distinguishable; }
    }
  std::printf("  first-order indices for Max throw F1: F1 %.3f, F2 %.3f\n", s1, s2);
  CHECK(s1 > 0.9);
  CHECK(!d2);
  double swing1 = 0, swing2 = 0;
  for (auto& o : r.oat)
    if (o.key == "Max throw F1") (o.parameter == "faults.F1.throw.max" ? swing1 : swing2) = std::fabs(o.high - o.low);
  std::printf("  one-at-a-time swings: F1 %.3f m, F2 %.3f m\n", swing1, swing2);
  CHECK_NEAR(swing1, 32.0, 0.5);  // P90 - P10 of U(80, 120)
  CHECK(swing1 > swing2);
}

TEST_CASE("U6 scenarios: option frequencies follow their weights; each member's option is recorded") {
  if (!full_ensembles("U6")) return;
  std::string sc = "  scenarios:\n    - name: slip\n      options:\n"
                   "        - {name: dip_slip, weight: " + q("0.3 fraction") + ", set: {faults.F1.slip_mode: dip_slip}}\n"
                   "        - {name: vertical, weight: " + q("0.7 fraction") + ", set: {faults.F1.slip_mode: vertical_offset}}\n";
  auto f = write_doc("scen", model_yaml(faults("100 m"), "", {1000, 1200, 1400, 1600}, block(200, sc)));
  auto r = unc::run_ensemble(f, f.parent_path() / "out", quick());
  std::size_t first = 0, n = 0;
  for (auto& m : r.members) {
    if (m.sample.index == 0) continue;
    REQUIRE(m.sample.options.size() == 1);
    first += m.sample.options[0] == 0;
    ++n;
  }
  const double p = static_cast<double>(first) / static_cast<double>(n), se = std::sqrt(0.3 * 0.7 / static_cast<double>(n));
  std::printf("  option 1 chosen in %.3f of %zu members (weight 0.3)\n", p, n);
  CHECK(std::fabs(p - 0.3) <= 3 * se);
  CHECK(read(f.parent_path() / "out" / "ensemble" / "members.csv").find("vertical") != std::string::npos);
}

TEST_CASE("U7 failed members are reported with their parameters, never dropped") {
  if (!full_ensembles("U7")) return;
  auto f = write_doc("fail", model_yaml(faults("100 m"), "", {1000, 1200, 1400, 1600},
                                        block(100, "  parameters:\n" + uniform("wells.W1.tops.H2", 950, 1250))));  // H2 above H1 in a sixth of them
  auto r = unc::run_ensemble(f, f.parent_path() / "out", quick());
  std::size_t failed = 0, with_reason = 0;
  for (auto& m : r.members) {
    if (m.sample.index == 0 || m.status == "ok") continue;
    ++failed;
    with_reason += !m.reason.empty() && m.sample.values.size() == 1;
  }
  std::printf("  %zu ok, %zu failed\n", r.ok(), failed);
  CHECK(failed > 0);
  CHECK(with_reason == failed);
  CHECK(r.ok() + failed == 100);
  auto* k = key_stats(r, "Max throw F1");
  REQUIRE(k != nullptr);
  CHECK(k->failed == failed);
}

TEST_CASE("U8 probability map: a unit at a depth where the horizon is uniform about it") {
  if (!full_ensembles("U8")) return;
  std::string pm = "  probability_maps:\n    - {kind: unit, unit: U1, depth: " + q("1200 m") + ", spacing: " + q("1000 m") + "}\n";
  auto f = write_doc("pmap", model_yaml("", "", {1000, 1200, 1400, 1600},
                                        block(200, "  parameters:\n" + uniform("wells.W1.tops.H2", 1150, 1250) + pm)));
  auto r = unc::run_ensemble(f, f.parent_path() / "out", quick());
  REQUIRE(r.maps.size() == 1);
  double lo = 1, hi = 0;
  for (double p : r.maps[0].p) { lo = std::min(lo, p); hi = std::max(hi, p); }
  std::printf("  P(U1 at 1200 m): %.3f to %.3f over %zu nodes\n", lo, hi, r.maps[0].p.size());
  CHECK(lo >= 0.44 && hi <= 0.56);
  CHECK(fs::exists(f.parent_path() / "out" / "ensemble" / "maps"));
}

TEST_CASE("U9 the same seed gives byte-identical results with 1 and 4 threads") {
  std::string text = model_yaml(faults("100 m"), "", {1000, 1200, 1400, 1600},
                                block(40, "  parameters:\n" + normal("faults.F1.throw.max", 100, 10) + uniform("wells.W1.tops.H2", 1180, 1220)));
  auto f = write_doc("rerun", text);
  unc::run_ensemble(f, f.parent_path() / "a", quick(1));
  unc::run_ensemble(f, f.parent_path() / "b", quick(4));
  for (const char* t : {"members.csv", "statistics.csv"}) {
    auto a = read(f.parent_path() / "a" / "ensemble" / t), b = read(f.parent_path() / "b" / "ensemble" / t);
    CHECK(!a.empty());
    CHECK(a == b);
  }
}

TEST_CASE("U10 description: an uncertainty block validates; rule F33 refuses bad distributions") {
  auto ok = validate(model_yaml(faults("100 m"), "", {1000, 1200, 1400, 1600}, block(10, "  parameters:\n" + normal("faults.F1.throw.max", 100, 10))));
  if (!ok.ok()) std::cerr << ok.text();
  CHECK(ok.ok());
  const std::vector<std::string> bad = {
      "  parameters:\n" + normal("faults.F9.throw.max", 100, 10),                       // unknown path
      "  parameters:\n" + normal("faults.F1.throw.max", 100, 0),                        // sd <= 0
      "  parameters:\n" + normal("faults.F1.throw.max", 100, 10) + normal("wells.W1.tops.H2", 1200, 5) +
          "  correlations: [{a: faults.F1.throw.max, b: wells.W1.tops.H2, rho: " + q("1.2 ratio") + "}]\n",  // |rho| >= 1
      "  scenarios:\n    - name: s\n      options:\n        - {name: a, weight: " + q("0.5 fraction") + ", set: {faults.F1.slip_mode: dip_slip}}\n"
      "        - {name: b, weight: " + q("0.4 fraction") + ", set: {faults.F1.slip_mode: vertical_offset}}\n",  // weights sum 0.9
  };
  for (auto& b : bad) {
    auto rep = validate(model_yaml(faults("100 m"), "", {1000, 1200, 1400, 1600}, block(10, b)));
    CHECK(!rep.ok());
    CHECK(rep.text().find("F33") != std::string::npos);
  }
}

STX_TEST_MAIN()
