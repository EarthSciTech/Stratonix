// Stage 2C (Prompts 23 and 24): stratigraphic series and unconformities, onlap, thickness trends
// and pinch-outs, wells honouring the series, isopachs, growth strata and clinoforms (STX-TN-07).
#include "v2_fixture.hpp"
#include "stratonix/engines/structure.hpp"

using namespace v2fx;

namespace {

/// Description with the given horizon list, reference-well tops, units, series and structure.
std::string doc(const std::string& horizons, const std::string& tops, const std::string& units, const std::string& series,
                const std::string& structure = "structure:\n  regional: {type: flat}\n", const std::string& faults = "",
                const std::string& wells = "") {
  std::string s = "stratonix_schema: 1\nmodel: {name: t, version: \"1.0.0\", status: DESIGN}\n"
                  "crs: {name: local, vertical_datum: datum, length_unit: m}\n"
                  "domain:\n  origin: " + q("[0, 0] m") + "\n  size: " + q("[10000, 10000] m") +
                  "\n  depth_range: " + q("[0, 3000] m") + "\n  grid_spacing: " + q("50 m") + "\n";
  s += "wells:\n  - {name: W1, role: control, head: " + q("[1000, 1000] m") + ", tops: {" + tops + "}}\n" + wells;
  s += "stratigraphy:\n  reference: {well: W1}\n  horizons: [" + horizons + "]\n  units: [" + units + "]\n";
  if (!series.empty()) s += "  series:\n" + series;
  s += structure;
  if (!faults.empty()) s += "faults:\n" + faults + "frame: {fixed: {well: W1}}\n";
  return s;
}
std::string top(const std::string& h, double d) { return h + ": " + q(fmt_fixed(d, 1) + " m"); }
std::string units_for(const std::vector<std::string>& hs) {
  std::string u;
  for (std::size_t k = 0; k + 1 < hs.size(); ++k)
    u += (k ? ", " : "") + std::string("{name: Unit_") + std::to_string(k + 1) + ", top: " + hs[k] + ", base: " + hs[k + 1] + "}";
  return u;
}
const Structure& structure_of(const Model& m) { return *static_cast<const Structure*>(m.capability("structure")); }

// Angular unconformity: S1 (H3-H5) tilted 5 degrees west-down, truncated by the flat erosional H2.
const double kTan5 = std::tan(deg2rad(5.0));
double tilted(double d0, double x) { return d0 - kTan5 * (x - 1000); }  // depth, shallower eastwards
std::string angular() {
  return doc("{name: H1}, {name: H2}, {name: H3}, {name: H4}, {name: H5}",
             top("H1", 900) + ", " + top("H2", 1100) + ", " + top("H3", 1150) + ", " + top("H4", 1300) + ", " + top("H5", 1450),
             units_for({"H1", "H2", "H3", "H4", "H5"}),
             "    - {name: S2, horizons: [H1, H2], base: erosional}\n"
             "    - {name: S1, horizons: [H3, H4, H5], base: conformable, tilt: {dip: " + q("5 deg") + ", azimuth: " + q("270 deg") + "}}\n",
             "structure:\n  regional: {type: flat}\n", "",
             "  - {name: W2, role: planned, head: " + q("[6000, 5000] m") + ", total_depth: " + q("2500 m") + "}\n"
             "  - {name: W3, role: planned, head: " + q("[1200, 5000] m") + ", total_depth: " + q("2500 m") + "}\n");
}

// Onlap: S2 (H1, H2) laps onto H3, folded by an anticline restricted to S1.
double bump(double t, double w) { return std::fabs(t) < w ? 0.5 * (1 + std::cos(kPi * std::fabs(t) / w)) : 0.0; }
std::string onlap() {
  return doc("{name: H1}, {name: H2}, {name: H3}, {name: H4}",
             top("H1", 900) + ", " + top("H2", 1000) + ", " + top("H3", 1100) + ", " + top("H4", 1300), units_for({"H1", "H2", "H3", "H4"}),
             "    - {name: S2, horizons: [H1, H2], base: onlap}\n    - {name: S1, horizons: [H3, H4], base: conformable}\n",
             "structure:\n  regional: {type: flat}\n  primitives:\n    - {type: anticline, name: A, series: S1, crest: " + q("[5000, 5000] m") +
                 ", axis_azimuth: " + q("0 deg") + ", half_width_left: " + q("2000 m") + ", half_width_right: " + q("2000 m") +
                 ", relief: {H3: " + q("300 m") + "}}\n");
}
double onlap_h3(double x) { return 1100 - 300 * bump(x - 5000, 2000); }

}  // namespace

TEST_CASE("P23 one conformable series equals a v1 description") {
  std::string a = model_yaml("");
  std::string b = a;
  b.replace(b.find("  units: ["), 0, "  series:\n    - {name: S1, horizons: [H1, H2, H3, H4], base: conformable}\n");
  auto ma = build_text(a), mb = build_text(b);
  double worst = 0;
  for (int h = 0; h < 4; ++h)
    for (double x = 500; x < 10000; x += 1500) worst = std::max(worst, std::fabs(depth_faulted(*ma, h, x, 4000) - depth_faulted(*mb, h, x, 4000)));
  CHECK(worst == 0.0);
}

TEST_CASE("P23 angular unconformity: truncation and missing section") {
  desc::Description d;
  auto m = build_text(angular(), &d);
  const double d0[] = {1150, 1300, 1450};
  double worst = 0;
  for (double x = 200; x <= 9800; x += 200)
    for (int k = 0; k < 3; ++k) worst = std::max(worst, std::fabs(depth_faulted(*m, 2 + k, x, 5000) - std::max(tilted(d0[k], x), 1100.0)));
  std::cout << "  worst truncated-horizon misfit " << fmt_fixed(worst, 4) << " m\n";
  CHECK(worst <= 0.01);
  CHECK_NEAR(depth_faulted(*m, 1, 7000, 5000), 1100, 0.01);
  // missing section of Unit_3 (H3-H4): raw 150 m minus what is left below the unconformity
  const Structure& st = structure_of(*m);
  for (double x : {1000.0, 2500.0, 4000.0, 8000.0}) {
    double left = std::max(tilted(1300, x), 1100.0) - std::max(tilted(1150, x), 1100.0);
    CHECK_NEAR(st.missing_section(3, x, 5000), 150 - left, 0.01);
  }
  auto rep = checks_for(*m, d);
  auto c13 = find_check(rep, "C13"), c14 = find_check(rep, "C14");
  REQUIRE(c13 != nullptr);
  REQUIRE(c14 != nullptr);
  CHECK(c13->passed);
  CHECK(c14->passed);
  for (auto& r : rep.results)
    if (!r.passed && !r.warn_only) std::cout << "  failed: " << r.id << " " << r.subject << " " << r.value << "\n";
  CHECK(rep.passed());
}

TEST_CASE("P23 wells honour the series: an eroded horizon has no top") {
  desc::Description d;
  auto m = build_text(angular(), &d);
  auto paths = make_well_paths(d);
  for (auto& w : paths) {
    if (w.name() == "W1") continue;
    auto tops = repr::well_intersections(*m, w);
    std::string names;
    for (auto& t : tops) names += t.horizon + " ";
    std::cout << "  " << w.name() << ": " << names << "\n";
    if (w.name() == "W2") CHECK(names == "H1 H2 ");                // S1 eroded entirely at x = 6000
    if (w.name() == "W3") CHECK(names == "H1 H2 H3 H4 H5 ");       // nothing eroded at x = 1200
  }
}

TEST_CASE("P23 onlap onto a fold of the older series") {
  auto m = build_text(onlap());
  double worst = 0;
  for (double x = 500; x <= 9500; x += 100) {
    double h3 = onlap_h3(x);
    worst = std::max(worst, std::fabs(depth_faulted(*m, 2, x, 5000) - h3));
    worst = std::max(worst, std::fabs(depth_faulted(*m, 1, x, 5000) - std::min(1000.0, h3)));
    worst = std::max(worst, std::fabs(depth_faulted(*m, 0, x, 5000) - std::min(900.0, h3)));
    worst = std::max(worst, std::fabs(depth_faulted(*m, 3, x, 5000) - (h3 + 200)));  // H4: relief constant below the anchor
  }
  std::cout << "  worst onlap misfit " << fmt_fixed(worst, 4) << " m\n";
  CHECK(worst <= 0.01);
  CHECK_NEAR(depth_faulted(*m, 1, 5000, 5000) - depth_faulted(*m, 2, 5000, 5000), 0, 0.01);  // Unit_2 absent on the crest
  // isopach of Unit_2: zero on the crest, 100 m far from it
  auto t = repr::unit_thickness_grid(*m, 2, GridSpec{500, 5000, 100, 1, 90, 1});
  CHECK_NEAR(t.at(0, 0), 100, 0.01);
  CHECK_NEAR(t.at(45, 0), 0, 0.01);
  auto lines = repr::termination_lines(*m, 2, GridSpec::covering(m->info().domain, 50));
  REQUIRE(!lines.empty());
  CHECK(lines[0].kind == "onlap");
}

TEST_CASE("P23 pinch-out of a linear thickness trend") {
  std::string hs = "{name: H1}, {name: H2, isochore: {from: H1, thickness: " + q("200 m") + ", trend: {type: linear, gradient: " +
                   q("[-0.05, 0] ratio") + "}}}, {name: H3}";
  auto m = build_text(doc(hs, top("H1", 900) + ", " + top("H2", 1100) + ", " + top("H3", 1400), units_for({"H1", "H2", "H3"}), ""));
  double worst = 0;
  for (double x = 100; x <= 9900; x += 100) {
    double t = std::max(0.0, 200 - 0.05 * (x - 1000));
    worst = std::max(worst, std::fabs(depth_faulted(*m, 1, x, 3000) - depth_faulted(*m, 0, x, 3000) - t));
  }
  std::cout << "  worst thickness misfit " << fmt_fixed(worst, 4) << " m\n";
  CHECK(worst <= 0.01);
  auto lines = repr::termination_lines(*m, 1, GridSpec::covering(m->info().domain, 50));
  REQUIRE(!lines.empty());
  CHECK(lines[0].kind == "pinch-out");
  double xs = 0;
  int n = 0;
  for (auto& p : lines[0].line) { xs += p.x; ++n; }
  CHECK_NEAR(xs / n, 5000, 50);
}

TEST_CASE("P24 growth fault: throw decreases upwards through the growth package") {
  BuildLog log;
  auto m = build_text(model_yaml(planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", "max: " + q("200 m"),
                                              "dip_slip", "{type: layer_rigid}", "    growth: {onset: H3, end: H1}\n")),
                      nullptr, &log);
  const double d0[] = {1000, 1200, 1400, 1600}, off[] = {0, 100, 200, 200};
  for (int h = 0; h < 4; ++h) CHECK_NEAR(depth_faulted(*m, h, 8000, 5000) - d0[h], off[h], 0.5);
  for (int h = 0; h < 4; ++h) CHECK_NEAR(depth_faulted(*m, h, 2000, 5000), d0[h], 0.01);
  // sequential unwinding: removing f(H2) T = 100 m flattens H2 and leaves H3 offset by 100 m
  CHECK_NEAR((depth_faulted(*m, 2, 8000, 5000) - 100) - d0[2], 100, 0.5);
  CHECK_NEAR(log_value(log, "F1.U1.growth_index"), 1.5, 0.01);
  CHECK_NEAR(log_value(log, "F1.U2.growth_index"), 1.5, 0.01);
  CHECK_NEAR(log_value(log, "F1.U3.growth_index"), 1.0, 0.01);
}

TEST_CASE("P24 growth fold: relief decreases upwards") {
  auto m = build_text(model_yaml("", "structure:\n  regional: {type: flat}\n  primitives:\n    - {type: anticline, name: A, crest: " +
                                         q("[5000, 5000] m") + ", axis_azimuth: " + q("0 deg") + ", half_width_left: " + q("1500 m") +
                                         ", half_width_right: " + q("1500 m") + ", relief: {H4: " + q("300 m") +
                                         "}, growth: {onset: H3, end: H1}}\n"));
  const double d0[] = {1000, 1200, 1400, 1600}, rel[] = {0, 150, 300, 300};
  for (int h = 0; h < 4; ++h) CHECK_NEAR(depth_faulted(*m, h, 5000, 5000), d0[h] - rel[h], 0.01);
}

namespace {
std::string clino(const std::string& profile) {
  return doc("{name: H1}, {name: C3}, {name: C2}, {name: C1}, {name: H4}, {name: H5}",
             top("H1", 1000) + ", " + top("H4", 1600) + ", " + top("H5", 1800), units_for({"H1", "C3", "C2", "C1", "H4", "H5"}),
             "    - {name: S3, horizons: [H1], base: conformable}\n"
             "    - {name: S2, horizons: [C3, C2, C1], base: conformable, clinoforms: {profile: " + profile + ", height: " + q("300 m") +
                 ", width: " + q("2000 m") + ", direction: " + q("90 deg") + ", first_break: " + q("[3000, 5000] m") +
                 ", progradation: " + q("1500 m") + "}}\n"
             "    - {name: S1, horizons: [H4, H5], base: conformable}\n");
}
double sig(double u) { return 0.5 * (1 - std::tanh(2 * u)); }
double obl(double u) { return std::min(1.0, std::max(0.0, 0.5 - u)); }
}  // namespace

TEST_CASE("P24 clinoforms against their formulas") {
  for (std::string prof : {"sigmoidal", "oblique"}) {
    BuildLog log;
    auto m = build_text(clino(prof), nullptr, &log);
    double worst = 0;
    for (int i = 0; i < 3; ++i)
      for (double x = 100; x <= 9900; x += 100) {
        double u = (x - 3000 - i * 1500) / 2000, p = prof == "sigmoidal" ? sig(u) : obl(u);
        worst = std::max(worst, std::fabs(depth_faulted(*m, 3 - i, x, 7000) - (1600 - 300 * p)));
      }
    std::cout << "  " << prof << ": worst misfit " << fmt_fixed(worst, 5) << " m\n";
    CHECK(worst <= 0.01);
    CHECK_NEAR(log_value(log, "S2.clinoform_foreset_dip_max"), rad2deg(std::atan(300.0 / 2000.0)), 0.1);
  }
}

TEST_CASE("P23/P24 refusals (F20, F22, F23)") {
  std::string bad = angular();
  bad.replace(bad.find("horizons: [H1, H2]"), 18, "horizons: [H2, H1]");
  auto r1 = validate(bad);
  CHECK(!r1.ok());
  CHECK(r1.text().find("F20") != std::string::npos);
  auto r2 = validate(model_yaml(planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", "max: " + q("200 m"),
                                             "dip_slip", "{type: rigid}", "    growth: {onset: H3, end: H1}\n")));
  CHECK(!r2.ok());
  CHECK(r2.text().find("F22") != std::string::npos);
  std::string c = clino("sigmoidal");
  c.replace(c.find("    - {name: S1, horizons: [H4, H5], base: conformable}\n"), 55, "");
  c.replace(c.find("[C3, C2, C1]"), 12, "[C3, C2, C1, H4, H5]");
  auto r3 = validate(c);
  CHECK(!r3.ok());
  CHECK(r3.text().find("F23") != std::string::npos);
}

STX_TEST_MAIN()
