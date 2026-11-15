// Stage 2F (Prompts 32 and 33): juxtaposition, shale gouge ratio, clay smear potential, fault-zone
// thickness, seal capacity, stress, slip and dilation tendency, curvature (STX-TN-10).
#include "v2_fixture.hpp"
#include "stratonix/repr/faultrock.hpp"

using namespace v2fx;

namespace {

const double kSin60 = std::sin(deg2rad(60.0));

/// Layer cake U1 sand (0.1), U2 shale (0.8), U3 sand (0.2) cut by a 60-degree normal fault.
std::string stack(double throw_m, const std::string& top_level = "", bool clay = true) {
  auto u = [&](const std::string& n, const std::string& t, const std::string& b, const std::string& v) {
    return "{name: " + n + ", top: " + t + ", base: " + b + (clay ? ", clay_fraction: " + q(v + " fraction") : std::string()) + "}";
  };
  std::string s = "stratonix_schema: 1\nmodel: {name: t, version: \"1.0.0\", status: DESIGN}\n"
                  "crs: {name: local, vertical_datum: datum, length_unit: m}\n"
                  "domain:\n  origin: " + q("[0, 0] m") + "\n  size: " + q("[10000, 10000] m") +
                  "\n  depth_range: " + q("[0, 3000] m") + "\n  grid_spacing: " + q("50 m") + "\n";
  s += "wells:\n  - {name: W1, role: control, head: " + q("[1000, 1000] m") + ", tops: {H1: " + q("1000 m") + ", H2: " + q("1200 m") +
       ", H3: " + q("1400 m") + ", H4: " + q("1600 m") + "}}\n";
  s += "stratigraphy:\n  reference: {well: W1}\n  horizons: [{name: H1}, {name: H2}, {name: H3}, {name: H4}]\n  units: [" +
       u("U1", "H1", "H2", "0.1") + ", " + u("U2", "H2", "H3", "0.8") + ", " + u("U3", "H3", "H4", "0.2") + "]\n";
  s += "structure:\n  regional: {type: flat}\n";
  s += "faults:\n" + planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", "max: " + q(fmt_fixed(throw_m, 1) + " m"),
                                  "dip_slip") + "frame: {fixed: {well: W1}}\n";
  return s + top_level;
}
/// The fault point whose footwall (unmoved) depth is D.
Vec3 on_fault(double D) { return {5000 + (D - 1300) / std::tan(deg2rad(60.0)), 5000, -D}; }

repr::FaultRockPoint at(const Model& m, const desc::Description& d, double D) {
  return repr::fault_rock_at(m, 0, on_fault(D), fault_rock_input(d, m), stress_input(d));
}
const std::string kRock = "fault_rock:\n  shale_threshold: " + q("0.5 fraction") + "\n  thickness_displacement_ratio: " + q("0.01 ratio") +
                          "\n  seal_calibration: {type: sgr_power, a: " + q("27 ratio") + ", c: " + q("0.5 ratio") + "}\n";
std::string regime() {
  return "stress:\n  regime: {sv_gradient: " + q("25 MPa/km") + ", shmax_ratio: " + q("0.8 ratio") + ", shmin_ratio: " + q("0.6 ratio") +
         ", shmax_azimuth: " + q("0 deg") + "}\n  pore_pressure_gradient: " + q("10 MPa/km") + "\n  friction: " + q("0.6 ratio") + "\n";
}

}  // namespace

TEST_CASE("P32 shale gouge ratio of a known stack") {
  desc::Description d;
  auto m = build_text(stack(200, kRock), &d);
  auto p = at(*m, d, 1300);
  std::cout << "  separation " << fmt_fixed(p.separation, 3) << " m, SGR " << fmt_fixed(p.sgr, 4) << "\n";
  CHECK_NEAR(p.separation, 200, 0.5);
  CHECK_NEAR(p.sgr, 0.45, 0.005);
  CHECK_NEAR(p.thickness, 0.01 * 200 / kSin60, 0.01);
  CHECK_NEAR(p.seal_bar, std::pow(10.0, 45.0 / 27.0 - 0.5), 0.15);
}

TEST_CASE("P32 juxtaposition and clay smear potential") {
  desc::Description d;
  auto m = build_text(stack(100, kRock), &d);
  auto p = at(*m, d, 1250);
  CHECK(p.fw_unit == 2);  // U2 (model unit 2: between H2 and H3)
  CHECK(p.hw_unit == 1);  // U1
  auto c = at(*m, d, 1550);
  const double csp = 200.0 * 200.0 * kSin60 / 150.0 + 200.0 * 200.0 * kSin60 / 50.0;
  std::cout << "  CSP " << fmt_fixed(c.csp, 2) << " (hand " << fmt_fixed(csp, 2) << ")\n";
  CHECK_NEAR(c.csp, csp, 0.01 * csp);
}

TEST_CASE("P32 missing inputs are reported, never defaulted") {
  desc::Description d;
  auto m = build_text(stack(200, "", false), &d);
  auto p = at(*m, d, 1300);
  CHECK(std::isnan(p.sgr));
  CHECK(std::isnan(p.seal_bar));
  CHECK(std::isnan(p.thickness));
  CHECK(!p.missing.empty());
}

TEST_CASE("P33 slip and dilation tendency of a known plane") {
  desc::Description d;
  auto m = build_text(stack(200, regime()), &d);
  auto p = at(*m, d, 2000);
  std::cout << "  sigma_n " << fmt_fixed(p.sigma_n / 1e6, 3) << " MPa, tau " << fmt_fixed(p.tau / 1e6, 3) << " MPa, Ts " << fmt_fixed(p.ts, 4)
            << ", Td " << fmt_fixed(p.td, 4) << ", dP " << fmt_fixed(p.dp / 1e6, 4) << " MPa\n";
  CHECK_NEAR(p.sigma_n / 1e6, 35.0, 0.01);
  CHECK_NEAR(p.tau / 1e6, 8.660, 0.01);
  CHECK_NEAR(p.ts, 0.5774, 0.001);
  CHECK_NEAR(p.td, 0.75, 0.001);
  CHECK_NEAR(p.dp / 1e6, 15.0 - 8.660 / 0.6, 0.01);
}

TEST_CASE("P33 the same stress as a tensor") {
  std::string t = "stress:\n  tensor_gradient: " + q("[15, 20, 25, 0, 0, 0] MPa/km") + "\n  pore_pressure_gradient: " + q("10 MPa/km") +
                  "\n  friction: " + q("0.6 ratio") + "\n";
  desc::Description d;
  auto m = build_text(stack(200, t), &d);
  auto p = at(*m, d, 2000);
  CHECK_NEAR(p.ts, 0.5774, 0.001);
  CHECK_NEAR(p.td, 0.75, 0.001);
}

TEST_CASE("P33 curvature and fracture trend at an anticline crest") {
  auto m = build_text(model_yaml("", "structure:\n  regional: {type: flat}\n  primitives:\n    - {type: anticline, name: A, crest: " +
                                         q("[5000, 5000] m") + ", axis_azimuth: " + q("0 deg") + ", half_width_left: " + q("1500 m") +
                                         ", half_width_right: " + q("1500 m") + ", relief: {H1: " + q("300 m") + "}}\n"));
  auto c = repr::curvature(*m, 0, GridSpec{4900, 4900, 10, 10, 21, 21}, 1e-4);
  const double kexp = -300 * kPi * kPi / (2 * 1500.0 * 1500.0);
  std::cout << "  k_min " << c.kmin.at(10, 10) << " (expected " << kexp << "), Gaussian " << c.gaussian.at(10, 10) << ", trend "
            << c.trend.at(10, 10) << " deg\n";
  CHECK_NEAR(c.kmin.at(10, 10), kexp, 0.01 * std::fabs(kexp));
  CHECK(std::fabs(c.gaussian.at(10, 10)) < 1e-9);
  CHECK(std::fabs(std::fmod(c.trend.at(10, 10) + 180.0, 180.0)) < 1.0 || std::fabs(std::fmod(c.trend.at(10, 10) + 180.0, 180.0) - 180) < 1.0);
}

TEST_CASE("P32/P33 checks C19 and C20") {
  desc::Description d;
  auto m = build_text(stack(200, kRock + regime()), &d);
  auto rep = checks_for(*m, d);
  auto c19 = find_check(rep, "C19"), c20 = find_check(rep, "C20");
  REQUIRE(c19 != nullptr);
  REQUIRE(c20 != nullptr);
  std::cout << "  C19 " << c19->value << ", C20 " << c20->value << "\n";
  CHECK(c19->passed);
  CHECK(c20->passed);
}

TEST_CASE("P32/P33 refusals (F27, F28)") {
  std::string a = stack(200, kRock);
  a.replace(a.find("0.8 fraction"), 12, "1.2 fraction");
  auto r1 = validate(a);
  CHECK(!r1.ok());
  CHECK(r1.text().find("F27") != std::string::npos);
  std::string b = stack(200, regime() + "  tensor_gradient: " + q("[15, 20, 25, 0, 0, 0] MPa/km") + "\n");
  auto r2 = validate(b);
  CHECK(!r2.ok());
  CHECK(r2.text().find("F28") != std::string::npos);
}

STX_TEST_MAIN()
