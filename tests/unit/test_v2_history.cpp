// Stage 2E (Prompts 29 and 30): event timeline, phases and inversion, states at any event,
// restorability and balance, flexural-slip unfolding, backstripping (STX-TN-09).
#include "v2_fixture.hpp"
#include "stratonix/engines/history.hpp"
#include "stratonix/engines/structure.hpp"

using namespace v2fx;

namespace {

std::string age(double a) { return q(fmt_fixed(a, 1) + " Ma"); }
/// Deposition events for H1..H4 at the given ages, plus extra events.
std::string events(const std::vector<double>& ages, const std::string& extra) {
  std::string s = "events:\n";
  for (std::size_t h = 0; h < ages.size(); ++h)
    s += "  - {name: D" + std::to_string(h + 1) + ", type: deposition, horizon: H" + std::to_string(h + 1) + ", age: " + age(ages[h]) + "}\n";
  return s + extra;
}
std::string lr_fault(const std::string& extra = "") {
  return planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", "max: " + q("300 m"), "dip_slip",
                      "{type: layer_rigid}", extra);
}
/// Offset of horizon h across F1 (hanging wall east of the fault, regional depth d0).
double offset(const Model& m, int h, double d0) { return depth_faulted(m, h, 8000, 5000) - d0; }

ModelPtr inversion(desc::Description* d = nullptr) {
  // H4 100, H3 90 pre-rift; normal 300 m at 85; H2 80 syn-rift; reverse 150 m at 70; H1 60.
  std::string ev = events({60, 80, 90, 100},
                          "  - {name: Rift, type: fault, fault: F1, age: " + age(85) + ", throw: " + q("300 m") + ", sense: normal}\n"
                          "  - {name: Inversion, type: fault, fault: F1, age: " + age(70) + ", throw: " + q("150 m") + ", sense: reverse}\n");
  return build_text(model_yaml(lr_fault(), "", {1000, 1200, 1400, 1600}, ev), d);
}

}  // namespace

TEST_CASE("P29 events reproduce the stage 2C growth fault") {
  auto a = build_text(model_yaml(lr_fault("    growth: {onset: H4, end: H1}\n"), "", {1000, 1200, 1400, 1600}));
  // growth from the deposition of H4 (100 Ma) to that of H1 (70 Ma); horizons equally spaced in depth and time
  auto b = build_text(model_yaml(lr_fault(), "", {1000, 1200, 1400, 1600},
                                 events({70, 80, 90, 100}, "  - {name: Slip, type: fault, fault: F1, from: " + age(100) + ", to: " + age(70) + "}\n")));
  double worst = 0;
  for (int h = 0; h < 4; ++h)
    for (double x : {2000.0, 4800.0, 6000.0, 8000.0}) worst = std::max(worst, std::fabs(depth_faulted(*a, h, x, 5000) - depth_faulted(*b, h, x, 5000)));
  std::cout << "  worst difference " << fmt_fixed(worst, 4) << " m\n";
  CHECK(worst <= 0.01);
}

TEST_CASE("P29 inversion: net normal pre-rift, reverse syn-rift, none post-inversion") {
  auto m = inversion();
  const double d0[] = {1000, 1200, 1400, 1600}, off[] = {0, -150, 150, 150};
  for (int h = 0; h < 4; ++h) CHECK_NEAR(offset(*m, h, d0[h]), off[h], 0.5);
  for (int h = 0; h < 4; ++h) CHECK_NEAR(depth_faulted(*m, h, 2000, 5000), d0[h], 0.01);
}

TEST_CASE("P29 state between the phases and before faulting") {
  desc::Description d;
  auto m = inversion(&d);
  auto s75 = build_state(d, *m, 75.0);
  REQUIRE(s75->info().horizons.size() == 3);  // H1 not yet deposited
  CHECK_NEAR(offset(*s75, 0, 1200), 0, 0.5);    // H2 flat
  CHECK_NEAR(offset(*s75, 1, 1400), 300, 0.5);  // H3, H4 offset by the rift phase only
  CHECK_NEAR(offset(*s75, 2, 1600), 300, 0.5);
  auto s95 = build_state(d, *m, 88.0);  // after H3 (90 Ma), before the rift (85 Ma)
  REQUIRE(s95->info().horizons.size() == 2);
  CHECK_NEAR(offset(*s95, 0, 1400), 0, 0.01);
  CHECK_NEAR(offset(*s95, 1, 1600), 0, 0.01);
}

TEST_CASE("P29 fold growth by ages, and the state during it") {
  std::string st = "structure:\n  regional: {type: flat}\n  primitives:\n    - {type: anticline, name: A, crest: " + q("[5000, 5000] m") +
                   ", axis_azimuth: " + q("0 deg") + ", half_width_left: " + q("1500 m") + ", half_width_right: " + q("1500 m") +
                   ", relief: {H4: " + q("300 m") + "}}\n";
  desc::Description d;
  auto m = build_text(model_yaml("", st, {1000, 1200, 1400, 1600},
                                 events({60, 80, 90, 100}, "  - {name: Fold, type: fold, primitive: A, from: " + age(85) + ", to: " + age(75) + "}\n")),
                      &d);
  const double d0[] = {1000, 1200, 1400, 1600}, rel[] = {0, 150, 300, 300};
  for (int h = 0; h < 4; ++h) CHECK_NEAR(d0[h] - depth_faulted(*m, h, 5000, 5000), rel[h], 0.01);
  auto s80 = build_state(d, *m, 80.0);
  REQUIRE(s80->info().horizons.size() == 3);
  const double rel80[] = {0, 150, 150};
  for (int h = 0; h < 3; ++h) CHECK_NEAR(d0[h + 1] - depth_faulted(*s80, h, 5000, 5000), rel80[h], 0.01);
}

TEST_CASE("P30 restorability and balance") {
  desc::Description d;
  auto m = inversion(&d);
  CheckContext ctx(*m, d, make_well_paths(d), 100);
  auto rep = run_checks(ctx);
  auto c17 = find_check(rep, "C17");
  REQUIRE(c17 != nullptr);
  CHECK(c17->passed);
  auto c18 = find_check(rep, "C18", "H4");
  REQUIRE(c18 != nullptr);
  std::cout << "  C17 " << c17->value << " " << c17->unit << "; C18 H4 " << c18->value << " " << c18->unit << "\n";
  CHECK(c18->value <= 0.001);
}

TEST_CASE("P30 flexural-slip unfolding of a parallel anticline") {
  std::string a = "{type: anticline, name: A, crest: " + q("[5000, 5000] m") + ", axis_azimuth: " + q("0 deg") + ", half_width_left: " +
                  q("1500 m") + ", half_width_right: " + q("1500 m") + ", relief: {H1: " + q("300 m") + "}, mechanics: parallel}";
  auto m = build_text(model_yaml("", "structure:\n  regional: {type: flat}\n  primitives:\n    - " + a + "\n"));
  auto s = repr::make_section(*m, "S", {{1000, 5000, 0}, {9000, 5000, 0}}, 10, 2);
  auto u = repr::unfold_flexural(*m, s);
  REQUIRE(u.size() == 4);
  std::cout << "  bed lengths H1 " << fmt_fixed(u[0].traced_length, 1) << ", H2 " << fmt_fixed(u[1].traced_length, 1) << " m\n";
  CHECK(std::fabs(u[0].traced_length - u[1].traced_length) <= 0.005 * u[0].traced_length);
  double restored = 0;
  for (std::size_t k = 1; k < u[0].restored.size(); ++k) restored += norm(u[0].restored[k] - u[0].restored[k - 1]);
  CHECK_NEAR(restored, u[0].traced_length, 1e-6 * u[0].traced_length + 1e-6);
  CHECK(u[0].traced_length > 8000 + 10);  // folded: longer than the section
}

TEST_CASE("P30 backstripping: the state at a deposition decompacts the units below") {
  std::string comp = "  compaction:\n    geometry: depositional\n    surface_depth: " + q("0 m") + "\n    units:\n" +
                     "      - {unit: U1, law: exponential, surface_porosity: " + q("0.5 fraction") + ", decay: " + q("0.0005 1/m") + "}\n" +
                     "      - {unit: U2, law: exponential, surface_porosity: " + q("0.5 fraction") + ", decay: " + q("0.0005 1/m") + "}\n";
  std::string y = model_yaml("", "", {1000, 1200, 1400, 1600}, events({60, 80, 90, 100}, ""));
  y.replace(y.find("structure:"), 0, comp);
  desc::Description d;
  auto m = build_text(y, &d);
  auto hist = static_cast<const Structure*>(m->capability("structure"))->burial_history(1000, 1000);
  // hist steps: H3 (index 2) at the surface, then H2, then H1, then present (U3 rigid below H3)
  REQUIRE(hist.size() >= 3);
  auto s80 = build_state(d, *m, 80.0);  // H2 just deposited
  const double h2 = depth_faulted(*s80, 0, 1000, 1000), h3 = depth_faulted(*s80, 1, 1000, 1000);
  std::cout << "  state at 80 Ma: H2 " << fmt_fixed(h2, 2) << " m, H3 " << fmt_fixed(h3, 2) << " m\n";
  CHECK_NEAR(h3 - h2, hist[1].depth[2] - hist[1].depth[1], 0.1);  // U2 at its burial when H2 was at the surface
}

TEST_CASE("P29 refusals (F26)") {
  auto r1 = validate(model_yaml(lr_fault(), "", {1000, 1200, 1400, 1600},
                                events({60, 80, 90, 100}, "  - {name: X, type: fault, fault: F9, age: " + age(85) + "}\n")));
  CHECK(!r1.ok());
  CHECK(r1.text().find("F26") != std::string::npos);
  // stage 4B (I-030): phases on a rigid planar fault go level by level and are accepted
  const std::string phased = model_yaml(planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", "max: " + q("300 m"),
                                                     "dip_slip"),
                                        "", {1000, 1200, 1400, 1600},
                                        events({60, 80, 90, 100},
                                               "  - {name: A, type: fault, fault: F1, age: " + age(85) + ", throw: " + q("300 m") + ", sense: normal}\n"
                                               "  - {name: B, type: fault, fault: F1, age: " + age(70) + ", throw: " + q("150 m") + ", sense: reverse}\n"));
  auto r2 = validate(phased);
  CHECK(r2.ok());
  auto m2 = build_text(phased);
  CHECK(m2 != nullptr);
  auto r3 = validate(model_yaml(lr_fault(), "", {1000, 1200, 1400, 1600}, events({60, 80, 90}, "")));
  CHECK(!r3.ok());
  CHECK(r3.text().find("F26") != std::string::npos);
}

STX_TEST_MAIN()
