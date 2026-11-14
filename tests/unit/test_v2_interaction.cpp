// Stage 2A, Prompt 17: drag, splays with branch lines and throw transfer, relay ramps, and the
// displacement-length check (STX-TN-05 sections 3 and 4).
#include "v2_fixture.hpp"

using namespace v2fx;

namespace {

// Upward normal of a plane dipping `dip` towards strike + 90 (right-hand rule), as (east, north, up).
Vec3 up_normal(double strike_deg, double dip_deg) {
  double phi = deg2rad(strike_deg + 90), d = deg2rad(dip_deg);
  return {std::sin(d) * std::sin(phi), std::sin(d) * std::cos(phi), std::cos(d)};
}

const std::string kThrow100 = "max: " + q("100 m");

}  // namespace

TEST_CASE("P17 branch line of a splay on its master") {
  std::string f = planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", kThrow100, "dip_slip") +
                  planar_fault("F2", "2", "[5150, 5000] m", "1300 m", "20 deg", "75 deg", "100000 m", "max: " + q("60 m"),
                               "dip_slip", "{type: rigid}", "    abuts: {fault: F1, side: hanging_wall, transfer: true}\n");
  BuildLog log;
  auto m = build_text(model_yaml(f, "", {1000, 1200, 1400, 1600, 2400}), nullptr, &log);
  Vec3 d = cross(up_normal(0, 60), up_normal(20, 75));
  if (d.z > 0) d = d * -1.0;  // plunging down
  double trend = std::fmod(rad2deg(std::atan2(d.x, d.y)) + 360.0, 360.0);
  double plunge = rad2deg(std::asin(-d.z / norm(d)));
  CHECK_NEAR(log_value(log, "F2.branch.F1.trend"), trend, 0.1);
  CHECK_NEAR(log_value(log, "F2.branch.F1.plunge"), plunge, 0.1);
}

TEST_CASE("P17 throw transfer: below the branch line the master carries both throws") {
  std::string f = planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", kThrow100, "dip_slip") +
                  planar_fault("F2", "2", "[5150, 5000] m", "1300 m", "20 deg", "75 deg", "100000 m", "max: " + q("60 m"),
                               "dip_slip", "{type: rigid}", "    abuts: {fault: F1, side: hanging_wall, transfer: true}\n");
  desc::Description d;
  auto m = build_text(model_yaml(f, "", {1000, 1200, 1400, 1600, 2400}), &d);
  // At y = 5000 the branch line lies near 1800 m: H1 (1000 m) is cut by both faults, H5 (2400 m)
  // only by the master below the branch line.
  CHECK_NEAR(depth_faulted(*m, 0, 8000, 5000) - depth_faulted(*m, 0, 3000, 5000), 160, 1.0);  // both faults
  CHECK_NEAR(depth_faulted(*m, 0, 5000, 5000) - depth_faulted(*m, 0, 3000, 5000), 100, 1.0);  // between them
  CHECK_NEAR(depth_faulted(*m, 4, 8000, 5000) - depth_faulted(*m, 4, 3000, 5000), 160, 1.0);  // master alone
  d.wells.clear();
  auto rep = checks_for(*m, d);
  auto* c11 = find_check(rep, "C11");
  CHECK(c11 != nullptr);
  CHECK(c11 && c11->passed);
}

TEST_CASE("P17 relay ramp: summed throw through the overlap equals one fault's throw") {
  std::string f = planar_fault("F1", "1", "[4000, 0] m", "1300 m", "0 deg", "60 deg", "12000 m", kThrow100, "dip_slip",
                               "{type: rigid}", "    relay: {fault: F2}\n") +
                  planar_fault("F2", "2", "[5000, 10000] m", "1300 m", "0 deg", "60 deg", "12000 m", kThrow100, "dip_slip");
  auto m = build_text(model_yaml(f));
  double worst = 0;
  for (double y = 1000; y <= 9000; y += 200)
    worst = std::max(worst, std::fabs(depth_faulted(*m, 1, 8000, y) - depth_faulted(*m, 1, 2000, y) - 100));
  std::cout << "  worst summed-throw error along strike " << fmt_fixed(worst, 3) << " m\n";
  CHECK(worst <= 1.0);
  // Inside the overlap the ramp between the faults sits between footwall and hanging wall.
  double ramp = depth_faulted(*m, 1, 4600, 5000) - depth_faulted(*m, 1, 2000, 5000);
  CHECK(ramp > 5 && ramp < 95);
}

TEST_CASE("P17 normal drag on both walls against d(sigma) = f T (1 - sigma / W)^2") {
  std::string drag = "    drag: {type: normal, fraction: " + q("0.2 fraction") + ", width: " + q("500 m") + "}\n";
  std::string f = planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", kThrow100, "dip_slip",
                               "{type: rigid}", drag);
  auto m = build_text(model_yaml(f));
  double s60 = std::sin(deg2rad(60)), cot = 1 / std::tan(deg2rad(60));
  auto xf = [&](double z) { return 5000 + (z - 1300) * cot; };
  auto dragf = [](double sigma) { double a = std::fabs(sigma); return a < 500 ? 0.2 * 100 * std::pow(1 - a / 500, 2) : 0.0; };
  double worst = 0;
  // Hanging wall: H2 at 1300 m after faulting; moved up the slip vector by d.
  for (double x = xf(1300) + 1; x <= xf(1300) + 800; x += 10) {
    double dd = dragf((x - xf(1300)) * s60);
    worst = std::max(worst, std::fabs(depth_in_block(*m, 1, 1, x - dd * cot, 5000) - (1300 - dd)));
  }
  // Footwall: H2 at 1200 m; moved down the slip vector by d.
  for (double x = xf(1200) - 1; x >= xf(1200) - 800; x -= 10) {
    double dd = dragf((x - xf(1200)) * s60);
    worst = std::max(worst, std::fabs(depth_in_block(*m, 1, 0, x + dd * cot, 5000) - (1200 + dd)));
  }
  std::cout << "  worst drag misfit " << fmt_fixed(worst, 4) << " m\n";
  CHECK(worst <= 0.5);
}

TEST_CASE("P17 drag that would carry a horizon across the fault is refused (F13)") {
  std::string drag = "    drag: {type: normal, fraction: " + q("0.6 fraction") + ", width: " + q("500 m") + "}\n";
  std::string f = planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", kThrow100, "dip_slip",
                               "{type: rigid}", drag);
  auto r = validate(model_yaml(f));
  CHECK(!r.ok());
  CHECK(r.text().find("F13") != std::string::npos);
}

TEST_CASE("P17 C10 displacement-length scaling: warns outside the stated range") {
  for (auto [len, thr, ok] : {std::tuple{"400 m", "200 m", false}, std::tuple{"5000 m", "100 m", true}}) {
    std::string f = planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", len, "max: " + q(thr), "dip_slip",
                                 "{type: rigid}", "", flat_top("0.25"));
    desc::Description d;
    auto m = build_text(model_yaml(f), &d);
    d.wells.clear();
    auto rep = checks_for(*m, d);
    auto* c10 = find_check(rep, "C10", "F1");
    REQUIRE(c10 != nullptr);
    std::cout << "  C10 " << len << " / " << thr << ": ratio " << c10->value << (c10->passed ? " PASS" : " WARN") << "\n";
    CHECK(c10->passed == ok);
    if (!ok) CHECK(c10->warn_only);
  }
}

STX_TEST_MAIN()
