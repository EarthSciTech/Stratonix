// Stage 2A, Prompt 16: curved (listric and ramp-flat) faults with vertical- and inclined-shear
// hanging walls, against the inclined-shear construction of STX-TN-05 section 2.
#include "v2_fixture.hpp"

using namespace v2fx;

namespace {

// Fault profiles in the dip plane: depth (positive down) at horizontal distance u along the dip
// direction from the reference point.
struct Profile {
  std::function<double(double)> zf;
};
Profile listric(double z0, double dip_deg, double zd) {
  double t = std::tan(deg2rad(dip_deg));
  return {[=](double u) { return zd - (zd - z0) * std::exp(-u * t / (zd - z0)); }};
}
Profile ramp_flat(std::vector<std::pair<double, double>> pts) {
  return {[pts](double u) {
    auto seg = [&](std::size_t k, double uu) {
      auto [u0, z0] = pts[k];
      auto [u1, z1] = pts[k + 1];
      return z0 + (z1 - z0) * (uu - u0) / (u1 - u0);
    };
    if (u <= pts.front().first) return seg(0, u);
    for (std::size_t k = 0; k + 1 < pts.size(); ++k)
      if (u <= pts[k + 1].first) return seg(k, u);
    return seg(pts.size() - 2, u);
  }};
}

// The construction: shear lines u(z) = xi + s z tan(alpha); G(xi) = depth where the shear line meets
// the fault.
struct Construction {
  Profile p;
  double s = -1, ta = 0, heave = 0, x_ref = 3000;
  double G(double xi) const {
    if (ta == 0) return p.zf(xi);
    return bisect([&](double z) { return p.zf(xi + s * z * ta) - z; }, -10000, 10000);
  }
  // Forward: original (u, z) in the hanging wall -> present (x, z).
  void forward(double u, double z, double& x1, double& z1) const {
    double xi = u - s * z * ta;
    double xi1 = xi + heave;
    z1 = z + G(xi1) - G(xi);
    x1 = x_ref + xi1 + s * z1 * ta;
  }
  bool in_hanging_wall(double u, double z) const { return z < G(u - s * z * ta); }
};

std::string curved_fault(const std::string& geometry, const std::string& hw, const std::string& heave,
                         const std::string& upper = kFree) {
  return "  - name: F1\n    age: " + q("1 ordinal") + "\n    geometry: " + geometry + "\n    extent: {upper_tip: " + upper +
         ", lower_tip: " + kFree + ", length: " + q("100000 m") + "}\n    throw: {heave: " + q(heave) +
         ", profile_strike: " + flat_top("0.01") + "}\n    slip_mode: dip_slip\n    hanging_wall: " + hw + "\n";
}
const std::string kListric = "{type: listric, reference_point: " + q("[3000, 5000] m") + ", reference_depth: " + q("0 m") +
                             ", strike: " + q("0 deg") + ", dip: " + q("60 deg") + ", detachment_depth: " + q("2500 m") + "}";
const std::string kRampFlat = "{type: ramp_flat, reference_point: " + q("[3000, 5000] m") + ", reference_depth: " + q("0 m") +
                              ", strike: " + q("0 deg") + ", profile: " +
                              q("[[0, 0], [635.085, 1100], [2000, 1100], [2800, 1900], [6000, 1900]] m") + "}";
const std::string kVertical = "{type: vertical_shear}";
std::string inclined(const std::string& angle, const std::string& sense) {
  return "{type: inclined_shear, shear_angle: " + q(angle + " deg") + ", shear_sense: " + sense + "}";
}

// Worst misfit of every horizon in the hanging wall against the construction.
double worst_misfit(const Model& m, const Construction& c, const std::vector<double>& depths) {
  double worst = 0;
  int n = 0;
  for (std::size_t h = 0; h < depths.size(); ++h)
    for (double u = -1500; u <= 6500; u += 25) {
      if (!c.in_hanging_wall(u, depths[h])) continue;
      double x1, z1;
      c.forward(u, depths[h], x1, z1);
      if (x1 < 300 || x1 > 9700) continue;
      double zm = depth_in_block(m, static_cast<int>(h), 1, x1, 5000);
      worst = std::max(worst, std::fabs(zm - z1));
      ++n;
    }
  std::cout << "  " << n << " hanging-wall points, worst misfit " << fmt_fixed(worst, 4) << " m\n";
  return n > 50 ? worst : 1e9;
}

const std::vector<double> kDepths = {1000, 1200, 1400, 1600};

}  // namespace

TEST_CASE("P16 listric fault, vertical shear: rollover against the construction") {
  auto m = build_text(model_yaml(curved_fault(kListric, kVertical, "600 m")));
  Construction c{listric(0, 60, 2500), -1, 0, 600};
  CHECK(worst_misfit(*m, c, kDepths) <= 0.5);
}

TEST_CASE("P16 listric fault, 30 degree antithetic shear: rollover against the construction") {
  auto m = build_text(model_yaml(curved_fault(kListric, inclined("30", "antithetic"), "600 m")));
  Construction c{listric(0, 60, 2500), -1, std::tan(deg2rad(30)), 600};
  CHECK(worst_misfit(*m, c, kDepths) <= 0.5);
}

TEST_CASE("P16 ramp-flat fault, vertical shear: hanging wall across the kinks") {
  auto m = build_text(model_yaml(curved_fault(kRampFlat, kVertical, "500 m")));
  Construction c{ramp_flat({{0, 0}, {635.085, 1100}, {2000, 1100}, {2800, 1900}, {6000, 1900}}), -1, 0, 500};
  CHECK(worst_misfit(*m, c, kDepths) <= 0.5);
}

TEST_CASE("P16 round trip: forward by the construction, restored by the model") {
  auto m = build_text(model_yaml(curved_fault(kListric, inclined("30", "antithetic"), "600 m")));
  Construction c{listric(0, 60, 2500), -1, std::tan(deg2rad(30)), 600};
  std::vector<Vec3> present, original;
  for (double u = -500; u <= 5500; u += 250)
    for (double z = 200; z <= 2400; z += 100) {
      if (!c.in_hanging_wall(u, z) || z > c.G(u - c.s * z * c.ta) - 20) continue;
      double x1, z1;
      c.forward(u, z, x1, z1);
      if (x1 < 300 || x1 > 9700) continue;
      present.push_back({x1, 5000, -z1});
      original.push_back({3000 + u, 5000, -z});
    }
  REQUIRE(present.size() > 100);
  std::vector<Vec3> back(present.size());
  m->restore(present, back);
  double worst = 0;
  for (std::size_t i = 0; i < back.size(); ++i) worst = std::max(worst, norm(back[i] - original[i]));
  std::cout << "  " << back.size() << " points, worst round-trip error " << fmt_fixed(worst, 5) << " m\n";
  CHECK(worst <= 0.01);
}

TEST_CASE("P16 conservation: unit thickness along the shear direction is unchanged") {
  // Vertical shear: vertical thickness of every unit in the hanging wall stays 200 m.
  {
    auto m = build_text(model_yaml(curved_fault(kListric, kVertical, "600 m")));
    double worst = 0;
    for (double x = 4200; x <= 9500; x += 100)
      for (int h = 0; h + 1 < 4; ++h)
        worst = std::max(worst, std::fabs(depth_in_block(*m, h + 1, 1, x, 5000) - depth_in_block(*m, h, 1, x, 5000) - 200));
    CHECK(worst <= 0.01);
  }
  // 30 degree antithetic shear: thickness measured along a shear line stays 200 / cos(30).
  {
    auto m = build_text(model_yaml(curved_fault(kListric, inclined("30", "antithetic"), "600 m")));
    double ta = std::tan(deg2rad(30)), worst = 0;
    for (double xi1 = 2500; xi1 <= 6000; xi1 += 250) {
      // Present shear line through xi1: x(z) = 3000 + xi1 - z tan(30).
      auto meet = [&](int h) {
        return bisect([&](double z) { return depth_in_block(*m, h, 1, 3000 + xi1 - z * ta, 5000) - z; }, 0, 3000);
      };
      double z1 = meet(0), z2 = meet(1);
      worst = std::max(worst, std::fabs((z2 - z1) / std::cos(deg2rad(30)) - 200 / std::cos(deg2rad(30))));
    }
    CHECK(worst <= 0.01);
  }
}

TEST_CASE("P16 refusals outside the admissible range (F10, F12)") {
  // Synthetic shear at 40 degrees on a fault 60 degrees steep: shear lines cross the fault twice.
  auto r1 = validate(model_yaml(curved_fault(kListric, inclined("40", "synthetic"), "600 m")));
  CHECK(!r1.ok());
  CHECK(r1.text().find("F12") != std::string::npos);
  // An upper tip inside a unit: the construction needs the fault to reach the model top.
  std::string tip = "{rule: within_unit, unit: U1, fraction: " + q("0.5 fraction") + "}";
  auto r2 = validate(model_yaml(curved_fault(kListric, kVertical, "600 m", tip)));
  CHECK(!r2.ok());
  CHECK(r2.text().find("F10") != std::string::npos);
  // Synthetic shear steeper than the fault is admissible.
  auto r3 = validate(model_yaml(curved_fault(kListric, inclined("20", "synthetic"), "600 m")));
  if (!r3.ok()) std::cerr << r3.text();
  CHECK(r3.ok());
}

TEST_CASE("P16 report: heave, throw and rollover width per horizon") {
  BuildLog log;
  auto m = build_text(model_yaml(curved_fault(kListric, kVertical, "600 m")), nullptr, &log);
  Profile p = listric(0, 60, 2500);
  for (int h = 0; h < 4; ++h) {
    std::string key = "F1.H" + std::to_string(h + 1);
    double zh = kDepths[static_cast<std::size_t>(h)];
    // Footwall cut-off: the fault at the horizon's depth; it moves to xi + heave.
    double xc = bisect([&](double u) { return p.zf(u) - zh; }, -5000, 8000);
    double throw_expected = p.zf(xc + 600) - zh;
    double xr = bisect([&](double u) { return p.zf(u + 600) - p.zf(u) - 1.0; }, xc, 20000);
    CHECK_NEAR(log_value(log, key + ".heave"), 600, 0.01);
    CHECK_NEAR(log_value(log, key + ".throw"), throw_expected, 0.5);
    CHECK_NEAR(log_value(log, key + ".rollover_width"), xr - xc, 5.0);
  }
}

STX_TEST_MAIN()
