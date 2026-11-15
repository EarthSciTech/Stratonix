// Stage 2D (Prompts 26 and 27): salt bodies with volume-balanced withdrawal, the overburden's
// response, compaction, decompaction, burial history and porosity (STX-TN-08).
#include "v2_fixture.hpp"
#include "stratonix/engines/structure.hpp"

using namespace v2fx;

namespace {

std::string doc(const std::string& horizons, const std::string& tops, const std::string& units, const std::string& structure,
                const std::string& strat_extra = "") {
  std::string s = "stratonix_schema: 1\nmodel: {name: t, version: \"1.0.0\", status: DESIGN}\n"
                  "crs: {name: local, vertical_datum: datum, length_unit: m}\n"
                  "domain:\n  origin: " + q("[0, 0] m") + "\n  size: " + q("[10000, 10000] m") +
                  "\n  depth_range: " + q("[0, 3000] m") + "\n  grid_spacing: " + q("50 m") + "\n";
  s += "wells:\n  - {name: W1, role: control, head: " + q("[1000, 1000] m") + ", tops: {" + tops + "}}\n";
  s += "stratigraphy:\n  reference: {well: W1}\n  horizons: [" + horizons + "]\n  units: [" + units + "]\n" + strat_extra;
  s += structure;
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

// Overburden H1-H3, salt TS-BS (source t0), H5 below.
std::string salt_doc(double t0, const std::string& body) {
  return doc("{name: H1}, {name: H2}, {name: H3}, {name: TS}, {name: BS}, {name: H5}",
             top("H1", 800) + ", " + top("H2", 1000) + ", " + top("H3", 1200) + ", " + top("BS", 1600) + ", " + top("H5", 1800),
             units_for({"H1", "H2", "H3", "TS", "BS", "H5"}),
             "structure:\n  regional: {type: flat}\n  salt:\n    top: TS\n    base: BS\n    source_thickness: " + q(fmt_fixed(t0, 1) + " m") +
                 "\n    bodies:\n      - " + body + "\n");
}
std::string pillow(const std::string& extra = "") {
  return "{type: pillow, name: P1, centre: " + q("[5000, 5000] m") + ", radius: " + q("1500 m") + ", height: " + q("300 m") +
         ", withdrawal_width: " + q("2000 m") + extra + "}";
}
double raised(double r, double R) { return r < R ? 0.5 * (1 + std::cos(kPi * r / R)) : 0.0; }
double ring(double r, double R, double W) { return r > R && r < R + W ? std::pow(std::sin(kPi * (r - R) / W), 2) : 0.0; }
/// Withdrawal amplitude a balancing a pillow (R, A) over a ring of width W with source t0 (weld:
/// min(t0, a g)), by the test's own radial integration.
double balance(double R, double A, double W, double t0) {
  auto integ = [&](const std::function<double(double)>& f, double r0, double r1) {
    const int n = 20000;
    double s = 0, h = (r1 - r0) / n;
    for (int i = 0; i < n; ++i) { double r = r0 + (i + 0.5) * h; s += f(r) * 2 * kPi * r * h; }
    return s;
  };
  double vol = A * integ([&](double r) { return raised(r, R); }, 0, R);
  double lo = 0, hi = 1e6;
  for (int it = 0; it < 200; ++it) {
    double a = 0.5 * (lo + hi);
    double w = integ([&](double r) { return std::min(t0, a * ring(r, R, W)); }, R, R + W);
    (w < vol ? lo : hi) = a;
  }
  return 0.5 * (lo + hi);
}

// Exponential solid thickness between burial depths z1 < z2.
double solid(double z1, double z2, double p0, double c) { return (z2 - z1) - p0 / c * (std::exp(-c * z1) - std::exp(-c * z2)); }
// Base depth holding solid thickness S below top z1 (bisection).
double base_for(double z1, double S, double p0, double c) {
  double lo = z1, hi = z1 + 10 * S + 1000;
  for (int i = 0; i < 200; ++i) { double m = 0.5 * (lo + hi); (solid(z1, m, p0, c) < S ? lo : hi) = m; }
  return 0.5 * (lo + hi);
}
// Top depth holding solid thickness S above base z2 (bisection).
double top_for(double z2, double S, double p0, double c) {
  double lo = std::max(0.0, z2 - 10 * S - 1000), hi = z2;
  for (int i = 0; i < 200; ++i) { double m = 0.5 * (lo + hi); (solid(m, z2, p0, c) > S ? lo : hi) = m; }
  return 0.5 * (lo + hi);
}
std::string compaction(const std::string& geometry, const std::string& units) {
  return "  compaction:\n    geometry: " + geometry + "\n    surface_depth: " + q("0 m") + "\n    units:\n" + units;
}
std::string expo(const std::string& u) {
  return "      - {unit: " + u + ", law: exponential, surface_porosity: " + q("0.5 fraction") + ", decay: " + q("0.0005 1/m") + "}\n";
}

}  // namespace

TEST_CASE("P26 pillow: salt volume conserved, crest height") {
  BuildLog log;
  auto m = build_text(salt_doc(200, pillow()), nullptr, &log);
  CHECK_NEAR(depth_faulted(*m, 3, 5000, 5000), 1600 - 200 - 300, 0.01);
  CHECK_NEAR(structure_of(*m).salt_thickness(5000, 5000), 500, 0.01);
  double dv = 0, vol = 0;
  for (double x = 12.5; x < 10000; x += 25)
    for (double y = 12.5; y < 10000; y += 25) {
      double t = depth_faulted(*m, 4, x, y) - depth_faulted(*m, 3, x, y);
      dv += (t - 200) * 625;
      vol += 300 * raised(std::hypot(x - 5000, y - 5000), 1500) * 625;
    }
  std::cout << "  net salt volume change " << fmt_fixed(dv / vol * 100, 4) << " % of the pillow's volume\n";
  CHECK(std::fabs(dv) <= 0.01 * vol);
}

TEST_CASE("P26 withdrawal basin around a pillow") {
  auto m = build_text(salt_doc(200, pillow()));
  const double a = balance(1500, 300, 2000, 200);
  std::cout << "  withdrawal amplitude (test) " << fmt_fixed(a, 3) << " m\n";
  CHECK_NEAR(depth_faulted(*m, 1, 5000, 5000), 1000 - 300, 0.01);       // H2 pre-kinematic: rises with the crest
  CHECK_NEAR(depth_faulted(*m, 1, 5000 + 2500, 5000), 1000 + a, 1.0);   // ring middle: subsides by a
  CHECK_NEAR(depth_faulted(*m, 1, 500, 500), 1000, 0.01);               // far field unchanged
}

TEST_CASE("P26 weld: a thin source layer welds in the ring, volume still conserved") {
  BuildLog log;
  auto m = build_text(salt_doc(20, pillow()), nullptr, &log);
  const double a = balance(1500, 300, 2000, 20);
  CHECK(a > 20);
  CHECK_NEAR(depth_faulted(*m, 3, 7500, 5000), depth_faulted(*m, 4, 7500, 5000), 0.01);  // welded at the ring middle
  CHECK(log_value(log, "salt.weld_area") > 0);
  double dv = 0, vol = 0;
  for (double x = 12.5; x < 10000; x += 25)
    for (double y = 12.5; y < 10000; y += 25) {
      dv += (depth_faulted(*m, 4, x, y) - depth_faulted(*m, 3, x, y) - 20) * 625;
      vol += 300 * raised(std::hypot(x - 5000, y - 5000), 1500) * 625;
    }
  CHECK(std::fabs(dv) <= 0.01 * vol);
}

TEST_CASE("P26 syn-kinematic units thin over the crest") {
  auto m = build_text(salt_doc(200, pillow(", growth: {onset: H3, end: H1}")));
  // f(H1) = 0, f(H2) = 0.5, f(H3) = 1: each unit thinner by (f_base - f_top) * 300 = 150 m
  CHECK_NEAR(depth_faulted(*m, 1, 5000, 5000) - depth_faulted(*m, 0, 5000, 5000), 200 - 150, 0.01);
  CHECK_NEAR(depth_faulted(*m, 2, 5000, 5000) - depth_faulted(*m, 1, 5000, 5000), 200 - 150, 0.01);
  CHECK_NEAR(depth_faulted(*m, 0, 5000, 5000), 800, 0.01);  // post-kinematic H1 flat
}

TEST_CASE("P26 diapir pierces its syn-kinematic overburden") {
  std::string body = "{type: diapir, name: D1, centre: " + q("[5000, 5000] m") + ", radius: " + q("600 m") + ", flank_width: " +
                     q("100 m") + ", height: " + q("900 m") + ", withdrawal_width: " + q("3000 m") + ", growth: {onset: TS, end: H1}}";
  auto m = build_text(salt_doc(300, body));
  const double ts = depth_faulted(*m, 3, 5000, 5000);
  CHECK_NEAR(ts, 1600 - 300 - 900, 1.0);
  for (int h = 0; h < 3; ++h) CHECK_NEAR(depth_faulted(*m, h, 5000, 5000), ts, 0.01);  // clipped onto top salt
  CHECK_NEAR(depth_faulted(*m, 0, 500, 500), 800, 0.01);
  CHECK(m->info().units[1].termination.find("salt") != std::string::npos);
}

TEST_CASE("P27 decompaction against the solid-thickness formula") {
  auto m = build_text(doc("{name: H1}, {name: H2}, {name: H3}", top("H1", 1000) + ", " + top("H2", 1200) + ", " + top("H3", 1500),
                          units_for({"H1", "H2", "H3"}), "structure:\n  regional: {type: flat}\n",
                          compaction("present", expo("Unit_1"))));
  const double S = solid(1000, 1200, 0.5, 0.0005);
  const double t0 = base_for(0, S, 0.5, 0.0005);
  std::cout << "  solid thickness " << fmt_fixed(S, 3) << " m; at the surface " << fmt_fixed(t0, 3) << " m\n";
  const Structure& st = structure_of(*m);
  CHECK_NEAR(st.solid_thickness(1, 3000, 3000), S, 0.01);
  CHECK_NEAR(st.decompacted_thickness(1, 3000, 3000, 0), t0, 0.01);
  CHECK_NEAR(st.decompacted_thickness(1, 3000, 3000, 2500), base_for(2500, S, 0.5, 0.0005) - 2500, 0.01);
  CHECK_NEAR(st.porosity(1, 3000, 3000), (200 - S) / 200, 1e-4);
}

namespace {
std::string drape() {
  return doc("{name: H1}, {name: H2}, {name: B}, {name: B2}", top("H1", 1000) + ", " + top("H2", 1300) + ", " + top("B", 1600) + ", " + top("B2", 2000),
             units_for({"H1", "H2", "B", "B2"}),
             "structure:\n  regional: {type: flat}\n  primitives:\n    - {type: dome, name: High, series: S1, crest: " + q("[5000, 5000] m") +
                 ", radius_major: " + q("2500 m") + ", relief: {B: " + q("400 m") + "}}\n",
             "  series:\n    - {name: S2, horizons: [H1, H2], base: onlap}\n    - {name: S1, horizons: [B, B2], base: conformable}\n" +
                 compaction("depositional", expo("Unit_1") + expo("Unit_2")));
}
}  // namespace

TEST_CASE("P27 differential compaction: drape over a rigid high") {
  auto m = build_text(drape());
  // the test's own stacking from the rigid base B upwards
  auto present_h1 = [&](double dB, double dH2, double dH1) {
    double S2 = dB - dH2 > 0 ? solid(0, dB - dH2, 0.5, 0.0005) : 0, S1 = solid(0, dH2 - dH1, 0.5, 0.0005);
    double z2 = S2 > 0 ? top_for(dB, S2, 0.5, 0.0005) : dB;
    return top_for(z2, S1, 0.5, 0.0005);
  };
  const double off = present_h1(1600, 1300, 1000), crest = present_h1(1200, 1200, 1000);
  std::cout << "  H1 present: off the high " << fmt_fixed(off, 2) << " m, over the crest " << fmt_fixed(crest, 2) << " m\n";
  CHECK_NEAR(depth_faulted(*m, 0, 500, 500), off, 0.1);
  CHECK_NEAR(depth_faulted(*m, 0, 5000, 5000), crest, 0.1);
  CHECK(off > crest + 10);                                 // drape: the top subsides more off the high
  CHECK_NEAR(depth_faulted(*m, 2, 5000, 5000), 1200, 0.01);  // the rigid high stays where it is
}

TEST_CASE("P27 burial history against an independent backstripping") {
  auto m = build_text(drape());
  auto hist = structure_of(*m).burial_history(500, 500);
  REQUIRE(hist.size() == 4);  // B, H2 and H1 at the surface, then the present
  const double S2 = solid(0, 300, 0.5, 0.0005), S1 = solid(0, 300, 0.5, 0.0005);
  CHECK_NEAR(hist[0].depth[2], 0, 1e-9);
  CHECK_NEAR(hist[1].depth[2], base_for(0, S2, 0.5, 0.0005), 0.1);
  double z2 = base_for(0, S1, 0.5, 0.0005);
  CHECK_NEAR(hist[2].depth[1], z2, 0.1);
  CHECK_NEAR(hist[2].depth[2], base_for(z2, S2, 0.5, 0.0005), 0.1);
  CHECK(std::isnan(hist[0].depth[0]));
  CHECK_NEAR(hist[3].depth[0], depth_faulted(*m, 0, 500, 500), 0.1);
}

TEST_CASE("P26/P27 refusals (F24, F25)") {
  std::string bad = salt_doc(200, pillow());
  bad.replace(bad.find("    top: TS"), 11, "    top: H3");
  auto r1 = validate(bad);
  CHECK(!r1.ok());
  CHECK(r1.text().find("F24") != std::string::npos);
  std::string c = doc("{name: H1}, {name: H2}", top("H1", 1000) + ", " + top("H2", 1200), units_for({"H1", "H2"}),
                      "structure:\n  regional: {type: flat}\n", compaction("present", expo("Unit_1")));
  c.replace(c.find("0.5 fraction"), 12, "1.0 fraction");
  auto r2 = validate(c);
  CHECK(!r2.ok());
  CHECK(r2.text().find("F25") != std::string::npos);
}

STX_TEST_MAIN()
