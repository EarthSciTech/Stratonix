// Structure (Prompt 3) and fault (Prompt 4) benchmarks with analytically known answers.
#include <filesystem>
#include <map>

#include "../stx_test.hpp"
#include "stratonix/checks/checks.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/parametric.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/model/model_suite.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

using namespace stratonix;
namespace fs = std::filesystem;

static const char* P = " | user_instruction | test";
static std::string q(const std::string& v) { return "\"" + v + P + "\""; }

/// Build a description text: 4 horizons at 1000/1200/1400/1600 m, domain 10 x 10 km x 3 km.
static std::string model_yaml(const std::string& structure, const std::string& faults, const std::string& extra = "",
                              const std::string& wells_extra = "") {
  std::string s = "stratonix_schema: 1\nmodel: {name: t, version: \"1.0.0\", status: DESIGN}\n"
                  "crs: {name: local, vertical_datum: datum, length_unit: m}\n"
                  "domain:\n  origin: " + q("[0, 0] m") + "\n  size: " + q("[10000, 10000] m") +
                  "\n  depth_range: " + q("[0, 3000] m") + "\n  grid_spacing: " + q("50 m") + "\n";
  s += "wells:\n  - {name: W1, role: control, head: " + q("[1000, 1000] m") + ", tops: {H1: " + q("1000 m") + ", H2: " +
       q("1200 m") + ", H3: " + q("1400 m") + ", H4: " + q("1600 m") + "}}\n" + wells_extra;
  s += "stratigraphy:\n  reference: {well: W1}\n  horizons: [{name: H1}, {name: H2}, {name: H3}, {name: H4}]\n"
       "  units: [{name: U1, top: H1, base: H2}, {name: U2, top: H2, base: H3}, {name: U3, top: H3, base: H4}]\n";
  s += structure;
  if (!faults.empty()) s += "faults:\n" + faults + "frame: {fixed: {well: W1}}\n";
  s += extra;
  return s;
}

static ModelPtr build_text(const std::string& text, desc::Description* dout = nullptr, BuildLog* log = nullptr,
                           const fs::path& dir = stxtest::tmp_dir()) {
  register_all_plugins();
  desc::Description d;
  auto rep = desc::validate_description_text(text, dir, &d);
  if (!rep.ok()) std::cerr << rep.text() << text;
  if (!rep.ok()) throw std::runtime_error("invalid test description");
  BuildOptions o;
  o.builder = "parametric";
  auto r = build_only(d, o);
  if (dout) *dout = d;
  if (log) *log = r.log;
  return r.model;
}

static CheckReport checks_for(const Model& m, const desc::Description& d) {
  CheckContext ctx(m, d, make_well_paths(d), 100);
  return run_checks(ctx);
}

static std::string fault_yaml(const std::string& name, const std::string& age, const std::string& point, const std::string& depth,
                              const std::string& strike, const std::string& dip, const std::string& upper, const std::string& lower,
                              const std::string& length, const std::string& tmax, const std::string& pdip,
                              const std::string& pstrike, const std::string& mode, const std::string& hw,
                              const std::string& extra = "") {
  return "  - name: " + name + "\n    age: " + q(age + " ordinal") + "\n    geometry: {type: planar, reference_point: " + q(point) +
         ", reference_depth: " + q(depth) + ", strike: " + q(strike) + ", dip: " + q(dip) + "}\n    extent: {upper_tip: " + upper +
         ", lower_tip: " + lower + ", length: " + q(length) + "}\n    throw: {max: " + q(tmax) + ", profile_dip: " + pdip +
         ", profile_strike: " + pstrike + "}\n    slip_mode: " + mode + "\n    hanging_wall: " + hw + "\n" + extra;
}
static const std::string kFree = "{rule: free}";
static std::string flat_top(const std::string& taper) { return "{type: flat_topped, taper: " + q(taper + " fraction") + "}"; }

// ================================ Prompt 3: structure ======================================

TEST_CASE("P3 flat layers through the parametric engine") {
  auto m = build_text(model_yaml("", ""));
  CHECK(m->info().engine == "parametric");
  GridSpec g{0, 0, 500, 500, 21, 21};
  for (int h = 0; h < 4; ++h) {
    auto z = repr::horizon_grid_block(*m, h, 0, g);
    double worst = 0;
    for (double v : z.v) worst = std::max(worst, std::fabs(v - (-1000.0 - 200.0 * h)));
    CHECK_NEAR(worst, 0, 1e-6);
  }
  CHECK(run_model_suite(*m).passed());
}

TEST_CASE("P3 symmetric dome: known closure and area") {
  // Dome radius 2000 m, relief 100 m on every horizon except the flat base; reference well off the dome.
  std::string st = "structure:\n  regional: {type: flat}\n  base: {type: flat}\n  primitives:\n"
                   "    - {type: dome, name: D, crest: " + q("[5000, 5000] m") + ", radius_major: " + q("2000 m") +
                   ", relief: {H1: " + q("100 m") + "}}\n";
  desc::Description d;
  auto m = build_text(model_yaml(st, "", "checks:\n  closure: {horizon: H1, faults_seal: " + q("true flag") + "}\n"), &d);
  GridSpec g{0, 0, 20, 20, 501, 501};
  auto z = repr::horizon_grid_faulted(*m, 0, g);
  auto c = repr::closure_analysis(z.z, true);
  CHECK(c.found);
  CHECK_NEAR(c.closure, 100.0, 1e-6);           // crest -900 m, regional -1000 m
  CHECK_NEAR(-c.crest_z, 900.0, 1e-6);
  CHECK_NEAR(c.area, kPi * 2000.0 * 2000.0, 0.01 * kPi * 2000.0 * 2000.0);
  CHECK_NEAR(c.spill_z, -1000.0, 1e-6);
  // Relief interpolates linearly to zero at the flat base (H4): H2 has 2/3 of the relief.
  auto z2 = repr::horizon_grid_block(*m, 1, 0, GridSpec{5000, 5000, 1, 1, 1, 1});
  CHECK_NEAR(z2.v[0], -1200 + 100.0 * 2.0 / 3.0, 1e-6);
  CHECK(run_model_suite(*m).passed());
  auto rep = checks_for(*m, d);
  CHECK(rep.passed());
}

TEST_CASE("P3 asymmetric anticline: known limb dips") {
  // Cylindrical anticline, axis north, limb half-widths 1500 m (west) and 3000 m (east), relief 200 m.
  std::string st = "structure:\n  regional: {type: flat}\n  base: {type: flat}\n  primitives:\n"
                   "    - {type: anticline, name: A, crest: " + q("[5000, 5000] m") + ", axis_azimuth: " + q("0 deg") +
                   ", half_width_left: " + q("1500 m") + ", half_width_right: " + q("3000 m") + ", relief: {H1: " +
                   q("200 m") + "}}\n";
  auto m = build_text(model_yaml(st, ""));
  // Profile across the axis at 1 m spacing.
  GridSpec g{2000, 5000, 1, 1, 6001, 1};
  auto z = repr::horizon_grid_block(*m, 0, 0, g);
  double west = 0, east = 0;
  for (std::size_t i = 1; i < g.nx; ++i) {
    double slope = (z.v[i] - z.v[i - 1]) / g.dx;
    if (g.x(i) < 5000) west = std::max(west, slope); else east = std::max(east, -slope);
  }
  double dip_w = rad2deg(std::atan(west)), dip_e = rad2deg(std::atan(east));
  CHECK_NEAR(dip_w, rad2deg(std::atan(kPi * 200 / (2 * 1500))), 0.01);
  CHECK_NEAR(dip_e, rad2deg(std::atan(kPi * 200 / (2 * 3000))), 0.01);
  std::cout << "  limb dips: west " << fmt_fixed(dip_w, 3) << " deg, east " << fmt_fixed(dip_e, 3) << " deg\n";
}

TEST_CASE("P3 control wells: second well honoured by a reported correction") {
  std::string st = "structure:\n  regional: {type: flat}\n  correction_range: " + q("3000 m") + "\n  base: {type: flat}\n"
                   "  primitives:\n    - {type: dome, name: D, crest: " + q("[5000, 5000] m") + ", radius_major: " + q("3000 m") +
                   ", relief: {H1: " + q("100 m") + "}}\n";
  std::string w2 = "  - {name: W2, role: control, head: " + q("[5500, 5200] m") + ", tops: {H1: " + q("925 m") + "}}\n"
                   "  - {name: W3, role: check, head: " + q("[7000, 7000] m") + ", tops: {H1: " + q("990 m") + "}}\n";
  desc::Description d;
  BuildLog log;
  auto m = build_text(model_yaml(st, "", "", w2), &d, &log);
  auto paths = make_well_paths(d);
  auto t2 = repr::well_intersections(*m, paths[1]);
  REQUIRE(!t2.empty());
  CHECK_NEAR(t2[0].p.z, -925.0, 1e-4);
  auto t1 = repr::well_intersections(*m, paths[0]);
  CHECK_NEAR(t1[0].p.z, -1000.0, 1e-4);
  bool reported = false;
  for (auto& n : log.notes) reported = reported || n.find("correction field") != std::string::npos;
  CHECK(reported);
  auto rep = checks_for(*m, d);
  bool w3_warn = false;
  for (auto& r : rep.results)
    if (r.subject == "W3 / H1") w3_warn = r.warn_only;
  CHECK(w3_warn);  // check wells are compared, never honoured silently
}

TEST_CASE("P3 imported surface with isochore stacking") {
  fs::path dir = fs::path(stxtest::tmp_dir()) / "imported";
  fs::create_directories(dir);
  GridSpec g{0, 0, 250, 250, 41, 41};
  Grid2D top(g);
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i) top.at(i, j) = 1000 + 0.01 * g.x(i) - 0.005 * g.y(j);  // depth
  write_file(dir / "top.irap", io::write_irap_ascii(top));
  std::string text = "stratonix_schema: 1\nmodel: {name: imp, version: \"1.0.0\", status: DESIGN}\n"
                     "crs: {name: local, vertical_datum: datum, length_unit: m}\n"
                     "domain: {origin: " + q("[0, 0] m") + ", size: " + q("[10000, 10000] m") + ", depth_range: " +
                     q("[0, 3000] m") + ", grid_spacing: " + q("250 m") + "}\n"
                     "inputs: [{id: top_grid, path: top.irap, role: surface, sharing: internal, source: test}]\n"
                     "stratigraphy:\n  reference: {point: " + q("[5000, 5000] m") + "}\n"
                     "  horizons:\n    - {name: H1, surface: {input: top_grid, format: irap_ascii, z_is: depth}}\n"
                     "    - {name: H2, isochore: {from: H1, thickness: " + q("150 m") + "}}\n"
                     "  units: [{name: U1, top: H1, base: H2}]\n";
  write_file(dir / "imp.yaml", text);
  auto d = desc::load_description(dir / "imp.yaml");
  BuildOptions o;
  auto m = build_only(d, o).model;
  auto z1 = repr::horizon_grid_block(*m, 0, 0, g);
  auto z2 = repr::horizon_grid_block(*m, 1, 0, g);
  double e1 = 0, e2 = 0;
  for (std::size_t k = 0; k < g.size(); ++k) {
    e1 = std::max(e1, std::fabs(z1.v[k] + top.v[k]));
    e2 = std::max(e2, std::fabs(z1.v[k] - z2.v[k] - 150));
  }
  CHECK_NEAR(e1, 0, 1e-3);  // IRAP written with 4 decimals
  CHECK_NEAR(e2, 0, 1e-6);
}

TEST_CASE("P3 closure targets are reported, never tuned") {
  std::string st = "structure:\n  regional: {type: flat}\n  base: {type: flat}\n  primitives:\n"
                   "    - {type: dome, name: D, crest: " + q("[5000, 5000] m") + ", radius_major: " + q("2000 m") +
                   ", relief: {H1: " + q("100 m") + "}}\n";
  std::string chk = "checks:\n  closure:\n    horizon: H1\n    faults_seal: " + q("true flag") +
                    "\n    targets: {closure: " + q("120 m") + ", tolerance_fraction: " + q("0.05 fraction") + "}\n";
  desc::Description d;
  auto m = build_text(model_yaml(st, "", chk), &d);
  auto rep = checks_for(*m, d);
  bool found = false;
  for (auto& r : rep.results)
    if (r.name == "closure against target") {
      found = true;
      CHECK_NEAR(r.value, -20.0, 0.5);  // model 100 m vs target 120 m: reported as a difference
      CHECK(!r.passed);
    }
  CHECK(found);
}

// ================================ Prompt 4: faults =========================================

static double throw_of(const Model& m, int f, int h, double spacing = 50) {
  auto fs = repr::fault_surface(m, f, spacing);
  auto cl = repr::cutoff_lines(m, fs, h);
  return repr::throw_row(m, fs, cl, h).throw_max;
}

TEST_CASE("P4 one planar fault with constant throw") {
  for (std::string mode : {"dip_slip", "vertical_offset"}) {
    std::string f = fault_yaml("F", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", kFree, kFree, "100000 m", "50 m",
                               flat_top("0.25"), flat_top("0.01"), mode, "{type: rigid}");
    desc::Description d;
    auto m = build_text(model_yaml("", f), &d);
    for (int h = 0; h < 4; ++h) {
      if (mode == "dip_slip") CHECK_NEAR(throw_of(*m, 0, h), 50.0, 1e-3);
    }
    // Hanging wall (east, away from W1) is 50 m lower than the footwall everywhere.
    auto zf = repr::horizon_grid_faulted(*m, 1, GridSpec{1000, 1000, 8000, 8000, 2, 2});
    CHECK_NEAR(zf.z.at(0, 0), -1200, 1e-6);
    CHECK_NEAR(zf.z.at(1, 1), -1250, 1e-6);
    CHECK(run_model_suite(*m).passed());
    d.wells.clear();
    auto rep = checks_for(*m, d);
    for (auto& r : rep.results)
      if (!r.passed && !r.info && r.id != "C6") std::cerr << r.id << " " << r.subject << " " << r.value << "\n";
    bool ok = true;
    for (auto& r : rep.results) ok = ok && (r.passed || r.info || r.id == "C6");  // a free fault offsets top and base
    CHECK(ok);
  }
}

TEST_CASE("P4 blind fault with both tips inside units; unit holding a tip changes thickness by the throw") {
  // Tips halfway through U1 and U3: H1 and H4 unoffset; H2 and H3 offset.
  std::string f = fault_yaml("F", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg",
                             "{rule: within_unit, unit: U1, fraction: " + q("0.5 fraction") + "}",
                             "{rule: within_unit, unit: U3, fraction: " + q("0.5 fraction") + "}", "6000 m", "40 m",
                             "{type: sine}", "{type: elliptical}", "dip_slip", "{type: rigid}");
  desc::Description d;
  auto m = build_text(model_yaml("", f), &d);
  double t1 = throw_of(*m, 0, 0), t2 = throw_of(*m, 0, 1), t3 = throw_of(*m, 0, 2), t4 = throw_of(*m, 0, 3);
  std::cout << "  throws H1..H4: " << t1 << " " << t2 << " " << t3 << " " << t4 << "\n";
  CHECK_NEAR(t1, 0, 1e-6);
  CHECK_NEAR(t4, 0, 1e-6);
  CHECK(t2 > 10);
  CHECK(t3 > 10);
  // Sine profile between tips at 1100 and 1500 m depth: H2 (1200) at eta = 0.25 -> sin(pi/4) * 40 at the
  // centre; the throw is read at the restored position, so it does not depend on the frame.
  CHECK_NEAR(t2, 40 * std::sin(kPi * 0.25), 0.05);
  // U1 holds the upper tip: across the fault it thickens by the throw at H2 (its offset base).
  auto* net = static_cast<const FaultNetwork*>(m->capability("faults"));
  BlockCode hwb = (net->fixed_block() & 1u) ? 0 : 1;
  BlockCode fwb = (net->fixed_block() & 1u);
  GridSpec at{5000 + 200, 5000, 1, 1, 1, 1}, at_fw{5000 - 400, 5000, 1, 1, 1, 1};
  // compare unit thickness just inside each wall at the same strike position, at the H2 level
  auto th = [&](BlockCode b, const GridSpec& g) {
    return repr::horizon_grid_block(*m, 0, b, g).v[0] - repr::horizon_grid_block(*m, 1, b, g).v[0];
  };
  double d_thick = th(hwb, at) - th(fwb, at_fw);
  std::cout << "  U1 thickening across the fault: " << d_thick << " m (throw at H2 " << t2 << " m)\n";
  CHECK_NEAR(std::fabs(d_thick), t2, 1.0);
  CHECK(checks_for(*m, d).passed());
  CHECK(run_model_suite(*m).passed());
}

TEST_CASE("P4 two faults with different extents") {
  std::string f = fault_yaml("FA", "1", "[3500, 5000] m", "1300 m", "0 deg", "60 deg",
                             "{rule: horizon, horizon: H1, offset: " + q("50 m") + "}", "{rule: horizon, horizon: H4, offset: " + q("-50 m") + "}",
                             "5000 m", "30 m", flat_top("0.2"), flat_top("0.2"), "dip_slip", "{type: rigid}") +
                  fault_yaml("FB", "2", "[7000, 5000] m", "1300 m", "180 deg", "70 deg",
                             "{rule: horizon, horizon: H2, offset: " + q("50 m") + "}", "{rule: horizon, horizon: H3, offset: " + q("50 m") + "}",
                             "3000 m", "20 m", flat_top("0.2"), flat_top("0.2"), "dip_slip", "{type: decay, influence_distance: " + q("1500 m") + "}");
  desc::Description d;
  auto m = build_text(model_yaml("", f), &d);
  CHECK(throw_of(*m, 0, 0) == 0);
  CHECK(throw_of(*m, 0, 1) > 20);
  CHECK(throw_of(*m, 0, 2) > 20);
  CHECK(throw_of(*m, 0, 3) == 0);
  std::cout << "  FB throws H1..H4: " << throw_of(*m, 1, 0) << " " << throw_of(*m, 1, 1) << " " << throw_of(*m, 1, 2) << " "
            << throw_of(*m, 1, 3) << "\n";
  CHECK(throw_of(*m, 1, 1) == 0);   // FB: upper tip below H2
  CHECK(throw_of(*m, 1, 2) > 10);   // H3 inside FB
  CHECK(throw_of(*m, 1, 3) == 0);
  auto rep = checks_for(*m, d);
  if (!rep.passed()) std::cerr << rep.text();
  CHECK(rep.passed());
  CHECK(run_model_suite(*m).passed());
}

TEST_CASE("P4 two intersecting faults with abutting") {
  // FA strikes north; FB strikes east, is younger and abuts FA on FA's footwall side (west).
  std::string f = fault_yaml("FA", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", kFree, kFree, "9000 m", "40 m",
                             flat_top("0.25"), flat_top("0.2"), "dip_slip", "{type: rigid}") +
                  fault_yaml("FB", "2", "[3000, 5000] m", "1300 m", "90 deg", "60 deg", kFree, kFree, "6000 m", "30 m",
                             flat_top("0.25"), flat_top("0.2"), "dip_slip", "{type: rigid}",
                             "    abuts: {fault: FA, side: footwall}\n");
  desc::Description d;
  auto m = build_text(model_yaml("", f, "", ""), &d);
  // FB displaces nothing east of FA (FA's hanging wall).
  std::vector<Vec3> east;
  for (int k = 0; k < 50; ++k) east.push_back({6000.0 + 50 * k, 4000.0, -1300.0});
  FaultEval e;
  m->evaluate_fault(1, east, e);
  double leak = 0;
  for (double v : e.displacement) leak = std::max(leak, v);
  CHECK_NEAR(leak, 0, 1e-12);
  CHECK(throw_of(*m, 1, 1) > 20);
  auto rep = checks_for(*m, d);
  bool c8 = false;
  for (auto& r : rep.results)
    if (r.id == "C8") { c8 = true; CHECK(r.passed); }
  CHECK(c8);
  CHECK(run_model_suite(*m).passed());
}

TEST_CASE("P4 fault tip inside a unit: offsets continuous where the fault dies out") {
  std::string f = fault_yaml("F", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg",
                             "{rule: within_unit, unit: U2, fraction: " + q("0.3 fraction") + "}",
                             "{rule: horizon, horizon: H4, offset: " + q("-20 m") + "}", "8000 m", "30 m",
                             flat_top("0.3"), flat_top("0.3"), "dip_slip", "{type: rigid}");
  desc::Description d;
  auto m = build_text(model_yaml("", f), &d);
  CHECK_NEAR(throw_of(*m, 0, 1), 0, 1e-9);  // H2 above the tip
  CHECK(throw_of(*m, 0, 2) > 5);            // H3 below
  // Continuity: the throw of H3 along strike decreases smoothly to zero at the ends.
  auto fs = repr::fault_surface(*m, 0, 50);
  auto cl = repr::cutoff_lines(*m, fs, 2);
  double jump = 0, prev = kNaN;
  for (std::size_t i = 0; i < fs.na; ++i) {
    auto a = cl[0].by_column[i], b = cl[1].by_column[i];
    double t = (std::isnan(a.x) || std::isnan(b.x)) ? 0 : b.z - a.z;
    if (!std::isnan(prev)) jump = std::max(jump, std::fabs(t - prev));
    prev = t;
  }
  CHECK(jump < 2.0);  // no step larger than 2 m between 50 m strike samples
  auto rep = checks_for(*m, d);
  if (!rep.passed()) std::cerr << rep.text();
  CHECK(rep.passed());
}

TEST_CASE("P4 dip-slip versus vertical offset on the synthetic example") {
  register_all_plugins();
  auto d1 = desc::load_description(std::string(stxtest::source_dir()) + "/examples/synthetic_fold_two_faults.yaml");
  auto d2 = d1;
  for (auto& f : d2.faults) f.slip_mode = "vertical_offset";
  BuildOptions o;
  auto m1 = build_only(d1, o).model, m2 = build_only(d2, o).model;
  GridSpec g = GridSpec::covering(d1.domain(), 100);
  double worst = 0;
  std::string where;
  for (int h = 0; h < 5; ++h) {
    auto a = repr::horizon_grid_faulted(*m1, h, g), b = repr::horizon_grid_faulted(*m2, h, g);
    for (std::size_t k = 0; k < g.size(); ++k)
      if (!std::isnan(a.z.v[k]) && !std::isnan(b.z.v[k]) && std::fabs(a.z.v[k] - b.z.v[k]) > worst) {
        worst = std::fabs(a.z.v[k] - b.z.v[k]);
        where = d1.horizons[static_cast<std::size_t>(h)].name + " at (" + fmt_fixed(g.x(k % g.nx), 0) + ", " + fmt_fixed(g.y(k / g.nx), 0) + ")";
      }
  }
  std::cout << "  dip-slip vs vertical offset: max horizon difference " << fmt_fixed(worst, 2) << " m, " << where << "\n";
  CHECK(worst > 0);
  auto rep = checks_for(*m2, d2);
  if (!rep.passed()) std::cerr << rep.text();
  CHECK(rep.passed());
}

STX_TEST_MAIN()
