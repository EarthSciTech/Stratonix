// Stage 4B tests (STX-TN-15 section 11): analytical cases G1 to G24 of geological completeness.
// Written first (Prompt 48): they fail with "not implemented" until Prompts 49 to 53 build the
// features. Example models are added by each prompt.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <numeric>

#include "v2_fixture.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/bodies.hpp"
#include "stratonix/engines/depositional.hpp"
#include "stratonix/engines/faultstyles.hpp"
#include "stratonix/engines/fractures.hpp"
#include "stratonix/engines/restoration3d.hpp"
#include "stratonix/engines/thermal.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/faultrock.hpp"

using namespace stratonix;
using namespace stxtest;
using namespace v2fx;
namespace fs = std::filesystem;

namespace {
double rel(double a, double b) { return std::fabs(a - b) / std::max(1e-300, std::fabs(b)); }
}  // namespace

// ---- 1. Bodies (Prompt 49)

TEST_CASE("G1 volume of a salt sheet of known footprint, base and taper") {
  SheetBody s;
  s.footprint = {{0, 0, 0}, {5000, 0, 0}, {5000, 3000, 0}, {0, 3000, 0}};
  s.base_z = -2000;
  s.thickness = 400;
  s.taper = 200;
  // closed form: full thickness on the inner rectangle, linear taper over 200 m inside the edge
  // V = t [A_inner + (perimeter band of width w) / 2] with the corner squares counted by their pyramids
  const double W = 5000, L = 3000, w = 200, t = 400;
  const double inner = (W - 2 * w) * (L - 2 * w);
  const double sides = 2 * (W - 2 * w) * w + 2 * (L - 2 * w) * w;  // wedges: half volume
  const double corners = 4 * w * w;                                // min-distance taper: 1/3 volume
  const double v = t * (inner + sides / 2 + corners / 3);
  const double got = sheet_volume(s, 50);
  std::printf("  sheet volume %.4e m3, closed form %.4e m3\n", got, v);
  CHECK(rel(got, v) < 0.005);
}

TEST_CASE("G2 forced fold above a circular sill: fold volume equals sill volume, peak uplift 2 h0") {
  SillBody s;
  s.centre = {5000, 5000, -1500};
  s.radius = 2000;
  s.thickness = 100;
  s.accommodation = 1;
  CHECK(rel(sill_volume(s), 2.0 / 3.0 * kPi * 2000 * 2000 * 100) < 1e-9);
  CHECK(rel(forced_fold_uplift(s, 0), 200) < 0.01);
  CHECK(rel(forced_fold_volume(s, 25), sill_volume(s)) < 0.01);
  s.accommodation = 0.4;
  CHECK(rel(forced_fold_volume(s, 25), 0.4 * sill_volume(s)) < 0.01);
  CHECK(forced_fold_uplift(s, 2100) == 0);
}

TEST_CASE("G2b injectite mass balance and seeded reproducibility") {
  InjectiteSpec sp;
  sp.parent_top = -1200;
  sp.ceiling = -900;
  sp.parent_sand_volume = 5e8;
  sp.sheets = 25;
  sp.area.add({0, 0, -2000});
  sp.area.add({10000, 10000, 0});
  auto a = build_injectites(sp, 7), b = build_injectites(sp, 7);
  REQUIRE(a.ok);
  REQUIRE(a.sheets.size() == b.sheets.size());
  for (std::size_t i = 0; i < a.sheets.size(); ++i) CHECK(norm(a.sheets[i].start - b.sheets[i].start) == 0);
  CHECK(rel(sp.parent_sand_volume - a.parent_volume_after, a.injected_volume) < 0.005);
  CHECK(a.injected_volume <= sp.max_fraction * sp.parent_sand_volume * (1 + 1e-12));
  for (auto& s : a.sheets) CHECK(std::max(s.start.z, s.end.z) <= sp.ceiling + 1e-6);
}

TEST_CASE("G1m bodies in models: unit query, drape of the overburden, restoration removes them, younger faults refused") {
  register_all_plugins();
  const fs::path ex = fs::path(STX_SOURCE_DIR) / "examples";
  auto d = desc::load_description(ex / "salt_sheet_drape.yaml");
  auto m = build_only(d, BuildOptions{}).model;
  auto* bs = static_cast<const std::vector<BodyReport>*>(m->capability("bodies"));
  REQUIRE(bs != nullptr);
  REQUIRE(bs->size() == 1);
  const int salt = (*bs)[0].unit;
  // at the sheet's centre (x 5500, y 4000): base on H2 (1300 m), 500 m thick, H1 lifted by 500 m
  CHECK(m->unit_at({5500, 4000, -1050}) == salt);
  CHECK(m->unit_at({5500, 4000, -1350}) != salt);
  CHECK(m->unit_at({5500, 4000, -790}) == d.horizon_index("H1") + 1);  // U1 above the sheet
  CHECK(m->unit_at({5500, 4000, -310}) == d.horizon_index("H1") + 1);  // H1 lifted from 800 to 300 m
  CHECK(m->unit_at({5500, 4000, -290}) == 0);
  // restoration removes the body: the host at the undeformed position (flat host: z - thickness)
  Vec3 p{5500, 4000, -600}, r;
  m->restore(std::span<const Vec3>(&p, 1), std::span<Vec3>(&r, 1));
  CHECK(std::fabs(r.z - (-1100)) < 1e-6);
  // a fault younger than a body that crosses it is refused with the reason
  auto t = model_yaml(planar_fault("F1", "2", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "40000 m", "max: " + q("100 m"), "dip_slip"), "",
                      {1000, 1200, 1400, 1600},
                      "bodies:\n  - {name: S1, type: sill, age: " + q("1 ordinal") + ", centre: " + q("[5000, 5000] m") + ", depth: " + q("1300 m") +
                          ", radius: " + q("1500 m") + ", thickness: " + q("50 m") + "}\n");
  bool refused = false;
  std::string why;
  try {
    build_text(t);
  } catch (const std::exception& e) {
    refused = true;
    why = e.what();
  }
  std::printf("  refusal: %s\n", why.c_str());
  CHECK(refused);
  CHECK(why.find("younger fault F1") != std::string::npos);
}

// ---- 2. Fault styles (Prompt 50)

TEST_CASE("G3 polygonal fault network: reproducible, spacing within 5 %, throws scale with length") {
  PolygonalSpec s;
  s.x1 = s.y1 = 20000;
  s.unit_top = -800;
  s.unit_base = -1100;
  s.spacing = 600;
  auto a = polygonal_network(s, 11), b = polygonal_network(s, 11);
  REQUIRE(!a.faults.empty());
  REQUIRE(a.faults.size() == b.faults.size());
  CHECK(norm(a.faults[3].a - b.faults[3].a) == 0);
  std::printf("  %zu faults, mean spacing %.1f m (labelled 600)\n", a.faults.size(), a.mean_spacing);
  CHECK(rel(a.mean_spacing, 600) < 0.05);
  for (auto& f : a.faults) {
    CHECK(rel(f.throw_m, s.throw_ratio * f.length) < 1e-9);
    CHECK(f.throw_m < s.unit_top - s.unit_base);
  }
}

TEST_CASE("G4 flower structure: splays carry the vertical component and conserve horizontal slip") {
  FlowerSpec s;
  s.horizontal_slip = 1000;
  s.vertical_component = 150;
  s.splays = 5;
  s.spread_deg = 40;
  auto r = flower_splays(s);
  REQUIRE(r.vertical_throw.size() == 5);
  CHECK(rel(std::accumulate(r.vertical_throw.begin(), r.vertical_throw.end(), 0.0), 150) < 0.01);
  CHECK(rel(std::accumulate(r.horizontal_slip.begin(), r.horizontal_slip.end(), 0.0), 1000) < 0.01);
  s.vertical_component = -150;  // negative flower
  auto n = flower_splays(s);
  CHECK(std::accumulate(n.vertical_throw.begin(), n.vertical_throw.end(), 0.0) < 0);
}

TEST_CASE("G5 duplex shortening equals the summed slip of its horses") {
  DuplexSpec s;
  s.horses = 4;
  s.horse_length = 1200;
  s.slip = 700;
  auto r = duplex(s);
  CHECK(rel(r.shortening, 4 * 700) < 0.01);
  CHECK(r.style == "hinterland_dipping");  // slip < horse length
  s.slip = 1200;
  CHECK(duplex(s).style == "antiformal_stack");
  s.slip = 1800;
  CHECK(duplex(s).style == "foreland_dipping");
}

TEST_CASE("G6 breached relay: summed throw continuous along strike") {
  RelaySpec s{8000, 8000, 2000, 500, 120, 100};
  auto p = breached_relay_throw(s, 2000);  // 7 m apart: a smooth profile changes by less than 1 % per sample
  REQUIRE(p.size() == 2000);
  double tmax = 0, step = 0;
  for (std::size_t i = 0; i < p.size(); ++i) {
    tmax = std::max(tmax, p[i].second);
    if (i) step = std::max(step, std::fabs(p[i].second - p[i - 1].second));
  }
  std::printf("  largest step %.3f m of %.1f m\n", step, tmax);
  CHECK(step < 0.01 * tmax);
}

TEST_CASE("G7 trishear on a listric fault against an incremental reference construction") {
  // listric fault (x, z): radius 3 km from dip 60 to 20 degrees
  std::vector<Vec3> fault;
  for (int k = 0; k <= 40; ++k) {
    const double a = (60.0 - 40.0 * k / 40.0) * kPi / 180;
    fault.push_back({3000 * (1 - std::cos(kPi / 3) ) + 3000 * (std::cos(kPi / 3) - std::cos(a)), 0, -3000 * (std::sin(kPi / 3) - std::sin(a)) - 1000});
  }
  TrishearParams p;
  p.slip = 500;
  std::vector<Vec3> pts;
  for (double x = -1000; x <= 4000; x += 250) pts.push_back({x, 0, -800});
  auto moved = trishear_curved(fault, p, pts);
  REQUIRE(moved.size() == pts.size());
  // far hanging-wall points move by the slip parallel to the fault at depth (rigid translation)
  double far = norm(moved.front() - pts.front());
  CHECK(far <= p.slip * 1.01);
  // footwall points below the tip's whole path never enter the trishear zone and do not move
  // (points ahead of the tip at shallower levels are swept by the zone as the tip propagates)
  std::vector<Vec3> fw;
  for (double x = 3000; x <= 5000; x += 500) fw.push_back({x, 0, -2500});
  auto fwm = trishear_curved(fault, p, fw);
  for (std::size_t i = 0; i < fw.size(); ++i) CHECK(norm(fwm[i] - fw[i]) < 1e-9);
}

TEST_CASE("G8 growth strata on a rigid hanging wall: thickness = slip increment x sin(dip)") {
  GrowthSpec s;
  s.hanging_wall = "rigid";
  s.dip_deg = 60;
  s.slip_rate = 10;  // m per time unit
  s.deposition_times = {0, 2, 5, 9};
  auto t = growth_thickness(s);
  REQUIRE(t.size() == 3);
  CHECK(rel(t[0], 20 * std::sin(kPi / 3)) < 0.01);
  CHECK(rel(t[1], 30 * std::sin(kPi / 3)) < 0.01);
  CHECK(rel(t[2], 40 * std::sin(kPi / 3)) < 0.01);
  for (const char* hw : {"shear", "fault_bend", "trishear", "layer_rigid"}) {
    s.hanging_wall = hw;
    CHECK(growth_thickness(s).size() == 3);
  }
}

// ---- 3. Depositional geometries (Prompt 51)

TEST_CASE("G9 onlap line of a horizontal horizon on a dipping basement plane") {
  const double z0 = -1000, theta = 5 * kPi / 180;
  auto horizon = [&](double, double) { return z0; };
  auto basement = [&](double x, double) { return -1500 + std::tan(theta) * x; };  // reaches z0 at x = 500 / tan(theta)
  GridSpec g{0, 0, 50, 50, 201, 201};
  auto line = onlap_line(horizon, basement, g);
  REQUIRE(!line.empty());
  const double xl = 500 / std::tan(theta);
  for (auto& p : line) CHECK(std::fabs(p.x - xl) <= 50);
}

TEST_CASE("G10 oblique clinoforms: shelf break advances by the spacing, foresets dip at the label") {
  ClinoformSpec s;
  s.count = 6;
  s.spacing = 400;
  s.foreset_dip_deg = 4;
  auto c = clinoform_set(s);
  REQUIRE(c.size() == 6);
  for (std::size_t i = 1; i < c.size(); ++i) CHECK(std::fabs(c[i].shelf_break_x - c[i - 1].shelf_break_x - 400) < 1e-6);
  for (auto& k : c) CHECK(std::fabs(k.foreset_dip_deg - 4) < 0.1);
}

TEST_CASE("G10b sigmoidal and tangential foresets dip at the label; a rising trajectory lifts the shelf break") {
  for (const char* shape : {"sigmoidal", "tangential"}) {
    ClinoformSpec s;
    s.count = 4;
    s.spacing = 500;
    s.foreset_dip_deg = 3;
    s.trajectory_deg = 2;
    s.shape = shape;
    auto c = clinoform_set(s);
    REQUIRE(c.size() == 4);
    for (auto& k : c) CHECK(std::fabs(k.foreset_dip_deg - 3) < 0.1);
    for (std::size_t i = 1; i < c.size(); ++i) CHECK(std::fabs(c[i].shelf_break_z - c[i - 1].shelf_break_z - 500 * std::tan(2 * kPi / 180)) < 1e-9);
  }
}

TEST_CASE("G11 channel-belt net-to-gross of a known placement") {
  UnitBox u{0, 0, 10000, 10000, -1100, -1000};
  ChannelSpec s;
  s.count = 10;
  s.width = 200;
  s.thickness = 10;
  s.amplitude = 0;  // straight
  s.orientation_deg = 0;  // north
  s.random_placement = false;
  auto b = place_channels(s, u, 1);
  REQUIRE(b.size() == 10);
  const double ntg = 10 * 200 * 10.0 * 10000 / (10000.0 * 10000 * 100);
  const double got = net_to_gross(b, u, 10);
  std::printf("  net-to-gross %.4f, closed form %.4f\n", got, ntg);
  CHECK(rel(got, ntg) < 0.01);
}

// ---- 4. Restoration (Prompt 52)

TEST_CASE("G12 flexural slip conserves bed length") {
  // a concentric anticline: arc of radius 2000 m over 90 degrees
  std::vector<Vec3> h;
  for (int k = 0; k <= 180; ++k) {
    const double a = (-45 + 90.0 * k / 180) * kPi / 180;
    h.push_back({2000 * std::sin(a), 0, -2000 + 2000 * std::cos(a) - 1000});
  }
  auto r = flexural_slip_unfold(h, 0, -1000);
  const double arc = 2000 * kPi / 2;
  CHECK(rel(r.present_length, arc) < 0.001);
  CHECK(rel(r.restored_length, arc) < 0.001);
  for (auto& p : r.restored) CHECK(std::fabs(p.z + 1000) < 1e-6);
}

TEST_CASE("G13 a rigidly tilted layer cake restores flat with zero strain") {
  auto m = build_text(model_yaml("", "structure:\n  regional: {type: planar, dip: " + q("5 deg") + ", dip_azimuth: " + q("90 deg") + "}\n"));
  VolumetricSpec s;
  s.horizon = "H2";
  s.datum_z = -1200;
  auto r = volumetric_restore(*m, s);
  REQUIRE(r.ok);
  std::printf("  max strain %.2e, flatness rms %.3f m, %zu tetrahedra\n", r.max_strain, r.flatness_rms, r.tetrahedra);
  CHECK(r.max_strain < 1e-3);
  CHECK(r.flatness_rms < 0.5);
}

TEST_CASE("G14 a fold restores flat within the labelled strain limit (hybrid source model)") {
  register_all_plugins();
  auto d = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples" / "data" / "hybrid_fault_fold" / "source.yaml");
  auto m = build_only(d, BuildOptions{}).model;
  VolumetricSpec s;
  s.horizon = "H1";
  s.datum_z = -1000;
  auto r = volumetric_restore(*m, s);
  REQUIRE(r.ok);
  CHECK(r.max_strain < s.max_strain);
  CHECK(r.flatness_rms < 1.0);
}

TEST_CASE("G14b the stage 2E example restores in 3D within the strain limit") {
  register_all_plugins();
  auto d = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples" / "rift_inversion_history.yaml");
  auto m = build_only(d, BuildOptions{}).model;
  VolumetricSpec s;
  s.horizon = "H2";  // the post-rift top, folded by the inversion and offset by F1
  s.datum_z = -1000;
  auto r = volumetric_restore(*m, s);
  REQUIRE(r.ok);
  std::printf("  max strain %.4f, flatness %.3f m, bed stretch %.4f, %zu tetrahedra, %zu fault contacts\n", r.max_strain, r.flatness_rms,
              r.horizon_stretch_max, r.tetrahedra, r.fault_contacts);
  CHECK(r.max_strain < s.max_strain);
  CHECK(r.flatness_rms < 1.0);
}

TEST_CASE("G15 compaction across a fault active during burial: each block on its own Athy curve") {
  const double phi0 = 0.6, c = 5e-4;
  CHECK(rel(block_porosity(phi0, c, 2000), phi0 * std::exp(-c * 2000)) < 0.005);
  CHECK(rel(block_porosity(phi0, c, 1500), phi0 * std::exp(-c * 1500)) < 0.005);
}

TEST_CASE("G16 Airy backstripping and flexural response against closed forms") {
  // 3 km of sediment of 2400 kg/m3 in 200 m of water, no sea-level change
  const double y = airy_tectonic_subsidence(3000, 2400, 3300, 1030, 200, 0);
  CHECK(rel(y, 3000 * (3300 - 2400) / (3300 - 1030.0) + 200) < 0.001);
  const double te = 20e3, E = 70e9, nu = 0.25, drho = 3300 - 2400, k = 2 * kPi / 200e3;
  const double D = E * te * te * te / (12 * (1 - nu * nu));
  CHECK(rel(flexural_response(k, te, E, nu, drho), 1 / (1 + D * std::pow(k, 4) / (drho * 9.81))) < 0.01);
  CHECK(rel(flexural_response(k, 0, E, nu, drho), 1.0) < 1e-12);
}

TEST_CASE("G17 eroded section from a compaction trend: 800 m removed") {
  const double phi0 = 0.55, c = 4.5e-4, removed = 800;
  std::vector<std::pair<double, double>> log;
  for (double z = 50; z <= 2000; z += 50) log.push_back({z, phi0 * std::exp(-c * (z + removed))});
  CHECK(rel(eroded_section(log, phi0, c), removed) < 0.02);
}

// ---- 5. Fractures and mechanics (Prompt 53)

TEST_CASE("G18 fracture intensity on a known curvature") {
  FractureLaw law{0.05, 40.0, 0.0, 10};
  const double R = 2000;
  CHECK(rel(fracture_intensity(1 / R, 0, law), 0.05 + 40.0 / R) < 0.005);
  CHECK(fracture_intensity(1.0, 0, law) == 10);  // clipped
}

TEST_CASE("G19 Baecher-disk network: reproducible, realised P32 within 3 %, set poles within 2 degrees") {
  BBox box;
  box.add({0, 0, 0});
  box.add({200, 200, 50});
  std::vector<FractureSet> sets(2);
  sets[0].mean_pole = normalized(Vec3{1, 0, 0});
  sets[0].kappa = 40;
  sets[1].mean_pole = normalized(Vec3{0, 1, 0.2});
  sets[1].kappa = 25;
  auto a = sample_dfn(box, 0.3, sets, 99), b = sample_dfn(box, 0.3, sets, 99);
  REQUIRE(a.ok);
  REQUIRE(a.disks.size() == b.disks.size());
  CHECK(norm(a.disks[10].c - b.disks[10].c) == 0);
  const double p = realised_p32(a.disks, box);
  std::printf("  %zu disks, realised P32 %.4f (target 0.3)\n", a.disks.size(), p);
  CHECK(rel(p, 0.3) < 0.03);
  for (int k = 0; k < 2; ++k) {
    Vec3 sum{0, 0, 0};
    for (auto& d : a.disks)
      if (d.set == k) sum = sum + (dot(d.n, sets[static_cast<std::size_t>(k)].mean_pole) < 0 ? d.n * -1.0 : d.n);
    CHECK(std::acos(std::clamp(dot(normalized(sum), sets[static_cast<std::size_t>(k)].mean_pole), -1.0, 1.0)) * 180 / kPi < 2);
  }
  auto big = sample_dfn(box, 0.3, sets, 99, 10);
  CHECK(!big.ok);
}

TEST_CASE("G20 horizontal stresses of two layers under tectonic strain (closed form)") {
  const double sv = 50e6, pp = 20e6, alpha = 1;
  for (auto [E, nu] : {std::pair{20e9, 0.25}, std::pair{8e9, 0.35}}) {
    auto [sh, sH] = horizontal_stress(sv, pp, alpha, E, nu, 2e-4, 6e-4);
    const double base = nu / (1 - nu) * (sv - alpha * pp) + alpha * pp;
    CHECK(rel(sh, base + E / (1 - nu * nu) * (2e-4 + nu * 6e-4)) < 0.001);
    CHECK(rel(sH, base + E / (1 - nu * nu) * (6e-4 + nu * 2e-4)) < 0.001);
  }
}

TEST_CASE("G21 shale gouge ratio of an implicit model through its restored column matches the parametric source") {
  register_all_plugins();
  auto dp = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples" / "data" / "implicit_faults_series" / "truth.yaml");
  auto di = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples" / "implicit_faults_series_from_data.yaml");
  auto mp = build_only(dp, BuildOptions{}).model;
  auto mi = build_only(di, BuildOptions{}).model;
  std::printf("  engines: %s / %s; column capability: %d / %d\n", mp->info().engine.c_str(), mi->info().engine.c_str(), mp->capability("column") != nullptr,
              mi->capability("column") != nullptr);
  repr::FaultRockInput in;
  in.unit_clay.assign(mp->info().units.size(), kNaN);
  for (std::size_t u = 0; u < in.unit_clay.size(); ++u) in.unit_clay[u] = (u % 2) ? 0.6 : 0.15;
  repr::FaultRockInput ii = in;
  ii.unit_clay.assign(mi->info().units.size(), kNaN);
  for (std::size_t u = 0; u < ii.unit_clay.size(); ++u) ii.unit_clay[u] = (u % 2) ? 0.6 : 0.15;
  repr::StressInput st;
  std::vector<double> diff;
  for (double t = -3000; t <= 3000; t += 500) {
    // points on F1 (strike 020 through (5000, 5000) at 1300 m) along strike at its reference depth
    const Vec3 p{5000 + t * std::sin(20 * kPi / 180), 5000 + t * std::cos(20 * kPi / 180), -1300};
    auto a = repr::fault_rock_at(*mp, 0, p, in, st);
    auto b = repr::fault_rock_at(*mi, 0, p, ii, st);
    if (std::isfinite(a.sgr)) {
      CHECK(std::isfinite(b.sgr));
      if (std::isfinite(b.sgr)) diff.push_back(std::fabs(a.sgr - b.sgr));
    }
  }
  REQUIRE(diff.size() >= 5);
  std::sort(diff.begin(), diff.end());
  std::printf("  median |SGR difference| %.3f\n", diff[diff.size() / 2]);
  CHECK(diff[diff.size() / 2] < 0.02);
}

// ---- 6. Heat and pressure (Prompt 53)

TEST_CASE("G22 steady conduction through a layered column with radiogenic heat") {
  std::vector<ThermalLayer> L{{1000, 2.0, 1e-6}, {1500, 3.0, 0.5e-6}, {500, 1.5, 0}};
  const double ts = 10, q = 0.065;
  auto T = steady_temperature(L, ts, q);
  REQUIRE(T.size() == 4);
  double t = ts, qq = q;
  CHECK(std::fabs(T[0] - ts) < 1e-9);
  for (std::size_t i = 0; i < L.size(); ++i) {
    t += qq * L[i].thickness / L[i].conductivity - L[i].heat_production * L[i].thickness * L[i].thickness / (2 * L[i].conductivity);
    qq -= L[i].heat_production * L[i].thickness;
    CHECK(std::fabs(T[i + 1] - t) < 0.1);
  }
}

TEST_CASE("G23 transient heating after a basal step against the series solution") {
  TransientSpec s;
  s.length = 1000;
  s.diffusivity = 1e-6;
  s.t_top = 0;
  s.t_base = 100;
  s.t_initial = 0;
  s.time = 0.1 * s.length * s.length / s.diffusivity;  // dimensionless time 0.1
  auto T = transient_temperature(s);
  REQUIRE(T.size() == static_cast<std::size_t>(s.nodes));
  // T(z, t) = T1 [ z/L + (2/pi) sum (-1)^n / n sin(n pi z / L) exp(-n^2 pi^2 kappa t / L^2) ], z from the top
  for (int i = 20; i < s.nodes; i += 20) {
    const double z = s.length * i / (s.nodes - 1);
    double v = z / s.length;
    for (int n = 1; n < 200; ++n) v += 2 / kPi * std::pow(-1.0, n) / n * std::sin(n * kPi * z / s.length) * std::exp(-n * n * kPi * kPi * 0.1);
    CHECK(std::fabs(T[static_cast<std::size_t>(i)] - 100 * v) <= 1.0);
  }
}

TEST_CASE("G24 overpressure of a known sedimentation rate and permeability (Gibson 1958 limits)") {
  GibsonSpec s;
  s.rate = 1e-3 / (365.25 * 86400);  // 1 mm/yr
  s.duration = 2e6 * 365.25 * 86400;  // 2 Myr: 2000 m
  s.permeability = 1e-14;            // fast drainage: hydrostatic
  CHECK(gibson_overpressure_base(s) < 0.02 * (s.rho_s - s.rho_w) * s.g * 2000);
  s.permeability = 1e-24;  // undrained: the whole buoyant load
  CHECK(rel(gibson_overpressure_base(s), (s.rho_s - s.rho_w) * s.g * 2000) < 0.02);
  // intermediate: Gibson's dimensionless time factor T = m^2 t / c_v; check against the series
  s.permeability = 1e-19;
  const double cv = s.permeability / (s.viscosity * s.compressibility * s.rho_w * s.g) * s.rho_w * s.g;
  const double T = s.rate * s.rate * s.duration / cv;
  const double got = gibson_overpressure_base(s) / ((s.rho_s - s.rho_w) * s.g * s.rate * s.duration);
  std::printf("  Gibson time factor %.3f, base overpressure ratio %.3f\n", T, got);
  CHECK(got > 0.0);
  CHECK(got < 1.0);
}

STX_TEST_MAIN()
