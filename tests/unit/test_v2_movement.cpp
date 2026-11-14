// Stage 2A, Prompt 18: reverse, strike-slip and oblique movement (STX-TN-05 section 5).
#include "v2_fixture.hpp"
#include "stratonix/io/io.hpp"

using namespace v2fx;

namespace {
// Layers dipping 5 degrees to the north: depth increases with y.
const std::string kNorthDip = "structure:\n  regional: {type: planar, dip: " + q("5 deg") + ", dip_azimuth: " + q("0 deg") + "}\n";
}  // namespace

TEST_CASE("P18 reverse fault: repeated section equals the throw") {
  std::string f = planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "30 deg", "100000 m",
                               "max: " + q("150 m") + ", sense: reverse", "dip_slip");
  auto m = build_text(model_yaml(f));
  // At x = 5300 the fault is at 1473 m: the raised hanging wall holds H1 to H4 above it, the footwall
  // H4 below it.
  double hw_h1 = depth_in_block(*m, 0, 1, 5300, 5000), fw_h4 = depth_in_block(*m, 3, 0, 5300, 5000);
  CHECK_NEAR(hw_h1, 850, 0.5);
  CHECK_NEAR(fw_h4 - hw_h1 - 600, 150, 0.5);
  // Down the vertical line the unit sequence repeats U3 (index 3) below the fault.
  std::vector<int> seq;
  for (double z = 800; z <= 1700; z += 0.5) {
    int u = m->unit_at({5300, 5000, -z});
    if (seq.empty() || seq.back() != u) seq.push_back(u);
  }
  int u3 = 0;
  for (int u : seq) u3 += u == 3;
  std::cout << "  unit sequence:";
  for (int u : seq) std::cout << " " << u;
  std::cout << "\n";
  CHECK(u3 == 2);
}

TEST_CASE("P18 right-lateral strike-slip offsets dipping layers by slip x tan(dip)") {
  std::string f = planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "90 deg", "100000 m",
                               "slip: " + q("400 m") + ", sense: right_lateral", "strike_slip");
  auto m = build_text(model_yaml(f, kNorthDip));
  // East block moves south by 400 m relative to the west (fixed) block.
  double step = depth_faulted(*m, 0, 6000, 5000) - depth_faulted(*m, 0, 4000, 5000);
  CHECK_NEAR(step, 400 * std::tan(deg2rad(5)), 0.5);
}

TEST_CASE("P18 oblique slip at rake -45 degrees: dip and strike components") {
  std::string f = planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", "slip: " + q("300 m"),
                               "oblique", "{type: rigid}", "    rake: " + q("-45 deg") + "\n");
  auto m = build_text(model_yaml(f, kNorthDip));
  // Hanging wall moves down 300 sin45 sin60 and north 300 cos45 (towards shallower layers).
  double expected = 300 * std::sin(deg2rad(45)) * std::sin(deg2rad(60)) - 300 * std::cos(deg2rad(45)) * std::tan(deg2rad(5));
  double step = depth_faulted(*m, 0, 7000, 5000) - depth_faulted(*m, 0, 3000, 5000);
  CHECK_NEAR(step, expected, 0.5);
}

TEST_CASE("P18 thrust exports: overlapping hanging-wall and footwall grids load back") {
  namespace fs = std::filesystem;
  register_all_plugins();
  fs::path out = fs::path(stxtest::tmp_dir()) / "v2_thrust_export";
  fs::remove_all(out);
  BuildOptions o;
  o.figures = false;
  o.exports = {"petrel"};
  auto r = run_build(fs::path(stxtest::source_dir()) / "examples/thrust_repeated_section.yaml", out, o);
  REQUIRE(r.checks_passed);
  // Between the hanging-wall and footwall cut-offs of H2 (x 4567 to 4827 m at y = 5000 m) both
  // sheets exist: the raised hanging wall at 1050 m and the footwall at 1200 m.
  auto hw = io::read_irap_ascii(read_file(out / "exports/horizons/per_block/H2_block1.irap"));
  auto fw = io::read_irap_ascii(read_file(out / "exports/horizons/per_block/H2_block0.irap"));
  CHECK_NEAR(-hw.sample(4700, 5000), 1050, 0.5);
  CHECK_NEAR(-fw.sample(4700, 5000), 1200, 0.5);
}

STX_TEST_MAIN()
