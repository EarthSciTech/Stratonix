// Layer-rigid hanging walls, tips on the unfaulted stratigraphy, completions, unit labels,
// outer units and section sheets (key numbers and the figure).
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "../stx_test.hpp"
#include "stratonix/checks/checks.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/faults.hpp"
#include "stratonix/engines/structure.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/figures/figures.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

using namespace stratonix;
namespace fs = std::filesystem;

namespace {
const fs::path kExample = fs::path(STX_SOURCE_DIR) / "examples/layered_faults_sheet.yaml";

std::string slurp(const fs::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}
std::string replace_all(std::string s, const std::string& a, const std::string& b) {
  for (std::size_t p = 0; (p = s.find(a, p)) != std::string::npos; p += b.size()) s.replace(p, a.size(), b);
  return s;
}
/// Present elevation of horizon h at (x, y) in block b (root of the block's scalar).
double horizon_in_block(const Model& m, int h, BlockCode b, double x, double y) {
  GridSpec g;
  g.x0 = x; g.y0 = y; g.dx = 1; g.dy = 1; g.nx = 1; g.ny = 1;
  return repr::horizon_grid_block(m, h, b, g).v[0];
}
}  // namespace

TEST_CASE("L1 layer-rigid hanging walls move every layer as a rigid slab") {
  register_all_plugins();
  auto d = desc::load_description(kExample);
  BuildOptions opt;
  auto o = build_only(d, opt);
  const Model& m = *o.model;
  auto* net = static_cast<const FaultNetwork*>(m.capability("faults"));
  auto* st = static_cast<const Structure*>(m.capability("structure"));
  REQUIRE(net && st);
  CHECK(net->layered());
  const int sand = d.horizon_index("Top_Sand");
  // South block (hanging wall of FA, footwall of FB): the same displacement 1 km and 3 km away
  // from the faults, although the horizon lies at different depths there.
  Vec3 p1{4000, 1000, -1300}, p2{4000, 3000, -1300};
  EvalResult r;
  std::vector<Vec3> pts{p1, p2};
  m.evaluate(pts, kEvalBlock, r);
  CHECK(r.block[0] == r.block[1]);
  double z1 = horizon_in_block(m, sand, r.block[0], p1.x, p1.y), z2 = horizon_in_block(m, sand, r.block[0], p2.x, p2.y);
  double d1 = z1 - st->z(static_cast<std::size_t>(sand), p1.x, p1.y), d2 = z2 - st->z(static_cast<std::size_t>(sand), p2.x, p2.y);
  CHECK(std::fabs(d1 - d2) < 1e-6);
  // ... and it equals the throws read where Top_Sand meets each fault: FA lowers it, FB raises it
  // (the frame holds the block north of FA).
  const auto& fa = net->faults()[static_cast<std::size_t>(d.fault_index("FA"))];
  const auto& fb = net->faults()[static_cast<std::size_t>(d.fault_index("FB"))];
  double ta = fa.throw_at(0, fa.layer_b(sand, 0)), tb = fb.throw_at(0, fb.layer_b(sand, 0));
  CHECK(ta > 1 && tb > 1);
  CHECK_NEAR(d1, -ta + tb, 1e-6);
  // Restoring a present point on the horizon lands on the unfaulted horizon.
  Vec3 p{4000, 2000, z1};
  Vec3 q;
  m.restore(std::span<const Vec3>(&p, 1), std::span<Vec3>(&q, 1));
  CHECK_NEAR(q.z, st->z(static_cast<std::size_t>(sand), q.x, q.y), 1e-6);
}

TEST_CASE("L2 the layered example builds with every check passing; mixing hanging-wall types is refused") {
  register_all_plugins();
  fs::path out = fs::path(STX_TMP_DIR) / "sheet_example";
  fs::remove_all(out);
  BuildOptions opt;
  auto o = run_build(kExample, out, opt, "test");
  CHECK(o.checks_passed);
  // one fault layer_rigid, the other rigid: not allowed
  std::string t = slurp(kExample);
  std::size_t first = t.find("hanging_wall: {type: layer_rigid}");
  REQUIRE(first != std::string::npos);
  t.replace(first, std::string("hanging_wall: {type: layer_rigid}").size(), "hanging_wall: {type: rigid}");
  fs::path mixed = fs::path(STX_TMP_DIR) / "mixed_layered.yaml";
  write_file(mixed, t);
  bool threw = false;
  try {
    auto d = desc::load_description(mixed);
    build_only(d, opt);
  } catch (const Error& e) {
    threw = std::string(e.what()).find("layer_rigid") != std::string::npos;
  }
  CHECK(threw);
}

TEST_CASE("L3 a tip on the unfaulted stratigraphy ignores the older fault's offset") {
  register_all_plugins();
  BuildOptions opt;
  // FB (younger) tip within CHALK: in its own frame vs on the unfaulted stratigraphy
  std::string t = slurp(kExample);
  t = replace_all(t, "upper_tip: {rule: horizon, horizon: Top_Chalk, offset: \"0 m | user_instruction | worked example\"}",
                  "upper_tip: {rule: horizon, horizon: Top_Sand, offset: \"0 m | user_instruction | worked example\"}");
  fs::path a = fs::path(STX_TMP_DIR) / "tips_frame.yaml";
  write_file(a, t);
  std::string tu = replace_all(t, "upper_tip: {rule: horizon, horizon: Top_Sand,", "upper_tip: {rule: horizon, horizon: Top_Sand, stratigraphy: unfaulted,");
  fs::path b = fs::path(STX_TMP_DIR) / "tips_unfaulted.yaml";
  write_file(b, tu);
  auto ma = build_only(desc::load_description(a), opt).model;
  auto mb = build_only(desc::load_description(b), opt).model;
  auto tip = [&](const Model& m) {
    auto* net = static_cast<const FaultNetwork*>(m.capability("faults"));
    const auto& f = net->faults()[1];
    double bu, bl;
    f.tips(0, bu, bl);
    return -f.at(0, bu).z;
  };
  auto* st = static_cast<const Structure*>(mb->capability("structure"));
  const auto& fb = static_cast<const FaultNetwork*>(mb->capability("faults"))->faults()[1];
  double bu, bl;
  fb.tips(0, bu, bl);
  Vec3 p = fb.at(0, bu);
  // unfaulted: the tip sits exactly on the unfaulted Top_Sand at its own map position
  CHECK_NEAR(p.z, st->z(3, p.x, p.y), 1e-5);
  // in FB's frame FA has displaced Top_Sand there, so the two tips differ
  CHECK(std::fabs(tip(*ma) - tip(*mb)) > 1.0);
}

TEST_CASE("S1 section sheet: key numbers from the model and a DejaVu Sans figure") {
  register_all_plugins();
  fs::path out = fs::path(STX_TMP_DIR) / "sheet_example";
  if (!fs::exists(out / "figures")) {
    BuildOptions opt;
    run_build(kExample, out, opt, "test");
  }
  fs::path csv = out / "figures" / "example_layered_faults_v1.0.0_SN_sheet_summary.csv";
  REQUIRE(fs::exists(csv));
  std::string s = slurp(csv);
  CHECK(s.find("\"top Sand at W1 (m TVDSS)\",\"1200.0\"") != std::string::npos);  // the column is given at the well
  CHECK(s.find("n/a") == std::string::npos);
  for (const char* ext : {".png", ".svg", ".tif", ".pptx"})
    CHECK(fs::exists(out / "figures" / (std::string("STX_example_layered_faults_v1.0.0_sheet_SN_sheet") + ext)));
  std::string svg = slurp(out / "figures" / "STX_example_layered_faults_v1.0.0_sheet_SN_sheet.svg");
  CHECK(svg.find("'DejaVu Sans'") != std::string::npos);
  CHECK(svg.find("6 units / 7 horizons") != std::string::npos);  // title placeholders filled
  CHECK(svg.find("two planar faults") != std::string::npos);
  CHECK(svg.find("1210\xE2\x80\x93" "1270 m TVDSS") != std::string::npos);  // completion annotation
  // the computed throw equals the fault's layer throw at Top_Sand at the section (fault centre)
  auto d = desc::load_description(kExample);
  BuildOptions opt;
  auto o = build_only(d, opt);
  auto* net = static_cast<const FaultNetwork*>(o.model->capability("faults"));
  const auto& fa = net->faults()[static_cast<std::size_t>(d.fault_index("FA"))];
  double ta = fa.throw_at(0, fa.layer_b(d.horizon_index("Top_Sand"), 0));
  CHECK(s.find("\"FA: throw at top Sand (m)\",\"" + fmt_fixed(ta, 1) + "\"") != std::string::npos);
}

TEST_CASE("S2 DejaVu Sans metrics and UTF-8 text widths") {
  using fig::text_width;
  CHECK_NEAR(text_width("A", 10, false, "DejaVu Sans"), 6.84, 1e-9);
  CHECK_NEAR(text_width("\xC3\xA5", 10, false, "DejaVu Sans"), text_width("a", 10, false, "DejaVu Sans"), 1e-9);  // å
  CHECK(text_width("AB", 10, true, "DejaVu Sans") > text_width("AB", 10, false, "DejaVu Sans"));
  CHECK_NEAR(text_width("Gassum", 10, false), text_width("Gassum", 10, false, ""), 1e-12);  // default metrics unchanged
}

TEST_CASE("S3 description rules: completions, outer units, sheets") {
  std::string t = slurp(kExample);
  auto issues = [&](const std::string& text) {
    desc::Description d;
    auto r = desc::validate_description_text(text, kExample.parent_path(), &d);
    return r.text();
  };
  CHECK(issues(t).empty());
  // completion base above its top
  CHECK(issues(replace_all(t, "base: \"1270 m", "base: \"1150 m")).find("W4") != std::string::npos);
  // the model_base unit must sit below the last horizon
  CHECK(issues(replace_all(t, "top: Top_Salt, base: model_base", "top: Top_Base, base: model_base")).find("S2") != std::string::npos);
  // a tip within the outer unit is refused
  CHECK(issues(replace_all(t, "unit: CHALK, stratigraphy", "unit: SALT, stratigraphy")).find("F1") != std::string::npos);
  // unknown references in a sheet
  CHECK(issues(replace_all(t, "fault: FA, side: footwall", "fault: FX, side: footwall")).find("O1") != std::string::npos);
  CHECK(issues(replace_all(t, "section: SN\n", "section: NS\n")).find("O1") != std::string::npos);
  // the named outer unit carries its colour into the model
  desc::Description d;
  REQUIRE(desc::validate_description_text(t, kExample.parent_path(), &d).ok());
  CHECK(d.below_base_name == "SALT");
  CHECK(d.wells[0].completions.size() == 1);
  CHECK(d.units[0].label == "Shale Gp");
}

STX_TEST_MAIN()
