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

namespace {
/// Value column of the summary row whose quantity is q ("" if absent).
std::string summary_value(const std::string& csv, const std::string& q, std::string* basis = nullptr) {
  std::string key = "\"" + q + "\",\"";
  std::size_t p = csv.find(key);
  if (p == std::string::npos) return {};
  p += key.size();
  std::size_t e = csv.find("\",\"", p);
  if (basis) {
    std::size_t b = e + 3, be = csv.find('"', b);
    *basis = csv.substr(b, be - b);
  }
  return csv.substr(p, e - p);
}
}  // namespace

TEST_CASE("S4 summary table: every quantity with its value and basis; key numbers bottom right") {
  register_all_plugins();
  fs::path out = fs::path(STX_TMP_DIR) / "sheet_example";
  if (!fs::exists(out / "figures")) {
    BuildOptions opt;
    run_build(kExample, out, opt, "test");
  }
  std::string s = slurp(out / "figures" / "example_layered_faults_v1.0.0_SN_sheet_summary.csv");
  CHECK(s.rfind("quantity,value,basis\n", 0) == 0);
  CHECK(s.find('{') == std::string::npos);  // basis placeholders filled
  CHECK(s.find("n/a") == std::string::npos);
  CHECK(summary_value(s, "section") == "S-N through the crest");
  CHECK(summary_value(s, "W1 completion (m TVDSS)") == "1210-1270");
  CHECK(summary_value(s, "vertical exaggeration of the figure") == "2.0");
  CHECK(summary_value(s, "FA: dip") == "55 deg towards south");
  CHECK(summary_value(s, "FB: units holding the tips (upper / lower)").rfind("Chalk Fm / ", 0) == 0);  // tip on top Chalk
  std::string basis;
  summary_value(s, "FA: top Sand footwall / hanging wall at the fault (m TVDSS)", &basis);
  CHECK(basis.find("footwall block ") == 0);
  CHECK(basis != "footwall block ?, hanging wall block ?");
  // planes: FA at (Y 6500, 1100 m) dips S 55 deg; FB at (Y 5500, 1100 m) dips N 60 deg. At FA's lower tip (1800 m):
  const double kd = kPi / 180;
  double y_fa = 6500 - (1800 - 1100) / std::tan(55 * kd), y_fb = 5500 + (1800 - 1100) / std::tan(60 * kd);
  CHECK(summary_value(s, "FA: horizontal gap to FB at its lower tip (m)") == fmt_fixed(std::fabs(y_fa - y_fb), 0));

  // overlap and throw ratio: thickness = vertical distance between the two footwall cut-offs on the plane,
  // found here independently by marching the plane through the footwall block.
  auto d = desc::load_description(kExample);
  BuildOptions opt;
  auto o = build_only(d, opt);
  const Model& m = *o.model;
  auto plane_z = [&](double y) { return -(1100 + (6500 - y) * std::tan(55 * kd)); };  // FA elevation
  auto fw_block = [&](double y) {
    Vec3 p{4000, y + 5, plane_z(y) + 1};  // just north of FA: its footwall
    EvalResult r;
    m.evaluate(std::span<const Vec3>(&p, 1), kEvalBlock, r);
    return r.block[0];
  };
  auto cut = [&](int h) {
    double lo = 5000, hi = 7900;
    auto g = [&](double y) { return horizon_in_block(m, h, fw_block(y), 4000, y) - plane_z(y); };
    double glo = g(lo), ghi = g(hi);
    REQUIRE((glo > 0) != (ghi > 0));
    for (int i = 0; i < 60; ++i) {
      double mid = 0.5 * (lo + hi), gm = g(mid);
      if ((gm > 0) == (glo > 0)) { lo = mid; glo = gm; } else hi = mid;
    }
    return plane_z(0.5 * (lo + hi));
  };
  double z_sand = cut(d.horizon_index("Top_Sand")), z_base = cut(d.horizon_index("Top_Base")), z_seal = cut(d.horizon_index("Top_Seal"));
  double t = std::stod(summary_value(s, "FA: throw at top Sand (m)"));
  CHECK_NEAR(std::stod(summary_value(s, "FA: Sand-on-Sand overlap (fraction)")), 1 - t / (z_sand - z_base), 0.006);
  CHECK_NEAR(std::stod(summary_value(s, "FA: throw / seal thickness")), t / (z_seal - z_sand), 0.006);
  CHECK_NEAR(std::stod(summary_value(s, "FA: Y at top Sand, footwall (m)")),
             6500 - (-z_sand - 1100) / std::tan(55 * kd), 1.0);
  // fold relief and limb dip on the unfaulted structure along the section
  auto* st = static_cast<const Structure*>(m.capability("structure"));
  const std::size_t hs = static_cast<std::size_t>(d.horizon_index("Top_Sand"));
  double crest = -1e30, mx = 0, prev = st->z(hs, 4000, 0);
  std::size_t ic = 0;
  for (std::size_t i = 0; i <= 12000; ++i) {
    double z = st->z(hs, 4000, static_cast<double>(i));
    if (z > crest) { crest = z; ic = i; }
  }
  for (std::size_t i = 1; i <= ic; ++i) {
    double z = st->z(hs, 4000, static_cast<double>(i));
    mx = std::max(mx, std::fabs(z - prev));
    prev = z;
  }
  CHECK_NEAR(std::stod(summary_value(s, "top Sand, crest to northern end (m)")), crest - st->z(hs, 4000, 12000), 0.06);
  CHECK_NEAR(std::stod(summary_value(s, "south limb max dip at top Sand (deg)")), std::atan(mx) / kd, 0.1);

  // the sheet: "Key numbers" heading, rows on the sheet only, right-aligned in the bottom-right corner
  std::string svg = slurp(out / "figures" / "STX_example_layered_faults_v1.0.0_sheet_SN_sheet.svg");
  CHECK(svg.find(">Key numbers<") != std::string::npos);
  CHECK(svg.find("KEY NUMBERS") == std::string::npos);
  CHECK(svg.find("FA: dip") == std::string::npos);
  CHECK(svg.find("block: relief from FA") != std::string::npos);  // drawn (the apostrophe is escaped in SVG)
}

TEST_CASE("L4 layer-rigid fault shorter than the model: the throw tapers to zero at its ends") {
  register_all_plugins();
  std::string t = slurp(kExample);
  // FA: 6 km long, centred at X = 4000 m, cosine-bell strike profile
  std::size_t p = t.find("length: \"30000 m");
  t.replace(p, t.find('"', p + 9) - p + 1, "length: \"6000 m | user_instruction | test\"");
  p = t.find("taper: \"0.25 fraction");
  t.replace(p, t.find('"', p + 8) - p + 1, "taper: \"0.5 fraction | user_instruction | test\"");
  fs::path f = fs::path(STX_TMP_DIR) / "short_fault.yaml";
  { std::ofstream(f) << t; }
  auto d = desc::load_description(f);
  BuildOptions opt;
  auto o = build_only(d, opt);
  auto* net = static_cast<const FaultNetwork*>(o.model->capability("faults"));
  const auto& fa = net->faults()[static_cast<std::size_t>(d.fault_index("FA"))];
  const double b = fa.layer_b(d.horizon_index("Top_Sand"), 0);
  CHECK(fa.throw_at(0, b) > 20);
  CHECK(std::fabs(fa.throw_at(2990, fa.layer_b(d.horizon_index("Top_Sand"), 2990))) < 0.01);  // near the end
  CHECK(fa.throw_at(3010, b) == 0);  // beyond it
  bool noted = false;
  for (auto& n : o.log.notes) noted = noted || n.find("dies out") != std::string::npos;
  CHECK(noted);
  // the faulted top Sand is continuous across the fault end in FA's hanging wall
  const Model& m = *o.model;
  auto top_sand = [&](double x) {
    Vec3 q{x, 6000, -1150};
    EvalResult r;
    m.evaluate(std::span<const Vec3>(&q, 1), kEvalBlock, r);
    return horizon_in_block(m, d.horizon_index("Top_Sand"), r.block[0], x, 6000);
  };
  CHECK(std::fabs(top_sand(4000 + 2995) - top_sand(4000 + 3005)) < 0.05);
  CHECK(std::fabs(top_sand(4000) - top_sand(7500)) > 5);  // the offset is there at the centre
}

STX_TEST_MAIN()
