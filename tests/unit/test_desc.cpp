#include <filesystem>

#include "../stx_test.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/desc/description.hpp"

using namespace stratonix;
using namespace stratonix::desc;
namespace fs = std::filesystem;

static fs::path ex(const std::string& f) { return fs::path(stxtest::source_dir()) / "examples" / f; }

static bool has_rule(const ValidationReport& r, const std::string& rule, const std::string& path_part = {}) {
  for (auto& i : r.issues)
    if (i.error && i.rule == rule && (path_part.empty() || i.path.find(path_part) != std::string::npos)) return true;
  return false;
}

TEST_CASE("YAML and JSON parse to the same tree") {
  auto a = parse_document("a: 1\nb: [1, 2.5, x]\nc: {d: true, e: null}\n");
  auto b = parse_document(R"({"a": 1, "b": [1, 2.5, "x"], "c": {"d": true, "e": null}})");
  CHECK(to_json(a) == to_json(b));
  CHECK(a.get("b")->seq[1].number == 2.5);
  CHECK(a.get("c")->get("e")->is_null());
  CHECK_THROWS(parse_document("a: 1\na: 2\n"));
  CHECK_THROWS(parse_document("a: [1, 2\n"));
}

TEST_CASE("synthetic example validates") {
  Description d;
  auto r = validate_description(ex("synthetic_fold_two_faults.yaml"), &d);
  if (!r.ok()) std::cerr << r.text();
  REQUIRE(r.ok());
  CHECK(d.horizons.size() == 5);
  CHECK(d.faults.size() == 2);
  CHECK(d.faults[0].upper.rule == "within_unit");
  CHECK_NEAR(d.faults[0].dip.si(), 60.0 * 3.14159265358979323846 / 180, 1e-12);
  CHECK(d.faults[0].throw_max.evidence == Evidence::UserInstruction);
  CHECK(d.all_params().size() > 40);
  CHECK(d.closure_faults_seal.flag == true);
  CHECK(d.sections.size() == 2);
  CHECK(d.sections[0].polyline.rows() == 2);
}

TEST_CASE("flat example validates") {
  auto r = validate_description(ex("flat_three_layers.yaml"));
  if (!r.ok()) std::cerr << r.text();
  CHECK(r.ok());
}

TEST_CASE("broken: parameter without evidence class") {
  auto r = validate_description(ex("broken/broken_missing_evidence.yaml"));
  CHECK(!r.ok());
  CHECK(has_rule(r, "E1", "faults[0].throw.max"));
  CHECK(r.text().find("no evidence class") != std::string::npos);
}

TEST_CASE("broken: overlapping units") {
  auto r = validate_description(ex("broken/broken_overlapping_units.yaml"));
  CHECK(!r.ok());
  CHECK(has_rule(r, "S2"));
  CHECK(r.text().find("overlapping units") != std::string::npos);
}

TEST_CASE("broken: fault tip below the model base") {
  auto r = validate_description(ex("broken/broken_tip_below_base.yaml"));
  CHECK(!r.ok());
  CHECK(has_rule(r, "F2", "faults[0].extent.lower_tip"));
}

TEST_CASE("compact parameter forms") {
  const std::string base = R"(stratonix_schema: 1
model: {name: t, version: "1.0.0", status: DESIGN}
crs: {name: local, vertical_datum: datum, length_unit: m}
domain:
  origin: "[0, 0] m | assumed | test"
  size: {value: [1, 1], unit: km, evidence: assumed, source: test}
  depth_range: "[0, 1000] m | assumed | test"
  grid_spacing: "%SPACING%"
stratigraphy:
  reference: {point: "[500, 500] m | assumed | test"}
  horizons: [{name: A, depth: "100 m | assumed | test"}, {name: B, depth: "200 m | assumed | test"}]
  units: [{name: U, top: A, base: B}]
)";
  auto with = [&](const std::string& sp) {
    std::string t = base;
    t.replace(t.find("%SPACING%"), 9, sp);
    Description d;
    auto r = validate_description_text(t, ".", &d);
    return std::make_pair(r, d);
  };
  auto [r1, d1] = with("50 m | assumed | test");
  CHECK(r1.ok());
  CHECK_NEAR(d1.size.si(0), 1000, 0);
  auto [r2, d2] = with("50 m | test");
  CHECK(!r2.ok());
  auto [r3, d3] = with("50 deg | assumed | test");
  CHECK(has_rule(r3, "U1"));
  auto [r4, d4] = with("50 m | guessed | test");
  CHECK(has_rule(r4, "E1"));
  auto [r5, d5] = with("50 | assumed | test");
  CHECK(has_rule(r5, "E1"));
}

TEST_CASE("schema version and migration hook") {
  auto doc = parse_document("stratonix_schema: 99\n");
  CHECK_THROWS(migrate(doc));
  auto ok = parse_document("stratonix_schema: 1\n");
  migrate(ok);
  auto r = validate_description_text("model: {}\n", ".");
  CHECK(!r.ok());
}

TEST_CASE("embedded schema is valid JSON with $defs") {
  auto s = parse_document(schema_v1_text());
  CHECK(s.get("$defs") != nullptr);
  CHECK(s.get("$defs")->get("param") != nullptr);
}

STX_TEST_MAIN()
