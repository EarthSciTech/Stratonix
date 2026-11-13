#include "../stx_test.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/model/model_suite.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/registry/registry.hpp"

using namespace stratonix;

static ModelInfo flat_info() {
  ModelInfo mi;
  mi.name = "flat";
  mi.domain = {0, 0, 1000, 1000, 0, -1000, 50};
  mi.horizons = {{"A", 0}, {"B", 1}, {"C", 2}};
  mi.units = {{"above", -1, 0}, {"U1", 0, 1}, {"U2", 1, 2}, {"below", 2, -1}};
  return mi;
}

TEST_CASE("flat-layers reference model passes the interface suite") {
  FlatLayersModel m(flat_info(), {-100, -300, -600});
  auto rep = run_model_suite(m);
  std::cout << rep.text();
  CHECK(rep.passed());
  CHECK(m.unit_at({10, 10, -50}) == 0);
  CHECK(m.unit_at({10, 10, -200}) == 1);
  CHECK(m.unit_at({10, 10, -900}) == 3);
}

TEST_CASE("flat-layers rejects bad input") {
  CHECK_THROWS(FlatLayersModel(flat_info(), {-100, -50, -600}));
  CHECK_THROWS(FlatLayersModel(flat_info(), {-100}));
}

TEST_CASE("a broken model fails the suite") {
  // Wrong unit tiling must be detected (S1).
  ModelInfo mi = flat_info();
  mi.units.pop_back();
  FlatLayersModel m(mi, {-100, -300, -600});
  CHECK(!run_model_suite(m).passed());
}

TEST_CASE("registry lists the builders and builds the flat example") {
  register_all_plugins();
  auto all = registry::Registry::global().all();
  CHECK(!all.empty());
  CHECK(registry::Registry::global().find(registry::Kind::StructureBuilder, "flat-layers") != nullptr);
  auto d = desc::load_description(std::string(stxtest::source_dir()) + "/examples/flat_three_layers.yaml");
  BuildOptions o;
  auto r = build_only(d, o);
  CHECK(r.model->info().engine == "flat-layers");
  CHECK(run_model_suite(*r.model).passed());
  CHECK(r.model->unit_at({1000, 1000, -600}) == 1);
}

STX_TEST_MAIN()
