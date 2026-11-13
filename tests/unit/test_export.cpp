// Petrel export set: round trips against the model, determinism, locale independence, golden hashes.
#include <clocale>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <locale>
#include <map>
#include <set>
#include <sstream>

#include "../stx_test.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/io/export.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

using namespace stratonix;
namespace fs = std::filesystem;

namespace {
struct Fixture {
  desc::Description d;
  ModelPtr m;
  std::unique_ptr<io::ExportContext> c;
  fs::path out;
};

Fixture& fx() {
  static Fixture f = [] {
    register_all_plugins();
    Fixture x;
    x.d = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples/synthetic_fold_two_faults.yaml");
    BuildLog log;
    x.m = build_model(x.d, {}, log);
    x.c = std::make_unique<io::ExportContext>(*x.m, make_well_paths(x.d), GridSpec::covering(x.m->info().domain, 50), 25, 250);
    x.out = fs::path(STX_TMP_DIR) / "export_petrel_a";
    fs::remove_all(x.out);
    io::export_set(*x.c, io::petrel_set_ids(), x.out, "test");
    return x;
  }();
  return f;
}

std::map<std::string, std::string> folder_files(const fs::path& root) {
  std::map<std::string, std::string> r;
  for (auto& e : fs::recursive_directory_iterator(root))
    if (e.is_regular_file()) r[fs::relative(e.path(), root).generic_string()] = read_file(e.path());
  return r;
}
std::string normalised(std::string s) {
  const std::string v = build_info().version;
  for (std::size_t p; (p = s.find("Stratonix " + v)) != std::string::npos;) s.replace(p, 10 + v.size(), "Stratonix VERSION");
  return s;
}
}  // namespace

TEST_CASE("E1 IRAP and ZMAP faulted grids read back equal to the model, gaps preserved") {
  auto& f = fx();
  const auto& inf = f.m->info();
  for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
    const auto& ref = f.c->faulted(static_cast<int>(h)).z;
    auto a = io::read_irap_ascii(read_file(f.out / "horizons" / (inf.horizons[h].name + "_faulted.irap")));
    auto b = io::read_zmap(read_file(f.out / "horizons" / (inf.horizons[h].name + "_faulted.zmap")));
    REQUIRE(a.spec.nx == ref.spec.nx && a.spec.ny == ref.spec.ny);
    REQUIRE(b.spec.nx == ref.spec.nx && b.spec.ny == ref.spec.ny);
    CHECK_NEAR(a.spec.x0, ref.spec.x0, 1e-6);
    CHECK_NEAR(b.spec.dy, ref.spec.dy, 1e-6);
    double worst = 0;
    std::size_t mask_diff = 0, gaps = 0;
    for (std::size_t k = 0; k < ref.v.size(); ++k) {
      bool na = std::isnan(ref.v[k]);
      gaps += na;
      if (na != std::isnan(a.v[k]) || na != std::isnan(b.v[k])) ++mask_diff;
      if (!na) worst = std::max({worst, std::fabs(a.v[k] - ref.v[k]), std::fabs(b.v[k] - ref.v[k])});
    }
    CHECK(mask_diff == 0);
    CHECK(worst <= 0.5e-4 + 1e-9);
    if (h >= 1 && h <= 3) CHECK(gaps > 0);  // H2-H4 are cut by faults
  }
}

TEST_CASE("E2 per-block surfaces equal the faulted surface inside their block") {
  auto& f = fx();
  for (int h = 1; h <= 3; ++h) {
    const auto& fsu = f.c->faulted(h);
    for (BlockCode b : f.c->blocks_on(h)) {
      auto g = io::read_irap_ascii(read_file(f.out / "horizons/per_block" / (f.m->info().horizons[static_cast<std::size_t>(h)].name + "_block" +
                                                                               std::to_string(b) + ".irap")));
      double worst = 0;
      std::size_t inside = 0;
      for (std::size_t k = 0; k < g.v.size(); ++k)
        if (!std::isnan(fsu.z.v[k]) && fsu.block[k] == b) {
          ++inside;
          worst = std::max(worst, std::fabs(g.v[k] - fsu.z.v[k]));
        }
      CHECK(inside > 0);
      CHECK(worst <= 1e-3);
    }
  }
}

TEST_CASE("E3 fault sticks and polygons read back; sticks lie on the fault, polygons are closed") {
  auto& f = fx();
  const auto& inf = f.m->info();
  for (std::size_t ff = 0; ff < inf.faults.size(); ++ff) {
    auto ref = repr::fault_sticks(f.c->fault_surface(static_cast<int>(ff)), 250);
    auto got = io::read_irap_lines(read_file(f.out / "faults" / (inf.faults[ff].name + "_sticks.dat")));
    REQUIRE(got.size() == ref.size());
    std::vector<Vec3> pts;
    for (std::size_t s = 0; s < got.size(); ++s) {
      REQUIRE(got[s].size() == ref[s].size());
      for (std::size_t i = 0; i < got[s].size(); ++i) {
        CHECK(norm(got[s][i] - ref[s][i]) < 1e-4);
        pts.push_back(got[s][i]);
      }
    }
    FaultEval fe;
    f.m->evaluate_fault(static_cast<int>(ff), pts, fe);
    double worst = 0;
    for (double dd : fe.distance) worst = std::max(worst, std::fabs(dd));
    CHECK(worst < 0.01);
  }
  for (int h = 1; h <= 3; ++h) {
    auto polys = io::read_irap_lines(read_file(f.out / "horizons/fault_polygons" / (inf.horizons[static_cast<std::size_t>(h)].name + "_fault_polygons.dat")));
    std::size_t n = 0;
    for (std::size_t ff = 0; ff < inf.faults.size(); ++ff) n += f.c->fault_polygons(static_cast<int>(ff), h).size();
    CHECK(polys.size() == n);
    for (auto& p : polys) CHECK(norm(p.front() - p.back()) < 1e-9);
  }
}

TEST_CASE("E4 TSurf read back equals the meshes; fault vertices lie on the fault") {
  auto& f = fx();
  const auto& inf = f.m->info();
  for (std::size_t ff = 0; ff < inf.faults.size(); ++ff) {
    io::TSurfInfo ti;
    auto m = io::read_tsurf(read_file(f.out / "faults" / (inf.faults[ff].name + ".ts")), &ti);
    CHECK(ti.name == inf.faults[ff].name);
    CHECK(ti.z_positive == "Elevation");
    const auto& ref = f.c->fault_surface(static_cast<int>(ff)).mesh;
    CHECK(m.t.size() == ref.t.size());
    FaultEval fe;
    f.m->evaluate_fault(static_cast<int>(ff), m.v, fe);
    double worst = 0;
    for (double dd : fe.distance) worst = std::max(worst, std::fabs(dd));
    CHECK(worst < 0.01);
    for (std::size_t t = 0; t < m.t.size(); t += 97)
      for (int k = 0; k < 3; ++k) CHECK(norm(m.v[m.t[t][k]] - ref.v[ref.t[t][k]]) < 1e-4);
  }
  auto hm = io::read_tsurf(read_file(f.out / "horizons/H3_faulted.ts"));
  auto ref = repr::grid_to_mesh(f.c->faulted(2).z);
  CHECK(hm.t.size() == ref.t.size());
}

TEST_CASE("E5 Petrel points, well heads and well tops read back; tops lie on their horizons") {
  auto& f = fx();
  auto pts = io::read_petrel_points(read_file(f.out / "faults/fault_sticks_points.txt"));
  std::size_t n = 0;
  for (std::size_t ff = 0; ff < f.m->info().faults.size(); ++ff)
    for (auto& s : repr::fault_sticks(f.c->fault_surface(static_cast<int>(ff)), 250)) n += s.size();
  CHECK(pts.pts.size() == n);
  REQUIRE(pts.attr.size() == 2);
  CHECK(pts.attr[0].first == "Fault" && pts.attr[0].second == "String");
  CHECK(pts.attr[1].first == "Stick" && pts.attr[1].second == "Int");
  CHECK(pts.values.front()[0] == "F1");

  auto heads = io::read_well_heads(read_file(f.out / "wells/well_heads.txt"));
  REQUIRE(heads.size() == f.c->wells().size());
  for (std::size_t w = 0; w < heads.size(); ++w) {
    CHECK(heads[w].name == f.c->wells()[w].name());
    CHECK_NEAR(heads[w].x, f.c->wells()[w].x(), 0.005);
    CHECK_NEAR(heads[w].td_md, f.c->wells()[w].total_depth_md(), 0.005);
  }
  auto tops = io::read_well_tops(read_file(f.out / "wells/well_tops.txt"));
  CHECK(tops.size() == 15);
  std::vector<Vec3> p;
  for (auto& t : tops) p.push_back({t.x, t.y, t.z});
  EvalResult r;
  f.m->evaluate(p, kEvalScalars, r);
  const auto& inf = f.m->info();
  for (std::size_t i = 0; i < tops.size(); ++i) {
    int h = -1;
    for (std::size_t k = 0; k < inf.horizons.size(); ++k)
      if (inf.horizons[k].name == tops[i].surface) h = static_cast<int>(k);
    REQUIRE(h >= 0);
    CHECK(std::fabs(r.s(static_cast<std::size_t>(h), i)) < 0.01);  // within rounding of 2 decimals
    CHECK_NEAR(tops[i].md, -tops[i].z, 0.01);                      // vertical wells, KB at 0
  }
}

TEST_CASE("E6 files are deterministic, locale independent and use LF line endings only") {
  auto& f = fx();
  fs::path out2 = fs::path(STX_TMP_DIR) / "export_petrel_b";
  fs::remove_all(out2);
  // A locale with a decimal comma must not change any byte: a named one when installed (also
  // affects C printf), otherwise a C++ locale with a decimal-comma facet (affects streams).
  struct Comma : std::numpunct<char> {
    char do_decimal_point() const override { return ','; }
    char do_thousands_sep() const override { return '.'; }
    std::string do_grouping() const override { return "\3"; }
  };
  bool named = false;
  for (const char* name : {"de_DE.UTF-8", "de_DE.utf8", "fr_FR.UTF-8", "German_Germany.1252"}) {
    try {
      std::locale::global(std::locale(name));
      named = std::setlocale(LC_ALL, name) != nullptr;
      break;
    } catch (...) {
    }
  }
  if (!named) std::locale::global(std::locale(std::locale::classic(), new Comma));
  std::ostringstream probe;
  probe << 1234.5;
  CHECK(probe.str() == "1.234,5");  // the comma locale is really active for streams
  io::ExportContext c2(*f.m, make_well_paths(f.d), GridSpec::covering(f.m->info().domain, 50), 25, 250);
  io::export_set(c2, io::petrel_set_ids(), out2, "test");
  std::locale::global(std::locale::classic());
  std::setlocale(LC_ALL, "C");
  auto a = folder_files(f.out), b = folder_files(out2);
  CHECK(a.size() == b.size());
  for (auto& [k, v] : a) {
    CHECK(b.count(k) == 1);
    CHECK(b[k] == v);
    CHECK(v.find('\r') == std::string::npos);
    CHECK(v.find('\t') == std::string::npos);
    for (unsigned char ch : v) CHECK(ch < 128);
  }
}

TEST_CASE("E7 golden hashes: the same bytes on every platform (version string normalised)") {
  auto& f = fx();
  fs::path golden = fs::path(STX_SOURCE_DIR) / "tests/data/petrel_synthetic_sha256.txt";
  auto files = folder_files(f.out);
  std::string now;
  for (auto& [k, v] : files) now += sha256_hex(normalised(v)) + "  " + k + "\n";
  write_file(fs::path(STX_TMP_DIR) / "petrel_synthetic_sha256.txt", now);  // to refresh the golden file deliberately
  REQUIRE(fs::exists(golden));
  std::string want = read_file(golden);
  if (want != now) {
    std::istringstream a(want), b(now);
    std::string la, lb;
    while (std::getline(a, la) && std::getline(b, lb))
      if (la != lb) std::printf("  golden: %s\n  now:    %s\n", la.c_str(), lb.c_str());
  }
  CHECK(want == now);
}

TEST_CASE("E8 names are made safe for Petrel") {
  CHECK(io::petrel_name("Top Reservoir (2)") == "Top_Reservoir__2_");
  CHECK(io::petrel_name("F-01_a") == "F-01_a");
  CHECK(io::petrel_name("") == "_");
}

STX_TEST_MAIN()
