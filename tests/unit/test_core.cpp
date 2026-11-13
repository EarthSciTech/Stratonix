#include <atomic>
#include <filesystem>

#include "../stx_test.hpp"
#include "stratonix/core/grid.hpp"
#include "stratonix/core/mesh.hpp"
#include "stratonix/core/units.hpp"
#include "stratonix/core/util.hpp"

using namespace stratonix;

TEST_CASE("sha256 known vectors") {
  CHECK(sha256_hex("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
  CHECK(sha256_hex("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
  CHECK(sha256_hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
  std::string million(1000000, 'a');
  CHECK(sha256_hex(million) == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}

TEST_CASE("number formatting is locale independent") {
  CHECK(fmt_fixed(1234.5678, 2) == "1234.57");
  CHECK(fmt_fixed(-0.0001, 3) == "0.000");
  CHECK(fmt_fixed(-2.5, 1) == "-2.5");
  CHECK(fmt_shortest(0.1) == "0.1");
  CHECK_NEAR(parse_double("1.5e3"), 1500, 0);
  CHECK_THROWS(parse_double("1,5"));
}

TEST_CASE("units") {
  CHECK(find_unit("km")->to_si == 1000);
  CHECK(find_unit("deg")->dim == Dimension::Angle);
  CHECK(!find_unit("furlong"));
}

TEST_CASE("grid bilinear sampling") {
  GridSpec g{0, 0, 10, 10, 3, 3};
  Grid2D gr(g, 0);
  for (std::size_t j = 0; j < 3; ++j)
    for (std::size_t i = 0; i < 3; ++i) gr.at(i, j) = 2 * g.x(i) + 3 * g.y(j);
  CHECK_NEAR(gr.sample(5, 5), 25, 1e-12);
  CHECK_NEAR(gr.sample(20, 20), 100, 1e-12);
  CHECK(std::isnan(gr.sample(25, 5)));
}

TEST_CASE("mesh quality") {
  TriMesh m;
  m.v = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
  m.t = {{0, 1, 2}, {1, 3, 2}};
  auto q = mesh_quality(m);
  CHECK(q.triangles == 2);
  CHECK(q.degenerate == 0);
  CHECK(q.boundary_edges == 4);
  CHECK(!q.watertight());
  CHECK_NEAR(m.area(), 1.0, 1e-12);
}

TEST_CASE("parallel_for covers every index once for any thread count") {
  for (unsigned t : {1u, 2u, 3u, 8u}) {
    set_thread_count(t);
    std::vector<int> hits(1003, 0);
    parallel_for(hits.size(), 17, [&](std::size_t a, std::size_t b) {
      for (std::size_t i = a; i < b; ++i) hits[i]++;
    });
    bool ok = true;
    for (int h : hits) ok = ok && h == 1;
    CHECK(ok);
  }
  set_thread_count(0);
  CHECK_THROWS(parallel_for(10, 1, [](std::size_t a, std::size_t) { if (a == 5) throw std::runtime_error("x"); }));
}

TEST_CASE("files are written byte-exact") {
  auto p = std::filesystem::path(stxtest::tmp_dir()) / "core_file.txt";
  write_file(p, "a\nb\n");
  CHECK(read_file(p) == "a\nb\n");
  CHECK(sha256_file(p) == sha256_hex("a\nb\n"));
}

STX_TEST_MAIN()
