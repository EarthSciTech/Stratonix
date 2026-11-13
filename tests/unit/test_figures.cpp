// Figure layer: encoders (PNG, TIFF, PPTX zip), SVG text, figure set of the synthetic example.
#include <cmath>
#include <cstdio>
#include <filesystem>

#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/figures/canvas.hpp"
#include "stratonix/package/pipeline.hpp"
#include "../stx_test.hpp"

using namespace stratonix;
using namespace stratonix::fig;
namespace fs = std::filesystem;

static std::uint32_t be32(const std::string& s, std::size_t o) {
  return (static_cast<std::uint32_t>(static_cast<unsigned char>(s[o])) << 24) | (static_cast<unsigned char>(s[o + 1]) << 16) |
         (static_cast<unsigned char>(s[o + 2]) << 8) | static_cast<unsigned char>(s[o + 3]);
}

static Figure sample() {
  Figure f(200, 100);
  f.rect(10, 10, 80, 40, rgb(200, 30, 30), rgb(0, 0, 0), 1);
  f.text({20, 80}, "Depth & <units>", 10, Anchor::Start, true);
  f.line({0, 0}, {200, 100}, rgb(0, 0, 200), 2);
  return f;
}

TEST_CASE("F1 SVG keeps text as text and escapes it") {
  std::string svg = to_svg(sample(), "t");
  CHECK(svg.find("<text") != std::string::npos);
  CHECK(svg.find("Depth &amp; &lt;units&gt;") != std::string::npos);
  CHECK(svg.find("<path") == std::string::npos || svg.find("Depth") != std::string::npos);
}

TEST_CASE("F2 PNG at 300 dpi, TIFF LZW round trip equals the raster") {
  if (!raster_available()) { std::puts("  (no rasteriser: skipped)"); return; }
  Raster r = rasterise(sample(), 300);
  CHECK(r.w == static_cast<int>(std::lround(200.0 / 72 * 300)));
  std::string png = encode_png(r, 300);
  CHECK(png.substr(1, 3) == "PNG");
  auto p = png.find("pHYs");
  REQUIRE(p != std::string::npos);
  CHECK(be32(png, p + 4) == 11811);  // 300 dpi = 11811 px/m
  CHECK(png[p + 12] == 1);           // unit: metre
  std::string tif = encode_tiff_lzw(r, 300, "test");
  Raster back = decode_tiff_lzw_for_test(tif);
  CHECK(back.w == r.w);
  CHECK(back.h == r.h);
  CHECK(back.rgb == r.rgb);
}

TEST_CASE("F3 PPTX is a zip holding one slide with picture and notes") {
  if (!raster_available()) return;
  Raster r = rasterise(sample(), 150);
  std::string png = encode_png(r, 150);
  std::string z = make_pptx({Slide{"Title", "Caption", "Notes", png, r.w, r.h}}, "deck");
  CHECK(z.substr(0, 4) == std::string("PK\x03\x04", 4));
  for (const char* part : {"[Content_Types].xml", "ppt/presentation.xml", "ppt/slides/slide1.xml", "ppt/notesSlides/notesSlide1.xml",
                           "ppt/media/image1.png"})
    CHECK(z.find(part) != std::string::npos);
}

TEST_CASE("F4 figure set of the synthetic example in four formats with footer") {
  register_all_plugins();
  fs::path out = fs::path(STX_TMP_DIR) / "figures_synthetic";
  fs::remove_all(out);
  BuildOptions o;
  auto r = run_build(fs::path(STX_SOURCE_DIR) / "examples/synthetic_fold_two_faults.yaml", out, o);
  CHECK(r.checks_passed);
  fs::path fd = out / "figures";
  REQUIRE(fs::exists(fd));
  std::size_t svgs = 0;
  for (auto& e : fs::directory_iterator(fd)) {
    if (e.path().extension() != ".svg") continue;
    ++svgs;
    std::string s = read_file(e.path());
    CHECK(s.find("Coordinate system: local engineering frame") != std::string::npos);
    CHECK(s.find("Evidence basis: DESIGN") != std::string::npos);
    CHECK(e.path().filename().string().rfind("STX_synthetic_fold_two_faults_v1.0.0_", 0) == 0);
    if (raster_available()) {
      std::string stem = e.path().filename().string();
      stem.resize(stem.size() - 4);
      for (const char* ext : {".png", ".tif", ".pptx"}) CHECK(fs::exists(fd / (stem + ext)));
    }
  }
  CHECK(svgs == 10);  // 2 sections, 5 maps, 2 3D views, 1 fence
}

STX_TEST_MAIN()
