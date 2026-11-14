// Delivery package: contents, manifest, identifier scan, sharing status, determinism, golden hashes.
#include <filesystem>
#include <map>
#include <sstream>

#include <zlib.h>

#include "../stx_test.hpp"
#include "stratonix/checks/manifest.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/package/delivery.hpp"
#include "stratonix/package/pipeline.hpp"

using namespace stratonix;
namespace fs = std::filesystem;

namespace {
fs::path tmp(const std::string& n) {
  fs::path p = fs::path(STX_TMP_DIR) / n;
  fs::remove_all(p);
  fs::remove(p.string() + ".zip");
  fs::remove_all(p.string() + "_REJECTED_SCAN");
  return p;
}
fs::path example() { return fs::path(STX_SOURCE_DIR) / "examples/synthetic_fold_two_faults.yaml"; }
std::map<std::string, std::string> files_of(const fs::path& root) {
  std::map<std::string, std::string> r;
  for (auto& e : fs::recursive_directory_iterator(root))
    if (e.is_regular_file()) r[fs::relative(e.path(), root).generic_string()] = read_file(e.path());
  return r;
}
std::string without_timestamp(const std::string& m) {
  std::istringstream in(m);
  std::string l, o;
  while (std::getline(in, l))
    if (l.find("\"timestamp_utc\"") == std::string::npos) o += l + "\n";
  return o;
}
/// Raw pixel bytes of a PNG (IDAT chunks inflated): independent of the zlib version that compressed it.
std::string png_pixels(const std::string& png) {
  auto be32 = [&](std::size_t i) {
    return (std::uint32_t(std::uint8_t(png[i])) << 24) | (std::uint32_t(std::uint8_t(png[i + 1])) << 16) |
           (std::uint32_t(std::uint8_t(png[i + 2])) << 8) | std::uint32_t(std::uint8_t(png[i + 3]));
  };
  std::string idat, ihdr;
  for (std::size_t i = 8; i + 8 <= png.size();) {
    std::uint32_t len = be32(i);
    std::string type = png.substr(i + 4, 4);
    if (i + 12 + len > png.size()) break;
    if (type == "IHDR") ihdr = png.substr(i + 8, len);
    if (type == "IDAT") idat += png.substr(i + 8, len);
    i += 12 + len;
  }
  if (ihdr.size() < 13) return "png-undecodable";
  auto u32 = [&](std::size_t k) {
    return (std::uint32_t(std::uint8_t(ihdr[k])) << 24) | (std::uint32_t(std::uint8_t(ihdr[k + 1])) << 16) |
           (std::uint32_t(std::uint8_t(ihdr[k + 2])) << 8) | std::uint32_t(std::uint8_t(ihdr[k + 3]));
  };
  const std::size_t w = u32(0), h = u32(4), depth = std::uint8_t(ihdr[8]), ctype = std::uint8_t(ihdr[9]);
  const std::size_t ch = ctype == 6 ? 4 : ctype == 2 ? 3 : ctype == 4 ? 2 : 1;
  uLongf n = static_cast<uLongf>(h * (1 + (w * ch * depth + 7) / 8));
  std::string raw(n, '\0');
  if (uncompress(reinterpret_cast<Bytef*>(raw.data()), &n, reinterpret_cast<const Bytef*>(idat.data()),
                 static_cast<uLong>(idat.size())) != Z_OK)
    return "png-undecodable";
  raw.resize(n);
  return raw;
}
std::string base64_decode(const std::string& b) {
  static const std::string a = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  std::string o;
  int val = 0, bits = -8;
  for (char c : b) {
    auto k = a.find(c);
    if (k == std::string::npos) continue;
    val = (val << 6) + static_cast<int>(k);
    bits += 6;
    if (bits >= 0) {
      o.push_back(static_cast<char>((val >> bits) & 0xFF));
      bits -= 8;
    }
  }
  return o;
}
std::string normalised(std::string s) {
  const std::string v = build_info().version;
  for (std::size_t p; (p = s.find("Stratonix " + v)) != std::string::npos;) s.replace(p, 10 + v.size(), "Stratonix VERSION");
  // PNG images embedded in SVG: compare their pixels, not the zlib stream (differs between zlib versions)
  const std::string tag = "data:image/png;base64,";
  for (std::size_t p = s.find(tag); p != std::string::npos; p = s.find(tag, p + 1)) {
    std::size_t a = p + tag.size(), e = s.find_first_of("\"'", a);
    if (e == std::string::npos) break;
    std::string digest = "png-pixels:" + sha256_hex(png_pixels(base64_decode(s.substr(a, e - a))));
    s.replace(p, e - p, digest);
  }
  return s;
}
PackageConfig fast() {
  PackageConfig c;
  c.formats = {"petrel", "csv"};
  c.zip = true;
  return c;
}
}  // namespace

TEST_CASE("K1 the synthetic example packages; the manifest verifies; the scan passes") {
  register_all_plugins();
  auto out = tmp("package_a");
  PackageConfig cfg = fast();
  cfg.figures = true;
  auto r = run_package(example(), out, cfg, {});
  CHECK(r.checks_passed);
  CHECK(r.scan_passed);
  CHECK(r.hits.empty());
  CHECK(fs::exists(out.string() + ".zip"));
  CHECK(verify_manifest(out / "manifest.json").empty());
  std::string rd = read_file(out / "README.txt");
  for (const char* s : {"synthetic_fold_two_faults version 1.0.0", "DESIGN", "local engineering frame", "Vertical datum", "Z convention",
                        "Grid spacing:        50 m", "H1", "F2", "user_instruction", "Assumptions", "Known limitations", "Terms of use",
                        "Contact", "DRAFT"})
    CHECK(rd.find(s) != std::string::npos);
  auto d = desc::load_description(example());
  std::istringstream pc(read_file(out / "parameters.csv"));
  std::size_t lines = 0;
  for (std::string l; std::getline(pc, l);) ++lines;
  CHECK(lines == d.all_params().size() + 1);
  for (const char* f : {"checks_report.txt", "checks_report.csv", "summary_table.csv", "model_description.yaml", "exports/README.txt",
                        "exports/horizons/H2_faulted.irap", "exports/wells/well_tops.txt", "figures"})
    CHECK(fs::exists(out / f));
  // the internal note identifier in the exports README was translated at the boundary
  CHECK(read_file(out / "exports/README.txt").find("STX-TN") == std::string::npos);
}

TEST_CASE("K2 a contaminated package fails the scan and is not released") {
  register_all_plugins();
  fs::path dir = fs::path(STX_TMP_DIR) / "contaminated";
  fs::create_directories(dir);
  std::string y = read_file(example());
  const std::string orig = "synthetic example: faults treated as sealing";
  auto p = y.find(orig);
  REQUIRE(p != std::string::npos);
  y.replace(p, orig.size(), "see decision D-012 and project ACME-77");  // register id + a project code
  write_file(dir / "contaminated.yaml", y);
  auto out = tmp("package_contaminated");
  PackageConfig cfg = fast();
  cfg.figures = false;
  cfg.banned_patterns = {"ACME-[0-9]+"};
  auto r = run_package(dir / "contaminated.yaml", out, cfg, {});
  CHECK(!r.scan_passed);
  CHECK(!fs::exists(out));
  CHECK(fs::exists(out.string() + "_REJECTED_SCAN/SCAN_FAILED.txt"));
  CHECK(!fs::exists(out.string() + ".zip"));
  bool reg = false, proj = false;
  for (auto& h : r.hits) {
    reg |= h.excerpt.find("D-012") != std::string::npos;
    proj |= h.excerpt.find("ACME-77") != std::string::npos;
  }
  CHECK(reg);
  CHECK(proj);
  // translating the labels at the boundary makes the same model package cleanly
  auto out2 = tmp("package_translated");
  cfg.translations = {{"decision D-012", "the design basis"}, {"project ACME-77", "the client project"}};
  auto r2 = run_package(dir / "contaminated.yaml", out2, cfg, {});
  CHECK(r2.scan_passed);
  CHECK(read_file(out2 / "parameters.csv").find("the design basis") != std::string::npos);
}

TEST_CASE("K3 an input with unknown sharing status is refused before anything is written") {
  register_all_plugins();
  fs::path dir = fs::path(STX_TMP_DIR) / "unknown_sharing";
  fs::create_directories(dir);
  std::string y = read_file(example());
  auto p = y.find("wells:");
  y.insert(p, "inputs:\n  - {id: SRC1, path: notes.txt, role: reference, sharing: unknown, source: received by e-mail}\n");
  write_file(dir / "unknown.yaml", y);
  auto out = tmp("package_unknown");
  CHECK_THROWS(run_package(dir / "unknown.yaml", out, fast(), {}));
  CHECK(!fs::exists(out));
}

TEST_CASE("K4 the same package twice is identical apart from the manifest timestamp; golden hashes") {
  register_all_plugins();
  PackageConfig cfg = fast();
  cfg.zip = false;
  auto a = tmp("package_det_a"), b = tmp("package_det_b");
  set_thread_count(1);
  run_package(example(), a, cfg, {});
  set_thread_count(0);
  run_package(example(), b, cfg, {});
  auto fa = files_of(a), fb = files_of(b);
  CHECK(fa.size() == fb.size());
  for (auto& [k, v] : fa) {
    if (k == "manifest.json") CHECK(without_timestamp(v) == without_timestamp(fb[k]));
    else CHECK(v == fb[k]);
  }
  CHECK(fa["manifest.json"].find("compiler") == std::string::npos);  // portable manifest
  // Golden hashes of every text file except the manifest (raster figures depend on fonts).
  std::string now;
  for (auto& [k, v] : fa) {
    if (k == "manifest.json" || k.find(".png") != std::string::npos || k.find(".tif") != std::string::npos || k.find(".pptx") != std::string::npos)
      continue;
    now += sha256_hex(normalised(v)) + "  " + k + "\n";
  }
  write_file(fs::path(STX_TMP_DIR) / "package_synthetic_sha256.txt", now);
  fs::path golden = fs::path(STX_SOURCE_DIR) / "tests/data/package_synthetic_sha256.txt";
  REQUIRE(fs::exists(golden));
  std::string want = read_file(golden);
  if (want != now) {
    std::istringstream x(want), y(now);
    std::string lx, ly;
    while (std::getline(x, lx) && std::getline(y, ly))
      if (lx != ly) std::printf("  golden: %s\n  now:    %s\n", lx.c_str(), ly.c_str());
  }
  CHECK(want == now);
}

TEST_CASE("K5 package configuration file") {
  fs::path f = fs::path(STX_TMP_DIR) / "pkg.yaml";
  write_file(f, "formats: [petrel, vtk]\nfigures: false\nzip: true\nterms_of_use: \"Approved text.\"\ncontact: \"Model owner\"\n"
                "banned_patterns: ['PRJ-[0-9]{4}']\ntranslations: {PRJ-0001: the project}\n");
  auto c = PackageConfig::load(f);
  CHECK(c.formats.size() == 2);
  CHECK(!c.figures);
  CHECK(c.zip);
  CHECK(c.terms_of_use == "Approved text.");
  CHECK(c.banned_patterns.size() == 1);
  REQUIRE(c.translations.size() == 1);
  CHECK(c.translations[0].first == "PRJ-0001");
  write_file(f, "banned_patterns: ['(unclosed']\n");
  CHECK_THROWS(PackageConfig::load(f));
}

STX_TEST_MAIN()
