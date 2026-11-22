// Stage 4C tests (STX-TN-16 section 9): cases D1 to D18 of real data. Written first (Prompt 54):
// they fail with "not implemented" until Prompts 55 to 57 build the features. The DLIS and SEG-Y
// files are written here by small encoders made from the published standards, so no third-party
// file is committed; CI also reads files written by dliswriter and segyio (test tools only).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "../stx_test.hpp"
#include "stratonix/core/crs.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/ages.hpp"
#include "stratonix/engines/surveys.hpp"
#include "stratonix/engines/timedepth.hpp"
#include "stratonix/io/maps.hpp"
#include "stratonix/io/seismic.hpp"
#include "stratonix/io/velocity.hpp"
#include "stratonix/io/welldata.hpp"
#include "stratonix/package/uncertainty.hpp"
#include "stratonix/repr/upscale.hpp"

using namespace stratonix;
namespace fs = std::filesystem;

namespace {
const double kPiD = 3.14159265358979323846;
double deg(double d) { return d * kPiD / 180.0; }
fs::path tmp(const std::string& name) {
  fs::path p = fs::path(STX_TMP_DIR) / "v4_data";
  fs::create_directories(p);
  return p / name;
}
void write_text(const fs::path& p, const std::string& s) { std::ofstream(p, std::ios::binary) << s; }
bool contains(const std::vector<std::string>& v, const std::string& what) {
  for (auto& s : v)
    if (s.find(what) != std::string::npos) return true;
  return false;
}

// ---- a minimal RP66 V1 (DLIS) encoder: storage unit label, one visible record, logical record
// segments for FILE-HEADER, ORIGIN, CHANNEL and FRAME, and FDATA records (API RP66 V1, POSC 1997)
struct Bytes {
  std::vector<std::uint8_t> b;
  void u8(unsigned v) { b.push_back(static_cast<std::uint8_t>(v)); }
  void u16(unsigned v) { u8(v >> 8); u8(v & 0xFF); }
  void u32(std::uint32_t v) { u16(v >> 16); u16(v & 0xFFFF); }
  void uvari(unsigned v) {
    if (v < 128) u8(v);
    else if (v < 16384) u16(0x8000 | v);
    else u32(0xC0000000u | v);
  }
  void ident(const std::string& s) { u8(static_cast<unsigned>(s.size())); for (char c : s) u8(static_cast<unsigned char>(c)); }
  void ascii(const std::string& s) { uvari(static_cast<unsigned>(s.size())); for (char c : s) u8(static_cast<unsigned char>(c)); }
  void obname(unsigned origin, unsigned copy, const std::string& id) { uvari(origin); u8(copy); ident(id); }
  void fdoubl(double v) { std::uint64_t x; std::memcpy(&x, &v, 8); for (int i = 7; i >= 0; --i) u8(static_cast<unsigned>((x >> (8 * i)) & 0xFF)); }
  void fsingl(float v) { std::uint32_t x; std::memcpy(&x, &v, 4); u32(x); }
  void add(const Bytes& o) { b.insert(b.end(), o.b.begin(), o.b.end()); }
};
enum : unsigned { kUSHORT = 15, kUVARI = 18, kIDENT = 19, kASCII = 20, kOBNAME = 23, kUNITS = 27, kFSINGL = 2, kFDOUBL = 7 };
// component descriptors: role in the top three bits, format bits below
const unsigned kSetT = 0xF0;            // SET with type
const unsigned kAttribL = 0x30;         // ATTRIB with label (template)
const unsigned kObjectN = 0x70;         // OBJECT with name
const unsigned kAttrCRV = 0x2D;         // ATTRIB with count, repcode and value
const unsigned kAttrCRUV = 0x2F;        // ATTRIB with count, repcode, units and value
Bytes lrs(bool eflr, unsigned type, const Bytes& body) {
  Bytes s;
  std::size_t len = 4 + body.b.size();
  const bool pad = len % 2 != 0;
  if (pad) ++len;
  s.u16(static_cast<unsigned>(len));
  s.u8((eflr ? 0x80 : 0x00) | (pad ? 0x01 : 0x00));
  s.u8(type);
  s.add(body);
  if (pad) s.u8(1);  // pad count
  return s;
}
std::vector<std::uint8_t> make_dlis(const std::vector<double>& depth, const std::vector<float>& gr) {
  Bytes fh;  // FILE-HEADER
  fh.u8(kSetT); fh.ident("FILE-HEADER");
  fh.u8(kAttribL); fh.ident("SEQUENCE-NUMBER");
  fh.u8(kAttribL); fh.ident("ID");
  fh.u8(kObjectN); fh.obname(1, 0, "0");
  fh.u8(kAttrCRV); fh.uvari(1); fh.u8(kASCII); fh.ascii("1");
  fh.u8(kAttrCRV); fh.uvari(1); fh.u8(kASCII); fh.ascii("STRATONIX TEST");
  Bytes og;  // ORIGIN
  og.u8(kSetT); og.ident("ORIGIN");
  og.u8(kAttribL); og.ident("WELL-NAME");
  og.u8(kObjectN); og.obname(1, 0, "DEFINING_ORIGIN");
  og.u8(kAttrCRV); og.uvari(1); og.u8(kASCII); og.ascii("TEST-1");
  Bytes ch;  // CHANNEL
  ch.u8(kSetT); ch.ident("CHANNEL");
  ch.u8(kAttribL); ch.ident("REPRESENTATION-CODE");
  ch.u8(kAttribL); ch.ident("UNITS");
  ch.u8(kAttribL); ch.ident("DIMENSION");
  auto chan = [&](const std::string& name, unsigned rc, const std::string& unit) {
    ch.u8(kObjectN); ch.obname(1, 0, name);
    ch.u8(kAttrCRV); ch.uvari(1); ch.u8(kUSHORT); ch.u8(rc);
    ch.u8(kAttrCRV); ch.uvari(1); ch.u8(kUNITS); ch.ident(unit);
    ch.u8(kAttrCRV); ch.uvari(1); ch.u8(kUVARI); ch.uvari(1);
  };
  chan("DEPT", kFDOUBL, "m");
  chan("GR", kFSINGL, "gAPI");
  Bytes fr;  // FRAME
  fr.u8(kSetT); fr.ident("FRAME");
  fr.u8(kAttribL); fr.ident("CHANNELS");
  fr.u8(kAttribL); fr.ident("INDEX-TYPE");
  fr.u8(kObjectN); fr.obname(1, 0, "MAIN");
  fr.u8(kAttrCRV); fr.uvari(2); fr.u8(kOBNAME); fr.obname(1, 0, "DEPT"); fr.obname(1, 0, "GR");
  fr.u8(kAttrCRV); fr.uvari(1); fr.u8(kIDENT); fr.ident("BOREHOLE-DEPTH");
  Bytes body;
  body.add(lrs(true, 0, fh));
  body.add(lrs(true, 1, og));
  body.add(lrs(true, 3, ch));
  body.add(lrs(true, 4, fr));
  for (std::size_t i = 0; i < depth.size(); ++i) {
    Bytes fd;
    fd.obname(1, 0, "MAIN");
    fd.uvari(static_cast<unsigned>(i + 1));
    fd.fdoubl(depth[i]);
    fd.fsingl(gr[i]);
    body.add(lrs(false, 0, fd));
  }
  Bytes out;
  const std::string sul = std::string("   1") + "V1.00" + "RECORD" + "08192" + std::string(60, ' ');
  for (char c : sul) out.u8(static_cast<unsigned char>(c));
  out.u16(static_cast<unsigned>(4 + body.b.size()));
  out.u8(0xFF);
  out.u8(0x01);
  out.add(body);
  return out.b;
}

// ---- a minimal SEG-Y encoder (SEG-Y revision 1 and 2.0): textual and binary headers, traces
std::uint32_t ibm_float(float f) {  // IEEE single to IBM System/360 single (SEG-Y format 1)
  if (f == 0.0f) return 0;
  std::uint32_t sign = f < 0 ? 0x80000000u : 0u;
  double v = std::fabs(static_cast<double>(f));
  int e = 0;
  while (v >= 1.0) { v /= 16.0; ++e; }
  while (v < 1.0 / 16.0) { v *= 16.0; --e; }
  const std::uint32_t frac = static_cast<std::uint32_t>(std::llround(v * 16777216.0)) & 0xFFFFFF;
  return sign | (static_cast<std::uint32_t>(e + 64) << 24) | frac;
}
struct SegyOut {
  std::vector<std::uint8_t> b;
  bool little = false;
  void put(std::size_t off, std::uint64_t v, int n) {
    for (int i = 0; i < n; ++i) b[off + (little ? i : n - 1 - i)] = static_cast<std::uint8_t>((v >> (8 * i)) & 0xFF);
  }
};
std::vector<std::uint8_t> make_segy(int revision, bool little, int nil, int nxl, int ns, double dt_us,
                                    const std::function<float(int, int, int)>& amp) {
  SegyOut s;
  s.little = little;
  s.b.assign(3600, 0);
  const std::string text = "C 1 STRATONIX TEST VOLUME";
  for (std::size_t i = 0; i < 3200; ++i) s.b[i] = static_cast<std::uint8_t>(i < text.size() ? text[i] : ' ');  // ASCII textual header
  const int fmt = revision >= 2 ? 5 : 1;
  s.put(3200 + 16, static_cast<std::uint64_t>(dt_us), 2);
  s.put(3200 + 20, static_cast<std::uint64_t>(ns), 2);
  s.put(3200 + 24, static_cast<std::uint64_t>(fmt), 2);
  s.put(3200 + 300, revision >= 2 ? 0x0200 : 0x0100, 2);   // revision 2.0 or 1.0
  if (little) s.put(3200 + 96, 0x01020304, 4);              // byte-order marker (revision 2)
  for (int i = 0; i < nil; ++i)
    for (int j = 0; j < nxl; ++j) {
      const std::size_t t0 = s.b.size();
      s.b.resize(t0 + 240 + 4 * static_cast<std::size_t>(ns), 0);
      s.put(t0 + 70, static_cast<std::uint64_t>(static_cast<std::uint16_t>(-100)), 2);  // coordinate scalar: divide by 100
      s.put(t0 + 72, static_cast<std::uint64_t>((1000 + 25 * j) * 100), 4);              // x
      s.put(t0 + 76, static_cast<std::uint64_t>((2000 + 25 * i) * 100), 4);              // y
      s.put(t0 + 114, static_cast<std::uint64_t>(ns), 2);
      s.put(t0 + 116, static_cast<std::uint64_t>(dt_us), 2);
      s.put(t0 + 188, static_cast<std::uint64_t>(100 + i), 4);  // inline (revision 1 position)
      s.put(t0 + 192, static_cast<std::uint64_t>(200 + j), 4);  // crossline
      for (int k = 0; k < ns; ++k) {
        const float a = amp(i, j, k);
        std::uint32_t w;
        if (fmt == 1) w = ibm_float(a);
        else std::memcpy(&w, &a, 4);
        s.put(t0 + 240 + 4 * static_cast<std::size_t>(k), w, 4);
      }
    }
  return s.b;
}
void write_bytes(const fs::path& p, const std::vector<std::uint8_t>& b) {
  std::ofstream(p, std::ios::binary).write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
}
}  // namespace

// ---- 1. Wells (Prompt 55)

TEST_CASE("D1 LAS 2.0: header, curves, units, nulls and values") {
  const std::string las =
      "~VERSION INFORMATION\n"
      " VERS.                  2.0 :   CWLS LOG ASCII STANDARD -VERSION 2.0\n"
      " WRAP.                   NO :   ONE LINE PER DEPTH STEP\n"
      "~WELL INFORMATION\n"
      "#MNEM.UNIT              DATA                       DESCRIPTION\n"
      " STRT.M              1670.0000 : START DEPTH\n"
      " STOP.M              1669.7500 : STOP DEPTH\n"
      " STEP.M               -0.1250 : STEP\n"
      " NULL.              -999.2500 : NULL VALUE\n"
      " WELL.   ANY ET AL OIL WELL #12 : WELL\n"
      "~CURVE INFORMATION\n"
      " DEPT.M                      :  1  DEPTH\n"
      " DT  .US/M                   :  2  SONIC TRANSIT TIME\n"
      " RHOB.K/M3                   :  3  BULK DENSITY\n"
      "~PARAMETER INFORMATION\n"
      " BHT .DEGC           35.5000 : BOTTOM HOLE TEMPERATURE\n"
      "~A  DEPTH     DT       RHOB\n"
      "1670.000   123.450 2550.000\n"
      "1669.875   123.450 -999.25\n"
      "1669.750   123.450 2550.000\n";
  auto w = io::parse_las(las, "inline D1");
  CHECK(w.version == "2.0");
  CHECK(w.well == "ANY ET AL OIL WELL #12");
  CHECK(w.index.mnemonic == "DEPT" && w.index.unit == "M");
  REQUIRE(w.curves.size() == 2);
  const auto* dt = w.find("DT");
  const auto* rhob = w.find("RHOB");
  REQUIRE(dt && rhob);
  CHECK(dt->unit == "US/M");
  CHECK(rhob->unit == "K/M3");
  CHECK_NEAR(rhob->null_value, -999.25, 0);
  CHECK(rhob->nulls == 1);
  REQUIRE(rhob->values.size() == 3);
  CHECK(std::isnan(rhob->values[1]));
  CHECK_NEAR(rhob->values[2], 2550.0, 0);
  CHECK_NEAR(w.index.values[1], 1669.875, 1e-9);
  CHECK(w.parameters.count("BHT") && w.parameters.at("BHT").unit == "DEGC");
  CHECK(contains(w.report, "null"));
}

TEST_CASE("D2 LAS 2.0 wrapped data, and header errors reported with their line") {
  const std::string wrapped =
      "~V\n VERS. 2.0 : CWLS\n WRAP. YES : MULTIPLE LINES PER DEPTH STEP\n"
      "~W\n STRT.M 910.0 :\n STOP.M 909.5 :\n STEP.M -0.5 :\n NULL. -999.25 :\n WELL. W2 : WELL\n"
      "~C\n DEPT.M : DEPTH\n A.V : A\n B.V : B\n C.V : C\n D.V : D\n"
      "~A\n910.0\n 1.0 2.0\n 3.0 4.0\n909.5\n 5.0 6.0\n 7.0 -999.25\n";
  auto w = io::parse_las(wrapped, "inline D2a");
  REQUIRE(w.curves.size() == 4);
  CHECK_NEAR(w.find("C")->values[1], 7.0, 0);
  CHECK(std::isnan(w.find("D")->values[1]));
  const std::string bad =
      "~V\n VERS. 2.0 : CWLS\n WRAP. NO : ONE\n"
      "~W\n STRT.M 10.0 :\n this line has no colon\n STOP.M 11.0 :\n STEP.M 1.0 :\n NULL. -999.25 :\n"
      "~C\n DEPT.M : DEPTH\n GR.API : GR\n"
      "~A\n10.0 50.0\n11.0 60.0\n";
  auto b = io::parse_las(bad, "inline D2b");
  CHECK(b.find("GR") && b.find("GR")->values.size() == 2);
  CHECK(contains(b.report, "line 6"));
}

TEST_CASE("D3 LAS 3.0: log, tops and inclinometry sections with a comma delimiter") {
  const std::string las =
      "~Version\n VERS. 3.0 : CWLS LOG ASCII STANDARD - VERSION 3.0\n WRAP. NO : ONE LINE PER DEPTH STEP\n DLM . COMMA : DELIMITING CHARACTER\n"
      "~Well\n STRT.M 1500.0 : \n STOP.M 1500.2 : \n STEP.M 0.1 : \n NULL. -999.25 : \n WELL. W3 : WELL\n"
      "~Log_Definition\n DEPT.M : Depth {F}\n GR.gAPI : Gamma ray {F}\n"
      "~Log_Data | Log_Definition\n1500.0,45.5\n1500.1,47.0\n1500.2,-999.25\n"
      "~Tops_Definition\n TOPT. : Top name {S}\n TOPD.M : Top depth {F}\n"
      "~Tops_Data | Tops_Definition\nH1,1200.5\nH2,1480.0\n"
      "~Inclinometry_Definition\n MD.M : Measured depth {F}\n INC.deg : Inclination {F}\n AZI.deg : Azimuth {F}\n"
      "~Inclinometry_Data | Inclinometry_Definition\n0,0,0\n1000,10,45\n1500,20,50\n";
  auto w = io::parse_las(las, "inline D3");
  CHECK(w.version == "3.0");
  REQUIRE(w.find("GR"));
  CHECK(w.find("GR")->nulls == 1);
  CHECK_NEAR(w.find("GR")->values[1], 47.0, 0);
  REQUIRE(w.tops.size() == 2);
  CHECK(w.tops[0].text.size() == 2 && w.tops[0].text[1] == "H2");
  CHECK_NEAR(w.tops[1].values[0], 1200.5, 0);
  REQUIRE(w.inclinometry.size() == 3);
  CHECK(w.inclinometry[1].unit == "deg");
  CHECK_NEAR(w.inclinometry[2].values[2], 50.0, 0);
}

TEST_CASE("D4 DLIS: two channels in one frame") {
  std::vector<double> depth;
  std::vector<float> gr;
  for (int i = 0; i < 20; ++i) { depth.push_back(1000 + 0.5 * i); gr.push_back(40.0f + static_cast<float>(i)); }
  const auto p = tmp("d4.dlis");
  write_bytes(p, make_dlis(depth, gr));
  auto sets = io::read_dlis(p);
  REQUIRE(sets.size() == 1);
  const auto& s = sets[0];
  CHECK(s.version == "RP66 V1");
  CHECK(s.well == "TEST-1");
  CHECK(s.index.mnemonic == "DEPT" && s.index.unit == "m");
  CHECK(s.index_type == "MD");
  REQUIRE(s.find("GR"));
  CHECK(s.find("GR")->unit == "gAPI");
  REQUIRE(s.find("GR")->values.size() == 20);
  CHECK_NEAR(s.index.values[19], 1009.5, 1e-12);
  CHECK_NEAR(s.find("GR")->values[7], 47.0, 1e-6);
}

TEST_CASE("D5 minimum curvature reproduces a circular arc; dogleg severity") {
  const double B = deg(3.0) / 30.0;  // 3 degrees per 30 m
  std::vector<SurveyStation> st;
  for (double md = 0; md <= 900 + 1e-9; md += 30) st.push_back({md, B * md, deg(30)});
  SurveySpec sp;
  sp.method = SurveyMethod::MinimumCurvature;
  sp.reference = AzimuthReference::GridNorth;
  auto pos = survey_positions(st, sp);
  REQUIRE(pos.size() == st.size());
  double worst = 0;
  for (std::size_t i = 0; i < st.size(); ++i) {
    const double md = st[i].md, h = (1 - std::cos(B * md)) / B, tvd = std::sin(B * md) / B;
    worst = std::max({worst, std::fabs(pos[i].north - h * std::cos(deg(30))), std::fabs(pos[i].east - h * std::sin(deg(30))),
                      std::fabs(pos[i].tvd - tvd)});
  }
  std::printf("  circular arc: largest position error %.2e m\n", worst);
  CHECK(worst < 1e-6);
  CHECK_NEAR(pos[5].dls, B, 1e-12);
}

TEST_CASE("D6 a straight inclined hole: every method agrees with the closed form") {
  std::vector<SurveyStation> st;
  for (double md = 0; md <= 2000; md += 100) st.push_back({md, md == 0 ? 0.0 : deg(30), deg(45)});
  // the first course builds from vertical to 30 degrees; test from the second station on
  for (auto m : {SurveyMethod::MinimumCurvature, SurveyMethod::RadiusOfCurvature, SurveyMethod::BalancedTangential,
                 SurveyMethod::AverageAngle, SurveyMethod::Tangential}) {
    SurveySpec sp;
    sp.method = m;
    sp.reference = AzimuthReference::TrueNorth;
    auto pos = survey_positions(st, sp);
    REQUIRE(pos.size() == st.size());
    const double dmd = 1000 - 100, dn = dmd * std::sin(deg(30)) * std::cos(deg(45)), dz = dmd * std::cos(deg(30));
    CHECK_NEAR(pos[10].north - pos[1].north, dn, 1e-6);
    CHECK_NEAR(pos[10].east - pos[1].east, dn, 1e-6);
    CHECK_NEAR(pos[10].tvd - pos[1].tvd, dz, 1e-6);
  }
}

TEST_CASE("D7 a survey without its azimuth reference is refused") {
  std::vector<SurveyStation> st{{0, 0, 0}, {100, deg(5), deg(90)}};
  SurveySpec sp;  // reference left unstated
  bool refused = false;
  try {
    survey_positions(st, sp);
  } catch (const Error& e) {
    refused = std::string(e.what()).find("azimuth reference") != std::string::npos;
  }
  CHECK(refused);
}

TEST_CASE("D8 upscaling a known log: means, most frequent, missing below the valid fraction") {
  repr::CornerPointGrid g;
  g.nx = g.ny = 1;
  g.nz = 3;
  g.x0 = 0; g.y0 = 0; g.dx = 100; g.dy = 100;
  for (int j = 0; j <= 1; ++j)
    for (int i = 0; i <= 1; ++i) {
      const double x = 100.0 * i, y = 100.0 * j;
      g.coord.insert(g.coord.end(), {x, y, -900.0, x, y, -1100.0});
    }
  g.zcorn.assign(8 * 3, 0);
  for (int k = 0; k < 3; ++k)
    for (int c = 0; c < 8; ++c) g.zcorn[g.zcorn_index(0, 0, k, c)] = (c & 4) ? -1010.0 - 10 * k : -1000.0 - 10 * k;
  g.actnum.assign(3, 1);
  g.zone.assign(3, 0);
  repr::PathLog log;
  log.well = "W";
  log.curve = "X";
  for (int s = 0; s < 30; ++s) {
    const double z = -1000.5 - s;
    log.p.push_back({50, 50, z});
    log.md.push_back(1000.5 + s);
    double v;
    if (s < 10) v = s + 1;                          // cell 0: 1 .. 10
    else if (s < 20) v = (s - 10) < 6 ? std::nan("") : 2.0;  // cell 1: 6 nulls of 10
    else v = (s - 20) < 3 ? 1 : (s - 20) < 7 ? 2 : 3;  // cell 2: facies 1, 2, 3
    log.v.push_back(v);
  }
  auto a = repr::upscale_log(g, log, repr::UpscaleMethod::Arithmetic, 0.5);
  auto h = repr::upscale_log(g, log, repr::UpscaleMethod::Harmonic, 0.5);
  auto ge = repr::upscale_log(g, log, repr::UpscaleMethod::Geometric, 0.5);
  auto mf = repr::upscale_log(g, log, repr::UpscaleMethod::MostFrequent, 0.5);
  auto val = [](const repr::UpscaledCells& u, std::size_t cell) -> double {
    for (std::size_t q = 0; q < u.cells.size(); ++q)
      if (u.cells[q] == cell) return u.value[q];
    return std::nan("");
  };
  double hs = 0;
  for (int i = 1; i <= 10; ++i) hs += 1.0 / i;
  double lg = 0;
  for (int i = 1; i <= 10; ++i) lg += std::log(i);
  CHECK_NEAR(val(a, 0), 5.5, 1e-12);
  CHECK_NEAR(val(h, 0), 10.0 / hs, 1e-12);
  CHECK_NEAR(val(ge, 0), std::exp(lg / 10), 1e-12);
  CHECK(std::isnan(val(a, 1)));
  CHECK(std::find(a.missing.begin(), a.missing.end(), std::size_t{1}) != a.missing.end());
  CHECK_NEAR(val(mf, 2), 2.0, 0);
  CHECK(a.method == "arithmetic mean");
  CHECK_NEAR(a.md_top[0], 1000.0, 0.51);
}

TEST_CASE("D9 ages: age at depth, correlation of three wells, an inversion reported") {
  std::vector<AgeDatum> a{{"chronostratigraphic", "s1", 1000, -1, 10, 0.2, "GTS2020", "t"},
                          {"biostratigraphic", "s2", 1500, -1, 20, 0.5, "GTS2020", "t"},
                          {"chronostratigraphic", "s3", 2000, -1, 40, 0.3, "GTS2020", "t"}};
  CHECK_NEAR(age_at_depth(a, 1250, AgeRule::LinearInDepth), 15.0, 1e-9);
  CHECK_NEAR(depth_of_age(a, 30, AgeRule::LinearInDepth), 1750.0, 1e-9);
  CHECK(std::isnan(age_at_depth(a, 2500, AgeRule::LinearInDepth)));
  std::vector<AgeDatum> b = a, c = a;
  for (auto& x : b) x.depth += 100;
  c[1].depth = 1100;
  c[1].age = 9;  // younger than the sample above it: an inversion
  auto cor = correlate_by_age({{"A", a}, {"B", b}, {"C", c}}, {{"H15", 15}, {"H50", 50}}, AgeRule::LinearInDepth);
  REQUIRE(cor.depth.size() == 2 && cor.depth[0].size() == 3);
  CHECK_NEAR(cor.depth[0][0], 1250, 1e-9);
  CHECK_NEAR(cor.depth[0][1], 1350, 1e-9);
  CHECK(std::isnan(cor.depth[1][0]));
  bool inversion = false;
  for (auto& k : cor.conflicts) inversion = inversion || (k.well == "C" && k.what.find("inversion") != std::string::npos);
  CHECK(inversion);
}

TEST_CASE("D7b descriptions: a survey file without its method, and ages without a correlation rule, are refused (F44)") {
  const fs::path dir = tmp("d7b");
  fs::create_directories(dir);
  write_text(dir / "s.csv", "md [m],inc [deg],azi [deg]\n0,0,0\n500,10,45\n");
  write_text(dir / "bad.csv", "md,inc,azi\n0,0,0\n500,10,45\n");
  auto yaml = [](const std::string& well) {
    return std::string("stratonix_schema: 1\nmodel: {name: t, version: \"1.0.0\", status: DESIGN}\n") +
           "crs: {name: local, epsg: 0, vertical_datum: datum, length_unit: m}\n"
           "domain: {origin: \"[0, 0] m | user_instruction | t\", size: \"[2000, 2000] m | user_instruction | t\", "
           "depth_range: \"[0, 2000] m | user_instruction | t\", grid_spacing: \"50 m | user_instruction | t\"}\n"
           "wells:\n  - {name: W1, role: control, head: \"[1000, 1000] m | user_instruction | t\", tops: {H1: \"800 m | user_instruction | t\"}}\n" +
           well + "stratigraphy: {reference: {well: W1}, horizons: [{name: H1}, {name: H2, thickness_above: \"200 m | user_instruction | t\"}], units: [{name: U, top: H1, base: H2}]}\n"
           "structure: {regional: {type: flat}}\n";
  };
  auto has = [](const desc::ValidationReport& r, const std::string& rule) {
    for (auto& i : r.issues)
      if (i.rule == rule && i.error) return true;
    return false;
  };
  const std::string head = "  - {name: W2, role: check, head: \"[500, 500] m | user_instruction | t\", ";
  auto r1 = desc::validate_description_text(yaml(head + "trajectory: {file: s.csv}}\n"), dir, nullptr);
  CHECK(has(r1, "F44"));
  auto r2 = desc::validate_description_text(yaml(head + "trajectory: {file: bad.csv, method: minimum_curvature, azimuth_reference: grid_north}}\n"), dir, nullptr);
  CHECK(has(r2, "F44"));  // the file states no units
  auto r3 = desc::validate_description_text(yaml(head + "trajectory: {file: s.csv, method: minimum_curvature, azimuth_reference: grid_north}}\n"), dir, nullptr);
  CHECK(r3.ok());
  auto r4 = desc::validate_description_text(
      yaml(head + "ages: [{kind: biostratigraphic, name: z, depth: \"900 m | assumed | t\", age: \"10 Ma | assumed | t\"}]}\n"), dir, nullptr);
  CHECK(has(r4, "F44"));
}

// ---- 2. Seismic, depth conversion, coordinates and maps (Prompt 56)

TEST_CASE("D10 SEG-Y revision 1 (IBM float) and 2.0 (IEEE, little-endian) read back") {
  auto amp = [](int i, int j, int k) { return static_cast<float>(std::sin(0.1 * k) * (1 + 0.1 * i) - 0.01 * j); };
  for (int rev : {1, 2}) {
    const bool little = rev == 2;
    const auto p = tmp("d10_rev" + std::to_string(rev) + ".sgy");
    write_bytes(p, make_segy(rev, little, 3, 4, 50, 4000, amp));
    auto v = io::read_segy(p);
    CHECK(v.revision == rev);
    CHECK(v.little_endian == little);
    CHECK(v.samples == 50);
    CHECK_NEAR(v.dt, 0.004, 1e-12);
    REQUIRE(v.traces.size() == 12);
    CHECK(v.traces[5].il == 101 && v.traces[5].xl == 201);
    CHECK_NEAR(v.traces[5].x, 1025.0, 1e-9);
    CHECK_NEAR(v.traces[5].y, 2025.0, 1e-9);
    double worst = 0;
    for (std::size_t t = 0; t < v.traces.size(); ++t)
      for (int k = 0; k < 50; ++k)
        worst = std::max(worst, std::fabs(static_cast<double>(v.traces[t].s[static_cast<std::size_t>(k)]) -
                                          amp(static_cast<int>(t / 4), static_cast<int>(t % 4), k)));
    CHECK(worst < 1e-6);
    CHECK_NEAR(v.sample(1025.0, 2025.0, 0.040), amp(1, 1, 10), 1e-6);
  }
}

TEST_CASE("D11 horizon picks and fault sticks keep their domain; no domain is refused") {
  const auto hp = tmp("d11_h1.dat");
  write_text(hp, "INLINE :      100 XLINE :      200     1000.00     2000.00     1500.00\n"
                 "INLINE :      100 XLINE :      201     1025.00     2000.00     1502.00\n"
                 "INLINE :      101 XLINE :      200     1000.00     2025.00     1504.00\n");
  auto ps = io::read_picks(hp, "charisma", "time", "ms");
  REQUIRE(ps.size() == 1);
  CHECK(ps[0].domain == "time" && ps[0].unit == "ms");
  REQUIRE(ps[0].picks.size() == 3);
  CHECK(ps[0].picks[1].xl == 201);
  CHECK_NEAR(ps[0].picks[2].p.z, 1504.0, 0);
  CHECK_THROWS(io::read_picks(hp, "charisma", "", ""));
  const auto sp = tmp("d11_sticks.csv");
  write_text(sp, "x,y,z,fault,stick\n1000,2000,1400,F1,1\n1010,2000,1600,F1,1\n1000,2100,1400,F1,2\n1012,2100,1620,F1,2\n");
  auto ss = io::read_fault_sticks(sp, "csv", "time", "ms");
  CHECK(ss.domain == "time");
  REQUIRE(ss.sticks.size() == 2);
  CHECK(ss.sticks[1].p.size() == 2 && ss.sticks[1].fault == "F1");
}

namespace {
// independent integration of dz/dt = v0 + k z over one-way time (RK4)
double rk4_depth(double v0, double k, double zt, double dt) {
  const int n = 20000;
  const double h = dt / n;
  double z = zt;
  auto f = [&](double zz) { return v0 + k * zz; };
  for (int i = 0; i < n; ++i) {
    const double a = f(z), b = f(z + 0.5 * h * a), c = f(z + 0.5 * h * b), d = f(z + h * c);
    z += h / 6 * (a + 2 * b + 2 * c + d);
  }
  return z;
}
VelocityModel two_layer(double t1, double t2) {
  VelocityModel m;
  m.layer_names = {"L1", "L2"};
  m.base_twt = [=](double x, double) { return std::vector<double>{t1 + 1e-5 * x, t2 + 2e-5 * x}; };
  m.velocity = [](double x, double) { return std::vector<LayerVelocity>{{1800 + 0.01 * x, 0.6}, {2500, 0.3}}; };
  return m;
}
}  // namespace

TEST_CASE("D12 layer-cake time to depth and back within 0.01 m") {
  CHECK_NEAR(layer_thickness({1800, 0.6}, 0, 0.4), rk4_depth(1800, 0.6, 0, 0.4), 1e-6);
  CHECK_NEAR(layer_thickness({2500, 0}, 300, 0.2), 500.0, 1e-9);
  CHECK_NEAR(layer_time({2500, 0.3}, 1000, 1000 + layer_thickness({2500, 0.3}, 1000, 0.25)), 0.25, 1e-9);
  auto m = two_layer(0.8, 1.5);
  double worst = 0;
  for (double x : {0.0, 2500.0, 5000.0})
    for (double t = 0.05; t < 1.9; t += 0.05) {
      const double z = time_to_depth(m, x, 0, t);
      worst = std::max(worst, std::fabs(time_to_depth(m, x, 0, depth_to_time(m, x, 0, z)) - z));
    }
  std::printf("  round trip: largest depth error %.2e m\n", worst);
  CHECK(worst < 0.01);
  const auto zb = base_depths(m, 0, 0);
  REQUIRE(zb.size() == 2);
  CHECK_NEAR(zb[0], rk4_depth(1800, 0.6, 0, 0.4), 1e-4);
  CHECK_NEAR(zb[1], rk4_depth(2500, 0.3, zb[0], 0.35), 1e-4);
}

TEST_CASE("D13 calibration recovers k and v0 from wells; residuals zero at the wells used") {
  std::vector<WellTie> w;
  const double k1 = 0.55, k2 = 0.25;
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> u(-150, 150);
  for (int i = 0; i < 5; ++i) {
    WellTie t;
    t.well = "W" + std::to_string(i + 1);
    t.x = 1000.0 * i;
    const double v01 = 1800 + u(rng), v02 = 2400 + u(rng), t1 = 0.7 + 0.05 * i, t2 = 1.3 + 0.04 * i;
    const double z1 = layer_thickness({v01, k1}, 0, t1 / 2), z2 = z1 + layer_thickness({v02, k2}, z1, (t2 - t1) / 2);
    t.top_depth = {z1, z2};
    t.top_twt = {t1, t2};
    w.push_back(t);
  }
  auto c = calibrate(w, 2);
  REQUIRE(c.k.size() == 2);
  std::printf("  k = %.5f, %.5f (true %.2f, %.2f)\n", c.k[0], c.k[1], k1, k2);
  CHECK(c.residual.size() == 2);
  double worst = 0;
  for (auto& row : c.residual)
    for (double r : row) worst = std::max(worst, std::fabs(r));
  CHECK(worst < 0.01);
  // with v0 free per well, every k ties the wells: the fitted k must still be finite and reported
  CHECK(std::isfinite(c.k[0]) && std::isfinite(c.k_sigma[0]));
  auto ex = calibrate(w, 2, {"W3"});
  CHECK(ex.excluded.size() == 1 && ex.excluded[0].find("W3") != std::string::npos);
}

TEST_CASE("D14 coordinate transformation against PROJ's own test values") {
  if (!crs_reprojection_available()) {
    std::printf("  skipped: this build has no PROJ (STRATONIX_WITH_PROJ=OFF)\n");
    return;
  }
  CrsInfo geo, utm;
  geo.epsg = 4258;   // ETRS89 geographic (GRS80)
  utm.epsg = 25832;  // ETRS89 / UTM zone 32N
  CrsTransform t(geo, utm);
  std::vector<double> x{12.0}, y{55.0}, z{0.0};
  t.apply(x, y, z);
  // PROJ gie documentation example: +proj=utm +zone=32 +ellps=GRS80, accept 12 55
  CHECK_NEAR(x[0], 691875.63214, 0.001);
  CHECK_NEAR(y[0], 6098907.82501, 0.001);
  CHECK(!t.operation().name.empty());
  CrsTransform back(utm, geo);
  back.apply(x, y, z);
  CHECK_NEAR(x[0], 12.0, 1e-8);
  CHECK_NEAR(y[0], 55.0, 1e-8);
}

TEST_CASE("D15 shapefile and ESRI ASCII grid read back") {
  io::MapLayer l;
  l.name = "units";
  io::MapFeature a;
  a.id = 1;
  a.type = 5;
  a.parts = {{{0, 0, 0}, {0, 1000, 0}, {500, 1000, 0}, {500, 0, 0}, {0, 0, 0}}};
  a.attributes["UNIT"] = "Young";
  io::MapFeature b = a;
  b.id = 2;
  b.parts = {{{500, 0, 0}, {500, 1000, 0}, {1000, 1000, 0}, {1000, 0, 0}, {500, 0, 0}}};
  b.attributes["UNIT"] = "Old";
  l.features = {a, b};
  const auto shp = tmp("d15_units.shp");
  io::write_shapefile(shp, l);
  auto r = io::read_shapefile(shp);
  REQUIRE(r.features.size() == 2);
  CHECK(r.features[1].attributes.at("UNIT") == "Old");
  CHECK(r.features[0].parts[0].size() == 5);
  const auto asc = tmp("d15_dem.asc");
  write_text(asc, "ncols 3\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 500\nNODATA_value -9999\n110 120 -9999\n100 105 108\n");
  auto g = io::read_esri_ascii_grid(asc);
  CHECK(g.spec.nx == 3 && g.spec.ny == 2);
  CHECK_NEAR(g.at(0, 0), 100, 0);   // the bottom row of the file is the southern row
  CHECK_NEAR(g.at(1, 1), 120, 0);
  CHECK(std::isnan(g.at(2, 1)));
}

TEST_CASE("D16 a map gives contact points on the elevation model and orientations") {
  io::MapLayer l;
  for (int f = 0; f < 2; ++f) {
    io::MapFeature m;
    m.id = f + 1;
    m.type = 5;
    const double x0 = 500.0 * f;
    m.parts = {{{x0, 0, 0}, {x0, 1000, 0}, {x0 + 500, 1000, 0}, {x0 + 500, 0, 0}, {x0, 0, 0}}};
    m.attributes["UNIT"] = f == 0 ? "Young" : "Old";
    l.features.push_back(m);
  }
  Grid2D dem(GridSpec{0, 0, 50, 50, 21, 21}, 100.0);
  std::vector<io::Measurement> ms{{{250, 500, 100}, deg(20), deg(270), "bedding", "deg", "field book", "2026-05-01"}};
  auto d = io::map_to_model_data(l, "UNIT", {"Young", "Old"}, nullptr, ms, dem);
  REQUIRE(!d.contact_points.empty());
  for (auto& [h, p] : d.contact_points) {
    CHECK(h == "Young");  // the contact is the base of the younger unit
    CHECK(std::fabs(p.x - 500) <= 50);
    CHECK_NEAR(p.z, 100, 1e-9);
  }
  REQUIRE(d.orientations.size() == 1);
  CHECK(d.orientations[0].source == "field book");
}

// ---- 3. Interpretation uncertainty (Prompt 57)

TEST_CASE("D17 velocity uncertainty: depth spread of one layer against the analytical value") {
  const LayerVelocity v{2000, 0.5};
  const double dt = 0.4;  // one-way
  VelocityUncertainty u;
  u.v0_sigma = 60;
  auto s = depth_spread(v, u, 0, 2 * dt, 2000, 11);
  const double expect = u.v0_sigma * (std::exp(v.k * dt) - 1) / v.k;
  std::printf("  sd %.3f m, analytical %.3f m\n", s.velocity_sd, expect);
  CHECK(std::fabs(s.velocity_sd - expect) / expect < 0.03);
  CHECK(s.members == 2000);
}

TEST_CASE("D18 velocity and pick uncertainty reported separately") {
  const LayerVelocity v{2000, 0.5};
  VelocityUncertainty u;
  u.v0_sigma = 60;
  u.pick_sigma_ms = 8;
  auto s = depth_spread(v, u, 0, 0.8, 4000, 12);
  const double zb = layer_thickness(v, 0, 0.4), vb = v.v0 + v.k * zb;
  const double pick = vb * (u.pick_sigma_ms / 1000.0) / 2;  // two-way ms to one-way s
  std::printf("  velocity %.2f m, picks %.2f m (analytical %.2f), total %.2f m\n", s.velocity_sd, s.pick_sd, pick, s.total_sd);
  CHECK(std::fabs(s.pick_sd - pick) / pick < 0.05);
  CHECK(std::fabs(s.total_sd * s.total_sd - (s.velocity_sd * s.velocity_sd + s.pick_sd * s.pick_sd)) /
            (s.total_sd * s.total_sd) < 0.05);
}

TEST_CASE("D17b ensemble: depth spread from a known velocity error against the analytical value") {
  // one layer, v0 stated (no well ties) with sd 60 m/s over the whole model, base at 800 ms two-way
  const fs::path dir = tmp("d17b");
  fs::create_directories(dir);
  std::string irap = "-996 21 100 100\n0 2000 0 2000\n21 0 0 0\n0 0 0 0 0 0 0\n";
  for (int k = 0; k < 21 * 21; ++k) irap += "800" + std::string((k % 6 == 5) ? "\n" : " ");
  write_text(dir / "h1.irap", irap + "\n");
  const std::string yaml = std::string("stratonix_schema: 1\nmodel: {name: t, version: \"1.0.0\", status: DESIGN}\n") +
      "crs: {name: local, epsg: 0, vertical_datum: datum, length_unit: m}\n"
      "domain: {origin: \"[0, 0] m | user_instruction | t\", size: \"[2000, 2000] m | user_instruction | t\", "
      "depth_range: \"[0, 2500] m | user_instruction | t\", grid_spacing: \"100 m | user_instruction | t\"}\n"
      "inputs: [{id: H1_T, path: h1.irap, role: horizon in time, sharing: shareable, source: t, evidence: interpreted}]\n"
      "stratigraphy: {reference: {point: \"[1000, 1000] m | user_instruction | t\"}, horizons: [{name: H1, surface: {input: H1_T, format: irap_ascii, z_is: twt_ms}}, "
      "{name: H2, thickness_above: \"200 m | user_instruction | t\"}], units: [{name: U, top: H1, base: H2}]}\n"
      "structure: {regional: {type: flat}}\n"
      "velocity_model: {layers: [{name: L1, base: H1, v0: \"2000 m/s | assumed | t\", k: \"0.5 1/s | assumed | t\"}]}\n"
      "uncertainty: {members: \"2000 ordinal | user_instruction | t\", seed: \"17 ordinal | user_instruction | t\", interpretation: {velocity: [{layer: L1, v0_sigma: \"60 m/s | assumed | t\"}]}}\n";
  desc::Description d;
  auto rep = desc::validate_description_text(yaml, dir, &d);
  for (auto& i : rep.issues) std::printf("  %s %s %s\n", i.rule.c_str(), i.path.c_str(), i.message.c_str());
  REQUIRE(rep.ok());
  const auto spec = unc::ensemble_spec(d);
  CHECK(spec.interpretation.size() == 1);
  const auto members = unc::sample_members(spec);
  const desc::Node base = desc::parse_document(yaml);
  std::vector<double> z;
  for (auto& m : members) {
    if (m.index == 0) continue;
    desc::Description dm;
    auto r = desc::validate_description_text(desc::to_json(unc::member_document(base, spec, m)), dir, &dm);
    REQUIRE(r.ok());
    auto vb = io::velocity_build(dm);
    z.push_back(vb.depth_grids.at("H1_T").sample(1000, 1000));
  }
  double mu = 0, sd = 0;
  for (double x : z) mu += x;
  mu /= static_cast<double>(z.size());
  for (double x : z) sd += (x - mu) * (x - mu);
  sd = std::sqrt(sd / static_cast<double>(z.size() - 1));
  const double t = 0.4, expect = 60.0 * std::expm1(0.5 * t) / 0.5;  // one-way 0.4 s
  std::printf("  %zu members: sd %.3f m, analytical %.3f m, mean %.2f m (base %.2f m)\n", z.size(), sd, expect, mu, 2000.0 * std::expm1(0.5 * t) / 0.5);
  CHECK(z.size() == 2000);
  CHECK(std::fabs(sd - expect) / expect < 0.03);  // STX-TN-16 section 3: 3 % at 2000 members
}

TEST_CASE("D17c correlated error: unit variance, correlation at L, zero at conditioning wells") {
  // covariance exp(-(r/L)^2): at r = L the correlation is 1/e
  const std::size_t n = 120;
  const double h = 50, L = 1000;
  double v0 = 0, cL = 0;
  std::size_t cnt = 0;
  for (std::uint64_t s = 1; s <= 40; ++s) {
    auto f = correlated_field(n, n, h, h, 2.0, L, s);
    for (std::size_t j = 10; j < n - 30; j += 7)
      for (std::size_t i = 10; i < n - 30; i += 7) {
        const double a = f[j * n + i], b = f[j * n + i + 20];  // 20 nodes = 1000 m
        v0 += a * a;
        cL += a * b;
        ++cnt;
      }
  }
  v0 /= static_cast<double>(cnt);
  cL /= static_cast<double>(cnt);
  std::printf("  variance %.3f (expect 4), correlation at L %.3f (expect %.3f)\n", v0, cL / v0, std::exp(-1.0));
  CHECK(std::fabs(v0 - 4.0) / 4.0 < 0.1);
  CHECK(std::fabs(cL / v0 - std::exp(-1.0)) < 0.06);
  // conditioning: zero at the wells, unchanged far from them
  auto f = correlated_field(n, n, h, h, 2.0, L, 99), g = f;
  const std::vector<std::pair<double, double>> pts = {{30, 30}, {60.5, 45.25}, {90, 90}};
  condition_field(g, n, n, h, h, L, pts);
  for (auto& [u, v] : pts) {
    const std::size_t i = static_cast<std::size_t>(u), j = static_cast<std::size_t>(v);
    CHECK(std::fabs(g[j * n + i]) < 0.3);  // the node beside the well
  }
  CHECK(std::fabs(g[30 * n + 30]) < 1e-2);  // a well on a node
  CHECK(std::fabs(g[110 * n + 5] - f[110 * n + 5]) < 1e-3);  // 5 L away
  // whole-model error: conditioning by any well removes it
  auto c = correlated_field(n, n, h, h, 2.0, 0.0, 5);
  CHECK(std::fabs(c[0] - c[n * n - 1]) < 1e-12);
  condition_field(c, n, n, h, h, 0.0, {{40, 40}});
  CHECK(std::fabs(c[7]) < 1e-4);
  // the stream is the same on every platform for a seed
  CHECK(std::fabs(seeded_normal(1, 0) - seeded_normal(1, 0)) == 0);
}

STX_TEST_MAIN()
