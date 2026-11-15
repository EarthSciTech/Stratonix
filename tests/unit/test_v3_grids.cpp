// Stage 3B analytical tests (STX-TN-12 section 11), written before the features (Prompt 37) and
// made to pass by Prompts 38 and 39: corner-point grids conforming to faults and horizons, their
// layering, the grid checks, GRDECL export, porosity from compaction, fault transmissibility,
// missing values, tetrahedral grids, the upscaling interface and the gridded examples.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>

#include "v2_fixture.hpp"
#include "stratonix/io/export.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/repr/grid.hpp"

namespace fs = std::filesystem;
using namespace stratonix;
using namespace v2fx;
using repr::CornerPointGrid;
using repr::GridOptions;
using repr::Layering;
using repr::ZoneSpec;

namespace {

fs::path tmpdir(const std::string& name) {
  fs::path p = fs::path(stxtest::tmp_dir()) / ("v3g_" + name);
  fs::remove_all(p);
  fs::create_directories(p);
  return p;
}

// Description with free horizons and tops (as test_v2_salt), domain 10 x 10 km.
std::string doc(const std::string& horizons, const std::string& tops, const std::string& units, const std::string& structure,
                const std::string& strat_extra = "") {
  std::string s = "stratonix_schema: 1\nmodel: {name: t, version: \"1.0.0\", status: DESIGN}\n"
                  "crs: {name: local, vertical_datum: datum, length_unit: m}\n"
                  "domain:\n  origin: " + q("[0, 0] m") + "\n  size: " + q("[10000, 10000] m") +
                  "\n  depth_range: " + q("[0, 3000] m") + "\n  grid_spacing: " + q("50 m") + "\n";
  s += "wells:\n  - {name: W1, role: control, head: " + q("[1000, 1000] m") + ", tops: {" + tops + "}}\n";
  s += "stratigraphy:\n  reference: {well: W1}\n  horizons: [" + horizons + "]\n  units: [" + units + "]\n" + strat_extra;
  s += structure.empty() ? "structure:\n  regional: {type: flat}\n" : structure;
  return s;
}
std::string top(const std::string& h, double d) { return h + ": " + q(fmt_fixed(d, 1) + " m"); }
std::string units_for(const std::vector<std::string>& hs) {
  std::string u;
  for (std::size_t k = 0; k + 1 < hs.size(); ++k)
    u += (k ? ", " : "") + std::string("{name: Unit_") + std::to_string(k + 1) + ", top: " + hs[k] + ", base: " + hs[k + 1] + "}";
  return u;
}
std::string compaction(const std::string& units) {
  return "  compaction:\n    geometry: present\n    surface_depth: " + q("0 m") + "\n    units:\n" + units;
}
std::string expo(const std::string& u) {
  return "      - {unit: " + u + ", law: exponential, surface_porosity: " + q("0.5 fraction") + ", decay: " + q("0.0005 1/m") + "}\n";
}

/// Zones for every inner unit of a model (between its top and base horizons), nk layers each.
std::vector<ZoneSpec> all_zones(const Model& m, int nk) {
  std::vector<ZoneSpec> z;
  const int nh = static_cast<int>(m.info().horizons.size());
  for (int u = 1; u < nh; ++u) z.push_back({u, nk, Layering::Proportional, 0});
  return z;
}
GridOptions options(double dx, std::vector<ZoneSpec> zones, double margin = 0) {
  GridOptions o;
  o.dx = o.dy = dx;
  o.margin = margin;
  o.zones = std::move(zones);
  return o;
}
/// Thickness at corner c (0..3, top) of cell (i, j, k): top minus bottom elevation.
double corner_thickness(const CornerPointGrid& g, int i, int j, int k, int c) {
  return g.zcorn[g.zcorn_index(i, j, k, c)] - g.zcorn[g.zcorn_index(i, j, k, c + 4)];
}

const std::string kThrow100 = "max: " + q("100 m");
/// One planar normal fault, strike 0, dip 60 deg east, throw 100 m, effectively infinite.
std::string fault_doc() {
  return model_yaml(planar_fault("F1", "1", "[5000, 5000] m", "1300 m", "0 deg", "60 deg", "100000 m", kThrow100, "dip_slip"));
}

// ---- minimal GRDECL reader (test only) -------------------------------------------------------
std::map<std::string, std::vector<double>> read_grdecl(const std::string& text) {
  std::map<std::string, std::vector<double>> kw;
  std::istringstream in(text);
  std::string line, cur;
  while (std::getline(in, line)) {
    auto c = line.find("--");
    if (c != std::string::npos) line = line.substr(0, c);
    std::istringstream ls(line);
    std::string tok;
    while (ls >> tok) {
      if (cur.empty()) {
        if (std::isalpha(static_cast<unsigned char>(tok[0]))) { cur = tok; kw[cur]; }
        continue;
      }
      if (tok == "/") { cur.clear(); continue; }
      bool end = tok.back() == '/';
      if (end) tok.pop_back();
      if (!tok.empty() && cur != "FAULTS") {
        auto star = tok.find('*');
        if (star != std::string::npos) {
          int n = std::stoi(tok.substr(0, star));
          double v = std::stod(tok.substr(star + 1));
          kw[cur].insert(kw[cur].end(), static_cast<std::size_t>(n), v);
        } else if (std::isdigit(static_cast<unsigned char>(tok[0])) || tok[0] == '-' || tok[0] == '.') {
          kw[cur].push_back(std::stod(tok));
        }
      }
      if (end && cur != "FAULTS") cur.clear();
    }
  }
  return kw;
}

}  // namespace

TEST_CASE("G1 layer cake: cell volumes dx dy t / nk, zone volumes area x thickness, checks pass") {
  auto m = build_text(model_yaml(""));
  auto g = repr::corner_point_grid(*m, options(500, {{1, 3, Layering::Proportional, 0}, {2, 2, Layering::Proportional, 0},
                                                     {3, 1, Layering::Proportional, 0}}));
  REQUIRE(g.nx == 20 && g.ny == 20 && g.nz == 6);
  const int nk[] = {3, 3, 3, 2, 2, 1};
  std::map<int, double> zone_vol;
  double worst = 0;
  for (int k = 0; k < g.nz; ++k)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i) {
        auto id = g.cell(i, j, k);
        REQUIRE(g.actnum[id] == 1);
        double v = repr::cell_volume(g.corners(i, j, k));
        worst = std::max(worst, std::fabs(v - 500.0 * 500.0 * 200.0 / nk[k]));
        zone_vol[g.zone[id]] += v;
      }
  std::printf("  worst cell volume error %.3g m3\n", worst);
  CHECK(worst <= 1e-6 * 500.0 * 500.0 * 200.0);
  for (int u = 1; u <= 3; ++u) CHECK_NEAR(zone_vol[u], 1e8 * 200.0, 1e-6 * 1e8 * 200.0);
  auto d = repr::grid_defects(*m, g, 5.0);
  CHECK(d.nonpositive_volume == 0);
  CHECK(d.inverted == 0);
}

TEST_CASE("G2 fault-conforming pillars: on the fault, the throw between the sides, different blocks, C25 passes") {
  auto m = build_text(fault_doc());
  auto g = repr::corner_point_grid(*m, options(250, all_zones(*m, 2), 50));
  int snapped = 0, compared = 0;
  double worst_on = 0, worst_throw = 0;
  FaultEval fe;
  for (int j = 0; j <= g.ny; ++j)
    for (int i = 0; i <= g.nx; ++i) {
      auto p = g.pillar(i, j);
      if (g.pillar_fault[p] != 0) continue;
      ++snapped;
      std::vector<Vec3> ends = {{g.coord[6 * p], g.coord[6 * p + 1], g.coord[6 * p + 2]},
                                {g.coord[6 * p + 3], g.coord[6 * p + 4], g.coord[6 * p + 5]}};
      m->evaluate_fault(0, ends, fe);
      worst_on = std::max({worst_on, std::fabs(fe.distance[0]), std::fabs(fe.distance[1])});
      if (i == 0 || i == g.nx || j == g.ny) continue;
      // H2 is the top of zone U2 (unit 2): the first layer of that zone is k = 2
      const int k = 2;
      auto left = g.cell(i - 1, j, k), right = g.cell(i, j, k);
      CHECK(g.block[left] != g.block[right]);
      double zl = g.zcorn[g.zcorn_index(i - 1, j, k, 1)], zr = g.zcorn[g.zcorn_index(i, j, k, 0)];
      worst_throw = std::max(worst_throw, std::fabs((zl - zr) - 100.0));
      ++compared;
    }
  std::printf("  %d snapped pillars, worst distance to the fault %.3g m, worst throw error %.3g m (%d faces)\n", snapped, worst_on,
              worst_throw, compared);
  CHECK(snapped >= g.ny + 1);
  CHECK(worst_on <= 1e-6);
  CHECK(compared > 0);
  CHECK(worst_throw <= 0.1);
  auto d = repr::grid_defects(*m, g, 0.01 * 250);
  CHECK(d.fault_faces > 0);
  CHECK(d.fault_face_off == 0);
  CHECK(d.nonpositive_volume == 0 && d.inverted == 0);
}

TEST_CASE("G3 layering on a wedge: proportional, follow_top and follow_base") {
  // H2 = H1 + t(x), t = 100 + 0.02 x (100 m at x = 0 to 300 m at x = 10 km)
  std::string hs = "{name: H1}, {name: H2, isochore: {from: H1, thickness: " + q("120 m") + ", trend: {type: linear, gradient: " +
                   q("[0.02, 0] ratio") + "}}}, {name: H3}";
  auto m = build_text(doc(hs, top("H1", 900) + ", " + top("H2", 1020) + ", " + top("H3", 1500), units_for({"H1", "H2", "H3"}), ""));
  auto t_at = [](double x) { return 100.0 + 0.02 * x; };
  {
    auto g = repr::corner_point_grid(*m, options(500, {{1, 5, Layering::Proportional, 0}, {2, 1, Layering::Proportional, 0}}));
    double worst = 0;
    for (int k = 0; k < 5; ++k)
      for (int i = 0; i < g.nx; ++i) worst = std::max(worst, std::fabs(corner_thickness(g, i, 3, k, 0) - t_at(i * 500.0) / 5));
    std::printf("  proportional: worst layer error %.3g m\n", worst);
    CHECK(worst <= 1e-3);
  }
  for (auto lay : {Layering::FollowTop, Layering::FollowBase}) {
    auto g = repr::corner_point_grid(*m, options(500, {{1, 15, lay, 20.0}, {2, 1, Layering::Proportional, 0}}));
    double worst = 0;
    int inactive_ok = 0, inactive_bad = 0;
    for (int k = 0; k < 15; ++k)
      for (int i = 0; i < g.nx; ++i) {
        const int m_from = lay == Layering::FollowTop ? k : 14 - k;  // index counted from the followed boundary
        double t = t_at(i * 500.0);
        double expect = m_from < 14 ? std::min(20.0, std::max(0.0, t - 20.0 * m_from)) : std::max(0.0, t - 280.0);
        worst = std::max(worst, std::fabs(corner_thickness(g, i, 3, k, 0) - expect));
        bool all_zero = true;
        for (int c = 0; c < 4; ++c) all_zero = all_zero && corner_thickness(g, i, 3, k, c) < 1e-9;
        if (all_zero) (g.actnum[g.cell(i, 3, k)] == 0 && !g.inactive_reason[g.cell(i, 3, k)].empty() ? inactive_ok : inactive_bad)++;
      }
    std::printf("  %s: worst layer error %.3g m, %d empty cells inactive\n", repr::layering_name(lay), worst, inactive_ok);
    CHECK(worst <= 1e-3);
    CHECK(inactive_ok > 0);
    CHECK(inactive_bad == 0);
  }
}

TEST_CASE("G4 the checks catch an inverted cell") {
  CornerPointGrid g;
  g.nx = g.ny = g.nz = 1;
  g.dx = g.dy = 100;
  for (int j = 0; j <= 1; ++j)
    for (int i = 0; i <= 1; ++i) g.coord.insert(g.coord.end(), {i * 100.0, j * 100.0, 0.0, i * 100.0, j * 100.0, -1000.0});
  g.zcorn.assign(8, 0.0);
  for (int c = 0; c < 8; ++c) g.zcorn[g.zcorn_index(0, 0, 0, c)] = c < 4 ? -500.0 : -400.0;  // bottom above top
  g.actnum = {1};
  g.zone = {1};
  g.block = {0};
  g.layer_unit = {1};
  g.pillar_fault.assign(4, -1);
  g.inactive_reason = {""};
  auto m = build_text(model_yaml(""));
  auto d = repr::grid_defects(*m, g, 1.0);
  CHECK(d.inverted == 1);
  CHECK(!d.first_inverted.empty());
}

TEST_CASE("G5 GRDECL round trip: COORD, ZCORN and ACTNUM read back to the written decimals") {
  auto m = build_text(fault_doc());
  auto g = repr::corner_point_grid(*m, options(500, all_zones(*m, 2), 50));
  io::GrdeclInfo info;
  info.header = {"test grid"};
  info.fault_names = {"F1"};
  std::string text = io::write_grdecl(g, nullptr, nullptr, info, 3);
  auto dir = tmpdir("grdecl");
  std::ofstream(dir / "grid.grdecl") << text;
  auto kw = read_grdecl(text);
  REQUIRE(kw.count("SPECGRID") && kw["SPECGRID"].size() >= 3);
  CHECK(static_cast<int>(kw["SPECGRID"][0]) == g.nx && static_cast<int>(kw["SPECGRID"][1]) == g.ny &&
        static_cast<int>(kw["SPECGRID"][2]) == g.nz);
  REQUIRE(kw["COORD"].size() == g.coord.size());
  REQUIRE(kw["ZCORN"].size() == g.zcorn.size());
  REQUIRE(kw["ACTNUM"].size() == g.actnum.size());
  double wc = 0, wz = 0;
  for (std::size_t n = 0; n < g.coord.size(); ++n) {
    double expect = n % 3 == 2 ? -g.coord[n] : g.coord[n];  // depth positive down
    wc = std::max(wc, std::fabs(kw["COORD"][n] - expect));
  }
  for (std::size_t n = 0; n < g.zcorn.size(); ++n) wz = std::max(wz, std::fabs(kw["ZCORN"][n] + g.zcorn[n]));
  std::size_t act = 0;
  for (std::size_t n = 0; n < g.actnum.size(); ++n) act += static_cast<int>(kw["ACTNUM"][n]) == g.actnum[n];
  std::printf("  COORD %.3g m, ZCORN %.3g m, ACTNUM %zu of %zu equal\n", wc, wz, act, g.actnum.size());
  CHECK(wc <= 0.0005 + 1e-9);
  CHECK(wz <= 0.0005 + 1e-9);
  CHECK(act == g.actnum.size());
}

TEST_CASE("G6 porosity from a known compaction law") {
  auto m = build_text(doc("{name: H1}, {name: H2}, {name: H3}", top("H1", 1000) + ", " + top("H2", 1200) + ", " + top("H3", 1500),
                          units_for({"H1", "H2", "H3"}), "", compaction(expo("Unit_1"))));
  auto g = repr::corner_point_grid(*m, options(1000, {{1, 4, Layering::Proportional, 0}, {2, 1, Layering::Proportional, 0}}));
  repr::PropertyInputs in;
  in.porosity_from_compaction = true;
  auto ps = repr::grid_properties(*m, g, in);
  auto* poro = ps.find("PORO");
  REQUIRE(poro != nullptr);
  CHECK(!poro->evidence.empty() && !poro->source.empty());
  const double p0 = 0.5, c = 0.0005;
  double worst = 0;
  for (int k = 0; k < 4; ++k) {
    double z1 = 1000 + 50.0 * k, z2 = z1 + 50;
    double phi = 1 - ((z2 - z1) - p0 / c * (std::exp(-c * z1) - std::exp(-c * z2))) / (z2 - z1);
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i) worst = std::max(worst, std::fabs(poro->v[g.cell(i, j, k)] - phi));
  }
  std::printf("  worst porosity error %.3g\n", worst);
  CHECK(worst <= 1e-6);
}

TEST_CASE("G7 transmissibility multiplier from a known shale gouge ratio (Manzocchi et al. 1999)") {
  const double sgr = 0.3, D = 50, tf = 0.5, k = 100, L = 100;
  double kf = std::pow(10.0, -4 * sgr - 0.25 * std::log10(D) * std::pow(1 - sgr, 5));
  double T = 1.0 / (1 + tf * (2 / kf - 1 / k - 1 / k) / (L / k + L / k));
  double kf_s = repr::manzocchi_fault_permeability(sgr, D);
  double T_s = repr::manzocchi_multiplier(tf, kf_s, k, k, L, L);
  std::printf("  k_f %.6g mD (test %.6g), T %.6g (test %.6g)\n", kf_s, kf, T_s, T);
  CHECK(std::fabs(kf_s - kf) <= 1e-9 * kf);
  CHECK(std::fabs(T_s - T) <= 1e-9 * T);
}

TEST_CASE("G8 missing, not defaulted: porosity without a law, a transform of it, facies not given") {
  auto m = build_text(doc("{name: H1}, {name: H2}, {name: H3}", top("H1", 1000) + ", " + top("H2", 1200) + ", " + top("H3", 1500),
                          units_for({"H1", "H2", "H3"}), "", compaction(expo("Unit_1"))));
  auto g = repr::corner_point_grid(*m, options(1000, {{1, 2, Layering::Proportional, 0}, {2, 2, Layering::Proportional, 0}}));
  repr::PropertyInputs in;
  in.porosity_from_compaction = true;
  in.facies = {{1, 1, "sand", "assumed", "test"}};
  in.transforms = {{"PERMX", "mD", "PORO", "exponential", 0.1, 20.0, "assumed", "test"}};
  auto ps = repr::grid_properties(*m, g, in);
  auto *poro = ps.find("PORO"), *perm = ps.find("PERMX"), *fac = ps.find("FACIES");
  REQUIRE(poro && perm && fac);
  const std::size_t per_layer = static_cast<std::size_t>(g.nx) * g.ny;
  for (auto* p : {poro, perm, fac}) {
    std::printf("  %s: %zu missing (%s)\n", p->name.c_str(), p->missing, p->missing_reason.c_str());
    CHECK(p->missing == 2 * per_layer);          // the two layers of Unit_2
    CHECK(!p->missing_reason.empty());
    CHECK(std::isnan(p->v[g.cell(0, 0, 3)]));
    CHECK(!std::isnan(p->v[g.cell(0, 0, 0)]));
  }
  CHECK_NEAR(perm->v[g.cell(0, 0, 0)], 0.1 * std::exp(20.0 * poro->v[g.cell(0, 0, 0)]), 1e-9 * perm->v[g.cell(0, 0, 0)]);
}

TEST_CASE("G9 tetrahedral grid: positive tetrahedra, volumes sum to the cells', shared faces") {
  auto planar = model_yaml("", "structure:\n  regional: {type: planar, dip: " + q("5 deg") + ", dip_azimuth: " + q("90 deg") + "}\n");
  for (int pass = 0; pass < 2; ++pass) {
    auto m = build_text(pass == 0 ? planar : fault_doc());
    auto g = repr::corner_point_grid(*m, options(500, all_zones(*m, 2), 0));
    auto t = repr::tetrahedral_grid(g);
    std::vector<double> sum(g.ncell(), 0.0);
    std::size_t negative = 0;
    std::map<std::array<std::int32_t, 3>, int> faces;
    for (std::size_t n = 0; n < t.tets.size(); ++n) {
      auto& q4 = t.tets[n];
      double v = repr::tet_volume(t.nodes[q4[0]], t.nodes[q4[1]], t.nodes[q4[2]], t.nodes[q4[3]]);
      negative += v <= 0;
      sum[t.cell[n]] += v;
      for (int f = 0; f < 4; ++f) {
        std::array<std::int32_t, 3> tri;
        for (int a = 0, b = 0; a < 4; ++a)
          if (a != f) tri[b++] = q4[a];
        std::sort(tri.begin(), tri.end());
        faces[tri]++;
      }
    }
    double worst = 0;
    for (int k = 0; k < g.nz; ++k)
      for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
          auto id = g.cell(i, j, k);
          if (!g.actnum[id]) continue;
          double v = repr::cell_volume(g.corners(i, j, k));
          worst = std::max(worst, std::fabs(sum[id] - v) / v);
        }
    int over = 0, shared = 0;
    for (auto& [f, n] : faces) { over += n > 2; shared += n == 2; }
    std::printf("  %s: %zu tetrahedra, %zu not positive, worst volume difference %.3g, %d shared faces\n",
                pass == 0 ? "planar layers" : "faulted", t.tets.size(), negative, worst, shared);
    CHECK(negative == 0);
    CHECK(worst <= (pass == 0 ? 1e-9 : 1e-2));  // exact where the faces are planar
    CHECK(over == 0);
    CHECK(shared > 0);
  }
}

TEST_CASE("G10 upscaling interface: volume-weighted mean of uniform and linear properties") {
  auto m = build_text(model_yaml(""));
  auto fine = repr::corner_point_grid(*m, options(250, all_zones(*m, 2)));
  auto coarse = repr::corner_point_grid(*m, options(500, all_zones(*m, 1)));
  repr::CellProperty uni{"U", "fraction", "assumed", "test"}, lin{"X", "m", "assumed", "test"};
  uni.v.assign(fine.ncell(), 0.2);
  lin.v.resize(fine.ncell());
  for (int k = 0; k < fine.nz; ++k)
    for (int j = 0; j < fine.ny; ++j)
      for (int i = 0; i < fine.nx; ++i) lin.v[fine.cell(i, j, k)] = repr::cell_centre(fine.corners(i, j, k)).x;
  auto cu = repr::upscale_volume_weighted_mean(fine, uni, coarse);
  auto cl = repr::upscale_volume_weighted_mean(fine, lin, coarse);
  double wu = 0, wl = 0;
  for (int k = 0; k < coarse.nz; ++k)
    for (int j = 0; j < coarse.ny; ++j)
      for (int i = 0; i < coarse.nx; ++i) {
        auto id = coarse.cell(i, j, k);
        wu = std::max(wu, std::fabs(cu.v[id] - 0.2));
        wl = std::max(wl, std::fabs(cl.v[id] - repr::cell_centre(coarse.corners(i, j, k)).x));
      }
  CHECK(wu <= 1e-12);
  CHECK(wl <= 1e-6);
  CHECK(cu.evidence == "upscaled" && !cu.source.empty());
}

TEST_CASE("G11 the stage 2A fault-system and stage 2E inversion examples, gridded, pass C23 to C25") {
  register_all_plugins();
  for (const char* ex : {"fault_system_relay_splay", "rift_inversion_history"}) {
    auto d = desc::load_description(fs::path(STX_SOURCE_DIR) / "examples" / (std::string(ex) + ".yaml"));
    BuildInputs ins;
    io::load_build_inputs(d, ins);
    BuildLog log;
    auto m = build_model(d, ins, log);
    auto g = repr::corner_point_grid(*m, options(200, all_zones(*m, 2), 50));
    std::size_t act = 0;
    for (auto a : g.actnum) act += a;
    auto dfx = repr::grid_defects(*m, g, 0.01 * 200);
    std::printf("  %s: %d x %d x %d cells, %zu active; %zu fault faces; defects %zu / %zu / %zu (worst fault distance %.3g m)\n", ex,
                g.nx, g.ny, g.nz, act, dfx.fault_faces, dfx.nonpositive_volume, dfx.inverted, dfx.fault_face_off,
                dfx.worst_fault_distance);
    CHECK(act > 0);
    CHECK(dfx.fault_faces > 0);
    CHECK(dfx.nonpositive_volume == 0);
    CHECK(dfx.inverted == 0);
    CHECK(dfx.fault_face_off == 0);
  }
}

TEST_CASE("G12 RESQML export through FESAPI (read back by resqpy in the CI job 'resqml')") {
#ifdef STRATONIX_HAVE_FESAPI
  auto dir = tmpdir("resqml");
  BuildOptions o;
  auto r = run_export(fs::path(STX_SOURCE_DIR) / "examples" / "fault_system_relay_splay.yaml", dir / "out", {"grdecl", "resqml"}, 0, o);
  CHECK(r.checks_passed);
  CHECK(fs::exists(dir / "out" / "exports" / "grid" / "example_fault_system_grid.epc"));
  CHECK(fs::exists(dir / "out" / "exports" / "grid" / "example_fault_system_grid.h5"));
#else
  std::printf("  skipped: this build has no FESAPI (the RESQML case runs in the CI job 'resqml')\n");
#endif
}

TEST_CASE("G13 description grid block: valid, and rule F32 refuses an unknown unit and zero layers") {
  const std::string grid_ok = "grid:\n  spacing: " + q("[500, 500] m") + "\n  zones:\n    - {unit: U1, layers: " + q("2 ordinal") +
                              ", layering: proportional}\n";
  auto ok = validate(model_yaml("", "", {1000, 1200, 1400, 1600}, grid_ok));
  if (!ok.ok()) std::cerr << ok.text();
  CHECK(ok.ok());
  const std::string grid_bad = "grid:\n  spacing: " + q("[500, 500] m") + "\n  zones:\n    - {unit: NOPE, layers: " + q("0 ordinal") +
                               ", layering: proportional}\n";
  auto bad = validate(model_yaml("", "", {1000, 1200, 1400, 1600}, grid_bad));
  CHECK(!bad.ok());
  CHECK(bad.text().find("F32") != std::string::npos);
}

STX_TEST_MAIN()
