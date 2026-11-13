// Benchmarks: full build and kernels at production size, reference and optimised.
#include "stratonix/package/bench.hpp"

#include <chrono>
#include <random>
#include <sstream>

#include "stratonix/checks/checks.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

namespace stratonix {

std::string production_model_yaml(double spacing) {
  const std::string S = " | user_instruction | synthetic benchmark\"";
  auto q = [&](const std::string& v) { return "\"" + v + S; };
  std::ostringstream o;
  o << "# Synthetic production-size benchmark: 40 x 40 km, six horizons, two folds, ten faults.\n"
    << "stratonix_schema: 1\nmodel: {name: synthetic_production_10_faults, version: \"1.0.0\", status: DESIGN}\n"
    << "crs: {name: \"local engineering frame\", epsg: 0, vertical_datum: \"model datum, z = 0\", length_unit: m}\n"
    << "domain:\n  origin: " << q("[0, 0] m") << "\n  size: " << q("[40000, 40000] m") << "\n  depth_range: "
    << q("[0, 4000] m") << "\n  grid_spacing: " << q(fmt_shortest(spacing) + " m") << "\n";
  o << "wells:\n  - name: W1\n    role: control\n    head: " << q("[8000, 8000] m") << "\n    tops:\n";
  const int depths[6] = {1500, 1700, 1900, 2150, 2400, 2800};
  for (int i = 0; i < 6; ++i) o << "      H" << i + 1 << ": " << q(std::to_string(depths[i]) + " m") << "\n";
  o << "stratigraphy:\n  reference: {well: W1}\n  horizons: [{name: H1}, {name: H2}, {name: H3}, {name: H4}, {name: H5}, {name: H6}]\n  units:\n";
  for (int i = 1; i < 6; ++i) o << "    - {name: U" << i << ", top: H" << i << ", base: H" << i + 1 << "}\n";
  o << "structure:\n  regional: {type: flat}\n  base: {type: flat}\n  primitives:\n"
    << "    - {type: anticline, name: A1, crest: " << q("[20000, 20000] m") << ", axis_azimuth: " << q("40 deg")
    << ", half_length: " << q("15000 m") << ", half_width_left: " << q("6000 m") << ", half_width_right: " << q("9000 m")
    << ", relief: {H1: " << q("250 m") << ", H5: " << q("150 m") << "}}\n"
    << "    - {type: dome, name: D1, crest: " << q("[30000, 10000] m") << ", radius_major: " << q("5000 m")
    << ", radius_minor: " << q("3500 m") << ", axis_azimuth: " << q("110 deg") << ", relief: {H1: " << q("120 m") << "}}\n";
  o << "faults:\n";
  for (int k = 0; k < 10; ++k) {
    int x = 6000 + 3200 * k, y = 12000 + (k % 3) * 8000, strike = k % 2 == 0 ? 35 : 215;
    std::string upper = k % 3 ? "{rule: within_unit, unit: U1, fraction: " + q("0.5 fraction") + "}"
                              : "{rule: horizon, horizon: H2, offset: " + q("-30 m") + "}";
    char nm[8];
    std::snprintf(nm, sizeof nm, "F%02d", k + 1);
    o << "  - name: " << nm << "\n    age: " << q(std::to_string(k + 1) + " ordinal") << "\n"
      << "    geometry: {type: planar, reference_point: " << q("[" + std::to_string(x) + ", " + std::to_string(y) + "] m")
      << ", reference_depth: " << q("2000 m") << ", strike: " << q(std::to_string(strike) + " deg") << ", dip: " << q("62 deg") << "}\n"
      << "    extent: {upper_tip: " << upper << ", lower_tip: {rule: within_unit, unit: U5, fraction: " << q("0.6 fraction")
      << "}, length: " << q(std::to_string(9000 + 1000 * (k % 4)) + " m") << ", tip_line_shape: elliptical}\n"
      << "    throw: {max: " << q(std::to_string(30 + 8 * (k % 5)) + " m") << ", profile_dip: {type: flat_topped, taper: "
      << q("0.3 fraction") << "}, profile_strike: {type: elliptical}}\n"
      << "    slip_mode: dip_slip\n    hanging_wall: {type: decay, influence_distance: " << q("2500 m") << "}\n";
  }
  o << "frame: {fixed: {well: W1}}\n";
  return o.str();
}

namespace {
struct Timer {
  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  double s() const { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }
};
}  // namespace

std::vector<BenchRow> run_benchmark(const BenchOptions& o) {
  register_all_plugins();
  if (o.threads) set_thread_count(o.threads);
  std::vector<BenchRow> rows;
  std::vector<std::pair<std::string, desc::Description>> models;
  if (!o.description.empty()) {
    models.push_back({o.description.filename().string(), desc::load_description(o.description)});
  } else {
    for (double sp : o.spacings) {
      desc::Description d;
      auto rep = desc::validate_description_text(production_model_yaml(sp), ".", &d);
      if (!rep.ok()) throw Error(ErrorCode::Internal, "benchmark description invalid:\n" + rep.text());
      models.push_back({"40x40km@" + fmt_shortest(sp) + "m", d});
    }
  }
  for (auto& [size, d] : models) {
    for (int mode = 0; mode < 2; ++mode) {
      bool opt = mode == 1;
      if ((opt && !o.optimised) || (!opt && !o.reference)) continue;
      std::string mname = opt ? "optimised" : "reference";
      auto row = [&](const std::string& step, double sec, const std::string& det) {
        rows.push_back({size, mname, step, sec, peak_rss_bytes() / 1048576.0, det});
      };
      BuildOptions bo;
      bo.reference_kernels = !opt;
      Timer t;
      auto bo_out = build_only(d, bo);
      auto& m = *bo_out.model;
      row("build model", t.s(), std::to_string(m.info().faults.size()) + " faults");
      const auto& dom = m.info().domain;
      GridSpec g = GridSpec::covering(dom, d.grid_spacing.si());
      // point evaluation throughput
      {
        std::mt19937_64 rng(1);
        std::uniform_real_distribution<double> ux(dom.x0, dom.x0 + dom.sx), uy(dom.y0, dom.y0 + dom.sy), uz(dom.z_base, dom.z_top);
        std::vector<Vec3> p(1000000);
        for (auto& x : p) x = {ux(rng), uy(rng), uz(rng)};
        Timer te;
        EvalResult r;
        m.evaluate(p, kEvalAll, r);
        row("evaluate 1M points", te.s(), "unit, block and 6 scalars");
      }
      {
        Timer tb;
        auto bl = repr::list_blocks(m);
        row("list fault blocks", tb.s(), std::to_string(bl.codes.size()) + " blocks");
      }
      {
        Timer ts;
        std::size_t gaps = 0;
        for (std::size_t h = 0; h < m.info().horizons.size(); ++h) gaps += repr::horizon_grid_faulted(m, static_cast<int>(h), g).gap_nodes;
        row("faulted horizon grids (all)", ts.s(), std::to_string(g.nx) + " x " + std::to_string(g.ny) + " nodes x " +
                                                       std::to_string(m.info().horizons.size()) + " horizons; " + std::to_string(gaps) + " gap nodes");
      }
      {
        Timer ts;
        for (std::size_t h = 0; h < m.info().horizons.size(); ++h) repr::horizon_grid_block(m, static_cast<int>(h), 0, g);
        row("per-block horizon grids (one block, all)", ts.s(), "block 0");
      }
      std::vector<repr::FaultSurface> fss;
      {
        Timer tf;
        std::size_t tris = 0;
        for (std::size_t f = 0; f < m.info().faults.size(); ++f) {
          fss.push_back(repr::fault_surface(m, static_cast<int>(f), d.grid_spacing.si() / 2));
          tris += fss.back().mesh.t.size();
        }
        row("fault surfaces (all)", tf.s(), std::to_string(tris) + " triangles");
      }
      {
        Timer tc;
        std::size_t n = 0;
        for (auto& fs : fss)
          for (std::size_t h = 0; h < m.info().horizons.size(); ++h) n += repr::cutoff_lines(m, fs, static_cast<int>(h)).size();
        row("cut-off lines (all faults x horizons)", tc.s(), std::to_string(n) + " lines");
      }
      {
        Timer tv;
        auto v = repr::voxel_grid(m, o.voxel_xy, o.voxel_xy, o.voxel_z, true);
        row("voxel grid", tv.s(), std::to_string(v.nx) + " x " + std::to_string(v.ny) + " x " + std::to_string(v.nz) + " cells at " +
                                      fmt_shortest(o.voxel_xy) + "/" + fmt_shortest(o.voxel_z) + " m");
      }
      {
        Timer tk;
        CheckContext ctx(m, d, make_well_paths(d), d.grid_spacing.si());
        auto rep = run_checks(ctx);
        row("checks (all)", tk.s(), std::to_string(rep.results.size()) + " checks, " + std::to_string(rep.failures()) + " failed");
      }
      row("total", t.s(), "");
    }
  }
  return rows;
}

std::string bench_table(const std::vector<BenchRow>& rows) {
  std::ostringstream s;
  s << "Stratonix benchmark (" << build_info().version << ", " << build_info().compiler << ", " << build_info().platform
    << ", " << thread_count() << " threads)\n\n";
  for (auto& r : rows)
    s << r.size << " | " << r.mode << " | " << r.step << " | " << fmt_fixed(r.seconds, 3) << " s | peak " << fmt_fixed(r.peak_rss_mb, 0)
      << " MB | " << r.detail << "\n";
  return s.str();
}
std::string bench_csv(const std::vector<BenchRow>& rows) {
  std::string s = "size,mode,step,seconds,peak_rss_mb,detail\n";
  for (auto& r : rows)
    s += r.size + "," + r.mode + ",\"" + r.step + "\"," + fmt_fixed(r.seconds, 4) + "," + fmt_fixed(r.peak_rss_mb, 1) + ",\"" + r.detail + "\"\n";
  return s;
}

}  // namespace stratonix
