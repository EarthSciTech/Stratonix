// Stage 4C (STX-TN-16 section 2.3, decision 23): layer-cake depth conversion of the description's
// time surfaces, calibrated to the wells' tops and times. Every converted surface keeps its time
// grid; the velocity model and its calibration are reported.
#include "stratonix/io/velocity.hpp"

#include <algorithm>
#include <cmath>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/registry/registry.hpp"

namespace stratonix::io {

VelocityBuild velocity_build(const desc::Description& d) {
  VelocityBuild out;
  if (!d.velocity.present) return out;
  out.present = true;
  const auto& V = d.velocity;
  std::vector<Grid2D> base(V.layers.size());
  std::vector<std::string> ids(V.layers.size());
  for (std::size_t L = 0; L < V.layers.size(); ++L) {
    out.layers.push_back(V.layers[L].name);
    const auto& h = d.horizons[static_cast<std::size_t>(d.horizon_index(V.layers[L].base))];
    const desc::InputDesc* in = nullptr;
    for (auto& x : d.inputs)
      if (x.id == h.surface_input) in = &x;
    if (!in) throw Error(ErrorCode::Build, "input '" + h.surface_input + "' is not listed under inputs");
    auto data = read_input(h.surface_format, in->resolved);
    if (!data.grid) throw Error(ErrorCode::Build, "input '" + h.surface_input + "' did not contain a grid");
    Grid2D g = *data.grid;
    const double sc = h.surface_z_is == "twt_s" ? 1.0 : 1e-3;  // two-way time in seconds
    for (auto& v : g.v) v = std::isfinite(v) ? v * sc : v;
    // section 3: one realisation of the pick uncertainty (ensemble members)
    if (h.surface_perturbation.present) {
      const auto& q = h.surface_perturbation;
      const auto e = correlated_field(g.spec.nx, g.spec.ny, std::fabs(g.spec.dx), std::fabs(g.spec.dy), q.sigma.si(),
                                      q.correlation_length.present ? q.correlation_length.si() : 0.0, q.seed);
      for (std::size_t n = 0; n < g.v.size(); ++n)
        if (std::isfinite(g.v[n])) g.v[n] += e[n];
    }
    base[L] = g;
    ids[L] = h.surface_input;
    out.time_grids[h.surface_input] = g;
  }
  // well ties: tops (TVD below the datum) and their two-way times, at the path's position
  std::vector<WellTie> ties;
  for (auto& w : d.wells) {
    if (w.tops_twt.empty()) continue;
    const WellPath path = make_well_path(w, d);
    WellTie t;
    t.well = w.name;
    t.top_depth.assign(V.layers.size(), kNaN);
    t.top_twt.assign(V.layers.size(), kNaN);
    for (std::size_t L = 0; L < V.layers.size(); ++L) {
      for (auto& [h, p] : w.tops)
        if (h == V.layers[L].base) t.top_depth[L] = p.si();
      for (auto& [h, p] : w.tops_twt)
        if (h == V.layers[L].base) t.top_twt[L] = p.si();
    }
    const double z0 = std::isfinite(t.top_depth[0]) ? t.top_depth[0] : 0.0;
    const Vec3 p = path.at_elevation(-z0);
    t.x = p.x;
    t.y = p.y;
    ties.push_back(t);
    out.wells.push_back(w.name);
  }
  std::vector<std::string> excl;
  for (auto& [w, r] : V.exclude) excl.push_back(w);
  std::vector<double> fixed;
  for (auto& L : V.layers) fixed.push_back(L.k.present && !L.fit_k ? L.k.si() : kNaN);
  out.cal = calibrate(ties, V.layers.size(), excl, fixed);
  // v0 at any point: from the wells (inverse distance squared or their mean), else as stated
  std::vector<std::vector<std::pair<Vec3, double>>> v0s(V.layers.size());
  for (std::size_t L = 0; L < V.layers.size(); ++L)
    for (std::size_t w = 0; w < ties.size(); ++w)
      if (std::find(excl.begin(), excl.end(), ties[w].well) == excl.end() && std::isfinite(out.cal.v0[L][w]))
        v0s[L].push_back({{ties[w].x, ties[w].y, 0}, out.cal.v0[L][w]});
  const bool idw = V.interpolation != "constant";
  std::vector<double> stated, ks = out.cal.k;
  for (auto& L : V.layers) stated.push_back(L.v0.si());
  // section 3: one realisation of the velocity uncertainty (ensemble members): k shifted by one
  // draw; v0 by a correlated field on the layer base's lattice, conditioned to zero at the wells
  // that calibrate it (the wells tie whatever the error between them)
  std::vector<Grid2D> dv0(V.layers.size());
  for (std::size_t L = 0; L < V.layers.size(); ++L) {
    const auto& q = V.layers[L].perturbation;
    if (!q.present) continue;
    if (q.k_sigma.present) ks[L] += q.k_sigma.si() * seeded_normal(q.seed, 1);
    if (!q.v0_sigma.present) continue;
    const GridSpec& gs = base[L].spec;
    const double Lc = q.correlation_length.present ? q.correlation_length.si() : 0.0;
    Grid2D f(gs, 0.0);
    f.v = correlated_field(gs.nx, gs.ny, std::fabs(gs.dx), std::fabs(gs.dy), q.v0_sigma.si(), Lc, q.seed);
    std::vector<std::pair<double, double>> pts;
    for (auto& [p, val] : v0s[L]) pts.push_back({(p.x - gs.x0) / gs.dx, (p.y - gs.y0) / gs.dy});
    condition_field(f.v, gs.nx, gs.ny, std::fabs(gs.dx), std::fabs(gs.dy), Lc, pts);
    dv0[L] = std::move(f);
  }
  out.model.layer_names = out.layers;
  out.model.base_twt = [base](double x, double y) {
    std::vector<double> t;
    for (auto& g : base) t.push_back(g.sample(x, y));
    return t;
  };
  out.model.velocity = [v0s, stated, ks, idw, dv0](double x, double y) {
    std::vector<LayerVelocity> v;
    for (std::size_t L = 0; L < v0s.size(); ++L) {
      double v0 = stated[L];
      if (!v0s[L].empty()) {
        double num = 0, den = 0;
        for (auto& [p, val] : v0s[L]) {
          const double wgt = idw ? 1.0 / std::max(1.0, (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y)) : 1.0;
          num += wgt * val;
          den += wgt;
        }
        v0 = num / den;
      }
      if (!dv0[L].v.empty()) {
        const double e = dv0[L].sample(x, y);
        if (std::isfinite(e)) v0 += e;
      }
      v.push_back({v0, ks[L]});
    }
    return v;
  };
  // depth grids on the time grids' nodes
  for (std::size_t L = 0; L < V.layers.size(); ++L) {
    Grid2D g = base[L];
    for (std::size_t j = 0; j < g.spec.ny; ++j)
      for (std::size_t i = 0; i < g.spec.nx; ++i) {
        double& z = g.at(i, j);
        if (!std::isfinite(z)) continue;
        const auto zb = base_depths(out.model, g.spec.x(i), g.spec.y(j));
        z = L < zb.size() ? zb[L] : kNaN;
      }
    out.depth_grids[ids[L]] = g;
  }
  // ties after conversion: the model's depth at each well's top against the top
  std::string csv = "well,layer,base,top_depth_m,top_twt_ms,v0_m_per_s,k_per_s,tie_residual_m,leave_one_out_residual_m,excluded\n";
  for (std::size_t w = 0; w < ties.size(); ++w) {
    const bool ex = std::find(excl.begin(), excl.end(), ties[w].well) != excl.end();
    for (std::size_t L = 0; L < V.layers.size(); ++L) {
      if (!std::isfinite(ties[w].top_depth[L]) || !std::isfinite(ties[w].top_twt[L])) continue;
      const double zm = time_to_depth(out.model, ties[w].x, ties[w].y, ties[w].top_twt[L]);
      VelocityBuild::Tie t{ties[w].well, V.layers[L].name, zm - ties[w].top_depth[L], out.cal.loo_residual[L][w], ex};
      out.ties.push_back(t);
      csv += ties[w].well + "," + V.layers[L].name + "," + V.layers[L].base + "," + fmt_fixed(ties[w].top_depth[L], 2) + "," +
             fmt_fixed(1e3 * ties[w].top_twt[L], 2) + "," + fmt_fixed(out.cal.v0[L][w], 2) + "," + fmt_fixed(out.cal.k[L], 5) + "," +
             fmt_fixed(t.residual, 3) + "," + (std::isfinite(t.loo) ? fmt_fixed(t.loo, 3) : std::string()) + "," + (ex ? "yes" : "no") + "\n";
    }
  }
  out.csv = csv;
  std::string rep = "Stratonix velocity model (stage 4C, STX-TN-16 section 2.3; decision 23)\n";
  rep += "Layer-cake, v(z) = v0 + k z per layer (z TVD below the model datum, time zero at the datum); v0 between wells by " +
         std::string(idw ? "inverse distance squared" : "the wells' mean") + "\n";
  for (std::size_t L = 0; L < V.layers.size(); ++L) {
    rep += "  " + V.layers[L].name + " (base " + V.layers[L].base + "): k = " + fmt_fixed(out.cal.k[L], 5) + " 1/s" +
           (std::isfinite(fixed[L]) ? " (stated)" : " (fitted, sigma " + fmt_fixed(out.cal.k_sigma[L], 5) + ")") + "; v0 " +
           (v0s[L].empty() ? "as stated, " + fmt_fixed(stated[L], 1) + " m/s (no well ties)" : "from " + std::to_string(v0s[L].size()) + " wells") + "\n";
  }
  for (auto& [w, r] : V.exclude) rep += "  well " + w + " excluded from the calibration: " + r + " (a blind check)\n";
  for (auto& L : V.layers)
    if (L.perturbation.present)
      rep += "  realisation of " + L.name + "'s velocity uncertainty (seed " + std::to_string(L.perturbation.seed) + "): v0 sigma " +
             (L.perturbation.v0_sigma.present ? fmt_fixed(L.perturbation.v0_sigma.si(), 1) + " m/s" : std::string("none")) + ", k sigma " +
             (L.perturbation.k_sigma.present ? fmt_fixed(L.perturbation.k_sigma.si(), 4) + " 1/s" : std::string("none")) + ", correlation length " +
             (L.perturbation.correlation_length.present ? fmt_fixed(L.perturbation.correlation_length.si(), 0) + " m" : std::string("whole model")) + "\n";
  for (auto& h : d.horizons)
    if (h.surface_perturbation.present)
      rep += "  realisation of " + h.name + "'s pick uncertainty (seed " + std::to_string(h.surface_perturbation.seed) + "): sigma " +
             fmt_fixed(1e3 * h.surface_perturbation.sigma.si(), 1) + " ms two-way, correlation length " +
             (h.surface_perturbation.correlation_length.present ? fmt_fixed(h.surface_perturbation.correlation_length.si(), 0) + " m" : std::string("whole surface")) + "\n";
  out.report = rep;
  return out;
}

double depth_to_time_safe(const VelocityModel& m, double x, double y, double z) {
  for (double t : m.base_twt(x, y))
    if (!std::isfinite(t)) return kNaN;
  try {
    return depth_to_time(m, x, y, z);
  } catch (const Error&) {
    return kNaN;
  }
}

}  // namespace stratonix::io
