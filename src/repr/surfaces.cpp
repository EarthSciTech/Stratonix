// Horizon grids (per block and faulted with gaps) and block listing, by evaluation only.
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include "repr_internal.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/repr/repr.hpp"

namespace stratonix::repr {

BlockList list_blocks(const Model& m, std::size_t n, std::size_t nz) {
  const Domain& d = m.info().domain;
  std::vector<Vec3> pts;
  pts.reserve(n * n * nz);
  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t i = 0; i < n; ++i)
        pts.push_back({d.x0 + d.sx * (i + 0.5) / n, d.y0 + d.sy * (j + 0.5) / n, d.z_base + (d.z_top - d.z_base) * (k + 0.5) / nz});
  EvalResult r;
  m.evaluate(pts, kEvalBlock, r);
  std::map<BlockCode, std::size_t> c;
  for (auto b : r.block) c[b]++;
  BlockList bl;
  for (auto& [k, v] : c) { bl.codes.push_back(k); bl.counts.push_back(v); }
  return bl;
}

// ---- vertical roots ------------------------------------------------------------------------
void vertical_roots(const Model& m, int h, BlockCode block, const std::vector<double>& xs, const std::vector<double>& ys,
                    double z_hi, double z_lo, std::vector<double>& z_out, const std::vector<double>* guess) {
  const std::size_t n = xs.size();
  z_out.assign(n, kNaN);
  if (n == 0) return;
  std::vector<double> hi(n, z_hi), lo(n, z_lo), fhi(n), flo(n);
  std::vector<Vec3> p(n);
  if (guess && optimised_kernels()) {
    // Narrow brackets around the guesses; nodes whose narrow bracket fails use the full one.
    const double w = 40;
    std::vector<std::size_t> g;
    for (std::size_t i = 0; i < n; ++i)
      if (!std::isnan((*guess)[i])) g.push_back(i);
    std::vector<Vec3> ph(g.size()), pl(g.size());
    for (std::size_t k = 0; k < g.size(); ++k) {
      std::size_t i = g[k];
      double gz = std::clamp((*guess)[i], z_lo + w, z_hi - w);
      ph[k] = {xs[i], ys[i], gz + w};
      pl[k] = {xs[i], ys[i], gz - w};
    }
    std::vector<double> sh(g.size()), sl(g.size());
    m.scalar_in_block(ph, block, h, sh);
    m.scalar_in_block(pl, block, h, sl);
    std::vector<char> narrowed(n, 0);
    for (std::size_t k = 0; k < g.size(); ++k)
      if (sh[k] > 0 && sl[k] < 0) {
        std::size_t i = g[k];
        hi[i] = ph[k].z; fhi[i] = sh[k];
        lo[i] = pl[k].z; flo[i] = sl[k];
        narrowed[i] = 1;
      }
    std::vector<std::size_t> rest;
    for (std::size_t i = 0; i < n; ++i)
      if (!narrowed[i]) rest.push_back(i);
    std::vector<Vec3> q1(rest.size()), q2(rest.size());
    for (std::size_t k = 0; k < rest.size(); ++k) {
      q1[k] = {xs[rest[k]], ys[rest[k]], z_hi};
      q2[k] = {xs[rest[k]], ys[rest[k]], z_lo};
    }
    std::vector<double> s1(rest.size()), s2(rest.size());
    m.scalar_in_block(q1, block, h, s1);
    m.scalar_in_block(q2, block, h, s2);
    for (std::size_t k = 0; k < rest.size(); ++k) { fhi[rest[k]] = s1[k]; flo[rest[k]] = s2[k]; }
  } else {
    for (std::size_t i = 0; i < n; ++i) p[i] = {xs[i], ys[i], z_hi};
    m.scalar_in_block(p, block, h, fhi);
    for (std::size_t i = 0; i < n; ++i) p[i].z = z_lo;
    m.scalar_in_block(p, block, h, flo);
  }
  std::vector<std::size_t> act;
  for (std::size_t i = 0; i < n; ++i) {
    if (fhi[i] == 0) { z_out[i] = hi[i]; continue; }
    if (flo[i] == 0) { z_out[i] = lo[i]; continue; }
    if (fhi[i] > 0 && flo[i] < 0) act.push_back(i);
  }
  const bool fast = optimised_kernels();
  std::vector<Vec3> q;
  std::vector<double> f;
  std::vector<int> side(n, 0);  // Illinois bookkeeping
  for (int it = 0; it < (fast ? 60 : 64) && !act.empty(); ++it) {
    q.resize(act.size());
    for (std::size_t k = 0; k < act.size(); ++k) {
      std::size_t i = act[k];
      double z;
      if (fast) z = hi[i] - fhi[i] * (hi[i] - lo[i]) / (fhi[i] - flo[i]);  // regula falsi
      else z = 0.5 * (hi[i] + lo[i]);                                         // bisection
      if (!(z < hi[i] && z > lo[i])) z = 0.5 * (hi[i] + lo[i]);
      q[k] = {xs[i], ys[i], z};
    }
    f.resize(act.size());
    m.scalar_in_block(q, block, h, f);
    std::vector<std::size_t> next;
    for (std::size_t k = 0; k < act.size(); ++k) {
      std::size_t i = act[k];
      double z = q[k].z, fz = f[k];
      if (fz > 0) {
        hi[i] = z; fhi[i] = fz;
        if (fast && side[i] == 1) flo[i] *= 0.5;
        side[i] = 1;
      } else {
        lo[i] = z; flo[i] = fz;
        if (fast && side[i] == -1) fhi[i] *= 0.5;
        side[i] = -1;
      }
      bool done = std::fabs(fz) < 1e-9 || (hi[i] - lo[i]) < (fast ? 1e-7 : 1e-8);
      if (done) z_out[i] = fz == 0 ? z : (std::fabs(fz) < 1e-9 ? z : 0.5 * (hi[i] + lo[i]));
      else next.push_back(i);
    }
    act.swap(next);
  }
  for (auto i : act) z_out[i] = 0.5 * (hi[i] + lo[i]);
}

/// Bilinear guesses on grid g from a coarse grid with stride `st` (NaN where unavailable).
static std::vector<double> coarse_guess(const GridSpec& g, std::size_t st, const std::vector<double>& cz, const GridSpec& c) {
  std::vector<double> out(g.size(), kNaN);
  Grid2D cg(c);
  cg.v = cz;
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i) out[g.index(i, j)] = cg.sample(g.x(i), g.y(j));
  (void)st;
  return out;
}

static GridSpec coarse_of(const GridSpec& g, std::size_t st) {
  GridSpec c{g.x0, g.y0, g.dx * st, g.dy * st, (g.nx - 1) / st + 1, (g.ny - 1) / st + 1};
  return c;
}

Grid2D horizon_grid_block(const Model& m, int h, BlockCode block, const GridSpec& g, double margin) {
  const Domain& d = m.info().domain;
  std::vector<double> xs(g.size()), ys(g.size()), z;
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i) { xs[g.index(i, j)] = g.x(i); ys[g.index(i, j)] = g.y(j); }
  const std::size_t st = 4;
  if (optimised_kernels() && g.nx > 4 * st && g.ny > 4 * st) {
    GridSpec c = coarse_of(g, st);
    std::vector<double> cx(c.size()), cy(c.size()), cz;
    for (std::size_t j = 0; j < c.ny; ++j)
      for (std::size_t i = 0; i < c.nx; ++i) { cx[c.index(i, j)] = c.x(i); cy[c.index(i, j)] = c.y(j); }
    vertical_roots(m, h, block, cx, cy, d.z_top + margin, d.z_base - margin, cz);
    auto guess = coarse_guess(g, st, cz, c);
    vertical_roots(m, h, block, xs, ys, d.z_top + margin, d.z_base - margin, z, &guess);
  } else {
    vertical_roots(m, h, block, xs, ys, d.z_top + margin, d.z_base - margin, z);
  }
  Grid2D out(g);
  out.v = std::move(z);
  return out;
}

// ---- faulted horizon -----------------------------------------------------------------------
static void faulted_roots_chunk(const Model& m, int h, const std::vector<double>& xs, const std::vector<double>& ys,
                                std::vector<double>& z, std::vector<BlockCode>& blk, std::size_t* gaps,
                                std::size_t* overlaps, const std::vector<double>* guess);

void faulted_roots(const Model& m, int h, const std::vector<double>& xs, const std::vector<double>& ys,
                   std::vector<double>& z, std::vector<BlockCode>& blk, std::size_t* gaps, std::size_t* overlaps,
                   const std::vector<double>* guess) {
  const std::size_t n = xs.size(), C = 16384;
  if (n <= C) return faulted_roots_chunk(m, h, xs, ys, z, blk, gaps, overlaps, guess);
  z.assign(n, kNaN);
  blk.assign(n, 0);
  std::size_t g = 0, o = 0;
  for (std::size_t a = 0; a < n; a += C) {
    std::size_t b = std::min(n, a + C), gg = 0, oo = 0;
    std::vector<double> cx(xs.begin() + a, xs.begin() + b), cy(ys.begin() + a, ys.begin() + b), cz, cg;
    if (guess) cg.assign(guess->begin() + a, guess->begin() + b);
    std::vector<BlockCode> cb;
    faulted_roots_chunk(m, h, cx, cy, cz, cb, &gg, &oo, guess ? &cg : nullptr);
    std::copy(cz.begin(), cz.end(), z.begin() + a);
    std::copy(cb.begin(), cb.end(), blk.begin() + a);
    g += gg;
    o += oo;
  }
  if (gaps) *gaps = g;
  if (overlaps) *overlaps = o;
}

static void faulted_roots_chunk(const Model& m, int h, const std::vector<double>& xs, const std::vector<double>& ys,
                                std::vector<double>& z, std::vector<BlockCode>& blk, std::size_t* gaps,
                                std::size_t* overlaps, const std::vector<double>* guess) {
  const Domain& d = m.info().domain;
  const std::size_t n = xs.size();
  z.assign(n, kNaN);
  blk.assign(n, 0);
  std::size_t ngap = 0, nover = 0;
  const double zt = d.z_top, zb = d.z_base;
  if (m.info().faults.empty()) {
    vertical_roots(m, h, 0, xs, ys, zt + 2000, zb - 2000, z);
    for (std::size_t i = 0; i < n; ++i)
      if (!std::isnan(z[i]) && (z[i] > zt || z[i] < zb)) z[i] = kNaN;
    if (gaps) *gaps = 0;
    if (overlaps) *overlaps = 0;
    return;
  }
  auto blocks_at = [&](const std::vector<std::size_t>& idx, const std::vector<double>& zz) {
    std::vector<Vec3> p(idx.size());
    for (std::size_t k = 0; k < idx.size(); ++k) p[k] = {xs[idx[k]], ys[idx[k]], zz[k]};
    EvalResult r;
    m.evaluate(p, kEvalBlock, r);
    return r.block;
  };
  if (!optimised_kernels()) {
    // Reference: candidate blocks from 96 samples per column; root in each; keep valid ones.
    const int L = 96;
    std::vector<std::set<BlockCode>> cand(n);
    std::vector<Vec3> p;
    p.reserve(n * L);
    for (std::size_t i = 0; i < n; ++i)
      for (int k = 0; k < L; ++k) p.push_back({xs[i], ys[i], zt - (zt - zb) * (k + 0.5) / L});
    EvalResult r;
    m.evaluate(p, kEvalBlock, r);
    for (std::size_t i = 0; i < n; ++i)
      for (int k = 0; k < L; ++k) cand[i].insert(r.block[i * L + static_cast<std::size_t>(k)]);
    std::map<BlockCode, std::vector<std::size_t>> by;
    for (std::size_t i = 0; i < n; ++i)
      for (auto b : cand[i]) by[b].push_back(i);
    std::vector<std::vector<std::pair<double, BlockCode>>> valid(n);
    std::vector<char> inside_any(n, 0);
    for (auto& [b, idx] : by) {
      std::vector<double> bx, byy, bz;
      for (auto i : idx) { bx.push_back(xs[i]); byy.push_back(ys[i]); }
      vertical_roots(m, h, b, bx, byy, zt + 2000, zb - 2000, bz);
      std::vector<std::size_t> ok;
      std::vector<double> okz;
      for (std::size_t k = 0; k < idx.size(); ++k)
        if (!std::isnan(bz[k]) && bz[k] <= zt && bz[k] >= zb) { ok.push_back(idx[k]); okz.push_back(bz[k]); inside_any[idx[k]] = 1; }
      auto bb = blocks_at(ok, okz);
      for (std::size_t k = 0; k < ok.size(); ++k)
        if (bb[k] == b) valid[ok[k]].push_back({okz[k], b});
    }
    for (std::size_t i = 0; i < n; ++i) {
      auto& v = valid[i];
      if (v.empty()) { if (inside_any[i]) ++ngap; continue; }
      std::sort(v.begin(), v.end(), [](auto& a, auto& b) { return a.first > b.first; });
      z[i] = v[0].first;
      blk[i] = v[0].second;
      if (v.size() > 1 && v[0].first - v.back().first > 0.01) ++nover;
    }
  } else {
    // Optimised: iterate block -> root -> block of the root, grouped by block.
    std::vector<double> zm(n, 0.5 * (zt + zb));
    if (guess)
      for (std::size_t i = 0; i < n; ++i)
        if (!std::isnan((*guess)[i])) zm[i] = (*guess)[i];
    std::vector<std::size_t> all(n);
    for (std::size_t i = 0; i < n; ++i) all[i] = i;
    std::vector<BlockCode> cur = blocks_at(all, zm);
    std::vector<char> done(n, 0), inside_any(n, 0);
    std::vector<std::vector<BlockCode>> seen(n);
    for (int it = 0; it < 8; ++it) {
      std::map<BlockCode, std::vector<std::size_t>> by;
      for (std::size_t i = 0; i < n; ++i)
        if (!done[i]) by[cur[i]].push_back(i);
      if (by.empty()) break;
      for (auto& [b, idx] : by) {
        std::vector<double> bx, byy, bz, bg;
        for (auto i : idx) { bx.push_back(xs[i]); byy.push_back(ys[i]); seen[i].push_back(b); bg.push_back(guess ? (*guess)[i] : kNaN); }
        vertical_roots(m, h, b, bx, byy, zt + 2000, zb - 2000, bz, guess ? &bg : nullptr);
        std::vector<std::size_t> ok;
        std::vector<double> okz;
        for (std::size_t k = 0; k < idx.size(); ++k) {
          if (!std::isnan(bz[k]) && bz[k] <= zt && bz[k] >= zb) { ok.push_back(idx[k]); okz.push_back(bz[k]); inside_any[idx[k]] = 1; }
          else done[idx[k]] = 1;  // no root inside the domain in this block: stop (outside or gap)
        }
        auto bb = blocks_at(ok, okz);
        for (std::size_t k = 0; k < ok.size(); ++k) {
          std::size_t i = ok[k];
          if (bb[k] == b) { z[i] = okz[k]; blk[i] = b; done[i] = 1; }
          else if (std::find(seen[i].begin(), seen[i].end(), bb[k]) != seen[i].end()) done[i] = 1;  // cycle: gap
          else cur[i] = bb[k];
        }
      }
    }
    // Reverse faults: a vertical can meet the horizon in two blocks. Collect the root of every block
    // found on this horizon and keep the upper one, as the reference does.
    bool can_overlap = false;
    for (auto& fi : m.info().faults) can_overlap = can_overlap || fi.overlapping;
    if (can_overlap) {
      std::vector<BlockCode> codes;
      for (std::size_t i = 0; i < n; ++i)
        if (!std::isnan(z[i])) codes.push_back(blk[i]);
      std::sort(codes.begin(), codes.end());
      codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
      for (BlockCode b : codes) {
        std::vector<std::size_t> idx;
        for (std::size_t i = 0; i < n; ++i)
          if (std::isnan(z[i]) || blk[i] != b) idx.push_back(i);
        if (idx.empty()) continue;
        std::vector<double> bx, byy, bz;
        for (auto i : idx) { bx.push_back(xs[i]); byy.push_back(ys[i]); }
        vertical_roots(m, h, b, bx, byy, zt + 2000, zb - 2000, bz);
        std::vector<std::size_t> ok;
        std::vector<double> okz;
        for (std::size_t k = 0; k < idx.size(); ++k)
          if (!std::isnan(bz[k]) && bz[k] <= zt && bz[k] >= zb) { ok.push_back(idx[k]); okz.push_back(bz[k]); }
        auto bb = blocks_at(ok, okz);
        for (std::size_t k = 0; k < ok.size(); ++k) {
          if (bb[k] != b) continue;
          std::size_t i = ok[k];
          if (std::isnan(z[i])) { z[i] = okz[k]; blk[i] = b; continue; }
          if (std::fabs(okz[k] - z[i]) > 0.01) ++nover;
          if (okz[k] > z[i]) { z[i] = okz[k]; blk[i] = b; }  // upper kept
        }
      }
    }
    for (std::size_t i = 0; i < n; ++i)
      if (std::isnan(z[i]) && inside_any[i]) ++ngap;
  }
  if (gaps) *gaps = ngap;
  if (overlaps) *overlaps = nover;
}

FaultedSurface horizon_grid_faulted(const Model& m, int h, const GridSpec& g) {
  std::vector<double> xs(g.size()), ys(g.size());
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i) { xs[g.index(i, j)] = g.x(i); ys[g.index(i, j)] = g.y(j); }
  FaultedSurface fs;
  fs.z = Grid2D(g);
  const std::size_t st = 4;
  if (optimised_kernels() && g.nx > 4 * st && g.ny > 4 * st) {
    // Coarse pass first; its (gap-free) values guide the fine pass.
    GridSpec c = coarse_of(g, st);
    std::vector<double> cx(c.size()), cy(c.size()), cz;
    std::vector<BlockCode> cb;
    for (std::size_t j = 0; j < c.ny; ++j)
      for (std::size_t i = 0; i < c.nx; ++i) { cx[c.index(i, j)] = c.x(i); cy[c.index(i, j)] = c.y(j); }
    faulted_roots(m, h, cx, cy, cz, cb, nullptr, nullptr);
    auto guess = coarse_guess(g, st, cz, c);
    faulted_roots(m, h, xs, ys, fs.z.v, fs.block, &fs.gap_nodes, &fs.overlap_nodes, &guess);
  } else {
    faulted_roots(m, h, xs, ys, fs.z.v, fs.block, &fs.gap_nodes, &fs.overlap_nodes);
  }
  return fs;
}

// Stage 2C (STX-TN-07 section 2): implemented in Prompt 23.
Grid2D unit_thickness_grid(const Model&, int, const GridSpec& g) { return Grid2D(g, std::nan("")); }
std::vector<TerminationLine> termination_lines(const Model&, int, const GridSpec&) { return {}; }

}  // namespace stratonix::repr
