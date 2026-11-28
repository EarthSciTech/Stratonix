// Fault surfaces, sticks, cut-off lines, throw tables, well intersections, sections, voxels.
#include <algorithm>
#include <cmath>
#include <map>

#include "repr_internal.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/repr/repr.hpp"

namespace stratonix::repr {

static bool isnan3(const Vec3& p) { return std::isnan(p.x); }
static const Vec3 kNaN3{kNaN, kNaN, kNaN};

// ---- fault roots along verticals -----------------------------------------------------------
void fault_vertical_roots(const Model& m, int f, const std::vector<double>& xs, const std::vector<double>& ys,
                          double z_hi, double z_lo, std::vector<double>& z_out, std::vector<double>* extent) {
  const std::size_t n = xs.size();
  z_out.assign(n, kNaN);
  if (extent) extent->assign(n, kNaN);
  std::vector<Vec3> p(n);
  FaultEval a, b;
  for (std::size_t i = 0; i < n; ++i) p[i] = {xs[i], ys[i], z_hi};
  m.evaluate_fault_geometry(f, p, a);
  for (std::size_t i = 0; i < n; ++i) p[i].z = z_lo;
  m.evaluate_fault_geometry(f, p, b);
  std::vector<double> hi(n, z_hi), lo(n, z_lo), dhi = a.distance, dlo = b.distance;
  std::vector<std::size_t> act;
  for (std::size_t i = 0; i < n; ++i)
    if ((dhi[i] > 0) != (dlo[i] > 0)) act.push_back(i);
  for (int it = 0; it < 60 && !act.empty(); ++it) {
    std::vector<Vec3> q(act.size());
    for (std::size_t k = 0; k < act.size(); ++k) {
      std::size_t i = act[k];
      q[k] = {xs[i], ys[i], 0.5 * (hi[i] + lo[i])};
    }
    FaultEval e;
    m.evaluate_fault_geometry(f, q, e);
    std::vector<std::size_t> next;
    for (std::size_t k = 0; k < act.size(); ++k) {
      std::size_t i = act[k];
      if ((e.distance[k] > 0) == (dhi[i] > 0)) { hi[i] = q[k].z; dhi[i] = e.distance[k]; }
      else { lo[i] = q[k].z; dlo[i] = e.distance[k]; }
      if (hi[i] - lo[i] < 1e-6) {
        z_out[i] = 0.5 * (hi[i] + lo[i]);
        if (extent) (*extent)[i] = e.extent[k];
      } else {
        next.push_back(i);
      }
    }
    act.swap(next);
  }
}

// ---- fault surface -------------------------------------------------------------------------
FaultSurface fault_surface(const Model& m, int f, double spacing, double extend) {
  const auto& inf = m.info();
  const FaultInfo& fi = inf.faults.at(static_cast<std::size_t>(f));
  const Domain& d = inf.domain;
  FaultSurface fs;
  fs.fault = f;
  double hl = fi.half_length * (1 + extend), bspan = fi.b_max - fi.b_min;
  double b0 = fi.b_min - extend * bspan, b1 = fi.b_max + extend * bspan;
  fs.na = std::max<std::size_t>(2, static_cast<std::size_t>(std::ceil(2 * hl / spacing)) + 1);
  fs.nb = std::max<std::size_t>(2, static_cast<std::size_t>(std::ceil((b1 - b0) / spacing)) + 1);
  fs.na = std::min<std::size_t>(fs.na, 4001);
  fs.nb = std::min<std::size_t>(fs.nb, 4001);
  for (std::size_t i = 0; i < fs.na; ++i) fs.a.push_back(-hl + 2 * hl * i / (fs.na - 1));
  for (std::size_t j = 0; j < fs.nb; ++j) fs.b.push_back(b0 + (b1 - b0) * j / (fs.nb - 1));
  const std::size_t N = fs.na * fs.nb;
  std::vector<Vec3> q(N);
  for (std::size_t j = 0; j < fs.nb; ++j)
    for (std::size_t i = 0; i < fs.na; ++i) q[i + fs.na * j] = fi.ref_point + fi.strike_dir * fs.a[i] + fi.dip_dir * fs.b[j];
  // Secant iteration along the approximate normal (exact in one step for planar faults). Only the
  // distance is needed, and a node whose distance is below 1e-8 keeps its position from then on,
  // so each step evaluates only the nodes still moving (I-051; the same values as evaluating all).
  FaultEval e;
  m.evaluate_fault_geometry(f, q, e);
  std::vector<double> t0(N, 0.0), s0 = e.distance, t1(N), s1(N);
  for (std::size_t k = 0; k < N; ++k) t1[k] = -s0[k];
  std::vector<std::size_t> idx;
  std::vector<Vec3> p;
  FaultEval ee;
  auto eval_at = [&](const std::vector<double>& t, std::vector<double>& out, const std::vector<std::size_t>* only) {
    const std::size_t n = only ? only->size() : N;
    p.resize(n);
    for (std::size_t r = 0; r < n; ++r) {
      const std::size_t k = only ? (*only)[r] : r;
      p[r] = q[k] + fi.normal * t[k];
    }
    m.evaluate_fault_geometry(f, p, ee);
    for (std::size_t r = 0; r < n; ++r) out[only ? (*only)[r] : r] = ee.distance[r];
  };
  eval_at(t1, s1, nullptr);
  std::vector<double> t2(N), s2(N);
  for (int it = 0; it < 20; ++it) {
    double worst = 0;
    for (double v : s1) worst = std::max(worst, std::fabs(v));
    if (worst < 1e-8) break;
    idx.clear();
    for (std::size_t k = 0; k < N; ++k) {
      double den = s1[k] - s0[k];
      t2[k] = std::fabs(den) > 1e-14 ? t1[k] - s1[k] * (t1[k] - t0[k]) / den : t1[k];
      if (std::fabs(s1[k]) < 1e-8) t2[k] = t1[k];
      if (t2[k] == t1[k]) s2[k] = s1[k];  // same point, same distance
      else idx.push_back(k);
    }
    eval_at(t2, s2, &idx);
    t0.swap(t1); s0.swap(s1);
    t1 = t2; s1 = s2;
  }
  fs.nodes.resize(N);
  for (std::size_t k = 0; k < N; ++k) fs.nodes[k] = q[k] + fi.normal * t1[k];
  FaultEval fe;
  m.evaluate_fault_geometry(f, fs.nodes, fe);
  fs.extent.resize(N);
  // Combine the tip-line extent with the domain box (<= 1 inside both).
  double cx = d.x0 + d.sx / 2, cy = d.y0 + d.sy / 2, cz = 0.5 * (d.z_top + d.z_base);
  for (std::size_t k = 0; k < N; ++k) {
    const Vec3& p = fs.nodes[k];
    double dom = std::max({std::fabs(p.x - cx) / (d.sx / 2), std::fabs(p.y - cy) / (d.sy / 2),
                           std::fabs(p.z - cz) / (0.5 * (d.z_top - d.z_base))});
    fs.extent[k] = std::max(fe.extent[k], dom);
    if (std::fabs(fe.distance[k]) > 1e-3) fs.extent[k] = std::max(fs.extent[k], 2.0);  // not converged
  }
  // Triangulate the region extent <= 1 by marching squares on the lattice.
  std::map<std::pair<std::size_t, std::size_t>, std::uint32_t> vid;  // corner or edge -> vertex
  auto corner = [&](std::size_t i, std::size_t j) {
    std::size_t k = i + fs.na * j;
    auto key = std::make_pair(k, k);
    auto it = vid.find(key);
    if (it != vid.end()) return it->second;
    std::uint32_t id = static_cast<std::uint32_t>(fs.mesh.v.size());
    fs.mesh.v.push_back(fs.nodes[k]);
    vid[key] = id;
    return id;
  };
  auto edge = [&](std::size_t k1, std::size_t k2) {
    auto key = std::make_pair(std::min(k1, k2), std::max(k1, k2));
    auto it = vid.find(key);
    if (it != vid.end()) return it->second;
    double e1 = fs.extent[k1] - 1, e2 = fs.extent[k2] - 1;
    double tt = std::clamp(e1 / (e1 - e2), 0.0, 1.0);
    std::uint32_t id = static_cast<std::uint32_t>(fs.mesh.v.size());
    fs.mesh.v.push_back(fs.nodes[k1] + (fs.nodes[k2] - fs.nodes[k1]) * tt);
    vid[key] = id;
    return id;
  };
  for (std::size_t j = 0; j + 1 < fs.nb; ++j)
    for (std::size_t i = 0; i + 1 < fs.na; ++i) {
      std::size_t ks[4] = {i + fs.na * j, i + 1 + fs.na * j, i + 1 + fs.na * (j + 1), i + fs.na * (j + 1)};
      std::size_t ci[4] = {i, i + 1, i + 1, i}, cj[4] = {j, j, j + 1, j + 1};
      bool in[4];
      int nin = 0;
      for (int c = 0; c < 4; ++c) { in[c] = fs.extent[ks[c]] <= 1; nin += in[c]; }
      if (nin == 0) continue;
      std::vector<std::uint32_t> poly;
      for (int c = 0; c < 4; ++c) {
        int c2 = (c + 1) % 4;
        if (in[c]) poly.push_back(corner(ci[c], cj[c]));
        if (in[c] != in[c2]) poly.push_back(edge(ks[c], ks[c2]));
      }
      for (std::size_t k = 1; k + 1 < poly.size(); ++k) {
        std::array<std::uint32_t, 3> tri{poly[0], poly[k], poly[k + 1]};
        Vec3 a = fs.mesh.v[tri[0]], b = fs.mesh.v[tri[1]], c = fs.mesh.v[tri[2]];
        if (0.5 * norm(cross(b - a, c - a)) > 1e-6) fs.mesh.t.push_back(tri);
      }
    }
  return fs;
}

std::vector<Polyline> fault_sticks(const FaultSurface& fs, double spacing) {
  std::vector<Polyline> out;
  if (fs.na < 2) return out;
  double a0 = fs.a.front(), a1 = fs.a.back();
  int nst = std::max(1, static_cast<int>(std::floor((a1 - a0) / spacing)));
  for (int s = 0; s <= nst; ++s) {
    double a = a0 + (a1 - a0) * (nst == 0 ? 0.5 : static_cast<double>(s) / nst);
    // nearest lattice column
    std::size_t i = static_cast<std::size_t>(std::llround((a - a0) / (a1 - a0) * (fs.na - 1)));
    i = std::min(i, fs.na - 1);
    Polyline pl;
    for (std::size_t j = 0; j < fs.nb; ++j) {
      std::size_t k = i + fs.na * j;
      bool in = fs.extent[k] <= 1;
      if (j > 0) {
        std::size_t kp = i + fs.na * (j - 1);
        bool inp = fs.extent[kp] <= 1;
        if (in != inp) {
          double e1 = fs.extent[kp] - 1, e2 = fs.extent[k] - 1;
          double tt = std::clamp(e1 / (e1 - e2), 0.0, 1.0);
          pl.push_back(fs.nodes[kp] + (fs.nodes[k] - fs.nodes[kp]) * tt);
          if (!in) { if (pl.size() >= 2) out.push_back(pl); pl.clear(); }
        }
      }
      if (in) pl.push_back(fs.nodes[k]);
    }
    if (pl.size() >= 2) out.push_back(pl);
  }
  return out;
}

// ---- cut-off lines -------------------------------------------------------------------------
std::vector<CutoffLine> cutoff_lines(const Model& m, const FaultSurface& fs, int h, double eps) {
  const auto& fi = m.info().faults.at(static_cast<std::size_t>(fs.fault));
  const std::size_t N = fs.nodes.size(), nh = m.info().horizons.size();
  std::vector<CutoffLine> res;
  // Only lattice nodes inside the tip line and the domain are read below, so only they are
  // evaluated (I-051); pos maps a lattice node to its place among them.
  constexpr std::size_t npos = static_cast<std::size_t>(-1);
  std::vector<std::size_t> pos(N, npos), inside;
  for (std::size_t k = 0; k < N; ++k)
    if (fs.extent[k] <= 1) { pos[k] = inside.size(); inside.push_back(k); }
  for (int side = 0; side < 2; ++side) {
    double sg = side == 0 ? 1.0 : -1.0;  // 0: hanging wall
    std::vector<Vec3> p(inside.size());
    for (std::size_t r = 0; r < inside.size(); ++r) p[r] = fs.nodes[inside[r]] + fi.normal * (sg * eps);
    EvalResult r;
    m.evaluate(p, kEvalScalars | kEvalBlock, r);
    (void)nh;
    CutoffLine cl;
    cl.fault = fs.fault;
    cl.horizon = h;
    cl.hanging_wall = side == 0;
    cl.by_column.assign(fs.na, kNaN3);
    // Bracket the crossing between lattice nodes, then refine by bisection on the scalar.
    struct Seg { std::size_t col, k1, k2; double lo, hi; bool s1pos; };
    std::vector<Seg> segs;
    for (std::size_t i = 0; i < fs.na; ++i)
      for (std::size_t j = 0; j + 1 < fs.nb; ++j) {
        std::size_t k1 = i + fs.na * j, k2 = i + fs.na * (j + 1);
        if (fs.extent[k1] > 1 || fs.extent[k2] > 1) continue;
        const std::size_t r1 = pos[k1], r2 = pos[k2];
        // A segment whose ends lie in different blocks crosses another fault: the fault surface
        // is offset between them, so no cut-off is interpolated there.
        if (r.block[r1] != r.block[r2]) continue;
        // Both samples must really lie on the intended wall (near offsets of a displaced fault
        // surface the side samples can fall on the other wall).
        bool on_hw = ((r.block[r1] >> fs.fault) & 1u) != 0;
        if (on_hw != (side == 0)) continue;
        double s1 = r.s(static_cast<std::size_t>(h), r1), s2 = r.s(static_cast<std::size_t>(h), r2);
        if ((s1 >= 0) != (s2 >= 0)) { segs.push_back({i, k1, k2, 0.0, 1.0, s1 >= 0}); break; }
      }
    // The side's scalar on the plane itself is extrapolated from eps and 2*eps off the plane,
    // so the cut-off lies where that wall's horizon meets the fault (no offset bias).
    for (int it = 0; it < 40 && !segs.empty(); ++it) {
      const std::size_t ns = segs.size();
      std::vector<Vec3> mid(2 * ns);
      for (std::size_t k = 0; k < ns; ++k) {
        double t = 0.5 * (segs[k].lo + segs[k].hi);
        Vec3 on = fs.nodes[segs[k].k1] + (fs.nodes[segs[k].k2] - fs.nodes[segs[k].k1]) * t;
        mid[k] = on + fi.normal * (sg * eps);
        mid[ns + k] = on + fi.normal * (sg * 2 * eps);
      }
      EvalResult rm;
      m.evaluate(mid, kEvalScalars, rm);
      for (std::size_t k = 0; k < ns; ++k) {
        double t = 0.5 * (segs[k].lo + segs[k].hi);
        double v = 2 * rm.s(static_cast<std::size_t>(h), k) - rm.s(static_cast<std::size_t>(h), ns + k);
        ((v >= 0) == segs[k].s1pos ? segs[k].lo : segs[k].hi) = t;
      }
    }
    // Keep a cut-off only if its side samples lie in the same block as the bracketing nodes
    // (rejects points where the fault surface is stepped by another fault between the nodes).
    std::vector<Vec3> fin(segs.size());
    for (std::size_t k = 0; k < segs.size(); ++k) {
      double t = 0.5 * (segs[k].lo + segs[k].hi);
      fin[k] = fs.nodes[segs[k].k1] + (fs.nodes[segs[k].k2] - fs.nodes[segs[k].k1]) * t;
    }
    std::vector<Vec3> side_pts(fin.size());
    for (std::size_t k = 0; k < fin.size(); ++k) side_pts[k] = fin[k] + fi.normal * (sg * eps);
    EvalResult rf;
    m.evaluate(side_pts, kEvalBlock, rf);
    for (std::size_t k = 0; k < segs.size(); ++k)
      if (rf.block[k] == r.block[pos[segs[k].k1]]) cl.by_column[segs[k].col] = fin[k];
    Polyline cur;
    for (std::size_t i = 0; i < fs.na; ++i) {
      if (isnan3(cl.by_column[i])) { if (cur.size() >= 2) cl.parts.push_back(cur); cur.clear(); continue; }
      cur.push_back(cl.by_column[i]);
    }
    if (cur.size() >= 2) cl.parts.push_back(cur);
    res.push_back(std::move(cl));
  }
  return res;
}

std::vector<double> scalar_jumps(const Model& m, const std::vector<Vec3>& pts, const Vec3& nrm, double eps) {
  const std::size_t n = pts.size(), nh = m.info().horizons.size();
  std::vector<Vec3> q(4 * n);
  for (std::size_t i = 0; i < n; ++i) {
    q[i] = pts[i] + nrm * eps;
    q[n + i] = pts[i] + nrm * (2 * eps);
    q[2 * n + i] = pts[i] - nrm * eps;
    q[3 * n + i] = pts[i] - nrm * (2 * eps);
  }
  EvalResult r;
  m.evaluate(q, kEvalScalars, r);
  std::vector<double> out(nh * n);
  for (std::size_t h = 0; h < nh; ++h)
    for (std::size_t i = 0; i < n; ++i) {
      double hw = 2 * r.s(h, i) - r.s(h, n + i);
      double fw = 2 * r.s(h, 2 * n + i) - r.s(h, 3 * n + i);
      out[h * n + i] = hw - fw;
    }
  return out;
}

ThrowRow throw_row(const Model& m, const FaultSurface& fs, const std::vector<CutoffLine>& lines, int h) {
  ThrowRow t;
  t.fault = fs.fault;
  t.horizon = h;
  const CutoffLine* hw = nullptr;
  const CutoffLine* fw = nullptr;
  for (auto& l : lines) (l.hanging_wall ? hw : fw) = &l;
  if (!hw || !fw) return t;
  const auto& fi = m.info().faults.at(static_cast<std::size_t>(fs.fault));
  std::size_t imid = 0;
  for (std::size_t i = 0; i < fs.na; ++i)
    if (std::fabs(fs.a[i]) < std::fabs(fs.a[imid])) imid = i;
  std::vector<Vec3> probe;
  std::vector<std::size_t> cols;
  for (std::size_t i = 0; i < fs.na; ++i) {
    const Vec3 &a = hw->by_column[i], &b = fw->by_column[i];
    if (isnan3(a) || isnan3(b)) continue;
    double dz = b.z - a.z;
    double heave = std::hypot(b.x - a.x, b.y - a.y);
    ++t.samples;
    if (std::fabs(dz) > std::fabs(t.throw_max)) t.throw_max = dz;
    t.heave_max = std::max(t.heave_max, heave);
    if (i == imid) t.throw_mid = dz;
    probe.push_back(a - fi.normal * 0.25);
    cols.push_back(i);
  }
  if (!probe.empty()) {
    for (auto& p : probe) p = p + fi.normal * 0.25;  // back onto the hanging-wall cut-off
    auto j = scalar_jumps(m, probe, fi.normal);
    for (std::size_t k = 0; k < probe.size(); ++k) {
      double jump = j[static_cast<std::size_t>(h) * probe.size() + k];
      if (std::fabs(jump) > std::fabs(t.jump_max)) t.jump_max = jump;
    }
  }
  t.cut = std::fabs(t.throw_max) > 1e-3 || std::fabs(t.jump_max) > 1e-3;
  return t;
}

// ---- wells ---------------------------------------------------------------------------------
std::vector<WellTop> well_intersections(const Model& m, const WellPath& w, double step) {
  const auto& inf = m.info();
  const std::size_t nh = inf.horizons.size();
  // Dense samples along the path.
  std::vector<Vec3> p;
  const auto& pts = w.points();
  for (std::size_t s = 1; s < pts.size(); ++s) {
    Vec3 a = pts[s - 1], b = pts[s];
    double len = norm(b - a);
    int n = std::max(1, static_cast<int>(std::ceil(len / step)));
    for (int k = (s == 1 ? 0 : 1); k <= n; ++k) p.push_back(a + (b - a) * (static_cast<double>(k) / n));
  }
  EvalResult r;
  m.evaluate(p, kEvalScalars | kEvalBlock, r);
  std::vector<WellTop> tops;
  auto eval1 = [&](const Vec3& q, double& s, BlockCode& b, std::size_t h) {
    EvalResult e;
    m.evaluate(std::span<const Vec3>(&q, 1), kEvalScalars | kEvalBlock, e);
    s = e.scalar[h];
    b = e.block[0];
  };
  for (std::size_t h = 0; h < nh; ++h) {
    for (std::size_t k = 0; k + 1 < p.size(); ++k) {
      double s0 = r.s(h, k), s1 = r.s(h, k + 1);
      if ((s0 >= 0) == (s1 >= 0) && r.block[k] == r.block[k + 1]) continue;
      // Refine: split the segment at block changes; keep sign changes inside one block.
      std::vector<std::pair<Vec3, Vec3>> stack{{p[k], p[k + 1]}};
      int guard = 0;
      while (!stack.empty() && guard++ < 200) {
        auto [a, b] = stack.back();
        stack.pop_back();
        double sa, sb;
        BlockCode ba, bb;
        eval1(a, sa, ba, h);
        eval1(b, sb, bb, h);
        if ((sa >= 0) == (sb >= 0) && ba == bb) continue;
        if (ba == bb) {
          for (int it = 0; it < 60; ++it) {
            Vec3 c = (a + b) * 0.5;
            double sc;
            BlockCode bc;
            eval1(c, sc, bc, h);
            ((sc >= 0) == (sa >= 0) ? a : b) = c;
          }
          Vec3 c = (a + b) * 0.5;
          tops.push_back({w.name(), inf.horizons[h].name, static_cast<int>(h), c, w.md_at_elevation(c.z), ba});
        } else if (norm(b - a) > 1e-4) {
          Vec3 c = (a + b) * 0.5;
          stack.push_back({c, b});
          stack.push_back({a, c});
        }
      }
    }
  }
  std::sort(tops.begin(), tops.end(), [](const WellTop& a, const WellTop& b) { return a.p.z > b.p.z; });
  // Stage 2C: a horizon clamped by a series relation at the well (eroded, or lapping out) has no
  // top there (STX-TN-07 section 2).
  if (auto* st = static_cast<const StratigraphyQuery*>(m.capability("stratigraphy")); st && !tops.empty()) {
    std::vector<Vec3> at, rs(tops.size());
    for (auto& t : tops) at.push_back(t.p);
    m.restore(at, rs);
    std::vector<WellTop> kept;
    for (std::size_t q = 0; q < tops.size(); ++q)
      if (!st->horizon_absent(tops[q].horizon_index, rs[q].x, rs[q].y)) kept.push_back(tops[q]);
    tops.swap(kept);
  }
  return tops;
}

// ---- sections ------------------------------------------------------------------------------
Section make_section(const Model& m, const std::string& name, const Polyline& line, double ds, double dz) {
  const auto& inf = m.info();
  Section s;
  s.name = name;
  s.line = line;
  s.ds = ds;
  s.dz = dz;
  s.z_top = inf.domain.z_top;
  s.z_base = inf.domain.z_base;
  std::vector<double> cum{0};
  for (std::size_t i = 1; i < line.size(); ++i)
    cum.push_back(cum.back() + std::hypot(line[i].x - line[i - 1].x, line[i].y - line[i - 1].y));
  s.length = cum.back();
  s.ns = static_cast<std::size_t>(std::ceil(s.length / ds)) + 1;
  s.nz = static_cast<std::size_t>(std::ceil((s.z_top - s.z_base) / dz)) + 1;
  for (std::size_t c = 0; c < s.ns; ++c) {
    double d = std::min(s.length, c * ds);
    std::size_t seg = 1;
    while (seg < line.size() - 1 && cum[seg] < d) ++seg;
    double t = cum[seg] > cum[seg - 1] ? (d - cum[seg - 1]) / (cum[seg] - cum[seg - 1]) : 0;
    s.station.push_back(d);
    s.xy.push_back({line[seg - 1].x + t * (line[seg].x - line[seg - 1].x), line[seg - 1].y + t * (line[seg].y - line[seg - 1].y), 0});
  }
  std::vector<Vec3> pts(s.ns * s.nz);
  for (std::size_t r = 0; r < s.nz; ++r)
    for (std::size_t c = 0; c < s.ns; ++c) pts[r * s.ns + c] = {s.xy[c].x, s.xy[c].y, std::max(s.z_base, s.z_top - r * dz)};
  EvalResult er;
  m.evaluate(pts, kEvalUnit | kEvalBlock, er);
  s.unit = std::move(er.unit);
  s.block = std::move(er.block);
  std::vector<double> xs, ys;
  for (auto& p : s.xy) { xs.push_back(p.x); ys.push_back(p.y); }
  for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
    std::vector<double> z;
    std::vector<BlockCode> b;
    faulted_roots(m, static_cast<int>(h), xs, ys, z, b, nullptr, nullptr);
    s.horizon_z.push_back(z);
  }
  bool can_overlap = false;
  for (auto& fi : inf.faults) can_overlap = can_overlap || fi.overlapping;
  if (can_overlap) {
    // Trace each block's horizon on its own, keeping a root only where the model puts it in that
    // block; used for drawing where two blocks hold the same horizon in one column.
    s.horizon_block_z.resize(inf.horizons.size());
    std::vector<BlockCode> codes(s.block.begin(), s.block.end());
    std::sort(codes.begin(), codes.end());
    codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
    bool any_overlap = false;
    for (std::size_t h = 0; h < inf.horizons.size(); ++h)
      for (BlockCode bc : codes) {
        std::vector<double> z;
        vertical_roots(m, static_cast<int>(h), bc, xs, ys, s.z_top + 2000, s.z_base - 2000, z);
        std::vector<Vec3> p(z.size());
        for (std::size_t i = 0; i < z.size(); ++i) p[i] = {xs[i], ys[i], std::isnan(z[i]) ? 0.0 : z[i]};
        EvalResult r;
        m.evaluate(p, kEvalBlock, r);
        bool any = false;
        for (std::size_t i = 0; i < z.size(); ++i) {
          if (std::isnan(z[i]) || z[i] > s.z_top || z[i] < s.z_base || r.block[i] != bc) z[i] = kNaN;
          else any = true;
        }
        if (any) s.horizon_block_z[h].push_back({bc, std::move(z)});
      }
    for (std::size_t h = 0; h < inf.horizons.size() && !any_overlap; ++h)
      for (std::size_t i = 0; i < xs.size() && !any_overlap; ++i) {
        int n = 0;
        for (auto& [bc, zz] : s.horizon_block_z[h]) n += std::isnan(zz[i]) ? 0 : 1;
        any_overlap = n > 1;
      }
    if (!any_overlap) s.horizon_block_z.clear();
  }
  for (std::size_t f = 0; f < inf.faults.size(); ++f) {
    std::vector<double> z, e;
    fault_vertical_roots(m, static_cast<int>(f), xs, ys, s.z_top, s.z_base, z, &e);
    for (std::size_t i = 0; i < z.size(); ++i)
      if (std::isnan(e[i]) || e[i] > 1) z[i] = kNaN;
    s.fault_z.push_back(z);
  }
  return s;
}

// ---- voxels --------------------------------------------------------------------------------
VoxelGrid voxel_grid(const Model& m, double dx, double dy, double dz, bool with_block) {
  const Domain& d = m.info().domain;
  VoxelGrid v;
  v.dx = dx; v.dy = dy; v.dz = dz;
  v.nx = static_cast<std::size_t>(std::floor(d.sx / dx));
  v.ny = static_cast<std::size_t>(std::floor(d.sy / dy));
  v.nz = static_cast<std::size_t>(std::floor((d.z_top - d.z_base) / dz));
  v.x0 = d.x0 + dx / 2; v.y0 = d.y0 + dy / 2; v.z0 = d.z_base + dz / 2;
  const std::size_t slab = v.nx * v.ny;
  v.unit.resize(slab * v.nz);
  if (with_block) v.block.resize(slab * v.nz);
  std::vector<Vec3> p(slab);
  for (std::size_t k = 0; k < v.nz; ++k) {
    for (std::size_t j = 0; j < v.ny; ++j)
      for (std::size_t i = 0; i < v.nx; ++i) p[i + v.nx * j] = v.centre(i, j, k);
    EvalResult r;
    m.evaluate(p, kEvalUnit | (with_block ? kEvalBlock : 0u), r);
    std::copy(r.unit.begin(), r.unit.end(), v.unit.begin() + k * slab);
    if (with_block) std::copy(r.block.begin(), r.block.end(), v.block.begin() + k * slab);
  }
  return v;
}

}  // namespace stratonix::repr
