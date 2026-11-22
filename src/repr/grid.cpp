// Grids and properties (stage 3B, STX-TN-12): corner-point grids conforming to faults and
// horizons (Prompt 38); properties, fault transmissibility, the upscaling interface and
// tetrahedral grids (Prompt 39).
#include "stratonix/repr/grid.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <stdexcept>
#include <tuple>

#include "repr_internal.hpp"
#include "stratonix/core/util.hpp"

namespace stratonix::repr {

const char* layering_name(Layering l) {
  switch (l) {
    case Layering::Proportional: return "proportional";
    case Layering::FollowTop: return "follow_top";
    case Layering::FollowBase: return "follow_base";
  }
  return "?";
}

// ---- grid geometry ---------------------------------------------------------------------------

std::size_t CornerPointGrid::zcorn_index(int i, int j, int k, int c) const {
  const std::size_t ib = c & 1, jb = (c >> 1) & 1, bot = (c >> 2) & 1;
  const std::size_t NX = 2 * static_cast<std::size_t>(nx), NY = 2 * static_cast<std::size_t>(ny);
  return ((2 * static_cast<std::size_t>(k) + bot) * NY + (2 * static_cast<std::size_t>(j) + jb)) * NX + 2 * static_cast<std::size_t>(i) + ib;
}

namespace {
/// Point on pillar p at elevation z (linear along the pillar, extended beyond its ends).
Vec3 on_pillar(const std::vector<double>& coord, std::size_t p, double z) {
  const double* c = &coord[6 * p];
  double dz = c[5] - c[2];
  double t = std::fabs(dz) > 1e-12 ? (z - c[2]) / dz : 0.0;
  return {c[0] + t * (c[3] - c[0]), c[1] + t * (c[4] - c[1]), z};
}
Vec3 pillar_at_t(const std::vector<double>& coord, std::size_t p, double t) {
  const double* c = &coord[6 * p];
  return {c[0] + t * (c[3] - c[0]), c[1] + t * (c[4] - c[1]), c[2] + t * (c[5] - c[2])};
}
}  // namespace

Vec3 CornerPointGrid::corner(int i, int j, int k, int c) const {
  const int ib = c & 1, jb = (c >> 1) & 1;
  return on_pillar(coord, pillar(i + ib, j + jb), zcorn[zcorn_index(i, j, k, c)]);
}
std::array<Vec3, 8> CornerPointGrid::corners(int i, int j, int k) const {
  std::array<Vec3, 8> a;
  for (int c = 0; c < 8; ++c) a[c] = corner(i, j, k, c);
  return a;
}

double cell_volume(const std::array<Vec3, 8>& X) {
  // V = -integral of det J over the unit cube (k, and so the third coordinate, points down).
  const double g = 0.5 / std::sqrt(3.0), q[2] = {0.5 - g, 0.5 + g};
  double v = 0;
  for (double a : q)
    for (double b : q)
      for (double cz : q) {
        Vec3 du{0, 0, 0}, dv{0, 0, 0}, dw{0, 0, 0};
        for (int c = 0; c < 8; ++c) {
          const int ib = c & 1, jb = (c >> 1) & 1, kb = (c >> 2) & 1;
          double Na = ib ? a : 1 - a, Nb = jb ? b : 1 - b, Nc = kb ? cz : 1 - cz;
          double sa = ib ? 1 : -1, sb = jb ? 1 : -1, sc = kb ? 1 : -1;
          du = du + X[c] * (sa * Nb * Nc);
          dv = dv + X[c] * (sb * Na * Nc);
          dw = dw + X[c] * (sc * Na * Nb);
        }
        v -= dot(cross(du, dv), dw) * 0.125;
      }
  return v;
}
Vec3 cell_centre(const std::array<Vec3, 8>& c) {
  Vec3 s{0, 0, 0};
  for (auto& p : c) s = s + p;
  return s / 8.0;
}
double tet_volume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  return dot(b - a, cross(c - a, d - a)) / 6.0;
}
std::array<std::array<int, 4>, 6> hex_tets_fixed() {
  // Kuhn (Freudenthal) split along the 0-7 diagonal: 0, 0 + e_a, 0 + e_a + e_b, 7 for the six
  // orders of the axes, each oriented positive for a hexahedron whose k points down.
  static const std::array<std::array<int, 4>, 6> t = [] {
    std::array<std::array<int, 4>, 6> r{};
    const int perm[6][3] = {{1, 2, 4}, {1, 4, 2}, {2, 1, 4}, {2, 4, 1}, {4, 1, 2}, {4, 2, 1}};
    auto ref = [](int c) { return Vec3{static_cast<double>(c & 1), static_cast<double>((c >> 1) & 1), -static_cast<double>((c >> 2) & 1)}; };
    for (int n = 0; n < 6; ++n) {
      r[n] = {0, perm[n][0], perm[n][0] + perm[n][1], 7};
      if (tet_volume(ref(r[n][0]), ref(r[n][1]), ref(r[n][2]), ref(r[n][3])) < 0) std::swap(r[n][2], r[n][3]);
    }
    return r;
  }();
  return t;
}

namespace {

/// Fault crossings along lattice lines at one elevation: for `along_x`, lines y = y_j sampled in x.
struct Crossing { int line; double s; };  // line index, coordinate along the line

void fault_crossings(const Model& m, int f, const Domain& d, bool along_x, int nlines, double d_line, double z, double step,
                     std::vector<Crossing>& out) {
  out.clear();
  const double s0 = along_x ? d.x0 : d.y0, len = along_x ? d.sx : d.sy, l0 = along_x ? d.y0 : d.x0;
  const int ns = std::max(2, static_cast<int>(std::ceil(len / step)) + 1);
  const double ds = len / (ns - 1);
  std::vector<Vec3> p;
  p.reserve(static_cast<std::size_t>(nlines) * ns);
  for (int l = 0; l < nlines; ++l)
    for (int s = 0; s < ns; ++s) {
      double a = s0 + s * ds, b = l0 + l * d_line;
      p.push_back(along_x ? Vec3{a, b, z} : Vec3{b, a, z});
    }
  FaultEval e;
  m.evaluate_fault(f, p, e);
  std::vector<std::array<double, 4>> br;  // line, lo, hi, d(lo)
  for (int l = 0; l < nlines; ++l)
    for (int s = 0; s + 1 < ns; ++s) {
      double a = e.distance[static_cast<std::size_t>(l) * ns + s], b = e.distance[static_cast<std::size_t>(l) * ns + s + 1];
      if (std::isnan(a) || std::isnan(b)) continue;
      if ((a > 0) != (b > 0)) br.push_back({static_cast<double>(l), s0 + s * ds, s0 + (s + 1) * ds, a});
    }
  for (int it = 0; it < 60 && !br.empty(); ++it) {
    std::vector<Vec3> q(br.size());
    for (std::size_t k = 0; k < br.size(); ++k) {
      double a = 0.5 * (br[k][1] + br[k][2]), b = l0 + br[k][0] * d_line;
      q[k] = along_x ? Vec3{a, b, z} : Vec3{b, a, z};
    }
    FaultEval ee;
    m.evaluate_fault(f, q, ee);
    for (std::size_t k = 0; k < br.size(); ++k) {
      double mid = 0.5 * (br[k][1] + br[k][2]);
      if ((ee.distance[k] > 0) == (br[k][3] > 0)) { br[k][1] = mid; br[k][3] = ee.distance[k]; } else br[k][2] = mid;
    }
  }
  for (auto& b : br) out.push_back({static_cast<int>(b[0]), 0.5 * (b[1] + b[2])});
}

double nearest_crossing(const std::vector<Crossing>& cs, int line, double s) {
  double best = kNaN, bd = 1e300;
  for (auto& c : cs)
    if (c.line == line && std::fabs(c.s - s) < bd) { bd = std::fabs(c.s - s); best = c.s; }
  return best;
}

/// Elevation where horizon h's scalar in `block` crosses zero going down pillar segment
/// t in [-0.5, 1.5]; NaN when it does not.
struct RootTask { std::size_t pillar; BlockCode block; int h; double z = kNaN; };

void pillar_roots(const Model& m, const std::vector<double>& coord, std::vector<RootTask>& tasks) {
  // group by (block, horizon) so each scalar_in_block call shares them
  std::map<std::pair<BlockCode, int>, std::vector<std::size_t>> groups;
  for (std::size_t n = 0; n < tasks.size(); ++n) groups[{tasks[n].block, tasks[n].h}].push_back(n);
  std::vector<std::pair<std::pair<BlockCode, int>, std::vector<std::size_t>>> glist(groups.begin(), groups.end());
  std::vector<std::pair<std::size_t, std::size_t>> chunks;  // (group, offset)
  const std::size_t C = 512;
  for (std::size_t g = 0; g < glist.size(); ++g)
    for (std::size_t o = 0; o < glist[g].second.size(); o += C) chunks.push_back({g, o});
  parallel_for(chunks.size(), 1, [&](std::size_t a, std::size_t b) {
    for (std::size_t ci = a; ci < b; ++ci) {
      auto& [key, ids_all] = glist[chunks[ci].first];
      std::size_t o = chunks[ci].second, n = std::min(C, ids_all.size() - o);
      const int NS = 33;
      std::vector<Vec3> p(n * NS);
      for (std::size_t k = 0; k < n; ++k)
        for (int s = 0; s < NS; ++s) p[k * NS + s] = pillar_at_t(coord, tasks[ids_all[o + k]].pillar, -0.5 + 2.0 * s / (NS - 1));
      std::vector<double> v(p.size());
      m.scalar_in_block(p, key.first, key.second, v);
      std::vector<double> lo(n, kNaN), hi(n, kNaN);
      for (std::size_t k = 0; k < n; ++k)
        for (int s = 0; s + 1 < NS; ++s) {
          double va = v[k * NS + s], vb = v[k * NS + s + 1];
          if (va >= 0 && vb < 0) { lo[k] = -0.5 + 2.0 * s / (NS - 1); hi[k] = -0.5 + 2.0 * (s + 1) / (NS - 1); break; }
        }
      std::vector<std::size_t> act;
      for (std::size_t k = 0; k < n; ++k)
        if (!std::isnan(lo[k])) act.push_back(k);
      std::vector<Vec3> q(act.size());
      std::vector<double> qv(act.size());
      for (int it = 0; it < 48 && !act.empty(); ++it) {
        for (std::size_t a2 = 0; a2 < act.size(); ++a2) {
          std::size_t k = act[a2];
          q[a2] = pillar_at_t(coord, tasks[ids_all[o + k]].pillar, 0.5 * (lo[k] + hi[k]));
        }
        m.scalar_in_block(q, key.first, key.second, qv);
        for (std::size_t a2 = 0; a2 < act.size(); ++a2) {
          std::size_t k = act[a2];
          double t = 0.5 * (lo[k] + hi[k]);
          (qv[a2] >= 0 ? lo[k] : hi[k]) = t;
        }
      }
      for (std::size_t k = 0; k < n; ++k)
        if (!std::isnan(lo[k])) tasks[ids_all[o + k]].z = pillar_at_t(coord, tasks[ids_all[o + k]].pillar, 0.5 * (lo[k] + hi[k])).z;
    }
  });
}

}  // namespace

CornerPointGrid corner_point_grid(const Model& m, const GridOptions& o) {
  const ModelInfo& info = m.info();
  const Domain& d = info.domain;
  if (!(o.dx > 0) || !(o.dy > 0)) throw std::invalid_argument("grid: spacing must be positive");
  if (o.zones.empty()) throw std::invalid_argument("grid: no zones");
  const int nunits = static_cast<int>(info.units.size());
  for (std::size_t z = 0; z < o.zones.size(); ++z) {
    const auto& zs = o.zones[z];
    if (zs.unit < 0 || zs.unit >= nunits || info.units[zs.unit].top_horizon < 0 || info.units[zs.unit].base_horizon < 0)
      throw std::invalid_argument("grid: zone " + std::to_string(z) + " is not a unit between two horizons");
    if (zs.layers < 1) throw std::invalid_argument("grid: zone " + info.units[zs.unit].name + " needs at least one layer");
    if (zs.layering != Layering::Proportional && !(zs.layer_thickness > 0))
      throw std::invalid_argument("grid: zone " + info.units[zs.unit].name + " needs a positive layer thickness");
    if (z > 0 && info.units[zs.unit].top_horizon != info.units[o.zones[z - 1].unit].base_horizon)
      throw std::invalid_argument("grid: zones must be consecutive units, top to base");
  }
  CornerPointGrid g;
  g.nx = std::max(1, static_cast<int>(std::llround(d.sx / o.dx)));
  g.ny = std::max(1, static_cast<int>(std::llround(d.sy / o.dy)));
  g.dx = d.sx / g.nx;
  g.dy = d.sy / g.ny;
  g.x0 = d.x0;
  g.y0 = d.y0;
  for (auto& zs : o.zones) g.nz += zs.layers;
  const int NX = g.nx + 1, NY = g.ny + 1;
  const std::size_t np = static_cast<std::size_t>(NX) * NY;

  // vertical extent: the top zone's top horizon to the base zone's base horizon, plus margin
  const int h_top = info.units[o.zones.front().unit].top_horizon, h_base = info.units[o.zones.back().unit].base_horizon;
  double z_hi = -1e300, z_lo = 1e300;
  {
    std::vector<double> xs, ys, zz;
    std::vector<BlockCode> bb;
    for (int j = 0; j < NY; ++j)
      for (int i = 0; i < NX; ++i) { xs.push_back(g.x0 + i * g.dx); ys.push_back(g.y0 + j * g.dy); }
    faulted_roots(m, h_top, xs, ys, zz, bb, nullptr, nullptr);
    for (double z : zz)
      if (!std::isnan(z)) z_hi = std::max(z_hi, z);
    faulted_roots(m, h_base, xs, ys, zz, bb, nullptr, nullptr);
    for (double z : zz)
      if (!std::isnan(z)) z_lo = std::min(z_lo, z);
  }
  if (z_hi < -1e299) z_hi = d.z_top;
  if (z_lo > 1e299) z_lo = d.z_base;
  z_hi = std::min(d.z_top, z_hi + o.margin);
  z_lo = std::max(d.z_base, z_lo - o.margin);
  if (!(z_hi > z_lo)) throw std::runtime_error("grid: the zones have no vertical extent in the domain");
  const double z_mid = 0.5 * (z_hi + z_lo);

  // pillars: vertical at the lattice nodes, then fault pillars and the tilt of their rows/columns
  std::vector<double> px(np), py(np), off_top_x(np, 0), off_base_x(np, 0), off_top_y(np, 0), off_base_y(np, 0);
  std::vector<double> snap_dist(np, 1e300);
  std::vector<std::array<double, 6>> snapped(np);
  g.pillar_fault.assign(np, -1);
  for (int j = 0; j < NY; ++j)
    for (int i = 0; i < NX; ++i) { px[g.pillar(i, j)] = g.x0 + i * g.dx; py[g.pillar(i, j)] = g.y0 + j * g.dy; }
  if (o.snap_to_faults) {
    for (std::size_t f = 0; f < info.faults.size(); ++f) {
      const auto& fi = info.faults[f];
      const bool along_x = std::fabs(fi.normal.x) >= std::fabs(fi.normal.y);  // rows y = const
      const int nl = along_x ? NY : NX, nn = along_x ? NX : NY;
      const double dl = along_x ? g.dy : g.dx, dn = along_x ? g.dx : g.dy, s0 = along_x ? g.x0 : g.y0;
      std::vector<Crossing> cm, ct, cb;
      fault_crossings(m, static_cast<int>(f), d, along_x, nl, dl, z_mid, dn / 4, cm);
      fault_crossings(m, static_cast<int>(f), d, along_x, nl, dl, z_hi, dn / 4, ct);
      fault_crossings(m, static_cast<int>(f), d, along_x, nl, dl, z_lo, dn / 4, cb);
      // inside the tip line at mid-depth
      std::vector<Vec3> mp(cm.size());
      for (std::size_t k = 0; k < cm.size(); ++k) {
        double b = (along_x ? g.y0 : g.x0) + cm[k].line * dl;
        mp[k] = along_x ? Vec3{cm[k].s, b, z_mid} : Vec3{b, cm[k].s, z_mid};
      }
      FaultEval fe;
      if (!mp.empty()) m.evaluate_fault(static_cast<int>(f), mp, fe);
      for (std::size_t k = 0; k < cm.size(); ++k) {
        if (std::isnan(fe.extent[k]) || fe.extent[k] > 1.0) continue;
        int node = static_cast<int>(std::llround((cm[k].s - s0) / dn));
        if (node <= 0 || node >= nn - 1) continue;  // the domain's edge pillars stay on the edge
        double dist = std::fabs(cm[k].s - (s0 + node * dn));
        if (dist >= 0.5 * dn) continue;
        double st = nearest_crossing(ct, cm[k].line, cm[k].s), sb = nearest_crossing(cb, cm[k].line, cm[k].s);
        if (std::isnan(st) || std::isnan(sb)) continue;
        std::size_t p = along_x ? g.pillar(node, cm[k].line) : g.pillar(cm[k].line, node);
        if (dist >= snap_dist[p]) continue;
        snap_dist[p] = dist;
        g.pillar_fault[p] = static_cast<int>(f);
        double b = (along_x ? g.y0 : g.x0) + cm[k].line * dl;
        snapped[p] = along_x ? std::array<double, 6>{st, b, z_hi, sb, b, z_lo} : std::array<double, 6>{b, st, z_hi, b, sb, z_lo};
      }
    }
    // tilt: along each row (x faults) or column (y faults), the snapped pillars' top and base
    // offsets interpolated linearly to the next snapped pillar or to zero at the domain's edge
    for (int pass = 0; pass < 2; ++pass) {
      const bool along_x = pass == 0;
      const int nl = along_x ? NY : NX, nn = along_x ? NX : NY;
      for (int l = 0; l < nl; ++l) {
        std::vector<std::pair<int, std::pair<double, double>>> known = {{0, {0.0, 0.0}}};
        for (int n = 1; n < nn - 1; ++n) {
          std::size_t p = along_x ? g.pillar(n, l) : g.pillar(l, n);
          int f = g.pillar_fault[p];
          if (f < 0) continue;
          const auto& fi = info.faults[static_cast<std::size_t>(f)];
          if ((std::fabs(fi.normal.x) >= std::fabs(fi.normal.y)) != along_x) continue;
          double node = along_x ? px[p] : py[p];
          double ot = (along_x ? snapped[p][0] : snapped[p][1]) - node, ob = (along_x ? snapped[p][3] : snapped[p][4]) - node;
          known.push_back({n, {ot, ob}});
        }
        known.push_back({nn - 1, {0.0, 0.0}});
        if (known.size() == 2) continue;
        for (std::size_t a = 0; a + 1 < known.size(); ++a) {
          int n0 = known[a].first, n1 = known[a + 1].first;
          for (int n = n0 + 1; n < n1; ++n) {
            double w = static_cast<double>(n - n0) / (n1 - n0);
            std::size_t p = along_x ? g.pillar(n, l) : g.pillar(l, n);
            double ot = (1 - w) * known[a].second.first + w * known[a + 1].second.first;
            double ob = (1 - w) * known[a].second.second + w * known[a + 1].second.second;
            (along_x ? off_top_x : off_top_y)[p] = ot;
            (along_x ? off_base_x : off_base_y)[p] = ob;
          }
        }
      }
    }
  }
  g.coord.resize(6 * np);
  for (std::size_t p = 0; p < np; ++p) {
    if (g.pillar_fault[p] >= 0) {
      std::copy(snapped[p].begin(), snapped[p].end(), g.coord.begin() + static_cast<std::ptrdiff_t>(6 * p));
    } else {
      double c[6] = {px[p] + off_top_x[p], py[p] + off_top_y[p], z_hi, px[p] + off_base_x[p], py[p] + off_base_y[p], z_lo};
      std::copy(c, c + 6, g.coord.begin() + static_cast<std::ptrdiff_t>(6 * p));
    }
  }
  // a row or column whose pillars cross at the top or base keeps vertical pillars (its fault is
  // then followed by cell edges)
  {
    auto crossing_free = [&](bool along_x, int l) {
      const int nn = along_x ? NX : NY;
      for (int n = 0; n + 1 < nn; ++n) {
        std::size_t a = along_x ? g.pillar(n, l) : g.pillar(l, n), b = along_x ? g.pillar(n + 1, l) : g.pillar(l, n + 1);
        int ax = along_x ? 0 : 1;
        double gap = 0.2 * (along_x ? g.dx : g.dy);
        if (g.coord[6 * b + ax] - g.coord[6 * a + ax] < gap || g.coord[6 * b + 3 + ax] - g.coord[6 * a + 3 + ax] < gap) return false;
      }
      return true;
    };
    for (int pass = 0; pass < 2; ++pass) {
      const bool along_x = pass == 0;
      const int nl = along_x ? NY : NX, nn = along_x ? NX : NY;
      for (int l = 0; l < nl; ++l) {
        if (crossing_free(along_x, l)) continue;
        for (int n = 0; n < nn; ++n) {
          std::size_t p = along_x ? g.pillar(n, l) : g.pillar(l, n);
          double c[6] = {px[p], py[p], z_hi, px[p], py[p], z_lo};
          std::copy(c, c + 6, g.coord.begin() + static_cast<std::ptrdiff_t>(6 * p));
          g.pillar_fault[p] = -1;
        }
      }
    }
  }

  // blocks: per column and zone, the most frequent block among the column's samples in the unit
  const std::size_t ncol = static_cast<std::size_t>(g.nx) * g.ny, nzone = o.zones.size();
  std::vector<BlockCode> zblock(ncol * nzone, 0);
  {
    const int NS = 160;
    parallel_for(ncol, 64, [&](std::size_t a, std::size_t b) {
      std::vector<Vec3> p(NS);
      EvalResult r;
      for (std::size_t col = a; col < b; ++col) {
        int i = static_cast<int>(col % g.nx), j = static_cast<int>(col / g.nx);
        for (int s = 0; s < NS; ++s) {
          double t = (s + 0.5) / NS;
          Vec3 c{0, 0, 0};
          for (int q = 0; q < 4; ++q) c = c + pillar_at_t(g.coord, g.pillar(i + (q & 1), j + (q >> 1)), t);
          p[s] = c / 4.0;
        }
        m.evaluate(p, kEvalUnit | kEvalBlock, r);
        BlockCode last = r.block[NS / 2];
        for (std::size_t z = 0; z < nzone; ++z) {
          std::map<BlockCode, int> count;
          for (int s = 0; s < NS; ++s)
            if (r.unit[s] == o.zones[z].unit) count[r.block[s]]++;
          if (!count.empty()) {
            last = std::max_element(count.begin(), count.end(), [](auto& x, auto& y) { return x.second < y.second; })->first;
          }
          zblock[col * nzone + z] = last;
        }
      }
    });
  }

  // zone boundaries on every pillar of every column, in the column's block for that zone
  std::vector<RootTask> tasks;
  std::map<std::tuple<std::size_t, BlockCode, int>, std::size_t> task_of;
  auto task = [&](std::size_t p, BlockCode b, int h) {
    auto key = std::make_tuple(p, b, h);
    auto it = task_of.find(key);
    if (it != task_of.end()) return it->second;
    tasks.push_back({p, b, h});
    task_of[key] = tasks.size() - 1;
    return tasks.size() - 1;
  };
  std::vector<std::array<std::size_t, 8>> col_task(ncol * nzone);  // 4 tops then 4 bases
  for (std::size_t col = 0; col < ncol; ++col) {
    int i = static_cast<int>(col % g.nx), j = static_cast<int>(col / g.nx);
    for (std::size_t z = 0; z < nzone; ++z) {
      const auto& u = info.units[o.zones[z].unit];
      for (int q = 0; q < 4; ++q) {
        std::size_t p = g.pillar(i + (q & 1), j + (q >> 1));
        col_task[col * nzone + z][q] = task(p, zblock[col * nzone + z], u.top_horizon);
        col_task[col * nzone + z][4 + q] = task(p, zblock[col * nzone + z], u.base_horizon);
      }
    }
  }
  pillar_roots(m, g.coord, tasks);

  // cells
  const std::size_t nc = g.ncell();
  g.zcorn.assign(8 * nc, kNaN);
  g.actnum.assign(nc, 0);
  g.zone.assign(nc, -1);
  g.block.assign(nc, 0);
  g.inactive_reason.assign(nc, "");
  g.layer_unit.clear();
  for (auto& zs : o.zones)
    for (int l = 0; l < zs.layers; ++l) g.layer_unit.push_back(zs.unit);
  for (std::size_t col = 0; col < ncol; ++col) {
    int i = static_cast<int>(col % g.nx), j = static_cast<int>(col / g.nx);
    double prev_base[4] = {1e300, 1e300, 1e300, 1e300};
    int k0 = 0;
    for (std::size_t z = 0; z < nzone; ++z) {
      const auto& zs = o.zones[z];
      double zt[4], zb[4];
      bool found = true;
      for (int q = 0; q < 4; ++q) {
        zt[q] = tasks[col_task[col * nzone + z][q]].z;
        zb[q] = tasks[col_task[col * nzone + z][4 + q]].z;
        if (std::isnan(zt[q]) || std::isnan(zb[q])) { found = false; continue; }
        zt[q] = std::min(zt[q], prev_base[q]);  // no overlap with the zone above
        zb[q] = std::min(zb[q], zt[q]);         // no negative thickness
        prev_base[q] = zb[q];
      }
      for (int l = 0; l < zs.layers; ++l) {
        const int k = k0 + l;
        const std::size_t id = g.cell(i, j, k);
        g.zone[id] = zs.unit;
        g.block[id] = zblock[col * nzone + z];
        if (!found) {
          for (int c = 0; c < 8; ++c) g.zcorn[g.zcorn_index(i, j, k, c)] = std::isnan(zt[c & 3]) ? z_hi : zt[c & 3];
          g.inactive_reason[id] = "zone boundary not found on a pillar in the cell's fault block";
          continue;
        }
        double tsum = 0;
        for (int q = 0; q < 4; ++q) {
          const double T = zt[q], B = zb[q], th = T - B, L = zs.layer_thickness;
          auto boundary = [&](int n) {  // n = 0 .. layers, from the top
            switch (zs.layering) {
              case Layering::Proportional: return T - th * n / zs.layers;
              case Layering::FollowTop: return n == zs.layers ? B : std::max(B, T - n * L);
              case Layering::FollowBase: return n == 0 ? T : std::min(T, B + (zs.layers - n) * L);
            }
            return T;
          };
          double top = boundary(l), bot = boundary(l + 1);
          g.zcorn[g.zcorn_index(i, j, k, q)] = top;
          g.zcorn[g.zcorn_index(i, j, k, q + 4)] = bot;
          tsum += top - bot;
        }
        if (tsum / 4 >= o.min_thickness) {
          g.actnum[id] = 1;
        } else {
          g.inactive_reason[id] = "thinner than the minimum thickness (pinch-out, truncation or layering)";
        }
      }
      k0 += zs.layers;
    }
  }
  return g;
}

// ---- checks (C23 to C25) ---------------------------------------------------------------------

GridDefects grid_defects(const Model& m, const CornerPointGrid& g, double fault_tol) {
  GridDefects d;
  auto label = [](int i, int j, int k) {
    return "cell (" + std::to_string(i + 1) + ", " + std::to_string(j + 1) + ", " + std::to_string(k + 1) + ")";
  };
  for (int k = 0; k < g.nz; ++k)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i) {
        const auto id = g.cell(i, j, k);
        if (!g.actnum[id]) continue;
        auto c = g.corners(i, j, k);
        double v = cell_volume(c);
        if (!(v > 0)) {
          if (d.nonpositive_volume++ == 0) d.first_nonpositive = label(i, j, k) + ": volume " + fmt_fixed(v, 3) + " m3";
          d.worst_volume = std::min(d.worst_volume, v);
        }
        // inversion: a top corner below its bottom corner, or a negative corner Jacobian
        bool inv = false;
        for (int q = 0; q < 4; ++q)
          if (g.zcorn[g.zcorn_index(i, j, k, q)] < g.zcorn[g.zcorn_index(i, j, k, q + 4)] - 1e-9) inv = true;
        if (!inv) {
          const double scale = g.dx * g.dy * 1e-9;
          for (int q = 0; q < 8; ++q) {
            int ib = q & 1, jb = (q >> 1) & 1, kb = (q >> 2) & 1;
            Vec3 ei = c[q ^ 1] - c[q], ej = c[q ^ 2] - c[q], ek = c[q ^ 4] - c[q];
            double sgn = (ib ? -1 : 1) * (jb ? -1 : 1) * (kb ? -1 : 1);
            double det = -sgn * dot(cross(ei, ej), ek);
            double h = norm(ek);
            if (h > 1e-6 && det < -scale * std::max(1.0, h)) inv = true;
          }
        }
        if (inv && d.inverted++ == 0) d.first_inverted = label(i, j, k);
      }
  // C25: faces whose two pillars lie on the same fault
  FaultEval fe;
  for (int dir = 0; dir < 2; ++dir)
    for (int j = 0; j < g.ny + (dir == 0 ? 0 : 1); ++j)
      for (int i = 0; i < g.nx + (dir == 0 ? 1 : 0); ++i) {
        // dir 0: X face at pillar line i between pillars (i, j) and (i, j + 1)
        // dir 1: Y face at pillar line j between pillars (i, j) and (i + 1, j)
        std::size_t pa = g.pillar(i, j), pb = dir == 0 ? g.pillar(i, j + 1) : g.pillar(i + 1, j);
        if (dir == 0 && (i == 0 || i == g.nx)) continue;
        if (dir == 1 && (j == 0 || j == g.ny)) continue;
        int f = g.pillar_fault[pa];
        if (f < 0 || g.pillar_fault[pb] != f) continue;
        for (int k = 0; k < g.nz; ++k) {
          std::vector<Vec3> pts;
          // the cells on both sides, their four corners on the face
          int il = dir == 0 ? i - 1 : i, jl = dir == 0 ? j : j - 1;
          for (int side = 0; side < 2; ++side) {
            int ci = side == 0 ? il : i, cj = side == 0 ? jl : j;
            if (!g.actnum[g.cell(ci, cj, k)]) continue;
            for (int q = 0; q < 8; ++q) {
              int ib = q & 1, jb = (q >> 1) & 1;
              bool on = dir == 0 ? (ib == (side == 0 ? 1 : 0)) : (jb == (side == 0 ? 1 : 0));
              if (on) pts.push_back(g.corner(ci, cj, k, q));
            }
          }
          if (pts.empty()) continue;
          ++d.fault_faces;
          m.evaluate_fault(f, pts, fe);
          // allowance for the chord of a curved fault: its distance from the pillar at mid-depth
          double sag = 0;
          {
            std::vector<Vec3> mids = {pillar_at_t(g.coord, pa, 0.5), pillar_at_t(g.coord, pb, 0.5)};
            FaultEval fm;
            m.evaluate_fault(f, mids, fm);
            sag = std::max(std::fabs(fm.distance[0]), std::fabs(fm.distance[1]));
          }
          double worst = 0;
          for (double x : fe.distance) worst = std::max(worst, std::fabs(x));
          d.worst_fault_distance = std::max(d.worst_fault_distance, worst);
          if (worst > fault_tol + sag && d.fault_face_off++ == 0)
            d.first_fault_off = "fault " + m.info().faults[static_cast<std::size_t>(f)].name + " face at pillar (" + std::to_string(i + 1) + ", " +
                                std::to_string(j + 1) + "), layer " + std::to_string(k + 1) + ": " + fmt_fixed(worst, 3) + " m";
        }
      }
  return d;
}

// ---- fault faces ---------------------------------------------------------------------------

std::vector<FaultFace> grid_fault_faces(const Model& m, const CornerPointGrid& g) {
  std::vector<FaultFace> faces;
  const int nf = static_cast<int>(m.info().faults.size());
  for (int k = 0; k < g.nz; ++k)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i)
        for (int dir = 0; dir < 2; ++dir) {
          int i2 = dir == 0 ? i + 1 : i, j2 = dir == 0 ? j : j + 1;
          if (i2 >= g.nx || j2 >= g.ny) continue;
          auto a = g.cell(i, j, k), b = g.cell(i2, j2, k);
          if (!g.actnum[a] && !g.actnum[b]) continue;
          BlockCode x = g.block[a] ^ g.block[b];
          if (!x) continue;
          int f = 0;
          while (f < 64 && !((x >> f) & 1u)) ++f;
          if (f >= nf) continue;
          FaultFace ff;
          ff.i = i; ff.j = j; ff.k = k;
          ff.dir = dir == 0 ? 'X' : 'Y';
          ff.fault = f;
          Vec3 c{0, 0, 0};
          for (int q = 0; q < 8; ++q) {
            int ib = q & 1, jb = (q >> 1) & 1;
            if ((dir == 0 && ib == 1) || (dir == 1 && jb == 1)) c = c + g.corner(i, j, k, q);
          }
          ff.centre = c / 4.0;
          faces.push_back(ff);
        }
  std::vector<std::vector<std::size_t>> by_fault(static_cast<std::size_t>(nf));
  for (std::size_t n = 0; n < faces.size(); ++n) by_fault[static_cast<std::size_t>(faces[n].fault)].push_back(n);
  std::vector<char> keep(faces.size(), 0);
  for (int f = 0; f < nf; ++f) {
    auto& ids = by_fault[static_cast<std::size_t>(f)];
    if (ids.empty()) continue;
    std::vector<Vec3> p(ids.size());
    for (std::size_t q = 0; q < ids.size(); ++q) p[q] = faces[ids[q]].centre;
    FaultEval e;
    m.evaluate_fault(f, p, e);
    for (std::size_t q = 0; q < ids.size(); ++q) {
      keep[ids[q]] = !std::isnan(e.extent[q]) && e.extent[q] <= 1.0;
      faces[ids[q]].displacement = e.displacement[q];
    }
  }
  std::vector<FaultFace> out;
  for (std::size_t n = 0; n < faces.size(); ++n)
    if (keep[n]) out.push_back(faces[n]);
  return out;
}

// ---- properties (Prompt 39) ------------------------------------------------------------------


const CellProperty* PropertySet::find(const std::string& name) const {
  for (auto& p : props)
    if (p.name == name) return &p;
  return nullptr;
}
CellProperty* PropertySet::find(const std::string& name) {
  for (auto& p : props)
    if (p.name == name) return &p;
  return nullptr;
}
PropertySet grid_properties(const Model& m, const CornerPointGrid& g, const PropertyInputs& in) {
  PropertySet ps;
  const auto& info = m.info();
  const std::size_t nc = g.ncell();
  CellProperty zone{"ZONE", "code", "derived", "model unit index of the cell (ModelInfo units)"};
  zone.integer = true;
  zone.v.resize(nc);
  for (std::size_t c = 0; c < nc; ++c) zone.v[c] = g.zone[c];
  ps.props.push_back(zone);
  // series of a unit: that of its top horizon (stage 2C), when the model has series
  std::vector<std::string> names;
  for (auto& h : info.horizons)
    if (!h.series.empty() && std::find(names.begin(), names.end(), h.series) == names.end()) names.push_back(h.series);
  if (!names.empty()) {
    std::string src = "stratigraphic series of the cell's unit:";
    for (std::size_t k = 0; k < names.size(); ++k) src += " " + std::to_string(k + 1) + " " + names[k];
    CellProperty ser{"SERIES", "code", "derived", src};
    ser.integer = true;
    ser.v.assign(nc, kNaN);
    for (std::size_t c = 0; c < nc; ++c) {
      const int u = g.zone[c];
      if (u < 0) continue;
      const int ht = info.units[static_cast<std::size_t>(u)].top_horizon;
      if (ht < 0) continue;
      auto it = std::find(names.begin(), names.end(), info.horizons[static_cast<std::size_t>(ht)].series);
      if (it != names.end()) ser.v[c] = static_cast<double>(it - names.begin() + 1);
    }
    for (double v : ser.v) ser.missing += std::isnan(v);
    if (ser.missing) ser.missing_reason = "unit without a series";
    ps.props.push_back(ser);
  }
  CellProperty blk{"BLOCK", "code", "derived", "fault-block code of the cell: bit f set on the hanging wall of fault f"};
  blk.integer = true;
  blk.v.resize(nc);
  for (std::size_t c = 0; c < nc; ++c) blk.v[c] = static_cast<double>(g.block[c]);
  ps.props.push_back(blk);
  const std::size_t nunit = info.units.size();
  auto unit_name = [&](int u) { return u >= 0 && static_cast<std::size_t>(u) < nunit ? info.units[static_cast<std::size_t>(u)].name : std::string("?"); };
  // facies by zone: a labelled ordinal per unit; units without one are missing
  if (!in.facies.empty()) {
    std::map<int, const FaciesSpec*> by_unit;
    std::string ev, src = "facies by zone:";
    for (auto& f : in.facies) {
      by_unit[f.unit] = &f;
      if (ev.empty()) ev = f.evidence;
      else if (ev.find(f.evidence) == std::string::npos) ev += " / " + f.evidence;
      src += " " + unit_name(f.unit) + " = " + std::to_string(f.code) + (f.name.empty() ? "" : " (" + f.name + ")") + ";";
      if (!f.source.empty() && src.find(f.source) == std::string::npos) src += " source " + f.source + ";";
    }
    CellProperty fac{"FACIES", "code", ev, src};
    fac.integer = true;
    fac.v.assign(nc, kNaN);
    std::vector<std::string> without;
    for (std::size_t c = 0; c < nc; ++c) {
      auto it = by_unit.find(g.zone[c]);
      if (it != by_unit.end()) fac.v[c] = it->second->code;
      else {
        ++fac.missing;
        auto nm = unit_name(g.zone[c]);
        if (std::find(without.begin(), without.end(), nm) == without.end()) without.push_back(nm);
      }
    }
    for (std::size_t k = 0; k < without.size(); ++k) fac.missing_reason += (k ? ", " : "no facies given for unit ") + without[k];
    ps.props.push_back(fac);
  }
  // porosity from the compaction trend (stage 2D), mean over each cell's depth interval
  if (in.porosity_from_compaction) {
    auto* pq = static_cast<const PorosityQuery*>(m.capability("porosity"));
    CellProperty por{"PORO", "fraction", "derived", "mean porosity of the unit's compaction law (stage 2D) over the cell's restored depth interval"};
    por.v.assign(nc, kNaN);
    if (!pq) {
      por.missing = nc;
      por.missing_reason = "the model offers no compaction porosity (no compaction block, or an engine without one)";
    } else {
      std::vector<std::string> without;
      for (int k = 0; k < g.nz; ++k)
        for (int j = 0; j < g.ny; ++j)
          for (int i = 0; i < g.nx; ++i) {
            auto id = g.cell(i, j, k);
            if (!pq->has_law(g.zone[id])) {
              ++por.missing;
              auto nm = unit_name(g.zone[id]);
              if (std::find(without.begin(), without.end(), nm) == without.end()) without.push_back(nm);
              continue;
            }
            auto cc = g.corners(i, j, k);
            Vec3 top = (cc[0] + cc[1] + cc[2] + cc[3]) / 4.0, base = (cc[4] + cc[5] + cc[6] + cc[7]) / 4.0;
            por.v[id] = pq->mean_porosity(g.zone[id], top, base);
          }
      for (std::size_t k = 0; k < without.size(); ++k) por.missing_reason += (k ? ", " : "no compaction law for unit ") + without[k];
    }
    por.evidence = "derived (from the compaction parameters and their evidence classes)";
    ps.props.push_back(por);
  }
  // transforms, in the order given; each needs a property computed before it
  for (auto& t : in.transforms) {
    const CellProperty* from = ps.find(t.from);
    if (!from) throw std::invalid_argument("grid transform " + t.name + ": no property '" + t.from + "' computed before it");
    ps.props.push_back(apply_transform(*from, t));
  }
  return ps;
}
CellProperty apply_transform(const CellProperty& from, const TransformSpec& t) {
  CellProperty out{t.name, t.unit, t.evidence,
                   t.law + " transform of " + from.name + " (a " + fmt_shortest(t.a) + ", b " + fmt_shortest(t.b) + ")" + (t.source.empty() ? "" : ": " + t.source)};
  out.integer = false;
  out.v.assign(from.v.size(), kNaN);
  std::size_t undefined = 0;
  for (std::size_t c = 0; c < from.v.size(); ++c) {
    const double x = from.v[c];
    if (std::isnan(x)) { ++out.missing; continue; }
    double y = kNaN;
    if (t.law == "linear") y = t.a + t.b * x;
    else if (t.law == "exponential") y = t.a * std::exp(t.b * x);
    else if (t.law == "power") y = x > 0 || t.b == std::floor(t.b) ? t.a * std::pow(x, t.b) : kNaN;
    else throw std::invalid_argument("grid transform " + t.name + ": unknown law '" + t.law + "'");
    if (!std::isfinite(y)) { y = kNaN; ++out.missing; ++undefined; }
    out.v[c] = y;
  }
  if (out.missing) {
    out.missing_reason = from.missing ? from.name + " missing (" + from.missing_reason + ")" : std::string();
    if (undefined) out.missing_reason += (out.missing_reason.empty() ? "" : "; ") + std::to_string(undefined) + " cells where the law is undefined";
  }
  return out;
}

double manzocchi_fault_permeability(double sgr, double displacement_m) {
  if (std::isnan(sgr) || !(displacement_m > 0)) return kNaN;
  return std::pow(10.0, -4.0 * sgr - 0.25 * std::log10(displacement_m) * std::pow(1.0 - sgr, 5));
}
double manzocchi_multiplier(double tf, double kf, double ki, double kj, double Li, double Lj) {
  if (!(tf >= 0) || !(kf > 0) || !(ki > 0) || !(kj > 0) || !(Li > 0) || !(Lj > 0)) return kNaN;
  return 1.0 / (1.0 + tf * (2.0 / kf - 1.0 / ki - 1.0 / kj) / (Li / ki + Lj / kj));
}

TransmissibilityResult fault_transmissibility(const Model& m, const CornerPointGrid& g, const PropertySet& props, const FaultRockInput& rock) {
  TransmissibilityResult r;
  const CellProperty *kx = props.find("PERMX"), *ky = props.find("PERMY");
  if (!kx || !ky) throw std::invalid_argument("fault transmissibility needs PERMX and PERMY");
  const std::size_t nc = g.ncell();
  const std::string src = "Manzocchi et al. (1999): fault permeability from SGR and displacement, thickness from the fault-rock ratio";
  r.multx = {"MULTX", "fraction", "literature_derived calibration on the fault-rock inputs", src};
  r.multy = {"MULTY", "fraction", r.multx.evidence, src};
  r.multx.v.assign(nc, 1.0);
  r.multy.v.assign(nc, 1.0);
  r.faces = grid_fault_faces(m, g);
  StressInput none;
  auto length = [&](int i, int j, int k, char dir) {
    auto cc = g.corners(i, j, k);
    if (dir == 'X') return norm((cc[1] + cc[3] + cc[5] + cc[7]) / 4.0 - (cc[0] + cc[2] + cc[4] + cc[6]) / 4.0);
    return norm((cc[2] + cc[3] + cc[6] + cc[7]) / 4.0 - (cc[0] + cc[1] + cc[4] + cc[5]) / 4.0);
  };
  std::map<std::string, std::size_t> reasons;
  for (auto& f : r.faces) {
    const int i2 = f.dir == 'X' ? f.i + 1 : f.i, j2 = f.dir == 'X' ? f.j : f.j + 1;
    const auto a = g.cell(f.i, f.j, f.k), b = g.cell(i2, j2, f.k);
    const CellProperty& kp = f.dir == 'X' ? *kx : *ky;
    // fault rock where the face meets its fault: the centre moved onto the fault along the normal
    // (a stair-step face lies up to a cell away from it)
    Vec3 on = f.centre;
    for (int it = 0; it < 4; ++it) {
      const double e = 0.5;
      std::vector<Vec3> q{on, on + Vec3{e, 0, 0}, on - Vec3{e, 0, 0}, on + Vec3{0, e, 0}, on - Vec3{0, e, 0}, on + Vec3{0, 0, e}, on - Vec3{0, 0, e}};
      FaultEval fe;
      m.evaluate_fault(f.fault, q, fe);
      Vec3 gr{(fe.distance[1] - fe.distance[2]) / (2 * e), (fe.distance[3] - fe.distance[4]) / (2 * e), (fe.distance[5] - fe.distance[6]) / (2 * e)};
      const double gn = norm(gr);
      if (gn < 1e-9 || std::isnan(fe.distance[0])) break;
      on = on - gr * (fe.distance[0] / (gn * gn));
      if (std::fabs(fe.distance[0]) < 1e-3) break;
    }
    auto pt = fault_rock_at(m, f.fault, on, rock, none);
    f.sgr = pt.sgr;
    f.thickness = pt.thickness;
    f.k_fault = manzocchi_fault_permeability(f.sgr, f.displacement);
    f.k_i = kp.v[a];
    f.k_j = kp.v[b];
    f.L_i = length(f.i, f.j, f.k, f.dir);
    f.L_j = length(i2, j2, f.k, f.dir);
    if (!g.actnum[a] || !g.actnum[b]) f.missing = "an inactive cell on one side";
    else if (!(f.displacement >= 1.0)) {
      // below 1 m of displacement there is no gouge (stage 2F) and the fault-zone thickness
      // (ratio x displacement) tends to zero: the formula's limit, T = 1
      f.mult = 1.0;
      f.note = "displacement below 1 m: no fault rock, multiplier 1 (the formula's limit as the thickness goes to zero)";
    } else if (std::isnan(f.sgr)) f.missing = "no shale gouge ratio (" + (pt.missing.empty() ? std::string("clay fractions not given") : pt.missing) + ")";
    else if (std::isnan(f.thickness)) f.missing = "no fault-zone thickness (no thickness / displacement ratio)";
    else if (!(f.displacement > 0)) f.missing = "no displacement at the face";
    else if (std::isnan(f.k_i) || std::isnan(f.k_j)) f.missing = "cell permeability missing";
    if (f.missing.empty() && f.note.empty()) f.mult = manzocchi_multiplier(f.thickness, f.k_fault, f.k_i, f.k_j, f.L_i, f.L_j);
    if (std::isnan(f.mult) && f.missing.empty()) f.missing = "multiplier undefined for these inputs";
    if (!g.actnum[a] || !g.actnum[b]) continue;  // inactive neighbours: no flow across, nothing to report in MULT
    CellProperty& mp = f.dir == 'X' ? r.multx : r.multy;
    mp.v[a] = f.mult;
    if (!f.missing.empty()) reasons[f.missing]++;
  }
  for (auto* mp : {&r.multx, &r.multy}) {
    for (double v : mp->v) mp->missing += std::isnan(v);
  }
  std::string why;
  for (auto& [k, n] : reasons) why += (why.empty() ? "" : "; ") + std::to_string(n) + " fault faces: " + k;
  r.multx.missing_reason = r.multx.missing ? why : "";
  r.multy.missing_reason = r.multy.missing ? why : "";
  return r;
}

namespace {
/// Inverse bilinear map of (x, y) in the quad a (0,0), b (1,0), c (0,1), d (1,1); false if outside.
bool inverse_bilinear(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d, double x, double y, double& u, double& v) {
  u = 0.5;
  v = 0.5;
  for (int it = 0; it < 30; ++it) {
    double px = (1 - u) * (1 - v) * a.x + u * (1 - v) * b.x + (1 - u) * v * c.x + u * v * d.x - x;
    double py = (1 - u) * (1 - v) * a.y + u * (1 - v) * b.y + (1 - u) * v * c.y + u * v * d.y - y;
    double xu = (1 - v) * (b.x - a.x) + v * (d.x - c.x), xv = (1 - u) * (c.x - a.x) + u * (d.x - b.x);
    double yu = (1 - v) * (b.y - a.y) + v * (d.y - c.y), yv = (1 - u) * (c.y - a.y) + u * (d.y - b.y);
    double det = xu * yv - xv * yu;
    if (std::fabs(det) < 1e-12) return false;
    double du = (px * yv - py * xv) / det, dv = (py * xu - px * yu) / det;
    u -= du;
    v -= dv;
    if (std::fabs(du) + std::fabs(dv) < 1e-12) break;
  }
  return u >= -1e-9 && u <= 1 + 1e-9 && v >= -1e-9 && v <= 1 + 1e-9;
}
}  // namespace

CellProperty upscale_volume_weighted_mean(const CornerPointGrid& fine, const CellProperty& p, const CornerPointGrid& coarse) {
  CellProperty out{p.name, p.unit, "upscaled",
                   "volume_weighted_mean of " + p.name + " (" + p.evidence + ") from a " + std::to_string(fine.nx) + " x " + std::to_string(fine.ny) +
                       " x " + std::to_string(fine.nz) + " grid"};
  const std::size_t nc = coarse.ncell();
  std::vector<double> sum(nc, 0.0), vol(nc, 0.0);
  for (int k = 0; k < fine.nz; ++k)
    for (int j = 0; j < fine.ny; ++j)
      for (int i = 0; i < fine.nx; ++i) {
        auto id = fine.cell(i, j, k);
        if (!fine.actnum[id] || std::isnan(p.v[id])) continue;
        auto cc = fine.corners(i, j, k);
        const Vec3 c = cell_centre(cc);
        const double v = cell_volume(cc);
        const int ci0 = static_cast<int>(std::floor((c.x - coarse.x0) / coarse.dx)), cj0 = static_cast<int>(std::floor((c.y - coarse.y0) / coarse.dy));
        bool done = false;
        for (int dj = -1; dj <= 1 && !done; ++dj)
          for (int di = -1; di <= 1 && !done; ++di) {
            const int ci = ci0 + di, cj = cj0 + dj;
            if (ci < 0 || cj < 0 || ci >= coarse.nx || cj >= coarse.ny) continue;
            Vec3 q[4];
            for (int n = 0; n < 4; ++n) q[n] = on_pillar(coarse.coord, coarse.pillar(ci + (n & 1), cj + (n >> 1)), c.z);
            double u, w;
            if (!inverse_bilinear(q[0], q[1], q[2], q[3], c.x, c.y, u, w)) continue;
            for (int ck = 0; ck < coarse.nz; ++ck) {
              auto cid = coarse.cell(ci, cj, ck);
              auto z = [&](int corner) { return coarse.zcorn[coarse.zcorn_index(ci, cj, ck, corner)]; };
              const double zt = (1 - u) * (1 - w) * z(0) + u * (1 - w) * z(1) + (1 - u) * w * z(2) + u * w * z(3);
              const double zb = (1 - u) * (1 - w) * z(4) + u * (1 - w) * z(5) + (1 - u) * w * z(6) + u * w * z(7);
              if (c.z <= zt + 1e-9 && c.z >= zb - 1e-9) {
                sum[cid] += p.v[id] * v;
                vol[cid] += v;
                done = true;
                break;
              }
            }
          }
      }
  out.v.assign(nc, kNaN);
  for (std::size_t c = 0; c < nc; ++c) {
    if (vol[c] > 0) out.v[c] = sum[c] / vol[c];
    else ++out.missing;
  }
  if (out.missing) out.missing_reason = "no fine cell centre with a value in the coarse cell";
  return out;
}

TetGrid tetrahedral_grid(const CornerPointGrid& g) {
  TetGrid t;
  std::map<std::pair<std::size_t, long long>, std::int32_t> node;
  auto id = [&](int i, int j, int k, int c) {
    std::size_t p = g.pillar(i + (c & 1), j + ((c >> 1) & 1));
    auto key = std::make_pair(p, std::llround(g.zcorn[g.zcorn_index(i, j, k, c)] * 1000.0));
    auto it = node.find(key);
    if (it != node.end()) return it->second;
    t.nodes.push_back(g.corner(i, j, k, c));
    return node[key] = static_cast<std::int32_t>(t.nodes.size() - 1);
  };
  const auto split = hex_tets_fixed();
  for (int k = 0; k < g.nz; ++k)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i) {
        auto c = g.cell(i, j, k);
        if (!g.actnum[c]) continue;
        std::int32_t n[8];
        for (int q = 0; q < 8; ++q) n[q] = id(i, j, k, q);
        for (auto& tt : split) {
          t.tets.push_back({n[tt[0]], n[tt[1]], n[tt[2]], n[tt[3]]});
          t.cell.push_back(static_cast<std::int32_t>(c));
        }
      }
  return t;
}

// ---- plugins -----------------------------------------------------------------------------------
namespace {
class CornerPointBuiltin final : public CornerPointExtractor {
  CornerPointGrid extract(const Model& m, const GridOptions& o) const override { return corner_point_grid(m, o); }
};
class TetrahedralBuiltin final : public TetrahedralExtractor {
  TetGrid extract(const CornerPointGrid& g) const override { return tetrahedral_grid(g); }
};
class VolumeWeightedMean final : public UpscalerPlugin {
  CellProperty upscale(const CornerPointGrid& fine, const CellProperty& p, const CornerPointGrid& coarse) const override {
    return upscale_volume_weighted_mean(fine, p, coarse);
  }
};
}  // namespace
}  // namespace stratonix::repr

STRATONIX_REGISTER(ext_corner_point, ::stratonix::registry::Kind::Extractor, "corner_point", "1.0.0",
                   "Corner-point grid conforming to faults and horizons (stage 3B)", [] { return std::make_unique<stratonix::repr::CornerPointBuiltin>(); })
STRATONIX_REGISTER(ext_tetrahedral, ::stratonix::registry::Kind::Extractor, "tetrahedral", "1.0.0",
                   "Tetrahedral grid from a corner-point grid (six tetrahedra per cell)", [] { return std::make_unique<stratonix::repr::TetrahedralBuiltin>(); })
STRATONIX_REGISTER(ups_vwm, ::stratonix::registry::Kind::Upscaler, "volume_weighted_mean", "1.0.0",
                   "Volume-weighted mean of the fine cells whose centres fall in each coarse cell",
                   [] { return std::make_unique<stratonix::repr::VolumeWeightedMean>(); })

void stratonix::repr::register_grid_plugins() {
  STRATONIX_PLUGIN_ANCHOR(ext_corner_point);
  STRATONIX_PLUGIN_ANCHOR(ext_tetrahedral);
  STRATONIX_PLUGIN_ANCHOR(ups_vwm);
}
