// Stage 4A (STX-TN-14 section 6): the block-and-contact graph of a model (after Thiele et al.
// 2016: units in fault blocks as nodes, face-neighbour contacts on a lattice as edges) and the
// watertightness of its horizon surfaces against faults, other horizons and the domain boundary.
#include "stratonix/repr/topology.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <tuple>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/repr/repr.hpp"

namespace stratonix {

namespace {
std::string json_text(const std::string& s) {
  std::string o;
  for (char c : s) {
    if (c == '"' || c == '\\') o += '\\';
    o += c;
  }
  return o;
}
int lowest_bit(BlockCode x) {
  for (int b = 0; b < 64; ++b)
    if (x & (BlockCode{1} << b)) return b;
  return -1;
}
}  // namespace

std::string ContactGraph::json() const {
  std::ostringstream o;
  o << "{\n  \"stratonix_topology\": 1,\n  \"spacing_m\": " << fmt_shortest(spacing) << ",\n  \"nodes\": [";
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const auto& n = nodes[i];
    o << (i ? "," : "") << "\n    {\"id\": " << i << ", \"unit\": " << n.unit;
    if (n.unit >= 0 && static_cast<std::size_t>(n.unit) < unit_names.size()) o << ", \"unit_name\": \"" << json_text(unit_names[static_cast<std::size_t>(n.unit)]) << "\"";
    o << ", \"block\": " << n.block << ", \"cells\": " << n.cells << "}";
  }
  o << "\n  ],\n  \"edges\": [";
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    o << (i ? "," : "") << "\n    {\"a\": " << e.a << ", \"b\": " << e.b << ", \"surface\": \"" << json_text(e.surface) << "\", \"area_m2\": " << fmt_fixed(e.area, 0)
      << "}";
  }
  o << "\n  ]\n}\n";
  return o.str();
}

ContactGraph contact_graph(const Model& m, double spacing) {
  if (!(spacing > 0)) throw Error(ErrorCode::InvalidArgument, "contact graph: the spacing must be positive");
  const auto& inf = m.info();
  const auto& d = inf.domain;
  const std::size_t nx = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(d.sx / spacing)));
  const std::size_t ny = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(d.sy / spacing)));
  const std::size_t nz = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround((d.z_top - d.z_base) / spacing)));
  const double dx = d.sx / static_cast<double>(nx), dy = d.sy / static_cast<double>(ny), dz = (d.z_top - d.z_base) / static_cast<double>(nz);
  auto at = [&](std::size_t i, std::size_t j, std::size_t k) {
    return Vec3{d.x0 + (static_cast<double>(i) + 0.5) * dx, d.y0 + (static_cast<double>(j) + 0.5) * dy, d.z_base + (static_cast<double>(k) + 0.5) * dz};
  };
  const std::size_t N = nx * ny * nz;
  std::vector<std::int32_t> unit(N);
  std::vector<BlockCode> block(N);
  // columns in parallel; each column evaluated as one batch
  parallel_for(nx * ny, 64, [&](std::size_t b, std::size_t e) {
    std::vector<Vec3> pts;
    EvalResult er;
    for (std::size_t c = b; c < e; ++c) {
      const std::size_t i = c % nx, j = c / nx;
      pts.clear();
      for (std::size_t k = 0; k < nz; ++k) pts.push_back(at(i, j, k));
      m.evaluate(pts, kEvalUnit | kEvalBlock, er);
      for (std::size_t k = 0; k < nz; ++k) {
        const std::size_t id = (k * ny + j) * nx + i;
        unit[id] = er.unit[k];
        block[id] = er.block[k];
      }
    }
  });
  std::map<std::pair<BlockCode, int>, std::size_t> cells;
  for (std::size_t id = 0; id < N; ++id) ++cells[{block[id], unit[id]}];
  ContactGraph g;
  g.spacing = spacing;
  for (auto& u : inf.units) g.unit_names.push_back(u.name);
  std::map<std::pair<BlockCode, int>, int> index;
  for (auto& [key, n] : cells) {
    index[key] = static_cast<int>(g.nodes.size());
    g.nodes.push_back({key.second, static_cast<int>(key.first), n});
  }
  // surface between two samples of different nodes; non-adjacent units in one block are resolved
  // along the segment between them (a unit thinner than the spacing is not skipped)
  std::map<std::tuple<int, int, std::string>, double> area;
  auto label = [&](BlockCode ba, int ua, BlockCode bb, int ub) -> std::string {
    if (ba != bb) {
      int f = lowest_bit(ba ^ bb);
      std::string nm = (f >= 0 && static_cast<std::size_t>(f) < inf.faults.size()) ? inf.faults[static_cast<std::size_t>(f)].name : std::to_string(f);
      return "fault:" + nm;
    }
    const int top = std::min(ua, ub), bot = std::max(ua, ub);
    // the boundary is the base horizon of the upper unit
    std::string hn = (top >= 0 && static_cast<std::size_t>(top) < inf.horizons.size()) ? inf.horizons[static_cast<std::size_t>(top)].name : "?";
    if (bot - top == 1) return "horizon:" + hn;
    // non-adjacent units: an unconformity when the series change or a unit between them terminates
    bool explained = false;
    for (int h = top; h < bot && h + 1 < static_cast<int>(inf.horizons.size()); ++h)
      explained = explained || inf.horizons[static_cast<std::size_t>(h)].series != inf.horizons[static_cast<std::size_t>(h + 1)].series;
    for (int u = top + 1; u < bot && u < static_cast<int>(inf.units.size()); ++u) explained = explained || !inf.units[static_cast<std::size_t>(u)].termination.empty();
    return (explained ? "unconformity:" : "horizon:") + hn;
  };
  auto node = [&](BlockCode b, int u) {
    auto it = index.find({b, u});
    if (it != index.end()) return it->second;
    // a unit seen only beside a fault or between lattice rows: a node without lattice cells
    index[{b, u}] = static_cast<int>(g.nodes.size());
    g.nodes.push_back({u, static_cast<int>(b), 0});
    return static_cast<int>(g.nodes.size()) - 1;
  };
  auto add = [&](BlockCode ba, int ua, BlockCode bb, int ub, double a) {
    const int na = node(ba, ua), nb = node(bb, ub);
    if (na == nb) return;
    area[{std::min(na, nb), std::max(na, nb), label(ba, ua, bb, ub)}] += a;
  };
  // fault crossings: located by bisection along the face segment, the units taken just either side
  // of the fault; horizontal faces also at eight heights between the lattice rows, so a contact
  // strip thinner than the spacing is still sampled
  struct Seg { Vec3 P, Q; double a; };
  std::vector<Seg> segs;
  auto pair = [&](std::size_t p, std::size_t q, const Vec3& P, const Vec3& Q, double a) {
    if (unit[p] == unit[q] && block[p] == block[q]) return;
    if (block[p] != block[q]) {
      const bool horizontal = std::fabs(P.z - Q.z) < 1e-9;
      if (!horizontal) { segs.push_back({P, Q, a}); return; }
      for (int k = -4; k <= 4; ++k) {
        const Vec3 off{0, 0, dz * static_cast<double>(k) / 9.0};
        segs.push_back({P + off, Q + off, a / 9.0});
      }
      return;
    }
    if (block[p] == block[q] && std::abs(unit[p] - unit[q]) > 1) {
      // walk the segment: every unit crossed in order
      const int ns = 64;
      std::vector<Vec3> s;
      for (int k = 0; k <= ns; ++k) s.push_back(P + (Q - P) * (static_cast<double>(k) / ns));
      EvalResult er;
      m.evaluate(s, kEvalUnit | kEvalBlock, er);
      BlockCode pb = er.block[0];
      int pu = er.unit[0];
      for (int k = 1; k <= ns; ++k) {
        if (er.unit[static_cast<std::size_t>(k)] == pu && er.block[static_cast<std::size_t>(k)] == pb) continue;
        add(pb, pu, er.block[static_cast<std::size_t>(k)], er.unit[static_cast<std::size_t>(k)], a);
        pu = er.unit[static_cast<std::size_t>(k)];
        pb = er.block[static_cast<std::size_t>(k)];
      }
      return;
    }
    add(block[p], unit[p], block[q], unit[q], a);
  };
  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i < nx; ++i) {
        const std::size_t id = (k * ny + j) * nx + i;
        if (i + 1 < nx) pair(id, id + 1, at(i, j, k), at(i + 1, j, k), dy * dz);
        if (j + 1 < ny) pair(id, id + nx, at(i, j, k), at(i, j + 1, k), dx * dz);
        if (k + 1 < nz) pair(id, id + nx * ny, at(i, j, k), at(i, j, k + 1), dx * dy);
      }
  {
    struct Side { BlockCode bp = 0, bq = 0; int up = -1, uq = -1; bool ok = false; };
    std::vector<Side> sides(segs.size());
    parallel_for(segs.size(), 256, [&](std::size_t b0, std::size_t b1) {
      const std::size_t n = b1 - b0;
      std::vector<Vec3> lo(n), hi(n), pts(2 * n);
      std::vector<BlockCode> bl(n);
      EvalResult er;
      for (std::size_t i = 0; i < n; ++i) { lo[i] = segs[b0 + i].P; hi[i] = segs[b0 + i].Q; pts[2 * i] = lo[i]; pts[2 * i + 1] = hi[i]; }
      m.evaluate(pts, kEvalBlock, er);
      std::vector<bool> active(n);
      for (std::size_t i = 0; i < n; ++i) {
        bl[i] = er.block[2 * i];
        active[i] = er.block[2 * i] != er.block[2 * i + 1];
      }
      std::vector<Vec3> mid(n);
      for (int it = 0; it < 12; ++it) {
        for (std::size_t i = 0; i < n; ++i) mid[i] = (lo[i] + hi[i]) * 0.5;
        m.evaluate(mid, kEvalBlock, er);
        for (std::size_t i = 0; i < n; ++i) (er.block[i] == bl[i] ? lo[i] : hi[i]) = mid[i];
      }
      for (std::size_t i = 0; i < n; ++i) { pts[2 * i] = lo[i]; pts[2 * i + 1] = hi[i]; }
      m.evaluate(pts, kEvalUnit | kEvalBlock, er);
      for (std::size_t i = 0; i < n; ++i)
        if (active[i] && er.block[2 * i] != er.block[2 * i + 1])
          sides[b0 + i] = {er.block[2 * i], er.block[2 * i + 1], er.unit[2 * i], er.unit[2 * i + 1], true};
    });
    for (std::size_t i = 0; i < segs.size(); ++i)
      if (sides[i].ok) add(sides[i].bp, sides[i].up, sides[i].bq, sides[i].uq, segs[i].a);
  }
  for (auto& [key, a] : area) g.edges.push_back({std::get<0>(key), std::get<1>(key), std::get<2>(key), a});
  return g;
}

bool same_topology(const ContactGraph& a, const ContactGraph& b, std::vector<std::string>* diff) {
  auto nodes = [](const ContactGraph& g) {
    std::set<std::pair<int, int>> s;
    for (auto& n : g.nodes) s.insert({n.block, n.unit});
    return s;
  };
  auto edges = [](const ContactGraph& g) {
    std::set<std::string> s;
    for (auto& e : g.edges) {
      auto& x = g.nodes[static_cast<std::size_t>(e.a)];
      auto& y = g.nodes[static_cast<std::size_t>(e.b)];
      std::string p = "unit " + std::to_string(x.unit) + " block " + std::to_string(x.block), q = "unit " + std::to_string(y.unit) + " block " + std::to_string(y.block);
      if (q < p) std::swap(p, q);
      s.insert(p + " - " + q + " across " + e.surface);
    }
    return s;
  };
  bool same = true;
  auto na = nodes(a), nb = nodes(b);
  for (auto& n : na)
    if (!nb.count(n)) { same = false; if (diff) diff->push_back("node only in the first: unit " + std::to_string(n.second) + " block " + std::to_string(n.first)); }
  for (auto& n : nb)
    if (!na.count(n)) { same = false; if (diff) diff->push_back("node only in the second: unit " + std::to_string(n.second) + " block " + std::to_string(n.first)); }
  auto ea = edges(a), eb = edges(b);
  for (auto& e : ea)
    if (!eb.count(e)) { same = false; if (diff) diff->push_back("edge only in the first: " + e); }
  for (auto& e : eb)
    if (!ea.count(e)) { same = false; if (diff) diff->push_back("edge only in the second: " + e); }
  return same;
}

WatertightReport watertightness(const Model& m, double spacing, double tolerance) {
  if (!(spacing > 0)) throw Error(ErrorCode::InvalidArgument, "watertightness: the spacing must be positive");
  const auto& d = m.info().domain;
  GridSpec g;
  g.dx = g.dy = spacing;
  g.nx = std::max<std::size_t>(2, static_cast<std::size_t>(std::floor(d.sx / spacing)) + 1);
  g.ny = std::max<std::size_t>(2, static_cast<std::size_t>(std::floor(d.sy / spacing)) + 1);
  g.x0 = d.x0 + 0.5 * (d.sx - static_cast<double>(g.nx - 1) * spacing);
  g.y0 = d.y0 + 0.5 * (d.sy - static_cast<double>(g.ny - 1) * spacing);
  std::map<int, repr::FaultedSurface> cache;
  return watertightness(m, g, [&](int h) -> const repr::FaultedSurface& {
    auto it = cache.find(h);
    if (it == cache.end()) it = cache.emplace(h, repr::horizon_grid_faulted(m, h, g)).first;
    return it->second;
  }, tolerance);
}

WatertightReport watertightness(const Model& m, const GridSpec& g, const std::function<const repr::FaultedSurface&(int)>& surface, double tolerance) {
  WatertightReport rep;
  rep.tolerance = tolerance;
  const auto& inf = m.info();
  const auto& d = inf.domain;
  const std::size_t nh = inf.horizons.size(), nf = inf.faults.size();
  // the horizon of block b at (x, y), followed from a nearby elevation: the root of its scalar field
  // inside block b nearest zprev, within +-win (NaN when the horizon has no root there in block b)
  auto z_track = [&](int h, BlockCode b, double x, double y, double zprev, double win, int ns) {
    std::vector<Vec3> pts(ns + 1);
    for (int k = 0; k <= ns; ++k) pts[static_cast<std::size_t>(k)] = {x, y, zprev - win + 2 * win * k / ns};
    EvalResult er;
    m.evaluate(pts, kEvalScalars | kEvalBlock, er);
    double best = kNaN, bestd = 1e300;
    auto eval1 = [&](double z, BlockCode* blk) {
      const Vec3 q{x, y, z};
      EvalResult e1;
      m.evaluate(std::span<const Vec3>(&q, 1), kEvalScalars | kEvalBlock, e1);
      if (blk) *blk = e1.block[0];
      return e1.s(static_cast<std::size_t>(h), 0);
    };
    for (int k = 0; k < ns; ++k) {
      const std::size_t a = static_cast<std::size_t>(k), c = a + 1;
      const bool ina = er.block[a] == b, inc = er.block[c] == b;
      if (!ina && !inc) continue;
      double za = pts[a].z, zc = pts[c].z;
      double sa = er.s(static_cast<std::size_t>(h), a), sc = er.s(static_cast<std::size_t>(h), c);
      if (ina != inc) {
        // the block boundary between the samples: keep the part inside block b
        double zin = ina ? za : zc, zout = ina ? zc : za;
        for (int it = 0; it < 30; ++it) {
          const double zm = 0.5 * (zin + zout);
          BlockCode bm;
          eval1(zm, &bm);
          (bm == b ? zin : zout) = zm;
        }
        const double sb = eval1(zin, nullptr);
        if (ina) { zc = zin; sc = sb; } else { za = zin; sa = sb; }
      }
      if (!(std::isfinite(sa) && std::isfinite(sc)) || (sa > 0) == (sc > 0)) continue;
      double lo = za, hi = zc, flo = sa;
      for (int it = 0; it < 24; ++it) {
        const Vec3 q{x, y, 0.5 * (lo + hi)};
        EvalResult e1;
        m.evaluate(std::span<const Vec3>(&q, 1), kEvalScalars, e1);
        const double f = e1.s(static_cast<std::size_t>(h), 0);
        if ((f > 0) == (flo > 0)) { lo = q.z; flo = f; } else { hi = q.z; }
      }
      const double z = 0.5 * (lo + hi);
      if (std::fabs(z - zprev) < bestd) { bestd = std::fabs(z - zprev); best = z; }
    }
    return best;
  };
  const double win = std::max(3 * g.dx, 100.0);
  // distance from an edge point to what may end a horizon: a fault, another horizon, the boundary
  auto gap_at = [&](int h, const Vec3& e, std::string* what) {
    double best = std::min({e.x - d.x0, d.x0 + d.sx - e.x, e.y - d.y0, d.y0 + d.sy - e.y, e.z - d.z_base, d.z_top - e.z});
    *what = "domain boundary";
    best = std::max(0.0, best);
    for (std::size_t f = 0; f < nf; ++f) {
      FaultEval fe;
      m.evaluate_fault(static_cast<int>(f), std::span<const Vec3>(&e, 1), fe);
      double v = std::fabs(fe.distance[0]);
      if (v < best) { best = v; *what = "fault:" + inf.faults[f].name; }
    }
    EvalResult er;
    m.evaluate(std::span<const Vec3>(&e, 1), kEvalScalars, er);
    for (std::size_t k = 0; k < nh; ++k) {
      if (static_cast<int>(k) == h) continue;
      double v = std::fabs(er.s(k, 0));
      if (v < best) { best = v; *what = "horizon:" + inf.horizons[k].name; }
    }
    return best;
  };
  for (std::size_t h = 0; h < nh; ++h) {
    const auto& fsurf = surface(static_cast<int>(h));
    // edge nodes: defined, with a 4-neighbour in a gap or in another block
    struct Edge { std::size_t n, nb; };
    std::vector<Edge> edges;
    for (std::size_t j = 0; j < g.ny; ++j)
      for (std::size_t i = 0; i < g.nx; ++i) {
        const std::size_t n = j * g.nx + i;
        if (std::isnan(fsurf.z.v[n])) continue;
        const BlockCode bn = fsurf.block[n];
        const long di[] = {1, -1, 0, 0}, dj[] = {0, 0, 1, -1};
        for (int q = 0; q < 4; ++q) {
          long ii = static_cast<long>(i) + di[q], jj = static_cast<long>(j) + dj[q];
          if (ii < 0 || jj < 0 || ii >= static_cast<long>(g.nx) || jj >= static_cast<long>(g.ny)) continue;
          const std::size_t nb = static_cast<std::size_t>(jj) * g.nx + static_cast<std::size_t>(ii);
          if (std::isnan(fsurf.z.v[nb]) || fsurf.block[nb] != bn) edges.push_back({n, nb});
        }
      }
    std::vector<WatertightGap> found(edges.size());
    parallel_for(edges.size(), 16, [&](std::size_t b0, std::size_t b1) {
      for (std::size_t q = b0; q < b1; ++q) {
        const auto& E = edges[q];
        const BlockCode b = fsurf.block[E.n];
        const double xa = g.x(E.n % g.nx), ya = g.y(E.n / g.nx), xb = g.x(E.nb % g.nx), yb = g.y(E.nb / g.nx);
        // where the faulted surface keeps another block's sheet above (reverse-fault overlap), the
        // horizon of block b continues beneath it: not an edge
        {
          const double z1 = z_track(static_cast<int>(h), b, xb, yb, fsurf.z.v[E.n], win, 48);
          if (std::isfinite(z1)) {
            found[q] = {inf.horizons[h].name + " (continues under an overlapping block)", {xb, yb, z1}, 0.0};
            continue;
          }
        }
        // bisect along the horizon of block b for the point where it leaves the block
        double lo = 0, hi = 1;
        Vec3 last{xa, ya, fsurf.z.v[E.n]};
        const double seg = std::hypot(xb - xa, yb - ya);
        for (int it = 0; it < 12; ++it) {
          const double t = 0.5 * (lo + hi);
          // the horizon moves little over one bisection step: a window for slopes up to about 84 degrees
          const double w = std::max(5.0, 10.0 * seg * (t - lo));
          const double x = xa + t * (xb - xa), y = ya + t * (yb - ya), z = z_track(static_cast<int>(h), b, x, y, last.z, w, 8);
          if (std::isfinite(z)) {
            lo = t;
            last = {x, y, z};
          } else {
            hi = t;
          }
        }
        std::string what;
        const double gap = gap_at(static_cast<int>(h), last, &what);
        found[q] = {inf.horizons[h].name + " (nearest: " + what + ")", last, gap};
      }
    });
    for (auto& w : found) {
      rep.max_gap = std::max(rep.max_gap, w.gap);
      if (w.gap > tolerance) rep.gaps.push_back(w);
    }
  }
  std::sort(rep.gaps.begin(), rep.gaps.end(), [](const WatertightGap& a, const WatertightGap& b) { return a.gap > b.gap; });
  if (rep.gaps.size() > 20) rep.gaps.resize(20);
  return rep;
}

}  // namespace stratonix
