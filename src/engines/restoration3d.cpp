// Stage 4B (STX-TN-15 section 4): flexural-slip and volumetric restoration, compaction across
// moving faults, flexural backstripping and eroded-section estimates.
#include "stratonix/engines/restoration3d.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <cmath>
#include <map>
#include <numeric>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/structure.hpp"

namespace stratonix {

// ---- 4.1 flexural slip -------------------------------------------------------------------------
UnfoldResult flexural_slip_unfold(const std::vector<Vec3>& h, double pin_x, double datum_z) {
  if (h.size() < 2) throw Error(ErrorCode::InvalidArgument, "flexural slip: a horizon needs two points");
  // arc length along the horizon in the section plane (x, z); the pin is the point nearest pin_x
  std::vector<double> s{0};
  for (std::size_t i = 1; i < h.size(); ++i) s.push_back(s.back() + std::hypot(h[i].x - h[i - 1].x, h[i].z - h[i - 1].z));
  std::size_t ip = 0;
  for (std::size_t i = 1; i < h.size(); ++i)
    if (std::fabs(h[i].x - pin_x) < std::fabs(h[ip].x - pin_x)) ip = i;
  UnfoldResult r;
  r.present_length = s.back();
  // bed length measured from the pin is laid out flat on the datum (Dahlstrom 1969)
  for (std::size_t i = 0; i < h.size(); ++i) r.restored.push_back({h[ip].x + (s[i] - s[ip]), h[i].y, datum_z});
  r.restored_length = r.restored.back().x - r.restored.front().x;
  return r;
}

// ---- 4.2 volumetric restoration -----------------------------------------------------------------
namespace {
using M3 = std::array<std::array<double, 3>, 3>;
bool inv3(const M3& a, M3& o, double& det) {
  det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
  if (!(std::fabs(det) > 1e-300)) return false;
  const double id = 1 / det;
  o[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * id;
  o[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * id;
  o[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * id;
  o[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * id;
  o[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * id;
  o[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * id;
  o[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * id;
  o[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * id;
  o[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * id;
  return true;
}
/// Eigenvalues of a symmetric 3 x 3 matrix (Jacobi rotations).
std::array<double, 3> sym_eig(M3 a) {
  for (int sweep = 0; sweep < 50; ++sweep) {
    double off = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
    if (off < 1e-15) break;
    for (int p = 0; p < 2; ++p)
      for (int q = p + 1; q < 3; ++q) {
        if (std::fabs(a[p][q]) < 1e-300) continue;
        const double th = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
        const double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1));
        const double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < 3; ++k) {
          const double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (int k = 0; k < 3; ++k) {
          const double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
      }
  }
  std::array<double, 3> e{a[0][0], a[1][1], a[2][2]};
  std::sort(e.begin(), e.end());
  return e;
}

/// Symmetric block-sparse matrix (3 x 3 blocks per node pair) with a fixed pattern.
struct BlockCSR {
  std::vector<std::size_t> row;   ///< n + 1
  std::vector<std::size_t> col;   ///< node columns, sorted per row
  std::vector<double> val;        ///< 9 per entry
  std::size_t n = 0;
  void build(std::size_t nn, std::vector<std::vector<std::size_t>>& nb) {
    n = nn;
    row.assign(n + 1, 0);
    col.clear();
    for (std::size_t i = 0; i < n; ++i) {
      auto& v = nb[i];
      v.push_back(i);
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
      row[i + 1] = row[i] + v.size();
      col.insert(col.end(), v.begin(), v.end());
    }
    val.assign(9 * col.size(), 0.0);
  }
  double* at(std::size_t i, std::size_t j) {
    auto b = col.begin() + static_cast<std::ptrdiff_t>(row[i]), e = col.begin() + static_cast<std::ptrdiff_t>(row[i + 1]);
    auto it = std::lower_bound(b, e, j);
    return &val[9 * static_cast<std::size_t>(it - col.begin())];
  }
  void zero() { std::fill(val.begin(), val.end(), 0.0); }
  void mul(const std::vector<double>& x, std::vector<double>& y) const {
    parallel_for(n, 512, [&](std::size_t i0, std::size_t i1) {
      for (std::size_t i = i0; i < i1; ++i) {
        double a = 0, b = 0, c = 0;
        for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
          const double* m = &val[9 * k];
          const double* v = &x[3 * col[k]];
          a += m[0] * v[0] + m[1] * v[1] + m[2] * v[2];
          b += m[3] * v[0] + m[4] * v[1] + m[5] * v[2];
          c += m[6] * v[0] + m[7] * v[1] + m[8] * v[2];
        }
        y[3 * i] = a;
        y[3 * i + 1] = b;
        y[3 * i + 2] = c;
      }
    });
  }
};

/// Preconditioned conjugate gradients with 3 x 3 block-Jacobi.
std::size_t pcg(const BlockCSR& K, const std::vector<double>& f, std::vector<double>& x, double tol, std::size_t maxit) {
  const std::size_t n = K.n, N = 3 * n;
  std::vector<M3> pinv(n);
  for (std::size_t i = 0; i < n; ++i) {
    M3 d{};
    for (std::size_t k = K.row[i]; k < K.row[i + 1]; ++k)
      if (K.col[k] == i)
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) d[a][b] = K.val[9 * k + 3 * a + b];
    double det;
    if (!inv3(d, pinv[i], det)) pinv[i] = M3{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  }
  auto prec = [&](const std::vector<double>& r, std::vector<double>& z) {
    for (std::size_t i = 0; i < n; ++i)
      for (int a = 0; a < 3; ++a) z[3 * i + a] = pinv[i][a][0] * r[3 * i] + pinv[i][a][1] * r[3 * i + 1] + pinv[i][a][2] * r[3 * i + 2];
  };
  std::vector<double> r(N), z(N), p(N), q(N);
  K.mul(x, q);
  for (std::size_t i = 0; i < N; ++i) r[i] = f[i] - q[i];
  prec(r, z);
  p = z;
  double rz = std::inner_product(r.begin(), r.end(), z.begin(), 0.0);
  const double fn = std::sqrt(std::inner_product(f.begin(), f.end(), f.begin(), 0.0));
  if (fn == 0) return 0;
  std::size_t it = 0;
  for (; it < maxit; ++it) {
    if (std::sqrt(std::inner_product(r.begin(), r.end(), r.begin(), 0.0)) <= tol * fn) break;
    K.mul(p, q);
    const double a = rz / std::inner_product(p.begin(), p.end(), q.begin(), 0.0);
    for (std::size_t i = 0; i < N; ++i) { x[i] += a * p[i]; r[i] -= a * q[i]; }
    prec(r, z);
    const double rz2 = std::inner_product(r.begin(), r.end(), z.begin(), 0.0);
    const double b = rz2 / rz;
    rz = rz2;
    for (std::size_t i = 0; i < N; ++i) p[i] = z[i] + b * p[i];
  }
  return it;
}

int lowest_set_bit(BlockCode b) {
  for (int i = 0; i < 64; ++i)
    if ((b >> i) & 1u) return i;
  return -1;
}
}  // namespace

VolumetricResult volumetric_restore(const Model& m, const VolumetricSpec& s) {
  VolumetricResult R;
  const ModelInfo& info = m.info();
  int hz = -1;
  for (std::size_t h = 0; h < info.horizons.size(); ++h)
    if (info.horizons[h].name == s.horizon) hz = static_cast<int>(h);
  if (hz < 0) { R.reason = "unknown horizon '" + s.horizon + "'"; return R; }
  if (!(s.poisson >= 0 && s.poisson < 0.5) || !(s.young > 0)) { R.reason = "Poisson's ratio in [0, 0.5) and a positive Young's modulus"; return R; }
  for (auto& [E, nu] : s.unit_props)
    if (!(nu >= 0 && nu < 0.5) || !(E > 0)) { R.reason = "Poisson's ratio in [0, 0.5) and a positive Young's modulus for every unit"; return R; }
  const Domain& d = info.domain;
  const double h = s.mesh_size > 0 ? s.mesh_size : std::max(d.sx, d.sy) / 30;
  const std::size_t nx = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(d.sx / h))) + 1;
  const std::size_t ny = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(d.sy / h))) + 1;
  const std::size_t nz = std::max<std::size_t>(4, static_cast<std::size_t>(std::llround((d.z_top - d.z_base) / h))) + 1;
  const double dx = d.sx / (nx - 1), dy = d.sy / (ny - 1), dz = (d.z_top - d.z_base) / (nz - 1);
  auto lid = [&](std::size_t i, std::size_t j, std::size_t k) { return (k * ny + j) * nx + i; };
  auto lpos = [&](std::size_t id) {
    const std::size_t i = id % nx, j = (id / nx) % ny, k = id / (nx * ny);
    return Vec3{d.x0 + i * dx, d.y0 + j * dy, d.z_base + k * dz};
  };
  // tetrahedra: six per cell (Kuhn), each a path from corner 0 to corner 7
  struct Tet { std::array<std::size_t, 4> node; BlockCode block = 0; int unit = -1; };
  std::vector<Tet> tets;
  static const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  for (std::size_t k = 0; k + 1 < nz; ++k)
    for (std::size_t j = 0; j + 1 < ny; ++j)
      for (std::size_t i = 0; i + 1 < nx; ++i)
        for (auto& pm : perm) {
          std::array<std::size_t, 3> c{i, j, k};
          Tet t;
          t.node[0] = lid(c[0], c[1], c[2]);
          for (int a = 0; a < 3; ++a) {
            c[static_cast<std::size_t>(pm[a])] += 1;
            t.node[static_cast<std::size_t>(a) + 1] = lid(c[0], c[1], c[2]);
          }
          tets.push_back(t);
        }
  {
    std::vector<Vec3> cen(tets.size());
    for (std::size_t t = 0; t < tets.size(); ++t) {
      Vec3 c{0, 0, 0};
      for (auto n : tets[t].node) c = c + lpos(n);
      cen[t] = c * 0.25;
    }
    EvalResult er;
    m.evaluate(cen, kEvalBlock | kEvalUnit, er);
    for (std::size_t t = 0; t < tets.size(); ++t) { tets[t].block = er.block[t]; tets[t].unit = er.unit[t]; }
  }
  // node copies: one per (lattice node, fault block) among the tetrahedra that use the node
  std::map<std::pair<std::size_t, BlockCode>, std::size_t> copy_of;
  std::vector<std::size_t> copy_node;
  std::vector<BlockCode> copy_block;
  std::vector<std::array<std::size_t, 4>> tc(tets.size());
  for (std::size_t t = 0; t < tets.size(); ++t)
    for (int a = 0; a < 4; ++a) {
      auto key = std::make_pair(tets[t].node[static_cast<std::size_t>(a)], tets[t].block);
      auto it = copy_of.find(key);
      if (it == copy_of.end()) {
        it = copy_of.emplace(key, copy_node.size()).first;
        copy_node.push_back(key.first);
        copy_block.push_back(key.second);
      }
      tc[t][static_cast<std::size_t>(a)] = it->second;
    }
  const std::size_t nc = copy_node.size();
  R.nodes = nc;
  R.tetrahedra = tets.size();
  R.mesh_spacing = std::max(dx, dy);
  // horizon scalar at each copy, in its block
  std::vector<double> sc(nc, kNaN);
  {
    std::map<BlockCode, std::vector<std::size_t>> by_block;
    for (std::size_t c = 0; c < nc; ++c) by_block[copy_block[c]].push_back(c);
    for (auto& [b, ids] : by_block) {
      std::vector<Vec3> p(ids.size());
      for (std::size_t q = 0; q < ids.size(); ++q) p[q] = lpos(copy_node[ids[q]]);
      std::vector<double> v(ids.size());
      m.scalar_in_block(p, b, hz, v);
      for (std::size_t q = 0; q < ids.size(); ++q) sc[ids[q]] = v[q];
    }
  }
  // the model at the restored horizon's age: material wholly above the horizon is removed
  std::vector<char> keep(tets.size(), 1), used(nc, 0);
  for (std::size_t t = 0; t < tets.size(); ++t) {
    bool above = true;
    for (auto q : tc[t]) above = above && std::isfinite(sc[q]) && sc[q] > 0;
    keep[t] = above ? 0 : 1;
    if (keep[t])
      for (auto q : tc[t]) used[q] = 1;
  }
  // horizon crossings on tetrahedron edges: (copy a, copy b, t), each edge once
  struct Cross { std::size_t a, b; double t; };
  std::vector<Cross> cross;
  std::vector<char> at_h(tets.size(), 0);
  {
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> seen;
    static const int E[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    for (std::size_t t = 0; t < tets.size(); ++t)
      for (auto& e : E) {
        if (!keep[t]) break;
        std::size_t a = tc[t][static_cast<std::size_t>(e[0])], b = tc[t][static_cast<std::size_t>(e[1])];
        const double sa = sc[a], sb = sc[b];
        if (!std::isfinite(sa) || !std::isfinite(sb) || (sa > 0) == (sb > 0) || sa == sb) continue;
        at_h[t] = 1;
        if (a > b) std::swap(a, b);
        if (seen.count({a, b})) continue;
        seen[{a, b}] = cross.size();
        const double ta = sc[a] / (sc[a] - sc[b]);
        cross.push_back({a, b, ta});
      }
  }
  if (cross.empty()) { R.reason = "the horizon " + s.horizon + " does not cross the mesh"; return R; }
  // fault contacts: copies of one lattice node in different blocks (sliding: normal gap held)
  struct Contact { std::size_t a, b; Vec3 n; };
  std::vector<Contact> contacts;
  {
    std::map<std::size_t, std::vector<std::size_t>> by_node;
    for (std::size_t c = 0; c < nc; ++c) by_node[copy_node[c]].push_back(c);
    for (auto& [nd, cs] : by_node)
      for (std::size_t p = 0; p + 1 < cs.size(); ++p)
        for (std::size_t q = p + 1; q < cs.size(); ++q) {
          if (!used[cs[p]] || !used[cs[q]]) continue;
          const int f = lowest_set_bit(copy_block[cs[p]] ^ copy_block[cs[q]]);
          if (f < 0 || static_cast<std::size_t>(f) >= info.faults.size()) continue;
          contacts.push_back({cs[p], cs[q], info.faults[static_cast<std::size_t>(f)].normal});
        }
  }
  R.fault_contacts = contacts.size();
  // pattern
  BlockCSR K;
  {
    std::vector<std::vector<std::size_t>> nb(nc);
    for (std::size_t t = 0; t < tc.size(); ++t)
      for (int a = 0; a < 4 && keep[t]; ++a)
        for (int b = 0; b < 4; ++b) nb[tc[t][static_cast<std::size_t>(a)]].push_back(tc[t][static_cast<std::size_t>(b)]);
    for (auto& c : contacts) { nb[c.a].push_back(c.b); nb[c.b].push_back(c.a); }
    for (auto& c : cross) { nb[c.a].push_back(c.b); nb[c.b].push_back(c.a); }
    K.build(nc, nb);
  }
  // pin: the lattice column nearest the pin, held in plan
  const double px = s.pin_given ? s.pin_x : d.x0 + 0.5 * d.sx, py = s.pin_given ? s.pin_y : d.y0 + 0.5 * d.sy;
  std::vector<Vec3> X(nc), x(nc);
  for (std::size_t c = 0; c < nc; ++c) X[c] = x[c] = lpos(copy_node[c]);
  // the pin acts on the horizon point nearest the pin line (an interpolated crossing): nodes held
  // above and below it would forbid the rotation that unfolds a dipping limb
  std::size_t pin_cross = 0;
  {
    double best = 1e300;
    for (std::size_t q = 0; q < cross.size(); ++q) {
      const Vec3 p = X[cross[q].a] * (1 - cross[q].t) + X[cross[q].b] * cross[q].t;
      const double dd = std::hypot(p.x - px, p.y - py);
      if (dd < best) { best = dd; pin_cross = q; }
    }
  }
  // material per tetrahedron
  auto mat = [&](int unit) -> std::pair<double, double> {
    if (unit >= 0 && static_cast<std::size_t>(unit) < s.unit_props.size()) return s.unit_props[static_cast<std::size_t>(unit)];
    return {s.young, s.poisson};
  };
  double Emax = s.young;
  for (auto& [E, nu] : s.unit_props) Emax = std::max(Emax, E);
  const double kpen = 100 * Emax * h, kreg = 1e-7 * Emax * h;
  const int steps = std::max(1, s.steps);
  std::vector<double> u(3 * nc, 0.0), f(3 * nc, 0.0);
  for (int step = 1; step <= steps; ++step) {
    K.zero();
    std::fill(f.begin(), f.end(), 0.0);
    // elastic stiffness on the current geometry
    for (std::size_t t = 0; t < tets.size(); ++t) {
      if (!keep[t]) continue;
      const auto& n = tc[t];
      M3 A;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
          const Vec3 e = x[n[static_cast<std::size_t>(b) + 1]] - x[n[0]];
          A[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = a == 0 ? e.x : a == 1 ? e.y : e.z;
        }
      M3 Ai;
      double det;
      if (!inv3(A, Ai, det)) continue;
      const double V = std::fabs(det) / 6;
      // gradients of the shape functions: rows of A^-1 for nodes 1..3, node 0 minus their sum
      double g[4][3];
      for (int b = 0; b < 3; ++b)
        for (int a = 0; a < 3; ++a) g[b + 1][a] = Ai[static_cast<std::size_t>(b)][static_cast<std::size_t>(a)];
      for (int a = 0; a < 3; ++a) g[0][a] = -(g[1][a] + g[2][a] + g[3][a]);
      const auto [E, nu] = mat(tets[t].unit);
      const double lam = E * nu / ((1 + nu) * (1 - 2 * nu)), mu = E / (2 * (1 + nu));
      // K_ab (3 x 3) = V [ lam g_a g_b^T + mu (g_b g_a^T + (g_a . g_b) I) ]
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
          double* blk = K.at(n[static_cast<std::size_t>(a)], n[static_cast<std::size_t>(b)]);
          const double gg = g[a][0] * g[b][0] + g[a][1] * g[b][1] + g[a][2] * g[b][2];
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) blk[3 * i + j] += V * (lam * g[a][i] * g[b][j] + mu * (g[a][j] * g[b][i] + (i == j ? gg : 0.0)));
        }
    }
    // the horizon crossings move a share of the way to the datum (the last increment reaches it)
    for (auto& c : cross) {
      const double w[2] = {1 - c.t, c.t};
      const std::size_t id[2] = {c.a, c.b};
      const double zc = w[0] * x[c.a].z + w[1] * x[c.b].z;
      const double dzt = (s.datum_z - zc) / static_cast<double>(steps - step + 1);
      for (int p = 0; p < 2; ++p) {
        f[3 * id[p] + 2] += kpen * w[p] * dzt;
        for (int q = 0; q < 2; ++q) K.at(id[p], id[q])[8] += kpen * w[p] * w[q];
      }
    }
    // sliding fault contacts: the normal component of the copies' relative motion is held
    for (auto& c : contacts) {
      const double nn[3] = {c.n.x, c.n.y, c.n.z};
      for (int p = 0; p < 2; ++p)
        for (int q = 0; q < 2; ++q) {
          const double sgn = p == q ? 1.0 : -1.0;
          double* blk = K.at(p == 0 ? c.a : c.b, q == 0 ? c.a : c.b);
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) blk[3 * i + j] += sgn * kpen * nn[i] * nn[j];
        }
    }
    for (std::size_t c = 0; c < nc; ++c) {
      double* blk = K.at(c, c);
      for (int i = 0; i < 3; ++i) blk[4 * i] += kreg;
    }
    {
      // pin: the horizon point's horizontal displacement is held at zero
      const auto& c = cross[pin_cross];
      const double w[2] = {1 - c.t, c.t};
      const std::size_t id[2] = {c.a, c.b};
      for (int p = 0; p < 2; ++p)
        for (int q = 0; q < 2; ++q) {
          double* blk = K.at(id[p], id[q]);
          blk[0] += kpen * w[p] * w[q];
          blk[4] += kpen * w[p] * w[q];
        }
    }
    std::vector<double> du(3 * nc, 0.0);
    R.cg_iterations += pcg(K, f, du, 1e-10, 20000);
    for (std::size_t c = 0; c < nc; ++c) x[c] = x[c] + Vec3{du[3 * c], du[3 * c + 1], du[3 * c + 2]};
  }
  // flatness of the restored horizon
  double ss = 0;
  for (auto& c : cross) {
    const double z = (1 - c.t) * x[c.a].z + c.t * x[c.b].z;
    ss += (z - s.datum_z) * (z - s.datum_z);
  }
  R.flatness_rms = std::sqrt(ss / static_cast<double>(cross.size()));
  // Green-Lagrange strain of the total deformation, per tetrahedron
  double sum = 0;
  std::size_t cnt = 0;
  R.tetrahedra = 0;
  for (std::size_t t = 0; t < tets.size(); ++t) {
    if (!keep[t]) continue;
    ++R.tetrahedra;
    const auto& n = tc[t];
    M3 A, B;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        const Vec3 e0 = X[n[static_cast<std::size_t>(b) + 1]] - X[n[0]], e1 = x[n[static_cast<std::size_t>(b) + 1]] - x[n[0]];
        A[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = a == 0 ? e0.x : a == 1 ? e0.y : e0.z;
        B[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = a == 0 ? e1.x : a == 1 ? e1.y : e1.z;
      }
    M3 Ai;
    double det;
    if (!inv3(A, Ai, det)) continue;
    M3 F{}, C{};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) F[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] += B[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] * Ai[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) C[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] += F[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] * F[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
        C[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = 0.5 * (C[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] - (i == j ? 1.0 : 0.0));
      }
    const auto e = sym_eig(C);
    double detF;
    M3 Fi;
    inv3(F, Fi, detF);
    StrainSample sm;
    for (auto q : n) { sm.present = sm.present + X[q] * 0.25; sm.restored = sm.restored + x[q] * 0.25; }
    sm.max_abs = std::max(std::fabs(e[0]), std::fabs(e[2]));
    sm.max_shear = e[2] - e[0];
    sm.dilatation = detF - 1;
    sm.at_horizon = at_h[t] != 0;
    sm.unit = tets[t].unit;
    R.max_strain = std::max(R.max_strain, sm.max_abs);
    sum += sm.max_abs;
    ++cnt;
    if (sm.at_horizon) {
      // stretch in the bedding plane: the plane normal is the gradient of the horizon scalar in the
      // tetrahedron (present geometry); the principal stretches of F restricted to that plane
      Vec3 gr{0, 0, 0};
      double gN[4][3];
      for (int b = 0; b < 3; ++b)
        for (int a = 0; a < 3; ++a) gN[b + 1][a] = Ai[static_cast<std::size_t>(b)][static_cast<std::size_t>(a)];
      for (int a = 0; a < 3; ++a) gN[0][a] = -(gN[1][a] + gN[2][a] + gN[3][a]);
      for (int q = 0; q < 4; ++q) gr = gr + Vec3{gN[q][0], gN[q][1], gN[q][2]} * sc[n[static_cast<std::size_t>(q)]];
      if (norm(gr) > 0) {
        const Vec3 nn = normalized(gr);
        Vec3 t1 = normalized(std::fabs(nn.z) < 0.9 ? stratonix::cross(nn, Vec3{0, 0, 1}) : stratonix::cross(nn, Vec3{1, 0, 0}));
        Vec3 t2 = stratonix::cross(nn, t1);
        auto Fv = [&](const Vec3& v) {
          return Vec3{F[0][0] * v.x + F[0][1] * v.y + F[0][2] * v.z, F[1][0] * v.x + F[1][1] * v.y + F[1][2] * v.z, F[2][0] * v.x + F[2][1] * v.y + F[2][2] * v.z};
        };
        const Vec3 a1 = Fv(t1), a2 = Fv(t2);
        const double g11 = dot(a1, a1), g22 = dot(a2, a2), g12 = dot(a1, a2);
        const double tr = 0.5 * (g11 + g22), dd = std::sqrt(std::max(0.0, 0.25 * (g11 - g22) * (g11 - g22) + g12 * g12));
        const double l1 = std::sqrt(std::max(0.0, tr + dd)) - 1, l2 = std::sqrt(std::max(0.0, tr - dd)) - 1;
        sm.bed_stretch = std::max(std::fabs(l1), std::fabs(l2));
        R.horizon_stretch_max = std::max(R.horizon_stretch_max, sm.bed_stretch);
      }
    }
    R.samples.push_back(sm);
  }
  R.mean_strain = cnt ? sum / static_cast<double>(cnt) : 0.0;
  R.ok = true;
  return R;
}

// ---- 4.3 to 4.5 --------------------------------------------------------------------------------
double block_porosity(double phi0, double c, double burial_depth) {
  if (!(phi0 > 0 && phi0 < 1) || c < 0) throw Error(ErrorCode::InvalidArgument, "Athy law: surface porosity in (0, 1), decay >= 0");
  return phi0 * std::exp(-c * std::max(0.0, burial_depth));
}

double airy_tectonic_subsidence(double S, double rho_s, double rho_m, double rho_w, double Wd, double dSL) {
  if (!(rho_m > rho_s && rho_m > rho_w)) throw Error(ErrorCode::InvalidArgument, "Airy backstripping: the mantle must be denser than sediment and water");
  // Steckler and Watts (1978)
  return S * (rho_m - rho_s) / (rho_m - rho_w) + Wd - dSL * rho_m / (rho_m - rho_w);
}

double flexural_response(double k, double te, double young, double poisson, double delta_rho) {
  if (te < 0 || te > 100e3) throw Error(ErrorCode::InvalidArgument, "flexure: Te from 0 to 100 km");
  if (te == 0) return 1.0;  // Airy
  const double D = young * te * te * te / (12 * (1 - poisson * poisson));
  return 1 / (1 + D * std::pow(k, 4) / (delta_rho * 9.81));  // Watts (2001)
}

double eroded_section(const std::vector<std::pair<double, double>>& log, double phi0, double c) {
  if (log.empty() || !(c > 0) || !(phi0 > 0)) throw Error(ErrorCode::InvalidArgument, "eroded section: a porosity log, phi0 > 0, c > 0");
  // each sample's equivalent maximum burial on the normal trend minus its present depth (Magara 1976)
  double sum = 0;
  std::size_t n = 0;
  for (auto& [z, phi] : log) {
    if (!(phi > 0)) continue;
    sum += -std::log(phi / phi0) / c - z;
    ++n;
  }
  if (n == 0) throw Error(ErrorCode::InvalidArgument, "eroded section: no positive porosity in the log");
  return std::max(0.0, sum / static_cast<double>(n));
}

// ---- 4.4 backstripping in models ----------------------------------------------------------------
namespace {
double step_value(const std::vector<double>& per_h, int h, double fallback) {
  if (h >= 0 && static_cast<std::size_t>(h) < per_h.size() && std::isfinite(per_h[static_cast<std::size_t>(h)])) return per_h[static_cast<std::size_t>(h)];
  return fallback;
}
}  // namespace

std::vector<TectonicStep> backstrip_column(const Structure& st, const BackstripSpec& s, double x, double y) {
  std::vector<TectonicStep> out;
  auto hist = st.burial_history(x, y);
  if (hist.empty() || s.basement < 0) return out;
  const std::size_t b = static_cast<std::size_t>(s.basement);
  for (std::size_t k = 0; k < hist.size(); ++k) {
    const auto& stp = hist[k];
    if (b >= stp.depth.size() || !std::isfinite(stp.depth[b])) continue;
    TectonicStep t;
    t.label = stp.label;
    t.age = stp.age;
    // the horizon at the surface at this step (the present step: the top horizon's present depth)
    int top = -1;
    for (std::size_t h = 0; h <= b; ++h)
      if (std::isfinite(stp.depth[h])) { top = static_cast<int>(h); break; }
    if (top < 0) continue;
    const double d0 = stp.depth[static_cast<std::size_t>(top)];
    t.sediment = stp.depth[b];
    double mass = 0;
    if (d0 > 0) {
      // overburden above the top horizon (present step): compacted as the top unit
      double phi = st.mean_porosity(top + 1, 0.0, d0);
      if (!std::isfinite(phi)) phi = 0;
      mass += d0 * (s.rho_g * (1 - phi) + s.rho_w * phi);
    }
    for (std::size_t u = static_cast<std::size_t>(top) + 1; u <= b; ++u) {
      const double z1 = stp.depth[u - 1], z2 = stp.depth[u];
      if (!(z2 > z1)) continue;
      double phi = st.mean_porosity(static_cast<int>(u), z1, z2);
      if (!std::isfinite(phi)) phi = 0;  // a unit without a compaction law: no pore space
      mass += (z2 - z1) * (s.rho_g * (1 - phi) + s.rho_w * phi);
    }
    t.density = t.sediment > 0 ? mass / t.sediment : s.rho_g;
    const bool present = stp.label == "present";
    t.water_depth = present ? s.water_depth_present : step_value(s.water_depth, top, 0.0);
    t.sea_level = present ? s.sea_level_present : step_value(s.sea_level, top, 0.0);
    t.airy = airy_tectonic_subsidence(t.sediment, std::min(t.density, s.rho_m - 1), s.rho_m, s.rho_w, t.water_depth, t.sea_level);
    out.push_back(t);
  }
  return out;
}

namespace {
/// In-place 1D DFT of length n (naive, n small) along a strided line.
void dft_line(std::vector<std::complex<double>>& a, std::size_t off, std::size_t stride, std::size_t n, bool inverse) {
  std::vector<std::complex<double>> o(n);
  const double sg = inverse ? 1.0 : -1.0;
  for (std::size_t k = 0; k < n; ++k) {
    std::complex<double> acc = 0;
    for (std::size_t j = 0; j < n; ++j) acc += a[off + j * stride] * std::polar(1.0, sg * 2 * kPi * static_cast<double>(k * j % n) / static_cast<double>(n));
    o[k] = inverse ? acc / static_cast<double>(n) : acc;
  }
  for (std::size_t k = 0; k < n; ++k) a[off + k * stride] = o[k];
}
}  // namespace

BackstripMaps backstrip_maps(const Structure& st, const BackstripSpec& s, const Domain& dom, std::size_t n) {
  BackstripMaps M;
  M.n = n;
  std::vector<std::vector<TectonicStep>> cols(n * n);
  parallel_for(n * n, 8, [&](std::size_t a0, std::size_t a1) {
    for (std::size_t q = a0; q < a1; ++q) {
      const double x = dom.x0 + (static_cast<double>(q % n) + 0.5) * dom.sx / static_cast<double>(n);
      const double y = dom.y0 + (static_cast<double>(q / n) + 0.5) * dom.sy / static_cast<double>(n);
      cols[q] = backstrip_column(st, s, x, y);
    }
  });
  const std::size_t ns = cols[0].size();
  for (auto& c : cols)
    if (c.size() != ns) return M;  // the basement is missing somewhere: no maps
  for (std::size_t k = 0; k < ns; ++k) {
    M.labels.push_back(cols[0][k].label);
    M.ages.push_back(cols[0][k].age);
    std::vector<double> airy(n * n), flex(n * n, kNaN);
    for (std::size_t q = 0; q < n * n; ++q) airy[q] = cols[q][k].airy;
    if (s.te > 0) {
      // the sediment load's Airy deflection, filtered by the plate response (Watts 2001), on a grid
      // padded to twice its size by extending the edges
      const std::size_t N = 2 * n;
      std::vector<std::complex<double>> L(N * N);
      double rho_mean = 0;
      for (std::size_t q = 0; q < n * n; ++q) rho_mean += cols[q][k].density / static_cast<double>(n * n);
      for (std::size_t j = 0; j < N; ++j)
        for (std::size_t i = 0; i < N; ++i) {
          const std::size_t ii = std::min(n - 1, i < n + n / 2 ? std::min(i, n - 1) : 0), jj = std::min(n - 1, j < n + n / 2 ? std::min(j, n - 1) : 0);
          const auto& c = cols[jj * n + ii][k];
          L[j * N + i] = c.sediment * (c.density - s.rho_w) / (s.rho_m - s.rho_w);
        }
      for (std::size_t j = 0; j < N; ++j) dft_line(L, j * N, 1, N, false);
      for (std::size_t i = 0; i < N; ++i) dft_line(L, i, N, N, false);
      const double dx = dom.sx / static_cast<double>(n), dy = dom.sy / static_cast<double>(n);
      for (std::size_t j = 0; j < N; ++j)
        for (std::size_t i = 0; i < N; ++i) {
          const double fi = static_cast<double>(i <= N / 2 ? i : N - i) / (static_cast<double>(N) * dx);
          const double fj = static_cast<double>(j <= N / 2 ? j : N - j) / (static_cast<double>(N) * dy);
          const double kk = 2 * kPi * std::hypot(fi, fj);
          L[j * N + i] *= flexural_response(kk, s.te, s.young, s.poisson, s.rho_m - std::min(rho_mean, s.rho_m - 1));
        }
      for (std::size_t i = 0; i < N; ++i) dft_line(L, i, N, N, true);
      for (std::size_t j = 0; j < N; ++j) dft_line(L, j * N, 1, N, true);
      for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i) {
          const auto& c = cols[j * n + i][k];
          const double w = L[j * N + i].real();
          flex[j * n + i] = c.sediment - w + c.water_depth - c.sea_level * s.rho_m / (s.rho_m - s.rho_w);
        }
    }
    M.airy.push_back(std::move(airy));
    M.flexural.push_back(std::move(flex));
  }
  return M;
}

double backstrip_map_value(const BackstripMaps& m, const std::vector<double>& g, const Domain& dom, double x, double y) {
  const double n = static_cast<double>(m.n);
  const double fx = std::clamp((x - dom.x0) / dom.sx * n - 0.5, 0.0, n - 1), fy = std::clamp((y - dom.y0) / dom.sy * n - 0.5, 0.0, n - 1);
  const std::size_t i = std::min(m.n - 2, static_cast<std::size_t>(fx)), j = std::min(m.n - 2, static_cast<std::size_t>(fy));
  const double u = fx - static_cast<double>(i), v = fy - static_cast<double>(j);
  auto at = [&](std::size_t a, std::size_t b) { return g[b * m.n + a]; };
  return (1 - u) * (1 - v) * at(i, j) + u * (1 - v) * at(i + 1, j) + (1 - u) * v * at(i, j + 1) + u * v * at(i + 1, j + 1);
}

}  // namespace stratonix
