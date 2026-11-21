// Stage 4A (STX-TN-14 section 1): discrete implicit field on a regular grid with free horizon
// values. Unknowns: node values and one value per free horizon. Rows: interface points
// (s(p) - c_k = 0), orientations (two tangent rows and a unit-magnitude row, scaled by the
// vertical spacing) and second differences along x, y and z (regularisation). Solved by
// conjugate-gradient least squares (CGLS) with Jacobi column scaling, coarse to fine. Every
// reduction is summed in fixed chunks in a fixed order, so results do not depend on the number
// of threads.
#include "stratonix/engines/fields.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"

namespace stratonix {
namespace {

constexpr std::size_t kChunk = 16384;  // reduction chunk (fixed: results independent of threads)

struct Lattice {
  std::size_t nx = 0, ny = 0, nz = 0;
  Vec3 o{0, 0, 0}, d{1, 1, 1};
  std::size_t size() const { return nx * ny * nz; }
  std::size_t idx(std::size_t i, std::size_t j, std::size_t k) const { return i + nx * (j + ny * k); }
  /// Cell and local coordinates; the cell is clamped to the lattice, t is not (extrapolation).
  void cell(const Vec3& p, std::size_t c[3], double t[3]) const {
    const double f[3] = {(p.x - o.x) / d.x, (p.y - o.y) / d.y, (p.z - o.z) / d.z};
    const std::size_t n[3] = {nx, ny, nz};
    for (int a = 0; a < 3; ++a) {
      double fl = std::floor(f[a]);
      long ci = static_cast<long>(fl);
      ci = std::clamp(ci, 0L, static_cast<long>(n[a]) - 2);
      c[a] = static_cast<std::size_t>(ci);
      t[a] = f[a] - static_cast<double>(ci);
    }
  }
  /// The 8 node indices and trilinear weights of p.
  void weights(const Vec3& p, std::size_t id[8], double w[8]) const {
    std::size_t c[3];
    double t[3];
    cell(p, c, t);
    int q = 0;
    for (int kk = 0; kk < 2; ++kk)
      for (int jj = 0; jj < 2; ++jj)
        for (int ii = 0; ii < 2; ++ii, ++q) {
          id[q] = idx(c[0] + ii, c[1] + jj, c[2] + kk);
          w[q] = (ii ? t[0] : 1 - t[0]) * (jj ? t[1] : 1 - t[1]) * (kk ? t[2] : 1 - t[2]);
        }
  }
  /// The 8 node indices and gradient weights (per unit length) of p.
  void grad_weights(const Vec3& p, std::size_t id[8], Vec3 g[8]) const {
    std::size_t c[3];
    double t[3];
    cell(p, c, t);
    int q = 0;
    for (int kk = 0; kk < 2; ++kk)
      for (int jj = 0; jj < 2; ++jj)
        for (int ii = 0; ii < 2; ++ii, ++q) {
          id[q] = idx(c[0] + ii, c[1] + jj, c[2] + kk);
          double wx = ii ? t[0] : 1 - t[0], wy = jj ? t[1] : 1 - t[1], wz = kk ? t[2] : 1 - t[2];
          double sx = (ii ? 1.0 : -1.0) / d.x, sy = (jj ? 1.0 : -1.0) / d.y, sz = (kk ? 1.0 : -1.0) / d.z;
          g[q] = {sx * wy * wz, wx * sy * wz, wx * wy * sz};
        }
  }
};

Lattice make_lattice(const BBox& box, double dxy, double dz) {
  Lattice L;
  auto n = [](double ext, double h) { return static_cast<std::size_t>(std::max(2.0, std::ceil(ext / h - 1e-9) + 1)); };
  L.nx = n(box.hi.x - box.lo.x, dxy);
  L.ny = n(box.hi.y - box.lo.y, dxy);
  L.nz = n(box.hi.z - box.lo.z, dz);
  L.o = box.lo;
  L.d = {(box.hi.x - box.lo.x) / static_cast<double>(L.nx - 1), (box.hi.y - box.lo.y) / static_cast<double>(L.ny - 1),
         (box.hi.z - box.lo.z) / static_cast<double>(L.nz - 1)};
  return L;
}

/// Sparse rows (CSR) of the data part, and their transpose (CSC) for deterministic A^T y.
struct Csr {
  std::vector<std::size_t> ptr{0};
  std::vector<std::uint32_t> col;
  std::vector<double> val, rhs;
  void add(std::size_t c, double v) { col.push_back(static_cast<std::uint32_t>(c)); val.push_back(v); }
  void end_row(double b) { rhs.push_back(b); ptr.push_back(col.size()); }
  std::size_t rows() const { return rhs.size(); }
};
struct Csc {
  std::vector<std::size_t> ptr;
  std::vector<std::uint32_t> row;
  std::vector<double> val;
};
Csc transpose(const Csr& a, std::size_t ncols) {
  Csc t;
  t.ptr.assign(ncols + 1, 0);
  for (auto c : a.col) ++t.ptr[c + 1];
  for (std::size_t c = 0; c < ncols; ++c) t.ptr[c + 1] += t.ptr[c];
  t.row.resize(a.col.size());
  t.val.resize(a.col.size());
  std::vector<std::size_t> pos(t.ptr.begin(), t.ptr.end() - 1);
  for (std::size_t r = 0; r + 1 < a.ptr.size(); ++r)
    for (std::size_t q = a.ptr[r]; q < a.ptr[r + 1]; ++q) {
      std::size_t& p = pos[a.col[q]];
      t.row[p] = static_cast<std::uint32_t>(r);
      t.val[p] = a.val[q];
      ++p;
    }
  return t;
}

double dot_det(const std::vector<double>& a, const std::vector<double>& b) {
  const std::size_t n = a.size(), nc = (n + kChunk - 1) / kChunk;
  std::vector<double> part(nc, 0.0);
  parallel_for(nc, 1, [&](std::size_t c0, std::size_t c1) {
    for (std::size_t c = c0; c < c1; ++c) {
      double s = 0;
      for (std::size_t i = c * kChunk, e = std::min(n, (c + 1) * kChunk); i < e; ++i) s += a[i] * b[i];
      part[c] = s;
    }
  });
  double s = 0;
  for (double p : part) s += p;
  return s;
}

/// One level's least-squares system: data rows (CSR) and three regularisation row blocks of N rows
/// each (second differences along x, y and z at interior nodes; zero rows elsewhere).
struct System {
  Lattice L;
  std::size_t nfree = 0;
  Csr A;
  Csc At;
  double wr = 0.2;
  std::size_t ncols() const { return L.size() + nfree; }
  std::size_t nrows() const { return A.rows() + 3 * L.size(); }

  /// y = A x (y sized nrows)
  void apply(const std::vector<double>& x, std::vector<double>& y) const {
    const std::size_t m = A.rows(), N = L.size();
    y.assign(nrows(), 0.0);
    parallel_for(m, 4096, [&](std::size_t r0, std::size_t r1) {
      for (std::size_t r = r0; r < r1; ++r) {
        double s = 0;
        for (std::size_t q = A.ptr[r]; q < A.ptr[r + 1]; ++q) s += A.val[q] * x[A.col[q]];
        y[r] = s;
      }
    });
    const std::size_t sx = 1, sy = L.nx, sz = L.nx * L.ny;
    parallel_for(N, 4096, [&](std::size_t j0, std::size_t j1) {
      for (std::size_t j = j0; j < j1; ++j) {
        std::size_t i = j % L.nx, jj = (j / L.nx) % L.ny, k = j / sz;
        if (i > 0 && i + 1 < L.nx) y[m + j] = wr * (x[j - sx] - 2 * x[j] + x[j + sx]);
        if (jj > 0 && jj + 1 < L.ny) y[m + N + j] = wr * (x[j - sy] - 2 * x[j] + x[j + sy]);
        if (k > 0 && k + 1 < L.nz) y[m + 2 * N + j] = wr * (x[j - sz] - 2 * x[j] + x[j + sz]);
      }
    });
  }
  /// x = A^T y (x sized ncols)
  void apply_t(const std::vector<double>& y, std::vector<double>& x) const {
    const std::size_t m = A.rows(), N = L.size(), nc = ncols();
    x.assign(nc, 0.0);
    const std::size_t sx = 1, sy = L.nx, sz = L.nx * L.ny;
    parallel_for(nc, 4096, [&](std::size_t c0, std::size_t c1) {
      for (std::size_t c = c0; c < c1; ++c) {
        double s = 0;
        for (std::size_t q = At.ptr[c]; q < At.ptr[c + 1]; ++q) s += At.val[q] * y[At.row[q]];
        if (c < N) {
          std::size_t i = c % L.nx, jj = (c / L.nx) % L.ny, k = c / sz;
          // column c appears in the x rows of c - 1, c, c + 1 (where those rows are interior)
          const double* yx = y.data() + m;
          const double* yy = y.data() + m + N;
          const double* yz = y.data() + m + 2 * N;
          s += wr * (-2 * yx[c] + (i > 0 ? yx[c - sx] : 0.0) + (i + 1 < L.nx ? yx[c + sx] : 0.0));
          s += wr * (-2 * yy[c] + (jj > 0 ? yy[c - sy] : 0.0) + (jj + 1 < L.ny ? yy[c + sy] : 0.0));
          s += wr * (-2 * yz[c] + (k > 0 ? yz[c - sz] : 0.0) + (k + 1 < L.nz ? yz[c + sz] : 0.0));
        }
        x[c] = s;
      }
    });
  }
  /// Squared column norms (for the Jacobi scaling).
  std::vector<double> column_norms2() const {
    const std::size_t N = L.size(), nc = ncols();
    std::vector<double> c2(nc, 0.0);
    for (std::size_t c = 0; c < nc; ++c)
      for (std::size_t q = At.ptr[c]; q < At.ptr[c + 1]; ++q) c2[c] += At.val[q] * At.val[q];
    const std::size_t n3[3] = {L.nx, L.ny, L.nz};
    for (std::size_t j = 0; j < N; ++j) {
      std::size_t ijk[3] = {j % L.nx, (j / L.nx) % L.ny, j / (L.nx * L.ny)};
      for (int a = 0; a < 3; ++a) {
        // rows of this axis that contain node j: the row at j (coefficient -2) if interior, and the
        // rows at its neighbours (coefficient 1) if those are interior
        double s = 0;
        if (ijk[a] > 0 && ijk[a] + 1 < n3[a]) s += 4;
        if (ijk[a] >= 2) s += 1;            // the neighbour below is interior
        if (ijk[a] + 2 < n3[a]) s += 1;     // the neighbour above is interior
        c2[j] += wr * wr * s;
      }
    }
    return c2;
  }
};

struct Rows {
  const std::vector<HorizonDatum>* data;
  const std::vector<RbfField::Gradient>* grads;
  const std::vector<long>* free_col;     // per horizon: -1 or index among free values
  const std::vector<double>* known;      // per horizon: known value (anchor or fixed) or NaN
  double ow = 1.0;
};

System assemble(const Lattice& L, const Rows& r, std::size_t nfree, double wr) {
  System s;
  s.L = L;
  s.nfree = nfree;
  s.wr = wr;
  const std::size_t N = L.size();
  s.A.col.reserve(r.data->size() * 9 + r.grads->size() * 24);
  s.A.val.reserve(s.A.col.capacity());
  std::size_t id[8];
  double w[8];
  for (auto& p : *r.data) {
    const double wt = p.sigma > 0 ? 1.0 / std::max(p.sigma, 0.1) : 1.0;
    L.weights(p.p, id, w);
    for (int q = 0; q < 8; ++q) s.A.add(id[q], wt * w[q]);
    const auto h = static_cast<std::size_t>(p.horizon);
    long fc = (*r.free_col)[h];
    if (fc >= 0) {
      s.A.add(N + static_cast<std::size_t>(fc), -wt);
      s.A.end_row(0.0);
    } else {
      s.A.end_row(wt * (*r.known)[h]);
    }
  }
  const double sc = L.d.z * r.ow;  // gradient rows: change of the field across one vertical cell
  Vec3 gw[8];
  for (auto& g : *r.grads) {
    Vec3 n = normalized(g.g);
    Vec3 t1 = normalized(cross(n, std::fabs(n.z) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0}));
    Vec3 t2 = cross(n, t1);
    L.grad_weights(g.p, id, gw);
    const Vec3 dirs[3] = {t1, t2, n};
    for (int a = 0; a < 3; ++a) {
      for (int q = 0; q < 8; ++q) s.A.add(id[q], sc * dot(dirs[a], gw[q]));
      s.A.end_row(a == 2 ? sc : 0.0);
    }
  }
  s.At = transpose(s.A, s.ncols());
  return s;
}

struct CglsResult { int iterations = 0; double residual = 0; bool converged = false; };

/// CGLS on A D z = b with D = diag(1 / ||A_j||), x = D z; x holds the start and the result.
CglsResult cgls(const System& S, std::vector<double>& x, int maxit, double tol) {
  const std::size_t nc = S.ncols();
  auto c2 = S.column_norms2();
  std::vector<double> D(nc);
  for (std::size_t c = 0; c < nc; ++c) D[c] = c2[c] > 0 ? 1.0 / std::sqrt(c2[c]) : 0.0;
  std::vector<double> z(nc);
  for (std::size_t c = 0; c < nc; ++c) z[c] = D[c] > 0 ? x[c] / D[c] : 0.0;
  // b
  std::vector<double> b(S.nrows(), 0.0);
  std::copy(S.A.rhs.begin(), S.A.rhs.end(), b.begin());
  std::vector<double> r, q, s, p(nc), tmp(nc);
  S.apply(x, q);
  r.resize(b.size());
  for (std::size_t i = 0; i < b.size(); ++i) r[i] = b[i] - q[i];
  std::vector<double> atb;
  S.apply_t(b, atb);
  for (std::size_t c = 0; c < nc; ++c) atb[c] *= D[c];
  const double nb = std::sqrt(dot_det(atb, atb));
  S.apply_t(r, s);
  for (std::size_t c = 0; c < nc; ++c) s[c] *= D[c];
  p = s;
  double gamma = dot_det(s, s);
  CglsResult out;
  if (nb == 0) { out.converged = true; return out; }
  out.residual = std::sqrt(gamma) / nb;
  if (out.residual < tol) { out.converged = true; return out; }
  for (int it = 0; it < maxit; ++it) {
    for (std::size_t c = 0; c < nc; ++c) tmp[c] = D[c] * p[c];
    S.apply(tmp, q);
    double qq = dot_det(q, q);
    if (!(qq > 0)) break;
    double alpha = gamma / qq;
    for (std::size_t c = 0; c < nc; ++c) z[c] += alpha * p[c];
    for (std::size_t i = 0; i < r.size(); ++i) r[i] -= alpha * q[i];
    S.apply_t(r, s);
    for (std::size_t c = 0; c < nc; ++c) s[c] *= D[c];
    double gnew = dot_det(s, s);
    out.iterations = it + 1;
    out.residual = std::sqrt(gnew) / nb;
    if (out.residual < tol) { out.converged = true; gamma = gnew; break; }
    double beta = gnew / gamma;
    gamma = gnew;
    for (std::size_t c = 0; c < nc; ++c) p[c] = s[c] + beta * p[c];
  }
  for (std::size_t c = 0; c < nc; ++c) x[c] = D[c] * z[c];
  return out;
}

double eval(const Lattice& L, const std::vector<double>& u, const Vec3& p) {
  std::size_t id[8];
  double w[8];
  L.weights(p, id, w);
  double s = 0;
  for (int q = 0; q < 8; ++q) s += w[q] * u[id[q]];
  return s;
}

}  // namespace

DiscreteField::DiscreteField(const std::vector<HorizonDatum>& data, const std::vector<RbfField::Gradient>& grads, int anchor,
                             double anchor_value, const BBox& box, const DiscreteOptions& o, const std::vector<double>& fixed) {
  if (!box.valid() || !(box.hi.z > box.lo.z)) throw Error(ErrorCode::InvalidArgument, "discrete field: empty box");
  int nh = anchor + 1;
  for (auto& p : data) nh = std::max(nh, p.horizon + 1);
  nh = std::max(nh, static_cast<int>(fixed.size()));
  const auto NH = static_cast<std::size_t>(nh);
  std::vector<std::size_t> count(NH, 0);
  std::vector<double> zsum(NH, 0.0);
  for (auto& p : data) {
    if (p.horizon < 0) throw Error(ErrorCode::InvalidArgument, "discrete field: a datum without a horizon");
    ++count[static_cast<std::size_t>(p.horizon)];
    zsum[static_cast<std::size_t>(p.horizon)] += p.p.z;
  }
  if (anchor < 0 || !count[static_cast<std::size_t>(anchor)])
    throw Error(ErrorCode::InvalidArgument, "discrete field: the anchor horizon has no interface data");
  std::vector<double> known(NH, kNaN);
  std::vector<long> free_col(NH, -1);
  known[static_cast<std::size_t>(anchor)] = anchor_value;
  for (std::size_t h = 0; h < NH; ++h) {
    if (static_cast<int>(h) == anchor) continue;
    if (h < fixed.size() && std::isfinite(fixed[h])) known[h] = fixed[h];
    else if (count[h]) free_col[h] = static_cast<long>(free_count_++);
  }
  // spacings and levels
  const double ext = std::max(box.hi.x - box.lo.x, box.hi.y - box.lo.y), vext = box.hi.z - box.lo.z;
  const double dxy = o.spacing > 0 ? o.spacing : ext / 80.0, dz = o.vertical_spacing > 0 ? o.vertical_spacing : vext / 80.0;
  Lattice fine = make_lattice(box, dxy, dz);
  int levels = o.levels;
  if (levels <= 0) {
    levels = 1;
    double f = 2;
    while (levels < 6) {
      Lattice c = make_lattice(box, dxy * f, dz * f);
      if (std::min({c.nx, c.ny, c.nz}) < 6) break;
      ++levels;
      f *= 2;
    }
  }
  Rows rows{&data, &grads, &free_col, &known, o.orientation_weight};
  std::vector<double> x;
  Lattice prev;
  std::vector<double> prev_u;
  for (int lv = levels - 1; lv >= 0; --lv) {
    const double f = std::pow(2.0, lv);
    Lattice L = lv == 0 ? fine : make_lattice(box, dxy * f, dz * f);
    System S = assemble(L, rows, free_count_, o.regularisation);
    const std::size_t N = L.size();
    std::vector<double> start(S.ncols(), 0.0);
    if (prev_u.empty()) {
      // field ~ elevation shifted to the anchor: s = z - z_anchor_mean + anchor_value
      const double shift = anchor_value - zsum[static_cast<std::size_t>(anchor)] / static_cast<double>(count[static_cast<std::size_t>(anchor)]);
      for (std::size_t k = 0; k < L.nz; ++k)
        for (std::size_t j = 0; j < L.ny; ++j)
          for (std::size_t i = 0; i < L.nx; ++i) start[L.idx(i, j, k)] = L.o.z + static_cast<double>(k) * L.d.z + shift;
      for (std::size_t h = 0; h < NH; ++h)
        if (free_col[h] >= 0) start[N + static_cast<std::size_t>(free_col[h])] = zsum[h] / static_cast<double>(count[h]) + shift;
    } else {
      parallel_for(L.nz, 1, [&](std::size_t k0, std::size_t k1) {
        for (std::size_t k = k0; k < k1; ++k)
          for (std::size_t j = 0; j < L.ny; ++j)
            for (std::size_t i = 0; i < L.nx; ++i)
              start[L.idx(i, j, k)] = eval(prev, prev_u, {L.o.x + static_cast<double>(i) * L.d.x, L.o.y + static_cast<double>(j) * L.d.y,
                                                          L.o.z + static_cast<double>(k) * L.d.z});
      });
      for (std::size_t c = 0; c < free_count_; ++c) start[N + c] = x[prev.size() + c];
    }
    auto res = cgls(S, start, o.max_iterations, o.tolerance);
    iterations_ += res.iterations;
    residual_ = res.residual;
    converged_ = res.converged;
    x = std::move(start);
    prev = L;
    prev_u.assign(x.begin(), x.begin() + static_cast<long>(N));
  }
  levels_ = levels;
  nx_ = fine.nx;
  ny_ = fine.ny;
  nz_ = fine.nz;
  o_ = fine.o;
  d_ = fine.d;
  u_ = prev_u;
  values_.assign(NH, kNaN);
  for (std::size_t h = 0; h < NH; ++h) {
    if (std::isfinite(known[h])) values_[h] = known[h];
    if (free_col[h] >= 0) values_[h] = x[u_.size() + static_cast<std::size_t>(free_col[h])];
  }
}

double DiscreteField::value(const Vec3& p) const {
  Lattice L;
  L.nx = nx_; L.ny = ny_; L.nz = nz_; L.o = o_; L.d = d_;
  return eval(L, u_, p);
}

Vec3 DiscreteField::gradient(const Vec3& p) const {
  Lattice L;
  L.nx = nx_; L.ny = ny_; L.nz = nz_; L.o = o_; L.d = d_;
  std::size_t id[8];
  Vec3 g[8];
  L.grad_weights(p, id, g);
  Vec3 s{0, 0, 0};
  for (int q = 0; q < 8; ++q) s += g[q] * u_[id[q]];
  return s;
}

std::string DiscreteField::summary() const {
  return "discrete grid " + std::to_string(nx_) + " x " + std::to_string(ny_) + " x " + std::to_string(nz_) + " (spacing " +
         fmt_fixed(d_.x, 1) + " m, vertical " + fmt_fixed(d_.z, 1) + " m), " + std::to_string(levels_) + " levels, " +
         std::to_string(iterations_) + " CGLS iterations, relative residual " + fmt_shortest(residual_) + (converged_ ? "" : " (not converged)");
}

std::string horizon_order_problem(const std::vector<double>& values, const std::vector<std::string>& names) {
  for (std::size_t k = 1; k < values.size(); ++k) {
    if (!std::isfinite(values[k]) || !std::isfinite(values[k - 1])) continue;
    if (!(values[k] < values[k - 1])) {
      auto nm = [&](std::size_t i) { return i < names.size() ? names[i] : "horizon " + std::to_string(i); };
      return "the fitted values of " + nm(k - 1) + " (" + fmt_fixed(values[k - 1], 2) + ") and " + nm(k) + " (" + fmt_fixed(values[k], 2) +
             ") are not in stratigraphic order";
    }
  }
  return {};
}

}  // namespace stratonix
