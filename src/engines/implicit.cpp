// Implicit (data-conditioned) structure builder: Hermite RBF stratigraphic field, surfaces in the
// restored frame, faults from the v1 fault network, iterated restoration of the data.
#include "stratonix/engines/implicit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <tuple>

#include "stratonix/core/error.hpp"
#include "stratonix/core/linalg.hpp"
#include "stratonix/core/roots.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/conditioning.hpp"
#include "stratonix/engines/fields.hpp"
#include "stratonix/engines/foldframe.hpp"
#include "stratonix/engines/wells.hpp"

namespace stratonix {

// ---- RBF field ---------------------------------------------------------------------------------
namespace {
inline double phi(const Vec3& d) {
  double r = norm(d);
  return r * r * r;
}
inline Vec3 dphi(const Vec3& d) {  // gradient of |d|^3 = 3 |d| d
  double r = norm(d);
  return d * (3 * r);
}
inline void hphi(const Vec3& d, double h[9]) {  // Hessian of |d|^3 = 3 (|d| I + d d^T / |d|)
  double r = norm(d);
  double dv[3] = {d.x, d.y, d.z};
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) h[3 * a + b] = 3 * ((a == b ? r : 0.0) + (r > 0 ? dv[a] * dv[b] / r : 0.0));
}
inline double comp(const Vec3& v, int k) { return k == 0 ? v.x : (k == 1 ? v.y : v.z); }
}  // namespace

namespace {
/// Assemble the Hermite RBF system (STX-TN-04; per-point nugget from STX-TN-11 section 5).
struct HermiteSystem {
  std::vector<Vec3> xi, yj, ej;
  std::vector<double> A, rhs;
  std::size_t n = 0, ni = 0, nd = 0;
};
HermiteSystem assemble(const std::vector<RbfField::Value>& values, const std::vector<RbfField::Gradient>& grads, double nugget,
                       const std::function<Vec3(const Vec3&)>& local, double scale, double stretch) {
  // Orientations constrain the direction of the gradient only: two tangent constraints per
  // orientation (t . grad s = 0), so the magnitude of the gradient is free and the interface
  // values set the spacing of the horizons (unit-gradient constraints conflict with values
  // wherever layer thickness varies).
  HermiteSystem S;
  S.ni = values.size();
  if (S.ni < 4) throw Error(ErrorCode::Build, "implicit field: at least four interface points are needed (have " + std::to_string(S.ni) + ")");
  for (auto& v : values) S.xi.push_back(local(v.p));
  for (auto& g : grads) {
    Vec3 n = normalized(g.g);
    Vec3 t1 = normalized(cross(n, std::fabs(n.z) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0}));
    Vec3 t2 = cross(n, t1);
    for (const Vec3& t : {t1, t2}) {
      S.yj.push_back(local(g.p));
      S.ej.push_back(normalized(Vec3{t.x, t.y, t.z * stretch}));  // tangent in stretched coordinates
    }
  }
  S.nd = S.yj.size();
  const std::size_t ni = S.ni, nd = S.nd, n = ni + nd + 4;
  S.n = n;
  S.A.assign(n * n, 0.0);
  S.rhs.assign(n, 0.0);
  auto at = [&](std::size_t r, std::size_t c) -> double& { return S.A[r * n + c]; };
  const std::size_t od = ni, op = ni + nd;
  for (std::size_t r = 0; r < ni; ++r) {
    for (std::size_t c = 0; c < ni; ++c) at(r, c) = phi(S.xi[r] - S.xi[c]);
    // the point's own uncertainty, else the field's nugget; values are scaled with the coordinates
    at(r, r) += (values[r].sigma >= 0 ? values[r].sigma : nugget) / scale;
    for (std::size_t j = 0; j < nd; ++j) at(r, od + j) = at(od + j, r) = -dot(S.ej[j], dphi(S.xi[r] - S.yj[j]));
    at(r, op) = at(op, r) = 1;
    for (int k = 0; k < 3; ++k) at(r, op + 1 + k) = at(op + 1 + k, r) = comp(S.xi[r], k);
    S.rhs[r] = values[r].v / scale;
  }
  for (std::size_t l = 0; l < nd; ++l) {
    for (std::size_t j = 0; j < nd; ++j) {
      double h[9];
      hphi(S.yj[l] - S.yj[j], h);
      const Vec3 &el = S.ej[l], &ej = S.ej[j];
      double v = 0;
      for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y) v += comp(el, x) * h[3 * x + y] * comp(ej, y);
      at(od + l, od + j) = -v;
    }
    for (int k = 0; k < 3; ++k) at(od + l, op + 1 + k) = at(op + 1 + k, od + l) = comp(S.ej[l], k);
    S.rhs[od + l] = 0;
  }
  return S;
}
}  // namespace

RbfField::RbfField(std::vector<Value> values, std::vector<Gradient> grads, double nugget, Vec3 centre, double scale, double vertical_stretch)
    : centre_(centre), scale_(scale), stretch_(vertical_stretch) {
  auto S = assemble(values, grads, nugget, [this](const Vec3& p) { return local(p); }, scale_, stretch_);
  const std::size_t n = S.n, ni = S.ni, od = ni, op = ni + S.nd;
  xi_ = std::move(S.xi);
  yj_ = std::move(S.yj);
  ej_ = std::move(S.ej);
  solver_ = optimised_kernels() && n >= 256 ? "blocked LU (optimised)" : "Gaussian elimination (reference)";
  solve_dense(S.A, S.rhs, n);
  auto& rhs = S.rhs;
  a_.assign(rhs.begin(), rhs.begin() + static_cast<long>(ni));
  beta_.assign(rhs.begin() + static_cast<long>(od), rhs.begin() + static_cast<long>(op));
  c0_ = rhs[op];
  c_ = {rhs[op + 1], rhs[op + 2], rhs[op + 3]};
}

double RbfField::value(const Vec3& xw) const {
  Vec3 x = local(xw);
  double s = c0_ + dot(c_, x);
  for (std::size_t i = 0; i < xi_.size(); ++i) s += a_[i] * phi(x - xi_[i]);
  for (std::size_t j = 0; j < yj_.size(); ++j) s -= beta_[j] * dot(ej_[j], dphi(x - yj_[j]));
  return s * scale_;
}

Vec3 RbfField::gradient(const Vec3& xw) const {
  Vec3 x = local(xw);
  Vec3 g = c_;
  for (std::size_t i = 0; i < xi_.size(); ++i) g = g + dphi(x - xi_[i]) * a_[i];
  for (std::size_t j = 0; j < yj_.size(); ++j) {
    double h[9];
    hphi(x - yj_[j], h);
    const Vec3& e = ej_[j];
    g = g - Vec3{h[0] * e.x + h[1] * e.y + h[2] * e.z, h[3] * e.x + h[4] * e.y + h[5] * e.z, h[6] * e.x + h[7] * e.y + h[8] * e.z} * beta_[j];
  }
  g.z *= stretch_;  // chain rule for the stretched vertical coordinate
  return g;
}

// ---- leave-one-out errors (Rippa 1999; STX-TN-11 section 5) ---------------------------------------
std::vector<double> rippa_loo(const std::vector<RbfField::Value>& values, const std::vector<RbfField::Gradient>& grads, double nugget,
                              Vec3 centre, double scale, double stretch) {
  auto local = [&](const Vec3& x) {
    Vec3 l = (x - centre) * (1.0 / scale);
    l.z *= stretch;
    return l;
  };
  auto S = assemble(values, grads, nugget, local, scale, stretch);
  const std::size_t n = S.n, m = 1 + S.ni;
  // right-hand sides: the data and the unit vectors of the value rows (for diag(A^-1))
  std::vector<double> B(n * m, 0.0);
  for (std::size_t r = 0; r < n; ++r) B[r * m] = S.rhs[r];
  for (std::size_t i = 0; i < S.ni; ++i) B[i * m + 1 + i] = 1.0;
  solve_dense_multi(S.A, B, n, m);
  std::vector<double> e(S.ni);
  for (std::size_t i = 0; i < S.ni; ++i) e[i] = B[i * m] / B[i * m + 1 + i] * scale;
  return e;
}

// ---- partition of unity (STX-TN-11 section 1) -----------------------------------------------------
namespace {
inline double wendland(double t) { return t >= 1 ? 0.0 : std::pow(1 - t, 4) * (4 * t + 1); }
inline double dwendland(double t) { return t >= 1 ? 0.0 : -20.0 * t * std::pow(1 - t, 3); }
/// True when the points span three dimensions (the linear drift is then determined).
bool spans_3d(const std::vector<Vec3>& p) {
  if (p.size() < 4) return false;
  Vec3 c{0, 0, 0};
  for (auto& x : p) c = c + x;
  c = c * (1.0 / static_cast<double>(p.size()));
  double cov[9] = {0};
  for (auto& x : p) {
    Vec3 d = x - c;
    double v[3] = {d.x, d.y, d.z};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) cov[3 * a + b] += v[a] * v[b];
  }
  double val[3], vec[9];
  symmetric_eigen3(cov, val, vec);
  return val[0] > 1e-9 * std::max(1.0, val[2]);
}
}  // namespace

PuField::PuField(std::vector<RbfField::Value> values, std::vector<RbfField::Gradient> grads, double nugget, Vec3 centre, double scale,
                 double stretch, double x0, double x1, double y0, double y1, const SolverOptions& o)
    : values_(std::move(values)), grads_(std::move(grads)), nugget_(nugget), scale_(scale), stretch_(stretch), centre_(centre), x0_(x0), y0_(y0) {
  const std::size_t N = values_.size();
  if (N < 4) throw Error(ErrorCode::Build, "implicit field: at least four interface points are needed (have " + std::to_string(N) + ")");
  const double sx = std::max(x1 - x0, 1.0), sy = std::max(y1 - y0, 1.0);
  const std::size_t target = std::max<std::size_t>(o.points_per_patch, 1);
  const double npatch = std::max(1.0, std::round(static_cast<double>(N) / static_cast<double>(target)));
  h_ = std::sqrt(sx * sy / npatch);
  nx_ = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(sx / h_ - 1e-9)));
  ny_ = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(sy / h_ - 1e-9)));
  const double hx = sx / static_cast<double>(nx_), hy = sy / static_cast<double>(ny_);
  const double r0 = o.overlap * std::hypot(hx, hy) / 2;
  const double rmax = 2 * std::hypot(sx, sy);
  for (std::size_t j = 0; j < ny_; ++j)
    for (std::size_t i = 0; i < nx_; ++i) {
      Patch p;
      p.cx = x0 + (static_cast<double>(i) + 0.5) * hx;
      p.cy = y0 + (static_cast<double>(j) + 0.5) * hy;
      p.r = nx_ * ny_ == 1 ? rmax : r0;  // one patch holds every datum
      patches_.push_back(p);
    }
  bool distinct_values = false;  // the data hold at least two values
  for (std::size_t i = 1; i < N && !distinct_values; ++i) distinct_values = values_[i].v != values_[0].v;
  // members: every datum within the radius; grow a patch with too few or coplanar data
  parallel_for(patches_.size(), 1, [&](std::size_t a, std::size_t b) {
    for (std::size_t k = a; k < b; ++k) {
      Patch& p = patches_[k];
      for (;;) {
        p.vi.clear();
        p.gi.clear();
        std::vector<Vec3> pos;
        for (std::size_t i = 0; i < N; ++i)
          if (std::hypot(values_[i].p.x - p.cx, values_[i].p.y - p.cy) < p.r) {
            p.vi.push_back(i);
            pos.push_back(values_[i].p);
          }
        for (std::size_t i = 0; i < grads_.size(); ++i)
          if (std::hypot(grads_[i].p.x - p.cx, grads_[i].p.y - p.cy) < p.r) p.gi.push_back(i);
        // at least 30 data spanning three dimensions, of at least two horizons (values): with one
        // value the magnitude of the gradient, and so the spacing of the horizons, is undetermined
        bool two_values = !distinct_values;
        for (std::size_t m = 1; m < p.vi.size() && !two_values; ++m) two_values = values_[p.vi[m]].v != values_[p.vi[0]].v;
        if ((p.vi.size() >= std::min<std::size_t>(30, N) && spans_3d(pos) && two_values) || p.r >= rmax) break;
        p.r *= 1.25;
      }
    }
  });
  // local fields (each a dense Hermite RBF with the global centre and scale)
  parallel_for(patches_.size(), 1, [&](std::size_t a, std::size_t b) {
    for (std::size_t k = a; k < b; ++k) {
      Patch& p = patches_[k];
      std::vector<RbfField::Value> v;
      std::vector<RbfField::Gradient> g;
      for (auto i : p.vi) v.push_back(values_[i]);
      for (auto i : p.gi) g.push_back(grads_[i]);
      p.f = std::make_shared<RbfField>(std::move(v), std::move(g), nugget_, centre_, scale_, stretch_);
    }
  });
}

std::size_t PuField::largest_patch() const {
  std::size_t m = 0;
  for (auto& p : patches_) m = std::max(m, p.f ? p.f->unknowns() : 0);
  return m;
}

double PuField::value(const Vec3& x) const {
  if (patches_.size() == 1) return patches_[0].f->value(x);
  double num = 0, den = 0;
  const Patch* nearest = nullptr;
  double best = 1e300;
  for (auto& p : patches_) {
    double r = std::hypot(x.x - p.cx, x.y - p.cy);
    if (r < best) { best = r; nearest = &p; }
    if (r >= p.r) continue;
    double w = wendland(r / p.r);
    num += w * p.f->value(x);
    den += w;
  }
  if (den > 0) return num / den;
  return nearest->f->value(x);  // outside every patch (beyond the covered area): the nearest patch
}

Vec3 PuField::gradient(const Vec3& x) const {
  if (patches_.size() == 1) return patches_[0].f->gradient(x);
  double W = 0, S = 0;
  Vec3 G{0, 0, 0}, GW{0, 0, 0};
  const Patch* nearest = nullptr;
  double best = 1e300;
  for (auto& p : patches_) {
    double dx = x.x - p.cx, dy = x.y - p.cy, r = std::hypot(dx, dy);
    if (r < best) { best = r; nearest = &p; }
    if (r >= p.r) continue;
    double t = r / p.r, w = wendland(t), s = p.f->value(x);
    Vec3 gw = r > 0 ? Vec3{dx / r, dy / r, 0} * (dwendland(t) / p.r) : Vec3{0, 0, 0};
    W += w;
    S += w * s;
    G = G + p.f->gradient(x) * w + gw * s;
    GW = GW + gw;
  }
  if (W > 0) return (G - GW * (S / W)) * (1.0 / W);
  return nearest->f->gradient(x);
}

std::vector<double> PuField::leave_one_out() const {
  const std::size_t N = values_.size();
  std::vector<double> num(N, 0.0), den(N, 0.0);
  std::vector<std::vector<double>> per(patches_.size());
  parallel_for(patches_.size(), 1, [&](std::size_t a, std::size_t b) {
    for (std::size_t k = a; k < b; ++k) {
      const Patch& p = patches_[k];
      std::vector<RbfField::Value> v;
      std::vector<RbfField::Gradient> g;
      for (auto i : p.vi) v.push_back(values_[i]);
      for (auto i : p.gi) g.push_back(grads_[i]);
      per[k] = rippa_loo(v, g, nugget_, centre_, scale_, stretch_);
    }
  });
  for (std::size_t k = 0; k < patches_.size(); ++k) {
    const Patch& p = patches_[k];
    for (std::size_t m = 0; m < p.vi.size(); ++m) {
      const auto& x = values_[p.vi[m]].p;
      double w = patches_.size() == 1 ? 1.0 : wendland(std::hypot(x.x - p.cx, x.y - p.cy) / p.r);
      num[p.vi[m]] += w * per[k][m];
      den[p.vi[m]] += w;
    }
  }
  std::vector<double> e(N, kNaN);
  for (std::size_t i = 0; i < N; ++i)
    if (den[i] > 0) e[i] = num[i] / den[i];
  return e;
}

// ---- plane fit -------------------------------------------------------------------------------------
FittedPlane fit_plane(const std::vector<Polyline>& sticks) {
  FittedPlane fp;
  Vec3 c{0, 0, 0};
  for (auto& s : sticks)
    for (auto& p : s) { c = c + p; ++fp.points; }
  if (fp.points < 3) throw Error(ErrorCode::Build, "fault sticks: at least three points are needed to fit a plane");
  c = c * (1.0 / static_cast<double>(fp.points));
  double cov[9] = {0};
  for (auto& s : sticks)
    for (auto& p : s) {
      Vec3 d = p - c;
      double v[3] = {d.x, d.y, d.z};
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) cov[3 * a + b] += v[a] * v[b];
    }
  double val[3], vec[9];
  symmetric_eigen3(cov, val, vec);
  Vec3 n{vec[0], vec[3], vec[6]};  // smallest eigenvalue
  if (n.z < 0) n = n * -1.0;
  fp.centroid = c;
  fp.normal = normalized(n);
  fp.dip = std::acos(std::clamp(fp.normal.z, -1.0, 1.0));
  double dipdir = std::atan2(fp.normal.x, fp.normal.y);  // up-normal leans towards the dip direction
  fp.strike = std::fmod(dipdir - kPi / 2 + 4 * kPi, 2 * kPi);
  double ss = 0;
  for (auto& s : sticks)
    for (auto& p : s) ss += std::pow(dot(p - c, fp.normal), 2);
  fp.rms = std::sqrt(ss / static_cast<double>(fp.points));
  if (fp.dip < 1e-3) throw Error(ErrorCode::Build, "fault sticks define a horizontal plane; a fault must dip");
  return fp;
}

// ---- model -------------------------------------------------------------------------------------
ImplicitModel::ImplicitModel(std::shared_ptr<ParametricModel> inner, ModelInfo info, std::vector<DataMisfit> misfit)
    : inner_(std::move(inner)), info_(std::move(info)), misfit_(std::move(misfit)) {}

const void* ImplicitModel::capability(std::string_view name) const {
  if (name == "data_misfit") return &misfit_;
  return inner_->capability(name);
}

// ---- builder -------------------------------------------------------------------------------------
namespace {
desc::Param make_param(const std::string& path, std::vector<double> values, const std::string& unit, Dimension dim, double to_si,
                       desc::Evidence ev, const std::string& source) {
  desc::Param p;
  p.present = true;
  p.path = path;
  p.values = std::move(values);
  if (p.values.size() > 1) p.shape = {p.values.size()};
  p.unit = unit;
  p.dim = dim;
  p.to_si = to_si;
  p.evidence = ev;
  p.source = source;
  p.text = "";
  for (std::size_t i = 0; i < p.values.size(); ++i) p.text += (i ? ", " : "") + fmt_shortest(p.values[i]);
  if (p.values.size() > 1) p.text = "[" + p.text + "]";
  p.text += " " + unit + " | " + desc::evidence_name(ev) + " | " + source;
  return p;
}

desc::Evidence evidence_of(const std::string& s) {
  static const std::map<std::string, desc::Evidence> m = {{"measured", desc::Evidence::Measured}, {"interpreted", desc::Evidence::Interpreted},
                                                         {"calibrated", desc::Evidence::Calibrated}, {"upscaled", desc::Evidence::Upscaled},
                                                         {"literature_derived", desc::Evidence::LiteratureDerived}, {"assumed", desc::Evidence::Assumed},
                                                         {"user_instruction", desc::Evidence::UserInstruction}};
  auto it = m.find(s);
  return it == m.end() ? desc::Evidence::Assumed : it->second;
}

struct Datum {
  Vec3 present, restored;
  Vec3 normal{0, 0, 0};           ///< present-day normal (data)
  Vec3 normal_restored{0, 0, 0};  ///< normal carried through the restoration
  bool oriented = false;
  int horizon = -1;
  double sigma = -1;              ///< stated uncertainty (m; degrees for orientations); < 0 = none
  std::string input, evidence, kind;
};

/// The field of one stratigraphic series: dense or partition of unity (STX-TN-11 sections 1, 2).
struct SeriesField {
  std::unique_ptr<RbfField> dense;
  std::unique_ptr<PuField> pu;
  std::unique_ptr<DiscreteField> disc;  ///< stage 4A (STX-TN-14 section 1)
  double value(const Vec3& x) const { return dense ? dense->value(x) : (pu ? pu->value(x) : disc->value(x)); }
};

/// Surfaces of the field (s = value_h) on a grid in the restored frame, by scanning each column
/// from the top and refining the first crossing (s decreases downwards).
template <class F>
std::vector<Grid2D> extract_surfaces(const F& f, const std::vector<double>& vals, const GridSpec& g, double ztop, double zbot,
                                     std::size_t& clamped) {
  const std::size_t nh = vals.size(), K = 96;
  std::vector<Grid2D> out(nh, Grid2D(g));
  std::vector<std::size_t> clamp_rows(g.ny, 0);
  parallel_for(g.ny, 1, [&](std::size_t j0, std::size_t j1) {
    std::vector<double> sz(K + 1);
    for (std::size_t j = j0; j < j1; ++j)
      for (std::size_t i = 0; i < g.nx; ++i) {
        double x = g.x(i), y = g.y(j);
        for (std::size_t k = 0; k <= K; ++k) sz[k] = f.value({x, y, ztop + (zbot - ztop) * static_cast<double>(k) / K});
        for (std::size_t h = 0; h < nh; ++h) {
          double v = vals[h], z = kNaN;
          for (std::size_t k = 0; k < K && std::isnan(z); ++k) {
            double fa = sz[k] - v, fb = sz[k + 1] - v;
            if (fa >= 0 && fb < 0) {
              double za = ztop + (zbot - ztop) * static_cast<double>(k) / K, zb = ztop + (zbot - ztop) * static_cast<double>(k + 1) / K;
              z = illinois([&](double zz) { return f.value({x, y, zz}) - v; }, za, zb, fa, fb, 1e-6);
            }
          }
          if (std::isnan(z)) {  // no crossing inside the domain: hold the surface at the nearer boundary
            z = sz[0] < v ? ztop : zbot;
            ++clamp_rows[j];
          }
          out[h].at(i, j) = z;
        }
      }
  });
  clamped = 0;
  for (auto c : clamp_rows) clamped += c;
  // Horizons of one series must not cross: a deeper horizon is never above a shallower one.
  for (std::size_t h = 1; h < nh; ++h)
    for (std::size_t k = 0; k < out[h].v.size(); ++k) out[h].v[k] = std::min(out[h].v[k], out[h - 1].v[k]);
  return out;
}

/// 2D thin-plate spline (Duchon 1977) with linear drift: thickness of an isochore horizon from its
/// data. More than 3000 data are first averaged on a lattice (stated in the build log).
struct Tps2D {
  std::vector<double> x, y, w;
  double c0 = 0, cx = 0, cy = 0, ox = 0, oy = 0, sc = 1;
  static double k(double r2) { return r2 > 0 ? 0.5 * r2 * std::log(r2) : 0.0; }  // r^2 log r
  Tps2D(std::vector<std::array<double, 3>> p, double lattice) {
    if (p.size() > 3000) {  // average on a lattice
      std::map<std::pair<long, long>, std::array<double, 4>> b;
      for (auto& q : p) {
        auto& a = b[{std::lround(q[0] / lattice), std::lround(q[1] / lattice)}];
        a[0] += q[0]; a[1] += q[1]; a[2] += q[2]; a[3] += 1;
      }
      p.clear();
      for (auto& [key, a] : b) p.push_back({a[0] / a[3], a[1] / a[3], a[2] / a[3]});
    }
    const std::size_t n = p.size();
    if (n < 3) throw Error(ErrorCode::Build, "isochore horizon: fewer than three thickness data");
    for (auto& q : p) { ox += q[0]; oy += q[1]; }
    ox /= static_cast<double>(n);
    oy /= static_cast<double>(n);
    for (auto& q : p) sc = std::max(sc, std::max(std::fabs(q[0] - ox), std::fabs(q[1] - oy)));
    for (auto& q : p) { x.push_back((q[0] - ox) / sc); y.push_back((q[1] - oy) / sc); }
    const std::size_t m = n + 3;
    std::vector<double> A(m * m, 0.0), r(m, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) A[i * m + j] = k((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
      A[i * m + n] = A[n * m + i] = 1;
      A[i * m + n + 1] = A[(n + 1) * m + i] = x[i];
      A[i * m + n + 2] = A[(n + 2) * m + i] = y[i];
      r[i] = p[i][2];
    }
    solve_dense(A, r, m);
    w.assign(r.begin(), r.begin() + static_cast<long>(n));
    c0 = r[n]; cx = r[n + 1]; cy = r[n + 2];
  }
  double operator()(double xw, double yw) const {
    double u = (xw - ox) / sc, v = (yw - oy) / sc, s = c0 + cx * u + cy * v;
    for (std::size_t i = 0; i < w.size(); ++i) s += w[i] * k((u - x[i]) * (u - x[i]) + (v - y[i]) * (v - y[i]));
    return s;
  }
};

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
}  // namespace

ModelPtr build_implicit_model(const desc::Description& d, const BuildInputs& in, BuildLog& log) {
  const std::size_t nh = d.horizons.size();
  std::map<std::string, const desc::InputDesc*> inputs;
  for (auto& x : d.inputs) inputs[x.id] = &x;
  // 1. Data.
  std::vector<Datum> data;
  std::set<std::string> ignored_labels;
  std::vector<std::string> provenance;
  for (auto& e : d.implicit.data) {
    auto it = in.data.find(e.input);
    if (it == in.data.end()) throw Error(ErrorCode::Build, "implicit data input '" + e.input + "' was not loaded");
    // stage 4C: a geological map's data carry the evidence of the map input
    const std::string src_id = e.input.rfind("geological_map:", 0) == 0 ? (e.input == "geological_map:orientations" ? d.geological_map.measurements.input
                                                                                                                       : d.geological_map.units.input)
                                                                       : e.input;
    auto inp_it = inputs.find(src_id);
    if (inp_it == inputs.end()) throw Error(ErrorCode::Build, "input '" + src_id + "' is not listed under inputs");
    const auto* inp = inp_it->second;
    std::size_t used = 0, with_sigma = 0;
    const bool orient = e.kind == "orientations";
    // input-level uncertainty (m; orientations in degrees)
    const double sig_in = e.uncertainty.present ? (orient ? rad2deg(e.uncertainty.si()) : e.uncertainty.si()) : -1.0;
    for (auto& dp : it->second) {
      int h = d.horizon_index(dp.label);
      if (h < 0) { ignored_labels.insert(dp.label); continue; }
      Datum dt;
      dt.present = dt.restored = dp.p;
      dt.horizon = h;
      dt.oriented = dp.oriented;
      dt.normal = dt.normal_restored = dp.n;
      dt.input = e.input;
      dt.evidence = inp->evidence.empty() ? "assumed" : inp->evidence;
      dt.kind = dp.oriented ? "orientation" : "point";
      dt.sigma = dp.sigma >= 0 ? dp.sigma : sig_in;
      with_sigma += dt.sigma >= 0;
      data.push_back(dt);
      ++used;
    }
    provenance.push_back("Data " + e.input + " (" + inp->path + "): " + std::to_string(used) + " " + (orient ? "orientations" : "points") +
                         ", evidence " + (inp->evidence.empty() ? std::string("not stated (treated as assumed)") : inp->evidence) +
                         (with_sigma ? ", " + std::to_string(with_sigma) + " with a stated uncertainty" : std::string()) + ", source: " + inp->source);
  }
  if (!ignored_labels.empty()) {
    std::string s;
    for (auto& l : ignored_labels) s += (s.empty() ? "" : ", ") + (l.empty() ? std::string("(no label)") : l);
    log.note("implicit: data labelled " + s + " are not horizons and were not used for the stratigraphic field");
  }
  if (d.implicit.use_well_tops) {
    auto wells = make_well_paths(d);
    std::size_t nt = 0;
    for (std::size_t w = 0; w < d.wells.size(); ++w) {
      auto& wd = d.wells[w];
      if (wd.role != "control") continue;
      for (auto& [hn, p] : wd.tops) {
        int h = d.horizon_index(hn);
        if (h < 0) continue;
        double z = -p.si();  // tops are depths below z = 0 (as in the parametric engine)
        Datum dt;
        dt.present = dt.restored = wells[w].vertical() ? Vec3{wells[w].x(), wells[w].y(), z} : wells[w].at_elevation(z);
        dt.horizon = h;
        dt.input = "well:" + wd.name;
        dt.evidence = desc::evidence_name(p.evidence);
        dt.kind = "well_top";
        dt.sigma = wd.top_uncertainty.present ? wd.top_uncertainty.si() : -1.0;
        data.push_back(dt);
        ++nt;
      }
    }
    if (nt) provenance.push_back("Data: " + std::to_string(nt) + " control-well tops (evidence as labelled in the description)");
  }
  // Interface data of different horizons at one position contradict each other (exact
  // interpolation is impossible); they are refused rather than altered.
  {
    std::map<std::tuple<long, long, long>, std::size_t> seen;
    for (std::size_t i = 0; i < data.size(); ++i) {
      auto& dt = data[i];
      if (dt.oriented || dt.sigma > 0) continue;
      auto key = std::make_tuple(std::lround(dt.present.x * 1000), std::lround(dt.present.y * 1000), std::lround(dt.present.z * 1000));
      auto [it, fresh] = seen.emplace(key, i);
      if (!fresh && data[it->second].horizon != dt.horizon)
        throw Error(ErrorCode::Build, "implicit: data of " + d.horizons[static_cast<std::size_t>(data[it->second].horizon)].name + " (" +
                                          data[it->second].input + ") and " + d.horizons[static_cast<std::size_t>(dt.horizon)].name + " (" + dt.input +
                                          ") coincide at (" + fmt_fixed(dt.present.x, 2) + ", " + fmt_fixed(dt.present.y, 2) + ", " +
                                          fmt_fixed(dt.present.z, 2) + "): give them an uncertainty or remove the duplicate pick");
    }
  }
  // 2. Series: one field per stratigraphic series (STX-TN-11 section 2); without series, one field.
  std::vector<std::vector<int>> groups;
  std::vector<std::string> group_names;
  {
    std::map<int, std::size_t> of;
    for (std::size_t h = 0; h < nh; ++h) {
      int sr = d.series.empty() ? -1 : d.series_of(static_cast<int>(h));
      if (!of.count(sr)) {
        of[sr] = groups.size();
        groups.emplace_back();
        group_names.push_back(sr < 0 ? std::string("all horizons") : d.series[static_cast<std::size_t>(sr)].name);
      }
      groups[of[sr]].push_back(static_cast<int>(h));
    }
  }
  std::vector<int> group_of(nh, 0);
  for (std::size_t g = 0; g < groups.size(); ++g)
    for (int h : groups[g]) group_of[static_cast<std::size_t>(h)] = static_cast<int>(g);
  // Isochore horizons (STX-TN-11 section 2): not in their series' field; built from the nearest
  // shallower field horizon of the series minus a thickness interpolated from their own data.
  std::vector<bool> iso(nh, false);
  std::vector<int> iso_ref(nh, -1);
  for (auto& hn : d.implicit.isochore_horizons) iso[static_cast<std::size_t>(d.horizon_index(hn))] = true;
  for (std::size_t h = 0; h < nh; ++h)
    if (iso[h])
      for (int a = static_cast<int>(h) - 1; a >= 0; --a)
        if (group_of[static_cast<std::size_t>(a)] == group_of[h] && !iso[static_cast<std::size_t>(a)]) { iso_ref[h] = a; break; }
  // 3. Horizon values: elevation at the reference well when it has the top, else the mean of the data.
  std::vector<double> vals(nh, kNaN), sum(nh, 0.0);
  std::vector<std::size_t> cnt(nh, 0);
  for (auto& dt : data)
    if (!dt.oriented) { sum[static_cast<std::size_t>(dt.horizon)] += dt.present.z; ++cnt[static_cast<std::size_t>(dt.horizon)]; }
  const desc::WellDesc* rw = d.ref_well.empty() ? nullptr : d.well(d.ref_well);
  std::vector<bool> from_well(nh, false);
  for (std::size_t h = 0; h < nh; ++h) {
    if (rw)
      for (auto& [hn, p] : rw->tops)
        if (hn == d.horizons[h].name) { vals[h] = -p.si(); from_well[h] = true; }
  }
  // A horizon without a reference-well top: the value of the nearest horizon above it in its series
  // with a value, minus their mean thickness where both have data (points within 150 m laterally);
  // else the mean elevation of its data.
  auto mean_thickness = [&](std::size_t ha, std::size_t hb, double& t) {
    const double cell = 150.0;
    std::map<std::pair<long, long>, std::vector<Vec3>> ga;
    for (auto& dt : data)
      if (!dt.oriented && dt.horizon == static_cast<int>(ha))
        ga[{std::lround(dt.present.x / cell), std::lround(dt.present.y / cell)}].push_back(dt.present);
    double s = 0;
    std::size_t n = 0;
    for (auto& dt : data) {
      if (dt.oriented || dt.horizon != static_cast<int>(hb)) continue;
      long ci = std::lround(dt.present.x / cell), cj = std::lround(dt.present.y / cell);
      double best = cell, za = kNaN;
      for (long a = -1; a <= 1; ++a)
        for (long b = -1; b <= 1; ++b) {
          auto it = ga.find({ci + a, cj + b});
          if (it == ga.end()) continue;
          for (auto& q : it->second) {
            double r = std::hypot(q.x - dt.present.x, q.y - dt.present.y);
            if (r < best) { best = r; za = q.z; }
          }
        }
      if (!std::isnan(za)) { s += za - dt.present.z; ++n; }
    }
    if (n < 3) return false;
    t = s / static_cast<double>(n);
    return true;
  };
  for (std::size_t h = 0; h < nh; ++h) {
    if (!std::isnan(vals[h])) continue;
    // stage 4C: a stated stratigraphic column (thickness_above) below a horizon with a value
    if (h > 0 && d.horizons[h].thickness_above.present && !std::isnan(vals[h - 1])) {
      vals[h] = vals[h - 1] - d.horizons[h].thickness_above.si();
      continue;
    }
    int sr = d.series.empty() ? -1 : d.series_of(static_cast<int>(h));
    for (int a = static_cast<int>(h) - 1; a >= 0 && std::isnan(vals[h]); --a) {
      if ((d.series.empty() ? -1 : d.series_of(a)) != sr || std::isnan(vals[static_cast<std::size_t>(a)])) continue;
      double t;
      if (mean_thickness(static_cast<std::size_t>(a), h, t) && t > 0) vals[h] = vals[static_cast<std::size_t>(a)] - t;
      break;
    }
    if (std::isnan(vals[h]) && cnt[h]) vals[h] = sum[h] / static_cast<double>(cnt[h]);
    if (std::isnan(vals[h])) throw Error(ErrorCode::Build, "implicit: horizon " + d.horizons[h].name + " has no data point and no reference-well top");
  }
  const bool free_values = d.implicit.horizon_values == "free";  // stage 4A: fitted on the discrete grid
  for (auto g : groups) {
    if (free_values) break;
    g.erase(std::remove_if(g.begin(), g.end(), [&](int h) { return iso[static_cast<std::size_t>(h)]; }), g.end());
    for (std::size_t k = 1; k < g.size(); ++k) {
      auto h0 = static_cast<std::size_t>(g[k - 1]), h1 = static_cast<std::size_t>(g[k]);
      if (!(vals[h1] < vals[h0]))
        throw Error(ErrorCode::Build, "implicit: the field values of " + d.horizons[h0].name + " and " + d.horizons[h1].name +
                                          " are not in stratigraphic order (check the data labels)");
    }
  }
  // 4. Faults given by sticks: plane or dip profile fitted to the sticks; throws from data.
  desc::Description d2 = d;
  for (auto& f : d2.faults) {
    if (f.geom_type != "sticks") continue;
    auto it = in.lines.find(f.geom_input);
    if (it == in.lines.end()) throw Error(ErrorCode::Build, "fault " + f.name + ": sticks input '" + f.geom_input + "' was not loaded");
    const auto* inp = inputs.at(f.geom_input);
    auto ev = evidence_of(inp->evidence);
    std::string b = "faults." + f.name + ".geometry";
    if (f.sticks_fit == "profile") {
      // listric or ramp-flat: dip profile from the sticks (STX-TN-11 section 3.2), ramp_flat in the curved mechanism
      StickProfile pr = profile_from_sticks(it->second);
      std::vector<double> rows;
      const double u0 = pr.u.front();
      for (std::size_t k = 0; k < pr.depth.size(); k += 2) {  // 100 m nodes
        rows.push_back(pr.u[k] - u0);
        rows.push_back(pr.depth[k]);
      }
      if (rows.size() < 4) throw Error(ErrorCode::Build, "fault " + f.name + ": the sticks give fewer than two profile nodes");
      std::string src = "dip profile fitted to the sticks of input " + f.geom_input + " (" + std::to_string(rows.size() / 2) + " nodes, RMS " +
                        fmt_fixed(pr.rms, 3) + " m)";
      f.geom_type = "ramp_flat";
      f.reference_point = make_param(b + ".reference_point", {pr.top_point.x, pr.top_point.y}, "m", Dimension::Length, 1, ev, src);
      f.reference_depth = make_param(b + ".reference_depth", {pr.depth.front()}, "m", Dimension::Length, 1, ev, src);
      f.strike = make_param(b + ".strike", {rad2deg(pr.strike)}, "deg", Dimension::Angle, kPi / 180, ev, src);
      f.profile = make_param(b + ".profile", rows, "m", Dimension::Length, 1, ev, src);
      f.profile.shape = {rows.size() / 2, 2};
      f.upper.rule = "free";
      f.lower.rule = "free";
      log.value("implicit.fault." + f.name + ".fitted_profile", "strike " + fmt_fixed(rad2deg(pr.strike), 2) + " deg, dip " +
                                                               fmt_fixed(rad2deg(pr.dip.front()), 2) + " deg at " + fmt_fixed(pr.depth.front(), 0) +
                                                               " m to " + fmt_fixed(rad2deg(pr.dip.back()), 2) + " deg at " +
                                                               fmt_fixed(pr.depth.back(), 0) + " m, RMS " + fmt_fixed(pr.rms, 3) + " m");
      provenance.push_back("Fault " + f.name + ": " + src + ", evidence " + desc::evidence_name(ev));
      continue;
    }
    FittedPlane fp = fit_plane(it->second);
    std::string src = "plane fitted to " + std::to_string(fp.points) + " stick points of input " + f.geom_input + " (RMS " + fmt_fixed(fp.rms, 3) + " m)";
    f.geom_type = "planar";
    f.reference_point = make_param(b + ".reference_point", {fp.centroid.x, fp.centroid.y}, "m", Dimension::Length, 1, ev, src);
    f.reference_depth = make_param(b + ".reference_depth", {-fp.centroid.z}, "m", Dimension::Length, 1, ev, src);
    f.strike = make_param(b + ".strike", {rad2deg(fp.strike)}, "deg", Dimension::Angle, kPi / 180, ev, src);
    f.dip = make_param(b + ".dip", {rad2deg(fp.dip)}, "deg", Dimension::Angle, kPi / 180, ev, src);
    log.value("implicit.fault." + f.name + ".fitted_plane", "strike " + fmt_fixed(rad2deg(fp.strike), 2) + " deg, dip " + fmt_fixed(rad2deg(fp.dip), 2) +
                                                            " deg, centre (" + fmt_fixed(fp.centroid.x, 1) + ", " + fmt_fixed(fp.centroid.y, 1) + ", " +
                                                            fmt_fixed(fp.centroid.z, 1) + "), RMS " + fmt_fixed(fp.rms, 3) + " m");
    provenance.push_back("Fault " + f.name + ": " + src + ", evidence " + desc::evidence_name(ev));
  }
  // Throws, length and tips from horizon cut-offs (STX-TN-11 section 3.1).
  {
    std::vector<HorizonPoint> hp;
    std::set<std::string> data_inputs;
    for (auto& dt : data)
      if (!dt.oriented) { hp.push_back({dt.present, dt.horizon}); data_inputs.insert(dt.input); }
    std::string inputs_txt;
    for (auto& x : data_inputs) inputs_txt += (inputs_txt.empty() ? "" : ", ") + x;
    for (auto& f : d2.faults) {
      if (!f.throw_estimate) continue;
      FittedPlane pl;
      pl.centroid = {f.reference_point.si(0), f.reference_point.si(1), -f.reference_depth.si()};
      pl.strike = f.strike.si();
      pl.dip = f.dip.si();
      pl.normal = normalized(Vec3{std::sin(pl.dip) * std::cos(pl.strike), -std::sin(pl.dip) * std::sin(pl.strike), std::cos(pl.dip)});
      ThrowOptions o;
      o.band = f.estimate_band.si();
      o.station_spacing = f.estimate_spacing.si();
      o.profile = f.profile_strike.type;
      o.taper = f.profile_strike.taper.present ? f.profile_strike.taper.si() : 0.25;
      auto est = estimate_throws(pl, hp, nh, o);
      std::vector<double> tm, len, cen;
      int hmin = -1, hmax = -1;
      std::size_t stations = 0;
      double rms = 0;
      std::string per_horizon;
      for (auto& e : est) {
        const std::string hname = d.horizons[static_cast<std::size_t>(e.horizon)].name;
        std::string key = "implicit.fault." + f.name + ".throw_" + hname;
        if (!e.ok) {
          log.value(key, "not estimated: " + e.reason);
          per_horizon += "; " + hname + ": " + e.reason;
          continue;
        }
        const bool offset = e.max_throw > std::max(1.0, 2 * e.rms);
        std::string txt = fmt_fixed(e.max_throw, 2) + " m over " + fmt_fixed(e.length, 0) + " m, centre " + fmt_fixed(e.centre_s, 0) + " m along strike, RMS " +
                          fmt_fixed(e.rms, 2) + " m, " + std::to_string(e.stations.size()) + " stations" + (offset ? "" : " (not offset)");
        log.value(key, txt);
        per_horizon += "; " + hname + ": " + txt;
        if (!offset) continue;
        tm.push_back(e.max_throw);
        len.push_back(e.length);
        cen.push_back(e.centre_s);
        stations += e.stations.size();
        rms = std::max(rms, e.rms);
        if (hmin < 0) hmin = e.horizon;
        hmax = e.horizon;
      }
      if (tm.empty())
        throw Error(ErrorCode::Build, "fault " + f.name + ": no horizon is offset by more than its station scatter; the throw cannot be estimated from the data" +
                                          per_horizon);
      double T = median(tm), L = median(len), C = median(cen);
      std::string src = "estimated from horizon cut-offs of " + inputs_txt + " (" + std::to_string(stations) + " stations on " +
                        std::to_string(tm.size()) + " horizons, RMS " + fmt_fixed(rms, 2) + " m)";
      // stage 4A (STX-TN-14 section 3.1): one displacement field over the fault surface
      DisplacementField dfield;
      if (f.estimate_model == "displacement_field") {
        DisplacementOptions dop;
        dop.shape_strike = f.profile_strike.type.empty() ? "flat_topped" : f.profile_strike.type;
        dop.taper = f.profile_strike.taper.present ? f.profile_strike.taper.si() : 0.25;
        dop.shape_dip = "sine";  // finite throw gradient at the tips (the stage 2A admissibility check)
        auto dobs = displacement_observations(pl, est, hp);
        // stations within the band of another fault see that fault's offset too: excluded
        {
          const Vec3 sd{std::sin(pl.strike), std::cos(pl.strike), 0};
          const Vec3 dn = normalized(Vec3{std::cos(pl.strike) * std::cos(pl.dip), -std::sin(pl.strike) * std::cos(pl.dip), -std::sin(pl.dip)});
          std::size_t before = dobs.size();
          dobs.erase(std::remove_if(dobs.begin(), dobs.end(),
                                    [&](const DisplacementObservation& q) {
                                      Vec3 X = pl.centroid + sd * q.s + dn * q.d;
                                      for (auto& g : d2.faults) {
                                        if (&g == &f || g.geom_type != "planar" || !g.reference_point.present || !g.strike.present || !g.dip.present)
                                          continue;
                                        Vec3 cg{g.reference_point.si(0), g.reference_point.si(1), -g.reference_depth.si()};
                                        double st = g.strike.si(), dp = g.dip.si();
                                        Vec3 ng{std::sin(dp) * std::cos(st), -std::sin(dp) * std::sin(st), std::cos(dp)};
                                        Vec3 sg{std::sin(st), std::cos(st), 0};
                                        if (g.length.present && std::fabs(dot(X - cg, sg)) > 0.5 * g.length.si() + o.band) continue;
                                        if (std::fabs(dot(X - cg, ng)) / std::max(1e-6, std::sin(dp)) <= o.band) return true;
                                      }
                                      return false;
                                    }),
                     dobs.end());
          if (dobs.size() < before)
            log.value("implicit.fault." + f.name + ".displacement_field_excluded",
                      std::to_string(before - dobs.size()) + " of " + std::to_string(before) +
                          " observations within the band of another fault (they see its offset too)");
        }
        dfield = fit_displacement_field(dobs, dop);
        if (!dfield.ok) throw Error(ErrorCode::Build, "fault " + f.name + ": the displacement field cannot be fitted: " + dfield.reason);
        T = dfield.max_throw;
        L = dfield.length;
        C = dfield.s0;
        src = "displacement field fitted to " + std::to_string(dfield.observations) + " cut-off separations of " + inputs_txt + " (RMS " +
              fmt_fixed(dfield.rms, 2) + " m)";
        log.value("implicit.fault." + f.name + ".displacement_field",
                  "maximum " + fmt_fixed(T, 2) + " m at " + fmt_fixed(C, 0) + " m along strike" +
                      (std::isfinite(dfield.height) ? " and " + fmt_fixed(dfield.d0, 0) + " m down dip; height " + fmt_fixed(dfield.height, 0) + " m"
                                                    : "; height beyond the data (tips free down dip)") +
                      ", length " + fmt_fixed(L, 0) + " m, RMS " + fmt_fixed(dfield.rms, 2) + " m, " + std::to_string(dfield.observations) +
                      " observations (calibrated)");
      }
      std::string b = "faults." + f.name;
      f.throw_max = make_param(b + ".throw.max", {T}, "m", Dimension::Length, 1, desc::Evidence::Calibrated, src);
      f.length = make_param(b + ".extent.length", {L}, "m", Dimension::Length, 1, desc::Evidence::Calibrated, src);
      const Vec3 sdir{std::sin(pl.strike), std::cos(pl.strike), 0};
      f.reference_point = make_param(b + ".geometry.reference_point", {pl.centroid.x + sdir.x * C, pl.centroid.y + sdir.y * C}, "m", Dimension::Length, 1,
                                     desc::Evidence::Calibrated, src);
      // vertical tips: in the unit between the outermost offset horizons and their unoffset neighbours
      auto unit_between = [&](int ha, int hb) -> std::string {
        for (auto& u : d.units)
          if (u.top == d.horizons[static_cast<std::size_t>(ha)].name && u.base == d.horizons[static_cast<std::size_t>(hb)].name) return u.name;
        return "";
      };
      auto set_tip = [&](desc::TipDesc& t, int ha, int hb, bool at_edge) {
        if (!t.rule.empty()) return;  // stated tips are kept
        std::string un = at_edge ? "" : unit_between(ha, hb);
        if (un.empty()) { t.rule = "free"; return; }
        t.rule = "within_unit";
        t.unit = un;
        t.fraction = make_param(b + ".extent.tip.fraction", {0.5}, "fraction", Dimension::Dimensionless, 1, desc::Evidence::Calibrated,
                                "midway between the unoffset and the offset horizon (" + src + ")");
      };
      if (dfield.ok && std::isfinite(dfield.height) && kPi * T / dfield.height >= 0.8) {
        log.note("fault " + f.name + ": the displacement field's height (" + fmt_fixed(dfield.height, 0) +
                 " m) would give a throw gradient of " + fmt_fixed(kPi * T / dfield.height, 2) +
                 " near the tips; the tips are set from the offset horizons instead");
        dfield.height = std::numeric_limits<double>::infinity();
      }
      if (dfield.ok && std::isfinite(dfield.height)) {
        // tips where the displacement field ends down dip (depths), unless stated
        const double sind = std::sin(pl.dip);
        auto depth_tip = [&](desc::TipDesc& t, double dd, const std::string& which) {
          if (!t.rule.empty()) return;
          t.rule = "depth";
          t.depth = make_param(b + ".extent." + which + ".depth", {-pl.centroid.z + dd * sind}, "m", Dimension::Length, 1, desc::Evidence::Calibrated, src);
        };
        depth_tip(f.upper, dfield.d0 - 0.5 * dfield.height, "upper_tip");
        depth_tip(f.lower, dfield.d0 + 0.5 * dfield.height, "lower_tip");
        f.profile_dip.type = "sine";
      }
      set_tip(f.upper, hmin - 1, hmin, hmin == 0);
      set_tip(f.lower, hmax, hmax + 1, hmax + 1 >= static_cast<int>(nh));
      f.throw_estimate = false;
      log.value("implicit.fault." + f.name + ".estimated_throw", fmt_fixed(T, 2) + " m, length " + fmt_fixed(L, 0) + " m, centre " + fmt_fixed(C, 0) +
                                                                 " m along strike from the plane's reference point (calibrated: " + src + ")");
      provenance.push_back("Fault " + f.name + ": throw, length and tips " + src);
    }
  }
  // Fault relations inferred from data and rules (stage 4A, STX-TN-14 section 3.2): logged; applied to
  // faults without a stated relation when implicit.fault_relations is infer.
  {
    std::vector<FaultDataSet> fds;
    for (std::size_t k = 0; k < d.faults.size(); ++k) {
      if (d.faults[k].geom_type != "sticks") continue;
      auto it = in.lines.find(d.faults[k].geom_input);
      if (it == in.lines.end()) continue;
      FaultDataSet fs;
      fs.name = d.faults[k].name;
      fs.plane = fit_plane(it->second);
      for (auto& st : it->second)
        for (auto& p : st) fs.points.push_back(p);
      if (d2.faults[k].length.present) fs.tip_extent = 0.5 * d2.faults[k].length.si();
      fds.push_back(std::move(fs));
    }
    const double rd = d.implicit.relation_distance.present ? d.implicit.relation_distance.si() : 150.0;  // Stratonix default
    for (auto& r : infer_fault_relations(fds, rd)) {
      log.value("implicit.fault_relation." + r.a + "." + r.b, r.kind + ": " + r.rule);
      if (d.implicit.fault_relations != "infer" || r.kind == "crosses") continue;
      int ka = d2.fault_index(r.a), kb = d2.fault_index(r.b);
      auto& fa = d2.faults[static_cast<std::size_t>(ka)];
      if (!fa.abuts_fault.empty()) {
        log.note("implicit: the stated relation of " + r.a + " is kept (inferred: " + r.kind + " " + r.b + ")");
        continue;
      }
      const auto& fb = d2.faults[static_cast<std::size_t>(kb)];
      if (fa.age.present && fb.age.present && fb.age.si() > fa.age.si()) {
        log.note("implicit: " + r.a + " " + r.kind + " " + r.b + " was inferred but " + r.b + " is younger: not applied");
        continue;
      }
      fa.abuts_fault = r.b;
      fa.abuts_side = r.side;
      fa.abuts_transfer = r.kind == "splay";
      provenance.push_back("Fault " + r.a + ": " + r.kind + " " + r.b + " inferred from data (" + r.rule + ")");
    }
  }
  // 5. Iterate: fields in the restored frame -> surfaces -> parametric model with faults -> restore data.
  const Domain dom = d.domain();
  double sp = d.implicit.surface_spacing.present ? d.implicit.surface_spacing.si() : 2 * d.grid_spacing.si();  // Stratonix default: 2 x grid
  GridSpec g = GridSpec::covering(dom, sp);
  // Maximum number of restoration passes; stop earlier when no datum moves by more than 1 cm.
  int iters = d.implicit.iterations.present ? static_cast<int>(d.implicit.iterations.si()) : 6;  // Stratonix default
  const double converged = 0.01;  // m, Stratonix default
  if (d.faults.empty()) iters = 1;
  double nugget = d.implicit.nugget.present ? d.implicit.nugget.si() : 0.0;
  double stretch = d.implicit.vertical_stretch.present ? d.implicit.vertical_stretch.si() : 5.0;  // Stratonix default
  Vec3 centre{dom.x0 + dom.sx / 2, dom.y0 + dom.sy / 2, 0.5 * (dom.z_top + dom.z_base)};
  double scale = 0.5 * std::max({dom.sx, dom.sy, dom.z_top - dom.z_base});
  SolverOptions so;
  so.type = d.implicit.solver_type;
  if (d.implicit.points_per_patch.present) so.points_per_patch = static_cast<std::size_t>(d.implicit.points_per_patch.si());
  if (d.implicit.overlap.present) so.overlap = d.implicit.overlap.si();
  for (std::size_t h = 0; h < nh; ++h) {
    d2.horizons[h].surface_input = "__implicit_" + d.horizons[h].name;
    d2.horizons[h].surface_format = "implicit";
    d2.horizons[h].surface_z_is = "elevation";
    d2.horizons[h].iso_from.clear();
    d2.horizons[h].iso_input.clear();
  }
  // series relations stay (truncation, onlap); design-only series settings do not apply to data
  for (auto& sr : d2.series) {
    sr.tilt_dip = desc::Param{};
    sr.tilt_azimuth = desc::Param{};
    sr.clino_profile.clear();
  }
  d2.primitives.clear();
  d2.implicit = desc::ImplicitDesc{};
  // The implicit field already interpolates the control-well tops; the parametric engine's
  // correction field must not adjust the surfaces further (misfits are reported, not hidden).
  for (auto& w : d2.wells)
    if (w.role == "control") w.role = "check";
  std::shared_ptr<ParametricModel> inner;
  double moved = 0;
  std::vector<std::string> solver_txt(groups.size());
  std::size_t unknowns = 0, clamped = 0;
  int passes = 1;
  std::vector<SeriesField> fields(groups.size());
  for (int it = 0; it < iters; ++it) {
    BuildInputs in2 = in;
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
      std::vector<RbfField::Value> vs;
      std::vector<RbfField::Gradient> gs;
      // isochore horizons: each datum takes the value of its reference horizon minus its thickness
      // below the nearest reference pick (within two surface spacings), so it constrains the field's
      // gradient without forcing a constant value on a unit that thins (STX-TN-11 section 2)
      std::map<int, std::map<std::pair<long, long>, std::vector<Vec3>>> refpicks;
      const double cell = 2 * sp;
      for (auto& dt : data)
        if (!dt.oriented && group_of[static_cast<std::size_t>(dt.horizon)] == static_cast<int>(gi))
          for (std::size_t h = 0; h < nh; ++h)
            if (iso[h] && iso_ref[h] == dt.horizon)
              refpicks[dt.horizon][{std::lround(dt.restored.x / cell), std::lround(dt.restored.y / cell)}].push_back(dt.restored);
      std::size_t iso_used = 0;
      for (auto& dt : data) {
        if (group_of[static_cast<std::size_t>(dt.horizon)] != static_cast<int>(gi)) continue;
        if (dt.oriented) {
          gs.push_back({dt.restored, dt.normal_restored});
          continue;
        }
        if (!iso[static_cast<std::size_t>(dt.horizon)]) {
          vs.push_back({dt.restored, vals[static_cast<std::size_t>(dt.horizon)], dt.sigma});
          continue;
        }
        const int ref = iso_ref[static_cast<std::size_t>(dt.horizon)];
        auto& picks = refpicks[ref];
        long ci = std::lround(dt.restored.x / cell), cj = std::lround(dt.restored.y / cell);
        double best = cell, zr = kNaN;
        for (long a2 = -1; a2 <= 1; ++a2)
          for (long b2 = -1; b2 <= 1; ++b2) {
            auto itp = picks.find({ci + a2, cj + b2});
            if (itp == picks.end()) continue;
            for (auto& q : itp->second) {
              double r = std::hypot(q.x - dt.restored.x, q.y - dt.restored.y);
              if (r < best) { best = r; zr = q.z; }
            }
          }
        const double t = zr - dt.restored.z;
        if (std::isnan(zr) || !(t > 0.5)) continue;  // no reference pick near, or (near) zero thickness
        vs.push_back({dt.restored, vals[static_cast<std::size_t>(ref)] - t, dt.sigma});
        ++iso_used;
      }
      if (it == 0 && iso_used) log.value("implicit.isochore_data_in_field." + group_names[gi], std::to_string(iso_used));
      // fold frames (stage 4A, STX-TN-14 section 2): predicted bedding normals at lattice points
      // without measured orientations constrain the series' field
      for (auto& fd : d.implicit.folds) {
        if (!(fd.series.empty() ? groups.size() == 1 || group_names[gi] == "all horizons" : group_names[gi] == fd.series)) continue;
        std::vector<OrientationDatum> bed, fol;
        BBox bb;
        double zs = 0;
        std::size_t nz = 0;
        for (auto& dt : data) {
          if (group_of[static_cast<std::size_t>(dt.horizon)] != static_cast<int>(gi)) continue;
          bb.add(dt.restored);
          if (dt.oriented) bed.push_back({dt.restored, dt.normal_restored});
          else { zs += dt.restored.z; ++nz; }
        }
        if (!fd.foliation.empty()) {
          auto itf = in.data.find(fd.foliation);
          if (itf != in.data.end())
            for (auto& dp : itf->second)
              if (dp.oriented) fol.push_back({dp.p, dp.n});
        }
        FoldFrameOptions fo;
        if (fd.wavelength.present) fo.wavelength = fd.wavelength.si();
        if (fd.order.present) fo.max_order = static_cast<int>(fd.order.si());
        FoldFrame ff = fit_fold_frame(bed, fol, fo);
        if (!ff.ok) throw Error(ErrorCode::Build, "implicit: fold frame " + fd.name + ": " + ff.reason);
        const double step = fd.constraint_spacing.present ? fd.constraint_spacing.si() : std::max(2 * sp, ff.wavelength / 16);
        // constraint levels: the mean elevation of each horizon of the series with data
        std::map<int, std::pair<double, std::size_t>> hz;
        for (auto& dt : data)
          if (!dt.oriented && group_of[static_cast<std::size_t>(dt.horizon)] == static_cast<int>(gi)) {
            hz[dt.horizon].first += dt.restored.z;
            ++hz[dt.horizon].second;
          }
        std::vector<double> levels;
        for (auto& [hh, v] : hz) levels.push_back(v.first / static_cast<double>(v.second));
        if (levels.empty()) levels.push_back(nz ? zs / static_cast<double>(nz) : 0.5 * (bb.lo.z + bb.hi.z));
        std::map<std::pair<long, long>, bool> measured;
        for (auto& o2 : bed) measured[{std::lround(o2.p.x / step), std::lround(o2.p.y / step)}] = true;
        std::size_t added = 0;
        for (double y = bb.lo.y; y <= bb.hi.y + 1e-9; y += step)
          for (double x = bb.lo.x; x <= bb.hi.x + 1e-9; x += step) {
            if (measured.count({std::lround(x / step), std::lround(y / step)})) continue;
            for (double zl : levels) {
              Vec3 q{x, y, zl};
              gs.push_back({q, ff.bedding_normal(q)});
              ++added;
            }
          }
        if (it == 0)
          log.value("implicit.fold." + fd.name,
                    "axis trend " + fmt_fixed(rad2deg(ff.trend), 1) + " deg, plunge " + fmt_fixed(rad2deg(ff.plunge), 1) + " deg; interlimb " +
                        fmt_fixed(rad2deg(ff.interlimb), 1) + " deg; wavelength " + fmt_fixed(ff.wavelength, 0) + " m; Fourier order " +
                        std::to_string(ff.order) + "; RMS " + fmt_fixed(rad2deg(ff.rms), 2) + " deg on " + std::to_string(ff.data) +
                        " bedding orientations; cylindricity " + fmt_fixed(ff.cylindricity, 3) + "; " + std::to_string(added) +
                        " constraints at " + fmt_fixed(step, 0) + " m (calibrated)" + (ff.reason.empty() ? "" : "; " + ff.reason));
      }
      if (vs.size() < 4)
        throw Error(ErrorCode::Build, "implicit: series " + group_names[gi] + " has " + std::to_string(vs.size()) +
                                          " interface data; at least four are needed for its field");
      const bool use_disc = free_values || so.type == "discrete";
      const bool use_pu = !use_disc && (so.type == "partition_of_unity" || (so.type == "auto" && vs.size() > so.dense_limit));
      SeriesField sf;
      if (use_disc) {
        // stage 4A (STX-TN-14 section 1): the series' horizons indexed locally; the anchor is the
        // shallowest horizon with data; free values are fitted, fixed values are known rows
        std::vector<int> loc(nh, -1), hs;
        for (int h : groups[gi])
          if (!iso[static_cast<std::size_t>(h)]) { loc[static_cast<std::size_t>(h)] = static_cast<int>(hs.size()); hs.push_back(h); }
        std::vector<HorizonDatum> hd;
        std::vector<std::size_t> nper(hs.size(), 0);
        for (auto& dt : data) {
          if (dt.oriented || group_of[static_cast<std::size_t>(dt.horizon)] != static_cast<int>(gi)) continue;
          int lk = loc[static_cast<std::size_t>(dt.horizon)];
          if (lk < 0) continue;
          hd.push_back({dt.restored, lk, dt.sigma});
          ++nper[static_cast<std::size_t>(lk)];
        }
        int anchor = -1;
        for (std::size_t k = 0; k < hs.size() && anchor < 0; ++k)
          if (nper[k]) anchor = static_cast<int>(k);
        if (anchor < 0) throw Error(ErrorCode::Build, "implicit: series " + group_names[gi] + " has no interface data for its discrete field");
        std::vector<double> fixed(hs.size(), kNaN);
        if (!free_values)
          for (std::size_t k = 0; k < hs.size(); ++k) fixed[k] = vals[static_cast<std::size_t>(hs[k])];
        DiscreteOptions dopt;
        if (d.implicit.discrete_spacing.present) dopt.spacing = d.implicit.discrete_spacing.si();
        if (d.implicit.discrete_vertical_spacing.present) dopt.vertical_spacing = d.implicit.discrete_vertical_spacing.si();
        if (d.implicit.discrete_regularisation.present) dopt.regularisation = d.implicit.discrete_regularisation.si();
        if (d.implicit.discrete_max_iterations.present) dopt.max_iterations = static_cast<int>(d.implicit.discrete_max_iterations.si());
        if (d.implicit.discrete_orientation_weight.present) dopt.orientation_weight = d.implicit.discrete_orientation_weight.si();
        const double vext = dom.z_top - dom.z_base;
        BBox box;
        box.add({dom.x0, dom.y0, dom.z_base - 0.05 * vext});
        box.add({dom.x0 + dom.sx, dom.y0 + dom.sy, dom.z_top + 0.05 * vext});
        sf.disc = std::make_unique<DiscreteField>(hd, gs, anchor, vals[static_cast<std::size_t>(hs[static_cast<std::size_t>(anchor)])], box, dopt, fixed);
        if (free_values) {
          auto& fv = sf.disc->horizon_values();
          std::vector<std::string> names;
          for (std::size_t k = 0; k < hs.size(); ++k) {
            names.push_back(d.horizons[static_cast<std::size_t>(hs[k])].name);
            if (!std::isfinite(fv[k]))
              throw Error(ErrorCode::Build, "implicit: free horizon values: " + names.back() + " has no interface data in series " + group_names[gi]);
            vals[static_cast<std::size_t>(hs[k])] = fv[k];
          }
          auto problem = horizon_order_problem(fv, names);
          if (!problem.empty()) throw Error(ErrorCode::Build, "implicit (rule F34): " + problem + " (check the data labels)");
        }
        solver_txt[gi] = group_names[gi] + ": " + sf.disc->summary() + ", " + std::to_string(hd.size()) + " interface data" +
                         (free_values ? ", free horizon values" : "");
        unknowns = std::max(unknowns, sf.disc->unknowns());
      } else if (use_pu) {
        sf.pu = std::make_unique<PuField>(vs, gs, nugget, centre, scale, stretch, dom.x0, dom.x0 + dom.sx, dom.y0, dom.y0 + dom.sy, so);
        solver_txt[gi] = group_names[gi] + ": partition of unity, " + std::to_string(sf.pu->patches()) + " patches, largest system " +
                         std::to_string(sf.pu->largest_patch()) + " unknowns, " + std::to_string(vs.size()) + " interface data";
        unknowns = std::max(unknowns, sf.pu->largest_patch());
      } else {
        sf.dense = std::make_unique<RbfField>(vs, gs, nugget, centre, scale, stretch);
        solver_txt[gi] = group_names[gi] + ": dense " + sf.dense->solver() + ", " + std::to_string(sf.dense->unknowns()) + " unknowns";
        unknowns = std::max(unknowns, sf.dense->unknowns());
      }
      std::vector<double> gv;
      std::vector<int> fh;  // the series' horizons from the field
      for (int h : groups[gi])
        if (!iso[static_cast<std::size_t>(h)]) { gv.push_back(vals[static_cast<std::size_t>(h)]); fh.push_back(h); }
      std::size_t cl = 0;
      auto surfs = extract_surfaces(sf, gv, g, dom.z_top, dom.z_base, cl);
      clamped += cl;
      if (sf.disc) {
        // stage 4A (STX-TN-14 section 1): the discrete field is a least-squares fit; its surfaces are
        // corrected by a thin-plate spline of their residuals at the horizon's own data, so the data
        // are honoured as by the stage 3A interpolants, and the series order is kept
        for (std::size_t k = 0; k < fh.size(); ++k) {
          std::vector<std::array<double, 3>> res;
          double ss = 0, mx = 0;
          for (auto& dt : data)
            if (!dt.oriented && dt.horizon == fh[k]) {
              double r = dt.restored.z - surfs[k].sample(dt.restored.x, dt.restored.y);
              if (!std::isfinite(r)) continue;
              res.push_back({dt.restored.x, dt.restored.y, r});
              ss += r * r;
              mx = std::max(mx, std::fabs(r));
            }
          if (res.size() < 3) continue;
          Tps2D tps(res, 2 * sp);
          for (std::size_t j = 0; j < g.ny; ++j)
            for (std::size_t i = 0; i < g.nx; ++i) surfs[k].at(i, j) += tps(g.x(i), g.y(j));
          log.value("implicit.discrete_correction." + d.horizons[static_cast<std::size_t>(fh[k])].name,
                    "residuals of the discrete field at " + std::to_string(res.size()) + " data: rms " +
                        fmt_fixed(std::sqrt(ss / static_cast<double>(res.size())), 2) + " m, max " + fmt_fixed(mx, 2) +
                        " m; corrected by a thin-plate spline");
        }
        for (std::size_t k = 1; k < fh.size(); ++k)
          for (std::size_t q = 0; q < surfs[k].v.size(); ++q) surfs[k].v[q] = std::min(surfs[k].v[q], surfs[k - 1].v[q]);
      }
      std::map<int, Grid2D> made;
      for (std::size_t k = 0; k < fh.size(); ++k) made.emplace(fh[k], std::move(surfs[k]));
      // isochore horizons: reference surface minus the interpolated thickness, never above it
      for (int h : groups[gi]) {
        if (!iso[static_cast<std::size_t>(h)]) continue;
        const Grid2D& ref = made.at(iso_ref[static_cast<std::size_t>(h)]);
        std::vector<std::array<double, 3>> th;
        for (auto& dt : data)
          if (!dt.oriented && dt.horizon == h) th.push_back({dt.restored.x, dt.restored.y, ref.sample(dt.restored.x, dt.restored.y) - dt.restored.z});
        th.erase(std::remove_if(th.begin(), th.end(), [](const std::array<double, 3>& q) { return std::isnan(q[2]); }), th.end());
        Tps2D tps(th, 2 * sp);
        Grid2D sg(g);
        for (std::size_t j = 0; j < g.ny; ++j)
          for (std::size_t i = 0; i < g.nx; ++i) sg.at(i, j) = ref.at(i, j) - std::max(0.0, tps(g.x(i), g.y(j)));
        // keep the series order with deeper field horizons
        for (auto& [hh, gg] : made)
          if (hh > h)
            for (std::size_t q = 0; q < gg.v.size(); ++q) gg.v[q] = std::min(gg.v[q], sg.v[q]);
        made.emplace(h, std::move(sg));
        if (it == 0) log.value("implicit.isochore." + d.horizons[static_cast<std::size_t>(h)].name,
                               "thickness below " + d.horizons[static_cast<std::size_t>(iso_ref[static_cast<std::size_t>(h)])].name + " from " +
                                   std::to_string(th.size()) + " data (thin-plate spline, clamped at zero)");
      }
      for (auto& [hh, gg] : made) in2.grids[d2.horizons[static_cast<std::size_t>(hh)].surface_input] = std::move(gg);
      fields[gi] = std::move(sf);
    }
    BuildLog l2;
    inner = std::make_shared<ParametricModel>(d2, in2, l2);
    if (it + 1 == iters || (it > 0 && moved < converged)) {
      for (auto& kv : l2.values) log.value(kv.first, kv.second);
      for (auto& n : l2.notes) log.note(n);
      passes = it + 1;
      break;
    }
    clamped = 0;
    // Restore the data; an orientation is carried by restoring two points 5 m away in its plane.
    std::vector<Vec3> pres;
    for (auto& dt : data) {
      pres.push_back(dt.present);
      if (dt.oriented) {
        Vec3 t1 = normalized(cross(dt.normal, std::fabs(dt.normal.z) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0}));
        Vec3 t2 = cross(dt.normal, t1);
        pres.push_back(dt.present + t1 * 5.0);
        pres.push_back(dt.present + t2 * 5.0);
      }
    }
    std::vector<Vec3> rest(pres.size());
    inner->restore(pres, rest);
    moved = 0;
    std::size_t k = 0;
    for (auto& dt : data) {
      moved = std::max(moved, norm(rest[k] - dt.restored));
      dt.restored = rest[k];
      if (dt.oriented) {
        Vec3 n = normalized(cross(rest[k + 1] - rest[k], rest[k + 2] - rest[k]));
        if (n.z < 0) n = n * -1.0;
        dt.normal_restored = n;
        k += 3;
      } else {
        ++k;
      }
    }
    log.value("implicit.iteration_" + std::to_string(it + 1) + ".max_data_shift_m", fmt_fixed(moved, 4));
  }
  // 6. Leave-one-out errors (data influence, STX-TN-11 section 5), in the restored frame of the last pass.
  std::vector<double> loo(data.size(), kNaN);
  if (d.implicit.data_influence) {
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
      std::vector<RbfField::Value> vs;
      std::vector<RbfField::Gradient> gs;
      std::vector<std::size_t> idx;
      for (std::size_t i = 0; i < data.size(); ++i) {
        auto& dt = data[i];
        if (group_of[static_cast<std::size_t>(dt.horizon)] != static_cast<int>(gi) || iso[static_cast<std::size_t>(dt.horizon)]) continue;
        if (dt.oriented) gs.push_back({dt.restored, dt.normal_restored});
        else { vs.push_back({dt.restored, vals[static_cast<std::size_t>(dt.horizon)], dt.sigma}); idx.push_back(i); }
      }
      if (fields[gi].disc) {  // a least-squares fit: leave-one-out errors are not computed (STX-TN-14 section 1)
        log.note("implicit: data influence (leave-one-out errors) is not computed for the discrete field of series " + group_names[gi]);
        continue;
      }
      auto e = fields[gi].pu ? fields[gi].pu->leave_one_out() : rippa_loo(vs, gs, nugget, centre, scale, stretch);
      for (std::size_t k = 0; k < idx.size() && k < e.size(); ++k) loo[idx[k]] = e[k];
    }
  }
  // 7. Misfit at every data point, through the final model (present-day positions).
  std::vector<DataMisfit> mis;
  std::vector<Vec3> pres;
  for (auto& dt : data) pres.push_back(dt.present);
  EvalResult er;
  inner->evaluate(pres, kEvalScalars, er);
  double maxp = 0, sump = 0, maxo = 0;
  std::size_t np = 0, no = 0;
  for (std::size_t i = 0; i < data.size(); ++i) {
    auto& dt = data[i];
    DataMisfit m{dt.input, dt.evidence, d.horizons[static_cast<std::size_t>(dt.horizon)].name, dt.kind, dt.present, 0};
    if (!dt.oriented) {
      m.misfit = er.s(static_cast<std::size_t>(dt.horizon), i);
      maxp = std::max(maxp, std::fabs(m.misfit));
      sump += std::fabs(m.misfit);
      ++np;
    } else {
      // horizon normal of the model at the point, from the scalar's gradient (finite differences)
      const double e = 2.0;
      std::vector<Vec3> q{dt.present + Vec3{e, 0, 0}, dt.present - Vec3{e, 0, 0}, dt.present + Vec3{0, e, 0}, dt.present - Vec3{0, e, 0},
                          dt.present + Vec3{0, 0, e}, dt.present - Vec3{0, 0, e}};
      EvalResult eq;
      inner->evaluate(q, kEvalScalars, eq);
      auto s = [&](std::size_t k) { return eq.s(static_cast<std::size_t>(dt.horizon), k); };
      Vec3 gr = normalized(Vec3{s(0) - s(1), s(2) - s(3), s(4) - s(5)});
      m.misfit = rad2deg(std::acos(std::clamp(dot(gr, dt.normal), -1.0, 1.0)));
      maxo = std::max(maxo, m.misfit);
      ++no;
    }
    if (dt.sigma > 0) {
      m.sigma = dt.sigma;
      m.normalised = m.misfit / dt.sigma;
    } else if (dt.sigma == 0) {
      m.sigma = 0;
    }
    m.loo = loo[i];
    mis.push_back(m);
  }
  log.value("implicit.data", std::to_string(np) + " points and tops, " + std::to_string(no) + " orientations");
  log.value("implicit.fields", std::to_string(groups.size()));
  std::string st;
  for (auto& t : solver_txt) st += (st.empty() ? "" : "; ") + t;
  log.value("implicit.solver", st + "; " + std::to_string(passes) + " pass(es)");
  if (free_values)
    for (std::size_t h = 0; h < nh; ++h)
      if (!iso[h]) log.value("implicit.horizon_value." + d.horizons[h].name, fmt_fixed(vals[h], 3) + " (fitted, free horizon values)");
  if (passes == iters && iters > 1 && moved >= converged)
    log.note("implicit: data restoration not converged after " + std::to_string(iters) + " passes (last shift " + fmt_fixed(moved, 3) + " m)");
  log.value("implicit.surface_spacing_m", fmt_shortest(sp));
  log.value("implicit.vertical_stretch", fmt_shortest(stretch));
  log.value("implicit.misfit_points_max_m", fmt_fixed(maxp, 4));
  log.value("implicit.misfit_points_mean_m", fmt_fixed(np ? sump / static_cast<double>(np) : 0.0, 4));
  log.value("implicit.misfit_orientations_max_deg", fmt_fixed(maxo, 3));
  if (d.implicit.data_influence) {
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < mis.size(); ++i)
      if (!std::isnan(mis[i].loo)) order.push_back(i);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return std::fabs(mis[a].loo) > std::fabs(mis[b].loo); });
    for (std::size_t k = 0; k < std::min<std::size_t>(10, order.size()); ++k) {
      auto& m = mis[order[k]];
      log.value("implicit.influence_" + std::to_string(k + 1), m.label + " " + m.kind + " of " + m.input + " at (" + fmt_fixed(m.p.x, 1) + ", " +
                                                                   fmt_fixed(m.p.y, 1) + ", " + fmt_fixed(m.p.z, 1) + "): leave-one-out error " +
                                                                   fmt_fixed(m.loo, 3) + " m");
    }
  }
  if (clamped) log.note("implicit: " + std::to_string(clamped) + " surface nodes had no crossing inside the domain and were held at its top or base");
  ModelInfo info = inner->info();
  info.engine = "implicit-rbf";
  info.provenance = provenance;
  return std::make_shared<ImplicitModel>(inner, info, std::move(mis));
}

namespace {
class ImplicitBuilder final : public StructureBuilder {
 public:
  bool accepts(const desc::Description& d) const override { return d.implicit.present; }
  ModelPtr build(const desc::Description& d, const BuildInputs& in, BuildLog& log) const override { return build_implicit_model(d, in, log); }
};
}  // namespace

STRATONIX_REGISTER(implicit_builder, registry::Kind::StructureBuilder, "implicit-rbf", "0.1.0",
                   "Data-conditioned: Hermite RBF stratigraphic field from points, orientations and tops; faults from the fault network",
                   [] { return std::make_unique<ImplicitBuilder>(); })

}  // namespace stratonix
