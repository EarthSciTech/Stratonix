// Stage 4C (STX-TN-16 sections 2.3 and 3, decision 23): layer-cake time-to-depth conversion with an
// instantaneous velocity linear in depth per layer, v(z) = v0 + k z, calibrated to wells.
#include "stratonix/engines/timedepth.hpp"

#include <algorithm>
#include <cmath>
#include <random>

#include "stratonix/core/error.hpp"
#include "stratonix/core/geometry.hpp"

namespace stratonix {

// dz/dt = v0 + k z  =>  z_b - z_t = (v0 + k z_t) (exp(k dt) - 1) / k  (dt one-way)
double layer_thickness(const LayerVelocity& v, double zt, double dt) {
  const double vt = v.v0 + v.k * zt;
  if (!(vt > 0)) throw Error(ErrorCode::Validation, "velocity at the layer top is not positive (v0 + k z = " + std::to_string(vt) + " m/s)");
  if (std::fabs(v.k) < 1e-12) return vt * dt;
  return vt * std::expm1(v.k * dt) / v.k;
}

double layer_time(const LayerVelocity& v, double zt, double zb) {
  const double vt = v.v0 + v.k * zt;
  if (!(vt > 0)) throw Error(ErrorCode::Validation, "velocity at the layer top is not positive");
  if (std::fabs(v.k) < 1e-12) return (zb - zt) / vt;
  const double r = v.k * (zb - zt) / vt;
  if (!(r > -1)) throw Error(ErrorCode::Validation, "depth beyond the velocity function's range (v0 + k z reaches zero)");
  return std::log1p(r) / v.k;
}

namespace {
struct Layers { std::vector<double> twt, z; std::vector<LayerVelocity> v; };
Layers layers_at(const VelocityModel& m, double x, double y) {
  if (!m.base_twt || !m.velocity) throw Error(ErrorCode::InvalidArgument, "velocity model without base times or velocities");
  Layers L;
  L.twt = m.base_twt(x, y);
  L.v = m.velocity(x, y);
  if (L.twt.size() != L.v.size() || L.twt.empty()) throw Error(ErrorCode::InvalidArgument, "velocity model: one velocity per layer");
  double z = 0, t = 0;
  for (std::size_t i = 0; i < L.twt.size(); ++i) {
    const double tb = std::max(t, L.twt[i]);  // a layer whose base time lies above its top is absent here
    z += layer_thickness(L.v[i], z, 0.5 * (tb - t));
    L.z.push_back(z);
    t = tb;
  }
  return L;
}
}  // namespace

std::vector<double> base_depths(const VelocityModel& m, double x, double y) { return layers_at(m, x, y).z; }

double time_to_depth(const VelocityModel& m, double x, double y, double twt) {
  const Layers L = layers_at(m, x, y);
  double zt = 0, tt = 0;
  for (std::size_t i = 0; i < L.twt.size(); ++i) {
    const double tb = std::max(tt, L.twt[i]);
    if (twt <= tb || i + 1 == L.twt.size()) return zt + layer_thickness(L.v[i], zt, 0.5 * (twt - tt));  // below the last base: its velocity
    zt = L.z[i];
    tt = tb;
  }
  return zt;
}

double depth_to_time(const VelocityModel& m, double x, double y, double z) {
  const Layers L = layers_at(m, x, y);
  double zt = 0, tt = 0;
  for (std::size_t i = 0; i < L.twt.size(); ++i) {
    if (z <= L.z[i] || i + 1 == L.twt.size()) return tt + 2 * layer_time(L.v[i], zt, z);
    zt = L.z[i];
    tt = std::max(tt, L.twt[i]);
  }
  return tt;
}

Calibration calibrate(const std::vector<WellTie>& wells, std::size_t layers, const std::vector<std::string>& exclude,
                      const std::vector<double>& fixed_k) {
  Calibration c;
  c.k.assign(layers, 0.0);
  c.k_sigma.assign(layers, 0.0);
  c.v0.assign(layers, std::vector<double>(wells.size(), kNaN));
  c.residual.assign(layers, std::vector<double>(wells.size(), kNaN));
  c.loo_residual.assign(layers, std::vector<double>(wells.size(), kNaN));
  std::vector<bool> use(wells.size(), true);
  for (std::size_t w = 0; w < wells.size(); ++w)
    for (auto& e : exclude)
      if (e == wells[w].well) {
        use[w] = false;
        c.excluded.push_back(wells[w].well + ": excluded by the user");
      }
  for (std::size_t L = 0; L < layers; ++L) {
    // k: weighted least squares of the interval velocity against the mid-depth of the layer
    // (Marsden 1989); weights are the user's, never re-weighted
    double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    std::size_t n = 0;
    struct Row { double zt, zb, dt; bool ok; };
    std::vector<Row> rows(wells.size());
    for (std::size_t w = 0; w < wells.size(); ++w) {
      const auto& t = wells[w];
      if (L >= t.top_depth.size() || L >= t.top_twt.size()) { rows[w].ok = false; continue; }
      const double zt = L ? t.top_depth[L - 1] : 0.0, tt = L ? t.top_twt[L - 1] : 0.0;
      rows[w] = {zt, t.top_depth[L], 0.5 * (t.top_twt[L] - tt),
                 std::isfinite(zt) && std::isfinite(tt) && std::isfinite(t.top_depth[L]) && std::isfinite(t.top_twt[L]) && t.top_depth[L] > zt && t.top_twt[L] > tt};
      if (!use[w] || !rows[w].ok) continue;
      const double vi = (rows[w].zb - zt) / rows[w].dt, zm = 0.5 * (zt + rows[w].zb), wt = t.weight;
      sw += wt; sx += wt * zm; sy += wt * vi; sxx += wt * zm * zm; sxy += wt * zm * vi;
      ++n;
    }
    const bool fixed = L < fixed_k.size() && std::isfinite(fixed_k[L]);
    if (fixed) c.k[L] = fixed_k[L];
    if (!fixed && n >= 2) {
      const double den = sw * sxx - sx * sx;
      if (std::fabs(den) > 1e-9 * std::max(1.0, sw * sxx)) {
        c.k[L] = (sw * sxy - sx * sy) / den;
        const double a = (sy - c.k[L] * sx) / sw;
        double ss = 0;
        for (std::size_t w = 0; w < wells.size(); ++w) {
          if (!use[w] || !rows[w].ok) continue;
          const double vi = (rows[w].zb - rows[w].zt) / rows[w].dt, zm = 0.5 * (rows[w].zt + rows[w].zb);
          ss += wells[w].weight * (vi - a - c.k[L] * zm) * (vi - a - c.k[L] * zm);
        }
        c.k_sigma[L] = n > 2 ? std::sqrt(ss / static_cast<double>(n - 2) * sw / den) : 0.0;
      }
    }
    // v0 at each well so that the layer's base ties its top
    for (std::size_t w = 0; w < wells.size(); ++w) {
      if (!rows[w].ok) continue;
      const double k = c.k[L], h = rows[w].zb - rows[w].zt, dt = rows[w].dt;
      const double vt = std::fabs(k) < 1e-12 ? h / dt : h * k / std::expm1(k * dt);
      c.v0[L][w] = vt - k * rows[w].zt;
      c.residual[L][w] = rows[w].zt + layer_thickness({c.v0[L][w], k}, rows[w].zt, dt) - rows[w].zb;
    }
    // leave-one-out: v0 from the other wells by inverse distance squared, the residual at this well
    for (std::size_t w = 0; w < wells.size(); ++w) {
      if (!rows[w].ok) continue;
      double num = 0, den = 0;
      for (std::size_t o = 0; o < wells.size(); ++o) {
        if (o == w || !use[o] || !std::isfinite(c.v0[L][o])) continue;
        const double d2 = std::max(1.0, (wells[o].x - wells[w].x) * (wells[o].x - wells[w].x) + (wells[o].y - wells[w].y) * (wells[o].y - wells[w].y));
        num += c.v0[L][o] / d2;
        den += 1.0 / d2;
      }
      if (den > 0) c.loo_residual[L][w] = rows[w].zt + layer_thickness({num / den, c.k[L]}, rows[w].zt, rows[w].dt) - rows[w].zb;
    }
  }
  return c;
}

DepthSpread depth_spread(const LayerVelocity& v, const VelocityUncertainty& u, double zt, double twt_layer, std::size_t members, std::uint64_t seed) {
  // seeded sampling: velocity parameters and the base pick, separately and together
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> g(0.0, 1.0);
  const double dt = 0.5 * twt_layer;
  auto sd = [](const std::vector<double>& a) {
    double m = 0, s = 0;
    for (double x : a) m += x;
    m /= static_cast<double>(a.size());
    for (double x : a) s += (x - m) * (x - m);
    return std::sqrt(s / static_cast<double>(a.size() - 1));
  };
  std::vector<double> zv, zp, zz;
  for (std::size_t i = 0; i < members; ++i) {
    const double ev = g(rng), ek = g(rng), ep = g(rng);
    const LayerVelocity vv{v.v0 + u.v0_sigma * ev, v.k + u.k_sigma * ek};
    const double dp = 0.5 * u.pick_sigma_ms / 1000.0 * ep;
    zv.push_back(layer_thickness(vv, zt, dt));
    zp.push_back(layer_thickness(v, zt, dt + dp));
    zz.push_back(layer_thickness(vv, zt, dt + dp));
  }
  DepthSpread s;
  s.members = members;
  s.velocity_sd = sd(zv);
  s.pick_sd = sd(zp);
  s.total_sd = sd(zz);
  return s;
}

// ---- section 3: correlated errors -------------------------------------------------------------
namespace {
std::uint64_t splitmix(std::uint64_t& x) {
  std::uint64_t z = (x += 0x9E3779B97F4A7C15ull);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}
struct Normal {
  std::uint64_t state;
  bool have = false;
  double spare = 0;
  explicit Normal(std::uint64_t seed) : state(seed) {}
  double uniform() { return (static_cast<double>(splitmix(state) >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
  double next() {
    if (have) { have = false; return spare; }
    const double u1 = uniform(), u2 = uniform(), r = std::sqrt(-2.0 * std::log(u1)), a = 2.0 * kPi * u2;
    spare = r * std::sin(a);
    have = true;
    return r * std::cos(a);
  }
};
/// Gaussian kernel taps for covariance exp(-(r/L)^2): standard deviation L/2 in lattice units of h,
/// normalised to unit sum of squares.
std::vector<double> taps(double L, double h) {
  const double s = 0.5 * L / h;
  const int m = static_cast<int>(std::ceil(3.0 * s));
  std::vector<double> g(static_cast<std::size_t>(2 * m + 1));
  double ss = 0;
  for (int i = -m; i <= m; ++i) {
    const double v = std::exp(-0.5 * (i / s) * (i / s));
    g[static_cast<std::size_t>(i + m)] = v;
    ss += v * v;
  }
  for (auto& v : g) v /= std::sqrt(ss);
  return g;
}
}  // namespace

double seeded_normal(std::uint64_t seed, std::uint64_t index) {
  std::uint64_t x = seed ^ (index * 0xD1B54A32D192ED03ull);
  Normal n(splitmix(x));
  return n.next();
}

std::vector<double> correlated_field(std::size_t nx, std::size_t ny, double dx, double dy, double sigma, double L, std::uint64_t seed) {
  std::vector<double> f(nx * ny, 0.0);
  if (f.empty() || !(sigma > 0)) return f;
  Normal g(seed);
  const double ex = dx * static_cast<double>(nx > 1 ? nx - 1 : 1), ey = dy * static_cast<double>(ny > 1 ? ny - 1 : 1);
  if (!(L > 0) || L > 20.0 * std::max(ex, ey)) {  // one value
    const double v = sigma * g.next();
    std::fill(f.begin(), f.end(), v);
    return f;
  }
  if (L < std::min(dx, dy)) {  // independent nodes
    for (auto& v : f) v = sigma * g.next();
    return f;
  }
  // coarse lattice of at most L/8 spacing over the extent, then bilinear to the nodes
  const double hx = std::max(dx, L / 8.0), hy = std::max(dy, L / 8.0);
  const std::size_t cx = static_cast<std::size_t>(std::ceil(ex / hx)) + 1, cy = static_cast<std::size_t>(std::ceil(ey / hy)) + 1;
  const auto tx = taps(L, hx), ty = taps(L, hy);
  const std::size_t mx = tx.size() / 2, my = ty.size() / 2, wx = cx + 2 * mx, wy = cy + 2 * my;
  std::vector<double> w(wx * wy);
  for (auto& v : w) v = g.next();
  std::vector<double> a(cx * wy, 0.0), c(cx * cy, 0.0);  // along x, then along y
  for (std::size_t j = 0; j < wy; ++j)
    for (std::size_t i = 0; i < cx; ++i) {
      double s = 0;
      for (std::size_t q = 0; q < tx.size(); ++q) s += tx[q] * w[j * wx + i + q];
      a[j * cx + i] = s;
    }
  for (std::size_t j = 0; j < cy; ++j)
    for (std::size_t i = 0; i < cx; ++i) {
      double s = 0;
      for (std::size_t q = 0; q < ty.size(); ++q) s += ty[q] * a[(j + q) * cx + i];
      c[j * cx + i] = s;
    }
  for (std::size_t j = 0; j < ny; ++j)
    for (std::size_t i = 0; i < nx; ++i) {
      const double u = static_cast<double>(i) * dx / hx, v = static_cast<double>(j) * dy / hy;
      const std::size_t i0 = std::min(static_cast<std::size_t>(u), cx > 1 ? cx - 2 : 0), j0 = std::min(static_cast<std::size_t>(v), cy > 1 ? cy - 2 : 0);
      const std::size_t i1 = std::min(i0 + 1, cx - 1), j1 = std::min(j0 + 1, cy - 1);
      const double fu = std::clamp(u - static_cast<double>(i0), 0.0, 1.0), fv = std::clamp(v - static_cast<double>(j0), 0.0, 1.0);
      f[j * nx + i] = sigma * ((1 - fu) * (1 - fv) * c[j0 * cx + i0] + fu * (1 - fv) * c[j0 * cx + i1] + (1 - fu) * fv * c[j1 * cx + i0] + fu * fv * c[j1 * cx + i1]);
    }
  return f;
}

void condition_field(std::vector<double>& f, std::size_t nx, std::size_t ny, double dx, double dy, double L,
                     const std::vector<std::pair<double, double>>& pts) {
  if (pts.empty() || f.empty()) return;
  const double ex = dx * static_cast<double>(nx > 1 ? nx - 1 : 1), ey = dy * static_cast<double>(ny > 1 ? ny - 1 : 1);
  const bool one = !(L > 0) || L > 20.0 * std::max(ex, ey);
  auto cov = [&](double ax, double ay, double bx, double by) {
    if (one) return 1.0;
    const double r2 = ((ax - bx) * dx) * ((ax - bx) * dx) + ((ay - by) * dy) * ((ay - by) * dy);
    return std::exp(-r2 / (L * L));
  };
  auto at = [&](double u, double v) {
    u = std::clamp(u, 0.0, static_cast<double>(nx - 1));
    v = std::clamp(v, 0.0, static_cast<double>(ny - 1));
    const std::size_t i0 = std::min(static_cast<std::size_t>(u), nx > 1 ? nx - 2 : 0), j0 = std::min(static_cast<std::size_t>(v), ny > 1 ? ny - 2 : 0);
    const std::size_t i1 = std::min(i0 + 1, nx - 1), j1 = std::min(j0 + 1, ny - 1);
    const double fu = u - static_cast<double>(i0), fv = v - static_cast<double>(j0);
    return (1 - fu) * (1 - fv) * f[j0 * nx + i0] + fu * (1 - fv) * f[j0 * nx + i1] + (1 - fu) * fv * f[j1 * nx + i0] + fu * fv * f[j1 * nx + i1];
  };
  const std::size_t n = pts.size();
  std::vector<double> C(n * n), b(n);
  for (std::size_t i = 0; i < n; ++i) {
    b[i] = at(pts[i].first, pts[i].second);
    for (std::size_t j = 0; j < n; ++j) C[i * n + j] = cov(pts[i].first, pts[i].second, pts[j].first, pts[j].second) + (i == j ? 1e-6 : 0.0);
  }
  // Cholesky solve C w = b
  for (std::size_t j = 0; j < n; ++j) {
    double d = C[j * n + j];
    for (std::size_t k = 0; k < j; ++k) d -= C[j * n + k] * C[j * n + k];
    d = std::sqrt(std::max(d, 1e-12));
    C[j * n + j] = d;
    for (std::size_t i = j + 1; i < n; ++i) {
      double s = C[i * n + j];
      for (std::size_t k = 0; k < j; ++k) s -= C[i * n + k] * C[j * n + k];
      C[i * n + j] = s / d;
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    double s = b[i];
    for (std::size_t k = 0; k < i; ++k) s -= C[i * n + k] * b[k];
    b[i] = s / C[i * n + i];
  }
  for (std::size_t i = n; i-- > 0;) {
    double s = b[i];
    for (std::size_t k = i + 1; k < n; ++k) s -= C[k * n + i] * b[k];
    b[i] = s / C[i * n + i];
  }
  // the correction is the field's continuous simple-kriging estimate from its values at the points
  // (removing exactly the values interpolated at the points, not their bilinear neighbours)
  for (std::size_t j = 0; j < ny; ++j)
    for (std::size_t i = 0; i < nx; ++i) {
      double e = 0;
      for (std::size_t q = 0; q < n; ++q) e += cov(static_cast<double>(i), static_cast<double>(j), pts[q].first, pts[q].second) * b[q];
      f[j * nx + i] -= e;
    }
}

}  // namespace stratonix
