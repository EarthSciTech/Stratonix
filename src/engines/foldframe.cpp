// Stage 4A (STX-TN-14 section 2): fold frames fitted to bedding orientations. The fold axis is the
// smallest eigenvector of the orientation tensor of the bedding normals (Ramsay 1967); bedding
// inclinations in the profile plane are fitted by a Fourier series of the profile coordinate
// (wavelength by search, order 1 to 3 by the Akaike information criterion; Laurent et al. 2016,
// Grose et al. 2017); the axial surface bisects the limbs and passes through a hinge.
#include "stratonix/engines/foldframe.hpp"

#include <algorithm>
#include <cmath>

#include "stratonix/core/error.hpp"
#include "stratonix/core/linalg.hpp"

namespace stratonix {
namespace {

struct Fourier {
  double W = 0;
  int order = 0;
  std::vector<double> c;  // a0, a1, b1, a2, b2, ...
  double rss = 0;
  double operator()(double x) const {
    double s = c.empty() ? 0.0 : c[0];
    for (int k = 1; k <= order; ++k) {
      double w = 2 * kPi * k * x / W;
      s += c[static_cast<std::size_t>(2 * k - 1)] * std::cos(w) + c[static_cast<std::size_t>(2 * k)] * std::sin(w);
    }
    return s;
  }
};

/// Least-squares Fourier series of order K with wavelength W.
Fourier fit_fixed(const std::vector<double>& x, const std::vector<double>& y, double W, int K) {
  const std::size_t m = static_cast<std::size_t>(2 * K + 1);
  std::vector<double> A(m * m, 0.0), b(m, 0.0), f(m);
  for (std::size_t i = 0; i < x.size(); ++i) {
    f[0] = 1;
    for (int k = 1; k <= K; ++k) {
      double w = 2 * kPi * k * x[i] / W;
      f[static_cast<std::size_t>(2 * k - 1)] = std::cos(w);
      f[static_cast<std::size_t>(2 * k)] = std::sin(w);
    }
    for (std::size_t r = 0; r < m; ++r) {
      b[r] += f[r] * y[i];
      for (std::size_t c = 0; c < m; ++c) A[r * m + c] += f[r] * f[c];
    }
  }
  for (std::size_t r = 0; r < m; ++r) A[r * m + r] += 1e-12 * static_cast<double>(x.size());
  solve_dense(A, b, m);
  Fourier F;
  F.W = W;
  F.order = K;
  F.c = b;
  for (std::size_t i = 0; i < x.size(); ++i) {
    double e = F(x[i]) - y[i];
    F.rss += e * e;
  }
  return F;
}

/// Wavelength by a log-spaced search over [wmin, wmax], then golden-section refinement.
Fourier fit_search(const std::vector<double>& x, const std::vector<double>& y, double wmin, double wmax, int K) {
  const int steps = 240;
  double best = std::numeric_limits<double>::max(), wb = wmin;
  for (int s = 0; s <= steps; ++s) {
    double W = wmin * std::pow(wmax / wmin, static_cast<double>(s) / steps);
    double r = fit_fixed(x, y, W, K).rss;
    if (r < best) { best = r; wb = W; }
  }
  const double step = std::pow(wmax / wmin, 1.0 / steps);
  double a = wb / step, b = wb * step;
  const double g = 0.5 * (std::sqrt(5.0) - 1);
  double c = b - g * (b - a), d = a + g * (b - a);
  double fc = fit_fixed(x, y, c, K).rss, fd = fit_fixed(x, y, d, K).rss;
  for (int it = 0; it < 60; ++it) {
    if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = fit_fixed(x, y, c, K).rss; }
    else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = fit_fixed(x, y, d, K).rss; }
  }
  return fit_fixed(x, y, 0.5 * (a + b), K);
}

/// Best order by AIC (wavelength searched for each order around the order-1 wavelength).
Fourier fit_best(const std::vector<double>& x, const std::vector<double>& y, double wmin, double wmax, int maxK) {
  Fourier best = fit_search(x, y, wmin, wmax, 1);
  const double n = static_cast<double>(x.size());
  auto aic = [&](const Fourier& F) { return n * std::log(std::max(F.rss, 1e-300) / n) + 2.0 * (2 * F.order + 1); };
  double ab = aic(best);
  for (int K = 2; K <= maxK; ++K) {
    Fourier F = fit_search(x, y, best.W / 1.5, best.W * 1.5, K);
    if (aic(F) < ab - 1e-9) { best = F; ab = aic(F); }
  }
  return best;
}

}  // namespace

double FoldFrame::frame_coordinate(const Vec3& x) const { return dot(x - axial_point, axial_normal); }

double FoldFrame::rotation_angle(double c) const {
  Fourier F;
  F.W = wavelength;
  F.order = order;
  F.c = fourier;
  return F(c);
}

Vec3 FoldFrame::bedding_normal(const Vec3& x) const {
  Vec3 up = normalized(cross(axial_normal, axis));
  if (up.z < 0) up = up * -1.0;
  double a = rotation_angle(frame_coordinate(x));
  Vec3 n = normalized(up * std::cos(a) - axial_normal * std::sin(a));
  return n.z < 0 ? n * -1.0 : n;
}

FoldFrame fit_fold_frame(const std::vector<OrientationDatum>& bedding, const std::vector<OrientationDatum>& foliation,
                         const FoldFrameOptions& o) {
  FoldFrame F;
  F.data = bedding.size();
  if (bedding.size() < 10) { F.reason = "a fold frame needs at least 10 bedding orientations"; return F; }
  // 1. fold axis: smallest eigenvector of the orientation tensor
  double T[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  Vec3 c0{0, 0, 0};
  for (auto& b : bedding) {
    Vec3 n = normalized(b.n);
    double v[3] = {n.x, n.y, n.z};
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) T[r * 3 + c] += v[r] * v[c];
    c0 += b.p;
  }
  c0 = c0 / static_cast<double>(bedding.size());
  double ev[3], vec[9];
  symmetric_eigen3(T, ev, vec);
  Vec3 axis{vec[0], vec[3], vec[6]};
  axis = normalized(axis);
  F.cylindricity = ev[1] > 0 ? ev[0] / ev[1] : 1.0;
  if (axis.z > 0 || (std::fabs(axis.z) < 1e-12 && std::atan2(axis.x, axis.y) < 0)) axis = axis * -1.0;
  if (std::fabs(axis.z) > 0.99) { F.reason = "the fold axis is (near) vertical: the profile plane is horizontal"; return F; }
  F.axis = axis;
  F.plunge = std::asin(std::clamp(-axis.z, -1.0, 1.0));
  F.trend = std::atan2(axis.x, axis.y);
  if (F.trend < 0) F.trend += 2 * kPi;
  // 2. profile frame and bedding inclinations
  Vec3 across = normalized(cross(axis, Vec3{0, 0, 1}));
  Vec3 e2 = normalized(cross(across, axis));
  if (e2.z < 0) e2 = e2 * -1.0;
  // axial foliation, when measured, sets the axial-surface orientation directly
  double psi = kNaN;
  if (foliation.size() >= 3) {
    Vec3 s{0, 0, 0};
    for (auto& f : foliation) {
      Vec3 n = normalized(f.n);
      if (dot(n, across) < 0) n = n * -1.0;
      s += n;
    }
    s = normalized(s);
    psi = std::atan2(dot(s, e2), dot(s, across));  // axial normal = (cos psi, sin psi) in (across, e2)
  }
  std::vector<double> xp, th;
  double xmin = std::numeric_limits<double>::max(), xmax = -xmin;
  for (auto& b : bedding) {
    Vec3 n = normalized(b.n);
    xp.push_back(dot(b.p - c0, across));
    th.push_back(std::atan2(-dot(n, across), dot(n, e2)));
    xmin = std::min(xmin, xp.back());
    xmax = std::max(xmax, xp.back());
  }
  const double ext = std::max(xmax - xmin, 1.0);
  const double wmin = o.wavelength > 0 ? o.wavelength / 2 : ext / 10, wmax = o.wavelength > 0 ? o.wavelength * 2 : 4 * ext;
  Fourier P = fit_best(xp, th, wmin, wmax, std::max(1, o.max_order));
  double tmax = -1e9, tmin = 1e9;
  for (int s = 0; s <= 2000; ++s) {
    double v = P(xmin + ext * s / 2000.0);
    tmax = std::max(tmax, v);
    tmin = std::min(tmin, v);
  }
  if (std::isnan(psi)) psi = 0.5 * (tmax + tmin);
  F.interlimb = kPi - (tmax - tmin);
  // 3. frame rotated by psi: axial normal and the fold rotation angle against the frame coordinate
  Vec3 an = normalized(across * std::cos(psi) + e2 * std::sin(psi));
  std::vector<double> cc(xp.size()), al(xp.size());
  for (std::size_t i = 0; i < xp.size(); ++i) {
    cc[i] = dot(bedding[i].p - c0, an);
    al[i] = th[i] - psi;
  }
  Fourier A = fit_best(cc, al, P.W / 1.5, P.W * 1.5, std::max(1, o.max_order));
  // hinge (antiform crest: alpha falls through zero) nearest the centroid, within the data
  double cmin = *std::min_element(cc.begin(), cc.end()), cmax = *std::max_element(cc.begin(), cc.end());
  double ch = 0, bestd = std::numeric_limits<double>::max();
  const int ns = 4000;
  for (int s = 0; s < ns; ++s) {
    double a = cmin + (cmax - cmin) * s / ns, b = cmin + (cmax - cmin) * (s + 1) / ns;
    double fa = A(a), fb = A(b);
    if (fa > 0 && fb <= 0) {
      double z = a + (b - a) * fa / (fa - fb);
      if (std::fabs(z) < bestd) { bestd = std::fabs(z); ch = z; }
    }
  }
  F.axial_point = c0 + an * ch;
  F.axial_normal = an;
  for (auto& v : cc) v -= ch;
  A = fit_fixed(cc, al, A.W, A.order);
  F.fourier = A.c;
  F.order = A.order;
  F.wavelength = A.W;
  // 4. misfit of the predicted normals
  F.ok = true;
  double ss = 0;
  for (auto& b : bedding) {
    double a = std::acos(std::clamp(dot(F.bedding_normal(b.p), normalized(b.n.z < 0 ? b.n * -1.0 : b.n)), -1.0, 1.0));
    ss += a * a;
  }
  F.rms = std::sqrt(ss / static_cast<double>(bedding.size()));
  if (F.cylindricity > o.max_cylindricity)
    F.reason = "cylindricity " + std::to_string(F.cylindricity) + " above " + std::to_string(o.max_cylindricity) +
               ": the fold is not cylindrical; one axis is fitted (a varying axis is not yet built)";
  return F;
}

}  // namespace stratonix
