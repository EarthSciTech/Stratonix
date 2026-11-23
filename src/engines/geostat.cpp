// Stage 4D (STX-TN-17 section 1, decision 24): geostatistics from the published GSLIB algorithms
// (Deutsch and Journel 1998: gamv, sgsim, sisim, declus) and object-based facies, on points in
// geological-time coordinates.
//
// Reproducibility (STX-TN-17 section 1.3): the random stream (SplitMix64, Marsaglia's polar
// method), the path, the neighbour order (sorted by distance, then node index) and the arithmetic
// are fixed. The transcendental functions that enter the simulated values (exp, log) are evaluated
// here from + - * / and sqrt only, which IEEE 754 rounds identically on every platform, so a seed
// gives the same bytes on Windows and Linux.
#include "stratonix/engines/geostat.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <array>
#include <limits>
#include <map>
#include <numeric>
#include <unordered_map>

#include "stratonix/core/error.hpp"

namespace stratonix::geo {
namespace {

constexpr double kPiG = 3.14159265358979323846;
constexpr double kLn2 = 0.693147180559945309417232121458;

// ---- deterministic elementary functions --------------------------------------------------------
/// exp(x) by range reduction x = k ln2 + r, |r| <= ln2/2, and a degree-13 Taylor polynomial in
/// Horner form (relative error below 1e-15); 2^k by exact scaling.
double det_exp(double x) {
  if (x < -745.0) return 0.0;
  if (x > 709.0) return std::numeric_limits<double>::infinity();
  const double kf = std::floor(x / kLn2 + 0.5);
  const double r = x - kf * kLn2;
  double p = 1.0;
  for (int n = 13; n >= 1; --n) p = 1.0 + p * r / n;
  return std::ldexp(p, static_cast<int>(kf));
}
/// log(x) for x > 0: x = m 2^e with m in [sqrt(1/2), sqrt(2)), log m = 2 atanh(s), s = (m-1)/(m+1),
/// series to s^31 (|s| < 0.172: error below 1e-16).
double det_log(double x) {
  if (!(x > 0)) return -std::numeric_limits<double>::infinity();
  int e = 0;
  double m = std::frexp(x, &e);  // m in [0.5, 1)
  if (m < 0.70710678118654752440) { m *= 2.0; --e; }
  const double s = (m - 1.0) / (m + 1.0), s2 = s * s;
  double t = 0.0;
  for (int k = 31; k >= 1; k -= 2) t = 1.0 / k + s2 * t;
  return 2.0 * s * t + e * kLn2;
}

double norm_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }  // not used in simulated values
/// Acklam's rational approximation with one Halley step (as stage 3C).
double norm_quantile(double p) {
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
  p = std::clamp(p, 1e-12, 1.0 - 1e-12);
  double x;
  if (p < 0.02425) {
    const double q = std::sqrt(-2 * std::log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else if (p > 1 - 0.02425) {
    const double q = std::sqrt(-2 * std::log(1 - p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else {
    const double q = p - 0.5, r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }
  const double e = norm_cdf(x) - p, u = e * std::sqrt(2 * kPiG) * std::exp(x * x / 2);
  return x - u / (1 + x * u / 2);
}

// ---- random stream -------------------------------------------------------------------------------
struct Rng {
  std::uint64_t s;
  bool have = false;
  double spare = 0;
  explicit Rng(std::uint64_t seed) : s(seed) {}
  std::uint64_t next() {
    std::uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  double uniform() { return (static_cast<double>(next() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }  // (0, 1)
  std::size_t below(std::size_t n) { return static_cast<std::size_t>(uniform() * static_cast<double>(n)) % n; }
  /// Marsaglia's polar method: log and sqrt only (deterministic).
  double normal() {
    if (have) { have = false; return spare; }
    double u, v, q;
    do {
      u = 2 * uniform() - 1;
      v = 2 * uniform() - 1;
      q = u * u + v * v;
    } while (q >= 1 || q == 0);
    const double f = std::sqrt(-2 * det_log(q) / q);
    spare = v * f;
    have = true;
    return u * f;
  }
};

// ---- variogram geometry ----------------------------------------------------------------------------
double reduced(const VariogramStructure& s, const Vec3& h) {
  const double a = s.azimuth * kPiG / 180.0, sa = std::sin(a), ca = std::cos(a);
  const double hM = h.x * sa + h.y * ca, hm = h.x * ca - h.y * sa;
  const double aM = s.range_major > 0 ? s.range_major : 1e-12;
  const double am = s.range_minor > 0 ? s.range_minor : aM;
  const double av = s.range_vertical > 0 ? s.range_vertical : aM;
  return std::sqrt((hM / aM) * (hM / aM) + (hm / am) * (hm / am) + (h.z / av) * (h.z / av));
}
double shape(const std::string& type, double r) {
  if (type == "spherical") return r >= 1 ? 1.0 : r * (1.5 - 0.5 * r * r);
  if (type == "exponential") return 1.0 - det_exp(-3.0 * r);
  if (type == "gaussian") return 1.0 - det_exp(-3.0 * r * r);
  throw Error(ErrorCode::InvalidArgument, "variogram: unknown structure type '" + type + "' (spherical, exponential, gaussian)");
}

// ---- small dense solver --------------------------------------------------------------------------
/// Solves A x = b for a symmetric positive definite A (n x n, row-major) by Cholesky; false when
/// A is not positive definite.
bool cholesky_solve(std::vector<double>& A, std::vector<double>& b, std::size_t n) {
  for (std::size_t j = 0; j < n; ++j) {
    double d = A[j * n + j];
    for (std::size_t k = 0; k < j; ++k) d -= A[j * n + k] * A[j * n + k];
    if (!(d > 1e-14)) return false;
    d = std::sqrt(d);
    A[j * n + j] = d;
    for (std::size_t i = j + 1; i < n; ++i) {
      double s = A[i * n + j];
      for (std::size_t k = 0; k < j; ++k) s -= A[i * n + k] * A[j * n + k];
      A[i * n + j] = s / d;
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    double s = b[i];
    for (std::size_t k = 0; k < i; ++k) s -= A[i * n + k] * b[k];
    b[i] = s / A[i * n + i];
  }
  for (std::size_t i = n; i-- > 0;) {
    double s = b[i];
    for (std::size_t k = i + 1; k < n; ++k) s -= A[k * n + i] * b[k];
    b[i] = s / A[i * n + i];
  }
  return true;
}
/// General solve by Gaussian elimination with partial pivoting (ordinary kriging is indefinite).
bool gauss_solve(std::vector<double>& A, std::vector<double>& b, std::size_t n) {
  for (std::size_t c = 0; c < n; ++c) {
    std::size_t p = c;
    for (std::size_t r = c + 1; r < n; ++r)
      if (std::fabs(A[r * n + c]) > std::fabs(A[p * n + c])) p = r;
    if (std::fabs(A[p * n + c]) < 1e-14) return false;
    if (p != c) {
      for (std::size_t k = 0; k < n; ++k) std::swap(A[c * n + k], A[p * n + k]);
      std::swap(b[c], b[p]);
    }
    for (std::size_t r = c + 1; r < n; ++r) {
      const double f = A[r * n + c] / A[c * n + c];
      if (f == 0) continue;
      for (std::size_t k = c; k < n; ++k) A[r * n + k] -= f * A[c * n + k];
      b[r] -= f * b[c];
    }
  }
  for (std::size_t i = n; i-- > 0;) {
    double s = b[i];
    for (std::size_t k = i + 1; k < n; ++k) s -= A[i * n + k] * b[k];
    b[i] = s / A[i * n + i];
  }
  return true;
}

// ---- neighbour search ------------------------------------------------------------------------------
/// Buckets in coordinates scaled by the search radii (the ellipsoid becomes the unit sphere).
class Search {
 public:
  Search(const std::vector<Vec3>& nodes, const SearchOptions& s) : nodes_(nodes), s_(s) {
    const double a = s.azimuth * kPiG / 180.0;
    sa_ = std::sin(a);
    ca_ = std::cos(a);
    rM_ = s.radius_major > 0 ? s.radius_major : 1e30;
    rm_ = s.radius_minor > 0 ? s.radius_minor : rM_;
    rv_ = s.radius_vertical > 0 ? s.radius_vertical : rM_;
    scaled_.resize(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) scaled_[i] = scale(nodes[i]);
  }
  void insert(std::size_t node) { cells_[key(cell_of(scaled_[node]))].push_back(node); }
  /// Nearest inserted nodes within the ellipsoid, at most max_data (and max_per_octant per octant),
  /// sorted by scaled distance then node index.
  std::vector<std::size_t> find(std::size_t node) const {
    const Vec3 p = scaled_[node];
    const auto c = cell_of(p);
    std::vector<std::pair<double, std::size_t>> cand;
    for (long dk = -1; dk <= 1; ++dk)
      for (long dj = -1; dj <= 1; ++dj)
        for (long di = -1; di <= 1; ++di) {
          auto it = cells_.find(key({c[0] + di, c[1] + dj, c[2] + dk}));
          if (it == cells_.end()) continue;
          for (std::size_t q : it->second) {
            if (q == node) continue;
            const Vec3 d = scaled_[q] - p;
            const double r2 = d.x * d.x + d.y * d.y + d.z * d.z;
            if (r2 <= 1.0) cand.push_back({r2, q});
          }
        }
    // (distance, index) is a total order, so sorting only the first max_data gives the same nodes
    // in the same order as a full sort; an octant limit may skip candidates, so it sorts them all
    if (!s_.max_per_octant && s_.max_data > 0 && cand.size() > s_.max_data)
      std::partial_sort(cand.begin(), cand.begin() + static_cast<std::ptrdiff_t>(s_.max_data), cand.end());
    else
      std::sort(cand.begin(), cand.end());
    std::vector<std::size_t> out;
    std::size_t oct[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (auto& [r2, q] : cand) {
      if (out.size() >= s_.max_data) break;
      if (s_.max_per_octant) {
        const Vec3 d = nodes_[q] - nodes_[node];
        const int o = (d.x >= 0 ? 1 : 0) + (d.y >= 0 ? 2 : 0) + (d.z >= 0 ? 4 : 0);
        if (oct[o] >= s_.max_per_octant) continue;
        ++oct[o];
      }
      out.push_back(q);
    }
    return out;
  }

 private:
  Vec3 scale(const Vec3& x) const {
    const double xM = x.x * sa_ + x.y * ca_, xm = x.x * ca_ - x.y * sa_;
    return {xM / rM_, xm / rm_, x.z / rv_};
  }
  std::array<long, 3> cell_of(const Vec3& p) const {
    return {static_cast<long>(std::floor(p.x)), static_cast<long>(std::floor(p.y)), static_cast<long>(std::floor(p.z))};
  }
  static std::uint64_t key(const std::array<long, 3>& c) {
    return (static_cast<std::uint64_t>(c[0] + (1 << 20)) << 42) ^ (static_cast<std::uint64_t>(c[1] + (1 << 20)) << 21) ^
           static_cast<std::uint64_t>(c[2] + (1 << 20));
  }
  const std::vector<Vec3>& nodes_;
  SearchOptions s_;
  double sa_ = 0, ca_ = 1, rM_ = 1, rm_ = 1, rv_ = 1;
  std::vector<Vec3> scaled_;
  std::unordered_map<std::uint64_t, std::vector<std::size_t>> cells_;  // lookups only: order of iteration never used
};

/// Kriging weights of neighbours `nb` for node `x` (simple, or ordinary with a Lagrange term);
/// returns false when the system is singular. `var` receives the kriging variance.
bool krige(const std::vector<Vec3>& nodes, const std::vector<std::size_t>& nb, std::size_t x, const Variogram& v, bool ordinary,
           std::vector<double>& w, double& var) {
  const std::size_t n = nb.size(), m = ordinary ? n + 1 : n;
  std::vector<double> A(m * m, 0.0), b(m, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) A[i * m + j] = A[j * m + i] = v.covariance(nodes[nb[j]] - nodes[nb[i]]);
    b[i] = v.covariance(nodes[x] - nodes[nb[i]]);
  }
  if (ordinary) {
    for (std::size_t i = 0; i < n; ++i) A[i * m + n] = A[n * m + i] = 1.0;
    A[n * m + n] = 0.0;
    b[n] = 1.0;
  }
  const std::vector<double> c(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(n));
  const bool ok = ordinary ? gauss_solve(A, b, m) : cholesky_solve(A, b, m);
  if (!ok) return false;
  w.assign(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(n));
  var = v.sill();
  for (std::size_t i = 0; i < n; ++i) var -= w[i] * c[i];
  if (ordinary) var -= b[n];
  return true;
}

/// Random path through the nodes not fixed by data, visited coarse to fine (the multiple-grid
/// search of GSLIB sgsim, generalised to scattered nodes): at each level one node per cell of a
/// lattice in search-scaled coordinates (cells of 1/2, 1/4, 1/8 of the search radius), the rest
/// last; random order within a level (Fisher-Yates on the stream). Long-range structure is then
/// conditioned by nodes the limited neighbourhood would otherwise not reach.
std::vector<std::size_t> random_path(const std::vector<Vec3>& nodes, const std::vector<char>& fixed, const SearchOptions& s, Rng& r) {
  const std::size_t n = nodes.size();
  std::vector<std::size_t> order;
  for (std::size_t i = 0; i < n; ++i)
    if (!fixed[i]) order.push_back(i);
  for (std::size_t i = order.size(); i > 1; --i) std::swap(order[i - 1], order[r.below(i)]);
  const double a = s.azimuth * kPiG / 180.0, sa = std::sin(a), ca = std::cos(a);
  const double rM = s.radius_major > 0 ? s.radius_major : 0, rm = s.radius_minor > 0 ? s.radius_minor : rM,
               rv = s.radius_vertical > 0 ? s.radius_vertical : rM;
  if (!(rM > 0)) return order;
  std::vector<char> taken(n, 0);
  std::vector<std::size_t> path;
  path.reserve(order.size());
  for (double frac : {0.5, 0.25, 0.125}) {
    std::map<std::array<long, 3>, char> used;
    std::vector<std::size_t> level;
    for (std::size_t i : order) {
      if (taken[i]) continue;
      const Vec3& x = nodes[i];
      const double xM = (x.x * sa + x.y * ca) / (rM * frac), xm = (x.x * ca - x.y * sa) / (rm * frac), xz = x.z / (rv * frac);
      const std::array<long, 3> c = {static_cast<long>(std::floor(xM)), static_cast<long>(std::floor(xm)), static_cast<long>(std::floor(xz))};
      if (used.emplace(c, 1).second) {
        taken[i] = 1;
        level.push_back(i);
      }
    }
    path.insert(path.end(), level.begin(), level.end());
  }
  for (std::size_t i : order)
    if (!taken[i]) path.push_back(i);
  return path;
}

}  // namespace

// ---- variograms -------------------------------------------------------------------------------------
double Variogram::sill() const {
  double s = nugget;
  for (auto& st : structures) s += st.contribution;
  return s;
}
double Variogram::gamma(const Vec3& h) const {
  if (h.x == 0 && h.y == 0 && h.z == 0) return 0.0;
  double g = nugget;
  for (auto& st : structures) g += st.contribution * shape(st.type, reduced(st, h));
  return g;
}
double Variogram::covariance(const Vec3& h) const { return sill() - gamma(h); }

ExperimentalVariogram experimental_variogram(const std::vector<Vec3>& points, const std::vector<double>& values, const VariogramDirection& d) {
  if (points.size() != values.size()) throw Error(ErrorCode::InvalidArgument, "experimental variogram: one value per point");
  if (!(d.lag > 0) || d.lags == 0) throw Error(ErrorCode::InvalidArgument, "experimental variogram: lag and number of lags are positive");
  ExperimentalVariogram e;
  e.distance.assign(d.lags, 0.0);
  e.gamma.assign(d.lags, 0.0);
  e.pairs.assign(d.lags, 0);
  const double az = d.azimuth * kPiG / 180.0, dp = d.dip * kPiG / 180.0;
  const Vec3 u{std::sin(az) * std::cos(dp), std::cos(az) * std::cos(dp), std::sin(dp)};
  const bool omni = d.angle_tolerance >= 90.0;
  const double cos_tol = std::cos(d.angle_tolerance * kPiG / 180.0);
  const double tol = d.lag_tolerance > 0 ? d.lag_tolerance : 0.5 * d.lag;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (!std::isfinite(values[i])) continue;
    for (std::size_t j = i + 1; j < points.size(); ++j) {
      if (!std::isfinite(values[j])) continue;
      const Vec3 h = points[j] - points[i];
      const double dist = std::sqrt(h.x * h.x + h.y * h.y + h.z * h.z);
      if (dist == 0) continue;
      const double lf = dist / d.lag;
      const long l = std::lround(lf);
      if (l < 1 || static_cast<std::size_t>(l) > d.lags || std::fabs(dist - static_cast<double>(l) * d.lag) > tol) continue;
      if (!omni) {
        const double proj = h.x * u.x + h.y * u.y + h.z * u.z;
        if (std::fabs(proj) / dist < cos_tol) continue;
        const Vec3 perp{h.x - proj * u.x, h.y - proj * u.y, h.z - proj * u.z};
        if (std::sqrt(perp.x * perp.x + perp.y * perp.y + perp.z * perp.z) > d.bandwidth) continue;
      }
      const std::size_t k = static_cast<std::size_t>(l - 1);
      e.distance[k] += dist;
      e.gamma[k] += 0.5 * (values[j] - values[i]) * (values[j] - values[i]);
      ++e.pairs[k];
    }
  }
  for (std::size_t k = 0; k < d.lags; ++k)
    if (e.pairs[k]) {
      e.distance[k] /= static_cast<double>(e.pairs[k]);
      e.gamma[k] /= static_cast<double>(e.pairs[k]);
    } else {
      e.distance[k] = static_cast<double>(k + 1) * d.lag;
      e.gamma[k] = std::numeric_limits<double>::quiet_NaN();
    }
  return e;
}

Variogram fit_variogram(const ExperimentalVariogram& e, const Variogram& start, bool fit_nugget, std::size_t min_pairs) {
  // parameters: log of each contribution and range (major, minor, vertical), optionally the nugget
  std::vector<double> x;
  for (auto& s : start.structures) {
    x.push_back(std::log(std::max(1e-9, s.contribution)));
    x.push_back(std::log(std::max(1e-9, s.range_major)));
  }
  if (fit_nugget) x.push_back(std::log(std::max(1e-9, start.nugget + 1e-6)));
  auto make = [&](const std::vector<double>& p) {
    Variogram v = start;
    for (std::size_t k = 0; k < v.structures.size(); ++k) {
      auto& s = v.structures[k];
      const double r0 = start.structures[k].range_major;
      s.contribution = std::exp(p[2 * k]);
      const double f = std::exp(p[2 * k + 1]) / (r0 > 0 ? r0 : 1.0);
      s.range_major = std::exp(p[2 * k + 1]);
      s.range_minor = start.structures[k].range_minor * f;
      s.range_vertical = start.structures[k].range_vertical * f;
    }
    if (fit_nugget) v.nugget = std::exp(p.back());
    return v;
  };
  auto cost = [&](const std::vector<double>& p) {
    const Variogram v = make(p);
    double c = 0;
    for (std::size_t k = 0; k < e.gamma.size(); ++k) {
      if (e.pairs[k] < std::max<std::size_t>(1, min_pairs) || !std::isfinite(e.gamma[k])) continue;
      const double gm = v.gamma({e.distance[k], 0, 0});  // along the direction of the experimental variogram
      c += static_cast<double>(e.pairs[k]) * (e.gamma[k] - gm) * (e.gamma[k] - gm) / std::max(1e-12, gm * gm);
    }
    return c;
  };
  // Nelder-Mead (deterministic start simplex and iteration count)
  const std::size_t n = x.size();
  if (n == 0) return start;
  std::vector<std::vector<double>> s(n + 1, x);
  for (std::size_t i = 0; i < n; ++i) s[i + 1][i] += 0.5;
  std::vector<double> f(n + 1);
  for (std::size_t i = 0; i <= n; ++i) f[i] = cost(s[i]);
  for (int it = 0; it < 2000; ++it) {
    std::vector<std::size_t> o(n + 1);
    std::iota(o.begin(), o.end(), 0);
    std::sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return f[a] < f[b] || (f[a] == f[b] && a < b); });
    const std::size_t best = o[0], worst = o[n], second = o[n - 1];
    if (std::fabs(f[worst] - f[best]) < 1e-12 * (1 + std::fabs(f[best]))) break;
    std::vector<double> cen(n, 0.0);
    for (std::size_t i = 0; i <= n; ++i)
      if (i != worst)
        for (std::size_t k = 0; k < n; ++k) cen[k] += s[i][k] / static_cast<double>(n);
    auto along = [&](double t) {
      std::vector<double> p(n);
      for (std::size_t k = 0; k < n; ++k) p[k] = cen[k] + t * (s[worst][k] - cen[k]);
      return p;
    };
    auto xr = along(-1.0);
    const double fr = cost(xr);
    if (fr < f[best]) {
      auto xe = along(-2.0);
      const double fe = cost(xe);
      if (fe < fr) { s[worst] = xe; f[worst] = fe; } else { s[worst] = xr; f[worst] = fr; }
    } else if (fr < f[second]) {
      s[worst] = xr;
      f[worst] = fr;
    } else {
      auto xc = along(fr < f[worst] ? -0.5 : 0.5);
      const double fc = cost(xc);
      if (fc < std::min(fr, f[worst])) {
        s[worst] = xc;
        f[worst] = fc;
      } else {
        for (std::size_t i = 0; i <= n; ++i)
          if (i != best) {
            for (std::size_t k = 0; k < n; ++k) s[i][k] = s[best][k] + 0.5 * (s[i][k] - s[best][k]);
            f[i] = cost(s[i]);
          }
      }
    }
  }
  std::size_t b = 0;
  for (std::size_t i = 1; i <= n; ++i)
    if (f[i] < f[b]) b = i;
  return make(s[b]);
}

// ---- normal scores ------------------------------------------------------------------------------------
NormalScore::NormalScore(std::vector<double> values, std::vector<double> weights, double zmin, double zmax) : zmin_(zmin), zmax_(zmax) {
  if (values.empty()) throw Error(ErrorCode::InvalidArgument, "normal score: no values");
  if (weights.size() != values.size()) weights.assign(values.size(), 1.0);
  std::vector<std::size_t> o(values.size());
  std::iota(o.begin(), o.end(), 0);
  std::sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return values[a] < values[b] || (values[a] == values[b] && a < b); });
  double W = 0;
  for (double w : weights) W += w;
  // cumulative probability at the middle of each value's weight; equal values share one point
  double cum = 0;
  std::vector<double> z, p;
  for (std::size_t k = 0; k < o.size();) {
    std::size_t e = k;
    double w = 0;
    while (e < o.size() && values[o[e]] == values[o[k]]) w += weights[o[e++]];
    z.push_back(values[o[k]]);
    p.push_back((cum + 0.5 * w) / W);
    cum += w;
    k = e;
  }
  if (!(zmin_ < z.front())) zmin_ = z.front() - 1e-9 * (1 + std::fabs(z.front()));
  if (!(zmax_ > z.back())) zmax_ = z.back() + 1e-9 * (1 + std::fabs(z.back()));
  z_.push_back(zmin_);
  y_.push_back(0.0);  // cumulative probability 0 at zmin (linear tail)
  for (std::size_t k = 0; k < z.size(); ++k) { z_.push_back(z[k]); y_.push_back(p[k]); }
  z_.push_back(zmax_);
  y_.push_back(1.0);
}
double NormalScore::forward(double z) const {
  z = std::clamp(z, zmin_, zmax_);
  auto it = std::upper_bound(z_.begin(), z_.end(), z);
  std::size_t k = it == z_.end() ? z_.size() - 1 : static_cast<std::size_t>(it - z_.begin());
  if (k == 0) k = 1;
  const double t = (z - z_[k - 1]) / (z_[k] - z_[k - 1]);
  return norm_quantile(y_[k - 1] + t * (y_[k] - y_[k - 1]));
}
double NormalScore::back(double y) const {
  const double p = norm_cdf(y);
  auto it = std::upper_bound(y_.begin(), y_.end(), p);
  std::size_t k = it == y_.end() ? y_.size() - 1 : static_cast<std::size_t>(it - y_.begin());
  if (k == 0) k = 1;
  const double t = (p - y_[k - 1]) / (y_[k] - y_[k - 1]);
  return z_[k - 1] + t * (z_[k] - z_[k - 1]);
}

std::vector<double> decluster_weights(const std::vector<Vec3>& points, double cell) {
  if (points.empty()) return {};
  if (!(cell > 0)) return std::vector<double>(points.size(), 1.0);
  Vec3 lo = points[0];
  for (auto& p : points) {
    lo.x = std::min(lo.x, p.x);
    lo.y = std::min(lo.y, p.y);
    lo.z = std::min(lo.z, p.z);
  }
  std::map<std::array<long, 3>, std::size_t> count;
  std::vector<std::array<long, 3>> key(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    key[i] = {static_cast<long>(std::floor((points[i].x - lo.x) / cell)), static_cast<long>(std::floor((points[i].y - lo.y) / cell)),
              static_cast<long>(std::floor((points[i].z - lo.z) / cell))};
    ++count[key[i]];
  }
  std::vector<double> w(points.size());
  double s = 0;
  for (std::size_t i = 0; i < points.size(); ++i) s += (w[i] = 1.0 / static_cast<double>(count[key[i]]));
  for (auto& x : w) x *= static_cast<double>(points.size()) / s;
  return w;
}

// ---- sequential Gaussian simulation -------------------------------------------------------------
std::vector<double> sgs(const std::vector<Vec3>& nodes, const std::vector<ConditioningPoint>& data, const Variogram& v, const SearchOptions& s,
                        std::uint64_t seed, const std::vector<double>& local_mean) {
  const std::size_t n = nodes.size();
  if (!local_mean.empty() && local_mean.size() != n) throw Error(ErrorCode::InvalidArgument, "sgs: one local mean per node");
  std::vector<double> y(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<char> fixed(n, 0);
  Search search(nodes, s);
  for (auto& c : data) {
    if (c.node >= n) throw Error(ErrorCode::InvalidArgument, "sgs: conditioning node outside the nodes");
    if (fixed[c.node]) continue;  // the first value at a node is kept
    y[c.node] = c.value;
    fixed[c.node] = 1;
    search.insert(c.node);
  }
  Rng r(seed);
  const auto path = random_path(nodes, fixed, s, r);
  const double sill = v.sill();
  std::vector<double> w;
  for (std::size_t x : path) {
    const double mx = local_mean.empty() ? 0.0 : local_mean[x];
    double mean = mx, var = sill;
    const auto nb = search.find(x);
    if (!nb.empty()) {
      double kv = 0;
      if (krige(nodes, nb, x, v, s.ordinary, w, kv)) {
        mean = s.ordinary ? 0.0 : mx;
        for (std::size_t i = 0; i < nb.size(); ++i) mean += w[i] * (s.ordinary ? y[nb[i]] : y[nb[i]] - (local_mean.empty() ? 0.0 : local_mean[nb[i]]));
        var = std::max(0.0, kv);
      }
    }
    y[x] = mean + std::sqrt(var) * r.normal();
    search.insert(x);
  }
  return y;
}

// ---- sequential indicator simulation ------------------------------------------------------------
std::vector<int> sis(const std::vector<Vec3>& nodes, const std::vector<std::pair<std::size_t, int>>& data, const std::vector<double>& proportions,
                     const std::vector<Variogram>& v, const SearchOptions& s, std::uint64_t seed) {
  const std::size_t n = nodes.size(), K = proportions.size();
  if (K < 2 || v.size() != K) throw Error(ErrorCode::InvalidArgument, "sis: at least two facies, one variogram each");
  std::vector<int> f(n, -1);
  std::vector<char> fixed(n, 0);
  Search search(nodes, s);
  for (auto& [node, code] : data) {
    if (node >= n || code < 0 || static_cast<std::size_t>(code) >= K) throw Error(ErrorCode::InvalidArgument, "sis: conditioning node or code out of range");
    if (fixed[node]) continue;
    f[node] = code;
    fixed[node] = 1;
    search.insert(node);
  }
  Rng r(seed);
  const auto path = random_path(nodes, fixed, s, r);
  std::vector<double> w, p(K), count(K, 0.0);
  double done = 0;
  for (std::size_t x : path) {
    const auto nb = search.find(x);
    double tot = 0;
    for (std::size_t k = 0; k < K; ++k) {
      double pk = proportions[k], kv = 0;
      if (!nb.empty() && krige(nodes, nb, x, v[k], false, w, kv)) {
        for (std::size_t i = 0; i < nb.size(); ++i) pk += w[i] * ((f[nb[i]] == static_cast<int>(k) ? 1.0 : 0.0) - proportions[k]);
      }
      // servosystem (Strebelle 2002, factor 1/2): pulls the running proportions to the targets
      if (done > 0) pk += proportions[k] - count[k] / done;
      p[k] = std::clamp(pk, 0.0, 1.0);  // order relations: clip, then renormalise
      tot += p[k];
    }
    if (!(tot > 0)) {
      p = proportions;
      tot = std::accumulate(p.begin(), p.end(), 0.0);
    }
    const double u = r.uniform() * tot;
    double cum = 0;
    int code = static_cast<int>(K) - 1;
    for (std::size_t k = 0; k < K; ++k) {
      cum += p[k];
      if (u < cum) { code = static_cast<int>(k); break; }
    }
    f[x] = code;
    count[static_cast<std::size_t>(code)] += 1;
    done += 1;
    search.insert(x);
  }
  return f;
}

// ---- object-based facies ----------------------------------------------------------------------------
ObjectResult object_facies(const std::vector<Vec3>& nodes, const Vec3& lo, const Vec3& hi, const ObjectSpec& o,
                           const std::vector<std::pair<std::size_t, int>>& wells, std::uint64_t seed, std::size_t max_iterations) {
  if (o.kind != "channel" && o.kind != "lobe") throw Error(ErrorCode::InvalidArgument, "object facies: kind is channel or lobe");
  if (!(o.thickness > 0) || (o.kind == "channel" && !(o.width > 0)) || (o.kind == "lobe" && !(o.radius > 0)))
    throw Error(ErrorCode::InvalidArgument, "object facies: positive thickness and width (channel) or radius (lobe)");
  if (!(o.proportion > 0 && o.proportion < 1)) throw Error(ErrorCode::InvalidArgument, "object facies: the proportion is between 0 and 1");
  const std::size_t n = nodes.size();
  ObjectResult res;
  res.facies.assign(n, 0);
  std::vector<int> want(n, -1);
  for (auto& [node, code] : wells)
    if (node < n) want[node] = code == o.code ? 1 : 0;
  res.well_cells = 0;
  for (int w : want) res.well_cells += w >= 0;
  Rng r(seed);
  const double az = o.azimuth * kPiG / 180.0, sa = std::sin(az), ca = std::cos(az);
  const double amp = o.kind == "channel" && o.sinuosity > 1 && o.wavelength > 0 ? o.wavelength / kPiG * std::sqrt(o.sinuosity - 1) : 0.0;
  struct Obj { Vec3 c; double phase = 0; };
  // a node is inside an object: channel - lens of width w and thickness t about a sinuous
  // centreline through c along the azimuth; lobe - lens of radius R and thickness t about c
  auto inside = [&](const Obj& ob, const Vec3& p) {
    const double dx = p.x - ob.c.x, dy = p.y - ob.c.y;
    double lateral;
    if (o.kind == "channel") {
      const double along = dx * sa + dy * ca, across = dx * ca - dy * sa;
      const double centre = amp > 0 ? amp * std::sin(2 * kPiG * along / o.wavelength + ob.phase) : 0.0;
      lateral = std::fabs(across - centre) / (0.5 * o.width);
    } else {
      lateral = std::sqrt(dx * dx + dy * dy) / o.radius;
    }
    if (lateral >= 1) return false;
    const double t = o.thickness * (1 - lateral * lateral);  // lens: thickest at the axis
    return p.z <= ob.c.z && p.z >= ob.c.z - t;            // c.z is the object's top
  };
  auto covers_zero_well = [&](const Obj& ob) {
    for (std::size_t i = 0; i < n; ++i)
      if (want[i] == 0 && inside(ob, nodes[i])) return true;
    return false;
  };
  auto place = [&](const Obj& ob) {
    std::size_t added = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (res.facies[i] != o.code && inside(ob, nodes[i])) {
        res.facies[i] = o.code;
        ++added;
      }
    ++res.objects;
    return added;
  };
  std::size_t coded = 0, it = 0;
  const std::size_t target = static_cast<std::size_t>(std::llround(o.proportion * static_cast<double>(n)));
  // 1. objects through the wells' object cells first (conditioning by placement)
  for (std::size_t i = 0; i < n && it < max_iterations; ++i) {
    if (want[i] != 1 || res.facies[i] == o.code) continue;
    for (std::size_t tries = 0; tries < 50 && it < max_iterations; ++tries, ++it) {
      // the object's reference point moved along its axis (channel) or anywhere within half its
      // radius (lobe); a channel's phase puts its centreline through the well
      Obj ob;
      const double shift = (o.kind == "channel" ? o.wavelength > 0 ? o.wavelength : 4 * o.width : 0.5 * o.radius) * (2 * r.uniform() - 1);
      const double side = o.kind == "lobe" ? 0.5 * o.radius * (2 * r.uniform() - 1) : 0.0;
      ob.c = {nodes[i].x - shift * sa - side * ca, nodes[i].y - shift * ca + side * sa, nodes[i].z + 0.3 * o.thickness * r.uniform()};
      ob.phase = amp > 0 ? -2 * kPiG * shift / o.wavelength : 0.0;
      if (!inside(ob, nodes[i]) || covers_zero_well(ob)) continue;
      coded += place(ob);
      break;
    }
  }
  // 2. objects anywhere until the target proportion, rejecting those that cover a background well cell
  while (coded < target && it < max_iterations) {
    ++it;
    Obj ob;
    ob.phase = 2 * kPiG * r.uniform();
    ob.c = {lo.x + (hi.x - lo.x) * r.uniform(), lo.y + (hi.y - lo.y) * r.uniform(), lo.z + (hi.z - lo.z + o.thickness) * r.uniform()};
    if (covers_zero_well(ob)) continue;
    // stop before overshooting by more than half an object
    std::size_t would = 0;
    for (std::size_t i = 0; i < n; ++i) would += res.facies[i] != o.code && inside(ob, nodes[i]);
    if (would == 0) continue;
    if (coded + would > target && coded + would - target > target - coded) break;
    coded += place(ob);
  }
  res.proportion = n ? static_cast<double>(coded) / static_cast<double>(n) : 0.0;
  res.well_misfit = 0;
  for (std::size_t i = 0; i < n; ++i)
    if (want[i] >= 0 && (want[i] == 1) != (res.facies[i] == o.code)) ++res.well_misfit;
  return res;
}

}  // namespace stratonix::geo
