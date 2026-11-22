// Stage 4A (STX-TN-14 section 4, decision D-046): Bayesian conditioning. Priors from the
// description's uncertainty block or the admissible ranges, a Gaussian likelihood of the
// normalised data residuals, and the affine-invariant ensemble sampler of Goodman and Weare
// (2010) with stretch moves and differential-evolution moves (ter Braak 2006; Nelson et al.
// 2014), updated in two halves (Foreman-Mackey et al. 2013).
#include "stratonix/engines/posterior.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <thread>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/conditioning.hpp"

namespace stratonix {

namespace {
constexpr double kNegInf = -std::numeric_limits<double>::infinity();

// Engines do not use the package random generator: xoshiro256** seeded by SplitMix64, one stream
// per walker, so every walker's draws are fixed by (seed, walker) alone.
struct Xoshiro {
  std::uint64_t s[4];
  explicit Xoshiro(std::uint64_t seed) {
    for (auto& v : s) {
      seed += 0x9e3779b97f4a7c15ULL;
      std::uint64_t z = seed;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      v = z ^ (z >> 31);
    }
  }
  static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  std::uint64_t next() {
    const std::uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return r;
  }
  double uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }  // [0, 1)
  std::size_t below(std::size_t n) { return static_cast<std::size_t>(uniform() * static_cast<double>(n)); }
};

double quantile(std::vector<double> v, double q) {
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(v.begin(), v.end());
  double pos = q * static_cast<double>(v.size() - 1);
  std::size_t i = static_cast<std::size_t>(pos);
  double f = pos - static_cast<double>(i);
  return i + 1 < v.size() ? v[i] * (1 - f) + v[i + 1] * f : v[i];
}

// Run fn(i) for i in [0, n) on `threads` threads, static blocks: each index runs exactly once.
template <class F>
void run_blocks(std::size_t n, unsigned threads, F&& fn) {
  threads = std::max(1u, std::min<unsigned>(threads, static_cast<unsigned>(n)));
  if (threads <= 1) {
    for (std::size_t i = 0; i < n; ++i) fn(i);
    return;
  }
  std::exception_ptr err;
  std::atomic<bool> failed{false};
  std::vector<std::thread> ts;
  for (unsigned t = 0; t < threads; ++t)
    ts.emplace_back([&, t] {
      for (std::size_t i = t; i < n; i += threads) {
        try {
          fn(i);
        } catch (...) {
          if (!failed.exchange(true)) err = std::current_exception();
        }
      }
    });
  for (auto& th : ts) th.join();
  if (err) std::rethrow_exception(err);
}
}  // namespace

double PriorSpec::log_density(double x) const {
  if (!std::isfinite(x) || x < min || x > max) return kNegInf;
  if (distribution == "normal") {
    if (!(b > 0)) return kNegInf;
    double z = (x - a) / b;
    return -0.5 * z * z - std::log(b * std::sqrt(2 * kPi));
  }
  if (distribution == "triangular") {
    if (x < a || x > c || !(c > a)) return kNegInf;
    double p = x <= b ? (b > a ? 2 * (x - a) / ((c - a) * (b - a)) : 2 / (c - a)) : (c > b ? 2 * (c - x) / ((c - a) * (c - b)) : 2 / (c - a));
    return p > 0 ? std::log(p) : kNegInf;
  }
  // uniform
  if (x < a || x > b || !(b > a)) return kNegInf;
  return -std::log(b - a);
}

double split_rhat(const std::vector<std::vector<double>>& chains) {
  // split each chain into halves, then the potential scale reduction (Gelman et al. 2013, 11.4)
  std::vector<std::vector<double>> h;
  for (auto& c : chains) {
    std::size_t n = c.size() / 2;
    if (n < 2) continue;
    h.emplace_back(c.begin(), c.begin() + static_cast<std::ptrdiff_t>(n));
    h.emplace_back(c.end() - static_cast<std::ptrdiff_t>(n), c.end());
  }
  if (h.size() < 2) return std::numeric_limits<double>::quiet_NaN();
  const double m = static_cast<double>(h.size());
  std::size_t nn = h[0].size();
  for (auto& c : h) nn = std::min(nn, c.size());
  const double n = static_cast<double>(nn);
  std::vector<double> mean(h.size()), var(h.size());
  for (std::size_t j = 0; j < h.size(); ++j) {
    double s = 0;
    for (std::size_t i = 0; i < nn; ++i) s += h[j][i];
    mean[j] = s / n;
    double v = 0;
    for (std::size_t i = 0; i < nn; ++i) v += (h[j][i] - mean[j]) * (h[j][i] - mean[j]);
    var[j] = v / (n - 1);
  }
  double gm = std::accumulate(mean.begin(), mean.end(), 0.0) / m, B = 0, W = 0;
  for (std::size_t j = 0; j < h.size(); ++j) {
    B += (mean[j] - gm) * (mean[j] - gm);
    W += var[j];
  }
  B *= n / (m - 1);
  W /= m;
  if (!(W > 0)) return B > 0 ? std::numeric_limits<double>::infinity() : 1.0;
  double vplus = (n - 1) / n * W + B / n;
  return std::sqrt(vplus / W);
}

double effective_sample_size(const std::vector<std::vector<double>>& chains) {
  // integrated autocorrelation time from the walker-averaged autocorrelation function, with
  // Sokal's automatic window (smallest M >= 5 tau(M)); ESS = draws / tau
  if (chains.empty()) return 0;
  std::size_t n = chains[0].size();
  for (auto& c : chains) n = std::min(n, c.size());
  if (n < 4) return 0;
  const std::size_t maxlag = n / 2;
  std::vector<double> acf(maxlag + 1, 0.0);
  std::size_t used = 0;
  for (auto& c : chains) {
    double mu = 0;
    for (std::size_t i = 0; i < n; ++i) mu += c[i];
    mu /= static_cast<double>(n);
    double c0 = 0;
    for (std::size_t i = 0; i < n; ++i) c0 += (c[i] - mu) * (c[i] - mu);
    if (!(c0 > 0)) continue;
    ++used;
    for (std::size_t t = 0; t <= maxlag; ++t) {
      double s = 0;
      for (std::size_t i = 0; i + t < n; ++i) s += (c[i] - mu) * (c[i + t] - mu);
      acf[t] += s / c0;
    }
  }
  const double total = static_cast<double>(n) * static_cast<double>(chains.size());
  if (used == 0) return total;
  for (auto& v : acf) v /= static_cast<double>(used);
  double tau = 1;
  for (std::size_t M = 1; M <= maxlag; ++M) {
    tau += 2 * acf[M];
    if (static_cast<double>(M) >= 5 * tau) break;
  }
  tau = std::max(tau, 1.0);
  return total / tau;
}

PosteriorResult sample_ensemble(const LogDensity& f, const std::vector<std::vector<double>>& start, const PosteriorOptions& o) {
  PosteriorResult R;
  const std::size_t W = start.size();
  if (W < 4 || W % 2) { R.reason = "the ensemble needs an even number of walkers, at least 4 (got " + std::to_string(W) + ")"; return R; }
  const std::size_t np = start[0].size();
  for (auto& s : start)
    if (s.size() != np || np == 0) { R.reason = "every walker needs the same, non-empty parameter vector"; return R; }
  if (W < 2 * np) { R.reason = "the ensemble needs at least twice as many walkers as parameters"; return R; }
  if (o.steps < 2 || o.burn_in < 0 || o.thin < 1 || !(o.stretch > 1) || !(o.de_fraction >= 0 && o.de_fraction <= 1)) {
    R.reason = "invalid sampler options (steps >= 2, burn_in >= 0, thin >= 1, stretch > 1, de_fraction in [0, 1])";
    return R;
  }
  const unsigned threads = o.threads > 0 ? static_cast<unsigned>(o.threads) : thread_count();
  const double a = o.stretch;
  std::vector<Xoshiro> rng;
  rng.reserve(W);
  for (std::size_t k = 0; k < W; ++k) rng.emplace_back(o.seed * 0x100000001b3ULL + 0x632be59bd9b4e019ULL * (k + 1));
  std::vector<std::vector<double>> x = start;
  std::vector<double> lp(W);
  std::atomic<long> evals{0};
  run_blocks(W, threads, [&](std::size_t k) {
    double v = f(x[k]);
    ++evals;
    lp[k] = std::isfinite(v) ? v : kNegInf;
  });
  for (std::size_t k = 0; k < W; ++k)
    if (!std::isfinite(lp[k])) { R.reason = "walker " + std::to_string(k) + " starts where the posterior density is zero"; R.builds = evals; return R; }
  const std::size_t half = W / 2, total = static_cast<std::size_t>(o.burn_in + o.steps);
  R.chain.resize(W * static_cast<std::size_t>(o.steps));
  R.log_posterior.resize(W * static_cast<std::size_t>(o.steps));
  std::vector<long> accepted(W, 0);
  for (std::size_t it = 0; it < total; ++it) {
    for (int part = 0; part < 2; ++part) {
      const std::size_t off = part ? half : 0, other = part ? 0 : half;
      // proposals depend only on the complementary half, fixed during this update
      run_blocks(half, threads, [&](std::size_t i) {
        const std::size_t k = off + i;
        auto& g = rng[k];
        std::vector<double> y(np);
        double logq = 0;  // log proposal ratio
        if (g.uniform() < o.de_fraction) {
          // differential evolution: y = x_k + gamma (x_j1 - x_j2) + small noise; symmetric.
          // gamma = 2.38 / sqrt(2 d), and 1 in one move of ten for jumps between modes
          std::size_t j1 = other + g.below(half), j2 = other + g.below(half - 1);
          if (j2 >= j1) ++j2;
          const double gamma = g.uniform() < 0.1 ? 1.0 : 2.38 / std::sqrt(2.0 * static_cast<double>(np));
          for (std::size_t p = 0; p < np; ++p) {
            const double dlt = x[j1][p] - x[j2][p];
            y[p] = x[k][p] + gamma * dlt * (1 + 1e-4 * (2 * g.uniform() - 1));
          }
        } else {
          const std::size_t j = other + g.below(half);
          const double z = std::pow((a - 1) * g.uniform() + 1, 2) / a;
          for (std::size_t p = 0; p < np; ++p) y[p] = x[j][p] + z * (x[k][p] - x[j][p]);
          logq = static_cast<double>(np - 1) * std::log(z);
        }
        double ly = f(y);
        ++evals;
        if (!std::isfinite(ly)) ly = kNegInf;
        const double logr = logq + ly - lp[k];
        const double v = g.uniform();
        if (ly > kNegInf && std::log(std::max(v, 1e-300)) < logr) {
          x[k] = std::move(y);
          lp[k] = ly;
          if (it >= static_cast<std::size_t>(o.burn_in)) ++accepted[k];
        }
      });
    }
    if (it >= static_cast<std::size_t>(o.burn_in)) {
      const std::size_t s = it - static_cast<std::size_t>(o.burn_in);
      for (std::size_t k = 0; k < W; ++k) {
        R.chain[k * static_cast<std::size_t>(o.steps) + s] = x[k];
        R.log_posterior[k * static_cast<std::size_t>(o.steps) + s] = lp[k];
      }
    }
  }
  R.builds = evals;
  R.acceptance = static_cast<double>(std::accumulate(accepted.begin(), accepted.end(), 0L)) / static_cast<double>(W * static_cast<std::size_t>(o.steps));
  R.seed = o.seed;
  R.walkers = static_cast<int>(W);
  R.steps = o.steps;
  R.burn_in = o.burn_in;
  R.thin = o.thin;
  R.converged = true;
  for (std::size_t p = 0; p < np; ++p) {
    std::vector<std::vector<double>> ch(W, std::vector<double>(static_cast<std::size_t>(o.steps)));
    std::vector<double> all;
    all.reserve(R.chain.size());
    for (std::size_t k = 0; k < W; ++k)
      for (std::size_t s = 0; s < static_cast<std::size_t>(o.steps); ++s) {
        ch[k][s] = R.chain[k * static_cast<std::size_t>(o.steps) + s][p];
        all.push_back(ch[k][s]);
      }
    double mu = std::accumulate(all.begin(), all.end(), 0.0) / static_cast<double>(all.size()), v = 0;
    for (double q : all) v += (q - mu) * (q - mu);
    R.mean.push_back(mu);
    R.sd.push_back(std::sqrt(v / static_cast<double>(all.size() - 1)));
    R.q025.push_back(quantile(all, 0.025));
    R.q975.push_back(quantile(all, 0.975));
    R.rhat.push_back(split_rhat(ch));
    R.ess.push_back(effective_sample_size(ch));
    R.converged = R.converged && R.rhat.back() < 1.01 && R.ess.back() >= 100;
  }
  R.ok = true;
  return R;
}

PosteriorOptions posterior_options(const desc::Description& d) {
  PosteriorOptions o;
  const auto& c = d.conditioning;
  if (c.walkers.present) o.walkers = static_cast<int>(std::llround(c.walkers.si()));
  if (c.steps.present) o.steps = static_cast<int>(std::llround(c.steps.si()));
  if (c.burn_in.present) o.burn_in = static_cast<int>(std::llround(c.burn_in.si()));
  if (c.thin.present) o.thin = static_cast<int>(std::llround(c.thin.si()));
  if (c.seed.present) o.seed = static_cast<std::uint64_t>(std::llround(c.seed.si()));
  return o;
}

PosteriorResult sample_posterior(const desc::Description& d, const BuildInputs& in, const PosteriorOptions& opts, BuildLog& log) {
  PosteriorResult R;
  ConditioningProblem pb(d, in);
  if (!pb.ok()) { R.reason = pb.reason(); return R; }
  const std::size_t np = pb.parameters();
  // priors: a distribution of the uncertainty block for the same path, else uniform on the range
  std::vector<PriorSpec> pri(np);
  for (std::size_t k = 0; k < np; ++k) {
    auto& P = pri[k];
    P.path = pb.paths()[k];
    P.min = pb.lower()[k];
    P.max = pb.upper()[k];
    P.a = P.min;
    P.b = P.max;
    for (auto& u : d.uncertainty.parameters) {
      if (u.path != P.path) continue;
      if (u.distribution == "normal" && u.mean.present && u.sd.present) {
        P.distribution = "normal";
        P.a = u.mean.si();
        P.b = u.sd.si();
      } else if (u.distribution == "triangular" && u.min.present && u.mode.present && u.max.present) {
        P.distribution = "triangular";
        P.a = std::max(P.min, u.min.si());
        P.b = u.mode.si();
        P.c = std::min(P.max, u.max.si());
      } else if (u.distribution == "uniform" && u.min.present && u.max.present) {
        P.a = std::max(P.min, u.min.si());
        P.b = std::min(P.max, u.max.si());
      }
    }
    log.value("posterior.prior." + P.path, P.distribution + " (" + fmt_shortest(P.a / pb.to_si()[k]) + ", " + fmt_shortest(P.b / pb.to_si()[k]) +
                                                (P.distribution == "triangular" ? ", " + fmt_shortest(P.c / pb.to_si()[k]) : std::string()) + ")");
  }
  LogDensity f = [&](const std::vector<double>& x) {
    double l = 0;
    for (std::size_t k = 0; k < np; ++k) l += pri[k].log_density(x[k]);
    if (!std::isfinite(l)) return kNegInf;  // outside the prior: rejected without a build
    std::vector<double> r;
    if (!pb.residuals(x, r)) return kNegInf;
    double c2 = 0;
    for (double v : r) c2 += v * v;
    return l - 0.5 * c2;
  };
  // start: a small ball around the least-squares fit (or the start values when it fails)
  BuildLog lq;
  auto ls = condition_to_data(d, in, lq);
  std::vector<double> centre = pb.start(), scale(np);
  for (std::size_t k = 0; k < np; ++k) scale[k] = 1e-3 * (pb.upper()[k] - pb.lower()[k]);
  if (ls.converged) {
    for (std::size_t k = 0; k < np; ++k) {
      centre[k] = ls.parameters[k].value * pb.to_si()[k];
      if (ls.parameters[k].sigma > 0) scale[k] = 0.1 * ls.parameters[k].sigma * pb.to_si()[k];
    }
    std::string lap;
    for (std::size_t k = 0; k < np; ++k)
      lap += (k ? "; " : "") + ls.parameters[k].path + " " + fmt_fixed(ls.parameters[k].value, 3) + " +- " + fmt_fixed(ls.parameters[k].sigma, 3) + " " + ls.parameters[k].unit;
    log.value("posterior.start", "least-squares fit (" + std::to_string(ls.builds) + " model builds)");
    log.value("posterior.laplace", lap);
  } else {
    log.value("posterior.start", "start values (the least-squares fit failed: " + ls.reason + ")");
  }
  PosteriorOptions o = opts;
  if (o.walkers <= 0) o.walkers = static_cast<int>(std::max<std::size_t>(4 * np, 16));
  if (o.walkers % 2) ++o.walkers;
  Xoshiro g(o.seed ^ 0xa0761d6478bd642fULL);
  std::vector<std::vector<double>> start;
  for (int w = 0; w < o.walkers; ++w) {
    std::vector<double> x(np);
    for (std::size_t k = 0; k < np; ++k) {
      // Box-Muller normal, clamped inside the admissible range
      double u1 = std::max(g.uniform(), 1e-300), u2 = g.uniform();
      double z = std::sqrt(-2 * std::log(u1)) * std::cos(2 * kPi * u2);
      double lo = pb.lower()[k], hi = pb.upper()[k], m = 1e-6 * (hi - lo);
      x[k] = std::clamp(centre[k] + z * scale[k], lo + m, hi - m);
    }
    start.push_back(std::move(x));
  }
  R = sample_ensemble(f, start, o);
  R.paths = pb.paths();
  R.priors = pri;
  if (ls.converged)
    for (std::size_t k = 0; k < np; ++k) {
      R.lsq.push_back(ls.parameters[k].value * pb.to_si()[k]);
      R.lsq_sigma.push_back(ls.parameters[k].sigma * pb.to_si()[k]);
    }
  if (!R.ok) return R;
  // report in SI internally; the log in the parameters' units
  for (std::size_t k = 0; k < np; ++k) {
    const double s = pb.to_si()[k];
    const desc::Param* P = d.param_by_name(R.paths[k]);
    const std::string unit = P ? P->unit : "";
    log.value("posterior." + R.paths[k], fmt_fixed(R.mean[k] / s, 3) + " +- " + fmt_fixed(R.sd[k] / s, 3) + " " + unit + " (95 % " +
                                             fmt_fixed(R.q025[k] / s, 3) + " to " + fmt_fixed(R.q975[k] / s, 3) + "; R-hat " + fmt_fixed(R.rhat[k], 4) +
                                             ", ESS " + fmt_fixed(R.ess[k], 0) + ")");
  }
  log.value("posterior.sampler", std::to_string(R.walkers) + " walkers, " + std::to_string(R.burn_in) + " burn-in + " + std::to_string(R.steps) +
                                     " steps, thin " + std::to_string(R.thin) + ", seed " + std::to_string(R.seed) + ", acceptance " +
                                     fmt_fixed(R.acceptance, 3) + ", " + std::to_string(R.builds) + " posterior evaluations");
  log.value("posterior.likelihood", "independent Gaussian residuals of " + std::to_string(pb.size()) + " data (stated sigma per datum)");
  log.value("posterior.converged", R.converged ? "yes (every R-hat < 1.01, ESS >= 100)" : "no (an R-hat >= 1.01 or an ESS < 100: run longer)");
  return R;
}

}  // namespace stratonix
