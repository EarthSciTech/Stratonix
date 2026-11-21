/// @file posterior.hpp
/// @brief Stage 4A (STX-TN-14 section 4): Bayesian conditioning. Priors from the description's
/// distributions or ranges, Gaussian likelihood of the data residuals, an affine-invariant
/// ensemble sampler (stretch move), split R-hat and effective sample size.
#pragma once
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"

namespace stratonix {

struct PosteriorOptions {
  int walkers = 0;            ///< 0 = max(4 * parameters, 16), even
  int steps = 2000;           ///< per walker, after burn-in
  int burn_in = 500;
  int thin = 10;              ///< members written: every thin-th step of every walker
  double stretch = 2.0;       ///< a of the stretch move
  std::uint64_t seed = 0;
  int threads = 0;            ///< 0 = hardware concurrency
};

/// Prior of one parameter (SI units).
struct PriorSpec {
  std::string path, distribution = "uniform";  ///< uniform | normal | triangular
  double a = 0, b = 0, c = 0;                  ///< uniform: min, max; normal: mean, sd; triangular: min, mode, max
  double min = 0, max = 0;                     ///< admissible range (hard)
  double log_density(double x) const;
};

struct PosteriorResult {
  bool ok = false;
  std::string reason;
  std::vector<std::string> paths;
  std::vector<std::vector<double>> chain;   ///< [walker * steps + step][parameter] after burn-in, SI
  std::vector<double> log_posterior;        ///< same order as chain
  std::vector<double> mean, sd, q025, q975;  ///< per parameter, SI
  std::vector<double> rhat, ess;            ///< split R-hat and effective sample size per parameter
  double acceptance = 0;
  bool converged = false;                   ///< every rhat < 1.01 and ess >= 100
  long builds = 0;                          ///< model evaluations used
  std::uint64_t seed = 0;
  int walkers = 0, steps = 0, burn_in = 0, thin = 1;
};

/// A log-posterior over a parameter vector (SI): NaN or -inf rejects the point.
using LogDensity = std::function<double(const std::vector<double>&)>;

/// Sample a log density with the affine-invariant ensemble sampler, walkers started from `start`
/// (one vector per walker). Deterministic for a given seed whatever the number of threads.
PosteriorResult sample_ensemble(const LogDensity& f, const std::vector<std::vector<double>>& start, const PosteriorOptions& o);

/// Split R-hat (Vehtari et al. 2021, without rank normalisation) of one parameter over chains
/// [chain][draw].
double split_rhat(const std::vector<std::vector<double>>& chains);
/// Effective sample size from the integrated autocorrelation time of the walker-mean series
/// (Goodman and Weare 2010; Sokal's automatic window).
double effective_sample_size(const std::vector<std::vector<double>>& chains);

/// Bayesian conditioning of a description (conditioning.method: bayesian): priors, likelihood of
/// the data residuals, sampler started around the least-squares fit.
PosteriorResult sample_posterior(const desc::Description& d, const BuildInputs& in, const PosteriorOptions& o, BuildLog& log);

}  // namespace stratonix
