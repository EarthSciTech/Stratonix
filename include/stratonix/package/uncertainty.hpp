// Uncertainty (stage 3C, STX-TN-13): distributions on labelled parameters, scenarios, the seeded
// ensemble runner, statistics, sensitivity and probability maps. An ensemble is a set of
// descriptions, each built by the existing pipeline; the model interface is unchanged.
#pragma once
#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "stratonix/desc/description.hpp"
#include "stratonix/desc/document.hpp"

namespace stratonix::unc {

// ---- generator (section 3) ------------------------------------------------------------------

/// SplitMix64 of a seed and a member index (Steele et al. 2014): the member's stream seed.
std::uint64_t member_seed(std::uint64_t seed, std::uint64_t member);

/// xoshiro256** (Blackman and Vigna 2021), seeded through SplitMix64. Integer arithmetic only:
/// identical streams on every platform.
class Rng {
 public:
  explicit Rng(std::uint64_t seed);
  std::uint64_t next();
  /// Uniform in the open interval (0, 1), from the top 53 bits.
  double uniform();
 private:
  std::uint64_t s_[4];
};

// ---- distributions (section 1) --------------------------------------------------------------

/// Standard normal CDF and its inverse (Acklam 2003, one Halley step).
double normal_cdf(double x);
double normal_quantile(double p);

struct Distribution {
  std::string type;              ///< uniform, normal, triangular, discrete; posterior (stage 4A: sorted samples, equal weights)
  double a = 0, b = 0, c = 0;    ///< uniform: min, max; normal: mean, sd; triangular: min, mode, max
  double lo = -1e300, hi = 1e300;  ///< normal truncation (min, max)
  std::vector<double> values, weights;  ///< discrete
  double quantile(double u) const;      ///< inverse CDF
  double mean() const;
  double variance() const;
};

struct UncertainParameter {
  std::string path;     ///< description path, e.g. faults.F1.throw.max
  Distribution dist;    ///< in the base parameter's unit
  std::string unit;     ///< the base parameter's unit
  double base = 0;      ///< the description's value (base unit)
  std::string evidence, source;  ///< of the distribution's parameters
};
struct ScenarioOption {
  std::string name;
  double weight = 0;
  std::vector<std::pair<std::string, desc::Node>> set;  ///< path, replacement node
  std::vector<std::string> remove;                      ///< paths of list items to delete
};
struct Scenario {
  std::string name;
  std::vector<ScenarioOption> options;
};
struct ProbabilityMapSpec {
  std::string kind;      ///< unit or closure
  std::string unit, horizon;
  double depth = 0;      ///< m, positive down (kind unit)
  double spacing = 100;  ///< m
};
struct EnsembleSpec {
  int members = 0;
  std::uint64_t seed = 0;
  std::vector<UncertainParameter> parameters;
  std::vector<std::pair<std::pair<int, int>, double>> correlations;  ///< parameter indices and rho
  std::vector<Scenario> scenarios;
  bool one_at_a_time = false;
  std::vector<ProbabilityMapSpec> maps;
  std::vector<std::string> key_numbers;  ///< extra key numbers (besides the summary table's)
  bool posterior = false;  ///< stage 4A: members drawn from the posterior of conditioning.method: bayesian
};
/// The ensemble specification of a description's uncertainty block (rule F33 already applied).
EnsembleSpec ensemble_spec(const desc::Description& d);

// ---- sampling ---------------------------------------------------------------------------------

struct MemberSample {
  int index = 0;           ///< 1..N (0: the base case)
  std::uint64_t seed = 0;  ///< member_seed(spec.seed, index)
  std::vector<double> values;   ///< per parameter (base unit)
  std::vector<int> options;     ///< per scenario: chosen option index
};
/// Members 0 (base) to N, deterministic for the spec's seed.
std::vector<MemberSample> sample_members(const EnsembleSpec& s);
/// Gaussian-copula draws for one member: uniforms per parameter (correlated through Cholesky).
std::vector<double> correlated_uniforms(const EnsembleSpec& s, Rng& rng);
/// Cholesky factor (lower, row-major) of an n x n correlation matrix; throws if not positive definite.
std::vector<double> cholesky(const std::vector<double>& corr, int n);

/// The member's description document: sampled values written with their unit and provenance,
/// scenario options applied.
desc::Node member_document(const desc::Node& base, const EnsembleSpec& s, const MemberSample& m);
/// Set (or read) a scalar parameter in a description document by path ("value unit | evidence | source").
void set_parameter(desc::Node& doc, const std::string& path, double value, const std::string& unit, const std::string& evidence,
                   const std::string& source);
bool get_parameter(const desc::Node& doc, const std::string& path, double* value, std::string* unit);

// ---- results (section 4) ---------------------------------------------------------------------

struct Stats {
  std::size_t n = 0;
  double mean = 0, sd = 0, p10 = 0, p50 = 0, p90 = 0, min = 0, max = 0;
};
/// Percentiles by linear interpolation (Hyndman and Fan 1996, type 7).
Stats statistics(std::vector<double> v);
/// First-order index of y on x by equal-count binning (sqrt(N) bins).
double first_order_index(const std::vector<double>& x, const std::vector<double>& y);

struct MemberResult {
  MemberSample sample;
  std::string status;     ///< ok, refused, build_failed, checks_failed
  std::string reason;
  std::vector<std::string> failed_checks;
  std::vector<std::pair<std::string, double>> keys;  ///< key numbers
  std::string topology;   ///< stage 4A: "same" as the base case, or the contacts that differ ("" when not built)
};
struct KeyStats { std::string key, unit; Stats stats; std::size_t failed = 0; };
struct OatRow { std::string parameter, key; double low = 0, high = 0, base = 0; };
struct SobolRow { std::string parameter, key; double index = 0; bool distinguishable = false; };
struct ProbabilityMap { ProbabilityMapSpec spec; double x0 = 0, y0 = 0, dx = 0, dy = 0; std::size_t nx = 0, ny = 0; std::vector<double> p; };

struct EnsembleResult {
  EnsembleSpec spec;
  std::vector<MemberResult> members;  ///< index 0 is the base case
  std::vector<KeyStats> stats;
  std::vector<OatRow> oat;
  std::vector<SobolRow> sobol;
  std::vector<ProbabilityMap> maps;
  std::size_t topology_differs = 0;  ///< stage 4A: members whose contact graph differs from the base case's
  std::size_t ok() const;
};

struct EnsembleOptions {
  int members = 0;          ///< 0: from the description
  std::int64_t seed = -1;   ///< < 0: from the description
  unsigned threads = 0;     ///< 0: default
  bool figures = true;
  bool base_build = false;  ///< also write the base case's full build into out_dir/base
};
/// Run the ensemble of a description; writes into out_dir when it is not empty.
EnsembleResult run_ensemble(const std::filesystem::path& description, const std::filesystem::path& out_dir, const EnsembleOptions& o);

/// Numeric value of a summary key number's text ("167.4 m at H4" -> 167.4), and its unit.
bool parse_key_value(const std::string& text, double* value, std::string* unit);

}  // namespace stratonix::unc
