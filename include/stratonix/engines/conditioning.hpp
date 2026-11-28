/// @file conditioning.hpp
/// @brief Stage 3A (STX-TN-11): faults from data (throws and tips from horizon cut-offs, fault
/// profiles from sticks) and hybrid conditioning (parameters of the v2 mechanisms fitted to data
/// by weighted least squares, with uncertainties and the data that control them).
#pragma once
#include <atomic>
#include <memory>
#include <set>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/engines/implicit.hpp"
#include "stratonix/model/model.hpp"

namespace stratonix {

/// An interface point of a horizon (present-day position).
struct HorizonPoint {
  Vec3 p;
  int horizon = -1;
};

/// Throw estimation (STX-TN-11 section 3.1).
struct ThrowOptions {
  double band = 500;              ///< m: data within this distance of the fault
  double station_spacing = 250;   ///< m along strike
  std::string profile = "elliptical";  ///< along-strike shape: elliptical | flat_topped | bell
  double taper = 0.3;             ///< flat_topped taper (fraction)
};
struct ThrowStation {
  double s = 0;          ///< m along strike from the plane's centroid
  double throw_m = 0;    ///< vertical separation of the cut-offs (+ normal)
  std::size_t n_fw = 0, n_hw = 0;
  double z = kNaN;       ///< stage 4A: mean elevation of the two cut-offs
};
struct ThrowEstimate {
  int horizon = -1;
  std::vector<ThrowStation> stations;
  double max_throw = kNaN, centre_s = kNaN, length = kNaN, rms = kNaN;
  bool ok = false;
  std::string reason;    ///< why not ok
};
/// Station throws and the fitted profile of every horizon with data on both sides. s is measured
/// along strike, direction (sin strike, cos strike), from plane.centroid.
std::vector<ThrowEstimate> estimate_throws(const FittedPlane& plane, const std::vector<HorizonPoint>& pts, std::size_t n_horizons,
                                           const ThrowOptions& o);

/// Profile of dip against depth fitted to the sticks of a listric or ramp-flat fault
/// (STX-TN-11 section 3.2): depths (m, positive down) and dips (rad), the strike of the plane fit,
/// the shallowest point of the mean profile and the RMS distance of the stick points.
struct StickProfile {
  std::vector<double> depth, dip;
  std::vector<double> u;           ///< m along the dip direction from the plane's centroid, per depth
  Vec3 centroid{0, 0, 0};
  double strike = kNaN, rms = kNaN;
  Vec3 top_point{0, 0, 0};
};
StickProfile profile_from_sticks(const std::vector<Polyline>& sticks);

/// Displacement field over a fault surface (stage 4A, STX-TN-14 section 3.1):
/// T(s, d) = T_max f_s(2 (s - s0) / Ls) f_d(2 (d - d0) / Ld), fitted to the cut-off separations of
/// every horizon. s along strike from the plane's centroid, d down dip (m) from the centroid.
struct DisplacementObservation { double s = 0, d = 0, throw_m = 0; int horizon = -1; };
struct DisplacementOptions {
  std::string shape_strike = "elliptical", shape_dip = "elliptical";  ///< stage 2A shapes
  double taper = 0.3;
};
struct DisplacementField {
  bool ok = false;
  std::string reason;
  double max_throw = kNaN, s0 = kNaN, d0 = kNaN, length = kNaN, height = kNaN, rms = kNaN;
  std::size_t observations = 0;
  std::string shape_strike = "elliptical", shape_dip = "elliptical";
  double taper = 0.3;
  double at(double s, double d) const;  ///< T(s, d), 0 outside the tips
};
/// Observations from the stage 3A stations of every horizon (station s, the cut-off's down-dip
/// position from its depth and the plane's dip).
std::vector<DisplacementObservation> displacement_observations(const FittedPlane& plane, const std::vector<ThrowEstimate>& est,
                                                              const std::vector<HorizonPoint>& pts);
DisplacementField fit_displacement_field(const std::vector<DisplacementObservation>& obs, const DisplacementOptions& o);

/// Fault relations inferred from data and rules (STX-TN-14 section 3.2).
struct FaultDataSet {
  std::string name;
  FittedPlane plane;
  std::vector<Vec3> points;  ///< stick points and cut-off points
  double tip_extent = kNaN;  ///< m along strike either side of the centroid, from the fitted field
};
struct FaultRelation {
  std::string a, b;          ///< a <kind> b
  std::string kind;          ///< abuts | splay | crosses
  std::string rule;          ///< the rule and the data it rests on, in words
  double distance = kNaN;    ///< m: closest approach of a's data to b's surface
  std::string side;          ///< abuts and splay: hanging_wall | footwall (the side of b where a lies)
};
std::vector<FaultRelation> infer_fault_relations(const std::vector<FaultDataSet>& faults, double relation_distance);

/// Hybrid conditioning (STX-TN-11 section 4).
struct FittedParameter {
  std::string path, unit;
  double start = kNaN, value = kNaN, sigma = kNaN, min = kNaN, max = kNaN;  ///< in the description's unit
  /// The data (index into ConditioningResult::residuals) with the largest sensitivity, and
  /// their share of the parameter's variance.
  std::vector<std::pair<std::size_t, double>> controlling;
};
struct ConditioningResult {
  bool converged = false;
  std::string reason;          ///< why the fit was refused
  int iterations = 0;
  double chi2 = kNaN, chi2_reduced = kNaN;
  std::vector<FittedParameter> parameters;
  std::vector<DataMisfit> residuals;  ///< at the fitted values; sigma and normalised filled
  std::vector<double> leverage;
  desc::Description fitted;           ///< the description with the fitted values (calibrated)
  long builds = 0;                    ///< model evaluations used (stage 4A: Broyden updates)
  int jacobians = 0;                  ///< full finite-difference Jacobians computed
};
/// A conditioning datum: an interface point or orientation (or, from stage 4A, a gravity or
/// magnetic station) with its stated uncertainty.
struct ConditioningDatum {
  Vec3 p, n;
  int horizon = -1;
  bool oriented = false;
  double sigma = 0;     ///< m; degrees for orientations; mGal or nT for stations
  double observed = 0;  ///< stations: the observed anomaly
  std::string input, evidence, kind;  ///< kind: point | orientation | well_top | gravity | magnetics
};

/// The data, parameters and admissible ranges of d.conditioning, and the forward model: normalised
/// residuals at a parameter vector (SI). Thread-safe: each evaluation builds its own model from
/// its own copy of the description (stage 4A, STX-TN-14 section 4).
class ConditioningProblem {
 public:
  ConditioningProblem(const desc::Description& d, const BuildInputs& in);
  ConditioningProblem(const ConditioningProblem&) = delete;
  ConditioningProblem& operator=(const ConditioningProblem&) = delete;
  bool ok() const { return reason_.empty(); }
  const std::string& reason() const { return reason_; }
  std::size_t size() const { return data_.size(); }
  std::size_t parameters() const { return paths_.size(); }
  const std::vector<std::string>& paths() const { return paths_; }
  const std::vector<double>& start() const { return x0_; }
  const std::vector<double>& lower() const { return lo_; }
  const std::vector<double>& upper() const { return hi_; }
  const std::vector<double>& to_si() const { return to_si_; }
  const std::vector<ConditioningDatum>& data() const { return data_; }
  const std::set<std::string>& inputs_used() const { return used_; }
  /// Normalised residuals (value / sigma); raw values when asked. False when the model cannot be
  /// built at x.
  bool residuals(const std::vector<double>& x, std::vector<double>& r, std::vector<double>* raw = nullptr) const;
  long builds() const { return builds_.load(); }
  /// The description used for the builds (conditioning removed, control wells as check wells).
  const desc::Description& base() const { return dd_; }

 private:
  desc::Description dd_;
  const BuildInputs* in_ = nullptr;
  std::vector<ConditioningDatum> data_;
  std::vector<std::string> paths_;
  std::vector<double> x0_, lo_, hi_, to_si_;
  std::set<std::string> used_;
  std::string reason_;
  mutable std::atomic<long> builds_{0};
};

/// Fit the parameters named in d.conditioning to the data. Never throws for a failed fit: the
/// result says why (converged = false).
ConditioningResult condition_to_data(const desc::Description& d, const BuildInputs& in, BuildLog& log);

/// A parametric model built from the fitted description, with its conditioning result:
/// capability("conditioning") -> const ConditioningResult*, "data_misfit" -> the residuals.
class ParametricModel;
struct PosteriorResult;

class ConditionedModel final : public Model {
 public:
  /// `posterior`: the samples of conditioning.method: bayesian (capability "posterior").
  ConditionedModel(std::shared_ptr<ParametricModel> inner, ConditioningResult r, std::shared_ptr<const PosteriorResult> posterior = nullptr);
  const ModelInfo& info() const override;
  void evaluate(std::span<const Vec3> pts, unsigned fields, EvalResult& out) const override;
  void scalar_in_block(std::span<const Vec3> pts, BlockCode block, int horizon, std::span<double> out) const override;
  void evaluate_fault(int fault, std::span<const Vec3> pts, FaultEval& out) const override;
  void evaluate_fault_geometry(int fault, std::span<const Vec3> pts, FaultEval& out) const override;
  void restore(std::span<const Vec3> pts, std::span<Vec3> out) const override;
  const void* capability(std::string_view name) const override;
  const ConditioningResult& result() const { return r_; }

 private:
  std::shared_ptr<ParametricModel> inner_;
  ConditioningResult r_;
  std::shared_ptr<const PosteriorResult> post_;
  ModelInfo info_;
};

/// Condition and build (the parametric builder calls this for a description with conditioning).
/// With conditioning.method: bayesian the posterior is sampled too and the model is the
/// least-squares fit with the posterior attached. Throws Error(Build) with the reason when the
/// fit or the sampler is refused.
ModelPtr build_conditioned_model(const desc::Description& d, const BuildInputs& in, BuildLog& log);

}  // namespace stratonix
