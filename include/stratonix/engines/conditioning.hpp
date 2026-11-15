/// @file conditioning.hpp
/// @brief Stage 3A (STX-TN-11): faults from data (throws and tips from horizon cut-offs, fault
/// profiles from sticks) and hybrid conditioning (parameters of the v2 mechanisms fitted to data
/// by weighted least squares, with uncertainties and the data that control them).
#pragma once
#include <string>
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
  double strike = kNaN, rms = kNaN;
  Vec3 top_point{0, 0, 0};
};
StickProfile profile_from_sticks(const std::vector<Polyline>& sticks);

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
};
/// Fit the parameters named in d.conditioning to the data. Never throws for a failed fit: the
/// result says why (converged = false).
ConditioningResult condition_to_data(const desc::Description& d, const BuildInputs& in, BuildLog& log);

}  // namespace stratonix
