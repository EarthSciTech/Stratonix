/// @file surveys.hpp
/// @brief Stage 4C (STX-TN-16 section 1.3): deviation surveys by a stated calculation method.
#pragma once
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"

namespace stratonix {

enum class SurveyMethod { MinimumCurvature, RadiusOfCurvature, BalancedTangential, AverageAngle, Tangential };
enum class AzimuthReference { Unstated, TrueNorth, GridNorth, MagneticNorth };

struct SurveyStation { double md = 0, inc = 0, azi = 0; };  ///< metres, radians, radians
struct SurveySpec {
  SurveyMethod method = SurveyMethod::MinimumCurvature;
  AzimuthReference reference = AzimuthReference::Unstated;
  double declination = 0;   ///< radians, east positive (magnetic north)
  double convergence = 0;   ///< radians: grid north minus true north at the well
};
struct SurveyPoint { double md = 0, north = 0, east = 0, tvd = 0, dls = 0; };  ///< dls in rad/m
/// Positions relative to the head (north and east in the grid frame, TVD positive down).
/// Refused (Error) when the azimuth reference is unstated.
std::vector<SurveyPoint> survey_positions(const std::vector<SurveyStation>& stations, const SurveySpec& spec);
std::string survey_method_name(SurveyMethod m);

}  // namespace stratonix
