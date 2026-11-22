// Stage 4C stubs (STX-TN-16 sections 1.3 and 1.4): implemented by Prompt 55.
#include "stratonix/engines/ages.hpp"
#include "stratonix/engines/surveys.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

std::vector<SurveyPoint> survey_positions(const std::vector<SurveyStation>&, const SurveySpec&) {
  throw Error(ErrorCode::NotImplemented, "survey_positions: not implemented yet (Prompt 55)");
}
std::string survey_method_name(SurveyMethod m) {
  switch (m) {
    case SurveyMethod::MinimumCurvature: return "minimum curvature";
    case SurveyMethod::RadiusOfCurvature: return "radius of curvature";
    case SurveyMethod::BalancedTangential: return "balanced tangential";
    case SurveyMethod::AverageAngle: return "average angle";
    case SurveyMethod::Tangential: return "tangential";
  }
  return "";
}
double age_at_depth(const std::vector<AgeDatum>&, double, AgeRule) { throw Error(ErrorCode::NotImplemented, "age_at_depth: not implemented yet (Prompt 55)"); }
double depth_of_age(const std::vector<AgeDatum>&, double, AgeRule) { throw Error(ErrorCode::NotImplemented, "depth_of_age: not implemented yet (Prompt 55)"); }
Correlation correlate_by_age(const std::vector<WellAges>&, const std::vector<std::pair<std::string, double>>&, AgeRule) {
  throw Error(ErrorCode::NotImplemented, "correlate_by_age: not implemented yet (Prompt 55)");
}

}  // namespace stratonix
