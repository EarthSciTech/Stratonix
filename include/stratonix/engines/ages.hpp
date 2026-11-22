/// @file ages.hpp
/// @brief Stage 4C (STX-TN-16 section 1.4): age data in wells and correlation by age.
#pragma once
#include <string>
#include <vector>

namespace stratonix {

struct AgeDatum {
  std::string kind;         ///< "biostratigraphic" or "chronostratigraphic"
  std::string name;         ///< event, zone or sample
  double depth = 0;         ///< metres (TVD below the datum), or the top of a range
  double depth_base = -1;   ///< base of a depth range, -1 for a point
  double age = 0, sigma = 0;  ///< Ma
  std::string timescale, source;
};
enum class AgeRule { LinearInDepth, ConstantRatePerUnit };
struct WellAges { std::string well; std::vector<AgeDatum> data; };
struct AgeConflict { std::string well, what; double depth = 0; };
/// Age at a depth between the dated samples of one well (NaN outside the dated interval).
double age_at_depth(const std::vector<AgeDatum>& data, double depth, AgeRule rule);
/// Depth of an age in one well (NaN when the age lies outside the dated interval).
double depth_of_age(const std::vector<AgeDatum>& data, double age, AgeRule rule);
struct Correlation {
  std::vector<std::string> horizons, wells;
  std::vector<std::vector<double>> depth;   ///< [horizon][well], NaN where it cannot be placed
  std::vector<AgeConflict> conflicts;
};
Correlation correlate_by_age(const std::vector<WellAges>& wells, const std::vector<std::pair<std::string, double>>& horizon_ages,
                             AgeRule rule);

}  // namespace stratonix
