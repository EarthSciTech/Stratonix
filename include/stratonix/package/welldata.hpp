/// @file welldata.hpp
/// @brief Stage 4C (STX-TN-16 section 1): the well data of a build - logs on the well paths,
/// correlation by age, logs upscaled to grid cells, the C35/C38 results and the panels.
#pragma once
#include <map>
#include <string>
#include <vector>

#include "stratonix/checks/checks.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/figures/figures.hpp"
#include "stratonix/repr/grid.hpp"

namespace stratonix {

struct WellDataBuild {
  bool present = false;
  CheckContext::WellDataResults results;
  std::vector<fig::WellPanel> panels;        ///< one per log file (index MD)
  bool correlation = false;
  fig::CorrelationPanel corr;
  std::vector<repr::CellProperty> upscaled;  ///< grid properties from logs
  std::string report;                        ///< well_data_report.txt
  std::map<std::string, std::string> csv;    ///< file name -> content
};
/// Reads the well files of the description (as delivered) and assembles the well data; `grid` may
/// be null (no upscaling). Throws on a file that cannot be read or a curve that is not there.
WellDataBuild well_data(const desc::Description& d, const Model& m, const repr::CornerPointGrid* grid);

}  // namespace stratonix
