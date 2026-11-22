/// @file geostat_props.hpp
/// @brief Stage 4D (STX-TN-17 section 1): the description's geostatistical properties on the
/// corner-point grid - cells mapped to geological-time coordinates, conditioned to upscaled logs,
/// simulated by SGS, SIS or objects; variograms, histograms, conditioning and statistics reported.
#pragma once
#include <string>
#include <vector>

#include "stratonix/desc/description.hpp"
#include "stratonix/engines/geostat.hpp"
#include "stratonix/model/model.hpp"
#include "stratonix/repr/grid.hpp"

namespace stratonix {

struct GeostatBuild {
  bool present = false;
  std::vector<repr::CellProperty> props;
  std::string report, csv;  ///< geostatistics_report.txt, variograms.csv
  struct VariogramPlot {
    std::string property, unit, direction;  ///< direction: "horizontal (major)", "vertical", "wells (vertical)"
    std::vector<double> lag, gamma;
    std::vector<std::size_t> pairs;
    geo::Variogram model;
    double azimuth = 0;  ///< degrees: the model is drawn along this direction
    bool vertical = false;
  };
  std::vector<VariogramPlot> variograms;
  struct Histogram { std::string property, unit, value_unit; std::vector<double> simulated, wells; std::vector<std::string> codes; };
  std::vector<Histogram> histograms;
  struct Conditioning { std::string property; std::size_t cells = 0, mismatched = 0; double largest = 0; bool exact = true; };
  std::vector<Conditioning> conditioning;  ///< C39
  struct Statistic { std::string property, what; double target = 0, got = 0, tolerance = 0; };
  std::vector<Statistic> statistics;       ///< C40
  std::size_t cells = 0, fallback = 0;     ///< simulated cells; those simulated in present-day coordinates
};
/// Simulate every grid.geostatistics entry on the grid; `existing` holds the properties already
/// on the grid (upscaled logs for conditioning, transforms for trends).
GeostatBuild geostat_properties(const desc::Description& d, const Model& m, const repr::CornerPointGrid& g, const repr::PropertySet& existing);

}  // namespace stratonix
