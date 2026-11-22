/// @file velocity.hpp
/// @brief Stage 4C (STX-TN-16 section 2.3, decision 23): the description's velocity model - time
/// surfaces of the layer bases, wells' tops and times, calibration, and the depth-converted surfaces.
#pragma once
#include <map>
#include <string>
#include <vector>

#include "stratonix/core/grid.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/timedepth.hpp"

namespace stratonix::io {

struct VelocityBuild {
  bool present = false;
  VelocityModel model;                         ///< base times and calibrated velocities at any (x, y)
  Calibration cal;
  std::vector<std::string> layers, wells;
  std::map<std::string, Grid2D> time_grids;    ///< by input id (s, two-way)
  std::map<std::string, Grid2D> depth_grids;   ///< by input id: the same nodes in depth (m, positive down)
  struct Tie { std::string well, layer; double residual = 0, loo = 0; bool excluded = false; };
  std::vector<Tie> ties;                       ///< after conversion; excluded wells are blind checks
  std::string report, csv;                     ///< velocity_report.txt, velocity_calibration.csv
};
/// Loads the time surfaces, calibrates to the wells and converts the surfaces to depth.
VelocityBuild velocity_build(const desc::Description& d);
/// Two-way time at a depth; NaN where the velocity model is undefined (outside the time surfaces).
double depth_to_time_safe(const VelocityModel& m, double x, double y, double z);

}  // namespace stratonix::io
