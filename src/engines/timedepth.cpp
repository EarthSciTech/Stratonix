// Stage 4C stubs (STX-TN-16 sections 2.3 and 3): implemented by Prompts 56 and 57.
#include "stratonix/engines/timedepth.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

double layer_thickness(const LayerVelocity&, double, double) { throw Error(ErrorCode::NotImplemented, "layer_thickness: not implemented yet (Prompt 56)"); }
double layer_time(const LayerVelocity&, double, double) { throw Error(ErrorCode::NotImplemented, "layer_time: not implemented yet (Prompt 56)"); }
std::vector<double> base_depths(const VelocityModel&, double, double) { throw Error(ErrorCode::NotImplemented, "base_depths: not implemented yet (Prompt 56)"); }
double time_to_depth(const VelocityModel&, double, double, double) { throw Error(ErrorCode::NotImplemented, "time_to_depth: not implemented yet (Prompt 56)"); }
double depth_to_time(const VelocityModel&, double, double, double) { throw Error(ErrorCode::NotImplemented, "depth_to_time: not implemented yet (Prompt 56)"); }
Calibration calibrate(const std::vector<WellTie>&, std::size_t, const std::vector<std::string>&) {
  throw Error(ErrorCode::NotImplemented, "calibrate: not implemented yet (Prompt 56)");
}
DepthSpread depth_spread(const LayerVelocity&, const VelocityUncertainty&, double, double, std::size_t, std::uint64_t) {
  throw Error(ErrorCode::NotImplemented, "depth_spread: not implemented yet (Prompt 57)");
}

}  // namespace stratonix
