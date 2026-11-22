/// @file timedepth.hpp
/// @brief Stage 4C (STX-TN-16 section 2.3, decision 23): layer-cake time-to-depth conversion with
/// an instantaneous velocity linear in depth per layer, calibrated to wells.
#pragma once
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace stratonix {

/// v(z) = v0 + k z (z positive down from the datum, m; v in m/s; k in 1/s). k = 0: constant velocity.
struct LayerVelocity { double v0 = 2000, k = 0; };
/// Thickness of a layer whose top is at depth zt, crossed in one-way time dt (s).
double layer_thickness(const LayerVelocity& v, double zt, double dt);
/// One-way time to cross from zt to zb.
double layer_time(const LayerVelocity& v, double zt, double zb);

struct VelocityModel {
  /// two-way times (s) of the layer bases at (x, y); size = layers
  std::function<std::vector<double>(double, double)> base_twt;
  /// velocity of each layer at (x, y) (v0 may vary laterally)
  std::function<std::vector<LayerVelocity>(double, double)> velocity;
  std::vector<std::string> layer_names;
};
/// Depths of the layer bases at (x, y) (m, positive down).
std::vector<double> base_depths(const VelocityModel& m, double x, double y);
double time_to_depth(const VelocityModel& m, double x, double y, double twt);
double depth_to_time(const VelocityModel& m, double x, double y, double z);

/// A well's data for calibration: the depths of the layer bases (tops) and their two-way times.
struct WellTie { std::string well; double x = 0, y = 0; std::vector<double> top_depth, top_twt; double weight = 1; };
struct Calibration {
  std::vector<double> k, k_sigma;                    ///< per layer
  std::vector<std::vector<double>> v0;               ///< [layer][well]
  std::vector<std::vector<double>> residual;         ///< [layer][well], model minus top (m)
  std::vector<std::string> excluded;                 ///< wells excluded by the user, with reasons
};
/// k by weighted least squares per layer (user weights only), v0 solved per well so the layer ties.
Calibration calibrate(const std::vector<WellTie>& wells, std::size_t layers, const std::vector<std::string>& exclude = {});

/// Stage 4C (Prompt 57): depth spread of one layer base from velocity and pick uncertainty.
struct VelocityUncertainty { double v0_sigma = 0, k_sigma = 0, pick_sigma_ms = 0; };
struct DepthSpread { double total_sd = 0, velocity_sd = 0, pick_sd = 0; std::size_t members = 0; };
DepthSpread depth_spread(const LayerVelocity& v, const VelocityUncertainty& u, double zt, double twt_layer, std::size_t members,
                         std::uint64_t seed);

}  // namespace stratonix
