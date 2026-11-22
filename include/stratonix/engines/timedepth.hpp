/// @file timedepth.hpp
/// @brief Stage 4C (STX-TN-16 section 2.3, decision 23): layer-cake time-to-depth conversion with
/// an instantaneous velocity linear in depth per layer, calibrated to wells.
#pragma once
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
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
  std::vector<std::vector<double>> residual;         ///< [layer][well], model minus top (m), with the well's own v0 (ties)
  std::vector<std::vector<double>> loo_residual;     ///< [layer][well], with v0 from the other wells (inverse distance squared)
  std::vector<std::string> excluded;                 ///< wells excluded by the user, with reasons
};
/// k by weighted least squares per layer (user weights only), v0 solved per well so the layer ties.
/// `fixed_k`: per layer, a stated k (NaN or absent: fitted). Tops at NaN are not used.
Calibration calibrate(const std::vector<WellTie>& wells, std::size_t layers, const std::vector<std::string>& exclude = {},
                      const std::vector<double>& fixed_k = {});

/// Stage 4C (Prompt 57): depth spread of one layer base from velocity and pick uncertainty.
struct VelocityUncertainty { double v0_sigma = 0, k_sigma = 0, pick_sigma_ms = 0; };
struct DepthSpread { double total_sd = 0, velocity_sd = 0, pick_sd = 0; std::size_t members = 0; };
DepthSpread depth_spread(const LayerVelocity& v, const VelocityUncertainty& u, double zt, double twt_layer, std::size_t members,
                         std::uint64_t seed);

/// Section 3: one realisation of a spatially correlated error on a regular lattice (nx by ny nodes
/// at spacings dx, dy): zero mean, standard deviation `sigma`, covariance sigma^2 exp(-(r/L)^2),
/// made as a moving average of white noise with a Gaussian kernel (Oliver 1995). L not positive:
/// one value for the whole lattice (fully correlated); L below the spacing: independent nodes.
/// Portable stream (SplitMix64 and Box-Muller): identical on every platform for a seed.
std::vector<double> correlated_field(std::size_t nx, std::size_t ny, double dx, double dy, double sigma, double L, std::uint64_t seed);
/// Conditions a zero-mean field to zero at points by simple kriging of its values there with the
/// field's covariance (conditional simulation, Journel and Huijbregts 1978):
/// f(x) - c(x)^T C^-1 f(x_i). `pts` are lattice coordinates (i, j as real numbers).
void condition_field(std::vector<double>& f, std::size_t nx, std::size_t ny, double dx, double dy, double L,
                     const std::vector<std::pair<double, double>>& pts);
/// A standard normal value from a seed and an index (portable).
double seeded_normal(std::uint64_t seed, std::uint64_t index);

}  // namespace stratonix
