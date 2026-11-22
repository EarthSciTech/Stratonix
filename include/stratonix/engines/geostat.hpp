/// @file geostat.hpp
/// @brief Stage 4D (STX-TN-17 section 1, decision 24): variograms, normal scores, sequential
/// Gaussian and indicator simulation, and object-based facies, on points given in geological-time
/// coordinates (u, v, w). Seeded and reproducible byte for byte on every platform.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"

namespace stratonix::geo {

/// One nested structure: type "spherical", "exponential" or "gaussian" (practical ranges).
struct VariogramStructure {
  std::string type;
  double contribution = 0;               ///< partial sill
  double range_major = 0, range_minor = 0, range_vertical = 0;  ///< m
  double azimuth = 0;                    ///< degrees clockwise from north, of the major range
};
struct Variogram {
  double nugget = 0;
  std::vector<VariogramStructure> structures;
  double sill() const;
  /// gamma for a separation vector (du, dv, dw) in metres.
  double gamma(const Vec3& h) const;
  /// covariance C(h) = sill - gamma(h).
  double covariance(const Vec3& h) const;
};

struct VariogramDirection {
  double azimuth = 0, dip = 0;           ///< degrees
  double angle_tolerance = 90;           ///< degrees (90: omnidirectional)
  double bandwidth = 1e30;               ///< m
  double lag = 0, lag_tolerance = 0;     ///< m
  std::size_t lags = 0;
};
struct ExperimentalVariogram {
  std::vector<double> distance, gamma;   ///< mean separation and semivariance per lag
  std::vector<std::size_t> pairs;
};
ExperimentalVariogram experimental_variogram(const std::vector<Vec3>& points, const std::vector<double>& values,
                                             const VariogramDirection& d);
/// Weighted least squares (weights N(h) / gamma(h)^2, Cressie 1985) of the free parameters of
/// `start` (contributions and ranges; the nugget when `fit_nugget`).
Variogram fit_variogram(const ExperimentalVariogram& e, const Variogram& start, bool fit_nugget, std::size_t min_pairs);

/// Normal-score transform of a declustered distribution with linear tails to [zmin, zmax].
class NormalScore {
 public:
  NormalScore(std::vector<double> values, std::vector<double> weights, double zmin, double zmax);
  double forward(double z) const;
  double back(double y) const;
 private:
  std::vector<double> z_, y_;
  double zmin_ = 0, zmax_ = 0;
};
/// Cell-declustering weights (cell size `cell`, Deutsch and Journel 1998 `declus`).
std::vector<double> decluster_weights(const std::vector<Vec3>& points, double cell);

struct SearchOptions {
  double radius_major = 0, radius_minor = 0, radius_vertical = 0, azimuth = 0;  ///< m, degrees
  std::size_t max_data = 24, max_per_octant = 0;   ///< 0: no octant limit
  bool ordinary = false;                           ///< ordinary kriging (else simple, mean 0)
};
struct ConditioningPoint { std::size_t node = 0; double value = 0; };

/// Sequential Gaussian simulation of normal scores at `nodes` (geological-time coordinates):
/// conditioning values (already normal scores) at their nodes, optional local means; returns one
/// value per node.
std::vector<double> sgs(const std::vector<Vec3>& nodes, const std::vector<ConditioningPoint>& data, const Variogram& v,
                        const SearchOptions& s, std::uint64_t seed, const std::vector<double>& local_mean = {});
/// Sequential indicator simulation of facies codes 0..K-1 with target proportions and one
/// variogram per code; conditioning codes at their nodes.
std::vector<int> sis(const std::vector<Vec3>& nodes, const std::vector<std::pair<std::size_t, int>>& data,
                     const std::vector<double>& proportions, const std::vector<Variogram>& v, const SearchOptions& s,
                     std::uint64_t seed);

/// Object-based facies: channels (width, thickness, sinuosity, azimuth) or lobes (radius,
/// thickness) placed by a seeded marked point process in (u, v, w) until the target proportion;
/// wells' facies honoured by rejection and repositioning.
struct ObjectSpec {
  std::string kind;                      ///< "channel" or "lobe"
  double width = 0, thickness = 0, sinuosity = 1, wavelength = 0, azimuth = 0, radius = 0;  ///< m, degrees
  double proportion = 0;                 ///< target fraction of the volume
  int code = 1;                          ///< facies code of the objects (background 0)
};
struct ObjectResult { std::vector<int> facies; double proportion = 0; std::size_t well_cells = 0, well_misfit = 0, objects = 0; };
ObjectResult object_facies(const std::vector<Vec3>& nodes, const Vec3& lo, const Vec3& hi, const ObjectSpec& o,
                           const std::vector<std::pair<std::size_t, int>>& wells, std::uint64_t seed, std::size_t max_iterations);

}  // namespace stratonix::geo
