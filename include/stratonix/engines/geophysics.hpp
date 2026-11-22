/// @file geophysics.hpp
/// @brief Stage 4A (STX-TN-14 section 7): gravity and magnetic forward modelling of a model
/// discretised into right rectangular prisms.
#pragma once
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/model/model.hpp"

namespace stratonix {

/// Vertical gravity anomaly (mGal, positive down) at `obs` of the prism [x1,x2]x[y1,y2]x[z1,z2]
/// (elevations, z up) with density contrast rho (kg/m^3) (Nagy et al. 2000).
double prism_gz(const Vec3& obs, double x1, double x2, double y1, double y2, double z1, double z2, double rho);
/// Total-field anomaly (nT) at `obs` of a uniformly magnetised prism with magnetisation M (A/m,
/// vector, z up) for an inducing field direction F (unit vector, z up) (Bhattacharyya 1964;
/// Blakely 1995).
double prism_dt(const Vec3& obs, double x1, double x2, double y1, double y2, double z1, double z2, const Vec3& M, const Vec3& F);
/// Unit vector of a field with inclination and declination (rad; inclination positive down,
/// declination clockwise from north) in x east, y north, z up.
Vec3 field_direction(double inclination, double declination);

struct PotentialFieldSpec {
  double cell = 100;                   ///< m prism size (x, y, z)
  std::vector<double> density;         ///< kg/m^3 per unit (NaN = not given)
  std::vector<double> susceptibility;  ///< SI per unit (NaN = not given)
  double background_density = 2670, background_susceptibility = 0;
  double inclination = 0, declination = 0, intensity = 50000;  ///< rad, rad, nT
};

/// Anomalies at stations (mGal / nT) of the model's units between z_top and z_base of its domain.
std::vector<double> gravity_anomaly(const Model& m, const PotentialFieldSpec& s, const std::vector<Vec3>& stations);
std::vector<double> magnetic_anomaly(const Model& m, const PotentialFieldSpec& s, const std::vector<Vec3>& stations);

/// The spec of a description's geophysics block for a model (densities and susceptibilities by unit
/// name; cell default outputs.voxel_spacing x, else the largest horizontal extent / 64).
PotentialFieldSpec potential_field_spec(const desc::Description& d, const ModelInfo& info);

/// Observed against modelled at every station of the description's geophysics block.
struct StationMisfit {
  std::string kind, input;  ///< gravity | magnetics
  Vec3 p;
  double observed = 0, modelled = 0, sigma = kNaN;  ///< mGal or nT
};
struct PotentialFieldResult {
  std::vector<StationMisfit> stations;
  double chi2_gravity = kNaN, chi2_magnetics = kNaN;  ///< reduced chi-square (stations with a sigma)
  std::size_t n_gravity = 0, n_magnetics = 0;
};
PotentialFieldResult potential_field_misfit(const Model& m, const desc::Description& d, const BuildInputs& in);

}  // namespace stratonix
