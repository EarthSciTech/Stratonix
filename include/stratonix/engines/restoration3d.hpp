// Stage 4B (STX-TN-15 section 4): flexural-slip and volumetric restoration, compaction across
// moving faults, flexural backstripping and eroded-section estimates.
#pragma once
#include <string>
#include <utility>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/model/model.hpp"

namespace stratonix {

/// Flexural-slip unfolding of a horizon polyline (x, z) in a section to its datum elevation,
/// pinned at `pin_x`: bed length measured along the horizon is conserved.
struct UnfoldResult { std::vector<Vec3> restored; double present_length = 0, restored_length = 0; };
UnfoldResult flexural_slip_unfold(const std::vector<Vec3>& horizon, double pin_x, double datum_z);

/// 3D volumetric restoration (decision 18): linear elasticity on a tetrahedral mesh, the horizon
/// flattened to its datum, faults as sliding contacts.
struct VolumetricSpec {
  std::string horizon;
  double datum_z = 0;
  double young = 1e10, poisson = 0.25;
  double max_strain = 0.1;  ///< labelled limit (C31)
  double mesh_size = 0;     ///< 0: largest extent / 30
};
struct VolumetricResult {
  bool ok = false;
  std::string reason;
  double max_strain = 0, mean_strain = 0, flatness_rms = 0;
  std::size_t tetrahedra = 0;
};
VolumetricResult volumetric_restore(const Model& m, const VolumetricSpec& s);

/// Porosity of a unit in its own fault block at the block's burial depth (compaction across faults
/// moving during burial): Athy phi = phi0 exp(-c z).
double block_porosity(double phi0, double c, double burial_depth);

/// Airy tectonic subsidence (Steckler and Watts 1978).
double airy_tectonic_subsidence(double sediment, double rho_s, double rho_m, double rho_w, double water_depth, double sea_level_change);
/// Flexural response of an elastic plate to a periodic load of wavenumber k (fraction of Airy).
double flexural_response(double k, double te, double young, double poisson, double delta_rho);

/// Eroded section from a porosity log (depth below the unconformity, porosity) against the normal
/// Athy trend phi0 exp(-c z): the maximum burial minus the present depth (m).
double eroded_section(const std::vector<std::pair<double, double>>& log, double phi0, double c);

}  // namespace stratonix
