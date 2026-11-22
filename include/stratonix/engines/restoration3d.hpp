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
  bool pin_given = false;
  double pin_x = 0, pin_y = 0;  ///< the vertical line held in plan (default: the domain centre)
  int steps = 8;            ///< increments (finite rotations are not taken as strain)
  /// per model unit: Young's modulus and Poisson's ratio (empty: young and poisson everywhere)
  std::vector<std::pair<double, double>> unit_props;
};
/// Strain of one tetrahedron (Green-Lagrange, from the total deformation present -> restored).
struct StrainSample {
  Vec3 present, restored;       ///< centroid
  double dilatation = 0;        ///< det F - 1
  double max_shear = 0;         ///< largest minus smallest principal strain
  double max_abs = 0;           ///< largest absolute principal strain
  bool at_horizon = false;      ///< the tetrahedron is crossed by the restored horizon
  double bed_stretch = 0;       ///< at the horizon: largest |principal stretch - 1| in the bedding plane
  int unit = -1;
};
struct VolumetricResult {
  bool ok = false;
  std::string reason;
  double max_strain = 0, mean_strain = 0, flatness_rms = 0;
  double horizon_stretch_max = 0;  ///< largest |in-plane stretch - 1| of the restored horizon (bed-length change)
  std::size_t tetrahedra = 0, nodes = 0, fault_contacts = 0, cg_iterations = 0;
  double mesh_spacing = 0;         ///< m (lattice cell size in plan)
  std::vector<StrainSample> samples;
};
VolumetricResult volumetric_restore(const Model& m, const VolumetricSpec& s);

/// Porosity of a unit in its own fault block at the block's burial depth (compaction across faults
/// moving during burial): Athy phi = phi0 exp(-c z).
double block_porosity(double phi0, double c, double burial_depth);

/// Airy tectonic subsidence (Steckler and Watts 1978).
double airy_tectonic_subsidence(double sediment, double rho_s, double rho_m, double rho_w, double water_depth, double sea_level_change);
/// Flexural response of an elastic plate to a periodic load of wavenumber k (fraction of Airy).
double flexural_response(double k, double te, double young, double poisson, double delta_rho);

/// Backstripping of the burial history (STX-TN-15 section 4.4): per deposition step, the sediment
/// above the basement horizon, its mean density from the compaction laws, and the tectonic
/// subsidence by Airy isostasy and, with Te > 0, by an elastic plate (on a map grid).
struct BackstripSpec {
  int basement = -1;                    ///< model horizon index
  double rho_m = 3300, rho_w = 1030, rho_g = 2650;  ///< kg/m3
  double te = 0, young = 70e9, poisson = 0.25;      ///< m, Pa
  std::vector<double> water_depth, sea_level;       ///< per model horizon (when it was at the surface); NaN = not given
  double water_depth_present = 0, sea_level_present = 0;
};
struct TectonicStep {
  std::string label;                    ///< horizon at the surface, or "present"
  double age = kNaN;
  double sediment = 0, density = 0, water_depth = 0, sea_level = 0;
  double airy = 0, flexural = kNaN;     ///< tectonic subsidence (m)
};
class Structure;
/// At one location: Airy only.
std::vector<TectonicStep> backstrip_column(const Structure& st, const BackstripSpec& s, double x, double y);
/// On an n x n grid over the domain, with flexure when te > 0: per step, the Airy and flexural
/// tectonic subsidence grids (row-major, x fastest).
struct BackstripMaps { std::size_t n = 0; std::vector<std::string> labels; std::vector<double> ages; std::vector<std::vector<double>> airy, flexural; };
BackstripMaps backstrip_maps(const Structure& st, const BackstripSpec& s, const Domain& dom, std::size_t n);
/// Bilinear value of a map at (x, y).
double backstrip_map_value(const BackstripMaps& m, const std::vector<double>& grid, const Domain& dom, double x, double y);

/// Eroded section from a porosity log (depth below the unconformity, porosity) against the normal
/// Athy trend phi0 exp(-c z): the maximum burial minus the present depth (m).
double eroded_section(const std::vector<std::pair<double, double>>& log, double phi0, double c);

}  // namespace stratonix
