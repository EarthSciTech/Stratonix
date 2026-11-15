// Grids and properties (stage 3B, STX-TN-12): corner-point grids conforming to faults and
// horizons, tetrahedral grids derived from them, properties per cell with units and evidence
// classes, fault transmissibility multipliers and the upscaling interface. Inputs are plain
// values (no description types); the model interface is unchanged.
#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "stratonix/model/model.hpp"
#include "stratonix/repr/faultrock.hpp"

namespace stratonix::repr {

// ---- corner-point grids (STX-TN-12 section 1) ----------------------------------------------

enum class Layering { Proportional, FollowTop, FollowBase };
const char* layering_name(Layering l);

/// One zone of the grid: a model unit and its layering.
struct ZoneSpec {
  int unit = -1;                  ///< model unit index (ModelInfo::units)
  int layers = 1;                 ///< nk
  Layering layering = Layering::Proportional;
  double layer_thickness = 0;     ///< m, for FollowTop / FollowBase
};

struct GridOptions {
  double dx = 0, dy = 0;          ///< areal spacing (m)
  double margin = 0;              ///< m above the top and below the base horizon
  double min_thickness = 0.01;    ///< m; thinner cells are inactive
  std::vector<ZoneSpec> zones;    ///< consecutive model units, top to base
  bool snap_to_faults = true;     ///< fault-conforming pillars
};

/// Corner-point grid in the Eclipse convention: i fastest, then j, then k (k = 0 at the top).
/// Coordinates are model coordinates with z up (elevation); exporters convert to depth.
struct CornerPointGrid {
  int nx = 0, ny = 0, nz = 0;
  double x0 = 0, y0 = 0, dx = 0, dy = 0;
  /// (nx + 1) (ny + 1) pillars, 6 values each: top x, y, z then base x, y, z.
  std::vector<double> coord;
  /// 8 nx ny nz corner elevations in Eclipse ZCORN order.
  std::vector<double> zcorn;
  std::vector<std::int32_t> actnum;   ///< per cell, 1 active
  std::vector<std::int32_t> zone;     ///< per cell: model unit index
  std::vector<BlockCode> block;       ///< per cell: fault block
  std::vector<std::int32_t> layer_unit;    ///< per layer k: model unit index
  std::vector<std::int32_t> pillar_fault;  ///< per pillar: fault index it lies on, or -1
  std::vector<std::string> inactive_reason; ///< per cell: why it is inactive (empty if active)

  std::size_t ncell() const { return static_cast<std::size_t>(nx) * ny * nz; }
  std::size_t cell(int i, int j, int k) const {
    return static_cast<std::size_t>(i) + static_cast<std::size_t>(nx) * (static_cast<std::size_t>(j) + static_cast<std::size_t>(ny) * k);
  }
  std::size_t pillar(int i, int j) const { return static_cast<std::size_t>(i) + static_cast<std::size_t>(nx + 1) * j; }
  /// Index into zcorn of corner c (bit 0: i + 1, bit 1: j + 1, bit 2: bottom) of cell (i, j, k).
  std::size_t zcorn_index(int i, int j, int k, int c) const;
  /// Corner c of cell (i, j, k): the point on its pillar at the corner's elevation.
  Vec3 corner(int i, int j, int k, int c) const;
  std::array<Vec3, 8> corners(int i, int j, int k) const;
};

/// Extract a corner-point grid from a model by evaluation (never edited by hand).
CornerPointGrid corner_point_grid(const Model& m, const GridOptions& o);

/// Exact volume of the trilinear hexahedron with these corners (bit order as above): the
/// Jacobian integrated by 2 x 2 x 2 Gauss quadrature, exact for trilinear maps.
double cell_volume(const std::array<Vec3, 8>& c);
Vec3 cell_centre(const std::array<Vec3, 8>& c);
/// The six tetrahedra of the hexahedron split along its 0-7 diagonal (used by C24).
std::array<std::array<int, 4>, 6> hex_tets_fixed();
double tet_volume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);

// ---- properties (section 2) -----------------------------------------------------------------

/// One value per cell (NaN where missing), with its unit, evidence class and source.
struct CellProperty {
  std::string name, unit, evidence, source;
  bool integer = false;
  std::vector<double> v;
  std::size_t missing = 0;
  std::string missing_reason;  ///< why the missing cells are missing ("" when none are)
};
struct PropertySet {
  std::vector<CellProperty> props;
  const CellProperty* find(const std::string& name) const;
  CellProperty* find(const std::string& name);
};

struct FaciesSpec {
  int unit = -1;
  int code = 0;
  std::string name, evidence, source;
};
struct TransformSpec {
  std::string name, unit, from, law;  ///< law: linear, exponential, power
  double a = kNaN, b = kNaN;
  std::string evidence, source;
};
struct PropertyInputs {
  std::vector<FaciesSpec> facies;
  bool porosity_from_compaction = false;
  std::vector<TransformSpec> transforms;
};

/// ZONE, SERIES (when the model has series), BLOCK, FACIES, PORO and the transforms.
PropertySet grid_properties(const Model& m, const CornerPointGrid& g, const PropertyInputs& in);
/// Apply one transform to `from` (cells missing in `from` stay missing, with its reason).
CellProperty apply_transform(const CellProperty& from, const TransformSpec& t);

// ---- fault transmissibility (section 3) -----------------------------------------------------

/// Manzocchi et al. (1999) fault permeability, mD: log10 k = -4 SGR - 1/4 log10(D) (1 - SGR)^5.
double manzocchi_fault_permeability(double sgr, double displacement_m);
/// Manzocchi et al. (1999) eq. 4: transmissibility multiplier of a fault of thickness tf (m)
/// and permeability kf between cells of permeabilities ki, kj and lengths Li, Lj normal to the face.
double manzocchi_multiplier(double tf, double kf, double ki, double kj, double Li, double Lj);

struct FaultFace {
  int i = 0, j = 0, k = 0;  ///< the cell on the low side
  char dir = 'X';           ///< 'X': face between (i, j) and (i + 1, j); 'Y': (i, j) and (i, j + 1)
  int fault = -1;
  Vec3 centre;
  double sgr = kNaN, thickness = kNaN, displacement = kNaN, k_fault = kNaN;
  double k_i = kNaN, k_j = kNaN, L_i = kNaN, L_j = kNaN, mult = kNaN;
  std::string missing;      ///< why mult is missing ("" when computed)
};
struct TransmissibilityResult {
  std::vector<FaultFace> faces;
  CellProperty multx, multy;
};
/// Fault faces of the grid and their multipliers; needs PERMX and PERMY in `props`.
TransmissibilityResult fault_transmissibility(const Model& m, const CornerPointGrid& g, const PropertySet& props,
                                              const FaultRockInput& rock);

// ---- upscaling interface (section 4) --------------------------------------------------------

/// Map a property from a fine grid to a coarse grid of the same model. Implementations are
/// registry plugins of kind Upscaler.
struct Upscaler {
  virtual ~Upscaler() = default;
  virtual CellProperty upscale(const CornerPointGrid& fine, const CellProperty& p, const CornerPointGrid& coarse) const = 0;
};
/// Reference implementation: volume-weighted mean of the fine cells whose centres fall in each
/// coarse cell.
CellProperty upscale_volume_weighted_mean(const CornerPointGrid& fine, const CellProperty& p, const CornerPointGrid& coarse);

// ---- tetrahedral grids (section 5) ----------------------------------------------------------

struct TetGrid {
  std::vector<Vec3> nodes;
  std::vector<std::array<std::int32_t, 4>> tets;
  std::vector<std::int32_t> cell;   ///< per tetrahedron: the corner-point cell it comes from
};
/// Split every active cell into six tetrahedra by the consistent rule of Dompierre et al. (1999).
TetGrid tetrahedral_grid(const CornerPointGrid& g);

// ---- checks support (section 7) -------------------------------------------------------------

struct GridDefects {
  std::size_t nonpositive_volume = 0, inverted = 0, fault_face_off = 0, fault_faces = 0;
  double worst_volume = 0, worst_fault_distance = 0;
  std::string first_nonpositive, first_inverted, first_fault_off;
};
/// Defects of C23 (volume), C24 (inversion) and C25 (fault faces within `fault_tol` m).
GridDefects grid_defects(const Model& m, const CornerPointGrid& g, double fault_tol);

}  // namespace stratonix::repr
