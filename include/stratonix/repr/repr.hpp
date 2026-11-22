/// @file repr.hpp
/// @brief Representations derived from any Model by evaluation only: horizon grids (per block
/// and faulted with gaps), fault meshes, cut-off lines, throw tables, well intersections,
/// sections, voxel grids and closure analysis.
///
/// Each kernel has a reference implementation (plain, robust, slow). Where an optimised version
/// exists it is selected by optimised_kernels() and must match the reference in parity tests.
#pragma once
#include <array>
#include <cstdint>
#include <map>
#include <mutex>
#include <string>
#include <tuple>
#include <vector>

#include "stratonix/core/grid.hpp"
#include "stratonix/core/mesh.hpp"
#include "stratonix/core/well.hpp"
#include "stratonix/model/model.hpp"

namespace stratonix::repr {

// ---- blocks --------------------------------------------------------------------------------
struct BlockList {
  std::vector<BlockCode> codes;    ///< distinct block codes found in the domain, ascending
  std::vector<std::size_t> counts; ///< samples per code
  /// per code: the map extent of its samples, x0, y0, x1, y1 (m), widened by one lattice cell
  std::vector<std::array<double, 4>> extent;
  double cell = 0;                 ///< lattice cell size (m), the larger of the two directions
};
/// Sample the domain on an n x n x nz lattice and list the block codes present.
BlockList list_blocks(const Model& m, std::size_t n = 48, std::size_t nz = 48);

// ---- horizon grids -------------------------------------------------------------------------
/// Root of a horizon's scalar along verticals, as seen from one block (extended beyond it).
/// Reference: bisection. Optimised: bracketed secant (Illinois). NaN if no root within
/// the domain's elevation range extended by `margin`.
Grid2D horizon_grid_block(const Model& m, int horizon, BlockCode block, const GridSpec& g, double margin = 2000);

struct FaultedSurface {
  Grid2D z;                         ///< elevation where the horizon is present; NaN in gaps
  std::vector<BlockCode> block;     ///< block of each defined node
  std::size_t gap_nodes = 0;        ///< nodes without the horizon inside the domain (fault gaps)
  std::size_t overlap_nodes = 0;    ///< nodes where the horizon occurs twice (reverse faults); upper kept
};
/// Faulted horizon surface with fault gaps.
/// Reference: candidate blocks from dense sampling of each column. Optimised: block iteration
/// from the per-block root. Parity tested.
FaultedSurface horizon_grid_faulted(const Model& m, int horizon, const GridSpec& g);

/// Faulted horizon surfaces of one model, each computed once per horizon and lattice and shared by
/// the checks, figures and exports of a build (stage 4B, I-020). Entries never move once made.
class FaultedSurfaceCache {
 public:
  explicit FaultedSurfaceCache(const Model& m) : m_(&m) {}
  const Model& model() const { return *m_; }
  const FaultedSurface& get(int horizon, const GridSpec& g);

 private:
  const Model* m_;
  std::mutex mu_;
  std::map<std::tuple<int, long long, long long, long long, long long, std::size_t, std::size_t>, FaultedSurface> map_;
};
/// The cached surface when `cache` is of `m`; otherwise computed into `scratch`, which is returned.
const FaultedSurface& faulted_surface(const Model& m, int horizon, const GridSpec& g, FaultedSurfaceCache* cache, FaultedSurface& scratch);

// ---- stratigraphic terminations (stage 2C, STX-TN-07 section 2) ----------------------------
/// Vertical thickness of model unit `unit` (between horizons unit-1 and unit) in the present
/// geometry, per node where both horizons are in one fault block; NaN in fault gaps.
Grid2D unit_thickness_grid(const Model& m, int unit, const GridSpec& g);
/// Where a unit's thickness falls to zero, with the kind of termination from ModelInfo
/// ("onlap", "truncation", "pinch-out", "clinoform"); empty for a unit that does not terminate.
struct TerminationLine {
  std::string unit, kind;
  Polyline line;
};
std::vector<TerminationLine> termination_lines(const Model& m, int unit, const GridSpec& g);

// ---- restoration (stage 2E, STX-TN-09 section 3) --------------------------------------------
struct Section;
/// Flexural-slip unfolding of a section: each horizon's trace straightened along its own length
/// from a pin at the section's start, fault gaps closed, at its depth at the pin.
struct UnfoldedHorizon {
  int horizon = -1;
  double traced_length = 0;  ///< bed length along the section, fault gaps closed
  Polyline restored;         ///< (along-section distance, 0, elevation) points of the straightened bed
};
std::vector<UnfoldedHorizon> unfold_flexural(const Model& m, const Section& s);
/// A horizon's trace along a section split into runs of one fault block (fault gaps between runs).
std::vector<std::vector<Vec3>> horizon_segments(const Model& m, const Section& s, int horizon);

// ---- faults --------------------------------------------------------------------------------
struct FaultSurface {
  int fault = -1;
  TriMesh mesh;                     ///< clipped to the tip line and the domain
  std::size_t na = 0, nb = 0;       ///< lattice size on the approximate plane
  std::vector<Vec3> nodes;          ///< lattice nodes on the fault surface (na x nb, a fastest)
  std::vector<double> extent;       ///< extent coordinate per node (<= 1 inside)
  std::vector<double> a, b;         ///< plane coordinates per lattice column / row
};
/// Fault surface located by root-finding the fault distance along the approximate normal.
/// `extend` enlarges the lattice beyond the tip line (fraction of the extent) so callers can
/// sample just outside it; the mesh always covers only the region inside the tip line.
FaultSurface fault_surface(const Model& m, int fault, double spacing, double extend = 0.0);

/// Fault sticks: polylines down dip at a strike spacing, inside the tip line.
std::vector<Polyline> fault_sticks(const FaultSurface& fs, double spacing);

struct CutoffLine {
  int fault = -1, horizon = -1;
  bool hanging_wall = true;
  std::vector<Polyline> parts;      ///< continuous parts
  std::vector<Vec3> by_column;      ///< crossing per lattice column of the fault surface (NaN if none)
};
/// Hanging-wall and footwall cut-off lines of a horizon on a fault (eps: offset from the plane).
std::vector<CutoffLine> cutoff_lines(const Model& m, const FaultSurface& fs, int horizon, double eps = 0.25);

struct ThrowRow {
  int fault = -1, horizon = -1;
  bool cut = false;                 ///< the horizon is offset by the fault
  double throw_max = 0;             ///< max vertical separation of cut-offs along strike (m)
  double throw_mid = 0;             ///< at the strike position nearest the reference point
  double heave_max = 0;             ///< max horizontal separation (m)
  double jump_max = 0;              ///< max scalar jump across the fault at the cut-off (m)
  std::size_t samples = 0;
};
/// Jump of each horizon's scalar across a surface at points with the given normal: the
/// hanging-wall value minus the footwall value, each extrapolated linearly to the point from
/// samples at eps and 2*eps on its side (exact for smooth linear fields). Result: [h * n + i].
std::vector<double> scalar_jumps(const Model& m, const std::vector<Vec3>& pts, const Vec3& normal, double eps = 0.25);

ThrowRow throw_row(const Model& m, const FaultSurface& fs, const std::vector<CutoffLine>& lines, int horizon);

// ---- wells ---------------------------------------------------------------------------------
struct WellTop {
  std::string well, horizon;
  int horizon_index = -1;
  Vec3 p;
  double md = 0;
  BlockCode block = 0;
};
/// Horizon crossings along a well path (handles fault crossings: a horizon top is recorded only
/// where the scalar changes sign inside one block).
std::vector<WellTop> well_intersections(const Model& m, const WellPath& w, double step = 2.0);

// ---- sections ------------------------------------------------------------------------------
struct Section {
  std::string name;
  Polyline line;                    ///< polyline vertices (z ignored)
  double length = 0, ds = 0, dz = 0, z_top = 0, z_base = 0;
  std::size_t ns = 0, nz = 0;       ///< columns along the line, rows from top
  std::vector<double> station;      ///< distance along the line per column
  std::vector<Vec3> xy;             ///< map position per column
  std::vector<std::int32_t> unit;   ///< nz x ns, row-major from the top
  std::vector<BlockCode> block;     ///< nz x ns
  std::vector<std::vector<double>> horizon_z;  ///< per horizon, per column: elevation or NaN (gap)
  std::vector<std::vector<double>> fault_z;    ///< per fault, per column: elevation of the fault or NaN
  /// Sections where fault blocks overlap (reverse faults): per horizon, the horizon in each fault
  /// block along the section (NaN where that block does not hold it), so both sheets are drawn.
  /// Empty when no block overlaps another along the section.
  std::vector<std::vector<std::pair<BlockCode, std::vector<double>>>> horizon_block_z;
};
Section make_section(const Model& m, const std::string& name, const Polyline& line, double ds, double dz);

// ---- voxels --------------------------------------------------------------------------------
struct VoxelGrid {
  double x0 = 0, y0 = 0, z0 = 0, dx = 0, dy = 0, dz = 0;  ///< centre of cell (0,0,0); z0 lowest
  std::size_t nx = 0, ny = 0, nz = 0;
  std::vector<std::int32_t> unit;   ///< i + nx*(j + ny*k)
  std::vector<BlockCode> block;
  Vec3 centre(std::size_t i, std::size_t j, std::size_t k) const {
    return {x0 + i * dx, y0 + j * dy, z0 + k * dz};
  }
};
VoxelGrid voxel_grid(const Model& m, double dx, double dy, double dz, bool with_block = true);

// ---- closure -------------------------------------------------------------------------------
struct Closure {
  bool found = false;
  double crest_z = 0, spill_z = 0, closure = 0, area = 0;
  Vec3 crest, spill_point;
  Grid2D mask;                      ///< 1 inside the closure, 0 outside, NaN undefined
  bool spill_at_boundary = false;   ///< spills over the domain edge (not over a saddle)
};
/// Structural closure of a surface (elevation grid, NaN = gap). With faults_seal, gaps are
/// barriers; otherwise reaching a gap counts as spilling.
Closure closure_analysis(const Grid2D& z, bool faults_seal);

// ---- meshes from grids ---------------------------------------------------------------------
/// Triangulated surface from a grid: two triangles per fully defined cell, one where exactly
/// three corners are defined; nothing across gaps. Vertices only at defined nodes.
TriMesh grid_to_mesh(const Grid2D& g);

// ---- contours ------------------------------------------------------------------------------
/// Contour lines of a grid at a level (marching squares; NaN cells skipped).
std::vector<Polyline> contour_lines(const Grid2D& g, double level);

}  // namespace stratonix::repr
