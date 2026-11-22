#pragma once
// Exporters: file formats (writers and readers) and the exporter plugins that write a built model.
//
// Conventions for every file written here:
// - z is elevation (negative below the datum), metres, unless a format states otherwise;
// - '\n' line endings and locale-independent number formatting (std::to_chars), so files are
//   byte-identical on every platform for the same model;
// - names are restricted to letters, digits, '_' and '-' (petrel_name()).
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "stratonix/core/grid.hpp"
#include "stratonix/core/mesh.hpp"
#include "stratonix/core/well.hpp"
#include "stratonix/model/model.hpp"
#include "stratonix/registry/registry.hpp"
#include "stratonix/repr/faultrock.hpp"
#include "stratonix/repr/grid.hpp"
#include "stratonix/repr/repr.hpp"

namespace stratonix::io {

/// Name safe for Petrel and file systems: letters, digits, '_' and '-'; other characters become '_'.
std::string petrel_name(const std::string& s);

// ---- IRAP classic lines (fault sticks, fault polygons) -------------------------------------
/// X Y Z per line; consecutive lines are separated by "999.000000 999.000000 999.000000".
constexpr double kIrapLineSeparator = 999.0;
std::string write_irap_lines(const std::vector<Polyline>& lines, int decimals = 4);
std::vector<Polyline> read_irap_lines(const std::string& text);

// ---- GOCAD TSurf --------------------------------------------------------------------------
struct TSurfInfo {
  std::string name;
  std::string z_positive = "Elevation";  ///< "Elevation" or "Depth"
  std::string crs_name = "Default";      ///< GOCAD_ORIGINAL_COORDINATE_SYSTEM NAME
  double r = 0.7, g = 0.7, b = 0.7;
  /// Optional vertex property (written with PROPERTIES / PVRTX); empty name = none.
  std::string property;
  std::vector<double> values;            ///< one per mesh vertex when property is set
};
std::string write_tsurf(const TriMesh& m, const TSurfInfo& info, int decimals = 4);
TriMesh read_tsurf(const std::string& text, TSurfInfo* info = nullptr);

// ---- Petrel points with attributes ----------------------------------------------------------
struct PointsTable {
  std::vector<Vec3> pts;
  /// Attribute columns: name and type ("String", "Int" or "Float"), values as text per point.
  std::vector<std::pair<std::string, std::string>> attr;
  std::vector<std::vector<std::string>> values;  ///< values[point][attribute]
};
std::string write_petrel_points(const PointsTable& t, const std::vector<std::string>& comments, int decimals = 4);
PointsTable read_petrel_points(const std::string& text);

// ---- Wells ---------------------------------------------------------------------------------
struct WellHead { std::string name; double x = 0, y = 0, kb = 0, td_md = 0; };
std::string write_well_heads(const std::vector<WellHead>& w, const std::vector<std::string>& comments, int decimals = 2);
std::vector<WellHead> read_well_heads(const std::string& text);
struct WellTopRow { std::string well, surface; double x = 0, y = 0, z = 0, md = 0; };
std::string write_well_tops(const std::vector<WellTopRow>& t, const std::vector<std::string>& comments, int decimals = 2);
std::vector<WellTopRow> read_well_tops(const std::string& text);
/// Deviation survey as MD X Y Z (elevation) rows with a commented header.
std::string write_well_path(const WellPath& w, int decimals = 2);

// ---- GRDECL (Eclipse keyword format, stage 3B) ----------------------------------------------
/// Grid keywords SPECGRID, COORD, ZCORN, ACTNUM (depth positive down), the integer properties
/// ZONE, BLOCK and FACIES, real properties, MULTX / MULTY and FAULTS; every keyword has a comment
/// header with provenance, unit and evidence class.
struct GrdeclInfo {
  std::vector<std::string> header;        ///< comment lines at the top of the file
  std::vector<std::string> fault_names;   ///< by fault index, for FAULTS
  const std::vector<repr::FaultFace>* faces = nullptr;  ///< for FAULTS (none written without)
};
/// Missing values of real properties in GRDECL (stated in the file's header).
constexpr double kGrdeclMissing = -999.25;
std::string write_grdecl(const repr::CornerPointGrid& g, const repr::PropertySet* props,
                         const repr::TransmissibilityResult* trans, const GrdeclInfo& info, int decimals = 3);

// ---- Neutral formats -------------------------------------------------------------------------
/// VTK XML unstructured grid (ASCII) of triangles or polylines with optional point and cell arrays.
struct VtuData {
  std::vector<Vec3> pts;
  std::vector<std::vector<std::int64_t>> cells;   ///< vertex ids per cell
  std::vector<std::uint8_t> types;                ///< VTK cell types (5 triangle, 4 polyline)
  std::vector<std::pair<std::string, std::vector<double>>> point_data;
  std::vector<std::pair<std::string, std::vector<std::int32_t>>> cell_data;
  std::vector<std::pair<std::string, std::vector<double>>> cell_data_real;  ///< Float64 cell arrays
};
std::string write_vtu(const VtuData& d, int decimals = 4);
VtuData read_vtu(const std::string& text);
/// VTK unstructured grids of the active hexahedra (cell type 12) and of a tetrahedral grid
/// (type 10), with ZONE, BLOCK and the properties as cell data.
VtuData grid_vtu(const repr::CornerPointGrid& g, const repr::PropertySet* props);
VtuData tet_vtu(const repr::TetGrid& t, const repr::CornerPointGrid& g, const repr::PropertySet* props);
/// VTK XML image data of a voxel grid (cell arrays unit and block).
std::string write_vti(const repr::VoxelGrid& v, int decimals = 4);
repr::VoxelGrid read_vti(const std::string& text);
std::string write_obj(const TriMesh& m, const std::string& name, const std::vector<std::string>& comments, int decimals = 4);
TriMesh read_obj(const std::string& text);
std::string write_stl(const TriMesh& m, const std::string& name, int decimals = 4);
TriMesh read_stl(const std::string& text);
/// Comma-separated table with a header line; numbers and plain strings (no commas or quotes).
struct CsvTable { std::vector<std::string> header; std::vector<std::vector<std::string>> rows; };
std::string write_csv(const CsvTable& t, const std::vector<std::string>& comments);
CsvTable read_csv(const std::string& text);
/// Points "x y z" per line (comment lines start with '#').
std::vector<Vec3> read_xyz_points(const std::string& text);

// ---- Export context ------------------------------------------------------------------------
/// What the exporters write from: the model, wells, grid and spacings, with lazily computed and
/// cached representations so several formats share one evaluation.
class ExportContext {
 public:
  ExportContext(const Model& m, std::vector<WellPath> wells, const GridSpec& grid, double fault_spacing, double stick_spacing);
  const Model& model() const { return m_; }
  const ModelInfo& info() const { return m_.info(); }
  const std::vector<WellPath>& wells() const { return wells_; }
  const GridSpec& grid() const { return grid_; }
  double stick_spacing() const { return stick_spacing_; }
  int decimals = 4;
  double vertical_exaggeration = 1;          ///< for viewer state files
  double voxel_dx = 0, voxel_dy = 0, voxel_dz = 0;  ///< voxel grid spacing (0: none)
  double sample_spacing = 250;               ///< point sampling for GemPy / LoopStructural sets
  /// Stage 2C series, youngest first (name, base relation, horizons), for input sets of other
  /// modellers (one feature per series); empty without series.
  struct SeriesExport { std::string name, base; std::vector<std::string> horizons; };
  std::vector<SeriesExport> series;
  /// Stage 3A: horizons built from the horizon above minus a thickness (implicit.isochore_horizons).
  std::vector<std::string> isochore_horizons;
  /// Stage 2F (STX-TN-10): fault-surface attributes and curvature grids, written only when given.
  bool fault_rock = false;
  repr::FaultRockInput rock;
  repr::StressInput stress;
  std::vector<int> curvature_horizons;
  double fracture_threshold = kNaN;
  /// Stage 3B (STX-TN-12): the grid of the description's grid block and what goes with it.
  struct CellGrid {
    const repr::CornerPointGrid* grid = nullptr;
    const repr::PropertySet* props = nullptr;
    const repr::TransmissibilityResult* trans = nullptr;
    const std::vector<repr::FaultFace>* faces = nullptr;
    const repr::TetGrid* tets = nullptr;
  };
  CellGrid cells;
  /// Lines written at the top of every file that takes comments (model, CRS, z convention).
  std::vector<std::string> comments() const;

  const repr::FaultedSurface& faulted(int h);
  /// Share the faulted surfaces already made by the checks of the same model (stage 4B, I-020).
  void share_surfaces(repr::FaultedSurfaceCache* c) { shared_ = c; }
  const repr::FaultSurface& fault_surface(int f);
  const std::vector<repr::CutoffLine>& cutoffs(int f, int h);
  const std::vector<repr::WellTop>& tops(std::size_t well);
  /// Per-block surface of horizon h for block b, kept on the block's footprint on that horizon
  /// dilated by `dilate` nodes (so it reaches across the fault gap); NaN elsewhere.
  Grid2D block_surface(int h, BlockCode b, int dilate = 2);
  /// Blocks present on horizon h (footprint not empty), ascending.
  std::vector<BlockCode> blocks_on(int h);
  /// "F1-HW_F2-FW" style label of a block (sides of every fault).
  std::string block_label(BlockCode b) const;
  /// Fault polygons of horizon h on fault f: footwall cut-off, then hanging-wall cut-off reversed.
  std::vector<Polyline> fault_polygons(int f, int h);

 private:
  const Model& m_;
  std::vector<WellPath> wells_;
  GridSpec grid_;
  double fault_spacing_, stick_spacing_;
  std::map<int, repr::FaultedSurface> faulted_;
  repr::FaultedSurfaceCache* shared_ = nullptr;
  std::map<int, repr::FaultSurface> fsurf_;
  std::map<std::pair<int, int>, std::vector<repr::CutoffLine>> cut_;
  std::map<std::size_t, std::vector<repr::WellTop>> tops_;
};

/// One written file and how to load it.
struct ExportedFile {
  std::string path;          ///< relative to the export folder, '/' separators
  std::string format;        ///< e.g. "IRAP classic ASCII grid"
  std::string petrel_import; ///< Petrel import dialog file type, or "-" when not for Petrel
  std::string object;        ///< what it holds (e.g. "H2 faulted surface, gaps at faults")
};

/// Exporter plugin: writes one format for a model into a folder.
class Exporter : public registry::PluginBase {
 public:
  virtual std::string id() const = 0;
  virtual bool petrel() const { return true; }
  virtual std::vector<ExportedFile> write(ExportContext& c, const std::filesystem::path& folder) const = 0;
};

/// The Petrel export set (exporter ids, in writing order).
std::vector<std::string> petrel_set_ids();
/// Write the listed exporters into folder (must be new or empty) plus README.txt describing every
/// file, the coordinate system, units, z convention and the Petrel import option.
std::vector<ExportedFile> export_set(ExportContext& c, const std::vector<std::string>& ids, const std::filesystem::path& folder,
                                     const std::string& evidence_note);
/// README text for an export folder.
std::string export_readme(const ExportContext& c, const std::vector<ExportedFile>& files, const std::string& evidence_note);

void register_exporters();
void register_importers();

}  // namespace stratonix::io
