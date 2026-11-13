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
  double r = 0.7, g = 0.7, b = 0.7;
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
  /// Lines written at the top of every file that takes comments (model, CRS, z convention).
  std::vector<std::string> comments() const;

  const repr::FaultedSurface& faulted(int h);
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

}  // namespace stratonix::io
