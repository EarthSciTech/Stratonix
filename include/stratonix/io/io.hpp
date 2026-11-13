/// @file io.hpp
/// @brief Input/output plugins: importers (data into Stratonix) and exporters (models out),
/// plus the plain format functions they use. All text is written with '\n' line endings and
/// locale-independent numbers, so files are byte-identical on Windows and Linux.
#pragma once
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "stratonix/core/grid.hpp"
#include "stratonix/core/mesh.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/registry/registry.hpp"

namespace stratonix::io {

struct NamedLine { std::string name; int id = 0; Polyline pts; };
struct NamedMesh { std::string name; TriMesh mesh; };
struct ImportedTop { std::string well, horizon; Vec3 p; double md = kNaN; };

/// Data read by an importer. Coordinates are as in the file (z as written; see z_is).
struct ImportedData {
  std::string kind;                 ///< "grid", "points", "lines", "meshes", "tops"
  std::optional<Grid2D> grid;
  std::vector<Vec3> points;
  std::vector<NamedLine> lines;
  std::vector<NamedMesh> meshes;
  std::vector<ImportedTop> tops;
  std::string z_is = "as_written";  ///< "elevation", "depth" or "as_written"
  std::vector<std::string> notes;
};

/// Importer plugin interface.
class Importer : public registry::PluginBase {
 public:
  virtual std::string data_kind() const = 0;
  virtual ImportedData read(const std::filesystem::path& file) const = 0;
};

// ---- grid formats --------------------------------------------------------------------------
constexpr double kIrapUndef = 9999900.0;
constexpr double kZmapNull = -99999.0;
/// IRAP classic ASCII grid. Header: -996 nrow xinc yinc / xmin xmax ymin ymax / ncol rot xrot yrot
/// / 0 0 0 0 0 0 0; values x fastest from (xmin, ymin), 6 per line, undefined = 9999900.0.
std::string write_irap_ascii(const Grid2D& g, int decimals = 4);
Grid2D read_irap_ascii(const std::string& text);
/// ZMAP+ grid: header between '@' lines; values column by column from the top-left corner.
std::string write_zmap(const Grid2D& g, const std::string& name, const std::vector<std::string>& comments, int decimals = 4);
Grid2D read_zmap(const std::string& text);
/// Regular-grid XYZ points (whitespace or comma separated; non-numeric lines skipped).
Grid2D read_xyz_grid(const std::string& text);

/// Read the input grids a description refers to (surfaces and isochores) with the registered
/// importers.
void load_build_inputs(const desc::Description& d, BuildInputs& in);

/// Register this layer's plugins.
void register_io_plugins();

}  // namespace stratonix::io
