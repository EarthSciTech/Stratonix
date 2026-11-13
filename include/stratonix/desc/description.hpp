/// @file description.hpp
/// @brief Typed, validated model description (schema version 1) with evidence on every parameter.
#pragma once
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "stratonix/core/crs.hpp"
#include "stratonix/core/geometry.hpp"
#include "stratonix/core/units.hpp"
#include "stratonix/desc/document.hpp"

namespace stratonix::desc {

constexpr int kSchemaVersion = 1;

/// Evidence class carried by every parameter.
enum class Evidence { Measured, Interpreted, Calibrated, Upscaled, LiteratureDerived, Assumed, UserInstruction };
const char* evidence_name(Evidence e);
/// Plain-words description used in README files ("assumed by the modeller", ...).
const char* evidence_plain(Evidence e);

/// A labelled parameter: value(s), unit, evidence class and source.
struct Param {
  bool present = false;
  std::string path;               ///< field path in the description, e.g. faults[0].throw.max
  std::vector<double> values;     ///< flattened values in the given unit
  std::vector<std::size_t> shape; ///< {} scalar, {n} list, {n, m} list of tuples
  bool flag = false;              ///< value for unit "flag"
  std::string unit;
  Dimension dim = Dimension::Any;
  double to_si = 1;
  Evidence evidence = Evidence::Assumed;
  std::string source;
  std::string text;               ///< value as written (for reports)

  double si(std::size_t i = 0) const { return values.at(i) * to_si; }
  std::vector<double> si_list() const;
  std::size_t rows() const { return shape.empty() ? 1 : shape[0]; }
  std::size_t cols() const { return shape.size() < 2 ? 1 : shape[1]; }
};

struct ModelMeta { std::string name, version, status, title, description, contact; };
struct InputDesc { std::string id, path, role, sharing, source, evidence; std::filesystem::path resolved; };
struct WellDesc {
  std::string name, role;
  Param head, datum_elevation, total_depth, stations;
  std::vector<std::pair<std::string, Param>> tops;  ///< horizon name -> depth (positive down)
};
struct HorizonDesc {
  std::string name;
  Param depth, thickness_above;
  std::string surface_input, surface_format, surface_z_is = "depth";
  std::string iso_from, iso_input, iso_format;
  Param iso_thickness;
  bool intended_offset = false;
};
struct UnitDesc { std::string name, top, base, role, color; };
struct PrimitiveDesc {
  std::string type, name;
  Param crest, axis_azimuth, plunge, half_length, half_width_left, half_width_right;
  Param radius_major, radius_minor, rim_depth, rim_width;
  std::vector<std::pair<std::string, Param>> relief;  ///< horizon name -> relief
};
struct TipDesc { std::string rule; std::string horizon, unit; Param offset, fraction, depth; };
struct ProfileDesc { std::string type = "flat_topped"; Param taper; };
struct FaultDesc {
  std::string name;
  Param age;
  std::string geom_type = "planar";
  Param reference_point, reference_depth, strike, dip, trace, trace_depth, decay_length;
  std::string dip_side;
  TipDesc upper, lower;
  Param length;
  std::string tip_line_shape = "rectangular";
  Param throw_max;
  ProfileDesc profile_dip, profile_strike;
  std::string slip_mode, hw_type;
  Param influence_distance;
  std::string abuts_fault, abuts_side;
};
struct SectionDesc { std::string name; Param polyline, detail_window; };

/// A complete, validated model description.
struct Description {
  std::filesystem::path file;
  std::string file_sha256;
  int schema_version = kSchemaVersion;
  Node resolved;  ///< the document after migration (written to the run manifest)

  ModelMeta model;
  CrsInfo crs;
  Param origin, size, depth_range, grid_spacing;
  std::vector<InputDesc> inputs;
  std::vector<WellDesc> wells;
  std::string ref_well;
  Param ref_point;
  std::vector<HorizonDesc> horizons;
  std::vector<UnitDesc> units;
  std::string above_top_name, below_base_name;
  std::string regional_type = "flat";
  Param regional_dip, regional_dip_azimuth;
  std::string base_type = "conformable";
  Param correction_range;  ///< range of the well-correction field (needed when >1 control well or imported surfaces)
  std::vector<PrimitiveDesc> primitives;
  std::vector<FaultDesc> faults;
  std::string frame_well;
  Param frame_point;
  // checks
  Param tol_thickness_min, tol_crossing, tol_offset_outside_tips, tol_throw, tol_wells, tol_unoffset;
  std::string closure_horizon;
  Param closure_faults_seal, target_closure, target_area, target_spill_depth, target_tolerance;
  std::vector<std::string> warn_only;
  // outputs
  std::vector<SectionDesc> sections;
  Param vertical_exaggeration, voxel_spacing, fault_stick_spacing;
  std::vector<std::string> exports;

  bool is_design() const { return model.status == "DESIGN"; }
  int horizon_index(const std::string& name) const;
  int unit_index(const std::string& name) const;
  int fault_index(const std::string& name) const;
  const WellDesc* well(const std::string& name) const;
  /// Every parameter present in the description (for the parameters table and manifest).
  std::vector<const Param*> all_params() const;
  Domain domain() const;
};

/// Result of validation: issues with field paths and rule identifiers.
struct ValidationReport {
  std::vector<Issue> issues;
  bool ok() const;
  std::string text() const;
};

/// The embedded JSON Schema (version 1) as text.
const std::string& schema_v1_text();

/// Upgrade a document to the current schema version in place. Throws on unknown versions.
void migrate(Node& doc);

/// Parse, migrate, validate (schema and semantic rules) and type a description.
/// Never throws for validation problems: they are returned in the report.
ValidationReport validate_description(const std::filesystem::path& file, Description* out = nullptr);
ValidationReport validate_description_text(std::string_view text, const std::filesystem::path& base_dir,
                                           Description* out = nullptr);

/// As validate_description, but throws Error(Validation) with the report text if invalid.
Description load_description(const std::filesystem::path& file);

}  // namespace stratonix::desc
