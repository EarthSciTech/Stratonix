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
struct CompletionDesc { std::string name; Param top, base; };  ///< depths positive down
struct WellDesc {
  std::string name, role;
  Param head, datum_elevation, total_depth, stations;
  std::vector<std::pair<std::string, Param>> tops;  ///< horizon name -> depth (positive down)
  std::vector<CompletionDesc> completions;
};
struct HorizonDesc {
  std::string name;
  Param depth, thickness_above;
  std::string surface_input, surface_format, surface_z_is = "depth";
  std::string iso_from, iso_input, iso_format;
  Param iso_thickness;
  bool intended_offset = false;
};
/// A unit between consecutive horizons; top "model_top" / base "model_base" name the units above
/// the first and below the last horizon.
struct UnitDesc { std::string name, top, base, role, color, label; };
struct PrimitiveDesc {
  std::string type, name;
  Param crest, axis_azimuth, plunge, half_length, half_width_left, half_width_right;
  Param radius_major, radius_minor, rim_depth, rim_width;
  std::vector<std::pair<std::string, Param>> relief;  ///< horizon name -> relief
};
/// Tip rule. stratigraphy: "fault_frame" (default) places the tip on the stratigraphy as displaced
/// by older faults; "unfaulted" on the stratigraphy before any faulting.
struct TipDesc { std::string rule; std::string horizon, unit; Param offset, fraction, depth; std::string stratigraphy = "fault_frame"; };
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
  std::string geom_input, geom_input_format;  ///< geometry type "sticks": fault data input
  Param throw_max;
  ProfileDesc profile_dip, profile_strike;
  std::string slip_mode, hw_type;
  Param influence_distance;
  std::string abuts_fault, abuts_side;
  // Curved faults (stage 2A, STX-TN-05 section 2): listric shape, ramp-flat profile, heave and
  // the shear hanging wall.
  Param detachment_depth;          ///< listric: depth of the horizontal detachment
  Param profile;                   ///< ramp_flat: [[u, z], ...] in the dip plane
  Param heave;                     ///< curved faults: heave (displacement parallel to the regional)
  Param shear_angle;               ///< inclined_shear: angle of the shear lines from vertical
  std::string shear_sense;         ///< inclined_shear: antithetic | synthetic
  bool curved() const { return geom_type == "listric" || geom_type == "ramp_flat"; }
};
struct SectionDesc { std::string name; Param polyline, detail_window; };
/// Section sheet (presentation page) settings.
struct SheetAnnotationDesc { std::string type, text, horizon, fault, side = "footwall", well, completion; Param at, text_at; };
struct SheetPanelDesc {
  bool present = false;
  std::string title;
  Param window, vertical_exaggeration, label_at;
  std::vector<std::string> label_units;
  std::vector<SheetAnnotationDesc> annotations;
};
struct SheetKeyDesc {
  std::string kind, label, horizon, unit, well, fault, end = "last";
  std::vector<std::string> faults;
  Param fraction;
  int decimals = -1;
  std::string basis, value, completion, side = "footwall", coordinate = "distance", from, to;
  bool on_sheet = true;
};
struct SheetDesc {
  std::string name, section, title, along_label, depth_label, end_left, end_right, well_legend, completion_legend;
  std::string key_title = "Key numbers";
  SheetPanelDesc main, detail;
  std::vector<std::string> fault_legend, notes;
  std::vector<SheetKeyDesc> keys;
};
/// Data-conditioned (implicit) structure: data inputs and interpolation settings.
struct ImplicitData { std::string input, kind, format, horizon; };
struct ImplicitDesc {
  bool present = false;
  std::vector<ImplicitData> data;
  bool use_well_tops = true;
  std::string kernel = "cubic";
  Param nugget, surface_spacing, iterations, vertical_stretch;
};

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
  Param tol_thickness_min, tol_crossing, tol_offset_outside_tips, tol_throw, tol_wells, tol_unoffset, tol_rollover;
  std::string closure_horizon;
  Param closure_faults_seal, target_closure, target_area, target_spill_depth, target_tolerance;
  std::vector<std::string> warn_only;
  ImplicitDesc implicit;
  // outputs
  std::vector<SectionDesc> sections;
  std::vector<SheetDesc> sheets;
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
