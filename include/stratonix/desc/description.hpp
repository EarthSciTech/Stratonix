/// @file description.hpp
/// @brief Typed, validated model description (schema version 1) with evidence on every parameter.
#pragma once
#include <cstdint>
#include <array>
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
struct InputDesc {
  std::string id, path, role, sharing, source, evidence;
  std::filesystem::path resolved;
  CrsInfo crs;            ///< stage 4C (F46): the input's CRS when it states one
  bool has_crs = false;
};
struct CompletionDesc { std::string name; Param top, base; };  ///< depths positive down
/// Stage 4C (STX-TN-16 section 1): log files of a well, read as delivered.
struct WellLogDesc { std::string file, frame, index; std::vector<std::string> curves; std::filesystem::path resolved; };
/// Stage 4C: an age datum in a well (TVD below the model datum, positive down, as the tops).
struct WellAgeDesc { std::string kind, name, timescale, source; Param depth, depth_base, age, sigma; };
struct WellDesc {
  std::string name, role;
  Param head, datum_elevation, total_depth, stations;
  /// stage 4C (rule F44): the survey's calculation method and azimuth reference; both empty for a
  /// trajectory written as in v3 (minimum curvature, grid azimuths)
  std::string survey_method, azimuth_reference, survey_file;
  Param declination, convergence;
  std::vector<std::array<double, 3>> survey;  ///< md (m), inclination and azimuth (rad), from stations or the file
  std::vector<WellLogDesc> logs;
  std::vector<WellAgeDesc> ages;
  std::vector<std::pair<std::string, Param>> tops_twt;  ///< stage 4C: two-way times of the tops (calibration)
  std::vector<std::pair<std::string, Param>> tops;  ///< horizon name -> depth (positive down)
  Param top_uncertainty;           ///< stage 3A: uncertainty of the tops (m)
  std::vector<CompletionDesc> completions;
};
struct HorizonDesc {
  std::string name;
  Param depth, thickness_above;
  std::string surface_input, surface_format, surface_z_is = "depth";
  /// stage 4C (STX-TN-16 section 3): one realisation of the pick uncertainty of a surface in time
  struct PickPerturbation { bool present = false; Param sigma, correlation_length; std::uint64_t seed = 0; } surface_perturbation;
  std::string iso_from, iso_input, iso_format;
  Param iso_thickness;
  Param age;                       ///< stage 2D: optional age (Ma) for burial-history axes
  std::string iso_trend;           ///< "" (constant), "linear" or "radial" (STX-TN-07 section 3)
  Param iso_gradient, iso_centre;
  bool intended_offset = false;
};
/// Stratigraphic series (STX-TN-07 section 1): consecutive horizons, youngest first.
struct SeriesDesc {
  std::string name;
  std::vector<std::string> horizons;
  std::string base;                ///< conformable, erosional or onlap
  Param tilt_dip, tilt_azimuth;
  std::string clino_profile;       ///< "" or sigmoidal / oblique (section 5)
  Param clino_height, clino_width, clino_direction, clino_first_break, clino_progradation;
  Param clino_trajectory;          ///< stage 4B: shelf-break trajectory (angle of rise, 0 to 15 degrees)
};
/// Growth strata (STX-TN-07 section 4): movement between onset and end.
struct GrowthDesc { std::string onset, end; bool present() const { return !onset.empty(); } };
/// Event of the timeline (STX-TN-09 section 1).
struct EventDesc {
  std::string name, type;          ///< deposition, fault, fold or salt
  std::string ref;                 ///< horizon, fault, primitive or body
  Param age, from, to, throw_;     ///< ages in Ma; per-phase throw (faults)
  std::string sense;               ///< per-phase sense (faults)
  double start() const { return from.present ? from.si() : age.si(); }  ///< older end (Ma)
  double end() const { return to.present ? to.si() : age.si(); }        ///< younger end (Ma)
};
/// Fault rock and stress (STX-TN-10).
struct FaultRockDesc {
  bool present = false;
  Param shale_threshold, thickness_ratio, seal_a, seal_c;
  std::string seal_type;           ///< "" or sgr_power
};
struct StressDesc {
  bool present = false;
  Param sv_gradient, shmax_ratio, shmin_ratio, shmax_azimuth;  ///< regime
  Param tensor_gradient;           ///< xx, yy, zz, xy, xz, yz (east, north, up)
  Param pore_pressure_gradient, friction, cohesion;
  bool regime() const { return sv_gradient.present || shmax_ratio.present; }
};
/// Salt (STX-TN-08 section 1).
struct SaltBodyDesc {
  std::string type, name;          ///< pillow, diapir or wall
  Param centre, radius, flank_width, height, withdrawal_width, axis_azimuth, half_length;
  GrowthDesc growth;
};
struct SaltDesc {
  std::string top, base;
  Param source_thickness;
  std::vector<SaltBodyDesc> bodies;
  bool present() const { return !top.empty(); }
};
/// Compaction (STX-TN-08 section 3).
struct CompactionUnitDesc {
  std::string unit, law;           ///< exponential or table
  Param surface_porosity, decay, table_depths, table_porosity;
};
struct CompactionDesc {
  std::string geometry;            ///< present or depositional ("" = no compaction)
  std::string burial = "block";    ///< stage 4B: block (maximum burial in its fault block) | restored (stage 2D, column-wise)
  Param surface_depth;
  std::vector<CompactionUnitDesc> units;
  bool present() const { return !geometry.empty(); }
};
/// A unit between consecutive horizons; top "model_top" / base "model_base" name the units above
/// the first and below the last horizon.
struct UnitDesc { std::string name, top, base, role, color, label; Param clay_fraction; };
struct PrimitiveDesc {
  std::string type, name;
  Param crest, axis_azimuth, plunge, half_length, half_width_left, half_width_right;
  Param radius_major, radius_minor, rim_depth, rim_width;
  std::vector<std::pair<std::string, Param>> relief;  ///< horizon name -> relief
  std::string series;              ///< stage 2C: youngest series folded ("" = all)
  GrowthDesc growth;               ///< stage 2C: growth strata
  // Stage 2B (STX-TN-06 sections 3 to 5)
  Param half_width, shortening;    ///< detachment_fold
  std::string detachment;          ///< detachment_fold: horizon of the detachment
  std::string mechanics;           ///< "" (interpolated relief), "parallel" or "similar"
  int array_count = 1;             ///< en echelon copies
  Param array_step_along, array_step_across;
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
  bool curved() const { return geom_type == "listric" || geom_type == "ramp_flat" || (geom_type == "sticks" && sticks_fit == "profile"); }
  // Faults from data (stage 3A, STX-TN-11 section 3).
  std::string sticks_fit = "plane";  ///< geometry.fit for sticks: plane | profile
  bool throw_estimate = false;       ///< throw.estimate: throw, length and tips from horizon cut-offs
  std::string estimate_model = "horizons";  ///< stage 4A: horizons | displacement_field
  Param estimate_band, estimate_spacing;
  // Interaction (STX-TN-05 sections 3 and 4).
  bool abuts_transfer = false;     ///< abuts.transfer: a linked splay of the abutted fault
  std::string drag_type;           ///< normal | reverse (empty: no drag)
  Param drag_fraction, drag_width;
  std::string relay_fault;         ///< relay.fault: soft-linked relay with this fault
  bool relay_breached = false;
  bool linking = false;            ///< stage 4B: connecting fault of a breached relay (no slip of its own)
  // Movement sense (STX-TN-05 section 5).
  std::string sense;               ///< throw.sense: normal | reverse | left_lateral | right_lateral
  Param net_slip;                  ///< throw.slip: net slip of oblique and strike-slip faults
  Param rake;                      ///< oblique: rake (Aki and Richards convention)
  bool strike_component() const { return slip_mode == "oblique" || slip_mode == "strike_slip"; }
  // Fault-related folding (STX-TN-06 section 1).
  GrowthDesc growth;               ///< stage 2C: growth strata (layer_rigid hanging walls)
  std::string fold_type;           ///< trishear (empty: none)
  Param fold_apical_angle, fold_ps_ratio;
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
struct ImplicitData { std::string input, kind, format, horizon; Param uncertainty; };
struct ImplicitDesc {
  bool present = false;
  std::vector<ImplicitData> data;
  bool use_well_tops = true;
  std::string kernel = "cubic";
  Param nugget, surface_spacing, iterations, vertical_stretch;
  // Stage 3A (STX-TN-11): large-data solver and data influence.
  std::string solver_type = "auto";  ///< auto | dense | partition_of_unity
  Param points_per_patch, overlap;
  bool data_influence = false;
  std::vector<std::string> isochore_horizons;  ///< built from the horizon above minus a thickness from their data
  // Stage 4A (STX-TN-14 section 1): free horizon values on the discrete solver.
  std::string horizon_values = "fixed";  ///< fixed | free
  Param discrete_spacing, discrete_vertical_spacing, discrete_regularisation, discrete_orientation_weight, discrete_max_iterations;
  // Stage 4A (STX-TN-14 sections 2 and 3): fold frames and fault relations from data.
  struct FoldFrameDesc { std::string name, series, foliation; Param wavelength, order, constraint_spacing; };
  std::vector<FoldFrameDesc> folds;
  std::string fault_relations = "stated";  ///< stated | infer
  Param relation_distance;
};

/// Hybrid conditioning (stage 3A, STX-TN-11 section 4): parameters of the description fitted to data.
struct FitParameterDesc { std::string path; Param min, max; };
struct ConditioningDesc {
  bool present = false;
  std::vector<ImplicitData> data;
  bool use_well_tops = true;
  Param uncertainty, orientation_uncertainty, max_iterations;
  std::vector<FitParameterDesc> parameters;
  std::string method = "least_squares";  ///< stage 4A: least_squares | bayesian
  Param walkers, steps, burn_in, thin, seed;  ///< stage 4A: conditioning.sampler
};

/// Grids and properties (stage 3B, STX-TN-12 section 10).
struct GridZoneDesc { std::string unit, layering = "proportional"; Param layers, layer_thickness; };
struct GridFaciesDesc { std::string unit, name; Param code; };
struct GridTransformDesc { std::string name, unit, from, law; Param a, b; };
struct GridDesc {
  bool present = false;
  Param spacing, margin, min_thickness;
  std::vector<GridZoneDesc> zones;          ///< empty: every unit between two horizons, one layer
  std::vector<GridFaciesDesc> facies;
  std::string porosity;                     ///< "" or "compaction"
  std::vector<GridTransformDesc> transforms;
  std::string transmissibility;             ///< "" or "manzocchi_1999"
  bool tetrahedral = false;
  std::string pillars = "fault_conforming";  ///< stage 4B: fault_conforming | vertical (dense intersecting networks)
  /// stage 4C: logs upscaled to the cells the wells pass through
  struct Upscaled { std::string curve, method, name; Param min_valid; };
  std::vector<Upscaled> upscaled_logs;
  /// stage 4D (STX-TN-17 section 1): properties simulated by geostatistics
  struct VariogramDesc {
    bool present = false, fit_vertical = false;
    Param nugget;
    struct Structure { std::string type; Param contribution, ranges, azimuth; };  ///< ranges [major, minor, vertical]
    std::vector<Structure> structures;
  };
  struct Geostat {
    std::string name, kind, unit, property_unit, condition_to, trend;
    VariogramDesc variogram;
    struct Code { int code = 0; std::string name; Param proportion; VariogramDesc variogram; };
    std::vector<Code> codes;
    Param min, max, mean, sd, decluster_cell;
    Param search_radii, search_azimuth, max_data;
    std::string object_kind, geobody;
    int object_code = 1;
    Param width, thickness, sinuosity, wavelength, azimuth, radius, proportion, max_iterations;
    Param seed;
  };
  std::vector<Geostat> geostatistics;
};

/// Restoration (stage 4B, STX-TN-15 section 4): volumetric restoration of the whole model to a
/// horizon (decision 18) and the bed-length limit of C31.
struct RestorationUnitDesc { std::string unit; Param young, poisson; };
struct RestorationDesc {
  bool present = false, volumetric = false;
  std::string horizon;
  Param datum, pin, young, poisson, mesh_size, max_strain, bed_length_limit;
  std::vector<RestorationUnitDesc> units;
};
/// Backstripping (STX-TN-15 section 4.4).
struct BackstripStepDesc { std::string horizon; Param water_depth, sea_level; };
struct BackstripDesc {
  bool present = false;
  std::string basement;
  Param mantle_density, water_density, grain_density, te, young, poisson;
  std::vector<BackstripStepDesc> steps;
};
/// Eroded-section estimates from porosity logs (STX-TN-15 section 4.5).
struct ErodedSectionDesc { std::string name, well, unconformity, log; Param surface_porosity, decay, vitrinite_estimate; };
/// Mechanical stratigraphy and tectonic strain (stage 4B, STX-TN-15 section 5.3, I-033).
struct MechanicsUnitDesc { std::string unit; Param young, poisson, ucs, biot; };
struct MechanicsDesc { bool present = false; std::vector<MechanicsUnitDesc> units; Param strain_min, strain_max; };
/// Fracture intensity per unit and discrete fracture networks (STX-TN-15 sections 5.1, 5.2).
struct FractureUnitDesc { std::string unit; Param p32_0, a_k, a_e, p32_max; };
struct FractureSetDesc { std::string name; Param pole, kappa, exponent, rmin, rmax, weight; };  ///< pole: [trend, plunge]
struct DfnDesc { bool present = false; std::string unit; Param area, seed, limit; std::vector<FractureSetDesc> sets; };
struct FracturesDesc { bool present = false; std::vector<FractureUnitDesc> units; DfnDesc dfn; };
/// Heat (STX-TN-15 section 6.1).
struct HeatUnitDesc { std::string unit; Param conductivity, heat_production; };
struct HeatDesc {
  bool present = false;
  Param surface_temperature, basal_heat_flow, water_conductivity, column_depth;
  std::vector<std::pair<Param, Param>> history;  ///< (age, basal heat flow)
  std::vector<HeatUnitDesc> units;
  bool vitrinite = false;
};
/// Overpressure from disequilibrium compaction (STX-TN-15 section 6.2).
struct PressureUnitDesc { std::string unit, law; Param k0, b, phi0; };  ///< exponential: k0 exp(b (phi - phi0)); kozeny_carman: k0 (phi/phi0)^3 ((1-phi0)/(1-phi))^2
struct PressureDesc { bool present = false; Param viscosity, rho_water, rho_grain; std::vector<PressureUnitDesc> units; };
/// Geobodies (stage 4B, STX-TN-15 section 3.3): channel belts and lobes placed inside a unit in
/// geological-time coordinates, so they follow the unit's deformation.
struct GeobodyDesc {
  std::string name, type, unit, color, shape = "rectangular";  ///< type: channels | lobes
  Param count, seed, width, thickness, wavelength, amplitude, orientation, radius, aspect;
};
/// Bodies (stage 4B, STX-TN-15 section 1): salt sheets, sills, saucer sills, dykes, injectites.
struct BodyDesc {
  std::string name, type;          ///< salt_sheet | sill | saucer_sill | dyke | injectite
  Param age;                       ///< emplacement, on the faults' ordinal scale (larger = younger)
  std::string color;
  // salt sheet
  Param footprint;                 ///< [[x, y], ...] m
  std::string base_horizon;        ///< base at a horizon (plus base_offset, positive up) or at base_depth
  Param base_offset, base_depth, thickness, taper;
  // sill and saucer sill (centre, depth of the mid-plane); also thickness
  Param centre, depth, radius, accommodation, uplift_decay, climb, outer_radius;
  // dyke
  Param trace, width, top_depth;   ///< base: base_depth
  // injectite
  std::string parent, ceiling;     ///< parent unit, ceiling horizon
  Param sheets, length, length_sd, thickness_sd, net_to_gross, max_fraction, seed;
};

/// Gravity and magnetics (stage 4A, STX-TN-14 section 7).
struct GeophysicsUnitDesc { std::string unit; Param density, susceptibility; };
struct StationSetDesc { std::string input, format; Param uncertainty; };
struct GeophysicsDesc {
  bool present = false;
  Param cell, background_density, background_susceptibility;
  Param inclination, declination, intensity;  ///< inducing field
  std::vector<GeophysicsUnitDesc> units;
  StationSetDesc gravity, magnetics;
};

/// Uncertainty (stage 3C, STX-TN-13 section 6).
struct UncertainParamDesc { std::string path, distribution; Param mean, sd, min, max, mode, values, weights; };
struct ScenarioOptionDesc { std::string name; Param weight; std::vector<std::pair<std::string, Node>> set; std::vector<std::string> remove; };
struct ScenarioDesc { std::string name; std::vector<ScenarioOptionDesc> options; };
struct ProbabilityMapDesc { std::string kind, unit, horizon; Param depth, spacing; };
struct CorrelationDesc { std::string a, b; Param rho; };
/// stage 4C (STX-TN-16 section 3): pick and velocity uncertainty carried into depth
struct PickUncertaintyDesc { std::string horizon; Param sigma, correlation_length; };
struct VelocityUncertaintyDesc { std::string layer; Param v0_sigma, k_sigma, correlation_length; };
struct UncertaintyDesc {
  bool present = false;
  std::vector<PickUncertaintyDesc> picks;
  std::vector<VelocityUncertaintyDesc> velocity;
  Param members, seed;
  std::vector<UncertainParamDesc> parameters;
  std::vector<CorrelationDesc> correlations;
  std::vector<ScenarioDesc> scenarios;
  bool one_at_a_time = false;
  std::vector<ProbabilityMapDesc> maps;
  std::vector<std::string> key_numbers;
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
  std::string correlation_rule;  ///< stage 4C: well_correlation.rule ("" when not given)
  /// stage 4C (STX-TN-16 section 2.3, decision 23): layer-cake velocity model
  /// one realisation of a layer's velocity uncertainty (STX-TN-16 section 3; written by the ensemble)
  struct VelocityPerturbation { bool present = false; Param v0_sigma, k_sigma, correlation_length; std::uint64_t seed = 0; };
  struct VelocityLayer { std::string name, base; Param v0, k; bool fit_k = false; VelocityPerturbation perturbation; };
  struct VelocityModelDesc {
    bool present = false;
    std::vector<VelocityLayer> layers;
    std::string interpolation = "inverse_distance";
    std::vector<std::pair<std::string, std::string>> exclude;  ///< well, reason
    Param tolerance;
  } velocity;
  /// stage 4C: a SEG-Y volume for sampling and display
  struct SeismicDesc { bool present = false; std::string input, domain; int il = 0, xl = 0, x = 0, y = 0; } seismic;
  /// stage 4C (section 2.5): a model from a geological map
  struct MapSource { std::string input, format, layer, field; };
  struct GeologicalMapDesc {
    bool present = false;
    MapSource units, faults, measurements;
    std::vector<std::pair<std::string, std::string>> unit_order;  ///< map unit -> base horizon, youngest first
    std::string dip_field, direction_field, convention, filter_field, filter_value, source, date;
    std::string dem_input, dem_format;
    Param sample_spacing;
  } geological_map;
  std::string ref_well;
  Param ref_point;
  std::vector<HorizonDesc> horizons;
  std::vector<UnitDesc> units;
  std::vector<SeriesDesc> series;  ///< empty: one conformable series
  CompactionDesc compaction;       ///< stage 2D
  std::vector<EventDesc> events;   ///< stage 2E: empty = implicit timeline
  ConditioningDesc conditioning;   ///< stage 3A
  FaultRockDesc fault_rock;        ///< stage 2F
  GridDesc grid;                   ///< stage 3B
  UncertaintyDesc uncertainty;     ///< stage 3C
  GeophysicsDesc geophysics;       ///< stage 4A
  std::vector<BodyDesc> bodies;    ///< stage 4B
  std::vector<GeobodyDesc> geobodies;  ///< stage 4B: channel belts and lobes inside units
  RestorationDesc restoration;     ///< stage 4B
  BackstripDesc backstripping;     ///< stage 4B
  std::vector<ErodedSectionDesc> eroded_sections;  ///< stage 4B
  MechanicsDesc mechanics;         ///< stage 4B
  FracturesDesc fractures;         ///< stage 4B
  HeatDesc heat;                   ///< stage 4B
  PressureDesc pressure;           ///< stage 4B
  StressDesc stress;               ///< stage 2F
  std::vector<std::string> curvature_horizons;  ///< stage 2F: outputs.curvature.horizons
  Param fracture_threshold;        ///< stage 2F: outputs.curvature.fracture_threshold (1/m)
  SaltDesc salt;                   ///< stage 2D
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
  Param tol_thickness_min, tol_crossing, tol_offset_outside_tips, tol_throw, tol_wells, tol_unoffset, tol_rollover, tol_area,
      tol_salt_volume, tol_solid_thickness, tol_restoration, tol_balance,
      tol_watertight;  ///< stage 4A (C27)
  Param scaling_min_ratio, scaling_max_ratio;  ///< checks.fault_scaling (C10)
  std::string closure_horizon;
  Param closure_faults_seal, target_closure, target_area, target_spill_depth, target_tolerance;
  std::vector<std::string> warn_only;
  ImplicitDesc implicit;
  // outputs
  std::vector<SectionDesc> sections;
  std::vector<SheetDesc> sheets;
  Param vertical_exaggeration, voxel_spacing, fault_stick_spacing;
  std::vector<std::string> exports;
  std::vector<std::string> isopachs;  ///< units with isopach maps besides those that terminate

  bool is_design() const { return model.status == "DESIGN"; }
  int horizon_index(const std::string& name) const;
  int unit_index(const std::string& name) const;
  int fault_index(const std::string& name) const;
  /// Events referring to a fault, primitive or body (stage 2E).
  std::vector<const EventDesc*> events_for(const std::string& type, const std::string& ref) const;
  /// Series of horizon h (index into series), -1 when no series are given.
  int series_of(int h) const;
  const WellDesc* well(const std::string& name) const;
  /// Every parameter present in the description (for the parameters table and manifest).
  std::vector<const Param*> all_params() const;
  /// The parameter named by a path with names (faults.F1.throw.max,
  /// structure.primitives.A1.relief.H1, stratigraphy.horizons.H2.depth, wells.W1.head), or null.
  const Param* param_by_name(const std::string& path) const;
  Param* param_by_name(const std::string& path) { return const_cast<Param*>(static_cast<const Description*>(this)->param_by_name(path)); }
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
