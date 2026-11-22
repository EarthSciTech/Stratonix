/// @file figures.hpp
/// @brief Figure set: sections along polylines, map views per horizon, 3D views, fence diagram.
/// Built only from the Model interface and representations. Every figure is written as PNG
/// (300 dpi), SVG (text as text), TIFF (300 dpi, LZW, RGB, no alpha) and PPTX (one 16:9 slide
/// with title, caption and notes), named by the figure naming rule, with a footer stating the
/// coordinate system, vertical datum and evidence basis.
#pragma once
#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "stratonix/core/well.hpp"
#include "stratonix/figures/canvas.hpp"
#include "stratonix/model/model.hpp"
#include "stratonix/repr/faultrock.hpp"
#include "stratonix/repr/geotime.hpp"
#include "stratonix/repr/grid.hpp"
#include "stratonix/repr/geotime.hpp"
#include "stratonix/repr/repr.hpp"

namespace stratonix::fig {

/// Information every figure needs besides the model.
struct FigureContext {
  const Model* model = nullptr;
  std::string model_name, model_version, title;
  std::vector<WellPath> wells;
  std::vector<std::pair<std::string, std::string>> key_numbers;  ///< read from the summary table
  std::vector<std::string> unit_roles;   ///< per model unit (may be empty)
  std::vector<std::string> unit_colors;  ///< per model unit "#rrggbb" or empty
  std::vector<std::string> unit_labels;  ///< per model unit display name (UTF-8) or empty
  double vertical_exaggeration = 1;
  /// Faulted surfaces of `model` shared with the checks and exports (may be null).
  repr::FaultedSurfaceCache* surfaces = nullptr;
  double grid_spacing = 50;
  int closure_horizon = -1;
  bool faults_seal = false;
  std::vector<int> isopach_units;  ///< units drawn as isopach maps besides those that terminate
  /// Burial history at a location (stage 2D, STX-TN-08 section 4): burial depth of each horizon
  /// per step (NaN before deposition); ages in Ma when every horizon has one.
  struct Burial {
    std::string location;
    std::vector<int> horizons;               ///< model horizon indices drawn
    std::vector<std::string> step_labels;
    std::vector<double> ages;                ///< per step; empty when ages are not given
    std::vector<std::vector<double>> depth;  ///< [step][horizon]
    /// stage 4B: tectonic subsidence per step (Airy; flexural when Te > 0), empty when not computed
    std::vector<double> tectonic_airy, tectonic_flexural;
    /// stage 4B: a value per step and horizon (temperature, overpressure) drawn as colour between
    /// the burial curves; NaN where it is not computed
    struct Field {
      std::string kind, label, unit;           ///< kind names the figure (temperature, overpressure)
      std::vector<std::vector<double>> v;      ///< [step][horizon], the steps of `depth`
    };
    std::vector<Field> fields;
  };
  std::vector<Burial> burial;
  /// Fault rock and stress (stage 2F, STX-TN-10): Allan diagram and fault-surface attributes.
  bool fault_rock = false;
  repr::FaultRockInput rock;
  repr::StressInput stress;
  std::vector<int> curvature_horizons;
  double fracture_threshold = kNaN;  ///< 1/m
};

/// Naming rule: STX_<model>_v<version>_<kind>_<name> (letters, digits, '-', '_' only).
std::string figure_basename(const FigureContext& c, const std::string& kind, const std::string& name);

Color unit_color(const FigureContext& c, std::size_t unit);

Figure section_figure(const FigureContext& c, const repr::Section& s, const std::vector<double>* detail_window);
Figure map_figure(const FigureContext& c, int horizon, const repr::FaultedSurface& surface,
                  const std::vector<std::vector<repr::CutoffLine>>& cutoffs, const repr::Closure* closure,
                  const std::vector<repr::TerminationLine>* terminations = nullptr);
/// Flexural-slip unfolded section (stage 2E): every bed straightened along its own length.
Figure unfolded_figure(const FigureContext& c, const std::string& section, const std::vector<repr::UnfoldedHorizon>& beds);
/// Burial-history curves at one location (stage 2D).
Figure burial_figure(const FigureContext& c, const FigureContext::Burial& b);

/// Stage 4C (STX-TN-16 section 1.6): a well log panel - tracks of curves by depth with tops, age
/// data and the model's units along the well.
struct WellPanel {
  std::string well, index_label;                    ///< e.g. "MD (m)"
  std::vector<double> index;                        ///< depth of each sample
  struct Track { std::string name, unit; std::vector<double> v; };
  std::vector<Track> tracks;
  std::vector<std::pair<std::string, double>> tops; ///< name, index depth
  struct Age { std::string name; double depth = 0, age = 0; };
  std::vector<Age> ages;
  struct Interval { int unit = -1; double top = 0, base = 0; };
  std::vector<Interval> units;                      ///< model units along the well
};
Figure well_log_figure(const FigureContext& c, const WellPanel& p);
/// Wells side by side by TVD below the model datum with the correlation lines, each labelled with
/// the age and rule (or the tops) it rests on.
struct CorrelationPanel {
  std::vector<WellPanel> wells;                     ///< index in TVD (m); one track drawn when present
  std::vector<std::string> horizons, basis;         ///< per line: name and what it rests on
  std::vector<std::vector<double>> depth;           ///< [line][well], NaN where it cannot be placed
  std::string note;                                 ///< the correlation rule, shown under the title
};
Figure correlation_figure(const FigureContext& c, const CorrelationPanel& p);
/// Stage 4C (STX-TN-16 section 2.2): seismic amplitude along a section, in depth through the stated
/// velocity model, with the model's horizons and faults drawn over it. `amp` is s.nz x s.ns, rows
/// from the section's top (NaN outside the volume).
Figure seismic_section_figure(const FigureContext& c, const repr::Section& s, const std::vector<float>& amp, const std::string& how);
/// The burial curves coloured by a field between consecutive horizons and steps (stage 4B).
Figure burial_field_figure(const FigureContext& c, const FigureContext::Burial& b, const FigureContext::Burial::Field& fd);
/// Isopach map of a unit (stage 2C): vertical thickness per fault block, with its termination lines.
Figure isopach_figure(const FigureContext& c, int unit, const Grid2D& thickness, const std::vector<repr::TerminationLine>& terminations);
Figure view3d_figure(const FigureContext& c, double azimuth_deg, double elevation_deg, const std::string& preset_name);
Figure fence_figure(const FigureContext& c, const std::vector<repr::Section>& legs);
/// Throw along strike of a fault system: points (strike distance, throw) per series.
struct StrikeThrowSeries {
  std::string label;
  Color color;
  bool summed = false;
  std::vector<std::pair<double, double>> points;
};
Figure strike_throw_figure(const FigureContext& c, const std::vector<StrikeThrowSeries>& series, const std::string& system_name);
/// Displacement on one fault surface (lattice of fault_surface) with the cut-offs of each horizon.
/// Stage 2F (STX-TN-10): a fault-surface attribute (SGR, CSP, thickness, seal capacity, slip or
/// dilation tendency, reactivation pressure) seen from the hanging wall, with the cut-offs.
/// vmin/vmax NaN: the range of the values.
Figure fault_attribute_figure(const FigureContext& c, int fault, const repr::FaultSurface& fs, const std::vector<double>& values,
                              const std::vector<std::vector<repr::CutoffLine>>& cutoffs, const std::string& attribute,
                              const std::string& unit, double vmin = kNaN, double vmax = kNaN);
/// Stage 3A (STX-TN-11): fitted parameters against their admissible ranges (start, fitted value
/// with +-2 sigma), one row per parameter.
struct FittedRow { std::string path, unit; double start, value, sigma, min, max; };
Figure conditioning_figure(const FigureContext& c, const std::vector<FittedRow>& rows, double chi2_reduced, std::size_t data);
/// Stage 4A (STX-TN-14 section 4): posterior samples for the trace and corner figures, in the
/// parameters' units. chain[walker * steps + step][parameter], after burn-in.
struct PosteriorSamples {
  std::vector<std::string> names, units;
  int walkers = 0, steps = 0, burn_in = 0, thin = 1;
  std::uint64_t seed = 0;
  std::vector<std::vector<double>> chain;
  std::vector<double> mean, q025, q975, rhat, ess;
  std::vector<double> lsq;  ///< least-squares values (NaN if none)
  double acceptance = 0;
  bool converged = false;
};
/// Stage 4A (STX-TN-14 section 7): gravity or magnetic stations, observed against modelled.
struct AnomalyStation { double x = 0, y = 0, observed = 0, modelled = 0, sigma = kNaN; };
/// Three maps: observed, modelled (same colour scale) and residual (modelled - observed).
Figure anomaly_map_figure(const FigureContext& c, const std::string& kind, const std::string& unit, const std::vector<AnomalyStation>& st, double chi2_reduced);
struct AnomalyProfile {
  std::string section;
  std::vector<double> s, modelled;                  ///< along the line
  std::vector<std::pair<double, double>> observed;  ///< (s, value) of stations within the band
  double band = 0;                                  ///< m either side of the line
};
Figure anomaly_profile_figure(const FigureContext& c, const std::string& kind, const std::string& unit, const AnomalyProfile& p);
/// Wheeler diagram (Wheeler 1958): distance along a section against geological time, coloured by
/// unit; hiatuses blank; horizon ages as dashed lines.
Figure wheeler_figure(const FigureContext& c, const WheelerDiagram& w);
/// Trace plot: every walker's path per parameter (at most 16 walkers drawn), with R-hat and ESS.
Figure posterior_trace_figure(const FigureContext& c, const PosteriorSamples& s);
/// Corner plot: marginal histograms with the 95 % interval, and pairwise scatter of the thinned
/// samples (at most 6 parameters).
Figure posterior_corner_figure(const FigureContext& c, const PosteriorSamples& s);
/// Data misfit map of one horizon: data positions coloured by normalised misfit (r / sigma, -3 to
/// 3) where an uncertainty is stated, else by misfit in metres.
struct MisfitPoint { double x, y, misfit, normalised; };
Figure misfit_map_figure(const FigureContext& c, const std::string& horizon, const std::vector<MisfitPoint>& pts);
/// Allan diagram: each cell of the fault surface split into its footwall (upper left) and
/// hanging-wall (lower right) unit colours, with the cut-offs of every horizon.
Figure allan_figure(const FigureContext& c, int fault, const repr::FaultSurface& fs, const std::vector<repr::FaultRockPoint>& pts,
                    const std::vector<std::vector<repr::CutoffLine>>& cutoffs);
/// Curvature map of a horizon: mean curvature (convex up red, concave blue) and fracture-trend ticks.
Figure curvature_figure(const FigureContext& c, int horizon, const repr::CurvatureGrids& cg);
/// Stage 3B (STX-TN-12): a section through the grid along row j (along_x) or column i, cells
/// coloured by zone (prop nullptr) or by a property; fault pillars in red.
Figure grid_section_figure(const FigureContext& c, const repr::CornerPointGrid& g, const repr::CellProperty* prop, bool along_x, int line);
/// A property (or the zone) on layer k, in map view.
Figure grid_layer_figure(const FigureContext& c, const repr::CornerPointGrid& g, const repr::CellProperty* prop, int k);
/// Stage 3C (STX-TN-13): ensemble figures. Each states the members it rests on and the seed.
struct EnsembleLabel { std::size_t members_ok = 0, members_failed = 0; std::uint64_t seed = 0; };
/// Tornado chart of one key number: per parameter, its value at the parameter's P10 and P90.
struct TornadoBar { std::string parameter; double low = 0, high = 0; };
Figure tornado_figure(const FigureContext& c, const std::string& key, double base, const std::vector<TornadoBar>& bars, const EnsembleLabel& l);
/// Histogram of a key number over the ok members, with P10, P50 and P90 marked.
/// Stage 4D (STX-TN-17 section 1.6): experimental variograms (points sized by their pairs) with
/// the model, per direction, and a property's histogram against the wells' upscaled values.
struct VariogramSeries {
  std::string label;
  std::vector<double> lag, gamma;
  std::vector<std::size_t> pairs;
  std::vector<double> model_lag, model_gamma;
};
Figure variogram_figure(const FigureContext& c, const std::string& property, const std::string& unit, double sill, const std::vector<VariogramSeries>& s);
Figure property_histogram_figure(const FigureContext& c, const std::string& property, const std::string& unit, const std::string& value_unit,
                                 const std::vector<double>& simulated, const std::vector<double>& wells, const std::vector<std::string>& codes);
Figure histogram_figure(const FigureContext& c, const std::string& key, const std::vector<double>& values, double p10, double p50, double p90,
                        const EnsembleLabel& l);
/// First-order variance-based indices of one key number, per parameter.
Figure sobol_figure(const FigureContext& c, const std::string& key, const std::vector<std::pair<std::string, double>>& indices, double threshold,
                    const EnsembleLabel& l);
/// Probability map (fraction of the ok members) on a regular lattice.
Figure probability_map_figure(const FigureContext& c, const std::string& what, const Grid2D& p, const EnsembleLabel& l);
/// Stage 4B: a scalar map over the domain (grey where undefined), with a colour bar from vmin to vmax.
Figure scalar_map_figure(const FigureContext& c, const std::string& title, const std::string& subtitle, const std::string& bar_label, const Grid2D& g,
                         double vmin, double vmax);
/// 3D view of the grid's outer faces coloured by zone.
Figure grid_view3d_figure(const FigureContext& c, const repr::CornerPointGrid& g, double azimuth_deg, double elevation_deg);
Figure fault_plane_figure(const FigureContext& c, int fault, const repr::FaultSurface& fs, const std::vector<double>& displacement,
                          const std::vector<std::vector<repr::CutoffLine>>& cutoffs);

struct WrittenFigure {
  std::string basename;
  std::vector<std::string> files;
  std::string title, caption;
};
/// Write a figure in the four formats into `folder`.
WrittenFigure write_figure(const Figure& f, const std::filesystem::path& folder, const std::string& basename,
                           const std::string& title, const std::string& caption, const std::string& notes);

/// Write the whole figure set for a model (sections, maps, 3D views, fence) and a combined PPTX.
std::vector<WrittenFigure> write_figure_set(const FigureContext& c, const std::vector<std::pair<std::string, Polyline>>& sections,
                                            const std::vector<std::vector<double>>& detail_windows,
                                            const std::filesystem::path& folder);

// ---- section sheet ---------------------------------------------------------------------------
/// A labelled annotation on a sheet panel. Text may hold placeholders: {value}, {value:N}
/// (computed depth), {top:N}, {base:N} (completion), N = decimals.
struct SheetAnnotation {
  std::string type;             ///< "leader", "horizon_at_fault" or "completion"
  std::string text;
  int horizon = -1, fault = -1;
  bool footwall = true;         ///< horizon_at_fault: side of the cut-off
  std::string well, completion;
  double at = 0;                ///< leader: distance along the section of the target on the horizon
  double text_s = 0, text_depth = 0;  ///< text anchor (left, first baseline): distance, depth (m, + down)
};
struct SheetPanel {
  bool present = false;
  double s0 = 0, s1 = 0;        ///< distance window along the section (s1 <= s0: whole section)
  double d0 = 0, d1 = 0;        ///< depth window (m, + down; d1 <= d0: model domain)
  double ve = 1;                ///< vertical exaggeration (sets the panel height)
  std::string title;            ///< placeholders {units}, {surfaces}, {faults}, {faults_word}, {ve}, {s0_km}, {s1_km}
  double label_at = -1;         ///< distance where unit labels sit (< 0: none)
  std::vector<int> label_units; ///< model unit indices labelled
  std::vector<SheetAnnotation> annotations;
};
struct SheetKey {
  /// horizon_at_well, horizon_highest, block_offset, fault_tips, fault_throw, well_fault_distance,
  /// block_highest_at_fault, relief_to_end, text, completion, vertical_exaggeration, fold_relief,
  /// limb_dip, fault_dip, fault_tip_units, fault_cutoff_position, fault_cutoff_depths,
  /// fault_overlap, throw_ratio, fault_gap
  std::string kind, label;
  int horizon = -1, unit = -1, fault = -1;
  std::vector<int> faults;
  std::string well, end;        ///< end: "first" or "last" (relief_to_end, fold_relief, limb_dip)
  double fraction = 0.5;        ///< well_fault_distance: level inside the unit
  int decimals = 1;
  std::string basis;            ///< third column of the summary; placeholders {x} {y} {s} {block} {fw_block} {hw_block}
  bool on_sheet = true;         ///< false: summary table only
  std::string value;            ///< text: the value as given
  std::string completion;       ///< completion: name of the well completion
  std::string side = "footwall";       ///< fault_cutoff_position: footwall or hanging_wall
  std::string coordinate = "distance"; ///< fault_cutoff_position: x, y or distance along the section
  int from = -1, to = -1;       ///< throw_ratio: horizons bounding the interval whose thickness divides the throw
};
struct SheetCompletion { std::string well, name; double top = 0, base = 0; };  ///< depths m, + down
/// A presentation sheet for one section: main panel, optional detail panel, legend by role,
/// computed key numbers and notes, laid out on a 15.5 x 9.6 inch page in DejaVu Sans.
struct SheetSpec {
  std::string name, title, along_label, depth_label, end_left, end_right;
  Polyline line;
  SheetPanel main, detail;
  std::string well_legend, completion_legend;
  std::vector<std::string> fault_legend, notes;
  std::vector<SheetKey> keys;
  std::string summary_name;     ///< file name of the summary table
  std::string key_title = "Key numbers";  ///< heading of the key-numbers block
  std::vector<SheetCompletion> completions;
  /// Unfaulted elevation of horizon h at (x, y), when the engine provides it (else empty).
  std::function<double(int, double, double)> unfaulted;
};
struct SheetValue { std::string label, text, basis; double value = 0; bool on_sheet = true; };
struct SheetResult {
  Figure figure{1116, 691.2};
  std::vector<SheetValue> keys;
};
/// Compute the sheet's key numbers and draw it.
SheetResult section_sheet(const FigureContext& c, const SheetSpec& spec);

/// Software renderer used by the 3D views (z-buffer, Lambert shading).
struct Mesh3 {
  TriMesh mesh;
  Color color;
  double alpha = 1.0;  ///< < 1: drawn see-through over the opaque meshes, in list order
};
Image render_meshes(const std::vector<Mesh3>& meshes, double azimuth_deg, double elevation_deg, double vexag, int w, int h,
                    const BBox& box);

}  // namespace stratonix::fig
