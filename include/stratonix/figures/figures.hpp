/// @file figures.hpp
/// @brief Figure set: sections along polylines, map views per horizon, 3D views, fence diagram.
/// Built only from the Model interface and representations. Every figure is written as PNG
/// (300 dpi), SVG (text as text), TIFF (300 dpi, LZW, RGB, no alpha) and PPTX (one 16:9 slide
/// with title, caption and notes), named by the figure naming rule, with a footer stating the
/// coordinate system, vertical datum and evidence basis.
#pragma once
#include <filesystem>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "stratonix/core/well.hpp"
#include "stratonix/figures/canvas.hpp"
#include "stratonix/model/model.hpp"
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
  double grid_spacing = 50;
  int closure_horizon = -1;
  bool faults_seal = false;
};

/// Naming rule: STX_<model>_v<version>_<kind>_<name> (letters, digits, '-', '_' only).
std::string figure_basename(const FigureContext& c, const std::string& kind, const std::string& name);

Color unit_color(const FigureContext& c, std::size_t unit);

Figure section_figure(const FigureContext& c, const repr::Section& s, const std::vector<double>* detail_window);
Figure map_figure(const FigureContext& c, int horizon, const repr::FaultedSurface& surface,
                  const std::vector<std::vector<repr::CutoffLine>>& cutoffs, const repr::Closure* closure);
Figure view3d_figure(const FigureContext& c, double azimuth_deg, double elevation_deg, const std::string& preset_name);
Figure fence_figure(const FigureContext& c, const std::vector<repr::Section>& legs);

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
  /// block_highest_at_fault, relief_to_end
  std::string kind, label;
  int horizon = -1, unit = -1, fault = -1;
  std::vector<int> faults;
  std::string well, end;        ///< end: "first" or "last" (relief_to_end)
  double fraction = 0.5;        ///< well_fault_distance: level inside the unit
  int decimals = 1;
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
  std::string summary_name;     ///< file name shown in the key-numbers heading
  std::vector<SheetCompletion> completions;
  /// Unfaulted elevation of horizon h at (x, y), when the engine provides it (else empty).
  std::function<double(int, double, double)> unfaulted;
};
struct SheetValue { std::string label, text; double value = 0; };
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
