/// @file figures.hpp
/// @brief Figure set: sections along polylines, map views per horizon, 3D views, fence diagram.
/// Built only from the Model interface and representations. Every figure is written as PNG
/// (300 dpi), SVG (text as text), TIFF (300 dpi, LZW, RGB, no alpha) and PPTX (one 16:9 slide
/// with title, caption and notes), named by the figure naming rule, with a footer stating the
/// coordinate system, vertical datum and evidence basis.
#pragma once
#include <filesystem>
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

/// Software renderer used by the 3D views (z-buffer, Lambert shading).
struct Mesh3 {
  TriMesh mesh;
  Color color;
  double alpha = 1.0;  ///< < 1: drawn see-through over the opaque meshes, in list order
};
Image render_meshes(const std::vector<Mesh3>& meshes, double azimuth_deg, double elevation_deg, double vexag, int w, int h,
                    const BBox& box);

}  // namespace stratonix::fig
