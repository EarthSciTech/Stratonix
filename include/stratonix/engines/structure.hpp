/// @file structure.hpp
/// @brief Unfaulted structure: stratigraphic column, fold primitives, imported surfaces with
/// isochore stacking, and the well-correction field. z_h(x, y) is the elevation of horizon h in
/// the restored (unfaulted) frame.
#pragma once
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "stratonix/core/grid.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"

namespace stratonix {

/// A fold primitive. shape() is 1 at the crest and 0 at the edge; rim() is the rim-syncline
/// shape (0..1).
struct FoldPrimitive {
  enum class Type { Anticline, Dome, Syncline, Basin };
  Type type = Type::Anticline;
  std::string name;
  double cx = 0, cy = 0;          ///< crest
  double azimuth = 0;             ///< axis azimuth (rad, clockwise from north)
  double plunge = 0;              ///< rad; tilt along the axis applied to every horizon
  double half_length = 0;         ///< 0 = cylindrical
  double w_left = 0, w_right = 0; ///< limb half-widths (anticline/syncline)
  double ra = 0, rb = 0;          ///< radii (dome/basin)
  double rim_depth = 0, rim_width = 0;
  std::vector<double> relief;     ///< per horizon (m), sign included
  double relief_max = 0;          ///< max |relief| (for rim scaling)
  /// Raised-cosine bump g(t; w) = (1 + cos(pi t / w)) / 2 for |t| < w.
  static double bump(double t, double w);
  void axis_coords(double x, double y, double& u, double& v) const;
  double shape(double x, double y) const;
  double rim(double x, double y) const;
  double tilt(double x, double y) const;  ///< -tan(plunge) * u
  /// Contribution to horizon h at (x, y).
  double contribution(std::size_t h, double x, double y) const;
};

/// Unfaulted horizon functions.
class Structure {
 public:
  /// Build from a validated description and its input grids. Honours the reference well by
  /// construction and other control wells by a correction field (reported in the log).
  Structure(const desc::Description& d, const BuildInputs& in, BuildLog& log);

  std::size_t horizons() const { return kind_.size(); }
  /// Elevation of horizon h at (x, y) in the restored frame.
  double z(std::size_t h, double x, double y) const;
  /// All horizons at (x, y) (out has horizons() entries).
  void z_all(double x, double y, double* out) const;
  /// Correction applied at a point (for reports).
  double correction(std::size_t h, double x, double y) const;
  const std::vector<FoldPrimitive>& primitives() const { return prims_; }
  /// Maximum |correction| at the wells' neighbourhood per horizon (reported).
  const std::vector<double>& max_correction() const { return max_corr_; }
  /// Reference depths (positive down) at the reference location.
  const std::vector<double>& reference_depths() const { return ref_depth_; }

 private:
  enum class Kind { Parametric, Surface, Isochore };
  double z_uncorrected(std::size_t h, double x, double y) const;
  double regional(double x, double y) const;

  std::vector<Kind> kind_;
  std::vector<double> zref_;           ///< constant term per horizon (parametric)
  std::vector<double> ref_depth_;
  std::vector<const Grid2D*> surf_;    ///< surface grid (Surface kind)
  std::vector<double> surf_sign_;      ///< -1 if grid holds depth, +1 if elevation
  std::vector<int> iso_from_;
  std::vector<double> iso_const_;
  std::vector<const Grid2D*> iso_grid_;
  std::vector<FoldPrimitive> prims_;
  double reg_gx_ = 0, reg_gy_ = 0, reg_x0_ = 0, reg_y0_ = 0;
  // correction field: per horizon, Gaussian RBF weights at well locations
  struct Corr { std::vector<double> x, y, w; };
  std::vector<Corr> corr_;
  double corr_range_ = 0;
  std::vector<double> max_corr_;
  std::map<std::string, Grid2D> owned_;  ///< copies of input grids
};

}  // namespace stratonix
