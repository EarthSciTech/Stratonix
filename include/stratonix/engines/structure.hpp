/// @file structure.hpp
/// @brief Unfaulted structure: stratigraphic column, fold primitives, imported surfaces with
/// isochore stacking, and the well-correction field. z_h(x, y) is the elevation of horizon h in
/// the restored (unfaulted) frame.
#pragma once
#include <functional>
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
  /// Parallel mechanics (STX-TN-06 section 4): per horizon, the offset of the anchor's surface at
  /// normal distance t below it, tabulated (1-D across a cylindrical axis, else a 2-D grid); null
  /// where the horizon takes relief * shape.
  struct Offset {
    double t = 0;                 ///< normal distance below the anchor
    bool cyl = true;
    double v0 = 0, dv = 1;        ///< cylindrical: table over v
    double x0 = 0, y0 = 0, dx = 1;///< grid origin and spacing
    std::size_t nx = 0, ny = 0;
    std::vector<double> val;      ///< contribution (offset elevation + t)
    double outside = 0;           ///< contribution beyond the table
    double eval(const FoldPrimitive& p, double x, double y) const;
  };
  std::vector<std::shared_ptr<const Offset>> offset;
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
class Structure final : public StratigraphyQuery, public ColumnQuery {
 public:
  void column(double x, double y, double* elevations) const override { z_all(x, y, elevations); }
  bool horizon_absent(int h, double x, double y) const override { return absent(static_cast<std::size_t>(h), x, y); }
  std::string termination_at(int unit, double x, double y, double radius) const override;
  /// Whether series relations or thickness trends are in force (stratigraphy queries apply).
  bool has_stratigraphy() const { return has_rel_ || has_trend_; }
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
  /// Raw (before the series relations) and effective horizons at (x, y), with a flag per
  /// horizon: 1 lowered by an erosional surface, 2 raised by onlap (STX-TN-07 section 1).
  void z_detail(double x, double y, double* raw, double* eff, unsigned char* flags) const;
  /// True when horizon h is clamped by a relation at (x, y): it has no top there.
  bool absent(std::size_t h, double x, double y) const;
  /// Whether series relations, tilts or clinoforms are in force.
  bool has_relations() const { return has_rel_; }
  /// Termination of model unit u ("" when it does not terminate in the model).
  const std::string& unit_termination(std::size_t u) const { return term_[u]; }
  /// Horizons that are the base of an erosional series (unconformity surfaces), youngest first.
  std::vector<std::size_t> erosional_surfaces() const {
    std::vector<std::size_t> out;
    for (auto& S : series_)
      if (S.base == 1) out.push_back(S.last);
    return out;
  }
  /// Series index of horizon h (0 when no series are given).
  int series_of(std::size_t h) const { return hseries_[h]; }
  /// Whether a raw thickness below zero is a designed pinch-out (thickness trend clamped).
  bool trend_clamped(std::size_t h, double x, double y) const;
  /// Age (Ma) of horizon h; NaN without a timeline.
  double horizon_age(std::size_t h) const { return hage_[h]; }
  /// Fraction of an event active from `from` to `to` (Ma) that followed the deposition of horizon
  /// h: 1 - g(a_h), g the completed fraction (STX-TN-09 section 1).
  double event_fraction(std::size_t h, double from, double to) const;
  /// Growth fraction of horizon h for a structure growing from `onset` to `end` (indices).
  double growth_fraction(std::size_t h, int onset, int end) const;
  // ---- stage 2D (STX-TN-08) -----------------------------------------------------------------
  /// Salt thickness at (x, y) in the restored frame (0 without salt).
  double salt_thickness(double x, double y) const;
  /// Solid (grain) thickness of model unit `unit` at (x, y); its present thickness for rigid units.
  double solid_thickness(int unit, double x, double y) const;
  /// Thickness of model unit `unit` at (x, y) moved to top burial depth `top_burial`.
  double decompacted_thickness(int unit, double x, double y, double top_burial) const;
  /// Present-day mean porosity of model unit `unit` at (x, y) (0 for rigid units).
  double porosity(int unit, double x, double y) const;
  /// Burial history at (x, y): one step per horizon of the compacting column reaching the surface
  /// (base first), then the present. Burial depths per horizon, NaN before deposition.
  struct BurialStep {
    double age = kNaN;  ///< Ma when horizon ages are given
    std::string label;  ///< horizon at the surface, or "present"
    std::vector<double> depth;
  };
  std::vector<BurialStep> burial_history(double x, double y) const;
  /// Missing (eroded) section of model unit `unit` (between horizons unit-1 and unit, both in one
  /// series) at (x, y) in the restored frame (STX-TN-07 section 1).
  double missing_section(int unit, double x, double y) const;
  /// Reference depths (positive down) at the reference location.
  const std::vector<double>& reference_depths() const { return ref_depth_; }

 private:
  void add_copies(const desc::PrimitiveDesc& pd, FoldPrimitive p, const desc::Description& d, BuildLog& log);
  std::shared_ptr<const FoldPrimitive::Offset> parallel_offset(const FoldPrimitive& p, std::size_t ha, double t,
                                                               const desc::Description& d);
  std::size_t offset_anchor_ = 0;
  enum class Kind { Parametric, Surface, Isochore, Clinoform, Salt };
  double z_uncorrected(std::size_t h, double x, double y) const;
  double z_raw(std::size_t h, double x, double y) const;
  void z_all_raw(double x, double y, double* out) const;
  void apply_relations(double x, double y, double* z, unsigned char* flags, bool erosion_only) const;
  double iso_thickness(std::size_t h, double x, double y, bool clamp) const;
  double clino_offset(std::size_t h, double x, double y) const;  ///< H P(...) above the series below
  void report_terminations(const desc::Description& d, BuildLog& log);
  // stage 2C (STX-TN-07)
  struct SeriesInfo {
    std::size_t first = 0, last = 0;  ///< youngest and oldest horizon
    int base = 0;                     ///< 0 conformable, 1 erosional, 2 onlap
    bool clino = false, sigmoid = true;
    double height = 0, width = 1, ax = 0, ay = 1, bx = 0, by = 0, rate = 0;
  };
  std::vector<SeriesInfo> series_;
  std::vector<int> hseries_;
  std::vector<double> tilt_gx_, tilt_gy_;
  bool has_tilt_ = false, has_rel_ = false, has_trend_ = false;
  std::vector<int> trend_;              ///< 0 none, 1 linear, 2 radial
  std::vector<double> trend_a_, trend_b_, trend_cx_, trend_cy_;
  std::vector<std::string> term_;
  // stage 2D (STX-TN-08)
  struct SaltBody {
    std::string name;
    int type = 0;                     ///< 0 pillow, 1 diapir, 2 wall
    double cx = 0, cy = 0, ax = 0, ay = 1, L = 0, R = 0, f = 1, A = 0, W = 0, a = 0;
    std::vector<double> fh;           ///< growth fraction of the overburden per horizon
    double dist(double x, double y) const;
    double shape(double d) const;
    double ring(double d) const;
    double contribution(double x, double y, double t0) const { double d = dist(x, y); return A * shape(d) - std::min(t0, a * ring(d)); }
  };
  std::vector<SaltBody> salt_;
  int salt_top_ = -1, salt_base_ = -1;
  bool salt_imported_ = false;
  double salt_t0_ = 0;
  double salt_t(double x, double y, double* contrib) const;  ///< thickness; per-body contributions
  double salt_net_ = 0, salt_vol_ = 0, weld_area_ = 0;
  int comp_mode_ = 0;                 ///< 0 off, 1 present, 2 depositional
  double surface_depth_ = 0;
  std::vector<int> law_;              ///< per model unit: 0 rigid, 1 exponential, 2 table
  std::vector<double> phi0_, decay_;
  std::vector<std::vector<double>> tz_, tp_;
  int comp_base_ = -1;                ///< deepest base horizon of a compacting unit
  std::vector<std::string> hname_;
  std::vector<double> hage_;
  bool column_ = false;               ///< a column pass is needed (relations, salt or compaction)
  bool from_top_ = false;             ///< reference depths from the reference well's (present) tops
  double phi(int u, double z) const;
  double solid(int u, double z1, double z2) const;
  double base_for(int u, double z1, double S) const;
  double top_for(int u, double z2, double S) const;
  void column_pass(double x, double y, double* z, unsigned char* flags) const;
  void report_salt_and_compaction(const desc::Description& d, BuildLog& log);
 public:
  /// Salt balance from the build (C15): net volume change and the bodies' volume (m3); weld area (m2).
  double salt_net_volume() const { return salt_net_; }
  double salt_bodies_volume() const { return salt_vol_; }
  bool salt_parametric() const { return !salt_.empty(); }
  /// Compaction mode: 0 off, 1 present, 2 depositional; the law of model unit u (0 rigid).
  int compaction_mode() const { return comp_mode_; }
  /// Mean porosity of unit u between burial depths d1 < d2 below the surface (stage 3B): 1 minus
  /// the grain thickness of its compaction law over the interval; NaN without a law.
  double mean_porosity(int u, double d1, double d2) const;
  double surface_depth() const { return surface_depth_; }
  int compaction_law(int u) const { return u >= 0 && static_cast<std::size_t>(u) < law_.size() ? law_[static_cast<std::size_t>(u)] : 0; }
  /// Depositional solid thickness of unit u at (x, y) (depositional mode), for C16.
  double depositional_solid(int u, double x, double y) const;
  /// Present thickness after decompaction to the surface and compaction back (C16 round trip).
  double compaction_round_trip(int u, double x, double y) const;
  const std::vector<std::string>& horizon_names() const { return hname_; }
 private:
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
