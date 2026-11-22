/// @file checks.hpp
/// @brief Checks that fail a build, each with a stated tolerance, plus closure reporting.
#pragma once
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "stratonix/core/well.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/geophysics.hpp"
#include "stratonix/model/model.hpp"
#include "stratonix/registry/registry.hpp"
#include "stratonix/repr/faultrock.hpp"
#include "stratonix/repr/grid.hpp"
#include "stratonix/repr/repr.hpp"
#include "stratonix/repr/topology.hpp"

namespace stratonix {

/// One check outcome.
struct CheckResult {
  std::string id;          ///< e.g. "C1"
  std::string name;        ///< plain-language name
  std::string subject;     ///< what was checked (horizon, fault, well, block)
  double value = 0;        ///< measured value
  double tolerance = 0;    ///< tolerance applied
  std::string unit = "m";
  std::string tolerance_source;  ///< description field or "Stratonix default"
  bool passed = true;
  bool warn_only = false;  ///< reported but does not fail the build
  bool info = false;       ///< report-only line (closure numbers without targets)
  std::string detail;
};

struct CheckReport {
  std::vector<CheckResult> results;
  bool passed() const;     ///< all non-warning checks passed
  std::string text() const;
  std::string csv() const;
  std::size_t failures() const;
};

/// Grids and properties of a description's grid block (stage 3B, STX-TN-12), computed once and
/// shared by the checks (C23 to C25), the exports and the figures.
struct GridProducts {
  repr::GridOptions options;
  repr::CornerPointGrid grid;
  repr::PropertySet props;
  bool has_trans = false;
  repr::TransmissibilityResult trans;
  bool has_tets = false;
  repr::TetGrid tets;
  repr::GridDefects defects;
  double fault_tolerance = 0;              ///< C25 tolerance (m)
  std::vector<std::string> fault_names;    ///< by fault index
};
/// Model unit index of a description unit (by its top and base horizons), or -1.
int model_unit_of(const desc::Description& d, const Model& m, const std::string& unit);
repr::GridOptions grid_options(const desc::Description& d, const Model& m);
repr::PropertyInputs property_inputs(const desc::Description& d, const Model& m);
/// Everything the grid block asks for; requires d.grid.present.
GridProducts make_grid_products(const desc::Description& d, const Model& m);

/// Shared, lazily computed representations used by checks (and reused by exports/figures).
class CheckContext {
 public:
  CheckContext(const Model& m, const desc::Description& d, std::vector<WellPath> wells, double grid_spacing);
  const Model& model() const { return m_; }
  const desc::Description& description() const { return d_; }
  /// False when checking a state of a timeline (C17 and C18 then belong to the full model).
  bool timeline_checks = true;
  const std::vector<WellPath>& wells() const { return wells_; }
  const GridSpec& grid() const { return grid_; }
  const repr::BlockList& blocks();
  const Grid2D& block_surface(int horizon, BlockCode b);
  const repr::FaultedSurface& faulted_surface(int horizon);
  const repr::FaultSurface& fault_lattice(int fault);           ///< extended lattice for checks
  const std::vector<repr::CutoffLine>& cutoffs(int fault, int horizon);
  const repr::ThrowRow& throw_row(int fault, int horizon);
  /// Tolerance from the description or the built-in default; `source` receives where it came from.
  double tolerance(const std::string& key, std::string* source) const;
  bool warn_only(const std::string& check_id) const;
  double fault_spacing() const { return fault_spacing_; }
  /// Grid products of the description's grid block (nullptr without one); computed once.
  const GridProducts* grid_products();
  /// The faulted surfaces of this model, shared with the figures and exports of the build.
  repr::FaultedSurfaceCache* surfaces() { return surfaces_.get(); }
  /// Stage 4B: add a cell property computed by the pipeline (P32, P32_DFN) before the checks run.
  void add_grid_property(repr::CellProperty p);
  /// Stage 4A: station misfits of the geophysics block (C29), set by the pipeline; nullptr = none.
  const PotentialFieldResult* potential_fields = nullptr;
  /// Stage 4B: the volumetric restoration of the restoration block (C31), set by the pipeline.
  const struct VolumetricResult* volumetric = nullptr;
  /// Stage 4B: heat, pressure and fracture-network results for C32 to C34, set by the pipeline.
  struct BasinResults {
    std::vector<std::pair<std::string, double>> heat_balance;  ///< location, relative error of the steady present column
    struct Bounds { std::string location; int below = 0, above = 0; double worst = 0; };  ///< Pa outside the bounds
    std::vector<Bounds> pressure;
    bool dfn = false, dfn_ok = false;
    std::string dfn_unit, dfn_reason;
    double dfn_target = 0, dfn_realised = 0;
    std::size_t dfn_count = 0;
  };
  const BasinResults* basin = nullptr;
  /// Stage 4C (STX-TN-16 section 6): well data (C35) and age consistency (C38), set by the pipeline.
  struct WellDataResults {
    struct Well {
      std::string well;
      std::string survey;                 ///< method and azimuth reference, or how the path was given
      std::size_t logs = 0, curves = 0;
      double beyond = 0;                  ///< m of log index beyond the well path's total depth
      double tops_off_path = 0;           ///< largest distance of a top below the path's total depth (m)
      std::vector<std::string> notes;     ///< nulls and units per curve
    };
    std::vector<Well> wells;
    struct Conflict { std::string well, what; };
    bool ages = false;
    std::size_t dated = 0;                ///< age data used
    std::vector<Conflict> conflicts;
  };
  const WellDataResults* well_data = nullptr;
  /// Stage 4C: time-depth ties after conversion (C36), set by the pipeline.
  struct TimeDepthResults {
    struct Tie { std::string well, layer; double residual = 0, loo = 0; bool excluded = false; };
    std::vector<Tie> ties;
    double tolerance = 5;
    std::string tolerance_source;
  };
  const TimeDepthResults* time_depth = nullptr;
  /// Stage 4D: geostatistical properties (C39 conditioning, C40 statistics), set by the pipeline.
  struct GeostatResults {
    struct Conditioning { std::string property; std::size_t cells = 0, mismatched = 0; double largest = 0; bool exact = true; };
    std::vector<Conditioning> conditioning;
    struct Statistic { std::string property, what; double target = 0, got = 0, tolerance = 0; };
    std::vector<Statistic> statistics;
    std::size_t cells = 0, fallback = 0;
  };
  const GeostatResults* geostat = nullptr;
  /// Stage 4A: the contact graph (C26, topology.json), computed once at topology_spacing().
  const ContactGraph& contact_graph();
  /// outputs.voxel_spacing (x) when given, else max(surface spacing, largest extent / 64).
  double topology_spacing() const;

 private:
  const Model& m_;
  const desc::Description& d_;
  std::vector<WellPath> wells_;
  GridSpec grid_;
  double fault_spacing_;
  std::unique_ptr<repr::BlockList> blocks_;
  std::map<std::pair<int, BlockCode>, Grid2D> block_surf_;
  std::shared_ptr<repr::FaultedSurfaceCache> surfaces_;
  std::map<int, repr::FaultSurface> lattice_;
  std::map<std::pair<int, int>, std::vector<repr::CutoffLine>> cutoffs_;
  std::map<std::pair<int, int>, repr::ThrowRow> throws_;
  std::unique_ptr<GridProducts> grid_products_;
  std::unique_ptr<ContactGraph> graph_;
};

/// Check plugin interface.
class Check : public registry::PluginBase {
 public:
  virtual std::vector<CheckResult> run(CheckContext& ctx) const = 0;
};

/// Run every registered check (in name order) and return the report.
CheckReport run_checks(CheckContext& ctx);

/// Line-length and area balance of a section (stage 2E, STX-TN-09 section 3): each horizon's bed
/// length and each unit's area in the present section against the same material restored (faults
/// undone by their own inverses). One row per horizon ("line length") and per unit ("area").
struct BalanceRow {
  std::string section, subject, method;
  double present = 0, restored = 0, mismatch = 0;  ///< mismatch: |restored - present| / present
  bool tested = true;  ///< false: reported only (a syn- or post-kinematic unit's area)
};
std::vector<BalanceRow> balance_rows(const Model& m, const desc::Description& d);

/// Register this layer's plugins.
void register_check_plugins();

/// Fault-rock and stress inputs from a description (stage 2F, STX-TN-10), in SI units.
repr::FaultRockInput fault_rock_input(const desc::Description& d, const Model& m);
repr::StressInput stress_input(const desc::Description& d);
/// With the description's mechanical stratigraphy (stage 4B) mapped onto the model's units.
repr::StressInput stress_input(const desc::Description& d, const Model& m);

}  // namespace stratonix
