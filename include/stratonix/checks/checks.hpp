/// @file checks.hpp
/// @brief Checks that fail a build, each with a stated tolerance, plus closure reporting.
#pragma once
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "stratonix/core/well.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/model/model.hpp"
#include "stratonix/registry/registry.hpp"
#include "stratonix/repr/faultrock.hpp"
#include "stratonix/repr/repr.hpp"

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

 private:
  const Model& m_;
  const desc::Description& d_;
  std::vector<WellPath> wells_;
  GridSpec grid_;
  double fault_spacing_;
  std::unique_ptr<repr::BlockList> blocks_;
  std::map<std::pair<int, BlockCode>, Grid2D> block_surf_;
  std::map<int, repr::FaultedSurface> faulted_;
  std::map<int, repr::FaultSurface> lattice_;
  std::map<std::pair<int, int>, std::vector<repr::CutoffLine>> cutoffs_;
  std::map<std::pair<int, int>, repr::ThrowRow> throws_;
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

}  // namespace stratonix
