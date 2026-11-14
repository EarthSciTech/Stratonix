/// @file model.hpp
/// @brief The model interface every engine implements, and the flat-layers reference model.
///
/// Every model, whatever builds it, answers the same questions: which unit and which fault block
/// is at a point, what is the value of each horizon's scalar field there, and (through the
/// representation layer) where its surfaces are. Sections, exports and checks use only this
/// interface, so they work unchanged for the parametric engine (v1) and an implicit engine (v2).
#pragma once
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "stratonix/core/crs.hpp"
#include "stratonix/core/geometry.hpp"

namespace stratonix {

/// Evidence basis of a whole model, written into every exported file.
enum class EvidenceBasis { Design, Interpreted };
inline const char* basis_label(EvidenceBasis b) { return b == EvidenceBasis::Design ? "DESIGN" : "INTERPRETED"; }

/// Fault-block code: bit f set when the point is on the hanging-wall side of fault f.
using BlockCode = std::uint64_t;

struct HorizonInfo {
  std::string name;
  int index = 0;
  bool intended_offset = false;  ///< true if a top/base horizon may be offset on purpose
};
struct UnitInfo {
  std::string name;
  int top_horizon = -1;   ///< -1 for the unit above the top horizon
  int base_horizon = -1;  ///< -1 for the unit below the base horizon
  std::string color;      ///< "#rrggbb" or empty
};
/// Approximate plane and extent of a fault, used by extractors to locate its surface.
struct FaultInfo {
  std::string name;
  int index = 0;
  int age = 0;              ///< larger = younger
  Vec3 ref_point;           ///< a point on the fault
  Vec3 strike_dir;          ///< unit vector along strike (horizontal)
  Vec3 dip_dir;             ///< unit vector down dip, in the plane
  Vec3 normal;              ///< unit normal pointing to the hanging wall
  double half_length = 0;   ///< extent along strike from ref_point (each way)
  double b_min = 0, b_max = 0;  ///< extent down dip from ref_point (b positive downwards)
  std::string slip_mode;    ///< "dip_slip" or "vertical_offset" (informational)
  std::string geometry = "planar";  ///< "planar", "listric" or "ramp_flat"
  /// Curved faults: the profile in the dip plane, (u, depth) with u the horizontal distance from
  /// the reference point along the dip direction; the plane fields above are then its chord.
  std::vector<std::pair<double, double>> profile;
};
/// Description of a model, independent of the engine.
struct ModelInfo {
  std::string name, version, title;
  EvidenceBasis basis = EvidenceBasis::Design;
  std::string engine;  ///< e.g. "parametric", "flat-layers", "implicit-rbf"
  CrsInfo crs;
  Domain domain;
  std::vector<HorizonInfo> horizons;  ///< ordered top to base
  std::vector<UnitInfo> units;        ///< unit k lies between horizons k-1 and k (k = 0: above top)
  std::vector<FaultInfo> faults;
  /// Data the model was conditioned on, one line per input with its evidence class and source
  /// (empty for purely parametric models). Carried into export READMEs and packages.
  std::vector<std::string> provenance;
};

/// Misfit of a data-conditioned model at one data point (capability "data_misfit" returns a
/// const std::vector<DataMisfit>*).
struct DataMisfit {
  std::string input;      ///< input id
  std::string evidence;   ///< evidence class of the input
  std::string label;      ///< horizon (or fault) the datum belongs to
  std::string kind;       ///< "point", "orientation" or "well_top"
  Vec3 p;                 ///< present-day position
  double misfit = 0;      ///< points: vertical distance to the horizon (m, + above); orientations: angle (deg)
};

/// Which fields evaluate() should fill.
enum EvalFields : unsigned { kEvalUnit = 1u, kEvalBlock = 2u, kEvalScalars = 4u, kEvalAll = 7u };

/// Results of evaluate(), structure of arrays. scalar is horizon-major: scalar[h * n + i].
struct EvalResult {
  std::vector<std::int32_t> unit;
  std::vector<BlockCode> block;
  std::vector<double> scalar;
  std::size_t n = 0;
  double s(std::size_t h, std::size_t i) const { return scalar[h * n + i]; }
};

/// Results of evaluate_fault().
struct FaultEval {
  std::vector<double> distance;      ///< signed distance to the present-day fault surface (+ hanging wall)
  std::vector<double> extent;        ///< <= 1 inside the tip line, > 1 outside
  std::vector<double> displacement;  ///< displacement magnitude at the point (m)
};

/// The model interface. Implementations are immutable after construction; all methods are
/// thread-safe and write only to caller-owned outputs.
class Model {
 public:
  virtual ~Model() = default;
  virtual const ModelInfo& info() const = 0;
  /// Evaluate unit index, block code and/or horizon scalar fields at points.
  /// Scalar convention: zero on the horizon, positive above it, approximately metres nearby.
  virtual void evaluate(std::span<const Vec3> pts, unsigned fields, EvalResult& out) const = 0;
  /// Scalar field of one horizon as seen from a given fault block (extended beyond the block).
  virtual void scalar_in_block(std::span<const Vec3> pts, BlockCode block, int horizon,
                               std::span<double> out) const = 0;
  /// Distance to, extent coordinate on, and displacement of fault `fault` at points.
  virtual void evaluate_fault(int fault, std::span<const Vec3> pts, FaultEval& out) const = 0;
  /// Positions in the restored (unfaulted) frame.
  virtual void restore(std::span<const Vec3> pts, std::span<Vec3> out) const = 0;
  /// Optional engine-specific capability by name (nullptr if not provided).
  virtual const void* capability(std::string_view name) const { (void)name; return nullptr; }

  // Convenience helpers built on the interface.
  int unit_at(const Vec3& p) const;
  BlockCode block_at(const Vec3& p) const;
};

using ModelPtr = std::shared_ptr<const Model>;

/// Reference implementation: flat, horizontal, unfaulted layers. It is the simplest model that
/// passes the interface suite, and documents the expected behaviour of every call.
class FlatLayersModel final : public Model {
 public:
  /// horizon_z: elevations of the horizons, top to base (strictly decreasing).
  FlatLayersModel(ModelInfo info, std::vector<double> horizon_z);
  const ModelInfo& info() const override { return info_; }
  void evaluate(std::span<const Vec3> pts, unsigned fields, EvalResult& out) const override;
  void scalar_in_block(std::span<const Vec3> pts, BlockCode block, int horizon,
                       std::span<double> out) const override;
  void evaluate_fault(int fault, std::span<const Vec3> pts, FaultEval& out) const override;
  void restore(std::span<const Vec3> pts, std::span<Vec3> out) const override;

 private:
  ModelInfo info_;
  std::vector<double> z_;
};

/// Unit index from horizon scalars at one point: the number of horizons the point lies below.
inline int unit_from_scalars(const double* s, std::size_t stride, std::size_t nh) {
  int u = 0;
  for (std::size_t h = 0; h < nh; ++h)
    if (s[h * stride] < 0) u = static_cast<int>(h) + 1;
  return u;
}

}  // namespace stratonix
