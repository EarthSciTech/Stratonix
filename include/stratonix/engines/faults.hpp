/// @file faults.hpp
/// @brief The fault network: planar faults with horizon-relative tips, throw profiles, slip modes
/// (kinematics), hanging-wall displacement fields, ages, abutting and the frame rule.
///
/// The network is its own object. Once its tips are resolved it no longer needs the horizons,
/// so another engine (v2 implicit) can reuse it unchanged.
#pragma once
#include <memory>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/model/model.hpp"
#include "stratonix/registry/registry.hpp"

namespace stratonix {

class Structure;

enum class ProfileType { Sine, FlatTopped, Elliptical, Linear };
/// Normalised throw profile on t in [0, 1]: zero at both ends, 1 at t = 0.5.
double profile_value(ProfileType p, double taper, double t);
ProfileType profile_from(const std::string& s);

/// Kinematics plugin: direction of hanging-wall displacement per unit (vertical) throw.
class Kinematics : public registry::PluginBase {
 public:
  /// dd: unit down-dip vector; dip in radians. Returned vector has vertical component -1.
  virtual Vec3 slip_vector(const Vec3& dd, double dip) const = 0;
  virtual std::string name() const = 0;
};

/// Hanging-wall field plugin: scales the displacement with distance from the fault.
class HangingWallField : public registry::PluginBase {
 public:
  /// sigma: distance from the plane on the moving side (>= 0); next: distance along the fault
  /// normal to the next fault plane (infinity if none); influence: influence distance.
  virtual double scale(double sigma, double next, double influence) const = 0;
  virtual bool needs_next_fault() const { return false; }
  virtual std::string name() const = 0;
};

/// Rule placing a fault tip relative to the stratigraphy.
struct TipRule {
  enum class Kind { Horizon, WithinUnit, Depth, Free };
  Kind kind = Kind::Free;
  int horizon = -1, top = -1, base = -1;
  double offset = 0;    ///< m, positive downwards (Horizon)
  double fraction = 0;  ///< WithinUnit
  double depth = 0;     ///< m, positive down (Depth)
};

/// One fault with resolved geometry.
struct Fault {
  std::string name;
  int index = 0;
  double age = 0;
  Vec3 p0;                   ///< reference point on the plane
  double strike = 0, dip = 0;
  Vec3 s, dh, dd, n;         ///< strike, horizontal dip direction, down-dip, normal (to hanging wall)
  double length = 0;
  bool elliptical = false;
  TipRule upper, lower;
  std::vector<double> a_tab, bu_tab, bl_tab;  ///< resolved tip tables (b of upper/lower tips)
  std::vector<unsigned char> clamped_tab;     ///< 1 where a tip was clamped to the domain
  double tmax = 0;
  ProfileType pdip = ProfileType::FlatTopped, pstrike = ProfileType::FlatTopped;
  double taper_dip = 0.25, taper_strike = 0.25;
  std::shared_ptr<const Kinematics> kin;
  std::shared_ptr<const HangingWallField> hw;
  double influence = 0;
  int abuts = -1;            ///< index of the older fault it abuts, or -1
  // Cached at build time (optimised kernels): slip vector, bounds of the active region.
  Vec3 slip;                 ///< kinematics slip vector (vertical component -1)
  double b_lo = 0, b_hi = 0; ///< extreme tips over the strike range (profile coordinate b)
  double slip_b = 0;         ///< max change of profile b caused by a displacement of tmax
  double slip_sigma = 0;     ///< max change of sigma caused by a displacement of tmax
  bool decays = false;       ///< hanging-wall field vanishes beyond the influence distance
  double max_slip_gradient = 0;  ///< max |d throw / d slip| over the fault (must be < 1)
  bool abut_hw_side = true;  ///< active on the hanging-wall side of the abutted fault

  double sigma(const Vec3& p) const { return dot(p - p0, n); }
  void ab(const Vec3& p, double& a, double& b) const { Vec3 d = p - p0; a = dot(d, s); b = dot(d, dd); }
  /// Profile coordinates of a point: strike position, and the down-dip position of the fault
  /// point at the same elevation. A hanging-wall point therefore takes the displacement of the
  /// fault where its own elevation meets the fault (equal to ab() on the plane).
  void profile_coords(const Vec3& p, double& a, double& b) const {
    a = dot(p - p0, s);
    b = (p0.z - p.z) / std::sin(dip);
  }
  Vec3 at(double a, double b) const { return p0 + s * a + dd * b; }
  /// Effective tip positions at strike position a (elliptical tip line applied).
  void tips(double a, double& bu, double& bl) const;
  /// Vertical throw (signed, + normal) at plane coordinates; 0 outside the tip line.
  double throw_at(double a, double b) const;
  /// Extent coordinate: <= 1 inside the tip line.
  double extent_coord(double a, double b) const;
  bool upper_free() const { return upper.kind == TipRule::Kind::Free; }
  bool lower_free() const { return lower.kind == TipRule::Kind::Free; }
};

/// Faults, their relations and the frame rule.
class FaultNetwork {
 public:
  FaultNetwork() = default;
  /// Resolve faults from a description against the restored structure.
  FaultNetwork(const desc::Description& d, const Structure& st, BuildLog& log);

  const std::vector<Fault>& faults() const { return faults_; }
  std::size_t size() const { return faults_.size(); }
  bool empty() const { return faults_.empty(); }
  const Vec3& fixed_point() const { return fixed_; }
  BlockCode fixed_block() const { return fixed_block_; }

  /// Restore a present-day point through all faults (youngest first). If `forced` is given,
  /// the side of each fault is taken from it instead of the point's position. Returns the
  /// restored position; writes the block code if `block` is not null.
  Vec3 restore(const Vec3& p, BlockCode* block, const BlockCode* forced = nullptr) const;
  /// Position of p in the frame of fault f (restored through faults younger than f).
  Vec3 to_fault_frame(const Vec3& p, int f) const;
  /// Relative hanging-wall displacement of fault f for material whose restored position is q
  /// (in f's frame), and whether the fault is active there (abutting).
  Vec3 relative_displacement(int f, const Vec3& q, bool* active = nullptr) const;
  /// Undo fault f for a present position q (in f's frame) on its moving side. sgn = +1 when the
  /// hanging wall moves, -1 when the footwall moves. Solves q0 = q - sgn * slip * tau(q0), where
  /// tau = throw x hanging-wall scale at the restored position q0 (bracketed, Illinois method),
  /// so the throw at a horizon is read where the horizon meets the fault before faulting and
  /// does not depend on which block the frame holds fixed. Returns q0; writes tau if asked.
  Vec3 undo(int f, const Vec3& q, int sgn, double* tau = nullptr) const;
  /// True when fault f provably cannot displace material whose present position is q (the
  /// restored position stays outside the active region for any displacement up to tmax). Then
  /// undo() would return q exactly, so skipping it changes nothing.
  bool cannot_move(int f, const Vec3& q) const;
  /// Restore a point given in fault f's frame, on the stated wall of f, through f and the
  /// faults older than it (used to measure throw in the fault's own frame).
  Vec3 restore_from_frame(const Vec3& q, int f, bool hanging_wall) const;
  /// Ages, youngest first.
  const std::vector<int>& young_to_old() const { return order_; }

 private:
  double next_fault_distance(int f, const Vec3& q, double sigma) const;
  double level_in_frame(const Structure& st, int h, double x, double y, int f, const Domain& dom) const;
  void resolve_tips(const desc::Description& d, const Structure& st, BuildLog& log);
  std::vector<Fault> faults_;
  std::vector<int> order_;       ///< youngest first
  std::vector<int> rank_;        ///< position of fault f in order_
  Vec3 fixed_;
  BlockCode fixed_block_ = 0;
};

/// Register built-in kinematics and hanging-wall plugins (part of register_engine_plugins()).
void register_fault_plugins();

}  // namespace stratonix
