/// @file faults.hpp
/// @brief The fault network: planar faults with horizon-relative tips, throw profiles, slip modes
/// (kinematics), hanging-wall displacement fields, ages, abutting and the frame rule.
///
/// The network is its own object. Once its tips are resolved it no longer needs the horizons,
/// so another engine (v2 implicit) can reuse it unchanged.
#pragma once
#include <cmath>
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
  /// Slip vector given the strike direction s and the rake (radians, Aki and Richards). Dip-slip
  /// kinematics ignore s and rake; oblique and strike-slip return a unit vector, and the fault's
  /// displacement magnitude is then its net slip.
  virtual Vec3 slip_vector(const Vec3& s, const Vec3& dd, double dip, double rake) const {
    (void)s;
    (void)rake;
    return slip_vector(dd, dip);
  }
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

/// Curved fault (listric or ramp-flat) with a vertical- or inclined-shear hanging wall
/// (STX-TN-05 section 2). Coordinates in the dip plane: u, horizontal distance from the reference
/// point along the dip direction; z, depth (positive down). Shear lines are u = xi + s z tan(alpha).
struct CurvedFault {
  enum class Shape { Listric, RampFlat };
  Shape shape = Shape::Listric;
  double z0 = 0, tan0 = 0, zd = 0;   ///< listric: reference depth, tan(dip at z0), detachment depth
  std::vector<double> pu, pz;        ///< ramp-flat profile points (u increasing, depth not decreasing)
  double s = -1;                     ///< shear sense: -1 antithetic, +1 synthetic
  double ta = 0;                     ///< tan(shear angle); 0 = vertical shear
  double heave = 0;                  ///< maximum heave (m); fault-parallel flow: signed slip (+ down dip)
  /// Fault-parallel flow (STX-TN-06 section 2): material moves along flow lines at constant
  /// offset d from the polyline (pu, pz), meeting on the bisectors of its bends.
  bool fpf = false;
  std::vector<double> fe_u, fe_z, fn_u, fn_z, fm_u, fm_z;  ///< unit tangents, normals, miters (per vertex)
  std::vector<double> fl0, fl1;      ///< segment length along a flow line: fl0 + d * fl1
  double fpf_dmax = 0;               ///< offset where neighbouring bisectors cross
  void setup_fpf();
  /// Flow-line coordinates (s, d) of (u, z), and the inverse.
  void fpf_coords(double u, double z, double& s, double& d) const;
  void fpf_at(double s, double d, double& u, double& z) const;
  /// Depth of the fault at u, and its slope dz/du (the dip is atan of it).
  double zf(double u) const;
  double dzdu(double u) const;
  /// u where the fault reaches depth z.
  double u_at(double z) const;
  /// Depth where the shear line through xi meets the fault.
  double G(double xi) const;
  /// Forward: material at (u, z) on the hanging-wall side, heave h -> present (u1, z1).
  void forward(double u, double z, double h, double& u1, double& z1) const;
  /// Inverse: present (u1, z1) on the hanging-wall side, heave h -> restored (u, z).
  void inverse(double u1, double z1, double h, double& u, double& z) const;
};

/// Trishear fold ahead of a propagating tip (STX-TN-06 section 1): symmetric linear velocity field
/// of Zehnder and Allmendinger (2000).
struct TrishearFold {
  double tanphi = 0;   ///< tan of the half apical angle
  double ps = 0;       ///< propagation-to-slip ratio
  double v0 = 1;       ///< +1: hanging wall moves towards the tip (reverse); -1: away from it (normal)
};

/// Rule placing a fault tip relative to the stratigraphy.
struct TipRule {
  enum class Kind { Horizon, WithinUnit, Depth, Free };
  Kind kind = Kind::Free;
  int horizon = -1, top = -1, base = -1;
  double offset = 0;    ///< m, positive downwards (Horizon)
  double fraction = 0;  ///< WithinUnit
  double depth = 0;     ///< m, positive down (Depth)
  bool unfaulted = false;  ///< level read on the unfaulted stratigraphy instead of the fault's frame
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
  /// Hanging wall "layer_rigid": every layer moves as a rigid slab by the throw where that layer
  /// meets the fault (in the fault's frame), so the displacement of a point depends on its
  /// stratigraphic level, not on its elevation.
  bool layered = false;
  std::vector<double> la_tab;              ///< strike positions of the layer table
  std::vector<std::size_t> l_start;        ///< per strike position: first entry in ls_tab / lb_tab
  std::vector<double> ls_tab, lb_tab;      ///< (level, down-dip position) pairs on the plane, level increasing
  /// Down-dip position where stratigraphic level s (horizon index, fractional inside units;
  /// -1 domain top, nh domain base) meets the plane at strike position a.
  double layer_b(double s, double a) const;
  /// Growth strata (STX-TN-07 section 4): fraction of the displacement each horizon carries
  /// (empty: 1 everywhere); between horizons linear in the level s.
  std::vector<double> growth_f;
  double growth_at(double s) const;

  /// Curved geometry (null for a planar fault). For a curved fault the down-dip coordinate b of
  /// at(), tips(), ab() and profile_coords() is the horizontal distance u along the dip direction.
  std::shared_ptr<const CurvedFault> curved;
  std::string hw_name;          ///< hanging-wall plugin name
  // Interaction (STX-TN-05 sections 3 and 4).
  int drag = 0;                 ///< 0 none, +1 normal drag, -1 reverse drag
  double drag_fraction = 0, drag_width = 0;
  bool transfer = false;        ///< linked splay of the abutted fault
  int relay_with = -1;          ///< soft-linked relay partner, or -1
  double rake = 0;              ///< oblique and strike-slip faults (radians)
  std::shared_ptr<const TrishearFold> tri;  ///< trishear fold kinematics (null: none)
  /// Slip along the fault at strike position a (trishear: constant behind the tip).
  double slip_at(double a) const;
  /// Relay window in this fault's strike coordinate: the strike factor of the throw is replaced by
  /// a linear ramp from 1 at ra_one to 0 at ra_zero (complementary on the partner).
  double ra_one = 0, ra_zero = 0;
  /// Drag shape at signed distance sigma: f (1 - |sigma| / W)^2 inside W, 0 beyond.
  double drag_shape(double sigma) const {
    double a = std::fabs(sigma);
    return a < drag_width ? drag_fraction * (1 - a / drag_width) * (1 - a / drag_width) : 0.0;
  }
  /// Heave of a curved fault at strike position a (tapered by the strike profile).
  double heave_at(double a) const;

  /// Signed distance to the fault (+ hanging wall). Exact for a plane; for a curved fault the
  /// vertical separation times the cosine of the local dip (exact in sign, first order in size).
  double sigma(const Vec3& p) const {
    if (!curved) return dot(p - p0, n);
    double u = dot(p - p0, dh), t = curved->dzdu(u);
    return (curved->zf(u) + p.z) / std::sqrt(1 + t * t);
  }
  /// Unit normal (towards the hanging wall) at plane coordinates (a, b).
  Vec3 normal_at(double b) const {
    if (!curved) return n;
    double t = curved->dzdu(b), c = 1 / std::sqrt(1 + t * t);
    return dh * (t * c) + Vec3{0, 0, c};
  }
  void ab(const Vec3& p, double& a, double& b) const {
    Vec3 d = p - p0;
    a = dot(d, s);
    b = curved ? dot(d, dh) : dot(d, dd);
  }
  /// Profile coordinates of a point: strike position, and the down-dip position of the fault
  /// point at the same elevation. A hanging-wall point therefore takes the displacement of the
  /// fault where its own elevation meets the fault (equal to ab() on the plane).
  void profile_coords(const Vec3& p, double& a, double& b) const {
    a = dot(p - p0, s);
    b = curved ? dot(p - p0, dh) : (p0.z - p.z) / std::sin(dip);
  }
  Vec3 at(double a, double b) const {
    if (!curved) return p0 + s * a + dd * b;
    Vec3 q = p0 + s * a + dh * b;
    q.z = -curved->zf(b);
    return q;
  }
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
  /// Drag (STX-TN-05 section 3). Forward drag displacement of material at its pre-drag position
  /// p1 on the given wall, and its inverse for a present position q.
  Vec3 drag_forward(int f, const Vec3& p1, bool hanging_wall) const;
  Vec3 drag_inverse(int f, const Vec3& q, bool hanging_wall) const;
  /// Trishear: material position mapped through the slip history of fault f (forward from the
  /// restored state to the present, or back).
  Vec3 trishear_map(int f, const Vec3& q, bool forward) const;
  /// One step of restoration through fault f: drag on either wall, then the fault itself when q is
  /// on the moving side.
  Vec3 step(int f, const Vec3& q, bool side, bool fast) const;
  /// Ages, youngest first.
  const std::vector<int>& young_to_old() const { return order_; }
  /// Layer-rigid networks: displacement magnitude of fault f for the material at present
  /// position p (its level found by restoring p through every fault).
  double layer_displacement(int f, const Vec3& p, bool* active = nullptr) const;
  /// True when the faults use layer-rigid hanging walls (all of them do, or none).
  bool layered() const { return layered_; }
  /// Stratigraphic coordinate of a restored point: horizon index, fractional inside a unit
  /// (-1 at the domain top, nh at the domain base, linear in metres / 1000 beyond them).
  double strat_coord(const Vec3& r) const;

 private:
  /// Layer-rigid restoration of q through the faults of ranks [r0, r1) for material at level s.
  /// first: 0 evaluate the side of the first fault; +1 / -1 move it as hanging wall / footwall
  /// (undo semantics); +2 / -2 take its side as hanging wall / footwall and respect the frame.
  Vec3 chain(const Vec3& q, std::size_t r0, std::size_t r1, double s, BlockCode* code, const BlockCode* forced,
             int first) const;
  /// Level s of the material at q such that restoring it through ranks [r0, end) lands on level s.
  double solve_layer(const Vec3& q, std::size_t r0, const BlockCode* forced, int first) const;
  double layer_throw(int f, double s, const Vec3& q) const;
  void resolve_layers(const Structure& st, BuildLog& log);
  bool layered_ = false;
  bool vertical_ = false;  ///< every slip vector is vertical (columns do not move when restored)
  const Structure* st_ = nullptr;
  double z_top_ = 0, z_base_ = 0;
  double next_fault_distance(int f, const Vec3& q, double sigma) const;
  double level_in_frame(const Structure& st, int h, double x, double y, int f, const Domain& dom) const;
  void resolve_tips(const desc::Description& d, const Structure& st, BuildLog& log);
  /// Curved fault k: tip tables spanning the model, and the per-horizon report (heave, throw,
  /// rollover width and dip) at the reference strike position.
  void resolve_curved(const desc::Description& d, const Structure& st, BuildLog& log, int k, std::size_t rank);
  std::vector<Fault> faults_;
  std::vector<int> order_;       ///< youngest first
  std::vector<int> rank_;        ///< position of fault f in order_
  Vec3 fixed_;
  BlockCode fixed_block_ = 0;
};

/// Register built-in kinematics and hanging-wall plugins (part of register_engine_plugins()).
void register_fault_plugins();

}  // namespace stratonix
