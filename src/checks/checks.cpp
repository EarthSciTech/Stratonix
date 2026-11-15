// Checks that fail a build. Each check states its tolerance and where the tolerance came from.
#include "stratonix/checks/checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/core/roots.hpp"
#include "stratonix/engines/faults.hpp"
#include "stratonix/engines/structure.hpp"
#include "stratonix/engines/history.hpp"
#include "stratonix/repr/faultrock.hpp"
#include "stratonix/repr/repr.hpp"

// Plugin anchors called below (C linkage, defined by STRATONIX_REGISTER in each plugin file).
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c1);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c2);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c3);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c4);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c5);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c6);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c7);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c8);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c9);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c10);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c11);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c12);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c13);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c14);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c15);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c16);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c17);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c18);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c19);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c20);

namespace stratonix {

// ---- report --------------------------------------------------------------------------------
bool CheckReport::passed() const { return failures() == 0; }
std::size_t CheckReport::failures() const {
  std::size_t n = 0;
  for (auto& r : results) n += (!r.passed && !r.warn_only && !r.info) ? 1 : 0;
  return n;
}
std::string CheckReport::text() const {
  std::ostringstream s;
  s << "Stratonix checks report\n";
  s << "result: " << (passed() ? "PASS" : "FAIL") << " (" << results.size() << " checks, " << failures() << " failed)\n\n";
  for (auto& r : results) {
    std::string st = r.info ? "INFO" : r.passed ? "PASS" : (r.warn_only ? "WARN" : "FAIL");
    s << st << "  " << r.id << " " << r.name << " [" << r.subject << "]: " << fmt_fixed(r.value, 4) << " " << r.unit;
    if (!r.info) s << " (tolerance " << fmt_fixed(r.tolerance, 4) << " " << r.unit << ", " << r.tolerance_source << ")";
    if (!r.detail.empty()) s << " - " << r.detail;
    s << "\n";
  }
  return s.str();
}
std::string CheckReport::csv() const {
  std::string s = "id,name,subject,status,value,unit,tolerance,tolerance_source,detail\n";
  auto q = [](const std::string& x) {
    std::string r = "\"";
    for (char c : x) r += c == '"' ? std::string("\"\"") : std::string(1, c);
    return r + "\"";
  };
  for (auto& r : results) {
    std::string st = r.info ? "INFO" : r.passed ? "PASS" : (r.warn_only ? "WARN" : "FAIL");
    s += r.id + "," + q(r.name) + "," + q(r.subject) + "," + st + "," + fmt_fixed(r.value, 4) + "," + r.unit + "," +
         (r.info ? "" : fmt_fixed(r.tolerance, 4)) + "," + q(r.tolerance_source) + "," + q(r.detail) + "\n";
  }
  return s;
}

// ---- context -------------------------------------------------------------------------------
CheckContext::CheckContext(const Model& m, const desc::Description& d, std::vector<WellPath> wells, double spacing)
    : m_(m), d_(d), wells_(std::move(wells)) {
  const Domain& dom = m.info().domain;
  // Limit check grids to ~160k nodes; the full-resolution grids are used for exports.
  double sp = spacing;
  while ((dom.sx / sp + 1) * (dom.sy / sp + 1) > 160000) sp *= 1.5;
  grid_ = GridSpec::covering(dom, sp);
  grid_.dx = dom.sx / (grid_.nx - 1);
  grid_.dy = dom.sy / (grid_.ny - 1);
  fault_spacing_ = std::max(sp / 2, 10.0);
}
const repr::BlockList& CheckContext::blocks() {
  if (!blocks_) blocks_ = std::make_unique<repr::BlockList>(repr::list_blocks(m_));
  return *blocks_;
}
const Grid2D& CheckContext::block_surface(int h, BlockCode b) {
  auto k = std::make_pair(h, b);
  auto it = block_surf_.find(k);
  if (it == block_surf_.end()) it = block_surf_.emplace(k, repr::horizon_grid_block(m_, h, b, grid_)).first;
  return it->second;
}
const repr::FaultedSurface& CheckContext::faulted_surface(int h) {
  auto it = faulted_.find(h);
  if (it == faulted_.end()) it = faulted_.emplace(h, repr::horizon_grid_faulted(m_, h, grid_)).first;
  return it->second;
}
const repr::FaultSurface& CheckContext::fault_lattice(int f) {
  auto it = lattice_.find(f);
  if (it == lattice_.end()) it = lattice_.emplace(f, repr::fault_surface(m_, f, fault_spacing_, 0.25)).first;
  return it->second;
}
const std::vector<repr::CutoffLine>& CheckContext::cutoffs(int f, int h) {
  auto k = std::make_pair(f, h);
  auto it = cutoffs_.find(k);
  if (it == cutoffs_.end()) it = cutoffs_.emplace(k, repr::cutoff_lines(m_, fault_lattice(f), h)).first;
  return it->second;
}
const repr::ThrowRow& CheckContext::throw_row(int f, int h) {
  auto k = std::make_pair(f, h);
  auto it = throws_.find(k);
  if (it == throws_.end()) it = throws_.emplace(k, repr::throw_row(m_, fault_lattice(f), cutoffs(f, h), h)).first;
  return it->second;
}
double CheckContext::tolerance(const std::string& key, std::string* source) const {
  struct T { const char* key; const desc::Param* p; double def; };
  const T table[] = {{"thickness_min", &d_.tol_thickness_min, 0.5}, {"crossing", &d_.tol_crossing, 0.01},
                     {"offset_outside_tips", &d_.tol_offset_outside_tips, 0.05}, {"throw", &d_.tol_throw, 0.5},
                     {"wells", &d_.tol_wells, 0.01}, {"unoffset", &d_.tol_unoffset, 0.05},
                     {"area", &d_.tol_area, 0.01}, {"salt_volume", &d_.tol_salt_volume, 0.01},
                     {"solid_thickness", &d_.tol_solid_thickness, 0.001}, {"restoration", &d_.tol_restoration, 0.05},
                     {"balance", &d_.tol_balance, 0.02}};
  for (auto& t : table)
    if (key == t.key) {
      if (t.p->present) {
        if (source) *source = "checks.tolerances." + key + " (" + desc::evidence_name(t.p->evidence) + ")";
        return t.p->si();
      }
      if (source) *source = "Stratonix default";
      return t.def;
    }
  throw Error(ErrorCode::Internal, "unknown tolerance " + key);
}
bool CheckContext::warn_only(const std::string& id) const {
  return std::find(d_.warn_only.begin(), d_.warn_only.end(), id) != d_.warn_only.end();
}

// ---- checks --------------------------------------------------------------------------------
namespace {

CheckResult base(CheckContext& c, const std::string& id, const std::string& name, const std::string& subject) {
  CheckResult r;
  r.id = id;
  r.name = name;
  r.subject = subject;
  r.warn_only = c.warn_only(id);
  return r;
}

/// C1/C2: unit thickness per fault block (only where the block is present at the unit).
std::vector<CheckResult> thickness_checks(CheckContext& c, bool crossing) {
  const auto& inf = c.model().info();
  const std::size_t nh = inf.horizons.size();
  std::vector<CheckResult> out;
  std::string src;
  double tol = c.tolerance(crossing ? "crossing" : "thickness_min", &src);
  for (std::size_t h = 0; h + 1 < nh; ++h) {
    double worst = std::numeric_limits<double>::infinity();
    std::string where;
    for (BlockCode b : c.blocks().codes) {
      const Grid2D& a = c.block_surface(static_cast<int>(h), b);
      const Grid2D& z = c.block_surface(static_cast<int>(h + 1), b);
      // nodes where the block is present at mid-unit
      std::vector<Vec3> p;
      std::vector<std::size_t> idx;
      const auto& g = c.grid();
      for (std::size_t k = 0; k < g.size(); ++k) {
        if (std::isnan(a.v[k]) || std::isnan(z.v[k])) continue;
        double zm = 0.5 * (a.v[k] + z.v[k]);
        if (zm > inf.domain.z_top || zm < inf.domain.z_base) continue;
        p.push_back({g.x(k % g.nx), g.y(k / g.nx), zm});
        idx.push_back(k);
      }
      EvalResult r;
      c.model().evaluate(p, kEvalBlock, r);
      const bool terminates = !crossing && !inf.units[h + 1].termination.empty();
      for (std::size_t q = 0; q < idx.size(); ++q) {
        if (r.block[q] != b) continue;
        double t = a.v[idx[q]] - z.v[idx[q]];
        if (terminates && t < 0.01) continue;  // absent by design (stage 2C): C13 checks it
        if (t < worst) {
          worst = t;
          where = "block " + std::to_string(b) + " at (" + fmt_fixed(p[q].x, 0) + ", " + fmt_fixed(p[q].y, 0) + ")";
        }
      }
    }
    if (!std::isfinite(worst)) continue;
    auto r = base(c, crossing ? "C2" : "C1", crossing ? "no horizon crossing" : "positive thickness in every fault block",
                  inf.units[h + 1].name);
    r.value = worst;
    r.tolerance = crossing ? -tol : tol;
    r.tolerance_source = src;
    r.passed = crossing ? worst >= -tol : worst >= tol;
    r.detail = "minimum thickness, " + where;
    if (!crossing && !inf.units[h + 1].termination.empty()) {
      // thins to zero towards its designed termination: positive wherever present (C13 checks the rest)
      r.tolerance = 0;
      r.passed = worst > 0;
      r.detail += "; where present, thinning to its " + inf.units[h + 1].termination + " by design";
    }
    out.push_back(r);
  }
  return out;
}

class C1 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override { return thickness_checks(c, false); }
};
class C2 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override { return thickness_checks(c, true); }
};

/// C3: no offset outside a fault's tips (scalar jump across the plane just outside the tip line).
class C3 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    const auto& inf = c.model().info();
    std::vector<CheckResult> out;
    std::string src;
    double tol = c.tolerance("offset_outside_tips", &src);
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      const auto& fs = c.fault_lattice(static_cast<int>(f));
      std::vector<Vec3> at;
      const Domain& d = inf.domain;
      for (std::size_t k = 0; k < fs.nodes.size(); ++k) {
        const Vec3& q = fs.nodes[k];
        bool in_dom = d.inside_xy(q.x, q.y) && q.z <= d.z_top && q.z >= d.z_base;
        if (!in_dom || fs.extent[k] <= 1.02) continue;  // only just outside the tip line
        at.push_back(q);
      }
      double worst = 0;
      std::string hn;
      if (!at.empty()) {
        auto j = repr::scalar_jumps(c.model(), at, inf.faults[f].normal);
        for (std::size_t h = 0; h < inf.horizons.size(); ++h)
          for (std::size_t k = 0; k < at.size(); ++k)
            if (std::fabs(j[h * at.size() + k]) > worst) { worst = std::fabs(j[h * at.size() + k]); hn = inf.horizons[h].name; }
      }
      auto r = base(c, "C3", "no offset outside the fault's tips", inf.faults[f].name);
      r.value = worst;
      r.tolerance = tol;
      r.tolerance_source = src;
      r.passed = worst <= tol;
      r.detail = std::to_string(at.size()) + " samples just outside the tip line" + (hn.empty() ? "" : "; largest at " + hn);
      out.push_back(r);
    }
    return out;
  }
};

/// C4 (parametric engine): throw at each horizon equals the fault's profile. Measured in the
/// fault's own frame (younger faults undone), where the fault is a plane: on each wall the
/// horizon's cut-off is found down dip from the scalar restored from that frame; the measured
/// throw is the vertical separation of the cut-offs (dip-slip) or the moving cut-off against the
/// fixed wall's horizon at the same map position (vertical offset). Expected: the profile at the
/// restored position of the moving cut-off.
class C4 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    const auto& inf = c.model().info();
    std::vector<CheckResult> out;
    auto* net = static_cast<const FaultNetwork*>(c.model().capability("faults"));
    auto* st = static_cast<const Structure*>(c.model().capability("structure"));
    if (!net || !st) return out;
    std::string src;
    double tol = c.tolerance("throw", &src);
    const double eps = 0.25, step = std::max(5.0, c.fault_spacing() / 4);
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      const Fault& F = net->faults()[f];
      const int fi = static_cast<int>(f);
      if (F.tri || (F.curved && F.curved->fpf)) continue;  // fold kinematics: no throw profile; area is checked by C12
      bool fixed_hw = ((net->fixed_block() >> f) & 1u) != 0;
      double bmin = 1e300, bmax = -1e300;
      for (std::size_t i = 0; i < F.a_tab.size(); ++i) {
        double bu, bl;
        F.tips(F.a_tab[i], bu, bl);
        bmin = std::min(bmin, bu);
        bmax = std::max(bmax, bl);
      }
      // wall scalar on the plane (extrapolated from eps and 2 eps off the plane)
      auto wall_s = [&](std::size_t h, double a, double b, bool hw) {
        double sg = hw ? 1 : -1;
        Vec3 p = F.at(a, b);
        const Vec3 nrm = F.normal_at(b);
        Vec3 r1 = net->restore_from_frame(p + nrm * (sg * eps), fi, hw);
        Vec3 r2 = net->restore_from_frame(p + nrm * (sg * 2 * eps), fi, hw);
        double s1 = r1.z - st->z(h, r1.x, r1.y), s2 = r2.z - st->z(h, r2.x, r2.y);
        return 2 * s1 - s2;
      };
      auto cutoff = [&](std::size_t h, double a, bool hw, double& b_out) {
        double bu, bl;
        F.tips(a, bu, bl);
        if (F.abuts >= 0 && !F.curved && !net->faults()[static_cast<std::size_t>(F.abuts)].curved) {
          // only the part of the fault on its active side of the abutted fault
          const Fault& M = net->faults()[static_cast<std::size_t>(F.abuts)];
          double den = dot(M.n, F.dd);
          if (std::fabs(den) > 1e-12) {
            double bb = -M.sigma(F.at(a, 0)) / den;
            bool above_active = (M.sigma(F.at(a, bb - 1)) > 0) == F.abut_hw_side;
            if (above_active) bl = std::min(bl, bb - 0.5); else bu = std::max(bu, bb + 0.5);
          }
        }
        if (!(bl > bu)) return false;
        double b0 = bu, s0 = wall_s(h, a, b0, hw);
        for (double b1 = bu + step;; b1 += step) {
          b1 = std::min(b1, bl);
          double s1 = wall_s(h, a, b1, hw);
          if ((s0 >= 0) != (s1 >= 0)) {
            b_out = illinois([&](double b) { return wall_s(h, a, b, hw); }, b0, b1, s0, s1, 1e-6);
            return true;
          }
          if (b1 >= bl) return false;
          b0 = b1;
          s0 = s1;
        }
      };
      const int NA = 60;
      for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
        double worst = 0, tmax = 0;
        std::size_t n = 0;
        for (int k = 1; k < NA; ++k) {
          double a = -0.5 * F.length + F.length * k / NA;
          double bhw, bfw;
          if (!cutoff(h, a, true, bhw) || !cutoff(h, a, false, bfw)) continue;
          Vec3 chw = F.at(a, bhw), cfw = F.at(a, bfw);
          Vec3 cm = fixed_hw ? cfw : chw;
          double expect = 0;
          if (F.drag) {
            // drag moves both cut-offs along the slip (STX-TN-05 section 3): undo it first, then
            // add its vertical parts to the expected separation
            Vec3 p_hw = net->drag_inverse(fi, chw, true), p_fw = net->drag_inverse(fi, cfw, false);
            net->undo(fi, fixed_hw ? p_fw : p_hw, fixed_hw ? -1 : 1, &expect);
            expect += net->drag_forward(fi, p_fw, false).z - net->drag_forward(fi, p_hw, true).z;
          } else {
            net->undo(fi, cm, fixed_hw ? -1 : 1, &expect);
          }
          double meas;
          const bool along_dip = F.kin->name() == "dip_slip" && std::fabs(dot(F.slip, F.s)) < 1e-9;
          if (along_dip) {
            meas = cfw.z - chw.z;
          } else {
            // fixed wall's horizon on the vertical through the position the moving cut-off came
            // from (the cut-off itself for vertical offsets; slip with a strike component moves it)
            Vec3 from = net->undo(fi, cm, fixed_hw ? -1 : 1);
            bool hwfix = fixed_hw;
            auto sv = [&](double z) {
              Vec3 p{from.x, from.y, z};
              Vec3 r = net->restore_from_frame(p, fi, hwfix);
              return r.z - st->z(h, r.x, r.y);
            };
            double zf = expand_and_solve(sv, cm.z, 5.0, 1e4, 1e-7);
            meas = fixed_hw ? cm.z - zf : zf - cm.z;
            expect *= -F.slip.z;  // vertical part of the displacement (1 for dip-slip and vertical offsets)
          }
          worst = std::max(worst, std::fabs(meas - expect));
          tmax = std::max(tmax, std::fabs(expect));
          ++n;
        }
        if (n == 0) continue;
        auto r = base(c, "C4", "throw at the horizon equals the fault profile", inf.faults[f].name + " / " + inf.horizons[h].name);
        r.value = worst;
        r.tolerance = tol;
        r.tolerance_source = src;
        r.passed = worst <= tol;
        r.detail = std::to_string(n) + " strike positions (fault frame); profile maximum " + fmt_fixed(tmax, 2) + " m";
        out.push_back(r);
      }
      (void)bmin;
      (void)bmax;
    }
    return out;
  }
};

/// C10: displacement-length scaling (STX-TN-05 section 4.4). A warning: a stated throw is never
/// changed. Faults longer than the model diagonal are reported for information.
class C10 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    const auto& inf = c.model().info();
    const auto& d = c.description();
    std::vector<CheckResult> out;
    auto* net = static_cast<const FaultNetwork*>(c.model().capability("faults"));
    const double lo = d.scaling_min_ratio.present ? d.scaling_min_ratio.si() : 0.001;
    const double hi = d.scaling_max_ratio.present ? d.scaling_max_ratio.si() : 0.1;
    const std::string src = d.scaling_min_ratio.present || d.scaling_max_ratio.present
                                ? std::string("checks.fault_scaling")
                                : std::string("Stratonix default (Schultz et al. 2008; Lathrop et al. 2022)");
    const double diag = std::hypot(inf.domain.sx, inf.domain.sy);
    for (std::size_t f = 0; f < inf.faults.size() && f < d.faults.size(); ++f) {
      const auto& fd = d.faults[f];
      double L = fd.length.present ? fd.length.si() : 0.0;
      double D = net && f < net->size() ? std::fabs(net->faults()[f].tmax) : (fd.throw_max.present ? std::fabs(fd.throw_max.si()) : 0.0);
      if (!(L > 0)) continue;
      auto r = base(c, "C10", "displacement-length scaling", inf.faults[f].name);
      r.unit = "ratio";
      r.value = D / L;
      r.tolerance = hi;
      r.tolerance_source = src;
      r.warn_only = true;
      if (L > diag) {
        r.info = true;
        r.passed = true;
        r.detail = "length " + fmt_fixed(L, 0) + " m exceeds the model diagonal; reported for information";
      } else {
        r.passed = r.value >= lo && r.value <= hi;
        r.detail = "maximum displacement " + fmt_fixed(D, 1) + " m over length " + fmt_fixed(L, 0) + " m; range " +
                   fmt_fixed(lo, 4) + " to " + fmt_fixed(hi, 4);
      }
      out.push_back(r);
    }
    return out;
  }
};

/// C11: fault-system consistency (STX-TN-05 section 4.4). Splays with transfer: the splay's
/// displacement continues on the block beyond the branch line. Relays: the summed throw through
/// the overlap equals the equivalent single fault's (a linear blend of the two throws).
class C11 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    auto* net = static_cast<const FaultNetwork*>(c.model().capability("faults"));
    if (!net) return out;
    std::string src;
    const double tol = c.tolerance("throw", &src);
    const auto& fs = net->faults();
    const Domain& dom = c.model().info().domain;
    auto mid_b = [](const Fault& F, double a) {
      double bu, bl;
      F.tips(a, bu, bl);
      return 0.5 * (bu + bl);
    };
    for (const Fault& F : fs) {
      if (F.transfer && F.abuts >= 0 && !F.curved) {
        const Fault& M = fs[static_cast<std::size_t>(F.abuts)];
        Vec3 dir = cross(M.n, F.n);
        if (norm(dir) < 1e-9) continue;
        dir = normalized(dir);
        Vec3 away = normalized(M.dd - dir * dot(M.dd, dir));
        double worst = 0;
        int n = 0;
        // points of the branch line inside the domain: intersect the two planes along the splay's strike
        for (int k = 1; k < 40; ++k) {
          double a = -0.5 * F.length + F.length * k / 40.0;
          // the splay's down-dip line at strike position a meets the master
          Vec3 o = F.at(a, 0);
          double den = dot(M.n, F.dd);
          if (std::fabs(den) < 1e-12) continue;
          double b = -M.sigma(o) / den;
          Vec3 B = F.at(a, b);
          if (!dom.inside_xy(B.x, B.y) || B.z > dom.z_top || B.z < dom.z_base) continue;
          double am, bm;
          M.profile_coords(B, am, bm);
          if (std::fabs(am) >= 0.5 * M.length) continue;
          double ta, tb;
          F.profile_coords(B, ta, tb);
          double ts = F.throw_at(ta, tb);
          Vec3 beyond = B + away * 2.0 + M.n * 2.0;
          bool active = true;
          Vec3 r = net->relative_displacement(F.index, beyond, &active);
          double carried = active && F.sigma(beyond) > 0 ? -r.z : 0.0;
          worst = std::max(worst, std::fabs(ts - carried));
          ++n;
        }
        if (n == 0) continue;
        auto r = base(c, "C11", "fault-system consistency", F.name + " splay of " + M.name);
        r.value = worst;
        r.tolerance = tol;
        r.tolerance_source = src;
        r.passed = worst <= tol;
        r.detail = std::to_string(n) + " points of the branch line; splay throw against the throw carried beyond it";
        out.push_back(r);
      }
      if (F.relay_with > F.index) {
        const Fault& G = fs[static_cast<std::size_t>(F.relay_with)];
        // throws just outside the overlap on each fault's full side
        double span = F.ra_zero - F.ra_one, ext = 0.02 * std::fabs(span);
        double aF_full = F.ra_one - (span > 0 ? ext : -ext);
        double tF = F.throw_at(aF_full, mid_b(F, aF_full));
        double aG_full = G.ra_one - ((G.ra_zero - G.ra_one) > 0 ? ext : -ext);
        double tG = G.throw_at(aG_full, mid_b(G, aG_full));
        double worst = 0;
        for (int k = 0; k <= 40; ++k) {
          double w = k / 40.0;
          double aF = F.ra_one + w * span;
          Vec3 p = F.at(aF, mid_b(F, aF));
          double aG = dot(p - G.p0, G.s);
          double sum = F.throw_at(aF, mid_b(F, aF)) + G.throw_at(aG, mid_b(G, aG));
          worst = std::max(worst, std::fabs(sum - ((1 - w) * tF + w * tG)));
        }
        auto r = base(c, "C11", "fault-system consistency", F.name + " relays with " + G.name);
        r.value = worst;
        r.tolerance = tol;
        r.tolerance_source = src;
        r.passed = worst <= tol;
        r.detail = "summed throw through the overlap against the equivalent single fault (" + fmt_fixed(tF, 1) + " to " +
                   fmt_fixed(tG, 1) + " m)";
        out.push_back(r);
      }
    }
    return out;
  }
};

/// C12: area conserved by fold kinematics (STX-TN-06 section 6): Jacobian determinant of the
/// fault's own restoration in its dip plane at the reference strike position.
class C12 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    auto* net = static_cast<const FaultNetwork*>(c.model().capability("faults"));
    if (!net) return out;
    std::string src;
    const double tol = c.tolerance("area", &src);
    const double e = 0.5;
    for (const Fault& F : net->faults()) {
      if (!F.tri && !F.curved) continue;
      // dip-plane basis (along, across) and the region sampled
      Vec3 ua, uc;
      std::vector<std::pair<double, double>> pts;  // (along, across) offsets from the origin
      Vec3 o;
      double bu, bl;
      F.tips(0, bu, bl);
      if (F.tri) {
        double S = F.slip_at(0), ps = F.tri->ps, tp = F.tri->tanphi;
        double len = ps * S + std::fabs(S) + 600, wid = len * tp + 200;
        ua = F.dd;
        uc = F.n;
        o = F.at(0, bu);
        for (int i = 0; i <= 24; ++i)
          for (int j = 0; j <= 24; ++j) {
            double al = -len + (len + ps * S + 300) * i / 24.0, ac = -wid + 2 * wid * j / 24.0;
            if (al > 0 && std::fabs(ac) < 3 * e + 1) continue;  // the fault surface behind the present tip
            pts.push_back({al, ac});
          }
      } else {
        ua = F.dh;
        uc = Vec3{0, 0, 1};
        o = F.at(0, 0);
        for (int i = 0; i <= 24; ++i)
          for (int j = 1; j <= 24; ++j) {
            double u = -500 + 4000 * i / 24.0;
            double zf = -F.curved->zf(u);
            pts.push_back({u, zf - o.z + 60.0 * j});  // above the fault
          }
      }
      auto map = [&](const Vec3& q) {
        if (F.tri) return net->trishear_map(F.index, q, false);
        return F.sigma(q) > 0 ? net->undo(F.index, q, 1) : q;
      };
      double worst = 0;
      int n = 0;
      for (auto& [al, ac] : pts) {
        Vec3 p = o + ua * al + uc * ac;
        if (F.curved && F.sigma(p) < 3 * e + 1) continue;
        Vec3 ra = map(p + ua * e), rb = map(p - ua * e), rc = map(p + uc * e), rd = map(p - uc * e);
        double a11 = dot(ra - rb, ua) / (2 * e), a12 = dot(rc - rd, ua) / (2 * e);
        double a21 = dot(ra - rb, uc) / (2 * e), a22 = dot(rc - rd, uc) / (2 * e);
        worst = std::max(worst, std::fabs(a11 * a22 - a12 * a21 - 1));
        ++n;
      }
      auto r = base(c, "C12", "area conserved by fold kinematics", F.name);
      r.unit = "fraction";
      r.value = worst;
      r.tolerance = tol;
      r.tolerance_source = src;
      r.passed = worst <= tol;
      r.detail = std::to_string(n) + " points in the dip plane at the reference strike position; |Jacobian - 1|";
      out.push_back(r);
    }
    return out;
  }
};

/// C13: no negative thickness in the restored frame, after the series relations and, within each
/// series, before them (STX-TN-07 section 6). Designed pinch-outs (trends clamped) are exempt.
class C13 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    auto* st = static_cast<const Structure*>(c.model().capability("structure"));
    if (!st || !st->has_relations()) return out;
    std::string src;
    const double tol = c.tolerance("crossing", &src);
    const auto& inf = c.model().info();
    const std::size_t nh = st->horizons();
    const auto& g = c.grid();
    std::vector<double> raw(nh), eff(nh), worst(nh, 1e300);
    std::vector<unsigned char> fl(nh);
    std::vector<std::string> where(nh);
    for (std::size_t k = 0; k < g.size(); ++k) {
      const double x = g.x(k % g.nx), y = g.y(k / g.nx);
      st->z_detail(x, y, raw.data(), eff.data(), fl.data());
      for (std::size_t u = 1; u < nh; ++u) {
        double t = eff[u - 1] - eff[u];
        if (st->series_of(u - 1) == st->series_of(u) && !st->trend_clamped(u, x, y)) t = std::min(t, raw[u - 1] - raw[u]);
        if (t < worst[u]) { worst[u] = t; where[u] = "(" + fmt_fixed(x, 0) + ", " + fmt_fixed(y, 0) + ")"; }
      }
    }
    for (std::size_t u = 1; u < nh; ++u) {
      auto r = base(c, "C13", "no negative thickness (restored frame, series relations)", inf.units[u].name);
      r.value = worst[u];
      r.tolerance = -tol;
      r.tolerance_source = src;
      r.passed = worst[u] >= -tol;
      r.detail = "minimum thickness at " + where[u] + ", before the relations within a series and after them";
      out.push_back(r);
    }
    return out;
  }
};

/// C14: truncation consistent with age order (STX-TN-07 section 6): a lowered horizon lies on a
/// younger erosional surface; a raised horizon lies on an older surface.
class C14 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    auto* st = static_cast<const Structure*>(c.model().capability("structure"));
    if (!st || !st->has_relations()) return out;
    const auto& inf = c.model().info();
    const std::size_t nh = st->horizons();
    const auto& g = c.grid();
    std::vector<double> raw(nh), eff(nh);
    std::vector<unsigned char> fl(nh);
    std::vector<std::size_t> low(nh, 0), up(nh, 0), bad(nh, 0);
    for (std::size_t k = 0; k < g.size(); ++k) {
      st->z_detail(g.x(k % g.nx), g.y(k / g.nx), raw.data(), eff.data(), fl.data());
      for (std::size_t h = 0; h < nh; ++h) {
        if (fl[h] & 1) {
          ++low[h];
          bool ok = false;  // a younger horizon at the same level
          for (std::size_t j = 0; j < h && !ok; ++j) ok = st->series_of(j) < st->series_of(h) && std::fabs(eff[j] - eff[h]) < 1e-9;
          if (!ok) ++bad[h];
        }
        if (fl[h] & 2) {
          ++up[h];
          bool ok = false;  // an older horizon at the same level
          for (std::size_t j = h + 1; j < nh && !ok; ++j) ok = st->series_of(j) > st->series_of(h) && std::fabs(eff[j] - eff[h]) < 1e-9;
          if (!ok) ++bad[h];
        }
      }
    }
    for (std::size_t h = 0; h < nh; ++h) {
      if (!low[h] && !up[h]) continue;
      auto r = base(c, "C14", "truncation consistent with age order", inf.horizons[h].name);
      r.unit = "columns";
      r.value = static_cast<double>(bad[h]);
      r.tolerance = 0;
      r.tolerance_source = "exact";
      r.passed = bad[h] == 0;
      r.detail = std::to_string(low[h]) + " truncated and " + std::to_string(up[h]) + " onlapping columns of " +
                 std::to_string(g.size());
      out.push_back(r);
    }
    return out;
  }
};

/// C15: salt volume conserved (STX-TN-08 section 5).
class C15 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    auto* st = static_cast<const Structure*>(c.model().capability("structure"));
    if (!st || (!st->salt_parametric() && st->salt_bodies_volume() == 0 && st->salt_net_volume() == 0)) return out;
    std::string src;
    const double tol = c.tolerance("salt_volume", &src);
    auto r = base(c, "C15", "salt volume conserved", "salt");
    r.unit = "fraction";
    if (st->salt_parametric()) {
      r.value = std::fabs(st->salt_net_volume()) / std::max(1.0, st->salt_bodies_volume());
      r.tolerance = tol;
      r.tolerance_source = src;
      r.passed = r.value <= tol;
      r.detail = "net change " + fmt_fixed(st->salt_net_volume(), 0) + " m3 against the bodies' " + fmt_fixed(st->salt_bodies_volume(), 0) + " m3";
    } else {
      r.info = true;
      r.value = st->salt_net_volume();
      r.unit = "m3";
      r.detail = "imported top salt: volume relative to the source layer, reported only";
    }
    out.push_back(r);
    return out;
  }
};

/// C16: solid (grain) thickness conserved through compaction and decompaction (STX-TN-08 section 5).
class C16 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    auto* st = static_cast<const Structure*>(c.model().capability("structure"));
    if (!st || st->compaction_mode() == 0) return out;
    std::string src;
    const double tol = c.tolerance("solid_thickness", &src);
    const auto& inf = c.model().info();
    const auto& g = c.grid();
    const std::size_t step = std::max<std::size_t>(1, g.size() / 2500);  // about 2,500 columns
    for (std::size_t u = 1; u + 1 < inf.units.size(); ++u) {
      if (!st->compaction_law(static_cast<int>(u))) continue;
      double worst = 0;
      std::size_t n = 0;
      for (std::size_t k = 0; k < g.size(); k += step) {
        const double x = g.x(k % g.nx), y = g.y(k / g.nx);
        const double s = st->solid_thickness(static_cast<int>(u), x, y);
        if (!(s > 1e-6)) continue;
        if (st->compaction_mode() == 2) worst = std::max(worst, std::fabs(s - st->depositional_solid(static_cast<int>(u), x, y)) / s);
        std::vector<double> col(st->horizons());
        st->z_all(x, y, col.data());
        const double t = col[u - 1] - col[u];
        worst = std::max(worst, std::fabs(st->compaction_round_trip(static_cast<int>(u), x, y) - t) / std::max(t, 1e-9));
        ++n;
      }
      auto r = base(c, "C16", "solid thickness conserved through compaction", inf.units[u].name);
      r.unit = "fraction";
      r.value = worst;
      r.tolerance = tol;
      r.tolerance_source = src;
      r.passed = worst <= tol;
      r.detail = std::to_string(n) + " columns; " + (st->compaction_mode() == 2 ? "present against depositional solid thickness, and " : "") +
                 "decompaction to the surface and back";
      out.push_back(r);
    }
    return out;
  }
};

/// C17: restorability at every event (STX-TN-09 section 4).
class C17 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const auto& d = c.description();
    if (d.events.empty() || !c.timeline_checks || !c.model().capability("structure")) return out;
    std::string src;
    std::string src2;
    const double tol = c.tolerance("restoration", &src), cross = c.tolerance("crossing", &src2);
    auto ages = event_ages(d);
    ages.push_back(0.0);
    const auto& g = c.grid();
    const std::size_t step = std::max<std::size_t>(1, g.size() / 1500);
    double round_trip = 0, worst_thk = 1e300;
    std::string bad;
    std::size_t built = 0;
    for (double a : ages) {
      ModelPtr s;
      try {
        s = build_state(d, c.model(), a);
      } catch (const std::exception& e) {
        bad += (bad.empty() ? "" : "; ") + std::string("state at ") + fmt_shortest(a) + " Ma does not build: " + e.what();
        continue;
      }
      ++built;
      auto* st = static_cast<const Structure*>(s->capability("structure"));
      const std::size_t nh = st->horizons();
      std::vector<double> z(nh);
      for (std::size_t k = 0; k < g.size(); k += step) {
        st->z_all(g.x(k % g.nx), g.y(k / g.nx), z.data());
        for (std::size_t h = 1; h < nh; ++h) worst_thk = std::min(worst_thk, z[h - 1] - z[h]);
      }
      if (a == 0.0) {  // the present state reproduces the full model
        auto* full = static_cast<const Structure*>(c.model().capability("structure"));
        std::vector<double> zf(full->horizons());
        for (std::size_t k = 0; k < g.size(); k += step) {
          const double x = g.x(k % g.nx), y = g.y(k / g.nx);
          st->z_all(x, y, z.data());
          full->z_all(x, y, zf.data());
          for (std::size_t h = 0; h < nh && h < zf.size(); ++h) round_trip = std::max(round_trip, std::fabs(z[h] - zf[h]));
        }
      }
    }
    auto r = base(c, "C17", "restorability at every event", "timeline");
    r.value = round_trip;
    r.tolerance = tol;
    r.tolerance_source = src;
    r.passed = bad.empty() && round_trip <= tol && worst_thk >= -cross;
    r.detail = std::to_string(built) + " of " + std::to_string(ages.size()) + " states built; minimum thickness " +
               fmt_fixed(worst_thk, 3) + " m; present state against the model" + (bad.empty() ? "" : "; " + bad);
    out.push_back(r);
    return out;
  }
};

/// C18: line-length and area balance per section (STX-TN-09 section 4), a warning.
class C18 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const auto& d = c.description();
    if (d.events.empty() || !c.timeline_checks) return out;
    std::string src;
    const double tol = c.tolerance("balance", &src);
    for (auto& b : balance_rows(c.model(), d)) {
      auto r = base(c, "C18", "balance (" + b.method + ")", b.section + " / " + b.subject);
      r.unit = "fraction";
      r.value = b.mismatch;
      r.tolerance = tol;
      r.tolerance_source = src;
      r.passed = !b.tested || b.mismatch <= tol;
      r.info = !b.tested;
      r.warn_only = true;
      r.detail = "present " + fmt_fixed(b.present, 1) + ", restored " + fmt_fixed(b.restored, 1);
      out.push_back(r);
    }
    return out;
  }
};

/// C5: control wells honoured (check wells reported as warnings).
class C5 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    std::string src;
    double tol = c.tolerance("wells", &src);
    const auto& d = c.description();
    for (std::size_t w = 0; w < d.wells.size(); ++w) {
      const auto& wd = d.wells[w];
      if (wd.tops.empty()) continue;
      auto tops = repr::well_intersections(c.model(), c.wells()[w]);
      for (auto& [hn, p] : wd.tops) {
        double want = -p.si(), best = std::numeric_limits<double>::infinity();
        for (auto& t : tops)
          if (t.horizon == hn) best = std::min(best, std::fabs(t.p.z - want));
        auto r = base(c, "C5", wd.role == "control" ? "control well honoured" : "check well compared", wd.name + " / " + hn);
        r.value = std::isfinite(best) ? best : 1e9;
        r.tolerance = tol;
        r.tolerance_source = src;
        r.passed = r.value <= tol;
        r.warn_only = r.warn_only || wd.role != "control";
        r.detail = std::isfinite(best) ? "top at " + fmt_fixed(-want, 2) + " m depth" : "horizon not intersected by the well";
        out.push_back(r);
      }
    }
    return out;
  }
};

/// C6: top and base horizons unoffset unless intended.
class C6 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    const auto& inf = c.model().info();
    std::vector<CheckResult> out;
    if (inf.faults.empty()) return out;
    std::string src;
    double tol = c.tolerance("unoffset", &src);
    for (int h : {0, static_cast<int>(inf.horizons.size()) - 1}) {
      if (inf.horizons[static_cast<std::size_t>(h)].intended_offset) continue;
      double worst = 0;
      std::string fn;
      for (std::size_t f = 0; f < inf.faults.size(); ++f) {
        const auto& t = c.throw_row(static_cast<int>(f), h);
        double v = std::max(std::fabs(t.throw_max), std::fabs(t.jump_max));
        if (v > worst) { worst = v; fn = inf.faults[f].name; }
      }
      auto r = base(c, "C6", h == 0 ? "top horizon unoffset" : "base horizon unoffset", inf.horizons[static_cast<std::size_t>(h)].name);
      r.value = worst;
      r.tolerance = tol;
      r.tolerance_source = src;
      r.passed = worst <= tol;
      r.detail = fn.empty() || worst <= tol ? "not offset by any fault"
                                            : "offset by " + fn + " (set intended_offset: true if this is intended)";
      out.push_back(r);
    }
    return out;
  }
};

/// C7: faults inside the domain.
class C7 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    const auto& inf = c.model().info();
    std::vector<CheckResult> out;
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      const auto& fs = c.fault_lattice(static_cast<int>(f));
      double area = fs.mesh.area();
      auto r = base(c, "C7", "fault inside the domain", inf.faults[f].name);
      r.value = area / 1e6;
      r.unit = "km2";
      r.tolerance = 0;
      r.tolerance_source = "must be > 0";
      r.passed = area > 0;
      auto q = mesh_quality(fs.mesh);
      r.detail = "fault surface area inside the domain and tip line; " + std::to_string(q.triangles) + " triangles, " +
                 std::to_string(q.degenerate) + " degenerate";
      if (q.degenerate) r.passed = false;
      out.push_back(r);
    }
    return out;
  }
};

/// C8: fault-fault intersections resolved (ordered by age or abutting; abutting faults inactive
/// on the far side).
class C8 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    const auto& inf = c.model().info();
    std::vector<CheckResult> out;
    const auto& d = c.description();
    for (std::size_t f = 0; f < inf.faults.size(); ++f)
      for (std::size_t g = f + 1; g < inf.faults.size(); ++g) {
        const auto& ff = c.fault_lattice(static_cast<int>(f));
        std::vector<Vec3> inside;
        for (std::size_t k = 0; k < ff.nodes.size(); ++k)
          if (ff.extent[k] <= 1) inside.push_back(ff.nodes[k]);
        FaultEval eg;
        c.model().evaluate_fault(static_cast<int>(g), inside, eg);
        bool pos = false, neg = false;
        for (std::size_t k = 0; k < inside.size(); ++k)
          if (eg.extent[k] <= 1) {
            (eg.distance[k] > 0 ? pos : neg) = true;
          }
        bool intersect = pos && neg;
        if (!intersect) continue;
        const auto& fa = d.faults[f];
        const auto& ga = d.faults[g];
        bool ordered = fa.age.si() != ga.age.si();
        bool abut = ga.abuts_fault == fa.name || fa.abuts_fault == ga.name;
        double leak = 0;
        if (abut) {
          // the younger, abutting fault must not displace anything on the far side
          std::size_t y = ga.abuts_fault == fa.name ? g : f, o = y == g ? f : g;
          const auto& fy = c.fault_lattice(static_cast<int>(y));
          FaultEval eo, ey;
          c.model().evaluate_fault(static_cast<int>(o), fy.nodes, eo);
          c.model().evaluate_fault(static_cast<int>(y), fy.nodes, ey);
          bool hw_side = d.faults[y].abuts_side == "hanging_wall";
          for (std::size_t k = 0; k < fy.nodes.size(); ++k)
            if ((eo.distance[k] > 0) != hw_side && std::fabs(eo.distance[k]) > 1) leak = std::max(leak, ey.displacement[k]);
        }
        auto r = base(c, "C8", "fault-fault intersection resolved", inf.faults[f].name + " x " + inf.faults[g].name);
        r.value = leak;
        r.tolerance = 0.001;
        r.tolerance_source = "abutting fault must not move the far side";
        r.passed = (ordered || abut) && leak <= 0.001;
        r.detail = std::string("faults intersect inside their extents; ") + (ordered ? "ordered by age" : "same age") +
                   (abut ? ", abutting declared" : "");
        out.push_back(r);
      }
    return out;
  }
};

/// C9: closure, spill point and closure area, compared with targets when stated.
class C9 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const auto& d = c.description();
    if (d.closure_horizon.empty()) return out;
    int h = d.horizon_index(d.closure_horizon);
    const auto& fs = c.faulted_surface(h);
    auto cl = repr::closure_analysis(fs.z, d.closure_faults_seal.flag);
    auto info = [&](const std::string& name, double v, const std::string& unit, const std::string& detail) {
      auto r = base(c, "C9", name, d.closure_horizon);
      r.value = v;
      r.unit = unit;
      r.info = true;
      r.detail = detail;
      out.push_back(r);
    };
    info("crest depth", -cl.crest_z, "m", "at (" + fmt_fixed(cl.crest.x, 0) + ", " + fmt_fixed(cl.crest.y, 0) + ")");
    info("spill depth", -cl.spill_z, "m", "spill point (" + fmt_fixed(cl.spill_point.x, 0) + ", " + fmt_fixed(cl.spill_point.y, 0) + ")" +
                                         (cl.spill_at_boundary ? " at the domain edge" : ""));
    info("closure", cl.closure, "m", d.closure_faults_seal.flag ? "faults sealing" : "faults not sealing");
    info("closure area", cl.area / 1e6, "km2", "check grid spacing " + fmt_fixed(c.grid().dx, 1) + " m");
    auto target = [&](const desc::Param& t, double value, const std::string& name, double scale, const std::string& unit) {
      if (!t.present) return;
      auto r = base(c, "C9", name + " against target", d.closure_horizon);
      double want = t.si() * scale;
      r.value = value - want;
      r.unit = unit;
      if (d.target_tolerance.present) {
        r.tolerance = std::fabs(want) * d.target_tolerance.si();
        r.tolerance_source = "checks.closure.targets.tolerance_fraction";
        r.passed = std::fabs(r.value) <= r.tolerance;
      } else {
        r.info = true;
      }
      r.detail = "model " + fmt_fixed(value, 2) + " vs target " + fmt_fixed(want, 2) + " (difference reported; no parameter is tuned)";
      out.push_back(r);
    };
    target(d.target_closure, cl.closure, "closure", 1.0, "m");
    target(d.target_spill_depth, -cl.spill_z, "spill depth", 1.0, "m");
    // area target: unit m2 or km2 (SI value in m2), reported in km2
    target(d.target_area, cl.area / 1e6, "closure area", 1e-6, "km2");
    return out;
  }
};

}  // namespace

STRATONIX_REGISTER(check_c1, registry::Kind::Check, "C1_thickness", "1.0.0", "Positive thickness in every fault block",
                   [] { return std::make_unique<C1>(); })
STRATONIX_REGISTER(check_c2, registry::Kind::Check, "C2_crossing", "1.0.0", "No horizon crossing",
                   [] { return std::make_unique<C2>(); })
STRATONIX_REGISTER(check_c3, registry::Kind::Check, "C3_offset_outside_tips", "1.0.0", "No offset outside a fault's tips",
                   [] { return std::make_unique<C3>(); })
STRATONIX_REGISTER(check_c4, registry::Kind::Check, "C4_throw_profile", "1.0.0",
                   "Throw at each horizon equals the fault profile (parametric engine)", [] { return std::make_unique<C4>(); })
STRATONIX_REGISTER(check_c5, registry::Kind::Check, "C5_wells", "1.0.0", "Control wells honoured",
                   [] { return std::make_unique<C5>(); })
STRATONIX_REGISTER(check_c6, registry::Kind::Check, "C6_top_base_unoffset", "1.0.0", "Top and base unoffset unless intended",
                   [] { return std::make_unique<C6>(); })
STRATONIX_REGISTER(check_c7, registry::Kind::Check, "C7_faults_in_domain", "1.0.0", "Faults inside the domain; mesh quality",
                   [] { return std::make_unique<C7>(); })
STRATONIX_REGISTER(check_c8, registry::Kind::Check, "C8_fault_intersections", "1.0.0", "Fault-fault intersections resolved",
                   [] { return std::make_unique<C8>(); })
/// Fault attributes on each fault surface at the check spacing (stage 2F).
std::vector<std::vector<repr::FaultRockPoint>> fault_rock_nodes(CheckContext& c, double max_extent) {
  std::vector<std::vector<repr::FaultRockPoint>> out;
  const auto& d = c.description();
  const auto in = fault_rock_input(d, c.model());
  const auto st = stress_input(d);
  const double spacing = std::max(c.model().info().domain.spacing * 2, 50.0);
  for (std::size_t f = 0; f < c.model().info().faults.size(); ++f) {
    auto fs = repr::fault_surface(c.model(), static_cast<int>(f), spacing);
    out.push_back(repr::fault_rock_surface(c.model(), fs, in, st, max_extent));
  }
  return out;
}

/// C19: juxtaposition consistent with the cut-offs (STX-TN-10 section 6).
class C19 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const auto& d = c.description();
    if (!d.fault_rock.present || !c.model().capability("column") || c.model().info().faults.empty()) return out;
    auto nodes = fault_rock_nodes(c, 0.98);
    for (std::size_t f = 0; f < nodes.size(); ++f) {
      std::size_t n = 0, bad = 0;
      for (auto& p : nodes[f]) {
        if (!p.valid || p.near_horizon || p.fw_unit < 0) continue;
        ++n;
        if (p.fw_unit != p.fw_unit_col || p.hw_unit != p.hw_unit_col) ++bad;
      }
      auto r = base(c, "C19", "juxtaposition consistent with the cut-offs", c.model().info().faults[f].name);
      r.unit = "fraction";
      r.value = n ? static_cast<double>(bad) / static_cast<double>(n) : 0.0;
      r.tolerance = 0;
      r.tolerance_source = "STX-TN-10 section 6 (exact)";
      r.passed = bad == 0;
      r.detail = std::to_string(n) + " nodes; " + std::to_string(bad) + " disagree";
      out.push_back(r);
    }
    return out;
  }
};

/// C20: slip and dilation tendency within their bounds (STX-TN-10 section 6).
class C20 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const auto& d = c.description();
    if (!d.stress.present || c.model().info().faults.empty()) return out;
    auto nodes = fault_rock_nodes(c, 1.0);
    for (std::size_t f = 0; f < nodes.size(); ++f) {
      std::size_t n = 0, bad = 0;
      double worst = 0;
      for (auto& p : nodes[f]) {
        if (!p.valid || std::isnan(p.td)) continue;
        ++n;
        const double e = 1e-9;
        const bool ok = p.td >= -e && p.td <= 1 + e && (std::isnan(p.ts) || (p.ts >= -e && p.ts <= p.ts_max * (1 + 1e-9) + e));
        if (!ok) ++bad;
        if (!std::isnan(p.ts) && p.ts_max > 0) worst = std::max(worst, p.ts / p.ts_max);
      }
      auto r = base(c, "C20", "slip and dilation tendency bounded", c.model().info().faults[f].name);
      r.unit = "count";
      r.value = static_cast<double>(bad);
      r.tolerance = 0;
      r.tolerance_source = "STX-TN-10 section 6 (by definition)";
      r.passed = bad == 0;
      r.detail = std::to_string(n) + " nodes; largest Ts / Ts_max " + fmt_fixed(worst, 3);
      out.push_back(r);
    }
    return out;
  }
};

STRATONIX_REGISTER(check_c9, registry::Kind::Check, "C9_closure", "1.0.0", "Closure, spill point and area against targets",
                   [] { return std::make_unique<C9>(); })
STRATONIX_REGISTER(check_c10, registry::Kind::Check, "C10_fault_scaling", "1.0.0", "Displacement-length scaling (warning)",
                   [] { return std::make_unique<C10>(); })
STRATONIX_REGISTER(check_c12, registry::Kind::Check, "C12_area", "1.0.0", "Area conserved by fold kinematics",
                   [] { return std::make_unique<C12>(); })
STRATONIX_REGISTER(check_c13, registry::Kind::Check, "C13_thickness_restored", "1.0.0", "No negative thickness (series relations)",
                   [] { return std::make_unique<C13>(); })
STRATONIX_REGISTER(check_c14, registry::Kind::Check, "C14_truncation_order", "1.0.0", "Truncation consistent with age order",
                   [] { return std::make_unique<C14>(); })
STRATONIX_REGISTER(check_c15, registry::Kind::Check, "C15_salt_volume", "1.0.0", "Salt volume conserved",
                   [] { return std::make_unique<C15>(); })
STRATONIX_REGISTER(check_c16, registry::Kind::Check, "C16_solid_thickness", "1.0.0", "Solid thickness conserved through compaction",
                   [] { return std::make_unique<C16>(); })
STRATONIX_REGISTER(check_c17, registry::Kind::Check, "C17_restorability", "1.0.0", "Restorability at every event",
                   [] { return std::make_unique<C17>(); })
STRATONIX_REGISTER(check_c18, registry::Kind::Check, "C18_balance", "1.0.0", "Line-length and area balance (warning)",
                   [] { return std::make_unique<C18>(); })
STRATONIX_REGISTER(check_c19, registry::Kind::Check, "C19_juxtaposition", "1.0.0", "Juxtaposition consistent with the cut-offs",
                   [] { return std::make_unique<C19>(); })
STRATONIX_REGISTER(check_c20, registry::Kind::Check, "C20_tendencies", "1.0.0", "Slip and dilation tendency bounded",
                   [] { return std::make_unique<C20>(); })
STRATONIX_REGISTER(check_c11, registry::Kind::Check, "C11_fault_system", "1.0.0",
                   "Throw continuous across branch lines and relays", [] { return std::make_unique<C11>(); })

std::vector<BalanceRow> balance_rows(const Model& m, const desc::Description& d) {
  std::vector<BalanceRow> rows;
  std::vector<std::pair<std::string, Polyline>> secs;
  for (auto& s : d.sections) {
    Polyline pl;
    for (std::size_t r = 0; r < s.polyline.rows(); ++r) pl.push_back({s.polyline.si(2 * r), s.polyline.si(2 * r + 1), 0});
    secs.push_back({s.name, pl});
  }
  const Domain& dm = m.info().domain;
  if (secs.empty()) secs.push_back({"X_centre", {{dm.x0, dm.y0 + dm.sy / 2, 0}, {dm.x0 + dm.sx, dm.y0 + dm.sy / 2, 0}}});
  const auto& inf = m.info();
  for (auto& [name, line] : secs) {
    auto s = repr::make_section(m, name, line, std::max(10.0, d.grid_spacing.si() / 2), 5.0);
    std::vector<std::vector<std::vector<Vec3>>> segs(inf.horizons.size());
    // a horizon where it coincides with a younger one is absent there (eroded, or lapping out)
    auto coincident = [&](std::size_t h, const Vec3& p) {
      for (std::size_t j = 0; j < h; ++j) {
        double zj = kNaN;
        double best = 1e300;
        for (std::size_t i = 0; i < s.ns; ++i)
          if (std::fabs(s.xy[i].x - p.x) + std::fabs(s.xy[i].y - p.y) < best) { best = std::fabs(s.xy[i].x - p.x) + std::fabs(s.xy[i].y - p.y); zj = s.horizon_z[j][i]; }
        if (!std::isnan(zj) && std::fabs(zj - p.z) < 0.05) return true;
      }
      return false;
    };
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      segs[h] = repr::horizon_segments(m, s, static_cast<int>(h));
      double lp = 0, lr = 0;
      for (auto& run : segs[h]) {
        std::vector<Vec3> rr(run.size());
        m.restore(run, rr);
        for (std::size_t k = 1; k < run.size(); ++k) {
          if (coincident(h, run[k]) || coincident(h, run[k - 1])) continue;
          lp += norm(run[k] - run[k - 1]);
          lr += norm(rr[k] - rr[k - 1]);
        }
      }
      if (lp > 0) rows.push_back({name, inf.horizons[h].name, "line length, faults restored by their own inverses", lp, lr, std::fabs(lr - lp) / lp});
    }
    // area of each unit: quads between consecutive columns where top and base are in one block
    for (std::size_t u = 1; u < inf.horizons.size(); ++u) {
      const auto &zt = s.horizon_z[u - 1], &zb = s.horizon_z[u];
      std::vector<Vec3> pts;
      std::vector<std::size_t> cols;
      for (std::size_t i = 0; i < s.ns; ++i)
        if (!std::isnan(zt[i]) && !std::isnan(zb[i]) && zt[i] - zb[i] > 0.01) {
          pts.push_back({s.xy[i].x, s.xy[i].y, zt[i] - 0.01});
          pts.push_back({s.xy[i].x, s.xy[i].y, zb[i] + 0.01});
          cols.push_back(i);
        }
      if (cols.size() < 2) continue;
      EvalResult er;
      m.evaluate(pts, kEvalBlock, er);
      std::vector<Vec3> rr(pts.size());
      m.restore(pts, rr);
      auto quad = [](const Vec3& a, const Vec3& b, const Vec3& cc, const Vec3& dd) {
        return 0.5 * norm(cross(cc - a, dd - b));
      };
      double ap = 0, ar = 0;
      for (std::size_t k = 0; k + 1 < cols.size(); ++k) {
        if (cols[k + 1] != cols[k] + 1) continue;
        const std::size_t a0 = 2 * k, b0 = 2 * k + 1, a1 = 2 * k + 2, b1 = 2 * k + 3;
        if (er.block[a0] != er.block[b0] || er.block[a0] != er.block[a1] || er.block[a1] != er.block[b1]) continue;
        ap += quad(pts[a0], pts[a1], pts[b1], pts[b0]);
        ar += quad(rr[a0], rr[a1], rr[b1], rr[b0]);
      }
      // only units deposited before every deformation event keep their area in the restored frame;
      // syn- and post-kinematic units are thicker or thinner by deposition, and are reported only
      bool pre = true;
      if (!d.events.empty()) {
        const auto& top = d.horizons[u - 1];
        for (auto& e : d.events)
          if (e.type != "deposition" && top.age.present && top.age.si() < e.start()) pre = false;
      }
      if (ap > 0)
        rows.push_back({name, inf.units[u].name,
                        pre ? "area, faults restored by their own inverses" : "area, syn- or post-kinematic unit (reported only)", ap, ar,
                        std::fabs(ar - ap) / ap, pre});
    }
  }
  return rows;
}

repr::FaultRockInput fault_rock_input(const desc::Description& d, const Model& m) {
  repr::FaultRockInput in;
  const auto& units = m.info().units;
  in.unit_clay.assign(units.size(), kNaN);
  for (const auto& u : d.units) {
    if (!u.clay_fraction.present) continue;
    const int b = d.horizon_index(u.base), t = d.horizon_index(u.top);
    for (std::size_t k = 0; k < units.size(); ++k)
      if ((b >= 0 && units[k].base_horizon == b) || (b < 0 && units[k].name == u.name) ||
          (b < 0 && t >= 0 && units[k].top_horizon == t && units[k].base_horizon < 0))
        in.unit_clay[k] = u.clay_fraction.si();
  }
  const auto& fr = d.fault_rock;
  if (fr.shale_threshold.present) in.shale_threshold = fr.shale_threshold.si();
  if (fr.thickness_ratio.present) in.thickness_ratio = fr.thickness_ratio.si();
  if (fr.seal_type == "sgr_power" && fr.seal_a.present && fr.seal_c.present) {
    in.seal = true;
    in.seal_a = fr.seal_a.si();
    in.seal_c = fr.seal_c.si();
  }
  return in;
}
repr::StressInput stress_input(const desc::Description& d) {
  repr::StressInput s;
  const auto& st = d.stress;
  if (!st.present) return s;
  if (st.tensor_gradient.present && st.tensor_gradient.values.size() == 6) {
    for (std::size_t k = 0; k < 6; ++k) s.g[k] = st.tensor_gradient.si(k);
    s.present = true;
  } else if (st.regime() && st.sv_gradient.present && st.shmax_ratio.present && st.shmin_ratio.present && st.shmax_azimuth.present) {
    const double sv = st.sv_gradient.si(), sh = st.shmax_ratio.si() * sv, sl = st.shmin_ratio.si() * sv, th = st.shmax_azimuth.si();
    const double sn = std::sin(th), cs = std::cos(th);  // SHmax along (sin, cos), Shmin along (cos, -sin)
    s.g[0] = sh * sn * sn + sl * cs * cs;
    s.g[1] = sh * cs * cs + sl * sn * sn;
    s.g[2] = sv;
    s.g[3] = (sh - sl) * sn * cs;
    s.present = true;
  }
  if (st.pore_pressure_gradient.present) s.pp = st.pore_pressure_gradient.si();
  if (st.friction.present) s.friction = st.friction.si();
  if (st.cohesion.present) s.cohesion = st.cohesion.si();
  return s;
}

CheckReport run_checks(CheckContext& ctx) {
  CheckReport rep;
  for (auto* e : registry::Registry::global().list(registry::Kind::Check)) {
    auto chk = registry::Registry::global().create<Check>(registry::Kind::Check, e->name);
    auto rs = chk->run(ctx);
    rep.results.insert(rep.results.end(), rs.begin(), rs.end());
  }
  // Report in numeric order of the check ids (C1, C2, ..., C10, C11), stable within a check.
  auto num = [](const std::string& id) { return id.size() > 1 && id[0] == 'C' ? std::atoi(id.c_str() + 1) : 1 << 20; };
  std::stable_sort(rep.results.begin(), rep.results.end(),
                   [&](const CheckResult& a, const CheckResult& b) { return num(a.id) < num(b.id); });
  return rep;
}

void register_check_plugins() {
  STRATONIX_PLUGIN_ANCHOR(check_c1);
  STRATONIX_PLUGIN_ANCHOR(check_c2);
  STRATONIX_PLUGIN_ANCHOR(check_c3);
  STRATONIX_PLUGIN_ANCHOR(check_c4);
  STRATONIX_PLUGIN_ANCHOR(check_c5);
  STRATONIX_PLUGIN_ANCHOR(check_c6);
  STRATONIX_PLUGIN_ANCHOR(check_c7);
  STRATONIX_PLUGIN_ANCHOR(check_c8);
  STRATONIX_PLUGIN_ANCHOR(check_c9);
  STRATONIX_PLUGIN_ANCHOR(check_c10);
  STRATONIX_PLUGIN_ANCHOR(check_c11);
  STRATONIX_PLUGIN_ANCHOR(check_c12);
  STRATONIX_PLUGIN_ANCHOR(check_c13);
  STRATONIX_PLUGIN_ANCHOR(check_c14);
  STRATONIX_PLUGIN_ANCHOR(check_c15);
  STRATONIX_PLUGIN_ANCHOR(check_c16);
  STRATONIX_PLUGIN_ANCHOR(check_c17);
  STRATONIX_PLUGIN_ANCHOR(check_c18);
  STRATONIX_PLUGIN_ANCHOR(check_c19);
  STRATONIX_PLUGIN_ANCHOR(check_c20);
}

}  // namespace stratonix
