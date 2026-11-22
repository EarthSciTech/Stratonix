// Checks that fail a build. Each check states its tolerance and where the tolerance came from.
#include "stratonix/checks/checks.hpp"

#include <chrono>
#include <algorithm>
#include <map>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/core/roots.hpp"
#include "stratonix/engines/faults.hpp"
#include "stratonix/engines/structure.hpp"
#include "stratonix/engines/bodies.hpp"
#include "stratonix/engines/depositional.hpp"
#include "stratonix/engines/restoration3d.hpp"
#include "stratonix/engines/conditioning.hpp"
#include "stratonix/engines/posterior.hpp"
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
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c21);
STRATONIX_DECLARE_PLUGIN_ANCHOR(check_c22);

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
    : m_(m), d_(d), wells_(std::move(wells)), surfaces_(std::make_shared<repr::FaultedSurfaceCache>(m)) {
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
  if (it != block_surf_.end()) return it->second;
  // The block's surface over its sampled map extent widened by a further lattice cell (two in
  // all); NaN beyond, where the block is not present. A network of many faults has many small
  // blocks, so a full grid per block was most of the cost of C1 and C2 (stage 4B).
  const auto& bl = blocks();
  const auto pos = static_cast<std::size_t>(std::lower_bound(bl.codes.begin(), bl.codes.end(), b) - bl.codes.begin());
  Grid2D out(grid_);
  if (pos < bl.codes.size() && bl.codes[pos] == b && !bl.extent.empty()) {
    auto a = bl.extent[pos];
    a[0] -= bl.cell; a[1] -= bl.cell; a[2] += bl.cell; a[3] += bl.cell;
    auto idx = [](double v, double v0, double dv, std::size_t n) {
      return static_cast<std::size_t>(std::clamp(std::floor((v - v0) / dv), 0.0, static_cast<double>(n - 1)));
    };
    const std::size_t i0 = idx(a[0], grid_.x0, grid_.dx, grid_.nx), j0 = idx(a[1], grid_.y0, grid_.dy, grid_.ny);
    const std::size_t i1 = std::min(grid_.nx - 1, idx(a[2], grid_.x0, grid_.dx, grid_.nx) + 1), j1 = std::min(grid_.ny - 1, idx(a[3], grid_.y0, grid_.dy, grid_.ny) + 1);
    GridSpec w{grid_.x(i0), grid_.y(j0), grid_.dx, grid_.dy, i1 - i0 + 1, j1 - j0 + 1};
    const Grid2D s = repr::horizon_grid_block(m_, h, b, w);
    for (std::size_t j = 0; j < w.ny; ++j)
      for (std::size_t i = 0; i < w.nx; ++i) out.at(i0 + i, j0 + j) = s.at(i, j);
  } else {
    out = repr::horizon_grid_block(m_, h, b, grid_);
  }
  return block_surf_.emplace(k, std::move(out)).first->second;
}
const repr::FaultedSurface& CheckContext::faulted_surface(int h) {
  return surfaces_->get(h, grid_);
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
const GridProducts* CheckContext::grid_products() {
  if (!d_.grid.present || !timeline_checks) return nullptr;  // grids are of the present-day model, not of restored states
  if (!grid_products_) grid_products_ = std::make_unique<GridProducts>(make_grid_products(d_, m_));
  return grid_products_.get();
}
void CheckContext::add_grid_property(repr::CellProperty p) {
  if (grid_products() && grid_products_) grid_products_->props.props.push_back(std::move(p));
}
const ContactGraph& CheckContext::contact_graph() {
  if (!graph_) graph_ = std::make_unique<ContactGraph>(stratonix::contact_graph(m_, topology_spacing()));
  return *graph_;
}
double CheckContext::topology_spacing() const {
  if (d_.voxel_spacing.present) {
    auto v = d_.voxel_spacing.si_list();
    if (!v.empty() && v[0] > 0) return v[0];
  }
  const auto& dom = m_.info().domain;
  return std::max(grid_.dx, std::max({dom.sx, dom.sy, dom.z_top - dom.z_base}) / 64.0);
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
        // a unit whose horizons are both held at the domain's top or base lies outside the domain there
        const double eps = 1e-6 * (inf.domain.z_top - inf.domain.z_base);
        if (zm >= inf.domain.z_top - eps || zm <= inf.domain.z_base + eps) continue;
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
      // Stage 4B: a sample whose probes straddle another fault's plane (fault networks meeting at
      // junctions, polygonal systems) measures that fault, not this one: it is set aside.
      std::size_t straddled = 0;
      if (inf.faults.size() > 1 && !at.empty()) {
        const double eps = 0.25;  // the probe offsets of scalar_jumps (default eps)
        const std::size_t n = at.size();
        std::vector<Vec3> q(4 * n);
        for (std::size_t i = 0; i < n; ++i) {
          q[i] = at[i] + inf.faults[f].normal * eps;
          q[n + i] = at[i] + inf.faults[f].normal * (2 * eps);
          q[2 * n + i] = at[i] - inf.faults[f].normal * eps;
          q[3 * n + i] = at[i] - inf.faults[f].normal * (2 * eps);
        }
        EvalResult er;
        c.model().evaluate(q, kEvalBlock, er);
        const BlockCode other = ~(BlockCode(1) << f);
        std::vector<Vec3> kept;
        for (std::size_t i = 0; i < n; ++i) {
          const BlockCode b0 = er.block[i] & other;
          if ((er.block[n + i] & other) != b0 || (er.block[2 * n + i] & other) != b0 || (er.block[3 * n + i] & other) != b0) { ++straddled; continue; }
          kept.push_back(at[i]);
        }
        at.swap(kept);
      }

      double worst = 0;
      std::string hn;
      if (!at.empty()) {
        auto j = repr::scalar_jumps(c.model(), at, inf.faults[f].normal);
        for (std::size_t h = 0; h < inf.horizons.size(); ++h)
          for (std::size_t k = 0; k < at.size(); ++k)
            if (std::fabs(j[h * at.size() + k]) > worst) {
              worst = std::fabs(j[h * at.size() + k]);
              hn = inf.horizons[h].name;
            }
      }
      auto r = base(c, "C3", "no offset outside the fault's tips", inf.faults[f].name);
      r.value = worst;
      r.tolerance = tol;
      r.tolerance_source = src;
      r.passed = worst <= tol;
      r.detail = std::to_string(at.size()) + " samples just outside the tip line" + (hn.empty() ? "" : "; largest at " + hn) +
                 (straddled ? "; " + std::to_string(straddled) + " set aside on another fault's plane" : "");
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
      if (f < c.description().faults.size() && c.description().faults[f].linking) continue;  // offset is the linked faults' throw
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
        std::size_t n = 0, set_aside = 0;
        for (int k = 1; k < NA; ++k) {
          double a = -0.5 * F.length + F.length * k / NA;
          double bhw, bfw;
          if (!cutoff(h, a, true, bhw) || !cutoff(h, a, false, bfw)) continue;
          Vec3 chw = F.at(a, bhw), cfw = F.at(a, bfw);
          {
            // stage 4B: another fault's plane between the two cut-offs (networks, polygonal systems)
            // adds its own offset to the separation; such a position measures two faults: set aside
            bool crossed = false;
            const Vec3 mid = (chw + cfw) * 0.5;
            for (std::size_t g = 0; g < net->size() && !crossed; ++g) {
              if (g == f) continue;
              const Fault& G = net->faults()[g];
              if (G.curved || std::fabs(dot(mid - G.p0, G.s)) >= 0.5 * G.length) continue;
              // between the cut-offs, or close enough to one that its offset reaches the cut-off search
              // the neighbour's vertical offset (none for strike-slip) is what reaches a cut-off
              const double near = std::max(5.0, 2 * std::fabs(G.tmax * G.slip.z)), s1 = G.sigma(chw), s2 = G.sigma(cfw);
              crossed = (s1 > 0) != (s2 > 0) || std::fabs(s1) < near || std::fabs(s2) < near;
            }
            if (crossed) { ++set_aside; continue; }
          }
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
        r.detail = std::to_string(n) + " strike positions (fault frame); profile maximum " + fmt_fixed(tmax, 2) + " m" +
                   (set_aside ? "; " + std::to_string(set_aside) + " set aside where another fault crosses between or near the cut-offs" : "");
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
      if (fd.linking) continue;  // stage 4B: a connecting fault has no slip of its own
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
  const auto st = stress_input(d, c.model());
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
      r.tolerance_source = "Stratonix rule (exact)";
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
      r.tolerance_source = "Stratonix rule: bounded by definition";
      r.passed = bad == 0;
      r.detail = std::to_string(n) + " nodes; largest Ts / Ts_max " + fmt_fixed(worst, 3);
      out.push_back(r);
    }
    return out;
  }
};

/// C21: data misfit within the stated uncertainty (STX-TN-11 section 7).
class C21 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    auto* mis = static_cast<const std::vector<DataMisfit>*>(c.model().capability("data_misfit"));
    if (!mis || std::none_of(mis->begin(), mis->end(), [](const DataMisfit& r) { return !std::isnan(r.sigma); })) return out;
    std::size_t n = 0, within = 0, exact = 0, exact_bad = 0;
    double worst = 0;
    const double tol = c.description().tol_wells.present ? c.description().tol_wells.si() : 0.01;  // Stratonix default (C5)
    for (auto& r : *mis) {
      if (std::isnan(r.sigma)) continue;
      if (r.sigma > 0) {
        ++n;
        within += std::fabs(r.normalised) <= 3.0;
        worst = std::max(worst, std::fabs(r.normalised));
      } else if (r.kind != "orientation") {
        ++exact;
        exact_bad += std::fabs(r.misfit) > tol;
      }
    }
    auto r = base(c, "C21", "data misfit within uncertainty", "all data with a stated uncertainty");
    r.unit = "fraction";
    r.value = n ? static_cast<double>(within) / static_cast<double>(n) : 1.0;
    r.tolerance = 0.95;
    r.tolerance_source = "Stratonix default: at least 95 % within 3 sigma";
    r.passed = r.value >= 0.95 && exact_bad == 0;
    r.detail = std::to_string(within) + " of " + std::to_string(n) + " within 3 sigma (largest |r| / sigma " + fmt_fixed(worst, 2) + ")" +
               (exact ? "; " + std::to_string(exact_bad) + " of " + std::to_string(exact) + " exact data beyond " + fmt_fixed(tol, 3) + " m" : std::string());
    out.push_back(r);
    return out;
  }
};

/// C22: conditioning converged inside the admissible ranges (STX-TN-11 section 7).
class C22 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    auto* cr = static_cast<const ConditioningResult*>(c.model().capability("conditioning"));
    if (!cr) return out;
    for (auto& p : cr->parameters) {
      auto r = base(c, "C22", "conditioning converged inside the range", p.path);
      r.unit = p.unit;
      r.value = p.value;
      r.tolerance = 0;
      r.tolerance_source = "Stratonix rule: strictly inside [" + fmt_shortest(p.min) + ", " + fmt_shortest(p.max) + "] " + p.unit;
      r.passed = cr->converged && p.value > p.min && p.value < p.max;
      r.detail = "fitted " + fmt_fixed(p.value, 3) + " +- " + fmt_fixed(p.sigma, 3) + " " + p.unit + " after " + std::to_string(cr->iterations) +
                 " iterations; reduced chi2 " + fmt_fixed(cr->chi2_reduced, 3);
      out.push_back(r);
    }
    if (cr->chi2_reduced > 4) {
      auto r = base(c, "C22", "conditioning reduced chi2", "fit");
      r.warn_only = true;
      r.value = cr->chi2_reduced;
      r.tolerance = 4;
      r.tolerance_source = "Stratonix default: warning above 4";
      r.passed = false;
      r.detail = "the stated uncertainties are too small for the misfit, or the model lacks a mechanism";
      out.push_back(r);
    }
    return out;
  }
};

/// C28: posterior convergence (STX-TN-14 section 4): every parameter's split R-hat below 1.01 and
/// effective sample size at least 100.
class C28 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    auto* pr = static_cast<const PosteriorResult*>(c.model().capability("posterior"));
    if (!pr) return out;
    for (std::size_t k = 0; k < pr->paths.size(); ++k) {
      auto r = base(c, "C28", "posterior split R-hat", pr->paths[k]);
      r.unit = "ratio";
      r.value = pr->rhat[k];
      r.tolerance = 1.01;
      r.tolerance_source = "Stratonix default: split R-hat below 1.01 (Vehtari et al. 2021)";
      r.passed = pr->rhat[k] < 1.01;
      r.detail = std::to_string(pr->walkers) + " walkers x " + std::to_string(pr->steps) + " steps after " + std::to_string(pr->burn_in) + " burn-in, seed " +
                 std::to_string(pr->seed) + (r.passed ? "" : ": run the sampler longer (conditioning.sampler.steps)");
      out.push_back(r);
      auto e = base(c, "C28", "posterior effective sample size", pr->paths[k]);
      e.unit = "samples";
      e.value = pr->ess[k];
      e.tolerance = 100;
      e.tolerance_source = "Stratonix default: at least 100 effective samples per parameter";
      e.passed = pr->ess[k] >= 100;
      e.detail = "integrated autocorrelation time " + fmt_fixed(static_cast<double>(pr->walkers) * pr->steps / std::max(1.0, pr->ess[k]), 1) + " steps";
      out.push_back(e);
    }
    return out;
  }
};

/// C26: topology (STX-TN-14 section 6): units in contact across a horizon inside one fault block are
/// adjacent in the column, unless an unconformity (a change of series) or a terminating unit lies
/// between them. Every explained exception is listed.
class C26 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    const auto& g = c.contact_graph();
    const auto& inf = c.model().info();
    std::vector<std::string> bad, explained;
    for (auto& e : g.edges) {
      const auto& a = g.nodes[static_cast<std::size_t>(e.a)];
      const auto& b = g.nodes[static_cast<std::size_t>(e.b)];
      if (a.block != b.block || std::abs(a.unit - b.unit) <= 1 || e.surface.rfind("body:", 0) == 0) continue;  // bodies cut across by design
      auto nm = [&](int u) { return (u >= 0 && static_cast<std::size_t>(u) < inf.units.size()) ? inf.units[static_cast<std::size_t>(u)].name : std::to_string(u); };
      std::string t = nm(std::min(a.unit, b.unit)) + " on " + nm(std::max(a.unit, b.unit)) + " in block " + std::to_string(a.block) + " (" +
                      fmt_fixed(e.area / 1e6, 3) + " km2)";
      (e.surface.rfind("unconformity:", 0) == 0 ? explained : bad).push_back(t);
    }
    auto r = base(c, "C26", "topology: units in contact are adjacent in the column", std::to_string(g.nodes.size()) + " nodes, " + std::to_string(g.edges.size()) + " contacts");
    r.unit = "contacts";
    r.value = static_cast<double>(bad.size());
    r.tolerance = 0;
    r.tolerance_source = "Stratonix rule: no unexplained contact (lattice " + fmt_shortest(g.spacing) + " m)";
    r.passed = bad.empty();
    std::string det;
    for (std::size_t k = 0; k < bad.size() && k < 5; ++k) det += (k ? "; " : "") + bad[k];
    if (!bad.empty()) det = "unexplained: " + det + (bad.size() > 5 ? " and " + std::to_string(bad.size() - 5) + " more" : "");
    if (!explained.empty()) {
      std::string ex;
      for (std::size_t k = 0; k < explained.size() && k < 5; ++k) ex += (k ? "; " : "") + explained[k];
      det += std::string(det.empty() ? "" : " | ") + "unconformities or terminations: " + ex + (explained.size() > 5 ? " and more" : "");
    }
    r.detail = det.empty() ? "every contact inside a block is between adjacent units" : det;
    return {r};
  }
};

/// C27: watertight horizon surfaces (STX-TN-14 section 6): every horizon edge on the surface grid lies
/// on a fault, another horizon (unconformity, pinch-out) or the domain boundary within tolerance.
class C27 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    // the exported surfaces are those of the present-day model; restored states are checked by C17 and C26
    if (!c.timeline_checks) return {};
    std::string src;
    double tol = 0.5 * c.grid().dx;
    if (c.description().tol_watertight.present) {
      tol = c.description().tol_watertight.si();
      src = std::string("checks.tolerances.watertight (") + desc::evidence_name(c.description().tol_watertight.evidence) + ")";
    } else {
      src = "Stratonix default: half the surface spacing";
    }
    auto w = watertightness(c.model(), c.grid(), [&](int h) -> const repr::FaultedSurface& { return c.faulted_surface(h); }, tol);
    auto r = base(c, "C27", "watertight horizon surfaces", "largest gap");
    r.unit = "m";
    r.value = w.max_gap;
    r.tolerance = tol;
    r.tolerance_source = src;
    r.passed = w.ok();
    std::string det;
    for (std::size_t k = 0; k < w.gaps.size() && k < 3; ++k)
      det += (k ? "; " : "") + w.gaps[k].surface + " " + fmt_fixed(w.gaps[k].gap, 2) + " m at (" + fmt_fixed(w.gaps[k].at.x, 0) + ", " +
             fmt_fixed(w.gaps[k].at.y, 0) + ", " + fmt_fixed(w.gaps[k].at.z, 0) + ")";
    r.detail = det.empty() ? "every horizon edge ends on a fault, a horizon or the domain boundary" : det;
    return {r};
  }
};

/// C30: bodies (STX-TN-15 section 1): every body has a volume inside the model and an emplacement
/// age; bodies thinner than half the grid spacing are reported.
class C30 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    // stage 4B: geobody sets (channel belts, lobes) - their net-to-gross in the host unit
    if (auto* gs = static_cast<const std::vector<GeobodyReport>*>(c.model().capability("geobodies")))
      for (auto& g : *gs) {
        auto r = base(c, "C30", "geobody net-to-gross", g.name + " (" + g.type + " in " + g.host + ")");
        r.unit = "fraction";
        r.value = g.ntg;
        r.tolerance = 1;
        r.tolerance_source = "Stratonix rule: a set occupies part of its unit (0 < net-to-gross < 1)";
        r.passed = g.ntg > 0 && g.ntg < 1;
        r.detail = std::to_string(g.count) + " placed in the unit's geological-time box (" + fmt_fixed(g.box_thickness, 1) +
                   " m high); closed form without overlap or clipping " + fmt_fixed(g.ntg_closed, 4);
        out.push_back(r);
      }
    auto* bs = static_cast<const std::vector<BodyReport>*>(c.model().capability("bodies"));
    if (!bs) return out;
    for (auto& b : *bs) {
      auto r = base(c, "C30", "body volume", b.name + " (" + b.type + ", age " + fmt_shortest(b.age) + ")");
      r.unit = "Mm3";
      r.value = b.volume / 1e6;
      r.tolerance = 0;
      r.tolerance_source = "Stratonix rule: a body has a positive volume";
      r.passed = b.volume > 0 && std::isfinite(b.age);
      r.detail = b.note.empty() ? "placed in age order; the host is evaluated where the body is not" : b.note;
      out.push_back(r);
    }
    return out;
  }
};

/// C31: restoration strain and bed length (STX-TN-15 section 4). The volumetric restoration's
/// largest strain and the restored horizon's largest in-plane stretch against their limits; and,
/// in models with fault-parallel flow or a restoration block, the bed length of every horizon along
/// every section, present against restored (I-026).
class C31 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const auto& d = c.description();
    const double blim = d.restoration.bed_length_limit.present ? d.restoration.bed_length_limit.si() : 0.05;
    const std::string bsrc = d.restoration.bed_length_limit.present ? "restoration.bed_length_limit" : "Stratonix default (5 %)";
    if (const VolumetricResult* v = c.volumetric) {
      const double lim = d.restoration.max_strain.present ? d.restoration.max_strain.si() : 0.1;
      auto r = base(c, "C31", "restoration strain", d.restoration.horizon + " to its datum (volumetric)");
      r.unit = "fraction";
      r.value = v->max_strain;
      r.tolerance = lim;
      r.tolerance_source = d.restoration.max_strain.present ? "restoration.volumetric.max_strain" : "Stratonix default (10 %)";
      r.passed = v->ok && v->max_strain <= lim;
      r.detail = v->ok ? "Green-Lagrange, largest principal; mean " + fmt_fixed(v->mean_strain, 5) + "; restored horizon flat within " +
                             fmt_fixed(v->flatness_rms, 3) + " m rms; " + std::to_string(v->tetrahedra) + " tetrahedra, " +
                             std::to_string(v->fault_contacts) + " fault contacts"
                       : "not restored: " + v->reason;
      out.push_back(r);
      if (v->ok) {
        auto b = base(c, "C31", "restored bed length", d.restoration.horizon + " (volumetric)");
        b.unit = "fraction";
        b.value = v->horizon_stretch_max;
        b.tolerance = blim;
        b.tolerance_source = bsrc;
        b.passed = v->horizon_stretch_max <= blim;
        b.detail = "largest in-plane stretch of the restored horizon, present against restored";
        out.push_back(b);
      }
    }
    // bed length along sections (fault-parallel flow, I-026, or any restoration block)
    bool fpf = false;
    if (auto* net = static_cast<const FaultNetwork*>(c.model().capability("faults")))
      for (auto& f : net->faults()) fpf = fpf || (f.curved && f.curved->fpf);
    if (!fpf && !d.restoration.present) return out;
    const auto& inf = c.model().info();
    for (auto& sd : d.sections) {
      if (!sd.polyline.present || sd.polyline.rows() < 2) continue;
      Polyline pl;
      for (std::size_t q = 0; q < sd.polyline.rows(); ++q) pl.push_back({sd.polyline.si(2 * q), sd.polyline.si(2 * q + 1), 0});
      auto s = repr::make_section(c.model(), sd.name, pl, std::max(5.0, c.fault_spacing() / 4), 2.0);
      double worst = 0;
      std::string hn;
      for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
        double present = 0, restored = 0;
        for (auto& seg : repr::horizon_segments(c.model(), s, static_cast<int>(h))) {
          if (seg.size() < 2) continue;
          std::vector<Vec3> r(seg.size());
          c.model().restore(seg, r);
          for (std::size_t q = 1; q < seg.size(); ++q) {
            present += norm(seg[q] - seg[q - 1]);
            restored += norm(r[q] - r[q - 1]);
          }
        }
        if (!(restored > 0)) continue;
        const double ch = present / restored - 1;
        if (std::fabs(ch) > std::fabs(worst)) { worst = ch; hn = inf.horizons[h].name; }
      }
      if (hn.empty()) continue;
      auto r = base(c, "C31", "bed length along a section", sd.name);
      r.unit = "fraction";
      r.value = std::fabs(worst);
      r.tolerance = blim;
      r.tolerance_source = bsrc;
      r.passed = std::fabs(worst) <= blim;
      r.detail = "largest change, present against restored: " + fmt_fixed(100 * worst, 3) + " % at " + hn + (fpf ? " (fault-parallel flow)" : "");
      out.push_back(r);
    }
    return out;
  }
};

/// C32 heat balance, C33 pressure bounds, C34 fracture network (STX-TN-15 sections 5 and 6).
class C32 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    if (!c.basin) return out;
    for (auto& [loc, err] : c.basin->heat_balance) {
      auto r = base(c, "C32", "heat balance", loc);
      r.unit = "fraction";
      r.value = err;
      r.tolerance = 0.01;
      r.tolerance_source = "Stratonix rule: in steady state the surface heat flow is the basal flow plus the radiogenic heat (1 %)";
      r.passed = err <= 0.01;
      r.detail = "present column in steady state, finite volumes";
      out.push_back(r);
    }
    for (auto& b : c.basin->pressure) {
      auto r = base(c, "C33", "pore pressure between hydrostatic and lithostatic", b.location);
      r.unit = "MPa";
      r.value = b.worst / 1e6;
      r.tolerance = 0;
      r.tolerance_source = "Stratonix rule: hydrostatic <= pore pressure <= lithostatic";
      r.passed = b.below == 0 && b.above == 0;
      r.detail = b.below + b.above == 0 ? "present profile within its bounds"
                                         : std::to_string(b.below) + " nodes below hydrostatic, " + std::to_string(b.above) +
                                               " above lithostatic: inputs the model cannot honour";
      out.push_back(r);
    }
    if (c.basin->dfn) {
      auto r = base(c, "C34", "fracture network", c.basin->dfn_unit);
      r.unit = "fraction";
      r.tolerance = 0.1;
      r.tolerance_source = "Stratonix rule: realised P32 within 10 % of the target (warning)";
      if (!c.basin->dfn_ok) {
        r.passed = false;
        r.detail = "not sampled: " + c.basin->dfn_reason;
      } else {
        r.value = c.basin->dfn_target > 0 ? std::fabs(c.basin->dfn_realised / c.basin->dfn_target - 1) : 0.0;
        r.warn_only = true;
        r.passed = r.value <= 0.1;
        r.detail = std::to_string(c.basin->dfn_count) + " disks; realised P32 " + fmt_fixed(c.basin->dfn_realised, 4) + " 1/m against " +
                   fmt_fixed(c.basin->dfn_target, 4) + " 1/m";
      }
      out.push_back(r);
    }
    return out;
  }
};

/// C35 well data and C38 age consistency (STX-TN-16 section 6).
class C35 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const auto* wd = c.well_data;
    if (!wd) return out;
    for (auto& w : wd->wells) {
      auto r = base(c, "C35", "well data", w.well);
      r.unit = "m";
      r.value = std::max(w.beyond, w.tops_off_path);
      r.tolerance = 1.0;
      r.tolerance_source = "Stratonix rule: logs and tops within the well path (1 m)";
      r.passed = r.value <= r.tolerance;
      std::string det = w.survey + "; " + std::to_string(w.logs) + " log file" + (w.logs == 1 ? "" : "s") + ", " + std::to_string(w.curves) + " curves";
      if (w.beyond > 1.0) det += "; logs run " + fmt_fixed(w.beyond, 1) + " m beyond the path's total depth";
      if (w.tops_off_path > 1.0) det += "; a top lies " + fmt_fixed(w.tops_off_path, 1) + " m below the path's total depth";
      for (std::size_t k = 0; k < w.notes.size() && k < 4; ++k) det += "; " + w.notes[k];
      if (w.notes.size() > 4) det += " and " + std::to_string(w.notes.size() - 4) + " more (well_data_report.txt)";
      r.detail = det;
      out.push_back(r);
    }
    if (wd->ages) {
      auto r = base(c, "C38", "age consistency", "wells");
      r.unit = "conflicts";
      r.value = static_cast<double>(wd->conflicts.size());
      r.tolerance = 0;
      r.tolerance_source = "Stratonix rule: no age inversion in a well; tops agree with the age data within their uncertainty";
      r.passed = wd->conflicts.empty();
      std::string det = std::to_string(wd->dated) + " age data";
      for (std::size_t k = 0; k < wd->conflicts.size() && k < 4; ++k) det += "; " + wd->conflicts[k].well + ": " + wd->conflicts[k].what;
      if (wd->conflicts.size() > 4) det += " and " + std::to_string(wd->conflicts.size() - 4) + " more";
      r.detail = det;
      out.push_back(r);
    }
    return out;
  }
};

/// C39 property conditioning and C40 reproduced statistics (STX-TN-17 section 8).
class C39 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const auto* gs = c.geostat;
    if (!gs) return out;
    for (auto& k : gs->conditioning) {
      auto r = base(c, "C39", "property conditioning", k.property);
      r.unit = "cells";
      r.value = static_cast<double>(k.mismatched);
      if (k.exact) {
        r.tolerance = 0;
        r.tolerance_source = "Stratonix rule: simulated values equal the upscaled logs at the well cells";
        r.passed = k.mismatched == 0;
      } else {  // object facies: honoured by placement, misfit reported (STX-TN-17 section 1.5)
        r.tolerance = std::floor(0.05 * static_cast<double>(k.cells));
        r.tolerance_source = "Stratonix default: object facies honour at least 95 % of the well cells";
        r.passed = r.value <= r.tolerance;
      }
      r.detail = std::to_string(k.cells) + " well cells, " + std::to_string(k.mismatched) + " not honoured" +
                 (k.exact && k.cells ? "; largest difference " + fmt_shortest(k.largest) : std::string());
      out.push_back(r);
    }
    for (auto& st : gs->statistics) {
      auto r = base(c, "C40", "reproduced statistics", st.property + ": " + st.what);
      r.unit = st.what.rfind("proportion", 0) == 0 ? "fraction" : "distance";
      r.value = st.what.rfind("proportion", 0) == 0 ? std::fabs(st.got - st.target) : st.got;
      r.tolerance = st.tolerance;
      r.tolerance_source = "Stratonix default (one realisation; a warning, not a failure)";
      r.passed = r.value <= r.tolerance;
      r.warn_only = true;
      r.detail = st.what.rfind("proportion", 0) == 0 ? fmt_fixed(st.got, 3) + " against the target " + fmt_fixed(st.target, 3)
                                                   : "Kolmogorov-Smirnov distance " + fmt_fixed(st.got, 3);
      out.push_back(r);
    }
    if (gs->cells) {
      auto r = base(c, "C40", "reproduced statistics", "cells in geological time");
      r.unit = "fraction";
      r.value = static_cast<double>(gs->fallback) / static_cast<double>(gs->cells);
      r.tolerance = 0.1;
      r.tolerance_source = "Stratonix default: at most 10 % of the cells simulated in present-day coordinates";
      r.passed = r.value <= r.tolerance;
      r.warn_only = true;
      r.detail = std::to_string(gs->fallback) + " of " + std::to_string(gs->cells) + " cells outside the geological-time mapping";
      out.push_back(r);
    }
    return out;
  }
};

/// C36 time-depth ties (STX-TN-16 section 6): each well top after depth conversion, per layer.
class C36 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const auto* td = c.time_depth;
    if (!td) return out;
    std::map<std::string, std::vector<const CheckContext::TimeDepthResults::Tie*>> by;
    for (auto& t : td->ties) by[t.layer].push_back(&t);
    for (auto& [layer, ts] : by) {
      auto r = base(c, "C36", "time-depth ties", layer);
      r.unit = "m";
      double worst = 0, loo = 0;
      std::string where, blind;
      for (auto* t : ts) {
        if (std::fabs(t->residual) > std::fabs(worst)) { worst = t->residual; where = t->well; }
        if (std::isfinite(t->loo) && std::fabs(t->loo) > std::fabs(loo)) loo = t->loo;
        if (t->excluded) blind += std::string(blind.empty() ? "" : ", ") + t->well + " " + fmt_fixed(t->residual, 2) + " m";
      }
      r.value = std::fabs(worst);
      r.tolerance = td->tolerance;
      r.tolerance_source = td->tolerance_source;
      r.passed = r.value <= r.tolerance;
      r.detail = std::to_string(ts.size()) + " tops; largest residual " + fmt_fixed(worst, 2) + " m" + (where.empty() ? "" : " at " + where) +
                 "; leave-one-out (v0 from the other wells) up to " + fmt_fixed(std::fabs(loo), 2) + " m" + (blind.empty() ? "" : "; blind wells: " + blind);
      out.push_back(r);
    }
    return out;
  }
};

/// C29: potential-field misfit (STX-TN-14 section 7): reduced chi-square of the stations reported;
/// a warning above 2.
class C29 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const PotentialFieldResult* pf = c.potential_fields;
    if (!pf) return out;
    for (int kind = 0; kind < 2; ++kind) {
      const std::size_t n = kind == 0 ? pf->n_gravity : pf->n_magnetics;
      if (n == 0) continue;
      const double chi = kind == 0 ? pf->chi2_gravity : pf->chi2_magnetics;
      auto r = base(c, "C29", "potential-field misfit", std::string(kind == 0 ? "gravity" : "magnetics") + ", " + std::to_string(n) + " stations");
      r.unit = "chi2";
      r.warn_only = true;
      r.value = chi;
      r.tolerance = 2;
      r.tolerance_source = "Stratonix default: warning above a reduced chi-square of 2";
      double rms = 0;
      std::size_t m = 0;
      for (auto& st : pf->stations)
        if (st.kind == (kind == 0 ? "gravity" : "magnetics")) { rms += std::pow(st.modelled - st.observed, 2); ++m; }
      rms = m ? std::sqrt(rms / static_cast<double>(m)) : 0;
      r.passed = !std::isfinite(chi) || chi <= 2;
      r.detail = std::string("rms misfit ") + fmt_fixed(rms, 4) + (kind == 0 ? " mGal" : " nT") +
                 (std::isfinite(chi) ? "" : "; no station uncertainty stated, so no chi-square");
      out.push_back(r);
    }
    return out;
  }
};

/// C23 to C25: grid checks (STX-TN-12 section 7), when the description asks for a grid.
class C23 final : public Check {
  std::vector<CheckResult> run(CheckContext& c) const override {
    std::vector<CheckResult> out;
    const GridProducts* gp = c.grid_products();
    if (!gp) return out;
    const auto& g = gp->grid;
    std::size_t act = 0;
    for (auto a : g.actnum) act += a;
    const std::string subj = std::to_string(g.nx) + " x " + std::to_string(g.ny) + " x " + std::to_string(g.nz) + " grid, " + std::to_string(act) + " active cells";
    auto r = base(c, "C23", "positive cell volumes", subj);
    r.unit = "cells";
    r.value = static_cast<double>(gp->defects.nonpositive_volume);
    r.tolerance = 0;
    r.tolerance_source = "Stratonix rule: none";
    r.passed = gp->defects.nonpositive_volume == 0;
    r.detail = r.passed ? "every active cell has a positive volume" : "first: " + gp->defects.first_nonpositive;
    out.push_back(r);
    auto r2 = base(c, "C24", "no inverted cells", subj);
    r2.unit = "cells";
    r2.value = static_cast<double>(gp->defects.inverted);
    r2.tolerance = 0;
    r2.tolerance_source = "Stratonix rule: none";
    r2.passed = gp->defects.inverted == 0;
    r2.detail = r2.passed ? "top corners at or above bottom corners; positive Jacobian at every corner" : "first: " + gp->defects.first_inverted;
    out.push_back(r2);
    auto r3 = base(c, "C25", "fault faces on the fault surfaces", std::to_string(gp->defects.fault_faces) + " fault faces on fault pillars");
    r3.unit = "m";
    r3.value = gp->defects.worst_fault_distance;
    r3.tolerance = gp->fault_tolerance;
    r3.tolerance_source = "Stratonix default: 1 % of the areal cell size (plus the chord of a curved fault)";
    r3.passed = gp->defects.fault_face_off == 0;
    r3.detail = r3.passed ? "largest distance of a fault-face corner from its fault " + fmt_fixed(gp->defects.worst_fault_distance, 4) + " m"
                          : std::to_string(gp->defects.fault_face_off) + " faces off; first: " + gp->defects.first_fault_off;
    out.push_back(r3);
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
STRATONIX_REGISTER(check_c21, registry::Kind::Check, "C21_data_misfit", "1.0.0", "Data misfit within the stated uncertainty",
                   [] { return std::make_unique<C21>(); })
STRATONIX_REGISTER(check_c22, registry::Kind::Check, "C22_conditioning", "1.0.0", "Conditioning converged inside the admissible ranges",
                   [] { return std::make_unique<C22>(); })
STRATONIX_REGISTER(check_c26, registry::Kind::Check, "C26_topology", "1.0.0", "Topology: units in contact are adjacent, or an unconformity separates them",
                   [] { return std::make_unique<C26>(); })
STRATONIX_REGISTER(check_c27, registry::Kind::Check, "C27_watertight", "1.0.0", "Watertight horizon surfaces",
                   [] { return std::make_unique<C27>(); })
STRATONIX_REGISTER(check_c30, registry::Kind::Check, "C30_bodies", "1.0.0", "Bodies: positive volume, emplacement age, thin bodies reported",
                   [] { return std::make_unique<C30>(); })
STRATONIX_REGISTER(check_c32, registry::Kind::Check, "C32_basin", "1.0.0", "Heat balance, pore-pressure bounds and fracture-network intensity",
                   [] { return std::make_unique<C32>(); })
STRATONIX_REGISTER(check_c35, registry::Kind::Check, "C35_well_data", "1.0.0", "Well data: survey method stated, logs and tops on the well path; age consistency (C38)",
                   [] { return std::make_unique<C35>(); })
STRATONIX_REGISTER(check_c36, registry::Kind::Check, "C36_time_depth", "1.0.0", "Time-depth ties: well tops after depth conversion within the tolerance",
                   [] { return std::make_unique<C36>(); })
STRATONIX_REGISTER(check_c39, registry::Kind::Check, "C39_geostatistics", "1.0.0", "Geostatistical properties: conditioning at the well cells (C39) and reproduced statistics (C40)",
                   [] { return std::make_unique<C39>(); })
STRATONIX_REGISTER(check_c31, registry::Kind::Check, "C31_restoration", "1.0.0", "Restoration strain and bed length within their limits",
                   [] { return std::make_unique<C31>(); })
STRATONIX_REGISTER(check_c29, registry::Kind::Check, "C29_potential_fields", "1.0.0", "Potential-field misfit: reduced chi-square of the stations (warning above 2)",
                   [] { return std::make_unique<C29>(); })
STRATONIX_REGISTER(check_c28, registry::Kind::Check, "C28_posterior", "1.0.0", "Posterior converged: split R-hat below 1.01, effective sample size at least 100",
                   [] { return std::make_unique<C28>(); })
STRATONIX_REGISTER(check_c23, registry::Kind::Check, "C23_grid", "1.0.0", "Grid: positive volumes, no inverted cells, fault faces on the faults",
                   [] { return std::make_unique<C23>(); })
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
int model_unit_of(const desc::Description& d, const Model& m, const std::string& unit) {
  const int u = d.unit_index(unit);
  if (u < 0) return -1;
  const int ht = d.horizon_index(d.units[static_cast<std::size_t>(u)].top), hb = d.horizon_index(d.units[static_cast<std::size_t>(u)].base);
  const auto& mu = m.info().units;
  for (std::size_t k = 0; k < mu.size(); ++k)
    if (mu[k].top_horizon == ht && mu[k].base_horizon == hb && ht >= 0 && hb >= 0) return static_cast<int>(k);
  return -1;
}
repr::GridOptions grid_options(const desc::Description& d, const Model& m) {
  repr::GridOptions o;
  const auto& g = d.grid;
  auto sp = g.spacing.si_list();
  o.dx = sp.empty() ? d.grid_spacing.si() : sp[0];
  o.dy = sp.size() > 1 ? sp[1] : o.dx;
  if (g.margin.present) o.margin = g.margin.si();
  if (g.min_thickness.present) o.min_thickness = g.min_thickness.si();
  o.snap_to_faults = g.pillars != "vertical";
  if (g.zones.empty()) {
    const auto& mu = m.info().units;
    for (std::size_t k = 0; k < mu.size(); ++k)
      if (mu[k].top_horizon >= 0 && mu[k].base_horizon >= 0) o.zones.push_back({static_cast<int>(k), 1, repr::Layering::Proportional, 0});
  }
  for (auto& z : g.zones) {
    repr::ZoneSpec s;
    s.unit = model_unit_of(d, m, z.unit);
    s.layers = z.layers.present ? static_cast<int>(std::llround(z.layers.si())) : 1;
    s.layering = z.layering == "follow_top" ? repr::Layering::FollowTop : z.layering == "follow_base" ? repr::Layering::FollowBase : repr::Layering::Proportional;
    if (z.layer_thickness.present) s.layer_thickness = z.layer_thickness.si();
    o.zones.push_back(s);
  }
  return o;
}
repr::PropertyInputs property_inputs(const desc::Description& d, const Model& m) {
  repr::PropertyInputs in;
  const auto& g = d.grid;
  for (auto& f : g.facies)
    in.facies.push_back({model_unit_of(d, m, f.unit), static_cast<int>(std::llround(f.code.si())), f.name, desc::evidence_name(f.code.evidence), f.code.source});
  in.porosity_from_compaction = g.porosity == "compaction";
  for (auto& t : g.transforms) {
    // a and b as written (their units are those of the law); evidence of the weaker of the two
    std::string ev = desc::evidence_name(t.a.evidence);
    if (t.b.evidence != t.a.evidence) ev += " / " + std::string(desc::evidence_name(t.b.evidence));
    in.transforms.push_back({t.name, t.unit, t.from, t.law, t.a.values.empty() ? kNaN : t.a.values[0], t.b.values.empty() ? kNaN : t.b.values[0], ev,
                             t.a.source == t.b.source ? t.a.source : t.a.source + "; " + t.b.source});
  }
  return in;
}
GridProducts make_grid_products(const desc::Description& d, const Model& m) {
  GridProducts gp;
  gp.options = grid_options(d, m);
  auto& reg = registry::Registry::global();
  gp.grid = reg.create<repr::CornerPointExtractor>(registry::Kind::Extractor, "corner_point")->extract(m, gp.options);
  gp.props = repr::grid_properties(m, gp.grid, property_inputs(d, m));
  if (!d.grid.transmissibility.empty()) {
    gp.trans = repr::fault_transmissibility(m, gp.grid, gp.props, fault_rock_input(d, m));
    gp.has_trans = true;
  }
  if (d.grid.tetrahedral) {
    gp.tets = reg.create<repr::TetrahedralExtractor>(registry::Kind::Extractor, "tetrahedral")->extract(gp.grid);
    gp.has_tets = true;
  }
  gp.fault_tolerance = 0.01 * std::min(gp.grid.dx, gp.grid.dy);  // Stratonix default
  gp.defects = repr::grid_defects(m, gp.grid, gp.fault_tolerance);
  for (auto& f : m.info().faults) gp.fault_names.push_back(f.name);
  return gp;
}
repr::StressInput stress_input(const desc::Description& d, const Model& m) {
  repr::StressInput s = stress_input(d);
  if (!s.present || !d.mechanics.present) return s;
  // stage 4B: mechanical stratigraphy per model unit
  s.lateral = true;
  s.unit_mech.assign(m.info().units.size(), {kNaN, kNaN, kNaN});
  for (auto& u : d.mechanics.units) {
    const int k = model_unit_of(d, m, u.unit);
    if (k >= 0) s.unit_mech[static_cast<std::size_t>(k)] = {u.young.si(), u.poisson.si(), u.biot.present ? u.biot.si() : 1.0};
  }
  if (d.mechanics.strain_min.present) {
    s.e_h = d.mechanics.strain_min.si();
    s.e_H = d.mechanics.strain_max.si();
  }
  if (d.stress.shmax_azimuth.present) s.azimuth = d.stress.shmax_azimuth.si();
  return s;
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
  static const bool timing = std::getenv("STRATONIX_TIMING") != nullptr;  // diagnostic: seconds per check on stderr
  for (auto* e : registry::Registry::global().list(registry::Kind::Check)) {
    auto chk = registry::Registry::global().create<Check>(registry::Kind::Check, e->name);
    const auto t0 = std::chrono::steady_clock::now();
    auto rs = chk->run(ctx);
    if (timing)
      std::fprintf(stderr, "[stratonix timing] %s %.3f s\n", e->name.c_str(), std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
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
  STRATONIX_PLUGIN_ANCHOR(check_c21);
  STRATONIX_PLUGIN_ANCHOR(check_c22);
}

}  // namespace stratonix
