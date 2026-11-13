// Checks that fail a build. Each check states its tolerance and where the tolerance came from.
#include "stratonix/checks/checks.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/core/roots.hpp"
#include "stratonix/engines/faults.hpp"
#include "stratonix/engines/structure.hpp"

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
                     {"wells", &d_.tol_wells, 0.01}, {"unoffset", &d_.tol_unoffset, 0.05}};
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
      for (std::size_t q = 0; q < idx.size(); ++q) {
        if (r.block[q] != b) continue;
        double t = a.v[idx[q]] - z.v[idx[q]];
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
        Vec3 r1 = net->restore_from_frame(p + F.n * (sg * eps), fi, hw);
        Vec3 r2 = net->restore_from_frame(p + F.n * (sg * 2 * eps), fi, hw);
        double s1 = r1.z - st->z(h, r1.x, r1.y), s2 = r2.z - st->z(h, r2.x, r2.y);
        return 2 * s1 - s2;
      };
      auto cutoff = [&](std::size_t h, double a, bool hw, double& b_out) {
        double bu, bl;
        F.tips(a, bu, bl);
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
          net->undo(fi, cm, fixed_hw ? -1 : 1, &expect);
          double meas;
          if (F.kin->name() == "dip_slip") {
            meas = cfw.z - chw.z;
          } else {
            // fixed wall's horizon on the vertical through the moving cut-off
            bool hwfix = fixed_hw;
            auto sv = [&](double z) {
              Vec3 p{cm.x, cm.y, z};
              Vec3 r = net->restore_from_frame(p, fi, hwfix);
              return r.z - st->z(h, r.x, r.y);
            };
            double zf = expand_and_solve(sv, cm.z, 5.0, 1e4, 1e-7);
            meas = fixed_hw ? cm.z - zf : zf - cm.z;
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
    // area target is given in m2 (length unit squared) or km2 via unit 'ratio'? use the value as m2
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
STRATONIX_REGISTER(check_c9, registry::Kind::Check, "C9_closure", "1.0.0", "Closure, spill point and area against targets",
                   [] { return std::make_unique<C9>(); })

CheckReport run_checks(CheckContext& ctx) {
  CheckReport rep;
  for (auto* e : registry::Registry::global().list(registry::Kind::Check)) {
    auto chk = registry::Registry::global().create<Check>(registry::Kind::Check, e->name);
    auto rs = chk->run(ctx);
    rep.results.insert(rep.results.end(), rs.begin(), rs.end());
  }
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
}

}  // namespace stratonix
