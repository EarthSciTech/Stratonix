// Deformation history (stage 2E, STX-TN-09): the event timeline and the model at any event.
#include "stratonix/engines/history.hpp"

#include <algorithm>
#include <cmath>
#include <set>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/structure.hpp"

namespace stratonix {

namespace {
desc::Param make_param(double v, const std::string& unit, Dimension dim, const std::string& source) {
  desc::Param p;
  p.present = true;
  p.values = {v};
  p.unit = unit;
  p.dim = dim;
  p.to_si = 1.0;
  p.evidence = desc::Evidence::Calibrated;
  p.source = source;
  p.text = fmt_shortest(v) + " " + unit;
  return p;
}
void scale(desc::Param& p, double k) {
  for (auto& v : p.values) v *= k;
  p.text = p.text.empty() ? p.text : p.text + " x " + fmt_fixed(k, 4);
}
/// Completed fraction at time t of an event from `from` to `to` (Ma).
double completed(double from, double to, double t) {
  if (!(from > to)) return t <= from ? 1.0 : 0.0;
  return std::clamp((from - t) / (from - to), 0.0, 1.0);
}
}  // namespace

std::vector<double> event_ages(const desc::Description& d) {
  std::set<double> a;
  for (auto& e : d.events) {
    a.insert(e.start());
    a.insert(e.end());
  }
  return std::vector<double>(a.rbegin(), a.rend());  // oldest first
}

desc::Description state_description(const desc::Description& d0, const Model& full, double t) {
  if (d0.events.empty()) throw Error(ErrorCode::InvalidArgument, "a model state needs a description with an events section (STX-TN-09)");
  auto* st = static_cast<const Structure*>(full.capability("structure"));
  if (!st) throw Error(ErrorCode::InvalidArgument, "model states are built by the parametric engine");
  desc::Description d = d0;
  const auto& rd = st->reference_depths();
  const std::size_t nh = d.horizons.size();
  const std::string src = "state at " + fmt_shortest(t) + " Ma: reference depth of the full model";
  std::vector<bool> keep(nh);
  std::set<std::string> gone;
  for (std::size_t h = 0; h < nh; ++h) {
    keep[h] = !(d.horizons[h].age.present && d.horizons[h].age.si() < t);
    if (!keep[h]) gone.insert(d.horizons[h].name);
  }
  // Horizons: the full model's described reference depths replace well tops (present depths)
  // and thicknesses; the state's own mechanisms act on them.
  const desc::WellDesc* rw = d.ref_well.empty() ? nullptr : d.well(d.ref_well);
  std::vector<desc::HorizonDesc> hs;
  for (std::size_t h = 0; h < nh; ++h) {
    if (!keep[h]) continue;
    auto hd = d.horizons[h];
    bool from_top = false;
    if (rw)
      for (auto& [n, p] : rw->tops) from_top = from_top || n == hd.name;
    const bool iso_lost = !hd.iso_from.empty() && gone.count(hd.iso_from);
    if (from_top || hd.thickness_above.present || iso_lost) {
      hd.depth = make_param(rd[h], "m", Dimension::Length, src);
      hd.thickness_above = {};
      if (iso_lost) { hd.iso_from.clear(); hd.iso_trend.clear(); hd.iso_thickness = {}; }
    }
    hs.push_back(hd);
  }
  if (!hs.empty()) hs.front().intended_offset = true;  // the state's depositional surface may carry fault scarps
  d.horizons = hs;
  for (auto& w : d.wells) w.tops.clear();  // present depths do not hold in a restored state
  // Units, series, outputs referring to removed horizons.
  std::vector<desc::UnitDesc> us;
  for (auto& u : d.units)
    if (!gone.count(u.top) && !gone.count(u.base) && !(u.top == "model_top" && gone.count(d0.horizons.front().name))) us.push_back(u);
  d.units = us;
  std::vector<desc::CompactionUnitDesc> cus;
  for (auto& cu : d.compaction.units)
    if (d.unit_index(cu.unit) >= 0) cus.push_back(cu);
  d.compaction.units = cus;
  if (d.compaction.units.empty()) d.compaction = {};
  std::vector<desc::SeriesDesc> ss;
  for (auto sr : d.series) {
    std::vector<std::string> kept;
    for (auto& n : sr.horizons)
      if (!gone.count(n)) kept.push_back(n);
    if (kept.empty()) continue;
    sr.horizons = kept;
    ss.push_back(sr);
  }
  d.series = ss;
  std::vector<std::string> iso;
  for (auto& n : d.isopachs)
    if (d.unit_index(n) >= 0) iso.push_back(n);
  d.isopachs = iso;
  d.sheets.clear();
  if (gone.count(d.closure_horizon)) d.closure_horizon.clear();
  // Events: those not begun by t removed; those under way scaled to their completed part.
  std::vector<desc::EventDesc> evs;
  std::set<std::string> active_faults, active_prims, active_bodies;
  for (auto e : d.events) {
    if (e.type == "deposition") {
      if (!gone.count(e.ref)) evs.push_back(e);
      continue;
    }
    if (e.start() < t) continue;  // begins after t
    const double g = completed(e.start(), e.end(), t);
    if (g <= 0) continue;
    if (e.end() < t) {  // under way at t: complete at t, with its completed part
      e.to = make_param(t, "Ma", Dimension::Time, src);
      if (!e.from.present) e.from = e.age;
      if (e.type == "fault") {
        const auto& fd = d.faults[static_cast<std::size_t>(d.fault_index(e.ref))];
        if (!e.throw_.present) e.throw_ = fd.throw_max;
        scale(e.throw_, g);
        if (e.sense.empty()) e.sense = fd.sense.empty() ? "normal" : fd.sense;
      } else if (e.type == "fold") {
        for (auto& p : d.primitives)
          if (p.name == e.ref) {
            for (auto& [hn, r] : p.relief) scale(r, g);
            if (p.shortening.present) scale(p.shortening, g);
          }
      } else if (e.type == "salt") {
        for (auto& b : d.salt.bodies)
          if (b.name == e.ref) scale(b.height, g);
      }
    }
    if (e.type == "fault") active_faults.insert(e.ref);
    if (e.type == "fold") active_prims.insert(e.ref);
    if (e.type == "salt") active_bodies.insert(e.ref);
    evs.push_back(e);
  }
  d.events = evs;
  // Mechanisms with no event by t are removed.
  std::vector<desc::FaultDesc> fs;
  for (auto& f : d.faults)
    if (active_faults.count(f.name)) fs.push_back(f);
  for (auto& f : fs) {
    if (!f.abuts_fault.empty() && !active_faults.count(f.abuts_fault)) { f.abuts_fault.clear(); f.abuts_transfer = false; }
    if (!f.relay_fault.empty() && !active_faults.count(f.relay_fault)) f.relay_fault.clear();
    for (auto* tip : {&f.upper, &f.lower})
      if (tip->rule == "horizon" && gone.count(tip->horizon)) tip->rule = "free";
  }
  d.faults = fs;
  std::vector<desc::PrimitiveDesc> ps;
  for (auto& p : d.primitives)
    if (active_prims.count(p.name)) ps.push_back(p);
  d.primitives = ps;
  std::vector<desc::SaltBodyDesc> bs;
  for (auto& b : d.salt.bodies)
    if (active_bodies.count(b.name)) bs.push_back(b);
  d.salt.bodies = bs;
  // Backstripping: the youngest deposited horizon is the depositional surface.
  if (d.compaction.present() && !d.horizons.empty()) {
    const int y = d0.horizon_index(d.horizons.front().name);
    d.compaction.surface_depth = make_param(rd[static_cast<std::size_t>(y)], "m", Dimension::Length, src);
  }
  d.model.title = (d.model.title.empty() ? d.model.name : d.model.title) + " at " + fmt_shortest(t) + " Ma";
  return d;
}

ModelPtr build_state(const desc::Description& d, const Model& full, double age, BuildLog* log) {
  BuildLog local;
  BuildLog& lg = log ? *log : local;
  auto sd = state_description(d, full, age);
  if (sd.horizons.empty()) throw Error(ErrorCode::InvalidArgument, "no horizon is deposited by " + fmt_shortest(age) + " Ma");
  return build_model(sd, BuildInputs{}, lg, "parametric");
}

}  // namespace stratonix
