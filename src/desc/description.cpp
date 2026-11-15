// Typed model description: parameter parsing, schema validation, migrations and semantic rules.
#include "stratonix/desc/description.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"

namespace stratonix::desc {

const char* evidence_name(Evidence e) {
  switch (e) {
    case Evidence::Measured: return "measured";
    case Evidence::Interpreted: return "interpreted";
    case Evidence::Calibrated: return "calibrated";
    case Evidence::Upscaled: return "upscaled";
    case Evidence::LiteratureDerived: return "literature_derived";
    case Evidence::Assumed: return "assumed";
    default: return "user_instruction";
  }
}
const char* evidence_plain(Evidence e) {
  switch (e) {
    case Evidence::Measured: return "measured directly";
    case Evidence::Interpreted: return "interpreted from data";
    case Evidence::Calibrated: return "calibrated against data";
    case Evidence::Upscaled: return "upscaled from finer data";
    case Evidence::LiteratureDerived: return "taken from published literature";
    case Evidence::Assumed: return "assumed by the modeller";
    default: return "set by instruction of the model owner";
  }
}
static bool parse_evidence(const std::string& s, Evidence& e) {
  static const std::pair<const char*, Evidence> t[] = {
      {"measured", Evidence::Measured}, {"interpreted", Evidence::Interpreted},
      {"calibrated", Evidence::Calibrated}, {"upscaled", Evidence::Upscaled},
      {"literature_derived", Evidence::LiteratureDerived}, {"assumed", Evidence::Assumed},
      {"user_instruction", Evidence::UserInstruction}};
  for (auto& [n, v] : t)
    if (s == n) { e = v; return true; }
  return false;
}

std::vector<double> Param::si_list() const {
  std::vector<double> r(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) r[i] = values[i] * to_si;
  return r;
}

bool ValidationReport::ok() const {
  for (auto& i : issues)
    if (i.error) return false;
  return true;
}
std::string ValidationReport::text() const {
  std::ostringstream s;
  for (auto& i : issues) {
    s << (i.error ? "error" : "warning") << ": " << i.path;
    if (i.line > 0) s << " (line " << i.line << ")";
    s << ": " << i.message << " [rule " << i.rule << "]\n";
  }
  return s.str();
}

int Description::horizon_index(const std::string& n) const {
  for (std::size_t i = 0; i < horizons.size(); ++i)
    if (horizons[i].name == n) return static_cast<int>(i);
  return -1;
}
int Description::unit_index(const std::string& n) const {
  for (std::size_t i = 0; i < units.size(); ++i)
    if (units[i].name == n) return static_cast<int>(i);
  return -1;
}
int Description::fault_index(const std::string& n) const {
  for (std::size_t i = 0; i < faults.size(); ++i)
    if (faults[i].name == n) return static_cast<int>(i);
  return -1;
}
std::vector<const EventDesc*> Description::events_for(const std::string& type, const std::string& ref) const {
  std::vector<const EventDesc*> out;
  for (auto& e : events)
    if (e.type == type && e.ref == ref) out.push_back(&e);
  return out;
}
int Description::series_of(int h) const {
  if (h < 0 || static_cast<std::size_t>(h) >= horizons.size()) return -1;
  for (std::size_t k = 0; k < series.size(); ++k)
    for (auto& n : series[k].horizons)
      if (n == horizons[static_cast<std::size_t>(h)].name) return static_cast<int>(k);
  return -1;
}
const WellDesc* Description::well(const std::string& n) const {
  for (auto& w : wells)
    if (w.name == n) return &w;
  return nullptr;
}
Domain Description::domain() const {
  Domain d;
  d.x0 = origin.si(0); d.y0 = origin.si(1);
  d.sx = size.si(0); d.sy = size.si(1);
  d.z_top = -depth_range.si(0); d.z_base = -depth_range.si(1);
  d.spacing = grid_spacing.si();
  return d;
}
std::vector<const Param*> Description::all_params() const {
  std::vector<const Param*> r;
  auto add = [&](const Param& p) { if (p.present) r.push_back(&p); };
  for (auto* p : {&origin, &size, &depth_range, &grid_spacing, &ref_point, &regional_dip, &correction_range,
                  &regional_dip_azimuth, &frame_point, &tol_thickness_min, &tol_crossing,
                  &tol_offset_outside_tips, &tol_throw, &tol_wells, &tol_unoffset, &tol_rollover, &tol_area, &tol_salt_volume, &tol_solid_thickness, &tol_restoration, &tol_balance, &scaling_min_ratio, &scaling_max_ratio,
                  &closure_faults_seal, &target_closure, &target_area, &target_spill_depth,
                  &target_tolerance, &vertical_exaggeration, &voxel_spacing, &fault_stick_spacing, &implicit.nugget,
                  &implicit.surface_spacing, &implicit.iterations, &implicit.vertical_stretch})
    add(*p);
  for (auto& w : wells) {
    add(w.head); add(w.datum_elevation); add(w.total_depth); add(w.stations);
    for (auto& [h, p] : w.tops) add(p);
    for (auto& c : w.completions) { add(c.top); add(c.base); }
  }
  for (auto& h : horizons) { add(h.depth); add(h.thickness_above); add(h.iso_thickness); add(h.iso_gradient); add(h.iso_centre); }
  for (auto& h : horizons) add(h.age);
  for (auto& e : events)
    for (auto* q : {&e.age, &e.from, &e.to, &e.throw_}) add(*q);
  add(salt.source_thickness);
  for (auto& b : salt.bodies)
    for (auto* q : {&b.centre, &b.radius, &b.flank_width, &b.height, &b.withdrawal_width, &b.axis_azimuth, &b.half_length}) add(*q);
  add(compaction.surface_depth);
  for (auto& u : compaction.units)
    for (auto* q : {&u.surface_porosity, &u.decay, &u.table_depths, &u.table_porosity}) add(*q);
  for (auto& sr : series)
    for (auto* q : {&sr.tilt_dip, &sr.tilt_azimuth, &sr.clino_height, &sr.clino_width, &sr.clino_direction, &sr.clino_first_break,
                    &sr.clino_progradation})
      add(*q);
  for (auto& p : primitives) {
    for (auto* q : {&p.crest, &p.axis_azimuth, &p.plunge, &p.half_length, &p.half_width_left,
                    &p.half_width_right, &p.radius_major, &p.radius_minor, &p.rim_depth, &p.rim_width,
                    &p.half_width, &p.shortening, &p.array_step_along, &p.array_step_across})
      add(*q);
    for (auto& [h, q] : p.relief) add(q);
  }
  for (auto& f : faults) {
    for (auto* q : {&f.age, &f.reference_point, &f.reference_depth, &f.strike, &f.dip, &f.trace,
                    &f.trace_depth, &f.decay_length, &f.upper.offset, &f.upper.fraction, &f.upper.depth,
                    &f.lower.offset, &f.lower.fraction, &f.lower.depth, &f.length, &f.throw_max,
                    &f.profile_dip.taper, &f.profile_strike.taper, &f.influence_distance, &f.detachment_depth,
                    &f.profile, &f.heave, &f.shear_angle, &f.drag_fraction, &f.drag_width, &f.net_slip,
                    &f.rake, &f.fold_apical_angle, &f.fold_ps_ratio})
      add(*q);
  }
  for (auto& s : sections) { add(s.polyline); add(s.detail_window); }
  for (auto& sh : sheets) {
    for (auto* pd : {&sh.main, &sh.detail}) {
      add(pd->window); add(pd->vertical_exaggeration); add(pd->label_at);
      for (auto& a : pd->annotations) { add(a.at); add(a.text_at); }
    }
    for (auto& k : sh.keys) add(k.fraction);
  }
  return r;
}

// ---- parameters ----------------------------------------------------------------------------
namespace {

void flatten(const Node& n, std::vector<double>& vals, std::vector<std::size_t>& shape, int depth,
             bool& ok) {
  if (n.is_number()) { vals.push_back(n.number); return; }
  if (!n.is_seq()) { ok = false; return; }
  if (static_cast<int>(shape.size()) <= depth) shape.push_back(n.seq.size());
  else if (shape[static_cast<std::size_t>(depth)] != n.seq.size()) { ok = false; return; }
  for (auto& c : n.seq) flatten(c, vals, shape, depth + 1, ok);
}

struct Ctx {
  std::vector<Issue>* issues;
  void err(const std::string& path, const std::string& rule, const std::string& msg, int line = 0) {
    issues->push_back({path, rule, msg, true, line});
  }
  void warn(const std::string& path, const std::string& rule, const std::string& msg, int line = 0) {
    issues->push_back({path, rule, msg, false, line});
  }
};

/// Parse a parameter node in compact or mapping form. Returns false with a message on failure.
bool parse_param(const Node& n, const std::string& path, Param& p, std::string& why) {
  p = Param{};
  p.path = path;
  std::string value_text, unit, ev, src;
  Node value_node;
  if (n.is_string()) {
    auto parts = split(n.str, '|');
    if (parts.size() < 3) {
      why = parts.size() == 1 ? "no evidence class and no source" : "no source";
      return false;
    }
    std::string vu = trim(parts[0]);
    ev = trim(parts[1]);
    src = trim(parts[2]);
    for (std::size_t i = 3; i < parts.size(); ++i) src += "|" + parts[i];
    src = trim(src);
    auto sp = vu.find_last_of(" \t");
    if (sp == std::string::npos) { why = "no unit (write the value followed by its unit, e.g. '60 m')"; return false; }
    value_text = trim(vu.substr(0, sp));
    unit = trim(vu.substr(sp + 1));
    try {
      value_node = parse_document(value_text, path);
    } catch (const Error&) {
      why = "value '" + value_text + "' is not a number or list";
      return false;
    }
  } else if (n.is_map()) {
    auto* v = n.get("value");
    auto* u = n.get("unit");
    auto* e = n.get("evidence");
    auto* s = n.get("source");
    std::vector<std::string> missing;
    if (!v) missing.push_back("value");
    if (!u || u->str.empty()) missing.push_back("unit");
    if (!e || e->str.empty()) missing.push_back("evidence class");
    if (!s || s->str.empty()) missing.push_back("source");
    if (!missing.empty()) {
      why = "no ";
      for (std::size_t i = 0; i < missing.size(); ++i) why += (i ? (i + 1 == missing.size() ? " and " : ", ") : "") + missing[i];
      return false;
    }
    value_node = *v;
    unit = u->str; ev = e->str; src = s->str;
    value_text = v->is_seq() ? to_json(*v, 0) : v->str;
    if (!value_text.empty() && value_text.back() == '\n') value_text.pop_back();
  } else if (n.is_number()) {
    why = "bare number without unit, evidence class and source";
    return false;
  } else {
    why = "not a parameter";
    return false;
  }
  if (ev.empty()) { why = "no evidence class"; return false; }
  if (src.empty()) { why = "no source"; return false; }
  if (!parse_evidence(ev, p.evidence)) {
    why = "unknown evidence class '" + ev + "' (allowed: measured, interpreted, calibrated, upscaled, "
          "literature_derived, assumed, user_instruction)";
    return false;
  }
  auto ud = find_unit(unit);
  if (!ud) { why = "unknown unit '" + unit + "'"; return false; }
  p.unit = unit;
  p.dim = ud->dim;
  p.to_si = ud->to_si;
  p.source = src;
  p.text = value_text;
  if (ud->dim == Dimension::Flag) {
    if (value_node.type != Node::Type::Bool) { why = "a flag must be true or false"; return false; }
    p.flag = value_node.boolean;
    p.values = {p.flag ? 1.0 : 0.0};
  } else {
    bool ok = true;
    flatten(value_node, p.values, p.shape, 0, ok);
    if (!ok || p.values.empty()) { why = "value '" + value_text + "' is not a number or a regular list of numbers"; return false; }
    for (double x : p.values)
      if (!std::isfinite(x)) { why = "value is not finite"; return false; }
  }
  p.present = true;
  return true;
}

struct Typer {
  Ctx c;
  const Node& root;
  Description& d;

  Param P(const Node* parent, const std::string& key, const std::string& base, Dimension dim,
          bool required, std::size_t n_values = 0, std::size_t tuple = 0) {
    Param p;
    std::string path = base.empty() ? key : base + "." + key;
    const Node* n = parent ? parent->get(key) : nullptr;
    if (!n) {
      if (required) c.err(path, "S0", "required parameter is missing");
      p.path = path;
      return p;
    }
    std::string why;
    if (!parse_param(*n, path, p, why)) {
      c.err(path, "E1", why + "; every parameter carries value, unit, evidence class and source", n->line);
      p.present = false;
      return p;
    }
    if (dim != Dimension::Any && p.dim != dim)
      c.err(path, "U1", std::string("unit '") + p.unit + "' is a " + dimension_name(p.dim) + " unit; this field needs a " +
                            dimension_name(dim) + " unit", n->line);
    if (n_values && p.values.size() != n_values && !tuple)
      c.err(path, "S0", "expected " + std::to_string(n_values) + " value(s), found " + std::to_string(p.values.size()), n->line);
    if (tuple && (p.shape.size() != 2 || p.shape[1] != tuple || p.shape[0] < n_values))
      c.err(path, "S0", "expected a list of at least " + std::to_string(n_values) + " points of " +
                            std::to_string(tuple) + " values each", n->line);
    return p;
  }
  static std::string S(const Node* parent, const std::string& key, const std::string& def = {}) {
    const Node* n = parent ? parent->get(key) : nullptr;
    return n && !n->is_null() ? n->str : def;
  }

  void run() {
    const Node* m = root.get("model");
    d.model = {S(m, "name"), S(m, "version"), S(m, "status"), S(m, "title"), S(m, "description"), S(m, "contact")};
    const Node* crs = root.get("crs");
    d.crs.name = S(crs, "name");
    if (crs && crs->get("epsg")) d.crs.epsg = static_cast<int>(crs->get("epsg")->number);
    d.crs.wkt = S(crs, "wkt");
    d.crs.vertical_datum = S(crs, "vertical_datum");
    d.crs.length_unit = S(crs, "length_unit", "m");

    const Node* dom = root.get("domain");
    d.origin = P(dom, "origin", "domain", Dimension::Length, true, 2);
    d.size = P(dom, "size", "domain", Dimension::Length, true, 2);
    d.depth_range = P(dom, "depth_range", "domain", Dimension::Length, true, 2);
    d.grid_spacing = P(dom, "grid_spacing", "domain", Dimension::Length, true, 1);

    if (const Node* in = root.get("inputs"))
      for (auto& i : in->seq)
        d.inputs.push_back({S(&i, "id"), S(&i, "path"), S(&i, "role"), S(&i, "sharing"), S(&i, "source"), S(&i, "evidence"), {}});

    if (const Node* ws = root.get("wells"))
      for (std::size_t k = 0; k < ws->seq.size(); ++k) {
        const Node& w = ws->seq[k];
        std::string b = "wells[" + std::to_string(k) + "]";
        WellDesc wd;
        wd.name = S(&w, "name");
        wd.role = S(&w, "role");
        wd.head = P(&w, "head", b, Dimension::Length, true, 2);
        wd.datum_elevation = P(&w, "datum_elevation", b, Dimension::Length, false, 1);
        wd.total_depth = P(&w, "total_depth", b, Dimension::Length, false, 1);
        if (const Node* tr = w.get("trajectory")) wd.stations = P(tr, "stations", b + ".trajectory", Dimension::Any, true, 2, 3);
        if (const Node* cs = w.get("completions"))
          for (std::size_t q = 0; q < cs->seq.size(); ++q) {
            std::string bc = b + ".completions[" + std::to_string(q) + "]";
            wd.completions.push_back({S(&cs->seq[q], "name"), P(&cs->seq[q], "top", bc, Dimension::Length, true, 1),
                                      P(&cs->seq[q], "base", bc, Dimension::Length, true, 1)});
          }
        if (const Node* tops = w.get("tops"))
          for (auto& [h, tn] : tops->map) wd.tops.emplace_back(h, P(tops, h, b + ".tops", Dimension::Length, true, 1));
        d.wells.push_back(std::move(wd));
      }

    const Node* st = root.get("stratigraphy");
    if (const Node* ref = st ? st->get("reference") : nullptr) {
      d.ref_well = S(ref, "well");
      d.ref_point = P(ref, "point", "stratigraphy.reference", Dimension::Length, false, 2);
    }
    if (const Node* hs = st ? st->get("horizons") : nullptr)
      for (std::size_t k = 0; k < hs->seq.size(); ++k) {
        const Node& h = hs->seq[k];
        std::string b = "stratigraphy.horizons[" + std::to_string(k) + "]";
        HorizonDesc hd;
        hd.name = S(&h, "name");
        hd.depth = P(&h, "depth", b, Dimension::Length, false, 1);
        hd.thickness_above = P(&h, "thickness_above", b, Dimension::Length, false, 1);
        hd.age = P(&h, "age", b, Dimension::Time, false, 1);
        if (const Node* s = h.get("surface")) {
          hd.surface_input = S(s, "input");
          hd.surface_format = S(s, "format");
          hd.surface_z_is = S(s, "z_is", "depth");
        }
        if (const Node* s = h.get("isochore")) {
          hd.iso_from = S(s, "from");
          hd.iso_input = S(s, "input");
          hd.iso_format = S(s, "format");
          hd.iso_thickness = P(s, "thickness", b + ".isochore", Dimension::Length, false, 1);
          if (const Node* tr = s->get("trend")) {
            hd.iso_trend = S(tr, "type");
            const std::string tb = b + ".isochore.trend";
            hd.iso_gradient = P(tr, "gradient", tb, Dimension::Dimensionless, true, hd.iso_trend == "radial" ? 1 : 2);
            hd.iso_centre = P(tr, "centre", tb, Dimension::Length, hd.iso_trend == "radial", 2);
          }
        }
        if (const Node* io = h.get("intended_offset")) hd.intended_offset = io->boolean;
        d.horizons.push_back(std::move(hd));
      }
    if (const Node* us = st ? st->get("units") : nullptr)
      for (auto& u : us->seq)
        d.units.push_back({S(&u, "name"), S(&u, "top"), S(&u, "base"), S(&u, "role"), S(&u, "color"), S(&u, "label")});
    if (const Node* ss = st ? st->get("series") : nullptr)
      for (std::size_t k = 0; k < ss->seq.size(); ++k) {
        const Node& sn = ss->seq[k];
        std::string b = "stratigraphy.series[" + std::to_string(k) + "]";
        SeriesDesc sd;
        sd.name = S(&sn, "name");
        sd.base = S(&sn, "base");
        if (const Node* hs2 = sn.get("horizons"))
          for (auto& x : hs2->seq) sd.horizons.push_back(x.str);
        if (const Node* t = sn.get("tilt")) {
          sd.tilt_dip = P(t, "dip", b + ".tilt", Dimension::Angle, true, 1);
          sd.tilt_azimuth = P(t, "azimuth", b + ".tilt", Dimension::Angle, true, 1);
        }
        if (const Node* cf = sn.get("clinoforms")) {
          const std::string cb = b + ".clinoforms";
          sd.clino_profile = S(cf, "profile");
          sd.clino_height = P(cf, "height", cb, Dimension::Length, true, 1);
          sd.clino_width = P(cf, "width", cb, Dimension::Length, true, 1);
          sd.clino_direction = P(cf, "direction", cb, Dimension::Angle, true, 1);
          sd.clino_first_break = P(cf, "first_break", cb, Dimension::Length, true, 2);
          sd.clino_progradation = P(cf, "progradation", cb, Dimension::Length, true, 1);
        }
        d.series.push_back(std::move(sd));
      }
    if (const Node* cp = st ? st->get("compaction") : nullptr) {
      const std::string cb = "stratigraphy.compaction";
      d.compaction.geometry = S(cp, "geometry");
      d.compaction.surface_depth = P(cp, "surface_depth", cb, Dimension::Length, true, 1);
      if (const Node* us = cp->get("units"))
        for (std::size_t k = 0; k < us->seq.size(); ++k) {
          const Node& un = us->seq[k];
          const std::string ub = cb + ".units[" + std::to_string(k) + "]";
          CompactionUnitDesc cu;
          cu.unit = S(&un, "unit");
          cu.law = S(&un, "law");
          const bool ex = cu.law == "exponential";
          cu.surface_porosity = P(&un, "surface_porosity", ub, Dimension::Dimensionless, ex, 1);
          cu.decay = P(&un, "decay", ub, Dimension::InverseLength, ex, 1);
          cu.table_depths = P(&un, "table_depths", ub, Dimension::Length, !ex, 0);
          cu.table_porosity = P(&un, "table_porosity", ub, Dimension::Dimensionless, !ex, 0);
          d.compaction.units.push_back(std::move(cu));
        }
    }
    d.above_top_name = S(st, "above_top_name");
    d.below_base_name = S(st, "below_base_name");
    for (auto& u : d.units) {
      if (u.base == "model_base" && u.top != "model_top" && d.below_base_name.empty()) d.below_base_name = u.name;
      if (u.top == "model_top" && u.base != "model_base" && d.above_top_name.empty()) d.above_top_name = u.name;
    }

    if (const Node* s = root.get("structure")) {
      if (const Node* r = s->get("regional")) {
        d.regional_type = S(r, "type", "flat");
        d.regional_dip = P(r, "dip", "structure.regional", Dimension::Angle, d.regional_type == "planar", 1);
        d.regional_dip_azimuth = P(r, "dip_azimuth", "structure.regional", Dimension::Angle, d.regional_type == "planar", 1);
      }
      if (const Node* bs = s->get("base")) d.base_type = S(bs, "type", "conformable");
      if (const Node* sl = s->get("salt")) {
        const std::string sb = "structure.salt";
        d.salt.top = S(sl, "top");
        d.salt.base = S(sl, "base");
        d.salt.source_thickness = P(sl, "source_thickness", sb, Dimension::Length, true, 1);
        if (const Node* bd = sl->get("bodies"))
          for (std::size_t k = 0; k < bd->seq.size(); ++k) {
            const Node& bn = bd->seq[k];
            const std::string bb = sb + ".bodies[" + std::to_string(k) + "]";
            SaltBodyDesc b;
            b.type = S(&bn, "type");
            b.name = S(&bn, "name");
            const bool steep = b.type == "diapir" || b.type == "wall", wall = b.type == "wall";
            b.centre = P(&bn, "centre", bb, Dimension::Length, true, 2);
            b.radius = P(&bn, "radius", bb, Dimension::Length, true, 1);
            b.flank_width = P(&bn, "flank_width", bb, Dimension::Length, steep, 1);
            b.height = P(&bn, "height", bb, Dimension::Length, true, 1);
            b.withdrawal_width = P(&bn, "withdrawal_width", bb, Dimension::Length, true, 1);
            b.axis_azimuth = P(&bn, "axis_azimuth", bb, Dimension::Angle, wall, 1);
            b.half_length = P(&bn, "half_length", bb, Dimension::Length, wall, 1);
            if (const Node* g = bn.get("growth")) b.growth = {S(g, "onset"), S(g, "end")};
            d.salt.bodies.push_back(std::move(b));
          }
      }
      d.correction_range = P(s, "correction_range", "structure", Dimension::Length, false, 1);
      if (const Node* ps = s->get("primitives"))
        for (std::size_t k = 0; k < ps->seq.size(); ++k) {
          const Node& p = ps->seq[k];
          std::string b = "structure.primitives[" + std::to_string(k) + "]";
          PrimitiveDesc pd;
          pd.type = S(&p, "type");
          pd.name = S(&p, "name");
          bool det = pd.type == "detachment_fold";
          bool fold = pd.type == "anticline" || pd.type == "syncline";
          pd.crest = P(&p, "crest", b, Dimension::Length, true, 2);
          pd.axis_azimuth = P(&p, "axis_azimuth", b, Dimension::Angle, fold || det, 1);
          pd.half_width = P(&p, "half_width", b, Dimension::Length, det, 1);
          pd.shortening = P(&p, "shortening", b, Dimension::Length, det, 1);
          pd.detachment = S(&p, "detachment");
          pd.mechanics = S(&p, "mechanics");
          pd.series = S(&p, "series");
          if (const Node* g = p.get("growth")) pd.growth = {S(g, "onset"), S(g, "end")};
          if (const Node* ar = p.get("array")) {
            std::string c0 = S(ar, "count", "1");
            pd.array_count = std::max(1, std::atoi(c0.c_str()));
            pd.array_step_along = P(ar, "step_along", b + ".array", Dimension::Length, true, 1);
            pd.array_step_across = P(ar, "step_across", b + ".array", Dimension::Length, true, 1);
          }
          if (det) {
            if (p.get("relief"))
              c.err(b + ".relief", "F19", "a detachment fold takes no relief: it follows from shortening x depth to the detachment");
            if (pd.detachment.empty()) c.err(b + ".detachment", "F19", "a detachment fold needs detachment: <horizon>");
            if (p.get("half_width_left") || p.get("half_width_right"))
              c.err(b + ".half_width", "F19", "a detachment fold is symmetric: give half_width");
          } else {
            if (!p.get("relief")) c.err(b + ".relief", "F7", "relief must be given (no hidden defaults)");
            for (const char* k : {"half_width", "shortening", "detachment"})
              if (p.get(k)) c.err(b + "." + k, "F19", std::string(k) + " belongs to detachment folds");
          }
          pd.plunge = P(&p, "plunge", b, Dimension::Angle, false, 1);
          pd.half_length = P(&p, "half_length", b, Dimension::Length, false, 1);
          pd.half_width_left = P(&p, "half_width_left", b, Dimension::Length, fold, 1);
          pd.half_width_right = P(&p, "half_width_right", b, Dimension::Length, fold, 1);
          pd.radius_major = P(&p, "radius_major", b, Dimension::Length, !fold && !det, 1);
          pd.radius_minor = P(&p, "radius_minor", b, Dimension::Length, false, 1);
          pd.rim_depth = P(&p, "rim_depth", b, Dimension::Length, false, 1);
          pd.rim_width = P(&p, "rim_width", b, Dimension::Length, pd.rim_depth.present, 1);
          if (const Node* rel = p.get("relief"))
            for (auto& [h, rn] : rel->map) pd.relief.emplace_back(h, P(rel, h, b + ".relief", Dimension::Length, true, 1));
          d.primitives.push_back(std::move(pd));
        }
    }

    if (const Node* fs = root.get("faults"))
      for (std::size_t k = 0; k < fs->seq.size(); ++k) {
        const Node& f = fs->seq[k];
        std::string b = "faults[" + std::to_string(k) + "]";
        FaultDesc fd;
        fd.name = S(&f, "name");
        fd.age = P(&f, "age", b, Dimension::Ordinal, true, 1);
        const Node* g = f.get("geometry");
        fd.geom_type = S(g, "type", "planar");
        bool by_trace = g && g->get("trace");
        bool by_sticks = fd.geom_type == "sticks";
        const bool curved = fd.curved(), ramp = fd.geom_type == "ramp_flat";
        fd.geom_input = S(g, "input");
        fd.geom_input_format = S(g, "input_format", "irap_lines");
        fd.reference_point = P(g, "reference_point", b + ".geometry", Dimension::Length, !by_trace && !by_sticks, 2);
        fd.reference_depth = P(g, "reference_depth", b + ".geometry", Dimension::Length, !by_trace && !by_sticks && !ramp, 1);
        fd.strike = P(g, "strike", b + ".geometry", Dimension::Angle, !by_trace && !by_sticks, 1);
        fd.dip = P(g, "dip", b + ".geometry", Dimension::Angle, !by_sticks && !ramp, 1);
        fd.detachment_depth = P(g, "detachment_depth", b + ".geometry", Dimension::Length, fd.geom_type == "listric", 1);
        fd.profile = P(g, "profile", b + ".geometry", Dimension::Length, ramp, 2, 2);
        fd.trace = P(g, "trace", b + ".geometry", Dimension::Length, false, 2, 2);
        fd.trace_depth = P(g, "trace_depth", b + ".geometry", Dimension::Length, by_trace, 1);
        fd.dip_side = S(g, "dip_side");
        fd.decay_length = P(g, "decay_length", b + ".geometry", Dimension::Length, false, 1);
        const Node* e = f.get("extent");
        auto tip = [&](const char* key, TipDesc& t) {
          const Node* tn = e ? e->get(key) : nullptr;
          std::string tb = b + ".extent." + key;
          t.rule = S(tn, "rule");
          t.horizon = S(tn, "horizon");
          t.unit = S(tn, "unit");
          t.offset = P(tn, "offset", tb, Dimension::Length, t.rule == "horizon", 1);
          t.fraction = P(tn, "fraction", tb, Dimension::Dimensionless, t.rule == "within_unit", 1);
          t.depth = P(tn, "depth", tb, Dimension::Length, t.rule == "depth", 1);
          t.stratigraphy = S(tn, "stratigraphy", "fault_frame");
        };
        tip("upper_tip", fd.upper);
        tip("lower_tip", fd.lower);
        fd.length = P(e, "length", b + ".extent", Dimension::Length, true, 1);
        fd.tip_line_shape = S(e, "tip_line_shape", "rectangular");
        const Node* th = f.get("throw");
        fd.slip_mode = S(&f, "slip_mode");
        const bool fpf_hw = f.get("hanging_wall") && S(f.get("hanging_wall"), "type") == "fault_parallel_flow";
        const bool by_slip = fd.strike_component() || (curved && fpf_hw);
        fd.throw_max = P(th, "max", b + ".throw", Dimension::Length, !curved && !by_slip, 1);
        fd.net_slip = P(th, "slip", b + ".throw", Dimension::Length, by_slip, 1);
        fd.sense = S(th, "sense");
        fd.rake = P(&f, "rake", b, Dimension::Angle, fd.slip_mode == "oblique", 1);
        fd.heave = P(th, "heave", b + ".throw", Dimension::Length, curved && !fpf_hw, 1);
        auto prof = [&](const char* key, ProfileDesc& pd) {
          const Node* pn = th ? th->get(key) : nullptr;
          pd.type = S(pn, "type", "flat_topped");
          pd.taper = P(pn, "taper", b + ".throw." + key, Dimension::Dimensionless, pn && pd.type == "flat_topped", 1);
        };
        prof("profile_dip", fd.profile_dip);
        prof("profile_strike", fd.profile_strike);
        if (const Node* g = f.get("growth")) fd.growth = {S(g, "onset"), S(g, "end")};
        const Node* fo = f.get("fold");
        if (fo) {
          fd.fold_type = S(fo, "type");
          fd.fold_apical_angle = P(fo, "apical_angle", b + ".fold", Dimension::Angle, true, 1);
          fd.fold_ps_ratio = P(fo, "propagation_slip_ratio", b + ".fold", Dimension::Dimensionless, true, 1);
        }
        if (fo && fd.fold_type == "trishear") {
          if (!th || !th->get("profile_strike"))
            c.err(b + ".throw", "F7", "profile_strike must be given (no hidden defaults)");
          if (th && th->get("profile_dip"))
            c.err(b + ".throw.profile_dip", "F17", "a trishear fault takes no down-dip profile: the slip is constant behind the "
                                                   "propagating tip");
        } else if (curved) {
          if (!th || !th->get("profile_strike"))
            c.err(b + ".throw", "F7", "profile_strike must be given (no hidden defaults)");
          if (th && th->get("profile_dip"))
            c.err(b + ".throw.profile_dip", "F14", "a curved fault takes no down-dip profile: the shear construction needs a "
                                                   "constant displacement down the fault (its throw follows from the shape)");
          if (fpf_hw && th && th->get("heave"))
            c.err(b + ".throw.heave", "F18", "fault-parallel flow takes the slip along the fault as throw.slip, not heave");
          if (th && th->get("max"))
            c.err(b + ".throw.max", "F14", std::string(fpf_hw ? "fault-parallel flow takes the slip along the fault as throw.slip; " :
                                           "a curved fault takes its displacement as throw.heave, constant down the fault; ") +
                                           "its throw at each horizon follows from the shape and is reported");
        } else if (!th || !th->get("profile_dip") || !th->get("profile_strike")) {
          c.err(b + ".throw", "F7", "both profile_dip and profile_strike must be given (no hidden defaults)");
        }
        if (!curved && th && th->get("heave"))
          c.err(b + ".throw.heave", "F14", "heave is the displacement of a curved (listric or ramp_flat) fault; a planar fault takes throw.max");
        fd.slip_mode = S(&f, "slip_mode");
        const Node* hw = f.get("hanging_wall");
        fd.hw_type = S(hw, "type");
        fd.influence_distance = P(hw, "influence_distance", b + ".hanging_wall", Dimension::Length, fd.hw_type == "decay", 1);
        fd.shear_angle = P(hw, "shear_angle", b + ".hanging_wall", Dimension::Angle, fd.hw_type == "inclined_shear", 1);
        fd.shear_sense = S(hw, "shear_sense");
        if (const Node* ab = f.get("abuts")) {
          fd.abuts_fault = S(ab, "fault");
          fd.abuts_side = S(ab, "side");
          if (const Node* t = ab->get("transfer")) fd.abuts_transfer = t->type == Node::Type::Bool && t->boolean;
        }
        if (const Node* dg = f.get("drag")) {
          fd.drag_type = S(dg, "type");
          fd.drag_fraction = P(dg, "fraction", b + ".drag", Dimension::Dimensionless, true, 1);
          fd.drag_width = P(dg, "width", b + ".drag", Dimension::Length, true, 1);
        }
        if (const Node* rl = f.get("relay")) {
          fd.relay_fault = S(rl, "fault");
          if (const Node* t = rl->get("breached")) fd.relay_breached = t->type == Node::Type::Bool && t->boolean;
        }
        d.faults.push_back(std::move(fd));
      }

    if (const Node* fr = root.get("frame"))
      if (const Node* fx = fr->get("fixed")) {
        d.frame_well = S(fx, "well");
        d.frame_point = P(fx, "point", "frame.fixed", Dimension::Length, false, 3);
      }

    if (const Node* ch = root.get("checks")) {
      if (const Node* t = ch->get("tolerances")) {
        std::string b = "checks.tolerances";
        d.tol_thickness_min = P(t, "thickness_min", b, Dimension::Length, false, 1);
        d.tol_crossing = P(t, "crossing", b, Dimension::Length, false, 1);
        d.tol_offset_outside_tips = P(t, "offset_outside_tips", b, Dimension::Length, false, 1);
        d.tol_throw = P(t, "throw", b, Dimension::Length, false, 1);
        d.tol_wells = P(t, "wells", b, Dimension::Length, false, 1);
        d.tol_unoffset = P(t, "unoffset", b, Dimension::Length, false, 1);
        d.tol_rollover = P(t, "rollover", b, Dimension::Length, false, 1);
        d.tol_area = P(t, "area", b, Dimension::Dimensionless, false, 1);
        d.tol_salt_volume = P(t, "salt_volume", b, Dimension::Dimensionless, false, 1);
        d.tol_solid_thickness = P(t, "solid_thickness", b, Dimension::Dimensionless, false, 1);
        d.tol_restoration = P(t, "restoration", b, Dimension::Length, false, 1);
        d.tol_balance = P(t, "balance", b, Dimension::Dimensionless, false, 1);
      }
      if (const Node* fsc = ch->get("fault_scaling")) {
        d.scaling_min_ratio = P(fsc, "min_ratio", "checks.fault_scaling", Dimension::Dimensionless, false, 1);
        d.scaling_max_ratio = P(fsc, "max_ratio", "checks.fault_scaling", Dimension::Dimensionless, false, 1);
      }
      if (const Node* cl = ch->get("closure")) {
        d.closure_horizon = S(cl, "horizon");
        d.closure_faults_seal = P(cl, "faults_seal", "checks.closure", Dimension::Flag, true, 1);
        if (const Node* tg = cl->get("targets")) {
          std::string b = "checks.closure.targets";
          d.target_closure = P(tg, "closure", b, Dimension::Length, false, 1);
          d.target_area = P(tg, "area", b, Dimension::Area, false, 1);
          d.target_spill_depth = P(tg, "spill_depth", b, Dimension::Length, false, 1);
          d.target_tolerance = P(tg, "tolerance_fraction", b, Dimension::Dimensionless, false, 1);
        }
      }
      if (const Node* wo = ch->get("warn_only"))
        for (auto& x : wo->seq) d.warn_only.push_back(x.str);
    }

    if (const Node* im = root.get("implicit")) {
      d.implicit.present = true;
      if (const Node* data = im->get("data"))
        for (auto& x : data->seq) d.implicit.data.push_back({S(&x, "input"), S(&x, "kind"), S(&x, "format"), S(&x, "horizon")});
      if (const Node* u = im->get("use_well_tops")) d.implicit.use_well_tops = u->boolean || u->str == "true";
      d.implicit.kernel = S(im, "kernel", "cubic");
      d.implicit.nugget = P(im, "nugget", "implicit", Dimension::Length, false, 1);
      d.implicit.surface_spacing = P(im, "surface_spacing", "implicit", Dimension::Length, false, 1);
      d.implicit.iterations = P(im, "iterations", "implicit", Dimension::Ordinal, false, 1);
      d.implicit.vertical_stretch = P(im, "vertical_stretch", "implicit", Dimension::Dimensionless, false, 1);
    }
    if (const Node* o = root.get("outputs")) {
      if (const Node* ss = o->get("sections"))
        for (std::size_t k = 0; k < ss->seq.size(); ++k) {
          std::string b = "outputs.sections[" + std::to_string(k) + "]";
          SectionDesc sd;
          sd.name = S(&ss->seq[k], "name");
          sd.polyline = P(&ss->seq[k], "polyline", b, Dimension::Length, true, 2, 2);
          sd.detail_window = P(&ss->seq[k], "detail_window", b, Dimension::Length, false, 4);
          d.sections.push_back(std::move(sd));
        }
      if (const Node* sh = o->get("section_sheets"))
        for (std::size_t k = 0; k < sh->seq.size(); ++k) {
          const Node& n = sh->seq[k];
          std::string b = "outputs.section_sheets[" + std::to_string(k) + "]";
          SheetDesc sd;
          sd.name = S(&n, "name");
          sd.section = S(&n, "section");
          sd.title = S(&n, "title");
          if (const Node* al = n.get("axis_labels")) { sd.along_label = S(al, "along"); sd.depth_label = S(al, "depth"); }
          if (const Node* el = n.get("end_labels"); el && el->seq.size() == 2) { sd.end_left = el->seq[0].str; sd.end_right = el->seq[1].str; }
          auto panel = [&](const char* key, SheetPanelDesc& pd) {
            const Node* pn = n.get(key);
            if (!pn) return;
            std::string bp = b + "." + key;
            pd.present = true;
            pd.title = S(pn, "title");
            pd.window = P(pn, "window", bp, Dimension::Length, false, 4);
            pd.vertical_exaggeration = P(pn, "vertical_exaggeration", bp, Dimension::Dimensionless, false, 1);
            if (const Node* ul = pn->get("unit_labels")) {
              pd.label_at = P(ul, "at", bp + ".unit_labels", Dimension::Length, true, 1);
              if (const Node* us = ul->get("units"))
                for (auto& u : us->seq) pd.label_units.push_back(u.str);
            }
            if (const Node* an = pn->get("annotations"))
              for (std::size_t q = 0; q < an->seq.size(); ++q) {
                const Node& a = an->seq[q];
                std::string ba = bp + ".annotations[" + std::to_string(q) + "]";
                SheetAnnotationDesc ad;
                ad.type = S(&a, "type");
                ad.text = S(&a, "text");
                ad.horizon = S(&a, "horizon");
                ad.fault = S(&a, "fault");
                ad.side = S(&a, "side", "footwall");
                ad.well = S(&a, "well");
                ad.completion = S(&a, "completion");
                ad.at = P(&a, "at", ba, Dimension::Length, ad.type == "leader", 1);
                ad.text_at = P(&a, "text_at", ba, Dimension::Length, true, 2);
                pd.annotations.push_back(std::move(ad));
              }
          };
          panel("main", sd.main);
          panel("detail", sd.detail);
          if (const Node* lg = n.get("legend")) {
            sd.well_legend = S(lg, "well");
            sd.completion_legend = S(lg, "completion");
            if (const Node* fl = lg->get("faults"))
              for (auto& x : fl->seq) sd.fault_legend.push_back(x.str);
          }
          sd.key_title = S(&n, "key_numbers_title", "Key numbers");
          if (const Node* kn = n.get("key_numbers"))
            for (std::size_t q = 0; q < kn->seq.size(); ++q) {
              const Node& a = kn->seq[q];
              SheetKeyDesc kd;
              kd.kind = S(&a, "kind");
              kd.label = S(&a, "label");
              kd.horizon = S(&a, "horizon");
              kd.unit = S(&a, "unit");
              kd.well = S(&a, "well");
              kd.fault = S(&a, "fault");
              kd.end = S(&a, "end", "last");
              if (const Node* fs = a.get("faults"))
                for (auto& x : fs->seq) kd.faults.push_back(x.str);
              kd.fraction = P(&a, "fraction", b + ".key_numbers[" + std::to_string(q) + "]", Dimension::Dimensionless, false, 1);
              if (const Node* dc = a.get("decimals")) kd.decimals = static_cast<int>(dc->number);
              kd.basis = S(&a, "basis");
              kd.value = S(&a, "value");
              kd.completion = S(&a, "completion");
              kd.side = S(&a, "side", "footwall");
              kd.coordinate = S(&a, "coordinate", "distance");
              kd.from = S(&a, "from");
              kd.to = S(&a, "to");
              if (const Node* os = a.get("on_sheet")) kd.on_sheet = os->type == Node::Type::Bool ? os->boolean : os->str != "false";
              sd.keys.push_back(std::move(kd));
            }
          if (const Node* nt = n.get("notes"))
            for (auto& x : nt->seq) sd.notes.push_back(x.str);
          d.sheets.push_back(std::move(sd));
        }
      d.vertical_exaggeration = P(o, "vertical_exaggeration", "outputs", Dimension::Dimensionless, false, 1);
      d.voxel_spacing = P(o, "voxel_spacing", "outputs", Dimension::Length, false, 3);
      d.fault_stick_spacing = P(o, "fault_stick_spacing", "outputs", Dimension::Length, false, 1);
      if (const Node* ex = o->get("exports"))
        for (auto& x : ex->seq) d.exports.push_back(x.str);
      if (const Node* ip = o->get("isopachs"))
        for (auto& x : ip->seq) d.isopachs.push_back(x.str);
    }
    if (const Node* ev = root.get("events"))
      for (std::size_t k = 0; k < ev->seq.size(); ++k) {
        const Node& en = ev->seq[k];
        const std::string eb = "events[" + std::to_string(k) + "]";
        EventDesc e;
        e.name = S(&en, "name");
        e.type = S(&en, "type");
        e.ref = e.type == "deposition" ? S(&en, "horizon") : e.type == "fault" ? S(&en, "fault") : e.type == "fold" ? S(&en, "primitive")
                                                                                                  : S(&en, "body");
        const bool span = en.get("from") || en.get("to");
        e.age = P(&en, "age", eb, Dimension::Time, !span, 1);
        e.from = P(&en, "from", eb, Dimension::Time, span, 1);
        e.to = P(&en, "to", eb, Dimension::Time, span, 1);
        e.throw_ = P(&en, "throw", eb, Dimension::Length, false, 1);
        e.sense = S(&en, "sense");
        if (e.type == "deposition" && e.age.present) {
          int h = d.horizon_index(e.ref);
          if (h >= 0 && !d.horizons[static_cast<std::size_t>(h)].age.present) d.horizons[static_cast<std::size_t>(h)].age = e.age;
        }
        d.events.push_back(std::move(e));
      }
  }

  // ---- semantic rules ----------------------------------------------------------------------
  /// Depth (positive down, metres) of each horizon at the reference location, NaN if not known
  /// before reading input files.
  std::vector<double> reference_depths() {
    std::vector<double> z(d.horizons.size(), kNaN);
    const WellDesc* rw = d.ref_well.empty() ? nullptr : d.well(d.ref_well);
    for (std::size_t i = 0; i < d.horizons.size(); ++i) {
      auto& h = d.horizons[i];
      if (rw)
        for (auto& [n, p] : rw->tops)
          if (n == h.name && p.present) z[i] = p.si();
      if (std::isnan(z[i]) && h.depth.present) z[i] = h.depth.si();
      if (std::isnan(z[i]) && h.thickness_above.present && i > 0 && !std::isnan(z[i - 1]))
        z[i] = z[i - 1] + h.thickness_above.si();
      if (std::isnan(z[i]) && !h.iso_from.empty() && h.iso_thickness.present) {
        int j = d.horizon_index(h.iso_from);
        if (j >= 0 && j < static_cast<int>(i) && !std::isnan(z[static_cast<std::size_t>(j)]))
          z[i] = z[static_cast<std::size_t>(j)] + h.iso_thickness.si();
      }
    }
    return z;
  }

  /// STX-TN-07: series (F20), trends (F21), growth (F22), clinoforms (F23).
  void series_rules() {
    const std::size_t nh = d.horizons.size();
    if (!d.series.empty()) {
      std::size_t next = 0;
      std::set<std::string> names;
      for (std::size_t k = 0; k < d.series.size(); ++k) {
        const auto& sr = d.series[k];
        const std::string b = "stratigraphy.series[" + std::to_string(k) + "]";
        if (!names.insert(sr.name).second) c.err(b + ".name", "F20", "duplicate series name '" + sr.name + "'");
        for (auto& h : sr.horizons) {
          if (next < nh && d.horizons[next].name == h) { ++next; continue; }
          c.err(b + ".horizons", "F20", "the series must list every horizon once, youngest first, in the order of stratigraphy.horizons: "
                                        "expected '" + (next < nh ? d.horizons[next].name : std::string("(none)")) + "', found '" + h + "'");
          next = nh + 1;
          break;
        }
        if (sr.base != "conformable" && sr.base != "erosional" && sr.base != "onlap")
          c.err(b + ".base", "F20", "base must be conformable, erosional or onlap");
        if (k + 1 == d.series.size() && sr.base != "conformable" && !sr.base.empty())
          c.err(b + ".base", "F20", "the oldest series rests on nothing: its base is conformable");
        if (sr.tilt_dip.present && !(sr.tilt_dip.si() >= 0 && sr.tilt_dip.si() < kPi / 4))
          c.err(b + ".tilt.dip", "F20", "tilt dip must be in [0, 45) degrees");
        if (!sr.clino_profile.empty()) {
          if (k + 1 == d.series.size())
            c.err(b + ".clinoforms", "F23", "clinoforms rest on the top of the series below; '" + sr.name + "' is the oldest series");
          if (sr.clino_height.present && !(sr.clino_height.si() > 0)) c.err(b + ".clinoforms.height", "F23", "height must be positive");
          if (sr.clino_width.present && !(sr.clino_width.si() > 0)) c.err(b + ".clinoforms.width", "F23", "width must be positive");
          if (sr.clino_progradation.present && sr.clino_progradation.si() < 0)
            c.err(b + ".clinoforms.progradation", "F23", "progradation must not be negative");
          if (sr.base == "onlap") c.err(b + ".base", "F23", "clinoforms downlap onto the series below: give base conformable or erosional");
        }
      }
      if (next == nh + 1) {
      } else if (next != nh)
        c.err("stratigraphy.series", "F20", "the series must list every horizon once; '" + d.horizons[next].name + "' is in none");
    }
    for (std::size_t i = 0; i < nh; ++i) {
      const auto& h = d.horizons[i];
      if (h.iso_trend.empty()) continue;
      const std::string b = "stratigraphy.horizons[" + std::to_string(i) + "].isochore.trend";
      if (h.iso_trend != "linear" && h.iso_trend != "radial") c.err(b + ".type", "F21", "trend type must be linear or radial");
      if (!h.iso_thickness.present) c.err(b, "F21", "a thickness trend needs the isochore's thickness at its origin");
      for (std::size_t q = 0; h.iso_gradient.present && q < h.iso_gradient.values.size(); ++q)
        if (!std::isfinite(h.iso_gradient.si(q))) c.err(b + ".gradient", "F21", "gradient must be finite");
      int j = d.horizon_index(h.iso_from);
      if (j >= 0 && d.series_of(j) != d.series_of(static_cast<int>(i)))
        c.err(b, "F21", "the isochore hangs from a horizon of another series");
    }
    auto growth = [&](const GrowthDesc& g, const std::string& b) {
      if (!g.present()) return;
      int on = d.horizon_index(g.onset), en = d.horizon_index(g.end);
      if (on < 0) c.err(b + ".growth.onset", "F22", "unknown horizon '" + g.onset + "'");
      if (en < 0) c.err(b + ".growth.end", "F22", "unknown horizon '" + g.end + "'");
      if (on >= 0 && en >= 0 && !(on > en))
        c.err(b + ".growth", "F22", "the onset horizon (pre-growth top) must be older than the end horizon (post-growth base)");
    };
    for (std::size_t k = 0; k < d.primitives.size(); ++k) {
      const auto& pd = d.primitives[k];
      const std::string b = "structure.primitives[" + std::to_string(k) + "]";
      if (!pd.series.empty() && std::none_of(d.series.begin(), d.series.end(), [&](const SeriesDesc& s) { return s.name == pd.series; }))
        c.err(b + ".series", "F20", "unknown series '" + pd.series + "'");
      growth(pd.growth, b);
      if (pd.growth.present() && pd.mechanics == "parallel")
        c.err(b + ".growth", "F22", "growth is not combined with parallel mechanics (offsets keep the anchor's thickness)");
    }
    // STX-TN-08 (F24, F25)
    if (d.salt.present()) {
      const std::string sb = "structure.salt";
      int ti = d.horizon_index(d.salt.top), bi = d.horizon_index(d.salt.base);
      if (ti < 0) c.err(sb + ".top", "F24", "unknown horizon '" + d.salt.top + "'");
      if (bi < 0) c.err(sb + ".base", "F24", "unknown horizon '" + d.salt.base + "'");
      if (ti >= 0 && bi >= 0 && ti != bi - 1)
        c.err(sb + ".top", "F24", "top salt must be the horizon directly above base salt (the salt is the unit between them)");
      if (d.salt.source_thickness.present && !(d.salt.source_thickness.si() > 0))
        c.err(sb + ".source_thickness", "F24", "the source-layer thickness must be positive");
      if (ti >= 0 && !d.horizons[static_cast<std::size_t>(ti)].surface_input.empty() && !d.salt.bodies.empty())
        c.err(sb + ".bodies", "F24", "top salt is an imported surface: give no parametric bodies");
      std::set<std::string> names;
      for (std::size_t k = 0; k < d.salt.bodies.size(); ++k) {
        const auto& b = d.salt.bodies[k];
        const std::string bb = sb + ".bodies[" + std::to_string(k) + "]";
        if (!names.insert(b.name).second) c.err(bb + ".name", "F24", "duplicate salt body '" + b.name + "'");
        for (auto [p, nm] : {std::pair{&b.radius, "radius"}, {&b.flank_width, "flank_width"}, {&b.height, "height"},
                             {&b.withdrawal_width, "withdrawal_width"}, {&b.half_length, "half_length"}})
          if (p->present && !(p->si() > 0)) c.err(bb + "." + nm, "F24", std::string(nm) + " must be positive");
        growth(b.growth, bb);
        if (b.centre.present && d.size.present && d.origin.present) {
          double x = b.centre.si(0), y = b.centre.si(1);
          if (x < d.origin.si(0) || y < d.origin.si(1) || x > d.origin.si(0) + d.size.si(0) || y > d.origin.si(1) + d.size.si(1))
            c.err(bb + ".centre", "F24", "a salt body lies inside the domain so its volume is counted");
        }
      }
    }
    if (d.compaction.present()) {
      const std::string cb = "stratigraphy.compaction";
      if (d.compaction.geometry != "present" && d.compaction.geometry != "depositional")
        c.err(cb + ".geometry", "F25", "geometry must be present or depositional");
      std::set<std::string> seen;
      for (std::size_t k = 0; k < d.compaction.units.size(); ++k) {
        const auto& u = d.compaction.units[k];
        const std::string ub = cb + ".units[" + std::to_string(k) + "]";
        int ui = d.unit_index(u.unit);
        if (ui < 0) c.err(ub + ".unit", "F25", "unknown unit '" + u.unit + "'");
        else if (d.horizon_index(d.units[static_cast<std::size_t>(ui)].top) < 0 || d.horizon_index(d.units[static_cast<std::size_t>(ui)].base) < 0)
          c.err(ub + ".unit", "F25", "a compacting unit lies between two horizons");
        if (!seen.insert(u.unit).second) c.err(ub + ".unit", "F25", "unit '" + u.unit + "' has two compaction laws");
        if (u.law == "exponential") {
          if (u.surface_porosity.present && !(u.surface_porosity.si() >= 0 && u.surface_porosity.si() < 1))
            c.err(ub + ".surface_porosity", "F25", "surface porosity must be in [0, 1)");
          if (u.decay.present && !(u.decay.si() > 0)) c.err(ub + ".decay", "F25", "the decay coefficient must be positive");
        } else if (u.law == "table") {
          const auto& zd = u.table_depths;
          const auto& ph = u.table_porosity;
          if (zd.present && ph.present && zd.values.size() != ph.values.size())
            c.err(ub + ".table_porosity", "F25", "table_depths and table_porosity need the same number of values");
          for (std::size_t q = 1; zd.present && q < zd.values.size(); ++q)
            if (!(zd.si(q) > zd.si(q - 1))) c.err(ub + ".table_depths", "F25", "table depths must increase");
          for (std::size_t q = 0; ph.present && q < ph.values.size(); ++q)
            if (!(ph.si(q) >= 0 && ph.si(q) < 1)) c.err(ub + ".table_porosity", "F25", "porosities must be in [0, 1)");
        } else {
          c.err(ub + ".law", "F25", "law must be exponential or table");
        }
      }
    }
    // STX-TN-09 (F26): the event timeline
    if (!d.events.empty()) {
      std::set<std::string> enames;
      for (std::size_t k = 0; k < d.events.size(); ++k) {
        const auto& e = d.events[k];
        const std::string eb = "events[" + std::to_string(k) + "]";
        if (!enames.insert(e.name).second) c.err(eb + ".name", "F26", "duplicate event name '" + e.name + "'");
        bool known = false;
        if (e.type == "deposition") known = d.horizon_index(e.ref) >= 0;
        else if (e.type == "fault") known = d.fault_index(e.ref) >= 0;
        else if (e.type == "fold") known = std::any_of(d.primitives.begin(), d.primitives.end(), [&](const PrimitiveDesc& p) { return p.name == e.ref; });
        else if (e.type == "salt") known = std::any_of(d.salt.bodies.begin(), d.salt.bodies.end(), [&](const SaltBodyDesc& b) { return b.name == e.ref; });
        if (!known) c.err(eb, "F26", "event '" + e.name + "' refers to an unknown " + e.type + " '" + e.ref + "'");
        if (e.from.present && e.to.present && !(e.from.si() > e.to.si())) c.err(eb + ".from", "F26", "from must be older (larger age) than to");
        if ((e.throw_.present || !e.sense.empty()) && e.type != "fault") c.err(eb, "F26", "throw and sense belong to fault events");
      }
      for (std::size_t h = 0; h < d.horizons.size(); ++h) {
        const auto& hd = d.horizons[h];
        if (!hd.age.present) c.err("stratigraphy.horizons[" + std::to_string(h) + "]", "F26",
                                   "horizon '" + hd.name + "' has no age: a timeline needs a deposition event or an age for every horizon");
        else if (h > 0 && d.horizons[h - 1].age.present && !(hd.age.si() > d.horizons[h - 1].age.si()))
          c.err("stratigraphy.horizons[" + std::to_string(h) + "]", "F26", "ages must increase downwards: '" + hd.name +
                                                                              "' is not older than '" + d.horizons[h - 1].name + "'");
      }
      auto timed = [&](const std::string& type, const std::string& ref, const std::string& b, bool has_growth) {
        auto ev = d.events_for(type, ref);
        if (ev.empty()) c.err(b, "F26", "in a timeline every " + type + " has an event; '" + ref + "' has none");
        if (has_growth && !ev.empty()) c.err(b + ".growth", "F26", "give the timing of '" + ref + "' by its events, not by growth");
        return ev;
      };
      for (std::size_t k = 0; k < d.primitives.size(); ++k)
        timed("fold", d.primitives[k].name, "structure.primitives[" + std::to_string(k) + "]", d.primitives[k].growth.present());
      for (std::size_t k = 0; k < d.salt.bodies.size(); ++k)
        timed("salt", d.salt.bodies[k].name, "structure.salt.bodies[" + std::to_string(k) + "]", d.salt.bodies[k].growth.present());
      for (std::size_t k = 0; k < d.faults.size(); ++k) {
        const auto& fd = d.faults[k];
        const std::string b = "faults[" + std::to_string(k) + "]";
        auto ev = timed("fault", fd.name, b, fd.growth.present());
        if (ev.size() > 1)
          for (auto* e : ev)
            if (!e->throw_.present || (e->sense != "normal" && e->sense != "reverse"))
              c.err(b, "F26", "fault '" + fd.name + "' moves in phases: event '" + e->name + "' needs its throw and sense (normal or reverse)");
        // level-dependent displacement needs layer-rigid hanging walls
        bool varies = ev.size() > 1;
        for (auto* e : ev)
          for (auto& hd : d.horizons) varies = varies || (hd.age.present && hd.age.si() < e->start());
        if (varies && fd.hw_type != "layer_rigid")
          c.err(b + ".hanging_wall.type", "F26", "fault '" + fd.name + "' displaces horizons by different amounts (phases, or deposition during or after "
                                                 "its events): it takes hanging_wall: layer_rigid");
      }
    }
    for (std::size_t k = 0; k < d.faults.size(); ++k) {
      const auto& fd = d.faults[k];
      const std::string b = "faults[" + std::to_string(k) + "]";
      growth(fd.growth, b);
      if (fd.growth.present()) {
        if (fd.hw_type != "layer_rigid")
          c.err(b + ".hanging_wall.type", "F22", "growth faults take hanging_wall: layer_rigid (each layer moves by its own throw)");
        int en = d.horizon_index(fd.growth.end);
        if (fd.upper.rule == "horizon" && en >= 0 && d.horizon_index(fd.upper.horizon) > en)
          c.err(b + ".extent.upper_tip", "F22", "the upper tip of a growth fault must be free or above the end horizon");
      }
    }
  }

  void semantic() {
    // domain
    if (d.size.present && d.size.values.size() == 2 && (d.size.si(0) <= 0 || d.size.si(1) <= 0))
      c.err("domain.size", "D1", "domain size must be positive");
    if (d.depth_range.present && d.depth_range.values.size() == 2 && d.depth_range.si(0) >= d.depth_range.si(1))
      c.err("domain.depth_range", "D1", "depth range must be [top depth, base depth] with top above base");
    if (d.grid_spacing.present && d.grid_spacing.si() <= 0) c.err("domain.grid_spacing", "D1", "grid spacing must be positive");
    // horizons S1, S3
    std::set<std::string> hn;
    for (std::size_t i = 0; i < d.horizons.size(); ++i) {
      auto& h = d.horizons[i];
      std::string b = "stratigraphy.horizons[" + std::to_string(i) + "]";
      if (!hn.insert(h.name).second) c.err(b + ".name", "S1", "duplicate horizon name '" + h.name + "'");
      int defs = (h.depth.present ? 1 : 0) + (h.thickness_above.present ? 1 : 0) + (!h.surface_input.empty() ? 1 : 0) +
                 (!h.iso_from.empty() ? 1 : 0);
      bool from_well = false;
      if (!d.ref_well.empty())
        if (auto* w = d.well(d.ref_well))
          for (auto& [n, p] : w->tops) from_well = from_well || n == h.name;
      if (from_well && (h.depth.present || h.thickness_above.present))
        c.err(b, "S3", "horizon '" + h.name + "' is defined twice: by the reference well top and by depth/thickness");
      const int si = d.series_of(static_cast<int>(i));
      const bool clino = (si >= 0 && !d.series[static_cast<std::size_t>(si)].clino_profile.empty()) ||
                         (d.salt.present() && h.name == d.salt.top && h.surface_input.empty());
      if (clino && (defs > 0 || from_well))
        c.err(b, h.name == d.salt.top ? "F24" : "F23",
              "horizon '" + h.name + "' follows from " + (h.name == d.salt.top ? "the salt thickness above base salt" : "the clinoform formula") +
                  ": give no depth, thickness, isochore or reference-well top");
      else if (clino) {
      } else if (!from_well && defs == 0)
        c.err(b, "S3", "horizon '" + h.name + "' has no definition (reference-well top, depth, thickness_above, surface or isochore)");
      else if (defs > 1)
        c.err(b, "S3", "horizon '" + h.name + "' is defined more than once");
      if (i == 0 && h.thickness_above.present) c.err(b + ".thickness_above", "S3", "the first horizon cannot be defined by a thickness");
      if (!h.iso_from.empty()) {
        int j = d.horizon_index(h.iso_from);
        if (j < 0) c.err(b + ".isochore.from", "R1", "unknown horizon '" + h.iso_from + "'");
        else if (j >= static_cast<int>(i)) c.err(b + ".isochore.from", "S3", "an isochore must hang from a shallower horizon");
        if (!h.iso_thickness.present && h.iso_input.empty())
          c.err(b + ".isochore", "S3", "isochore needs a thickness or an input grid");
      }
      for (auto* inp : {&h.surface_input, &h.iso_input})
        if (!inp->empty()) {
          bool found = false;
          for (auto& x : d.inputs) found = found || x.id == *inp;
          if (!found) c.err(b, "R1", "input '" + *inp + "' is not listed under inputs");
        }
    }
    auto zref = reference_depths();
    for (std::size_t i = 1; i < zref.size(); ++i)
      if (!std::isnan(zref[i]) && !std::isnan(zref[i - 1]) && zref[i] < zref[i - 1] &&
          d.series_of(static_cast<int>(i)) == d.series_of(static_cast<int>(i - 1)))
        c.err("stratigraphy.horizons[" + std::to_string(i) + "]", "S1",
              "horizon '" + d.horizons[i].name + "' is shallower than '" + d.horizons[i - 1].name + "' at the reference location");
    if (d.depth_range.present && d.depth_range.values.size() == 2)
      for (std::size_t i = 0; i < zref.size(); ++i)
        if (!std::isnan(zref[i]) && (zref[i] < d.depth_range.si(0) || zref[i] > d.depth_range.si(1)))
          c.err("stratigraphy.horizons[" + std::to_string(i) + "]", "S4",
                "horizon '" + d.horizons[i].name + "' lies outside the domain depth range at the reference location");
    series_rules();
    for (auto& n : d.isopachs)
      if (d.unit_index(n) < 0) c.err("outputs.isopachs", "O1", "unknown unit '" + n + "'");
    // reference
    if (d.ref_well.empty() && !d.ref_point.present)
      c.err("stratigraphy.reference", "R1", "give a reference well or a reference point");
    if (!d.ref_well.empty() && !d.well(d.ref_well))
      c.err("stratigraphy.reference.well", "R1", "unknown well '" + d.ref_well + "'");
    // units S2
    std::vector<int> cover(d.horizons.size() > 0 ? d.horizons.size() - 1 : 0, 0);
    std::set<std::string> un;
    for (std::size_t k = 0; k < d.units.size(); ++k) {
      auto& u = d.units[k];
      std::string b = "stratigraphy.units[" + std::to_string(k) + "]";
      if (!un.insert(u.name).second) c.err(b + ".name", "S2", "duplicate unit name '" + u.name + "'");
      if (u.base == "model_base" || u.top == "model_top") {
        bool below = u.base == "model_base";
        const std::string& edge = below ? u.top : u.base;
        int e = d.horizon_index(edge);
        if (u.base == "model_base" && u.top == "model_top")
          c.err(b, "S2", "unit '" + u.name + "' cannot span the whole model");
        else if (e != (below ? static_cast<int>(d.horizons.size()) - 1 : 0))
          c.err(b, "S2", std::string("a unit with ") + (below ? "base model_base must have the last horizon as its top"
                                                                 : "top model_top must have the first horizon as its base"));
        else if ((below ? d.below_base_name : d.above_top_name) != u.name)
          c.err(b, "S2", "the unit " + std::string(below ? "below the last" : "above the first") + " horizon is named twice");
        continue;
      }
      int t = d.horizon_index(u.top), bb = d.horizon_index(u.base);
      if (t < 0) { c.err(b + ".top", "R1", "unknown horizon '" + u.top + "'"); continue; }
      if (bb < 0) { c.err(b + ".base", "R1", "unknown horizon '" + u.base + "'"); continue; }
      if (bb <= t) { c.err(b, "S2", "unit '" + u.name + "' has its base above or at its top"); continue; }
      for (int i = t; i < bb; ++i) cover[static_cast<std::size_t>(i)]++;
      if (bb != t + 1)
        c.err(b, "S2", "unit '" + u.name + "' spans more than one interval; in v1 each unit lies between consecutive horizons");
    }
    for (std::size_t i = 0; i < cover.size(); ++i) {
      if (cover[i] > 1)
        c.err("stratigraphy.units", "S2", "overlapping units between '" + d.horizons[i].name + "' and '" +
                                              d.horizons[i + 1].name + "'");
      else if (cover[i] == 0)
        c.err("stratigraphy.units", "S2", "no unit between '" + d.horizons[i].name + "' and '" + d.horizons[i + 1].name + "'");
    }
    // wells W1
    for (std::size_t k = 0; k < d.wells.size(); ++k) {
      auto& w = d.wells[k];
      for (auto& [h, p] : w.tops)
        if (d.horizon_index(h) < 0)
          c.err("wells[" + std::to_string(k) + "].tops." + h, "W1", "unknown horizon '" + h + "'");
      if (w.role == "control" && w.tops.empty())
        c.err("wells[" + std::to_string(k) + "]", "W2", "a control well must have tops");
    }
    {
      int controls = 0;
      for (auto& w : d.wells) controls += (w.role == "control" && !w.tops.empty()) ? 1 : 0;
      bool imported = false;
      for (auto& h : d.horizons) imported = imported || !h.surface_input.empty() || !h.iso_input.empty();
      if ((controls > 1 || (imported && controls > 0)) && !d.correction_range.present)
        c.err("structure.correction_range", "W3",
              "control wells beyond the reference (or imported surfaces) are honoured by a correction field; "
              "give its range as a parameter");
    }
    // primitives
    for (std::size_t k = 0; k < d.primitives.size(); ++k) {
      const auto& pd = d.primitives[k];
      const std::string b = "structure.primitives[" + std::to_string(k) + "]";
      for (auto& [h, p] : pd.relief)
        if (d.horizon_index(h) < 0) c.err(b + ".relief." + h, "R1", "unknown horizon '" + h + "'");
      // STX-TN-06 sections 3 to 5: admissible range (F19)
      if (pd.type == "detachment_fold") {
        if (pd.shortening.present && pd.shortening.si() < 0) c.err(b + ".shortening", "F19", "shortening must not be negative");
        if (pd.half_width.present && !(pd.half_width.si() > 0)) c.err(b + ".half_width", "F19", "half_width must be positive");
        int hd = d.horizon_index(pd.detachment);
        if (!pd.detachment.empty() && hd < 0) c.err(b + ".detachment", "F19", "unknown horizon '" + pd.detachment + "'");
        else if (hd == 0)
          c.err(b + ".detachment", "F19", "the detachment must lie below the folded horizons; " + pd.detachment +
                                              " is the shallowest horizon, so nothing would fold");
        if (!pd.mechanics.empty()) c.err(b + ".mechanics", "F19", "a detachment fold's shape follows from the area balance; give no mechanics");
      }
      if (pd.mechanics == "parallel") {
        if (pd.relief.size() != 1)
          c.err(b + ".mechanics", "F19", "parallel mechanics needs exactly one anchor horizon with relief (the others are its offsets)");
        else if (int ha = d.horizon_index(pd.relief[0].first); ha >= 0)
          for (std::size_t h = static_cast<std::size_t>(ha) + 1; h < d.horizons.size(); ++h)
            if (!d.horizons[h].surface_input.empty() || !d.horizons[h].iso_from.empty())
              c.err(b + ".mechanics", "F19", "parallel mechanics is not combined with imported surfaces or isochores below the anchor (" +
                                                 d.horizons[h].name + ")");
      }
      if (pd.array_count > 1 && pd.array_step_along.present && pd.array_step_across.present &&
          pd.array_step_along.si() == 0 && pd.array_step_across.si() == 0)
        c.err(b + ".array", "F19", "an en echelon array needs a non-zero step along or across the axis");
    }
    // implicit (data-conditioned) builder
    if (d.implicit.present) {
      for (std::size_t k = 0; k < d.implicit.data.size(); ++k) {
        auto& x = d.implicit.data[k];
        std::string b = "implicit.data[" + std::to_string(k) + "]";
        if (std::none_of(d.inputs.begin(), d.inputs.end(), [&](const InputDesc& i) { return i.id == x.input; }))
          c.err(b + ".input", "I1", "input '" + x.input + "' is not listed under inputs");
        if (!x.horizon.empty() && d.horizon_index(x.horizon) < 0) c.err(b + ".horizon", "I2", "unknown horizon '" + x.horizon + "'");
      }
      if (d.implicit.nugget.present && d.implicit.nugget.si() < 0) c.err("implicit.nugget", "I3", "nugget must be zero or positive");
      if (d.implicit.surface_spacing.present && d.implicit.surface_spacing.si() <= 0)
        c.err("implicit.surface_spacing", "I3", "surface spacing must be positive");
      if (d.implicit.vertical_stretch.present && !(d.implicit.vertical_stretch.si() >= 1 && d.implicit.vertical_stretch.si() <= 100))
        c.err("implicit.vertical_stretch", "I3", "vertical stretch must be between 1 and 100");
      if (d.implicit.iterations.present && (d.implicit.iterations.si() < 1 || d.implicit.iterations.si() > 10))
        c.err("implicit.iterations", "I3", "iterations must be between 1 and 10");
    }
    // faults
    std::map<double, std::vector<std::size_t>> ages;
    double base_depth = d.depth_range.present && d.depth_range.values.size() == 2 ? d.depth_range.si(1) : kNaN;
    for (std::size_t k = 0; k < d.faults.size(); ++k) {
      auto& f = d.faults[k];
      std::string b = "faults[" + std::to_string(k) + "]";
      if (f.geom_type == "sticks") {
        if (f.geom_input.empty()) c.err(b + ".geometry.input", "F9", "a fault given by sticks needs geometry.input (an input id)");
        else if (std::none_of(d.inputs.begin(), d.inputs.end(), [&](const InputDesc& i) { return i.id == f.geom_input; }))
          c.err(b + ".geometry.input", "F9", "input '" + f.geom_input + "' is not listed under inputs");
        if (!d.implicit.present) c.err(b + ".geometry.type", "F9", "faults given by sticks need the implicit (data-conditioned) builder: add an implicit section");
      }
      const bool fpf = f.hw_type == "fault_parallel_flow";
      const bool shear_hw = f.hw_type == "vertical_shear" || f.hw_type == "inclined_shear" || fpf;
      if (f.curved()) {
        // STX-TN-05 section 2.3: the admissible range of the shear construction.
        if (d.implicit.present)
          c.err(b + ".geometry.type", "F14", "curved faults are built by the parametric engine; the implicit builder takes planar faults");
        if (!shear_hw)
          c.err(b + ".hanging_wall.type", "F14", "a curved fault needs a shear or flow hanging wall (vertical_shear, inclined_shear or "
                                                 "fault_parallel_flow): "
                                                 "a rigid hanging wall would open a void above the curved fault");
        if (f.slip_mode != "dip_slip")
          c.err(b + ".slip_mode", "F14", "a curved fault moves by dip_slip");
        if (f.upper.rule != "free" || f.lower.rule != "free")
          c.err(b + ".extent", "F10", "a curved fault cuts from the model top to its detachment: both tips must be free "
                                      "(an upper tip inside the section would need the displacement to fall to zero)");
        if (f.trace.present)
          c.err(b + ".geometry.trace", "F6", "a curved fault is placed by reference_point and strike, not by a trace");
        double theta_max = 0;  // steepest dip of the fault inside the model
        const double z_top = d.depth_range.present && !d.depth_range.values.empty() ? d.depth_range.si(0) : 0.0;
        if (f.geom_type == "listric" && f.dip.present && f.detachment_depth.present && f.reference_depth.present) {
          double z0 = f.reference_depth.si(), zd = f.detachment_depth.si(), th = f.dip.si();
          if (!(th > 0 && th < kPi / 2)) c.err(b + ".geometry.dip", "F6", "the dip of a listric fault must be in (0, 90) degrees");
          if (!(zd > z0)) c.err(b + ".geometry.detachment_depth", "F6", "the detachment must lie below the reference depth");
          else theta_max = std::atan(std::tan(th) * (zd - std::min(z0, z_top)) / (zd - z0));
        }
        if (f.geom_type == "ramp_flat" && f.profile.present && f.profile.cols() == 2) {
          const std::size_t n = f.profile.rows();
          if (!fpf && std::fabs(f.profile.si(0)) > 1e-9)
            c.err(b + ".geometry.profile", "F6", "the first profile point must be at u = 0 (the reference point)");
          if (fpf && !(f.profile.si(0) <= 0 && f.profile.si(2 * n - 2) >= 0))
            c.err(b + ".geometry.profile", "F18", "the profile must span u = 0 (the reference point)");
          if (!fpf && f.reference_depth.present && std::fabs(f.reference_depth.si() - f.profile.si(1)) > 1e-6)
            c.err(b + ".geometry.reference_depth", "F6", "reference_depth differs from the depth of the first profile point");
          for (std::size_t k = 1; k < n; ++k) {
            double du = f.profile.si(2 * k) - f.profile.si(2 * k - 2), dz = f.profile.si(2 * k + 1) - f.profile.si(2 * k - 1);
            if (!(du > 0)) c.err(b + ".geometry.profile", "F6", "profile u must increase from point to point");
            else if (dz < 0) c.err(b + ".geometry.profile", "F6", "the fault must not rise in its dip direction (depth must not decrease)");
            else theta_max = std::max(theta_max, std::atan2(dz, du));
          }
          if (!fpf && n >= 2 && !(f.profile.si(3) > f.profile.si(1)))
            c.err(b + ".geometry.profile", "F6", "the first profile segment must dip (it continues upwards to the model top)");
        }
        if (f.hw_type == "inclined_shear" && f.shear_angle.present) {
          double a = f.shear_angle.si();
          if (!(a >= 0 && a < kPi / 2)) c.err(b + ".hanging_wall.shear_angle", "F12", "the shear angle must be in [0, 90) degrees");
          if (f.shear_sense != "antithetic" && f.shear_sense != "synthetic")
            c.err(b + ".hanging_wall.shear_sense", "F12", "shear_sense must be antithetic or synthetic");
          if (f.shear_sense == "synthetic" && theta_max > 0 && a >= kPi / 2 - theta_max - 1e-12)
            c.err(b + ".hanging_wall.shear_angle", "F12",
                  "synthetic shear at " + fmt_fixed(a * 180 / kPi, 2) + " degrees is not steeper than the fault (steepest dip " +
                      fmt_fixed(theta_max * 180 / kPi, 2) + " degrees): shear lines would cross the fault twice. Use an angle below " +
                      fmt_fixed(90 - theta_max * 180 / kPi, 2) + " degrees or antithetic shear");
        }
        if (f.heave.present && f.heave.si() < 0) c.err(b + ".throw.heave", "F14", "heave must not be negative");
      } else if (shear_hw) {
        c.err(b + ".hanging_wall.type", fpf ? "F18" : "F14",
              fpf ? "fault-parallel flow needs a curved (ramp_flat or listric) fault: flow lines follow its bends"
                  : "vertical and inclined shear hanging walls belong to curved (listric or ramp_flat) faults");
      }
      if (f.trace.present && f.dip_side.empty())
        c.err(b + ".geometry.dip_side", "F6", "a fault given by its trace needs dip_side");
      if (f.dip.present && (f.dip.si() <= 0 || f.dip.si() > kPi / 2 + 1e-12))
        c.err(b + ".geometry.dip", "F6", "dip must be in (0, 90] degrees");
      if (f.length.present && f.length.si() <= 0) c.err(b + ".extent.length", "F1", "length must be positive");
      if (f.age.present) ages[f.age.si()].push_back(k);
      auto tip_depth = [&](const TipDesc& t, const std::string& tb) -> double {
        if (t.rule == "depth") return t.depth.present ? t.depth.si() : kNaN;
        if (t.rule == "horizon") {
          int h = d.horizon_index(t.horizon);
          if (h < 0) { c.err(tb + ".horizon", "F1", "unknown horizon '" + t.horizon + "'"); return kNaN; }
          return t.offset.present ? zref[static_cast<std::size_t>(h)] + t.offset.si() : kNaN;
        }
        if (t.rule == "within_unit") {
          int u = d.unit_index(t.unit);
          if (u < 0) { c.err(tb + ".unit", "F1", "unknown unit '" + t.unit + "'"); return kNaN; }
          if (t.fraction.present && (t.fraction.si() < 0 || t.fraction.si() > 1))
            c.err(tb + ".fraction", "F1", "fraction must be between 0 and 1");
          int ht = d.horizon_index(d.units[static_cast<std::size_t>(u)].top), hb = d.horizon_index(d.units[static_cast<std::size_t>(u)].base);
          if (ht < 0 || hb < 0) {
            c.err(tb + ".unit", "F1", "unit '" + t.unit + "' lies outside the horizons; a tip within a unit needs a unit between two horizons");
            return kNaN;
          }
          if (!t.fraction.present) return kNaN;
          double a = zref[static_cast<std::size_t>(ht)], bb = zref[static_cast<std::size_t>(hb)];
          return a + t.fraction.si() * (bb - a);
        }
        return kNaN;  // free
      };
      double du = tip_depth(f.upper, b + ".extent.upper_tip");
      double dl = tip_depth(f.lower, b + ".extent.lower_tip");
      for (auto [dd, nm] : {std::pair{du, "upper_tip"}, std::pair{dl, "lower_tip"}})
        if (!std::isnan(dd) && !std::isnan(base_depth) && dd > base_depth)
          c.err(b + ".extent." + nm, "F2", "fault tip at " + fmt_fixed(dd, 1) + " m depth lies below the model base (" +
                                                fmt_fixed(base_depth, 1) + " m)");
      if (!std::isnan(du) && !std::isnan(dl) && du >= dl)
        c.err(b + ".extent", "F3", "upper tip is not above the lower tip at the reference location");
      if (f.upper.rule == "free" && f.lower.rule == "free")
        c.warn(b + ".extent", "F3", "both tips free: the fault cuts the whole model");
      {
        // STX-TN-05 section 5: movement sense.
        const bool ss = f.slip_mode == "strike_slip", ob = f.slip_mode == "oblique";
        if ((ss || ob) && f.throw_max.present)
          c.err(b + ".throw.max", "F16", "an " + f.slip_mode + " fault takes its net slip as throw.slip, not throw.max");
        if (!ss && !ob && f.net_slip.present && !(fpf && f.curved()))
          c.err(b + ".throw.slip", "F16", "throw.slip is the net slip of oblique and strike_slip faults; dip-slip faults take throw.max");
        if (ss && f.sense != "left_lateral" && f.sense != "right_lateral")
          c.err(b + ".throw.sense", "F16", "a strike_slip fault needs throw.sense: left_lateral or right_lateral");
        if ((f.slip_mode == "dip_slip" || f.slip_mode == "vertical_offset") && !f.sense.empty() && f.sense != "normal" && f.sense != "reverse")
          c.err(b + ".throw.sense", "F16", "a " + f.slip_mode + " fault's sense is normal or reverse");
        if (ob && !f.sense.empty()) c.err(b + ".throw.sense", "F16", "an oblique fault's sense follows from its rake; give no throw.sense");
        if (!ob && f.rake.present) c.err(b + ".rake", "F16", "rake belongs to oblique faults");
        if (ob && f.rake.present && std::fabs(f.rake.si()) > kPi + 1e-12) c.err(b + ".rake", "F16", "rake must be in [-180, 180] degrees");
        if ((ss || ob) && f.hw_type != "rigid")
          c.err(b + ".hanging_wall.type", "F16", "a fault with a strike component of slip takes hanging_wall: rigid (the scaled fields assume dip-slip motion)");
        if (f.net_slip.present && f.net_slip.si() < 0) c.err(b + ".throw.slip", "F16", "net slip must not be negative (the sense gives the direction)");
        if (f.throw_max.present && f.throw_max.si() < 0 && f.sense == "reverse")
          c.err(b + ".throw.max", "F16", "give a positive throw with sense: reverse");
        if (f.curved() && !fpf && f.sense == "reverse")
          c.err(b + ".throw.sense", "F14", "shear hanging walls collapse onto normal curved faults; contractional curved faults take hanging_wall: fault_parallel_flow");
      }
      if (f.fold_type == "trishear") {
        // STX-TN-06 section 1: admissible range of trishear.
        if (f.fold_apical_angle.present && !(f.fold_apical_angle.si() > 0 && f.fold_apical_angle.si() < kPi))
          c.err(b + ".fold.apical_angle", "F17", "the apical angle must be in (0, 180) degrees");
        if (f.fold_ps_ratio.present && f.fold_ps_ratio.si() < 0)
          c.err(b + ".fold.propagation_slip_ratio", "F17", "the propagation-to-slip ratio must not be negative");
        if (f.upper.rule == "free")
          c.err(b + ".extent.upper_tip", "F17", "a trishear fault needs an upper tip inside the model: the fold forms ahead of it");
        if (f.geom_type != "planar" || f.trace.present)
          c.err(b + ".geometry", "F17", "trishear applies to planar faults placed by reference point, strike and dip");
        if (f.slip_mode != "dip_slip") c.err(b + ".slip_mode", "F17", "trishear needs dip_slip");
        if (f.hw_type != "rigid") c.err(b + ".hanging_wall.type", "F17", "trishear moves the hanging wall rigidly behind the zone: give hanging_wall: rigid");
        if (!f.drag_type.empty() || !f.relay_fault.empty() || !f.abuts_fault.empty())
          c.err(b + ".fold", "F17", "trishear is not combined with drag, relays or abutting in stage 2B");
      }
      if (!f.drag_type.empty()) {
        // STX-TN-05 section 3: drag moves material parallel to the slip, inside the fault surface.
        if (f.curved()) c.err(b + ".drag", "F13", "drag applies to planar faults (curved faults deform their hanging wall by shear)");
        if (f.slip_mode == "vertical_offset")
          c.err(b + ".drag", "F13", "drag needs the slip to lie in the fault surface (dip_slip, oblique or strike_slip)");
        if (f.hw_type == "layer_rigid") c.err(b + ".drag", "F13", "drag is not combined with layer_rigid hanging walls");
        if (f.drag_width.present && !(f.drag_width.si() > 0)) c.err(b + ".drag.width", "F13", "the drag width must be positive");
        if (f.drag_fraction.present) {
          double fr = f.drag_fraction.si();
          if (f.drag_type == "normal" && !(fr >= 0 && fr < 0.5))
            c.err(b + ".drag.fraction", "F13", "normal drag needs a fraction in [0, 0.5): at " + fmt_fixed(fr, 3) +
                                                   " a horizon would be carried across the fault (apparent throw T (1 - 2f))");
          if (f.drag_type == "reverse" && !(fr >= 0 && fr <= 1))
            c.err(b + ".drag.fraction", "F13", "reverse drag needs a fraction in [0, 1]");
        }
      }
      if (!f.relay_fault.empty()) {
        int o = d.fault_index(f.relay_fault);
        if (o < 0) c.err(b + ".relay.fault", "R1", "unknown fault '" + f.relay_fault + "'");
        else if (static_cast<std::size_t>(o) == k) c.err(b + ".relay.fault", "F15", "a fault cannot relay with itself");
        else {
          const auto& g = d.faults[static_cast<std::size_t>(o)];
          if (g.slip_mode != f.slip_mode) c.err(b + ".relay", "F15", "relay faults must share the slip mode");
          if (f.curved() || g.curved()) c.err(b + ".relay", "F15", "relays are between planar faults in stage 2A");
          if (!g.relay_fault.empty()) c.err(b + ".relay", "F15", "declare a relay on one of the two faults only");
          if (f.strike.present && g.strike.present && f.reference_point.present && g.reference_point.present &&
              f.length.present && g.length.present) {
            // overlap along F's strike between the two strike extents
            double st = f.strike.si(), sx = std::sin(st), sy = std::cos(st);
            double c2 = (g.reference_point.si(0) - f.reference_point.si(0)) * sx + (g.reference_point.si(1) - f.reference_point.si(1)) * sy;
            double lo = std::max(-0.5 * f.length.si(), c2 - 0.5 * g.length.si()), hi = std::min(0.5 * f.length.si(), c2 + 0.5 * g.length.si());
            if (!(hi > lo)) c.err(b + ".relay", "F15", "the two faults do not overlap along strike");
            else if (lo <= -0.5 * f.length.si() + 1e-6 && hi >= 0.5 * f.length.si() - 1e-6)
              c.err(b + ".relay", "F15", "one fault lies wholly inside the other's strike extent; a relay needs overlapping tips");
            double ds = std::fabs(std::remainder(st - g.strike.si(), 2 * kPi));
            if (ds > kPi / 6) c.err(b + ".relay", "F15", "relay faults must strike within 30 degrees of each other");
          }
          if (f.relay_breached)
            c.err(b + ".relay.breached", "F15", "breached relays are not implemented yet: declare the connecting fault "
                                                "explicitly, abutting both faults (open item in STX-TN-05)");
        }
      }
      if (!f.abuts_fault.empty()) {
        int a = d.fault_index(f.abuts_fault);
        if (a < 0) c.err(b + ".abuts.fault", "R1", "unknown fault '" + f.abuts_fault + "'");
        else if (static_cast<std::size_t>(a) == k) c.err(b + ".abuts.fault", "R1", "a fault cannot abut itself");
        else if (f.age.present && d.faults[static_cast<std::size_t>(a)].age.present &&
                 d.faults[static_cast<std::size_t>(a)].age.si() > f.age.si())
          c.err(b + ".abuts", "F4", "a fault can only abut an older (or equal-age) fault");
      }
    }
    for (auto& [age, ks] : ages)
      if (ks.size() > 1)
        for (auto k : ks)
          if (d.faults[k].abuts_fault.empty())
            c.err("faults[" + std::to_string(k) + "].age", "F4",
                  "fault '" + d.faults[k].name + "' shares its age with another fault and declares no abutting relation");
    if (d.faults.size() > 64) c.err("faults", "F8", "Stratonix v1 supports at most 64 faults per model");
    // frame
    if (!d.frame_well.empty() && !d.well(d.frame_well)) c.err("frame.fixed.well", "R1", "unknown well '" + d.frame_well + "'");
    if (!d.faults.empty() && d.frame_well.empty() && !d.frame_point.present && d.ref_well.empty())
      c.err("frame", "R2", "a faulted model needs a frame rule: a fixed well or a fixed point");
    if (!d.closure_horizon.empty() && d.horizon_index(d.closure_horizon) < 0)
      c.err("checks.closure.horizon", "R1", "unknown horizon '" + d.closure_horizon + "'");
    for (auto& s : d.sections)
      if (s.polyline.present && s.polyline.rows() < 2) c.err("outputs.sections", "S0", "a section needs at least two points");
    // completions W4
    for (std::size_t k = 0; k < d.wells.size(); ++k)
      for (std::size_t q = 0; q < d.wells[k].completions.size(); ++q) {
        auto& cp = d.wells[k].completions[q];
        if (cp.top.present && cp.base.present && cp.base.si() <= cp.top.si())
          c.err("wells[" + std::to_string(k) + "].completions[" + std::to_string(q) + "]", "W4",
                "a completion's base must lie below its top");
      }
    // section sheets O1 (references) and O2 (windows)
    auto unit_between = [&](const std::string& n) {
      int u = d.unit_index(n);
      return u >= 0;
    };
    for (std::size_t k = 0; k < d.sheets.size(); ++k) {
      auto& sh = d.sheets[k];
      std::string b = "outputs.section_sheets[" + std::to_string(k) + "]";
      bool sec = false;
      for (auto& s : d.sections) sec = sec || s.name == sh.section;
      if (!sec) c.err(b + ".section", "O1", "unknown section '" + sh.section + "' (give it under outputs.sections)");
      if (!sh.main.present || !sh.main.vertical_exaggeration.present)
        c.err(b + ".main", "O2", "a section sheet needs main.vertical_exaggeration");
      if (sh.detail.present && (!sh.detail.window.present || !sh.detail.vertical_exaggeration.present))
        c.err(b + ".detail", "O2", "a detail panel needs a window [distance from, distance to, depth top, depth base] and a vertical exaggeration");
      if (sh.detail.window.present && sh.detail.window.values.size() == 4 &&
          (sh.detail.window.si(1) <= sh.detail.window.si(0) || sh.detail.window.si(3) <= sh.detail.window.si(2)))
        c.err(b + ".detail.window", "O2", "window must be [distance from, distance to, depth top, depth base] in increasing order");
      for (auto* pd : {&sh.main, &sh.detail}) {
        for (auto& u : pd->label_units)
          if (!unit_between(u)) c.err(b + ".unit_labels", "O1", "unknown unit '" + u + "'");
        for (auto& a : pd->annotations) {
          if (!a.horizon.empty() && d.horizon_index(a.horizon) < 0) c.err(b + ".annotations", "O1", "unknown horizon '" + a.horizon + "'");
          if (!a.fault.empty() && d.fault_index(a.fault) < 0) c.err(b + ".annotations", "O1", "unknown fault '" + a.fault + "'");
          if (!a.well.empty() && !d.well(a.well)) c.err(b + ".annotations", "O1", "unknown well '" + a.well + "'");
          if (a.type == "leader" && a.horizon.empty()) c.err(b + ".annotations", "O1", "a leader names the horizon it points to");
          if (a.type == "horizon_at_fault" && (a.horizon.empty() || a.fault.empty()))
            c.err(b + ".annotations", "O1", "horizon_at_fault needs a horizon and a fault");
          if (a.type == "completion" && a.well.empty()) c.err(b + ".annotations", "O1", "a completion annotation names its well");
        }
      }
      for (std::size_t q = 0; q < sh.keys.size(); ++q) {
        auto& kd = sh.keys[q];
        std::string bk = b + ".key_numbers[" + std::to_string(q) + "]";
        if (!kd.horizon.empty() && d.horizon_index(kd.horizon) < 0) c.err(bk, "O1", "unknown horizon '" + kd.horizon + "'");
        if (!kd.unit.empty() && !unit_between(kd.unit)) c.err(bk, "O1", "unknown unit '" + kd.unit + "'");
        if (!kd.well.empty() && !d.well(kd.well)) c.err(bk, "O1", "unknown well '" + kd.well + "'");
        if (!kd.fault.empty() && d.fault_index(kd.fault) < 0) c.err(bk, "O1", "unknown fault '" + kd.fault + "'");
        for (auto& fn : kd.faults)
          if (d.fault_index(fn) < 0) c.err(bk, "O1", "unknown fault '" + fn + "'");
        static const std::set<std::string> kinds{"horizon_at_well", "horizon_highest", "block_offset", "fault_tips", "fault_throw",
                                                 "well_fault_distance", "block_highest_at_fault", "relief_to_end", "text",
                                                 "completion", "vertical_exaggeration", "fold_relief", "limb_dip", "fault_dip",
                                                 "fault_tip_units", "fault_cutoff_position", "fault_cutoff_depths", "fault_overlap",
                                                 "throw_ratio", "fault_gap", "growth_index"};
        if (!kinds.count(kd.kind)) { c.err(bk, "O1", "unknown key-number kind '" + kd.kind + "'"); continue; }
        if (!kd.from.empty() && d.horizon_index(kd.from) < 0) c.err(bk, "O1", "unknown horizon '" + kd.from + "'");
        if (!kd.to.empty() && d.horizon_index(kd.to) < 0) c.err(bk, "O1", "unknown horizon '" + kd.to + "'");
        static const std::set<std::string> with_h{"horizon_at_well", "horizon_highest", "block_offset", "fault_throw",
                                                  "block_highest_at_fault", "relief_to_end", "fold_relief", "limb_dip",
                                                  "fault_cutoff_position", "fault_cutoff_depths", "throw_ratio", "fault_overlap"};
        if (with_h.count(kd.kind) && kd.horizon.empty()) c.err(bk, "O1", kd.kind + " needs a horizon");
        static const std::set<std::string> with_f{"fault_tips", "fault_throw", "block_highest_at_fault", "relief_to_end", "fault_dip",
                                                  "fault_tip_units", "fault_cutoff_position", "fault_cutoff_depths", "fault_overlap",
                                                  "throw_ratio", "growth_index"};
        if (with_f.count(kd.kind) && kd.fault.empty()) c.err(bk, "O1", kd.kind + " needs a fault");
        if ((kd.kind == "block_offset" || kd.kind == "fault_gap") && kd.faults.size() != 2)
          c.err(bk, "O1", kd.kind + " needs two faults");
        if (kd.kind == "well_fault_distance" && kd.faults.empty()) c.err(bk, "O1", kd.kind + " needs faults");
        if ((kd.kind == "horizon_at_well" || kd.kind == "block_offset" || kd.kind == "well_fault_distance" ||
             kd.kind == "block_highest_at_fault" || kd.kind == "relief_to_end" || kd.kind == "completion") && kd.well.empty())
          c.err(bk, "O1", kd.kind + " needs a well");
        if (kd.kind == "completion" && kd.completion.empty()) c.err(bk, "O1", "completion needs a completion name");
        if (kd.kind == "text" && kd.value.empty()) c.err(bk, "O1", "text needs a value");
        if ((kd.kind == "fault_overlap" || kd.kind == "growth_index") && kd.unit.empty()) c.err(bk, "O1", kd.kind + " needs a unit");
        if (kd.kind == "throw_ratio" && (kd.from.empty() || kd.to.empty())) c.err(bk, "O1", "throw_ratio needs from and to horizons");
        if (kd.side != "footwall" && kd.side != "hanging_wall") c.err(bk, "O1", "side must be footwall or hanging_wall");
        if (kd.coordinate != "x" && kd.coordinate != "y" && kd.coordinate != "distance") c.err(bk, "O1", "coordinate must be x, y or distance");
        if (kd.kind == "well_fault_distance" && kd.unit.empty()) c.err(bk, "O1", "well_fault_distance needs a unit");
      }
    }
  }
};

}  // namespace

// The schema text is generated into schema_embed.cpp by CMake.
const std::string& schema_v1_text();

void migrate(Node& doc) {
  const Node* v = doc.get("stratonix_schema");
  if (!v || !v->is_number())
    throw Error(ErrorCode::Validation, "missing schema version (rule V1: write 'stratonix_schema: 1')", "stratonix_schema");
  int ver = static_cast<int>(v->number);
  if (ver > kSchemaVersion)
    throw Error(ErrorCode::Validation,
                "schema version " + std::to_string(ver) + " is newer than this Stratonix supports (" +
                    std::to_string(kSchemaVersion) + ")", "stratonix_schema");
  if (ver < 1) throw Error(ErrorCode::Validation, "unknown schema version " + std::to_string(ver), "stratonix_schema");
  // Migration chain: for (int k = ver; k < kSchemaVersion; ++k) migrations[k](doc);
  // Version 1 is the first version; there is nothing to migrate yet.
}

ValidationReport validate_description_text(std::string_view text, const std::filesystem::path& base_dir,
                                           Description* out) {
  ValidationReport rep;
  Node doc;
  try {
    doc = parse_document(text, "description");
    if (!doc.is_map()) {
      rep.issues.push_back({"(root)", "P1", "the description must be a mapping", true, 0});
      return rep;
    }
    migrate(doc);
  } catch (const Error& e) {
    rep.issues.push_back({e.field().empty() ? "(document)" : e.field(), "P1", e.what(), true, 0});
    return rep;
  }
  static const Node schema = parse_document(schema_v1_text(), "schema");
  RefHook hook = [](const std::string& def, const Node& v, const std::string& path) -> std::optional<Issue> {
    if (def != "param") return std::nullopt;
    Param p;
    std::string why;
    if (parse_param(v, path, p, why)) return std::nullopt;
    return Issue{path, "E1", why + "; every parameter carries value, unit, evidence class and source", true, v.line};
  };
  rep.issues = validate_schema(doc, schema, hook);
  if (!rep.ok()) return rep;
  Description d;
  Typer t{Ctx{&rep.issues}, doc, d};
  t.run();
  t.semantic();
  for (auto& in : d.inputs) {
    std::filesystem::path p(in.path);
    in.resolved = p.is_absolute() ? p : base_dir / p;
  }
  d.resolved = std::move(doc);
  d.file_sha256 = sha256_hex(text);
  if (out) *out = std::move(d);
  return rep;
}

ValidationReport validate_description(const std::filesystem::path& file, Description* out) {
  std::string text;
  try {
    text = read_file(file);
  } catch (const Error& e) {
    ValidationReport r;
    r.issues.push_back({file.generic_string(), "P0", e.what(), true, 0});
    return r;
  }
  auto r = validate_description_text(text, file.parent_path(), out);
  if (out) out->file = file;
  return r;
}

Description load_description(const std::filesystem::path& file) {
  Description d;
  auto r = validate_description(file, &d);
  if (!r.ok()) throw Error(ErrorCode::Validation, "invalid model description " + file.generic_string() + "\n" + r.text());
  for (auto& i : r.issues)
    if (!i.error) log(LogLevel::Warn, i.path + ": " + i.message);
  return d;
}

}  // namespace stratonix::desc
