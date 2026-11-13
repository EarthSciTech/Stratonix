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
                  &tol_offset_outside_tips, &tol_throw, &tol_wells, &tol_unoffset,
                  &closure_faults_seal, &target_closure, &target_area, &target_spill_depth,
                  &target_tolerance, &vertical_exaggeration, &voxel_spacing, &fault_stick_spacing, &implicit.nugget,
                  &implicit.surface_spacing, &implicit.iterations, &implicit.vertical_stretch})
    add(*p);
  for (auto& w : wells) {
    add(w.head); add(w.datum_elevation); add(w.total_depth); add(w.stations);
    for (auto& [h, p] : w.tops) add(p);
  }
  for (auto& h : horizons) { add(h.depth); add(h.thickness_above); add(h.iso_thickness); }
  for (auto& p : primitives) {
    for (auto* q : {&p.crest, &p.axis_azimuth, &p.plunge, &p.half_length, &p.half_width_left,
                    &p.half_width_right, &p.radius_major, &p.radius_minor, &p.rim_depth, &p.rim_width})
      add(*q);
    for (auto& [h, q] : p.relief) add(q);
  }
  for (auto& f : faults) {
    for (auto* q : {&f.age, &f.reference_point, &f.reference_depth, &f.strike, &f.dip, &f.trace,
                    &f.trace_depth, &f.decay_length, &f.upper.offset, &f.upper.fraction, &f.upper.depth,
                    &f.lower.offset, &f.lower.fraction, &f.lower.depth, &f.length, &f.throw_max,
                    &f.profile_dip.taper, &f.profile_strike.taper, &f.influence_distance})
      add(*q);
  }
  for (auto& s : sections) { add(s.polyline); add(s.detail_window); }
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
        }
        if (const Node* io = h.get("intended_offset")) hd.intended_offset = io->boolean;
        d.horizons.push_back(std::move(hd));
      }
    if (const Node* us = st ? st->get("units") : nullptr)
      for (auto& u : us->seq) d.units.push_back({S(&u, "name"), S(&u, "top"), S(&u, "base"), S(&u, "role"), S(&u, "color")});
    d.above_top_name = S(st, "above_top_name");
    d.below_base_name = S(st, "below_base_name");

    if (const Node* s = root.get("structure")) {
      if (const Node* r = s->get("regional")) {
        d.regional_type = S(r, "type", "flat");
        d.regional_dip = P(r, "dip", "structure.regional", Dimension::Angle, d.regional_type == "planar", 1);
        d.regional_dip_azimuth = P(r, "dip_azimuth", "structure.regional", Dimension::Angle, d.regional_type == "planar", 1);
      }
      if (const Node* bs = s->get("base")) d.base_type = S(bs, "type", "conformable");
      d.correction_range = P(s, "correction_range", "structure", Dimension::Length, false, 1);
      if (const Node* ps = s->get("primitives"))
        for (std::size_t k = 0; k < ps->seq.size(); ++k) {
          const Node& p = ps->seq[k];
          std::string b = "structure.primitives[" + std::to_string(k) + "]";
          PrimitiveDesc pd;
          pd.type = S(&p, "type");
          pd.name = S(&p, "name");
          bool fold = pd.type == "anticline" || pd.type == "syncline";
          pd.crest = P(&p, "crest", b, Dimension::Length, true, 2);
          pd.axis_azimuth = P(&p, "axis_azimuth", b, Dimension::Angle, fold, 1);
          pd.plunge = P(&p, "plunge", b, Dimension::Angle, false, 1);
          pd.half_length = P(&p, "half_length", b, Dimension::Length, false, 1);
          pd.half_width_left = P(&p, "half_width_left", b, Dimension::Length, fold, 1);
          pd.half_width_right = P(&p, "half_width_right", b, Dimension::Length, fold, 1);
          pd.radius_major = P(&p, "radius_major", b, Dimension::Length, !fold, 1);
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
        fd.geom_input = S(g, "input");
        fd.geom_input_format = S(g, "input_format", "irap_lines");
        fd.reference_point = P(g, "reference_point", b + ".geometry", Dimension::Length, !by_trace && !by_sticks, 2);
        fd.reference_depth = P(g, "reference_depth", b + ".geometry", Dimension::Length, !by_trace && !by_sticks, 1);
        fd.strike = P(g, "strike", b + ".geometry", Dimension::Angle, !by_trace && !by_sticks, 1);
        fd.dip = P(g, "dip", b + ".geometry", Dimension::Angle, !by_sticks, 1);
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
        };
        tip("upper_tip", fd.upper);
        tip("lower_tip", fd.lower);
        fd.length = P(e, "length", b + ".extent", Dimension::Length, true, 1);
        fd.tip_line_shape = S(e, "tip_line_shape", "rectangular");
        const Node* th = f.get("throw");
        fd.throw_max = P(th, "max", b + ".throw", Dimension::Length, true, 1);
        auto prof = [&](const char* key, ProfileDesc& pd) {
          const Node* pn = th ? th->get(key) : nullptr;
          pd.type = S(pn, "type", "flat_topped");
          pd.taper = P(pn, "taper", b + ".throw." + key, Dimension::Dimensionless, pn && pd.type == "flat_topped", 1);
        };
        prof("profile_dip", fd.profile_dip);
        prof("profile_strike", fd.profile_strike);
        if (!th || !th->get("profile_dip") || !th->get("profile_strike"))
          c.err(b + ".throw", "F7", "both profile_dip and profile_strike must be given (no hidden defaults)");
        fd.slip_mode = S(&f, "slip_mode");
        const Node* hw = f.get("hanging_wall");
        fd.hw_type = S(hw, "type");
        fd.influence_distance = P(hw, "influence_distance", b + ".hanging_wall", Dimension::Length, fd.hw_type == "decay", 1);
        if (const Node* ab = f.get("abuts")) { fd.abuts_fault = S(ab, "fault"); fd.abuts_side = S(ab, "side"); }
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
      d.vertical_exaggeration = P(o, "vertical_exaggeration", "outputs", Dimension::Dimensionless, false, 1);
      d.voxel_spacing = P(o, "voxel_spacing", "outputs", Dimension::Length, false, 3);
      d.fault_stick_spacing = P(o, "fault_stick_spacing", "outputs", Dimension::Length, false, 1);
      if (const Node* ex = o->get("exports"))
        for (auto& x : ex->seq) d.exports.push_back(x.str);
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
      else if (!from_well && defs == 0)
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
      if (!std::isnan(zref[i]) && !std::isnan(zref[i - 1]) && zref[i] < zref[i - 1])
        c.err("stratigraphy.horizons[" + std::to_string(i) + "]", "S1",
              "horizon '" + d.horizons[i].name + "' is shallower than '" + d.horizons[i - 1].name + "' at the reference location");
    if (d.depth_range.present && d.depth_range.values.size() == 2)
      for (std::size_t i = 0; i < zref.size(); ++i)
        if (!std::isnan(zref[i]) && (zref[i] < d.depth_range.si(0) || zref[i] > d.depth_range.si(1)))
          c.err("stratigraphy.horizons[" + std::to_string(i) + "]", "S4",
                "horizon '" + d.horizons[i].name + "' lies outside the domain depth range at the reference location");
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
    for (std::size_t k = 0; k < d.primitives.size(); ++k)
      for (auto& [h, p] : d.primitives[k].relief)
        if (d.horizon_index(h) < 0)
          c.err("structure.primitives[" + std::to_string(k) + "].relief." + h, "R1", "unknown horizon '" + h + "'");
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
      if (f.geom_type == "listric")
        c.err(b + ".geometry.type", "F5", "listric faults are reserved in schema v1 but not implemented in Stratonix v1");
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
          if (ht < 0 || hb < 0 || !t.fraction.present) return kNaN;
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
