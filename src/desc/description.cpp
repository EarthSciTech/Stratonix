// Typed model description: parameter parsing, schema validation, migrations and semantic rules.
#include "stratonix/desc/description.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/faultsystems.hpp"
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
                  &tol_offset_outside_tips, &tol_throw, &tol_wells, &tol_unoffset, &tol_rollover, &tol_area, &tol_salt_volume, &tol_solid_thickness, &tol_restoration, &tol_balance, &tol_watertight, &scaling_min_ratio, &scaling_max_ratio,
                  &closure_faults_seal, &target_closure, &target_area, &target_spill_depth,
                  &target_tolerance, &vertical_exaggeration, &voxel_spacing, &fault_stick_spacing, &implicit.nugget,
                  &implicit.surface_spacing, &implicit.iterations, &implicit.vertical_stretch, &implicit.points_per_patch, &implicit.overlap,
                  &implicit.discrete_spacing, &implicit.discrete_vertical_spacing, &implicit.discrete_regularisation, &implicit.discrete_orientation_weight, &implicit.discrete_max_iterations, &implicit.relation_distance})
    add(*p);
  for (auto& ff : implicit.folds) { add(ff.wavelength); add(ff.order); add(ff.constraint_spacing); }
  for (auto& w : wells) {
    add(w.head); add(w.datum_elevation); add(w.total_depth); add(w.stations); add(w.top_uncertainty);
    for (auto& [h, p] : w.tops) add(p);
    for (auto& c : w.completions) { add(c.top); add(c.base); }
  }
  for (auto& h : horizons) {
    add(h.depth); add(h.thickness_above); add(h.iso_thickness); add(h.iso_gradient); add(h.iso_centre);
    add(h.surface_perturbation.sigma); add(h.surface_perturbation.correlation_length);
  }
  for (auto& L : velocity.layers) {
    add(L.v0); add(L.k);
    add(L.perturbation.v0_sigma); add(L.perturbation.k_sigma); add(L.perturbation.correlation_length);
  }
  add(velocity.tolerance);
  for (auto& q : grid.geostatistics) {
    auto vario = [&](const GridDesc::VariogramDesc& v) {
      add(v.nugget);
      for (auto& st : v.structures) { add(st.contribution); add(st.ranges); add(st.azimuth); }
    };
    vario(q.variogram);
    for (auto& c : q.codes) { add(c.proportion); vario(c.variogram); }
    for (auto* p : {&q.min, &q.max, &q.mean, &q.sd, &q.decluster_cell, &q.search_radii, &q.search_azimuth, &q.max_data, &q.width, &q.thickness, &q.sinuosity,
                    &q.wavelength, &q.azimuth, &q.radius, &q.proportion, &q.max_iterations, &q.seed})
      add(*p);
  }
  add(geological_map.sample_spacing);
  for (auto& q : uncertainty.picks) { add(q.sigma); add(q.correlation_length); }
  add(uncertainty.sensitivity.base_samples); add(uncertainty.sensitivity.bootstrap); add(uncertainty.sensitivity.confidence);
  for (auto& q : uncertainty.velocity) { add(q.v0_sigma); add(q.k_sigma); add(q.correlation_length); }
  for (auto& h : horizons) add(h.age);
  for (auto& u : units) add(u.clay_fraction);
  for (auto& x : implicit.data) add(x.uncertainty);
  for (auto& x : conditioning.data) add(x.uncertainty);
  for (auto* q : {&conditioning.uncertainty, &conditioning.orientation_uncertainty, &conditioning.max_iterations, &conditioning.walkers,
                  &conditioning.steps, &conditioning.burn_in, &conditioning.thin, &conditioning.seed})
    add(*q);
  for (auto& fp : conditioning.parameters) { add(fp.min); add(fp.max); }
  for (auto* q : {&geophysics.cell, &geophysics.background_density, &geophysics.background_susceptibility, &geophysics.inclination,
                  &geophysics.declination, &geophysics.intensity, &geophysics.gravity.uncertainty, &geophysics.magnetics.uncertainty})
    add(*q);
  for (auto& u : geophysics.units) { add(u.density); add(u.susceptibility); }
  for (auto* q : {&restoration.datum, &restoration.pin, &restoration.young, &restoration.poisson, &restoration.mesh_size, &restoration.max_strain,
                  &restoration.bed_length_limit, &backstripping.mantle_density, &backstripping.water_density, &backstripping.grain_density,
                  &backstripping.te, &backstripping.young, &backstripping.poisson})
    add(*q);
  for (auto& u : restoration.units) { add(u.young); add(u.poisson); }
  for (auto& x : backstripping.steps) { add(x.water_depth); add(x.sea_level); }
  for (auto& e : eroded_sections) { add(e.surface_porosity); add(e.decay); add(e.vitrinite_estimate); }
  for (auto& u : mechanics.units) { add(u.young); add(u.poisson); add(u.ucs); add(u.biot); }
  add(mechanics.strain_min);
  add(mechanics.strain_max);
  for (auto& u : fractures.units) { add(u.p32_0); add(u.a_k); add(u.a_e); add(u.p32_max); }
  for (auto* q : {&fractures.dfn.area, &fractures.dfn.seed, &fractures.dfn.limit}) add(*q);
  for (auto& s : fractures.dfn.sets) for (auto* q : {&s.pole, &s.kappa, &s.exponent, &s.rmin, &s.rmax, &s.weight}) add(*q);
  for (auto* q : {&heat.surface_temperature, &heat.basal_heat_flow, &heat.water_conductivity, &heat.column_depth}) add(*q);
  for (auto& [a, q] : heat.history) { add(a); add(q); }
  for (auto& u : heat.units) { add(u.conductivity); add(u.heat_production); }
  for (auto* q : {&pressure.viscosity, &pressure.rho_water, &pressure.rho_grain}) add(*q);
  for (auto& u : pressure.units) { add(u.k0); add(u.b); add(u.phi0); }
  for (auto& g : geobodies)
    for (auto* q : {&g.count, &g.seed, &g.width, &g.thickness, &g.wavelength, &g.amplitude, &g.orientation, &g.radius, &g.aspect}) add(*q);
  for (auto& b : bodies)
    for (auto* q : {&b.age, &b.footprint, &b.base_offset, &b.base_depth, &b.thickness, &b.taper, &b.centre, &b.depth, &b.radius, &b.accommodation,
                    &b.uplift_decay, &b.climb, &b.outer_radius, &b.trace, &b.width, &b.top_depth, &b.sheets, &b.length, &b.length_sd, &b.thickness_sd,
                    &b.net_to_gross, &b.max_fraction, &b.seed})
      add(*q);
  for (auto* q : {&fault_rock.shale_threshold, &fault_rock.thickness_ratio, &fault_rock.seal_a, &fault_rock.seal_c, &stress.sv_gradient,
                  &stress.shmax_ratio, &stress.shmin_ratio, &stress.shmax_azimuth, &stress.tensor_gradient, &stress.pore_pressure_gradient,
                  &stress.friction, &stress.cohesion, &fracture_threshold})
    add(*q);
  for (auto& e : events)
    for (auto* q : {&e.age, &e.from, &e.to, &e.throw_}) add(*q);
  for (auto* q : {&grid.spacing, &grid.margin, &grid.min_thickness}) add(*q);
  for (auto& z : grid.zones) { add(z.layers); add(z.layer_thickness); }
  for (auto& f : grid.facies) add(f.code);
  for (auto& t : grid.transforms) { add(t.a); add(t.b); }
  add(uncertainty.members); add(uncertainty.seed);
  for (auto& u : uncertainty.parameters)
    for (auto* q : {&u.mean, &u.sd, &u.min, &u.max, &u.mode, &u.values, &u.weights}) add(*q);
  for (auto& c : uncertainty.correlations) add(c.rho);
  for (auto& sc : uncertainty.scenarios)
    for (auto& o : sc.options) add(o.weight);
  for (auto& m : uncertainty.maps) { add(m.depth); add(m.spacing); }
  add(salt.source_thickness);
  for (auto& b : salt.bodies)
    for (auto* q : {&b.centre, &b.radius, &b.flank_width, &b.height, &b.withdrawal_width, &b.axis_azimuth, &b.half_length}) add(*q);
  add(compaction.surface_depth);
  for (auto& u : compaction.units)
    for (auto* q : {&u.surface_porosity, &u.decay, &u.table_depths, &u.table_porosity}) add(*q);
  for (auto& sr : series)
    for (auto* q : {&sr.tilt_dip, &sr.tilt_azimuth, &sr.clino_height, &sr.clino_width, &sr.clino_direction, &sr.clino_first_break,
                    &sr.clino_progradation, &sr.clino_trajectory})
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
                    &f.rake, &f.fold_apical_angle, &f.fold_ps_ratio, &f.estimate_band, &f.estimate_spacing})
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

const Param* Description::param_by_name(const std::string& path) const {
  // map names to indices: faults.<name> -> faults[k], structure.primitives.<name> -> structure.primitives[k], ...
  std::string ip;
  auto take = [&](const std::string& prefix, const std::string& indexed, auto names) {
    if (path.rfind(prefix, 0) != 0) return false;
    std::string rest = path.substr(prefix.size());
    auto dot = rest.find('.');
    std::string nm = rest.substr(0, dot), tail = dot == std::string::npos ? "" : rest.substr(dot);
    for (std::size_t k = 0; k < names.size(); ++k)
      if (names[k] == nm) { ip = indexed + "[" + std::to_string(k) + "]" + tail; return true; }
    return false;
  };
  std::vector<std::string> fn, pn, hn, wn, vn;
  for (auto& L : velocity.layers) vn.push_back(L.name);
  std::vector<std::string> gn;
  for (auto& q : grid.geostatistics) gn.push_back(q.name);
  for (auto& f : faults) fn.push_back(f.name);
  for (auto& p : primitives) pn.push_back(p.name);
  for (auto& h : horizons) hn.push_back(h.name);
  for (auto& w : wells) wn.push_back(w.name);
  if (!take("faults.", "faults", fn) && !take("structure.primitives.", "structure.primitives", pn) &&
      !take("stratigraphy.horizons.", "stratigraphy.horizons", hn) && !take("wells.", "wells", wn) &&
      !take("velocity_model.layers.", "velocity_model.layers", vn) && !take("grid.geostatistics.", "grid.geostatistics", gn))
    ip = path;
  for (auto* p : all_params())
    if (p->path == ip) return p;
  return nullptr;
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
      for (auto& i : in->seq) {
        InputDesc x{S(&i, "id"), S(&i, "path"), S(&i, "role"), S(&i, "sharing"), S(&i, "source"), S(&i, "evidence"), {}, {}, false};
        if (const Node* cr = i.get("crs")) {
          x.has_crs = true;
          x.crs.name = S(cr, "name");
          if (cr->get("epsg")) x.crs.epsg = static_cast<int>(cr->get("epsg")->number);
          x.crs.wkt = S(cr, "wkt");
        }
        d.inputs.push_back(x);
      }

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
        wd.top_uncertainty = P(&w, "top_uncertainty", b, Dimension::Length, false, 1);
        if (const Node* tr = w.get("trajectory")) {
          wd.stations = P(tr, "stations", b + ".trajectory", Dimension::Any, false, 2, 3);
          wd.survey_file = S(tr, "file");
          wd.survey_method = S(tr, "method");
          wd.azimuth_reference = S(tr, "azimuth_reference");
          wd.declination = P(tr, "declination", b + ".trajectory", Dimension::Angle, false, 1);
          wd.convergence = P(tr, "convergence", b + ".trajectory", Dimension::Angle, false, 1);
          if (!wd.stations.present && wd.survey_file.empty())
            c.err(b + ".trajectory", "F44", "a trajectory needs stations or a survey file");
          for (std::size_t r = 0; wd.stations.present && r < wd.stations.rows(); ++r)
            wd.survey.push_back({wd.stations.values[3 * r] * wd.stations.to_si, deg2rad(wd.stations.values[3 * r + 1]),
                                 deg2rad(wd.stations.values[3 * r + 2])});
        }
        if (const Node* ls = w.get("logs"))
          for (auto& l : ls->seq) {
            WellLogDesc ld{S(&l, "file"), S(&l, "frame"), S(&l, "index"), {}, {}};
            if (const Node* cv = l.get("curves"))
              for (auto& x : cv->seq) ld.curves.push_back(x.str);
            wd.logs.push_back(ld);
          }
        if (const Node* as = w.get("ages"))
          for (std::size_t q = 0; q < as->seq.size(); ++q) {
            const Node& a = as->seq[q];
            const std::string ba = b + ".ages[" + std::to_string(q) + "]";
            wd.ages.push_back({S(&a, "kind"), S(&a, "name"), S(&a, "timescale"), S(&a, "source"), P(&a, "depth", ba, Dimension::Length, true, 1),
                               P(&a, "depth_base", ba, Dimension::Length, false, 1), P(&a, "age", ba, Dimension::Time, true, 1),
                               P(&a, "sigma", ba, Dimension::Time, false, 1)});
          }
        if (const Node* cs = w.get("completions"))
          for (std::size_t q = 0; q < cs->seq.size(); ++q) {
            std::string bc = b + ".completions[" + std::to_string(q) + "]";
            wd.completions.push_back({S(&cs->seq[q], "name"), P(&cs->seq[q], "top", bc, Dimension::Length, true, 1),
                                      P(&cs->seq[q], "base", bc, Dimension::Length, true, 1)});
          }
        if (const Node* tops = w.get("tops"))
          for (auto& [h, tn] : tops->map) wd.tops.emplace_back(h, P(tops, h, b + ".tops", Dimension::Length, true, 1));
        if (const Node* tt = w.get("tops_twt"))
          for (auto& [h, tn] : tt->map) wd.tops_twt.emplace_back(h, P(tt, h, b + ".tops_twt", Dimension::TravelTime, true, 1));
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
          if (const Node* pp = s->get("perturbation")) {
            auto& q = hd.surface_perturbation;
            q.present = true;
            q.sigma = P(pp, "sigma", b + ".surface.perturbation", Dimension::TravelTime, true, 1);
            q.correlation_length = P(pp, "correlation_length", b + ".surface.perturbation", Dimension::Length, false, 1);
            if (const Node* sd = pp->get("seed")) q.seed = static_cast<std::uint64_t>(sd->is_number() ? sd->number : std::stod(sd->str));
          }
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
      for (std::size_t k = 0; k < us->seq.size(); ++k) {
        const Node& u = us->seq[k];
        d.units.push_back({S(&u, "name"), S(&u, "top"), S(&u, "base"), S(&u, "role"), S(&u, "color"), S(&u, "label"),
                           P(&u, "clay_fraction", "stratigraphy.units[" + std::to_string(k) + "]", Dimension::Dimensionless, false, 1)});
      }
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
          sd.clino_trajectory = P(cf, "trajectory", cb, Dimension::Angle, false, 1);
        }
        d.series.push_back(std::move(sd));
      }
    if (const Node* cp = st ? st->get("compaction") : nullptr) {
      const std::string cb = "stratigraphy.compaction";
      d.compaction.geometry = S(cp, "geometry");
      d.compaction.burial = S(cp, "burial", "block");
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
        fd.sticks_fit = S(g, "fit", "plane");
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
        const Node* th0 = f.get("throw");
        const Node* est = th0 ? th0->get("estimate") : nullptr;
        fd.throw_estimate = est != nullptr;
        if (est) {
          fd.estimate_band = P(est, "band", b + ".throw.estimate", Dimension::Length, true, 1);
          fd.estimate_spacing = P(est, "station_spacing", b + ".throw.estimate", Dimension::Length, true, 1);
          fd.estimate_model = S(est, "model", "horizons");
        }
        fd.length = P(e, "length", b + ".extent", Dimension::Length, !fd.throw_estimate, 1);
        fd.tip_line_shape = S(e, "tip_line_shape", "rectangular");
        const Node* th = f.get("throw");
        fd.slip_mode = S(&f, "slip_mode");
        const bool fpf_hw = f.get("hanging_wall") && S(f.get("hanging_wall"), "type") == "fault_parallel_flow";
        const bool by_slip = fd.strike_component() || (curved && fpf_hw);
        fd.throw_max = P(th, "max", b + ".throw", Dimension::Length, !curved && !by_slip && !fd.throw_estimate, 1);
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
        if (const Node* lk = f.get("linking")) fd.linking = lk->type == Node::Type::Bool && lk->boolean;
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
        d.tol_watertight = P(t, "watertight", b, Dimension::Length, false, 1);
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
        for (std::size_t k = 0; k < data->seq.size(); ++k) {
          const Node& x = data->seq[k];
          d.implicit.data.push_back({S(&x, "input"), S(&x, "kind"), S(&x, "format"), S(&x, "horizon"),
                                     P(&x, "uncertainty", "implicit.data[" + std::to_string(k) + "]",
                                       S(&x, "kind") == "orientations" ? Dimension::Angle : Dimension::Length, false, 1)});
        }
      if (const Node* u = im->get("use_well_tops")) d.implicit.use_well_tops = u->boolean || u->str == "true";
      d.implicit.kernel = S(im, "kernel", "cubic");
      d.implicit.nugget = P(im, "nugget", "implicit", Dimension::Length, false, 1);
      d.implicit.surface_spacing = P(im, "surface_spacing", "implicit", Dimension::Length, false, 1);
      d.implicit.iterations = P(im, "iterations", "implicit", Dimension::Ordinal, false, 1);
      d.implicit.vertical_stretch = P(im, "vertical_stretch", "implicit", Dimension::Dimensionless, false, 1);
      if (const Node* so = im->get("solver")) {
        d.implicit.solver_type = S(so, "type", "auto");
        d.implicit.points_per_patch = P(so, "points_per_patch", "implicit.solver", Dimension::Ordinal, false, 1);
        d.implicit.overlap = P(so, "overlap", "implicit.solver", Dimension::Dimensionless, false, 1);
      }
      if (const Node* u = im->get("data_influence")) d.implicit.data_influence = u->boolean || u->str == "true";
      if (const Node* ih = im->get("isochore_horizons"))
        for (auto& x : ih->seq) d.implicit.isochore_horizons.push_back(x.str);
      d.implicit.horizon_values = S(im, "horizon_values", "fixed");
      d.implicit.fault_relations = S(im, "fault_relations", "stated");
      d.implicit.relation_distance = P(im, "relation_distance", "implicit", Dimension::Length, false, 1);
      if (const Node* fl = im->get("folds"))
        for (std::size_t k = 0; k < fl->seq.size(); ++k) {
          const Node& x = fl->seq[k];
          std::string bb = "implicit.folds[" + std::to_string(k) + "]";
          ImplicitDesc::FoldFrameDesc ff;
          ff.name = S(&x, "name");
          ff.series = S(&x, "series");
          ff.foliation = S(&x, "foliation");
          ff.wavelength = P(&x, "wavelength", bb, Dimension::Length, false, 1);
          ff.order = P(&x, "order", bb, Dimension::Ordinal, false, 1);
          ff.constraint_spacing = P(&x, "constraint_spacing", bb, Dimension::Length, false, 1);
          d.implicit.folds.push_back(ff);
        }
      if (const Node* dz = im->get("discrete")) {
        d.implicit.discrete_spacing = P(dz, "spacing", "implicit.discrete", Dimension::Length, false, 1);
        d.implicit.discrete_vertical_spacing = P(dz, "vertical_spacing", "implicit.discrete", Dimension::Length, false, 1);
        d.implicit.discrete_regularisation = P(dz, "regularisation", "implicit.discrete", Dimension::Dimensionless, false, 1);
        d.implicit.discrete_orientation_weight = P(dz, "orientation_weight", "implicit.discrete", Dimension::Dimensionless, false, 1);
        d.implicit.discrete_max_iterations = P(dz, "max_iterations", "implicit.discrete", Dimension::Ordinal, false, 1);
      }
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
    if (const Node* cn = root.get("conditioning")) {
      auto& cd = d.conditioning;
      cd.present = true;
      if (const Node* data = cn->get("data"))
        for (std::size_t k = 0; k < data->seq.size(); ++k) {
          const Node& x = data->seq[k];
          cd.data.push_back({S(&x, "input"), S(&x, "kind"), S(&x, "format"), S(&x, "horizon"),
                             P(&x, "uncertainty", "conditioning.data[" + std::to_string(k) + "]",
                               S(&x, "kind") == "orientations" ? Dimension::Angle
                               : S(&x, "kind") == "gravity"    ? Dimension::Acceleration
                               : S(&x, "kind") == "magnetics"  ? Dimension::MagneticField
                                                               : Dimension::Length,
                               false, 1)});
        }
      if (const Node* u = cn->get("use_well_tops")) cd.use_well_tops = u->boolean || u->str == "true";
      cd.uncertainty = P(cn, "uncertainty", "conditioning", Dimension::Length, false, 1);
      cd.orientation_uncertainty = P(cn, "orientation_uncertainty", "conditioning", Dimension::Angle, false, 1);
      cd.max_iterations = P(cn, "max_iterations", "conditioning", Dimension::Ordinal, false, 1);
      if (cn->get("method")) cd.method = S(cn, "method");
      if (const Node* sm = cn->get("sampler")) {
        cd.walkers = P(sm, "walkers", "conditioning.sampler", Dimension::Ordinal, false, 1);
        cd.steps = P(sm, "steps", "conditioning.sampler", Dimension::Ordinal, false, 1);
        cd.burn_in = P(sm, "burn_in", "conditioning.sampler", Dimension::Ordinal, false, 1);
        cd.thin = P(sm, "thin", "conditioning.sampler", Dimension::Ordinal, false, 1);
        cd.seed = P(sm, "seed", "conditioning.sampler", Dimension::Ordinal, false, 1);
      }
      if (const Node* ps = cn->get("parameters"))
        for (std::size_t k = 0; k < ps->seq.size(); ++k) {
          const Node& x = ps->seq[k];
          std::string b = "conditioning.parameters[" + std::to_string(k) + "]";
          cd.parameters.push_back({S(&x, "path"), P(&x, "min", b, Dimension::Any, true, 1), P(&x, "max", b, Dimension::Any, true, 1)});
        }
    }
    if (const Node* rs = root.get("restoration")) {
      auto& R = d.restoration;
      R.present = true;
      R.bed_length_limit = P(rs, "bed_length_limit", "restoration", Dimension::Dimensionless, false, 1);
      if (const Node* v = rs->get("volumetric")) {
        const std::string b = "restoration.volumetric";
        R.volumetric = true;
        R.horizon = S(v, "horizon");
        R.datum = P(v, "datum", b, Dimension::Length, true, 1);
        R.pin = P(v, "pin", b, Dimension::Length, true, 2);
        R.young = P(v, "young", b, Dimension::Pressure, true, 1);
        R.poisson = P(v, "poisson", b, Dimension::Dimensionless, true, 1);
        R.mesh_size = P(v, "mesh_size", b, Dimension::Length, false, 1);
        R.max_strain = P(v, "max_strain", b, Dimension::Dimensionless, false, 1);
        if (const Node* us = v->get("units"))
          for (std::size_t k = 0; k < us->seq.size(); ++k) {
            const Node& x = us->seq[k];
            const std::string ub = b + ".units[" + std::to_string(k) + "]";
            R.units.push_back({S(&x, "unit"), P(&x, "young", ub, Dimension::Pressure, true, 1), P(&x, "poisson", ub, Dimension::Dimensionless, true, 1)});
          }
      }
    }
    if (const Node* bs2 = root.get("backstripping")) {
      auto& B = d.backstripping;
      const std::string b = "backstripping";
      B.present = true;
      B.basement = S(bs2, "basement");
      B.mantle_density = P(bs2, "mantle_density", b, Dimension::Density, true, 1);
      B.water_density = P(bs2, "water_density", b, Dimension::Density, true, 1);
      B.grain_density = P(bs2, "grain_density", b, Dimension::Density, true, 1);
      B.te = P(bs2, "te", b, Dimension::Length, true, 1);
      B.young = P(bs2, "young", b, Dimension::Pressure, false, 1);
      B.poisson = P(bs2, "poisson", b, Dimension::Dimensionless, false, 1);
      if (const Node* ss = bs2->get("steps"))
        for (std::size_t k = 0; k < ss->seq.size(); ++k) {
          const Node& x = ss->seq[k];
          const std::string sb = b + ".steps[" + std::to_string(k) + "]";
          B.steps.push_back({S(&x, "horizon"), P(&x, "water_depth", sb, Dimension::Length, true, 1), P(&x, "sea_level", sb, Dimension::Length, false, 1)});
        }
    }
    if (const Node* es = root.get("eroded_section"))
      for (std::size_t k = 0; k < es->seq.size(); ++k) {
        const Node& x = es->seq[k];
        const std::string b = "eroded_section[" + std::to_string(k) + "]";
        d.eroded_sections.push_back({S(&x, "name"), S(&x, "well"), S(&x, "unconformity"), S(&x, "log"),
                                     P(&x, "surface_porosity", b, Dimension::Dimensionless, true, 1), P(&x, "decay", b, Dimension::InverseLength, true, 1),
                                     P(&x, "vitrinite_estimate", b, Dimension::Length, false, 1)});
      }
    if (const Node* mc = root.get("mechanics")) {
      auto& M = d.mechanics;
      M.present = true;
      if (const Node* us = mc->get("units"))
        for (std::size_t k = 0; k < us->seq.size(); ++k) {
          const Node& x = us->seq[k];
          const std::string b = "mechanics.units[" + std::to_string(k) + "]";
          M.units.push_back({S(&x, "unit"), P(&x, "young", b, Dimension::Pressure, true, 1), P(&x, "poisson", b, Dimension::Dimensionless, true, 1),
                             P(&x, "ucs", b, Dimension::Pressure, false, 1), P(&x, "biot", b, Dimension::Dimensionless, false, 1)});
        }
      if (const Node* ts = mc->get("tectonic_strain")) {
        M.strain_min = P(ts, "min", "mechanics.tectonic_strain", Dimension::Dimensionless, true, 1);
        M.strain_max = P(ts, "max", "mechanics.tectonic_strain", Dimension::Dimensionless, true, 1);
      }
    }
    if (const Node* fr = root.get("fractures")) {
      auto& F = d.fractures;
      F.present = true;
      if (const Node* us = fr->get("units"))
        for (std::size_t k = 0; k < us->seq.size(); ++k) {
          const Node& x = us->seq[k];
          const std::string b = "fractures.units[" + std::to_string(k) + "]";
          F.units.push_back({S(&x, "unit"), P(&x, "p32_0", b, Dimension::InverseLength, true, 1), P(&x, "a_k", b, Dimension::Dimensionless, true, 1),
                             P(&x, "a_e", b, Dimension::InverseLength, false, 1), P(&x, "p32_max", b, Dimension::InverseLength, true, 1)});
        }
      if (const Node* df = fr->get("dfn")) {
        auto& N = F.dfn;
        const std::string b = "fractures.dfn";
        N.present = true;
        N.unit = S(df, "unit");
        N.area = P(df, "area", b, Dimension::Length, true, 2, 2);
        N.seed = P(df, "seed", b, Dimension::Ordinal, true, 1);
        N.limit = P(df, "limit", b, Dimension::Ordinal, false, 1);
        if (const Node* ss = df->get("sets"))
          for (std::size_t k = 0; k < ss->seq.size(); ++k) {
            const Node& x = ss->seq[k];
            const std::string sb = b + ".sets[" + std::to_string(k) + "]";
            N.sets.push_back({S(&x, "name"), P(&x, "pole", sb, Dimension::Angle, true, 2), P(&x, "kappa", sb, Dimension::Dimensionless, true, 1),
                              P(&x, "exponent", sb, Dimension::Dimensionless, true, 1), P(&x, "rmin", sb, Dimension::Length, true, 1),
                              P(&x, "rmax", sb, Dimension::Length, true, 1), P(&x, "weight", sb, Dimension::Dimensionless, true, 1)});
          }
      }
    }
    if (const Node* ht = root.get("heat")) {
      auto& T = d.heat;
      const std::string b = "heat";
      T.present = true;
      T.surface_temperature = P(ht, "surface_temperature", b, Dimension::Temperature, true, 1);
      T.basal_heat_flow = P(ht, "basal_heat_flow", b, Dimension::HeatFlow, true, 1);
      T.water_conductivity = P(ht, "water_conductivity", b, Dimension::ThermalConductivity, true, 1);
      T.column_depth = P(ht, "column_depth", b, Dimension::Length, false, 1);
      if (const Node* v = ht->get("vitrinite")) T.vitrinite = v->type == Node::Type::Bool && v->boolean;
      if (const Node* hs = ht->get("history"))
        for (std::size_t k = 0; k < hs->seq.size(); ++k) {
          const Node& x = hs->seq[k];
          const std::string hb = b + ".history[" + std::to_string(k) + "]";
          T.history.push_back({P(&x, "age", hb, Dimension::Time, true, 1), P(&x, "basal_heat_flow", hb, Dimension::HeatFlow, true, 1)});
        }
      if (const Node* us = ht->get("units"))
        for (std::size_t k = 0; k < us->seq.size(); ++k) {
          const Node& x = us->seq[k];
          const std::string ub = b + ".units[" + std::to_string(k) + "]";
          T.units.push_back({S(&x, "unit"), P(&x, "conductivity", ub, Dimension::ThermalConductivity, true, 1),
                             P(&x, "heat_production", ub, Dimension::HeatProduction, true, 1)});
        }
    }
    if (const Node* pr = root.get("pressure")) {
      auto& Q = d.pressure;
      const std::string b = "pressure";
      Q.present = true;
      Q.viscosity = P(pr, "viscosity", b, Dimension::Viscosity, true, 1);
      Q.rho_water = P(pr, "water_density", b, Dimension::Density, true, 1);
      Q.rho_grain = P(pr, "grain_density", b, Dimension::Density, true, 1);
      if (const Node* us = pr->get("units"))
        for (std::size_t k = 0; k < us->seq.size(); ++k) {
          const Node& x = us->seq[k];
          const std::string ub = b + ".units[" + std::to_string(k) + "]";
          PressureUnitDesc u;
          u.unit = S(&x, "unit");
          u.law = S(&x, "law");
          u.k0 = P(&x, "k0", ub, Dimension::Permeability, true, 1);
          u.b = P(&x, "b", ub, Dimension::Dimensionless, u.law == "exponential", 1);
          u.phi0 = P(&x, "phi0", ub, Dimension::Dimensionless, true, 1);
          Q.units.push_back(std::move(u));
        }
    }
    if (const Node* gs = root.get("geobodies"))
      for (std::size_t k = 0; k < gs->seq.size(); ++k) {
        const Node& x = gs->seq[k];
        const std::string b = "geobodies[" + std::to_string(k) + "]";
        GeobodyDesc g;
        g.name = S(&x, "name");
        g.type = S(&x, "type");
        g.unit = S(&x, "unit");
        g.color = S(&x, "color");
        g.shape = S(&x, "shape", "rectangular");
        const bool ch = g.type == "channels", lb = g.type == "lobes";
        g.count = P(&x, "count", b, Dimension::Ordinal, true, 1);
        g.seed = P(&x, "seed", b, Dimension::Ordinal, true, 1);
        g.thickness = P(&x, "thickness", b, Dimension::Length, true, 1);
        g.orientation = P(&x, "orientation", b, Dimension::Angle, true, 1);
        g.width = P(&x, "width", b, Dimension::Length, ch, 1);
        g.wavelength = P(&x, "wavelength", b, Dimension::Length, false, 1);
        g.amplitude = P(&x, "amplitude", b, Dimension::Length, false, 1);
        g.radius = P(&x, "radius", b, Dimension::Length, lb, 1);
        g.aspect = P(&x, "aspect", b, Dimension::Dimensionless, lb, 1);
        d.geobodies.push_back(std::move(g));
      }
    if (const Node* bs = root.get("bodies"))
      for (std::size_t k = 0; k < bs->seq.size(); ++k) {
        const Node& x = bs->seq[k];
        const std::string b = "bodies[" + std::to_string(k) + "]";
        BodyDesc bd;
        bd.name = S(&x, "name");
        bd.type = S(&x, "type");
        bd.color = S(&x, "color");
        bd.age = P(&x, "age", b, Dimension::Ordinal, true, 1);
        const bool sheet = bd.type == "salt_sheet", sill = bd.type == "sill" || bd.type == "saucer_sill", dyke = bd.type == "dyke", inj = bd.type == "injectite";
        if (sheet) {
          bd.footprint = P(&x, "footprint", b, Dimension::Length, true, 2, 2);
          if (const Node* bn = x.get("base")) {
            bd.base_horizon = S(bn, "horizon");
            bd.base_offset = P(bn, "offset", b + ".base", Dimension::Length, false, 1);
            bd.base_depth = P(bn, "depth", b + ".base", Dimension::Length, false, 1);
          }
          bd.thickness = P(&x, "thickness", b, Dimension::Length, true, 1);
          bd.taper = P(&x, "taper", b, Dimension::Length, false, 1);
        }
        if (sill) {
          bd.centre = P(&x, "centre", b, Dimension::Length, true, 2);
          bd.depth = P(&x, "depth", b, Dimension::Length, true, 1);
          bd.radius = P(&x, "radius", b, Dimension::Length, true, 1);
          bd.thickness = P(&x, "thickness", b, Dimension::Length, true, 1);
          bd.accommodation = P(&x, "accommodation", b, Dimension::Dimensionless, false, 1);
          bd.uplift_decay = P(&x, "uplift_decay", b, Dimension::Dimensionless, false, 1);
          bd.climb = P(&x, "climb", b, Dimension::Angle, bd.type == "saucer_sill", 1);
          bd.outer_radius = P(&x, "outer_radius", b, Dimension::Length, bd.type == "saucer_sill", 1);
        }
        if (dyke) {
          bd.trace = P(&x, "trace", b, Dimension::Length, true, 2, 2);
          bd.width = P(&x, "width", b, Dimension::Length, true, 1);
          bd.top_depth = P(&x, "top_depth", b, Dimension::Length, true, 1);
          bd.base_depth = P(&x, "base_depth", b, Dimension::Length, true, 1);
        }
        if (inj) {
          bd.parent = S(&x, "parent");
          bd.ceiling = S(&x, "ceiling");
          bd.footprint = P(&x, "area", b, Dimension::Length, true, 2, 2);
          bd.sheets = P(&x, "sheets", b, Dimension::Ordinal, true, 1);
          bd.length = P(&x, "length", b, Dimension::Length, true, 1);
          bd.length_sd = P(&x, "length_sd", b, Dimension::Length, false, 1);
          bd.thickness = P(&x, "thickness", b, Dimension::Length, true, 1);
          bd.thickness_sd = P(&x, "thickness_sd", b, Dimension::Length, false, 1);
          bd.climb = P(&x, "climb", b, Dimension::Angle, true, 1);
          bd.net_to_gross = P(&x, "net_to_gross", b, Dimension::Dimensionless, true, 1);
          bd.max_fraction = P(&x, "max_fraction", b, Dimension::Dimensionless, false, 1);
          bd.seed = P(&x, "seed", b, Dimension::Ordinal, true, 1);
        }
        d.bodies.push_back(bd);
      }
    if (const Node* gp = root.get("geophysics")) {
      auto& g = d.geophysics;
      g.present = true;
      g.cell = P(gp, "cell", "geophysics", Dimension::Length, false, 1);
      g.background_density = P(gp, "background_density", "geophysics", Dimension::Density, false, 1);
      g.background_susceptibility = P(gp, "background_susceptibility", "geophysics", Dimension::Dimensionless, false, 1);
      if (const Node* f = gp->get("inducing_field")) {
        g.inclination = P(f, "inclination", "geophysics.inducing_field", Dimension::Angle, false, 1);
        g.declination = P(f, "declination", "geophysics.inducing_field", Dimension::Angle, false, 1);
        g.intensity = P(f, "intensity", "geophysics.inducing_field", Dimension::MagneticField, false, 1);
      }
      if (const Node* us = gp->get("units"))
        for (std::size_t k = 0; k < us->seq.size(); ++k) {
          const Node& x = us->seq[k];
          const std::string b = "geophysics.units[" + std::to_string(k) + "]";
          g.units.push_back({S(&x, "unit"), P(&x, "density", b, Dimension::Density, false, 1), P(&x, "susceptibility", b, Dimension::Dimensionless, false, 1)});
        }
      if (const Node* gr = gp->get("gravity"))
        g.gravity = {S(gr, "input"), S(gr, "format"), P(gr, "uncertainty", "geophysics.gravity", Dimension::Acceleration, false, 1)};
      if (const Node* mg = gp->get("magnetics"))
        g.magnetics = {S(mg, "input"), S(mg, "format"), P(mg, "uncertainty", "geophysics.magnetics", Dimension::MagneticField, false, 1)};
    }
    if (const Node* un = root.get("uncertainty")) {
      auto& u = d.uncertainty;
      u.present = true;
      u.members = P(un, "members", "uncertainty", Dimension::Ordinal, true, 1);
      u.seed = P(un, "seed", "uncertainty", Dimension::Ordinal, true, 1);
      if (const Node* ps = un->get("parameters"))
        for (std::size_t k = 0; k < ps->seq.size(); ++k) {
          const Node& x = ps->seq[k];
          const std::string b = "uncertainty.parameters[" + std::to_string(k) + "]";
          UncertainParamDesc q;
          q.path = S(&x, "path");
          q.distribution = S(&x, "distribution");
          q.mean = P(&x, "mean", b, Dimension::Any, false, 1);
          q.sd = P(&x, "sd", b, Dimension::Any, false, 1);
          q.min = P(&x, "min", b, Dimension::Any, false, 1);
          q.max = P(&x, "max", b, Dimension::Any, false, 1);
          q.mode = P(&x, "mode", b, Dimension::Any, false, 1);
          q.values = P(&x, "values", b, Dimension::Any, false);
          q.weights = P(&x, "weights", b, Dimension::Dimensionless, false);
          u.parameters.push_back(std::move(q));
        }
      if (const Node* cs = un->get("correlations"))
        for (std::size_t k = 0; k < cs->seq.size(); ++k) {
          const Node& x = cs->seq[k];
          u.correlations.push_back({S(&x, "a"), S(&x, "b"), P(&x, "rho", "uncertainty.correlations[" + std::to_string(k) + "]", Dimension::Dimensionless, true, 1)});
        }
      if (const Node* ss = un->get("scenarios"))
        for (std::size_t k = 0; k < ss->seq.size(); ++k) {
          const Node& x = ss->seq[k];
          ScenarioDesc sc;
          sc.name = S(&x, "name");
          if (const Node* os = x.get("options"))
            for (std::size_t q = 0; q < os->seq.size(); ++q) {
              const Node& o = os->seq[q];
              ScenarioOptionDesc od;
              od.name = S(&o, "name");
              od.weight = P(&o, "weight", "uncertainty.scenarios[" + std::to_string(k) + "].options[" + std::to_string(q) + "]", Dimension::Dimensionless, true, 1);
              if (const Node* st = o.get("set"))
                for (auto& [key, val] : st->map) od.set.emplace_back(key, val);
              if (const Node* rm = o.get("remove"))
                for (auto& r : rm->seq) od.remove.push_back(r.str);
              sc.options.push_back(std::move(od));
            }
          u.scenarios.push_back(std::move(sc));
        }
      if (const Node* oa = un->get("one_at_a_time")) u.one_at_a_time = oa->boolean || oa->str == "true";
      if (const Node* ms = un->get("probability_maps"))
        for (std::size_t k = 0; k < ms->seq.size(); ++k) {
          const Node& x = ms->seq[k];
          const std::string b = "uncertainty.probability_maps[" + std::to_string(k) + "]";
          u.maps.push_back({S(&x, "kind"), S(&x, "unit"), S(&x, "horizon"), P(&x, "depth", b, Dimension::Length, false, 1),
                            P(&x, "spacing", b, Dimension::Length, false, 1)});
        }
      if (const Node* kn = un->get("key_numbers"))
        for (auto& x : kn->seq) u.key_numbers.push_back(x.str);
      if (const Node* se = un->get("sensitivity")) {
        auto& q = u.sensitivity;
        q.present = true;
        q.method = S(se, "method");
        q.base_samples = P(se, "base_samples", "uncertainty.sensitivity", Dimension::Ordinal, true, 1);
        q.bootstrap = P(se, "bootstrap", "uncertainty.sensitivity", Dimension::Ordinal, true, 1);
        q.confidence = P(se, "confidence", "uncertainty.sensitivity", Dimension::Dimensionless, true, 1);
      }
      if (const Node* in = un->get("interpretation")) {
        if (const Node* ps = in->get("picks"))
          for (std::size_t k = 0; k < ps->seq.size(); ++k) {
            const Node& x = ps->seq[k];
            const std::string b = "uncertainty.interpretation.picks[" + std::to_string(k) + "]";
            u.picks.push_back({S(&x, "horizon"), P(&x, "sigma", b, Dimension::TravelTime, true, 1), P(&x, "correlation_length", b, Dimension::Length, false, 1)});
          }
        if (const Node* vs = in->get("velocity"))
          for (std::size_t k = 0; k < vs->seq.size(); ++k) {
            const Node& x = vs->seq[k];
            const std::string b = "uncertainty.interpretation.velocity[" + std::to_string(k) + "]";
            u.velocity.push_back({S(&x, "layer"), P(&x, "v0_sigma", b, Dimension::Velocity, false, 1), P(&x, "k_sigma", b, Dimension::VelocityGradient, false, 1),
                                  P(&x, "correlation_length", b, Dimension::Length, false, 1)});
          }
      }
    }
    if (const Node* gn = root.get("grid")) {
      auto& g = d.grid;
      g.present = true;
      g.spacing = P(gn, "spacing", "grid", Dimension::Length, true);
      g.margin = P(gn, "margin", "grid", Dimension::Length, false, 1);
      g.min_thickness = P(gn, "min_thickness", "grid", Dimension::Length, false, 1);
      if (const Node* zs = gn->get("zones"))
        for (std::size_t k = 0; k < zs->seq.size(); ++k) {
          const Node& x = zs->seq[k];
          std::string b = "grid.zones[" + std::to_string(k) + "]";
          GridZoneDesc z;
          z.unit = S(&x, "unit");
          if (!S(&x, "layering").empty()) z.layering = S(&x, "layering");
          z.layers = P(&x, "layers", b, Dimension::Ordinal, true, 1);
          z.layer_thickness = P(&x, "layer_thickness", b, Dimension::Length, false, 1);
          g.zones.push_back(std::move(z));
        }
      if (const Node* fs = gn->get("facies"))
        for (std::size_t k = 0; k < fs->seq.size(); ++k) {
          const Node& x = fs->seq[k];
          g.facies.push_back({S(&x, "unit"), S(&x, "name"), P(&x, "code", "grid.facies[" + std::to_string(k) + "]", Dimension::Ordinal, true, 1)});
        }
      g.porosity = S(gn, "porosity");
      if (const Node* ts = gn->get("transforms"))
        for (std::size_t k = 0; k < ts->seq.size(); ++k) {
          const Node& x = ts->seq[k];
          std::string b = "grid.transforms[" + std::to_string(k) + "]";
          g.transforms.push_back({S(&x, "name"), S(&x, "unit"), S(&x, "from"), S(&x, "law"), P(&x, "a", b, Dimension::Any, true, 1),
                                  P(&x, "b", b, Dimension::Any, true, 1)});
        }
      if (const Node* ft = gn->get("fault_transmissibility")) g.transmissibility = S(ft, "calibration");
      if (const Node* tt = gn->get("tetrahedral")) g.tetrahedral = tt->boolean || tt->str == "true";
      g.pillars = S(gn, "pillars", "fault_conforming");
      if (const Node* gs = gn->get("geostatistics"))
        for (std::size_t k = 0; k < gs->seq.size(); ++k) {
          const Node& x = gs->seq[k];
          const std::string b = "grid.geostatistics[" + std::to_string(k) + "]";
          GridDesc::Geostat q;
          q.name = S(&x, "name");
          q.kind = S(&x, "kind");
          q.unit = S(&x, "unit");
          q.property_unit = S(&x, "property_unit");
          q.condition_to = S(&x, "condition_to");
          q.trend = S(&x, "trend");
          auto vario = [&](const Node* v, const std::string& vb) {
            GridDesc::VariogramDesc vd;
            if (!v) return vd;
            vd.present = true;
            vd.nugget = P(v, "nugget", vb, Dimension::Dimensionless, false, 1);
            if (const Node* fv = v->get("fit_vertical")) vd.fit_vertical = fv->boolean || fv->str == "true";
            if (const Node* ss = v->get("structures"))
              for (std::size_t i = 0; i < ss->seq.size(); ++i) {
                const Node& st = ss->seq[i];
                const std::string sb = vb + ".structures[" + std::to_string(i) + "]";
                vd.structures.push_back({S(&st, "type"), P(&st, "contribution", sb, Dimension::Dimensionless, true, 1), P(&st, "ranges", sb, Dimension::Length, true, 3),
                                         P(&st, "azimuth", sb, Dimension::Angle, false, 1)});
              }
            return vd;
          };
          q.variogram = vario(x.get("variogram"), b + ".variogram");
          if (const Node* cs = x.get("codes"))
            for (std::size_t i = 0; i < cs->seq.size(); ++i) {
              const Node& c = cs->seq[i];
              const std::string cb = b + ".codes[" + std::to_string(i) + "]";
              GridDesc::Geostat::Code cd;
              if (const Node* cn = c.get("code")) cd.code = static_cast<int>(cn->is_number() ? cn->number : std::stod(cn->str));
              cd.name = S(&c, "name");
              cd.proportion = P(&c, "proportion", cb, Dimension::Dimensionless, true, 1);
              cd.variogram = vario(c.get("variogram"), cb + ".variogram");
              q.codes.push_back(std::move(cd));
            }
          if (const Node* ds = x.get("distribution")) {
            const std::string db = b + ".distribution";
            q.min = P(ds, "min", db, Dimension::Any, true, 1);
            q.max = P(ds, "max", db, Dimension::Any, true, 1);
            q.mean = P(ds, "mean", db, Dimension::Any, false, 1);
            q.sd = P(ds, "sd", db, Dimension::Any, false, 1);
            q.decluster_cell = P(ds, "decluster_cell", db, Dimension::Length, false, 1);
          }
          if (const Node* se = x.get("search")) {
            q.search_radii = P(se, "radii", b + ".search", Dimension::Length, true, 3);
            q.search_azimuth = P(se, "azimuth", b + ".search", Dimension::Angle, false, 1);
            q.max_data = P(se, "max_data", b + ".search", Dimension::Ordinal, true, 1);
          }
          if (const Node* ob = x.get("object")) {
            const std::string obb = b + ".object";
            q.object_kind = S(ob, "kind");
            q.geobody = S(ob, "geobody");
            if (const Node* oc = ob->get("code")) q.object_code = static_cast<int>(oc->is_number() ? oc->number : std::stod(oc->str));
            q.width = P(ob, "width", obb, Dimension::Length, false, 1);
            q.thickness = P(ob, "thickness", obb, Dimension::Length, true, 1);
            q.sinuosity = P(ob, "sinuosity", obb, Dimension::Dimensionless, false, 1);
            q.wavelength = P(ob, "wavelength", obb, Dimension::Length, false, 1);
            q.azimuth = P(ob, "azimuth", obb, Dimension::Angle, false, 1);
            q.radius = P(ob, "radius", obb, Dimension::Length, false, 1);
            q.proportion = P(ob, "proportion", obb, Dimension::Dimensionless, true, 1);
            q.max_iterations = P(ob, "max_iterations", obb, Dimension::Ordinal, false, 1);
          }
          q.seed = P(&x, "seed", b, Dimension::Ordinal, true, 1);
          g.geostatistics.push_back(std::move(q));
        }
      if (const Node* ul = gn->get("upscaled_logs"))
        for (std::size_t k = 0; k < ul->seq.size(); ++k) {
          const Node& x = ul->seq[k];
          g.upscaled_logs.push_back({S(&x, "curve"), S(&x, "method"), S(&x, "name"),
                                     P(&x, "min_valid", "grid.upscaled_logs[" + std::to_string(k) + "]", Dimension::Dimensionless, false, 1)});
        }
    }
    if (const Node* wc = root.get("well_correlation")) d.correlation_rule = S(wc, "rule");
    if (const Node* vm = root.get("velocity_model")) {
      auto& V = d.velocity;
      V.present = true;
      if (const Node* ls = vm->get("layers"))
        for (std::size_t k = 0; k < ls->seq.size(); ++k) {
          const Node& l = ls->seq[k];
          const std::string b = "velocity_model.layers[" + std::to_string(k) + "]";
          const Node* fk = l.get("fit_k");
          V.layers.push_back({S(&l, "name"), S(&l, "base"), P(&l, "v0", b, Dimension::Velocity, true, 1), P(&l, "k", b, Dimension::VelocityGradient, false, 1),
                              fk && (fk->boolean || fk->str == "true"), {}});
          if (const Node* pp = l.get("perturbation")) {
            auto& q = V.layers.back().perturbation;
            q.present = true;
            q.v0_sigma = P(pp, "v0_sigma", b + ".perturbation", Dimension::Velocity, false, 1);
            q.k_sigma = P(pp, "k_sigma", b + ".perturbation", Dimension::VelocityGradient, false, 1);
            q.correlation_length = P(pp, "correlation_length", b + ".perturbation", Dimension::Length, false, 1);
            if (const Node* sd = pp->get("seed")) q.seed = static_cast<std::uint64_t>(sd->is_number() ? sd->number : std::stod(sd->str));
          }
        }
      if (const Node* ca = vm->get("calibration")) {
        V.interpolation = S(ca, "interpolation", "inverse_distance");
        if (const Node* ex = ca->get("exclude"))
          for (auto& e : ex->seq) V.exclude.push_back({S(&e, "well"), S(&e, "reason")});
      }
      V.tolerance = P(vm, "tolerance", "velocity_model", Dimension::Length, false, 1);
    }
    if (const Node* sm = root.get("seismic")) {
      d.seismic.present = true;
      d.seismic.input = S(sm, "input");
      d.seismic.domain = S(sm, "domain");
      if (const Node* lo = sm->get("layout")) {
        auto gi = [&](const char* k) { const Node* n = lo->get(k); return n ? static_cast<int>(n->number) : 0; };
        d.seismic.il = gi("inline_byte");
        d.seismic.xl = gi("crossline_byte");
        d.seismic.x = gi("x_byte");
        d.seismic.y = gi("y_byte");
      }
    }
    if (const Node* gm = root.get("geological_map")) {
      auto& G = d.geological_map;
      G.present = true;
      auto src = [&](const Node* n) { return Description::MapSource{S(n, "input"), S(n, "format"), S(n, "layer"), S(n, "field")}; };
      G.units = src(gm->get("units"));
      if (gm->get("faults")) G.faults = src(gm->get("faults"));
      if (const Node* ms = gm->get("measurements")) {
        G.measurements = src(ms);
        G.dip_field = S(ms, "dip_field");
        G.direction_field = S(ms, "direction_field");
        G.convention = S(ms, "convention");
        G.filter_field = S(ms, "filter_field");
        G.filter_value = S(ms, "filter_value");
        G.source = S(ms, "source");
        G.date = S(ms, "date");
      }
      if (const Node* uo = gm->get("unit_order"))
        for (auto& u : uo->seq) G.unit_order.push_back({S(&u, "map_unit"), S(&u, "base")});
      if (const Node* em = gm->get("elevation_model")) { G.dem_input = S(em, "input"); G.dem_format = S(em, "format"); }
      G.sample_spacing = P(gm, "sample_spacing", "geological_map", Dimension::Length, false, 1);
      // the map's contacts and measurements enter the implicit builder as data
      d.implicit.data.push_back({"geological_map:contacts", "points", "geological_map", "", Param{}});
      if (!G.measurements.input.empty()) d.implicit.data.push_back({"geological_map:orientations", "orientations", "geological_map", "", Param{}});
    }
    if (const Node* fr = root.get("fault_rock")) {
      d.fault_rock.present = true;
      d.fault_rock.shale_threshold = P(fr, "shale_threshold", "fault_rock", Dimension::Dimensionless, false, 1);
      d.fault_rock.thickness_ratio = P(fr, "thickness_displacement_ratio", "fault_rock", Dimension::Dimensionless, false, 1);
      if (const Node* sc = fr->get("seal_calibration")) {
        d.fault_rock.seal_type = S(sc, "type");
        d.fault_rock.seal_a = P(sc, "a", "fault_rock.seal_calibration", Dimension::Dimensionless, true, 1);
        d.fault_rock.seal_c = P(sc, "c", "fault_rock.seal_calibration", Dimension::Dimensionless, true, 1);
      }
    }
    if (const Node* sn = root.get("stress")) {
      d.stress.present = true;
      if (const Node* rg = sn->get("regime")) {
        d.stress.sv_gradient = P(rg, "sv_gradient", "stress.regime", Dimension::PressureGradient, true, 1);
        d.stress.shmax_ratio = P(rg, "shmax_ratio", "stress.regime", Dimension::Dimensionless, true, 1);
        d.stress.shmin_ratio = P(rg, "shmin_ratio", "stress.regime", Dimension::Dimensionless, true, 1);
        d.stress.shmax_azimuth = P(rg, "shmax_azimuth", "stress.regime", Dimension::Angle, true, 1);
      }
      d.stress.tensor_gradient = P(sn, "tensor_gradient", "stress", Dimension::PressureGradient, false, 6);
      d.stress.pore_pressure_gradient = P(sn, "pore_pressure_gradient", "stress", Dimension::PressureGradient, true, 1);
      d.stress.friction = P(sn, "friction", "stress", Dimension::Dimensionless, false, 1);
      d.stress.cohesion = P(sn, "cohesion", "stress", Dimension::Pressure, false, 1);
    }
    if (const Node* o2 = root.get("outputs"))
      if (const Node* cv = o2->get("curvature")) {
        if (const Node* hs2 = cv->get("horizons"))
          for (auto& x : hs2->seq) d.curvature_horizons.push_back(x.str);
        d.fracture_threshold = P(cv, "fracture_threshold", "outputs.curvature", Dimension::InverseLength, true, 1);
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

  /// Stage 4B: hanging walls that move level by level (growth strata and phases).
  static bool level_hw(const FaultDesc& fd) {
    if (fd.curved() || !fd.fold_type.empty()) return false;
    return fd.hw_type == "layer_rigid" || fd.hw_type == "rigid" || fd.hw_type == "decay" || fd.hw_type == "linear_between_faults";
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
          if (sr.clino_trajectory.present && !(sr.clino_trajectory.si() >= 0 && sr.clino_trajectory.si() <= deg2rad(15)))
            c.err(b + ".clinoforms.trajectory", "F23", "the shelf-break trajectory rises at 0 to 15 degrees (a falling trajectory would cut the "
                                                       "older topsets: not built)");
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
      if (d.compaction.burial != "block" && d.compaction.burial != "restored")
        c.err("stratigraphy.compaction.burial", "F25", "burial is block (maximum burial in its fault block) or restored (column-wise)");
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
    // Stage 4B: restoration, backstripping and eroded-section estimates (rule F40)
    if (d.restoration.volumetric) {
      const auto& R = d.restoration;
      const std::string b = "restoration.volumetric";
      if (d.horizon_index(R.horizon) < 0) c.err(b + ".horizon", "F40", "unknown horizon '" + R.horizon + "'");
      auto nu_ok = [&](const Param& p, const std::string& path) {
        if (p.present && !(p.si() >= 0 && p.si() < 0.5)) c.err(path, "F40", "Poisson's ratio from 0 to 0.5 (exclusive)");
      };
      auto e_ok = [&](const Param& p, const std::string& path) {
        if (p.present && !(p.si() > 0)) c.err(path, "F40", "Young's modulus must be positive");
      };
      nu_ok(R.poisson, b + ".poisson");
      e_ok(R.young, b + ".young");
      for (std::size_t k = 0; k < R.units.size(); ++k) {
        const std::string ub = b + ".units[" + std::to_string(k) + "]";
        if (d.unit_index(R.units[k].unit) < 0) c.err(ub + ".unit", "F40", "unknown unit '" + R.units[k].unit + "'");
        nu_ok(R.units[k].poisson, ub + ".poisson");
        e_ok(R.units[k].young, ub + ".young");
      }
      if (R.max_strain.present && !(R.max_strain.si() > 0 && R.max_strain.si() < 1)) c.err(b + ".max_strain", "F40", "the strain limit lies in (0, 1)");
      if (R.mesh_size.present && !(R.mesh_size.si() > 0)) c.err(b + ".mesh_size", "F40", "the mesh size must be positive");
    }
    if (d.restoration.bed_length_limit.present && !(d.restoration.bed_length_limit.si() > 0 && d.restoration.bed_length_limit.si() < 1))
      c.err("restoration.bed_length_limit", "F40", "the bed-length limit lies in (0, 1)");
    if (d.backstripping.present) {
      const auto& B = d.backstripping;
      const std::string b = "backstripping";
      if (!d.compaction.present()) c.err(b, "F40", "backstripping needs the compaction laws of stratigraphy.compaction");
      const int bh = d.horizon_index(B.basement);
      if (bh < 0) c.err(b + ".basement", "F40", "unknown basement horizon '" + B.basement + "'");
      if (B.te.present && !(B.te.si() >= 0 && B.te.si() <= 100e3)) c.err(b + ".te", "F40", "Te from 0 (Airy) to 100 km");
      if (B.te.present && B.te.si() > 0 && (!B.young.present || !B.poisson.present)) c.err(b, "F40", "flexure (Te > 0) needs young and poisson");
      if (B.mantle_density.present && B.water_density.present && B.grain_density.present &&
          !(B.mantle_density.si() > B.grain_density.si() && B.grain_density.si() > B.water_density.si()))
        c.err(b, "F40", "densities: mantle > grain > water");
      std::set<std::string> given;
      for (std::size_t k = 0; k < B.steps.size(); ++k) {
        const auto& s = B.steps[k];
        if (s.horizon != "present" && d.horizon_index(s.horizon) < 0)
          c.err(b + ".steps[" + std::to_string(k) + "].horizon", "F40", "unknown horizon '" + s.horizon + "' (or 'present')");
        given.insert(s.horizon);
      }
      // a paleo-water depth for every step (zero allowed, with its evidence class)
      for (int h = 0; bh >= 0 && h < bh; ++h)
        if (!given.count(d.horizons[static_cast<std::size_t>(h)].name))
          c.err(b + ".steps", "F40", "no water depth for the step when " + d.horizons[static_cast<std::size_t>(h)].name + " was at the surface");
    }
    for (std::size_t k = 0; k < d.eroded_sections.size(); ++k) {
      const auto& e = d.eroded_sections[k];
      const std::string b = "eroded_section[" + std::to_string(k) + "]";
      if (!d.well(e.well)) c.err(b + ".well", "F40", "unknown well '" + e.well + "'");
      if (d.horizon_index(e.unconformity) < 0) c.err(b + ".unconformity", "F40", "unknown horizon '" + e.unconformity + "'");
      if (std::none_of(d.inputs.begin(), d.inputs.end(), [&](const InputDesc& i) { return i.id == e.log; }))
        c.err(b + ".log", "F40", "unknown input '" + e.log + "' (a porosity log: depth, porosity)");
      if (e.surface_porosity.present && !(e.surface_porosity.si() > 0 && e.surface_porosity.si() < 1)) c.err(b + ".surface_porosity", "F40", "surface porosity in (0, 1)");
      if (e.decay.present && !(e.decay.si() > 0)) c.err(b + ".decay", "F40", "the decay must be positive");
    }
    // Stage 4B: mechanics and fractures (rule F42), heat and pressure (rule F41)
    {
      auto unit_ok = [&](const std::string& u, const std::string& path, const char* rule) {
        if (d.unit_index(u) < 0) c.err(path, rule, "unknown unit '" + u + "'");
      };
      for (std::size_t k = 0; k < d.mechanics.units.size(); ++k) {
        const auto& u = d.mechanics.units[k];
        const std::string b = "mechanics.units[" + std::to_string(k) + "]";
        unit_ok(u.unit, b + ".unit", "F42");
        if (u.young.present && !(u.young.si() > 0)) c.err(b + ".young", "F42", "Young's modulus must be positive");
        if (u.poisson.present && !(u.poisson.si() >= 0 && u.poisson.si() < 0.5)) c.err(b + ".poisson", "F42", "Poisson's ratio from 0 to 0.5 (exclusive)");
        if (u.biot.present && !(u.biot.si() >= 0 && u.biot.si() <= 1)) c.err(b + ".biot", "F42", "Biot's coefficient from 0 to 1");
      }
      for (const Param* p : {&d.mechanics.strain_min, &d.mechanics.strain_max})
        if (p->present && !(p->si() >= -0.005 && p->si() <= 0.005)) c.err("mechanics.tectonic_strain", "F42", "tectonic strains from -0.005 to 0.005");
      if (d.mechanics.strain_min.present && !d.stress.shmax_azimuth.present)
        c.err("mechanics.tectonic_strain", "F42", "tectonic strain needs stress.regime.shmax_azimuth (the direction of the larger strain)");
      for (std::size_t k = 0; k < d.fractures.units.size(); ++k) {
        const auto& u = d.fractures.units[k];
        const std::string b = "fractures.units[" + std::to_string(k) + "]";
        unit_ok(u.unit, b + ".unit", "F42");
        for (const Param* p : {&u.p32_0, &u.a_k, &u.a_e})
          if (p->present && p->si() < 0) c.err(b, "F42", "fracture-law coefficients must not be negative");
        if (u.p32_max.present && !(u.p32_max.si() > 0 && u.p32_max.si() <= 10)) c.err(b + ".p32_max", "F42", "the maximum P32 lies in (0, 10] 1/m");
      }
      if (d.fractures.dfn.present) {
        const auto& N = d.fractures.dfn;
        unit_ok(N.unit, "fractures.dfn.unit", "F42");
        if (std::none_of(d.fractures.units.begin(), d.fractures.units.end(), [&](const FractureUnitDesc& u) { return u.unit == N.unit; }))
          c.err("fractures.dfn.unit", "F42", "the network's unit needs a fracture law under fractures.units");
        if (N.sets.empty()) c.err("fractures.dfn.sets", "F42", "a network needs at least one set");
        for (std::size_t k = 0; k < N.sets.size(); ++k) {
          const auto& s = N.sets[k];
          const std::string b = "fractures.dfn.sets[" + std::to_string(k) + "]";
          if (s.kappa.present && s.kappa.si() < 0) c.err(b + ".kappa", "F42", "Fisher kappa must not be negative");
          if (s.exponent.present && !(s.exponent.si() > 2)) c.err(b + ".exponent", "F42", "the power-law exponent must exceed 2");
          if (s.rmin.present && s.rmax.present && !(s.rmin.si() > 0 && s.rmin.si() < s.rmax.si())) c.err(b, "F42", "0 < rmin < rmax");
          if (s.weight.present && !(s.weight.si() > 0)) c.err(b + ".weight", "F42", "a set's weight must be positive");
        }
      }
      if (d.heat.present) {
        const auto& T = d.heat;
        auto q_ok = [&](const Param& q, const std::string& path) {
          if (q.present && !(q.si() >= 0.020 && q.si() <= 0.200)) c.err(path, "F41", "heat flow from 20 to 200 mW/m2");
        };
        q_ok(T.basal_heat_flow, "heat.basal_heat_flow");
        for (std::size_t k = 0; k < T.history.size(); ++k) q_ok(T.history[k].second, "heat.history[" + std::to_string(k) + "].basal_heat_flow");
        if (!d.compaction.present()) c.err("heat", "F41", "heat needs the burial history of stratigraphy.compaction");
        for (std::size_t k = 0; k < T.units.size(); ++k) {
          const auto& u = T.units[k];
          const std::string b = "heat.units[" + std::to_string(k) + "]";
          unit_ok(u.unit, b + ".unit", "F41");
          if (u.conductivity.present && !(u.conductivity.si() >= 0.5 && u.conductivity.si() <= 6)) c.err(b + ".conductivity", "F41", "conductivity from 0.5 to 6 W/m/K");
          if (u.heat_production.present && !(u.heat_production.si() >= 0 && u.heat_production.si() <= 10e-6))
            c.err(b + ".heat_production", "F41", "radiogenic heat from 0 to 10 uW/m3");
        }
        // every compacting unit needs its thermal properties (no hidden defaults)
        for (auto& cu : d.compaction.units)
          if (std::none_of(T.units.begin(), T.units.end(), [&](const HeatUnitDesc& u) { return u.unit == cu.unit; }))
            c.err("heat.units", "F41", "unit '" + cu.unit + "' has a compaction law but no conductivity and heat production");
      }
      if (d.pressure.present) {
        const auto& Q = d.pressure;
        if (!d.compaction.present()) c.err("pressure", "F41", "pressure needs the compaction laws of stratigraphy.compaction");
        for (std::size_t k = 0; k < Q.units.size(); ++k) {
          const auto& u = Q.units[k];
          const std::string b = "pressure.units[" + std::to_string(k) + "]";
          unit_ok(u.unit, b + ".unit", "F41");
          if (u.law != "exponential" && u.law != "kozeny_carman") c.err(b + ".law", "F41", "law must be exponential or kozeny_carman");
          if (u.k0.present && !(u.k0.si() >= 1e-23 && u.k0.si() <= 1e-12)) c.err(b + ".k0", "F41", "permeability from 1e-23 to 1e-12 m2 (1e-8 to 1e3 mD)");
          if (u.phi0.present && !(u.phi0.si() > 0 && u.phi0.si() < 1)) c.err(b + ".phi0", "F41", "reference porosity in (0, 1)");
        }
        for (auto& cu : d.compaction.units)
          if (std::none_of(Q.units.begin(), Q.units.end(), [&](const PressureUnitDesc& u) { return u.unit == cu.unit; }))
            c.err("pressure.units", "F41", "unit '" + cu.unit + "' has a compaction law but no permeability law");
      }
    }
    // Stage 4B: geobodies (rule F39)
    {
      std::set<std::string> gnames;
      for (std::size_t k = 0; k < d.geobodies.size(); ++k) {
        const auto& g = d.geobodies[k];
        const std::string b = "geobodies[" + std::to_string(k) + "]";
        if (g.name.empty() || !gnames.insert(g.name).second || d.unit_index(g.name) >= 0)
          c.err(b + ".name", "F39", "every geobody set needs a unique name that no unit has");
        if (g.type != "channels" && g.type != "lobes") c.err(b + ".type", "F39", "type must be channels or lobes (got '" + g.type + "')");
        const int u = d.unit_index(g.unit);
        if (u < 0) c.err(b + ".unit", "F39", "unknown unit '" + g.unit + "'");
        else if (d.horizon_index(d.units[static_cast<std::size_t>(u)].top) < 0 || d.horizon_index(d.units[static_cast<std::size_t>(u)].base) < 0)
          c.err(b + ".unit", "F39", "geobodies lie in a unit between two horizons");
        if (g.shape != "rectangular" && g.shape != "parabolic") c.err(b + ".shape", "F39", "shape must be rectangular or parabolic");
        if (g.count.present && !(g.count.si() >= 1 && g.count.si() <= 500)) c.err(b + ".count", "F39", "from 1 to 500 geobodies per set");
        if (g.thickness.present && !(g.thickness.si() > 0)) c.err(b + ".thickness", "F39", "thickness must be positive");
        const double w = g.type == "lobes" && g.radius.present && g.aspect.present ? 2 * g.radius.si() / g.aspect.si() : (g.width.present ? g.width.si() : kNaN);
        if (std::isfinite(w) && g.thickness.present && g.thickness.si() > 0 && !(w / g.thickness.si() >= 5 && w / g.thickness.si() <= 500))
          c.err(b, "F39", "width-to-thickness ratio from 5 to 500 (got " + fmt_fixed(w / g.thickness.si(), 1) + ")");
        if (g.type == "lobes" && g.aspect.present && !(g.aspect.si() >= 1)) c.err(b + ".aspect", "F39", "aspect (major / minor axis) of at least 1");
        if (g.amplitude.present && g.amplitude.si() > 0 && !g.wavelength.present) c.err(b + ".wavelength", "F39", "a sinuous channel belt needs a wavelength");
      }
    }
    // STX-TN-15 (F38): bodies
    {
      std::set<std::string> names;
      for (std::size_t k = 0; k < d.bodies.size(); ++k) {
        const auto& bd = d.bodies[k];
        const std::string b = "bodies[" + std::to_string(k) + "]";
        if (bd.name.empty() || !names.insert(bd.name).second) c.err(b + ".name", "F38", "every body needs a unique name");
        if (std::any_of(d.units.begin(), d.units.end(), [&](const UnitDesc& u) { return u.name == bd.name; }))
          c.err(b + ".name", "F38", "a body cannot share its name with a unit");
        if (bd.type != "salt_sheet" && bd.type != "sill" && bd.type != "saucer_sill" && bd.type != "dyke" && bd.type != "injectite")
          c.err(b + ".type", "F38", "type must be salt_sheet, sill, saucer_sill, dyke or injectite (got '" + bd.type + "')");
        auto pos = [&](const Param& p, const std::string& key) {
          if (p.present && !(p.si() > 0)) c.err(b + "." + key, "F38", key + " must be positive");
        };
        pos(bd.thickness, "thickness");
        pos(bd.radius, "radius");
        pos(bd.width, "width");
        pos(bd.length, "length");
        if (bd.type == "salt_sheet") {
          if (bd.footprint.present && bd.footprint.rows() < 3) c.err(b + ".footprint", "F38", "a footprint needs at least three vertices");
          if (bd.base_horizon.empty() == !bd.base_depth.present) c.err(b + ".base", "F38", "give the sheet's base as a horizon (with an offset) or as a depth, not both");
          if (!bd.base_horizon.empty() && d.horizon_index(bd.base_horizon) < 0) c.err(b + ".base.horizon", "F38", "unknown horizon '" + bd.base_horizon + "'");
          if (bd.taper.present && !(bd.taper.si() >= 0)) c.err(b + ".taper", "F38", "the taper cannot be negative");
        }
        if (bd.type == "sill" || bd.type == "saucer_sill") {
          if (bd.thickness.present && (bd.thickness.si() < 1 || bd.thickness.si() > 300)) c.err(b + ".thickness", "F38", "sill thickness from 1 to 300 m");
          if (bd.accommodation.present && !(bd.accommodation.si() >= 0 && bd.accommodation.si() <= 1)) c.err(b + ".accommodation", "F38", "accommodation from 0 to 1");
          if (bd.uplift_decay.present && !(bd.uplift_decay.si() >= 0 && bd.uplift_decay.si() <= 1)) c.err(b + ".uplift_decay", "F38", "uplift decay from 0 to 1 per km");
          if (bd.type == "saucer_sill") {
            if (bd.climb.present && (bd.climb.si() < 5 * kPi / 180 || bd.climb.si() > 45 * kPi / 180)) c.err(b + ".climb", "F38", "saucer climb from 5 to 45 degrees");
            if (bd.outer_radius.present && bd.radius.present && !(bd.outer_radius.si() > bd.radius.si())) c.err(b + ".outer_radius", "F38", "the outer radius must exceed the inner sill's radius");
          }
        }
        if (bd.type == "dyke" && bd.top_depth.present && bd.base_depth.present && !(bd.base_depth.si() > bd.top_depth.si()))
          c.err(b + ".base_depth", "F38", "a dyke's base must be deeper than its top");
        if (bd.type == "injectite") {
          int pu = -1;
          for (std::size_t u = 0; u < d.units.size(); ++u)
            if (d.units[u].name == bd.parent) pu = static_cast<int>(u);
          if (pu < 0) c.err(b + ".parent", "F38", "unknown parent unit '" + bd.parent + "'");
          const int ch = d.horizon_index(bd.ceiling);
          if (ch < 0) c.err(b + ".ceiling", "F38", "unknown ceiling horizon '" + bd.ceiling + "'");
          if (pu >= 0 && ch >= 0 && !(ch < d.horizon_index(d.units[static_cast<std::size_t>(pu)].top)))
            c.err(b + ".ceiling", "F38", "the ceiling must be younger than (above) the parent unit's top");
          if (bd.max_fraction.present && !(bd.max_fraction.si() > 0 && bd.max_fraction.si() <= 1)) c.err(b + ".max_fraction", "F38", "the fraction must lie in (0, 1]");
          if (bd.net_to_gross.present && !(bd.net_to_gross.si() > 0 && bd.net_to_gross.si() <= 1)) c.err(b + ".net_to_gross", "F38", "net-to-gross must lie in (0, 1]");
          if (bd.footprint.present && bd.footprint.rows() < 3) c.err(b + ".area", "F38", "an area needs at least three vertices");
        }
      }
    }
    // STX-TN-14 (F37): gravity and magnetics
    {
      const auto& g = d.geophysics;
      bool cond_grav = false, cond_mag = false;
      for (auto& x : d.conditioning.data) {
        cond_grav = cond_grav || x.kind == "gravity";
        cond_mag = cond_mag || x.kind == "magnetics";
      }
      if ((cond_grav || cond_mag) && !g.present) c.err("conditioning.data", "F37", "gravity or magnetic data need a geophysics block (unit densities, susceptibilities)");
      if (g.present) {
        if (g.cell.present && !(g.cell.si() > 0)) c.err("geophysics.cell", "F37", "the prism size must be positive");
        bool any_density = false, any_susc = false;
        for (std::size_t k = 0; k < g.units.size(); ++k) {
          const auto& u = g.units[k];
          const std::string b = "geophysics.units[" + std::to_string(k) + "]";
          if (std::none_of(d.units.begin(), d.units.end(), [&](const UnitDesc& x) { return x.name == u.unit; }))
            c.err(b + ".unit", "F37", "unit '" + u.unit + "' is not a unit of the stratigraphy");
          if (u.density.present) {
            any_density = true;
            if (!(u.density.si() > 0)) c.err(b + ".density", "F37", "a density must be positive");
          }
          if (u.susceptibility.present) {
            any_susc = true;
            if (!(u.susceptibility.si() > -1)) c.err(b + ".susceptibility", "F37", "a susceptibility must exceed -1 SI");
          }
        }
        const bool grav = !g.gravity.input.empty() || cond_grav, mag = !g.magnetics.input.empty() || cond_mag;
        if (!grav && !mag) c.err("geophysics", "F37", "give gravity or magnetics stations (geophysics.gravity, geophysics.magnetics or conditioning data)");
        if (grav && !g.background_density.present) c.err("geophysics.background_density", "F37", "gravity needs the background density (anomalies are contrasts against it)");
        if (grav && !any_density) c.err("geophysics.units", "F37", "gravity needs the density of at least one unit");
        if (mag && !any_susc) c.err("geophysics.units", "F37", "magnetics needs the susceptibility of at least one unit");
        if (mag && !(g.inclination.present && g.declination.present && g.intensity.present))
          c.err("geophysics.inducing_field", "F37", "magnetics needs the inducing field: inclination, declination and intensity");
        if (g.background_density.present && !(g.background_density.si() > 0)) c.err("geophysics.background_density", "F37", "a density must be positive");
        for (auto* st : {&g.gravity, &g.magnetics})
          if (!st->input.empty() && std::none_of(d.inputs.begin(), d.inputs.end(), [&](const InputDesc& i) { return i.id == st->input; }))
            c.err("geophysics", "F37", "input '" + st->input + "' is not listed under inputs");
      }
    }
    // STX-TN-11 (F31): hybrid conditioning
    if (d.conditioning.present) {
      const auto& cd = d.conditioning;
      if (d.implicit.present) c.err("conditioning", "F31", "conditioning fits the parameters of a parametric model; the implicit builder is conditioned by its data");
      for (std::size_t k = 0; k < cd.data.size(); ++k) {
        auto& x = cd.data[k];
        std::string b = "conditioning.data[" + std::to_string(k) + "]";
        if (std::none_of(d.inputs.begin(), d.inputs.end(), [&](const InputDesc& i) { return i.id == x.input; }))
          c.err(b + ".input", "F31", "input '" + x.input + "' is not listed under inputs");
        if (x.uncertainty.present && !(x.uncertainty.si() > 0)) c.err(b + ".uncertainty", "F31", "an uncertainty must be positive");
      }
      if (cd.uncertainty.present && !(cd.uncertainty.si() > 0)) c.err("conditioning.uncertainty", "F31", "the uncertainty must be positive");
      if (cd.orientation_uncertainty.present && !(cd.orientation_uncertainty.si() > 0))
        c.err("conditioning.orientation_uncertainty", "F31", "the uncertainty must be positive");
      if (cd.max_iterations.present && !(cd.max_iterations.si() >= 1)) c.err("conditioning.max_iterations", "F31", "at least one iteration");
      // STX-TN-14 (F36): the Bayesian method needs a prior per parameter and a valid sampler
      if (cd.method != "least_squares" && cd.method != "bayesian")
        c.err("conditioning.method", "F36", "method must be least_squares or bayesian (got '" + cd.method + "')");
      if (cd.method != "bayesian" && (cd.walkers.present || cd.steps.present || cd.burn_in.present || cd.thin.present || cd.seed.present))
        c.err("conditioning.sampler", "F36", "a sampler is used only by conditioning.method: bayesian");
      if (cd.method == "bayesian") {
        const double np = static_cast<double>(cd.parameters.size());
        if (cd.walkers.present) {
          double w = cd.walkers.si();
          if (!(w >= std::max(4.0, 2 * np)) || std::fmod(w, 2.0) != 0)
            c.err("conditioning.sampler.walkers", "F36", "an even number of walkers, at least 4 and at least twice the number of parameters");
        }
        if (cd.steps.present && !(cd.steps.si() >= 2)) c.err("conditioning.sampler.steps", "F36", "at least 2 steps");
        if (cd.burn_in.present && !(cd.burn_in.si() >= 0)) c.err("conditioning.sampler.burn_in", "F36", "burn-in cannot be negative");
        if (cd.thin.present && !(cd.thin.si() >= 1)) c.err("conditioning.sampler.thin", "F36", "thin is at least 1");
        if (cd.seed.present && !(cd.seed.si() >= 0)) c.err("conditioning.sampler.seed", "F36", "the seed cannot be negative");
        for (std::size_t k = 0; k < cd.parameters.size(); ++k) {
          auto& fp = cd.parameters[k];
          std::string b = "conditioning.parameters[" + std::to_string(k) + "]";
          if (!fp.min.present || !fp.max.present) c.err(b, "F36", "a Bayesian fit needs the range [min, max] of '" + fp.path + "' (its prior when no distribution is given)");
          for (auto& u : d.uncertainty.parameters) {
            if (u.path != fp.path) continue;
            if (u.distribution != "normal" && u.distribution != "uniform" && u.distribution != "triangular")
              c.err(b, "F36", "the prior of '" + fp.path + "' must be normal, uniform or triangular (got '" + u.distribution + "')");
          }
        }
      }
      for (std::size_t k = 0; k < cd.parameters.size(); ++k) {
        auto& fp = cd.parameters[k];
        std::string b = "conditioning.parameters[" + std::to_string(k) + "]";
        const Param* p = d.param_by_name(fp.path);
        if (!p) { c.err(b + ".path", "F31", "'" + fp.path + "' names no parameter of this description"); continue; }
        if (p->values.size() != 1) { c.err(b + ".path", "F31", "'" + fp.path + "' is not a scalar parameter"); continue; }
        if (!fp.min.present || !fp.max.present) continue;
        if (fp.min.dim != p->dim || fp.max.dim != p->dim) { c.err(b, "F31", "min and max must have the dimension of '" + fp.path + "'"); continue; }
        if (!(fp.min.si() < p->si() && p->si() < fp.max.si()))
          c.err(b, "F31", "the start value of '" + fp.path + "' (" + p->text + ") must lie strictly inside [min, max]");
      }
    }
    // STX-TN-10 (F27, F28): fault rock and stress
    for (std::size_t k = 0; k < d.units.size(); ++k)
      if (d.units[k].clay_fraction.present && !(d.units[k].clay_fraction.si() >= 0 && d.units[k].clay_fraction.si() <= 1))
        c.err("stratigraphy.units[" + std::to_string(k) + "].clay_fraction", "F27", "clay fraction must be in [0, 1]");
    if (d.fault_rock.present) {
      const auto& fr = d.fault_rock;
      if (fr.shale_threshold.present && !(fr.shale_threshold.si() > 0 && fr.shale_threshold.si() <= 1))
        c.err("fault_rock.shale_threshold", "F27", "the shale threshold is a clay fraction in (0, 1]");
      if (fr.thickness_ratio.present && !(fr.thickness_ratio.si() > 0)) c.err("fault_rock.thickness_displacement_ratio", "F27", "must be positive");
      if (!fr.seal_type.empty() && fr.seal_type != "sgr_power") c.err("fault_rock.seal_calibration.type", "F27", "calibration type must be sgr_power");
      if (fr.seal_a.present && !(fr.seal_a.si() > 0)) c.err("fault_rock.seal_calibration.a", "F27", "a must be positive");
      if (!fr.seal_type.empty() && std::none_of(d.units.begin(), d.units.end(), [](const UnitDesc& u) { return u.clay_fraction.present; }))
        c.err("fault_rock.seal_calibration", "F27", "a seal calibration needs clay fractions on the units");
    }
    // STX-TN-17 section 3 (F50): sensitivity designs
    if (d.uncertainty.sensitivity.present) {
      const auto& q = d.uncertainty.sensitivity;
      const std::string b = "uncertainty.sensitivity";
      if (q.method != "saltelli") c.err(b + ".method", "F50", "the sensitivity design is 'saltelli' (Saltelli et al. 2010)");
      if (q.base_samples.present && !(q.base_samples.si() >= 8)) c.err(b + ".base_samples", "F50", "at least 8 base samples");
      if (q.bootstrap.present && !(q.bootstrap.si() >= 100)) c.err(b + ".bootstrap", "F50", "at least 100 bootstrap resamples for an interval");
      if (q.confidence.present && !(q.confidence.si() > 0 && q.confidence.si() < 1)) c.err(b + ".confidence", "F50", "the confidence level is between 0 and 1");
      if (!d.uncertainty.correlations.empty()) c.err(b, "F50", "a Saltelli design needs independent parameters: remove the correlations or use a random ensemble");
      if (d.uncertainty.parameters.empty() && d.uncertainty.scenarios.empty()) c.err(b, "F50", "a sensitivity design needs uncertain parameters or scenarios");
      if (d.conditioning.method == "bayesian") c.err(b, "F50", "posterior members are not a Saltelli design");
    }
    // STX-TN-17 section 1 (F49): geostatistical properties
    {
      std::set<std::string> names;
      for (std::size_t k = 0; k < d.grid.geostatistics.size(); ++k) {
        const auto& q = d.grid.geostatistics[k];
        const std::string b = "grid.geostatistics[" + std::to_string(k) + "]";
        if (!names.insert(q.name).second) c.err(b + ".name", "F49", "'" + q.name + "' is simulated twice");
        bool unit_ok = false;
        for (auto& u : d.units) unit_ok = unit_ok || u.name == q.unit;
        if (!unit_ok) c.err(b + ".unit", "F49", "unknown unit '" + q.unit + "'");
        auto vario_ok = [&](const GridDesc::VariogramDesc& v, const std::string& vb) {
          if (!v.present || v.structures.empty()) { c.err(vb, "F49", "a variogram with at least one structure is stated"); return; }
          double s = v.nugget.present ? v.nugget.si() : 0.0;
          for (auto& st : v.structures) s += st.contribution.si();
          if (std::fabs(s - 1.0) > 1e-6) c.err(vb, "F49", "a standardised variogram: nugget and contributions add up to 1 (they add up to " + fmt_shortest(s) + ")");
          for (auto& st : v.structures)
            for (std::size_t i = 0; i < st.ranges.values.size(); ++i)
              if (!(st.ranges.si(i) > 0)) c.err(vb, "F49", "variogram ranges are positive");
        };
        bool cond_ok = q.condition_to.empty();
        for (auto& up : d.grid.upscaled_logs) cond_ok = cond_ok || (up.name.empty() ? up.curve : up.name) == q.condition_to;
        if (!cond_ok) c.err(b + ".condition_to", "F49", "'" + q.condition_to + "' is not one of grid.upscaled_logs (conditioning data are upscaled logs)");
        if (q.kind == "sgs" || q.kind == "sis") {
          if (!q.search_radii.present) c.err(b + ".search", "F49", "the search neighbourhood (radii, max_data) is stated");
        }
        if (q.kind == "sgs") {
          vario_ok(q.variogram, b + ".variogram");
          if (!q.min.present || !q.max.present) c.err(b + ".distribution", "F49", "the property's range (distribution min and max) is stated for the normal-score tails");
          else if (!(q.min.si() < q.max.si())) c.err(b + ".distribution", "F49", "min is below max");
          if (q.condition_to.empty() && (!q.mean.present || !q.sd.present))
            c.err(b + ".distribution", "F49", "without conditioning data the target distribution (mean and sd) is stated");
        } else if (q.kind == "sis") {
          if (q.codes.size() < 2) c.err(b + ".codes", "F49", "indicator simulation needs at least two facies codes");
          double s = 0;
          for (std::size_t i = 0; i < q.codes.size(); ++i) {
            s += q.codes[i].proportion.si();
            vario_ok(q.codes[i].variogram, b + ".codes[" + std::to_string(i) + "].variogram");
          }
          if (!q.codes.empty() && std::fabs(s - 1.0) > 1e-6) c.err(b + ".codes", "F49", "facies proportions add up to 1 (they add up to " + fmt_shortest(s) + ")");
        } else if (q.kind == "objects") {
          if (q.object_kind.empty()) c.err(b + ".object", "F49", "object facies state their object (channel or lobe)");
          if (q.object_kind == "channel" && !q.width.present) c.err(b + ".object.width", "F49", "a channel states its width");
          if (q.object_kind == "lobe" && !q.radius.present) c.err(b + ".object.radius", "F49", "a lobe states its radius");
          if (q.sinuosity.present && !(q.sinuosity.si() >= 1)) c.err(b + ".object.sinuosity", "F49", "sinuosity is 1 or more");
          if (q.proportion.present && !(q.proportion.si() > 0 && q.proportion.si() < 1)) c.err(b + ".object.proportion", "F49", "the target proportion is between 0 and 1");
          if (!q.geobody.empty()) {
            bool found = false;
            for (auto& gb : d.geobodies) found = found || gb.name == q.geobody;
            if (!found) c.err(b + ".object.geobody", "F49", "unknown geobody '" + q.geobody + "'");
          }
        }
      }
    }
    // STX-TN-16 section 3 (F48): pick and velocity uncertainty
    {
      const auto& u = d.uncertainty;
      for (std::size_t k = 0; k < u.picks.size(); ++k) {
        const auto& q = u.picks[k];
        const std::string b = "uncertainty.interpretation.picks[" + std::to_string(k) + "]";
        const int hi = d.horizon_index(q.horizon);
        if (hi < 0) { c.err(b + ".horizon", "F48", "unknown horizon '" + q.horizon + "'"); continue; }
        const auto& z = d.horizons[static_cast<std::size_t>(hi)].surface_z_is;
        if (z != "twt_ms" && z != "twt_s") c.err(b + ".horizon", "F48", "pick uncertainty is in time: '" + q.horizon + "' is not a surface in two-way time");
        if (q.sigma.present && !(q.sigma.si() > 0)) c.err(b + ".sigma", "F48", "the pick uncertainty is a positive standard deviation");
        if (q.correlation_length.present && !(q.correlation_length.si() > 0)) c.err(b + ".correlation_length", "F48", "the correlation length is positive");
      }
      for (std::size_t k = 0; k < u.velocity.size(); ++k) {
        const auto& q = u.velocity[k];
        const std::string b = "uncertainty.interpretation.velocity[" + std::to_string(k) + "]";
        bool found = false;
        for (auto& L : d.velocity.layers) found = found || L.name == q.layer;
        if (!found) c.err(b + ".layer", "F48", "unknown velocity_model layer '" + q.layer + "'");
        if (!q.v0_sigma.present && !q.k_sigma.present) c.err(b, "F48", "a velocity uncertainty states v0_sigma, k_sigma or both");
        if (q.v0_sigma.present && !(q.v0_sigma.si() > 0)) c.err(b + ".v0_sigma", "F48", "v0_sigma is a positive standard deviation");
        if (q.k_sigma.present && !(q.k_sigma.si() > 0)) c.err(b + ".k_sigma", "F48", "k_sigma is a positive standard deviation");
        if (q.correlation_length.present && !(q.correlation_length.si() > 0)) c.err(b + ".correlation_length", "F48", "the correlation length is positive");
      }
      for (std::size_t h = 0; h < d.horizons.size(); ++h)
        if (d.horizons[h].surface_perturbation.present && d.horizons[h].surface_z_is.rfind("twt", 0) != 0)
          c.err("stratigraphy.horizons." + d.horizons[h].name + ".surface.perturbation", "F48", "a pick perturbation applies to a surface in two-way time");
    }
    // STX-TN-13 (F33): the uncertainty block
    if (d.uncertainty.present) {
      const auto& u = d.uncertainty;
      if (u.members.present && !(u.members.si() >= 2)) c.err("uncertainty.members", "F33", "an ensemble needs at least 2 members");
      if (u.seed.present && !(u.seed.si() >= 0 && u.seed.si() == std::floor(u.seed.si()))) c.err("uncertainty.seed", "F33", "the seed is a whole number, 0 or more");
      std::set<std::string> seen;
      auto scalar = [&](const std::string& path) {
        const Param* p = d.param_by_name(path);
        return p && p->present && p->values.size() == 1 ? p : nullptr;
      };
      for (std::size_t k = 0; k < u.parameters.size(); ++k) {
        const auto& q = u.parameters[k];
        const std::string b = "uncertainty.parameters[" + std::to_string(k) + "]";
        const Param* base = scalar(q.path);
        if (!base) { c.err(b + ".path", "F33", "'" + q.path + "' is not a scalar parameter of this description"); continue; }
        if (!seen.insert(q.path).second) c.err(b + ".path", "F33", "'" + q.path + "' has two distributions");
        auto same_dim = [&](const Param& p, const char* nm) {
          if (p.present && base->dim != Dimension::Any && p.dim != base->dim && p.dim != Dimension::Any)
            c.err(b + "." + nm, "F33", std::string(nm) + " must be in a unit of the parameter's dimension (" + dimension_name(base->dim) + ")");
        };
        for (auto [p, nm] : {std::pair{&q.mean, "mean"}, {&q.sd, "sd"}, {&q.min, "min"}, {&q.max, "max"}, {&q.mode, "mode"}, {&q.values, "values"}})
          same_dim(*p, nm);
        if (q.distribution == "uniform") {
          if (!q.min.present || !q.max.present || !(q.min.si() < q.max.si())) c.err(b, "F33", "uniform needs min < max");
        } else if (q.distribution == "normal") {
          if (!q.mean.present || !q.sd.present || !(q.sd.si() > 0)) c.err(b, "F33", "normal needs a mean and sd > 0");
          if (q.min.present && q.max.present && !(q.min.si() < q.max.si())) c.err(b, "F33", "truncation needs min < max");
        } else if (q.distribution == "triangular") {
          if (!q.min.present || !q.mode.present || !q.max.present || !(q.min.si() <= q.mode.si() && q.mode.si() <= q.max.si() && q.min.si() < q.max.si()))
            c.err(b, "F33", "triangular needs min <= mode <= max, min < max");
        } else if (q.distribution == "discrete") {
          double sw = 0;
          bool pos = true;
          for (std::size_t i = 0; i < q.weights.values.size(); ++i) { sw += q.weights.si(i); pos = pos && q.weights.si(i) > 0; }
          if (!q.values.present || !q.weights.present || q.values.values.size() != q.weights.values.size() || q.values.values.empty())
            c.err(b, "F33", "discrete needs values and as many weights");
          else if (!pos || std::fabs(sw - 1) > 1e-6) c.err(b + ".weights", "F33", "weights must be positive and sum to 1");
        } else {
          c.err(b + ".distribution", "F33", "distribution is uniform, normal, triangular or discrete");
        }
      }
      for (std::size_t k = 0; k < u.correlations.size(); ++k) {
        const auto& cr = u.correlations[k];
        const std::string b = "uncertainty.correlations[" + std::to_string(k) + "]";
        if (!seen.count(cr.a) || !seen.count(cr.b) || cr.a == cr.b) c.err(b, "F33", "a correlation joins two different parameters that have distributions");
        if (cr.rho.present && !(std::fabs(cr.rho.si()) < 1)) c.err(b + ".rho", "F33", "|rho| must be below 1");
      }
      for (std::size_t k = 0; k < u.scenarios.size(); ++k) {
        const auto& sc = u.scenarios[k];
        double sw = 0;
        bool pos = true;
        for (auto& o : sc.options) { sw += o.weight.present ? o.weight.si() : 0; pos = pos && o.weight.present && o.weight.si() > 0; }
        if (sc.options.size() < 2) c.err("uncertainty.scenarios[" + std::to_string(k) + "]", "F33", "a scenario needs at least two options");
        else if (!pos || std::fabs(sw - 1) > 1e-6) c.err("uncertainty.scenarios[" + std::to_string(k) + "]", "F33", "option weights must be positive and sum to 1");
      }
      for (std::size_t k = 0; k < u.maps.size(); ++k) {
        const auto& m = u.maps[k];
        const std::string b = "uncertainty.probability_maps[" + std::to_string(k) + "]";
        if (m.kind == "unit") {
          if (d.unit_index(m.unit) < 0) c.err(b + ".unit", "F33", "unknown unit '" + m.unit + "'");
          if (!m.depth.present) c.err(b + ".depth", "F33", "a unit map needs a depth");
        } else if (m.kind == "closure") {
          if (d.horizon_index(m.horizon) < 0) c.err(b + ".horizon", "F33", "unknown horizon '" + m.horizon + "'");
        } else {
          c.err(b + ".kind", "F33", "kind is unit or closure");
        }
        if (m.spacing.present && !(m.spacing.si() > 0)) c.err(b + ".spacing", "F33", "spacing must be positive");
      }
    }
    // STX-TN-12 (F32): the grid block
    if (d.grid.present) {
      const auto& g = d.grid;
      if (g.spacing.present && (g.spacing.values.empty() || g.spacing.values.size() > 2))
        c.err("grid.spacing", "F32", "the grid spacing is one value or [dx, dy]");
      for (double v : g.spacing.si_list())
        if (!(v > 0)) c.err("grid.spacing", "F32", "the grid spacing must be positive");
      if (g.margin.present && g.margin.si() < 0) c.err("grid.margin", "F32", "the margin must not be negative");
      if (g.min_thickness.present && !(g.min_thickness.si() >= 0)) c.err("grid.min_thickness", "F32", "the minimum thickness must not be negative");
      std::set<std::string> seen;
      int prev_base = -2;
      for (std::size_t k = 0; k < g.zones.size(); ++k) {
        const auto& z = g.zones[k];
        const std::string b = "grid.zones[" + std::to_string(k) + "]";
        const int u = d.unit_index(z.unit);
        if (u < 0) { c.err(b + ".unit", "F32", "unknown unit '" + z.unit + "'"); continue; }
        if (!seen.insert(z.unit).second) c.err(b + ".unit", "F32", "unit '" + z.unit + "' is gridded twice");
        const int ht = d.horizon_index(d.units[static_cast<std::size_t>(u)].top), hb = d.horizon_index(d.units[static_cast<std::size_t>(u)].base);
        if (ht < 0 || hb < 0) c.err(b + ".unit", "F32", "unit '" + z.unit + "' is not between two horizons");
        if (prev_base != -2 && ht != prev_base) c.err(b + ".unit", "F32", "zones must be consecutive units, listed top to base");
        prev_base = hb;
        if (z.layers.present && !(z.layers.si() >= 1 && std::fabs(z.layers.si() - std::round(z.layers.si())) < 1e-9))
          c.err(b + ".layers", "F32", "a zone needs a whole number of layers, 1 or more");
        if (z.layering != "proportional" && z.layering != "follow_top" && z.layering != "follow_base")
          c.err(b + ".layering", "F32", "layering is proportional, follow_top or follow_base");
        if (z.layering != "proportional" && !(z.layer_thickness.present && z.layer_thickness.si() > 0))
          c.err(b + ".layer_thickness", "F32", "follow_top and follow_base need a positive layer thickness");
      }
      for (std::size_t k = 0; k < g.facies.size(); ++k)
        if (d.unit_index(g.facies[k].unit) < 0) c.err("grid.facies[" + std::to_string(k) + "].unit", "F32", "unknown unit '" + g.facies[k].unit + "'");
      if (!g.porosity.empty() && g.porosity != "compaction") c.err("grid.porosity", "F32", "porosity comes from compaction (stage 2D)");
      if (g.porosity == "compaction" && !d.compaction.present())
        c.err("grid.porosity", "F32", "porosity from compaction needs a compaction block (stage 2D)");
      std::set<std::string> known = {"ZONE", "BLOCK"};
      if (!g.facies.empty()) known.insert("FACIES");
      if (g.porosity == "compaction") known.insert("PORO");
      for (std::size_t k = 0; k < g.transforms.size(); ++k) {
        const auto& t = g.transforms[k];
        const std::string b = "grid.transforms[" + std::to_string(k) + "]";
        if (!known.count(t.from)) c.err(b + ".from", "F32", "'" + t.from + "' is not a property computed before this transform");
        if (t.law != "linear" && t.law != "exponential" && t.law != "power") c.err(b + ".law", "F32", "law is linear, exponential or power");
        if (known.count(t.name)) c.err(b + ".name", "F32", "property '" + t.name + "' is defined twice");
        known.insert(t.name);
      }
      if (!g.transmissibility.empty()) {
        if (g.transmissibility != "manzocchi_1999") c.err("grid.fault_transmissibility.calibration", "F32", "calibration is manzocchi_1999");
        if (!d.fault_rock.present) c.err("grid.fault_transmissibility", "F32", "fault transmissibility needs a fault_rock block (stage 2F)");
        if (!known.count("PERMX") || !known.count("PERMY"))
          c.err("grid.fault_transmissibility", "F32", "fault transmissibility needs PERMX and PERMY (transforms)");
      }
    }
    if (d.stress.present) {
      const auto& st2 = d.stress;
      if (st2.regime() == st2.tensor_gradient.present)
        c.err("stress", "F28", "give the stress as a regime or as a tensor gradient, one of the two");
      for (auto [p, nm] : {std::pair{&st2.sv_gradient, "regime.sv_gradient"}, {&st2.shmax_ratio, "regime.shmax_ratio"},
                           {&st2.shmin_ratio, "regime.shmin_ratio"}, {&st2.pore_pressure_gradient, "pore_pressure_gradient"},
                           {&st2.friction, "friction"}})
        if (p->present && !(p->si() > 0)) c.err(std::string("stress.") + nm, "F28", std::string(nm) + " must be positive");
      if (st2.cohesion.present && st2.cohesion.si() < 0) c.err("stress.cohesion", "F28", "cohesion must not be negative");
    }
    for (auto& h : d.curvature_horizons)
      if (d.horizon_index(h) < 0) c.err("outputs.curvature.horizons", "O1", "unknown horizon '" + h + "'");
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
        if (varies && !level_hw(fd))
          c.err(b + ".hanging_wall.type", "F26", "fault '" + fd.name + "' displaces horizons by different amounts (phases, or deposition during or after "
                                                 "its events): it takes a planar geometry with hanging_wall layer_rigid, rigid, decay or "
                                                 "linear_between_faults");
      }
    }
    for (std::size_t k = 0; k < d.faults.size(); ++k) {
      const auto& fd = d.faults[k];
      const std::string b = "faults[" + std::to_string(k) + "]";
      growth(fd.growth, b);
      if (fd.growth.present()) {
        if (!level_hw(fd))
          c.err(b + ".hanging_wall.type", "F22", "growth faults are planar with hanging_wall layer_rigid, rigid, decay or linear_between_faults "
                                                 "(each layer moves by its own throw); curved and fold-kinematic hanging walls are not "
                                                 "combined with growth strata yet");
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
      } else if (!from_well && defs == 0 && !d.implicit.present)  // implicit models: the data define the horizons (STX-TN-11)
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
      // F44 (stage 4C): a survey read from a file, or given with a method, states its method and
      // azimuth reference; a trajectory written as in v3 keeps its documented meaning
      const std::string bw = "wells[" + std::to_string(k) + "].trajectory";
      if ((!w.survey_file.empty() || !w.survey_method.empty() || !w.azimuth_reference.empty()) &&
          (w.survey_method.empty() || w.azimuth_reference.empty()))
        c.err(bw, "F44", "a survey states its calculation method and azimuth reference (they are never assumed)");
      if (w.azimuth_reference == "magnetic_north" && !w.declination.present)
        c.err(bw + ".declination", "F44", "a magnetic azimuth reference needs the declination");
      if (!w.ages.empty() && d.correlation_rule.empty())
        c.err("wells[" + std::to_string(k) + "].ages", "F44", "age data need well_correlation.rule (the rule a correlation rests on is stated, never assumed)");
      for (std::size_t q = 0; q < w.ages.size(); ++q) {
        const auto& a = w.ages[q];
        if (a.depth_base.present && a.depth_base.si() < a.depth.si())
          c.err("wells[" + std::to_string(k) + "].ages[" + std::to_string(q) + "]", "F44", "depth_base lies above depth");
        if (a.sigma.present && a.sigma.si() < 0) c.err("wells[" + std::to_string(k) + "].ages[" + std::to_string(q) + "].sigma", "F44", "an age uncertainty is zero or positive");
      }
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
    // stage 4C (STX-TN-16): time inputs (F45), coordinate systems (F46), geological maps (F47)
    {
      auto has_input = [&](const std::string& id) { return std::any_of(d.inputs.begin(), d.inputs.end(), [&](const InputDesc& i) { return i.id == id; }); };
      bool time_surfaces = false;
      for (auto& h : d.horizons) time_surfaces = time_surfaces || h.surface_z_is.rfind("twt", 0) == 0;
      if (time_surfaces && !d.velocity.present)
        c.err("velocity_model", "F45", "horizons given in two-way time need a velocity_model (time and depth are never mixed without a stated conversion)");
      const auto& V = d.velocity;
      int prev = -1;
      for (std::size_t k = 0; k < V.layers.size(); ++k) {
        const auto& L = V.layers[k];
        const std::string b = "velocity_model.layers[" + std::to_string(k) + "]";
        const int h = d.horizon_index(L.base);
        if (h < 0) { c.err(b + ".base", "F45", "unknown horizon '" + L.base + "'"); continue; }
        if (h <= prev) c.err(b + ".base", "F45", "layers are listed from the top down: '" + L.base + "' lies above the previous base");
        prev = h;
        if (d.horizons[static_cast<std::size_t>(h)].surface_z_is.rfind("twt", 0) != 0)
          c.err(b + ".base", "F45", "a layer's base horizon is given in two-way time (stratigraphy.horizons[].surface z_is: twt_ms or twt_s)");
        if (L.v0.present && !(L.v0.si() > 0)) c.err(b + ".v0", "F45", "v0 must be positive");
        if (!L.k.present && !L.fit_k) c.err(b + ".k", "F45", "state k (0 for a constant velocity) or fit_k: true");
      }
      for (auto& h : d.horizons)
        if (h.surface_z_is.rfind("twt", 0) == 0 && V.present && std::none_of(V.layers.begin(), V.layers.end(), [&](const Description::VelocityLayer& L) { return L.base == h.name; }))
          c.err("stratigraphy.horizons." + h.name, "F45", "a horizon in two-way time is the base of a velocity_model layer");
      for (auto& [w, r] : V.exclude)
        if (std::none_of(d.wells.begin(), d.wells.end(), [&](const WellDesc& x) { return x.name == w; }))
          c.err("velocity_model.calibration.exclude", "F45", "unknown well '" + w + "'");
        else if (r.empty())
          c.err("velocity_model.calibration.exclude", "F45", "a well is excluded with its reason");
      for (std::size_t k = 0; k < d.wells.size(); ++k)
        for (auto& [h, p] : d.wells[k].tops_twt)
          if (d.horizon_index(h) < 0) c.err("wells[" + std::to_string(k) + "].tops_twt." + h, "F45", "unknown horizon '" + h + "'");
      if (d.seismic.present) {
        if (!has_input(d.seismic.input)) c.err("seismic.input", "I1", "input '" + d.seismic.input + "' is not listed under inputs");
        if (d.seismic.domain == "time" && !V.present) c.err("seismic.domain", "F45", "a time volume is drawn in depth through a velocity_model");
      }
      // F46: an input in another CRS is transformed (PROJ); a local frame is never mixed with a CRS
      for (auto& in : d.inputs)
        if (in.has_crs) {
          if (in.crs.epsg <= 0 && in.crs.wkt.empty()) c.err("inputs." + in.id + ".crs", "F46", "an input's crs gives an EPSG code or WKT");
          else if (d.crs.is_local()) c.err("inputs." + in.id + ".crs", "F46", "the model is in a local frame (epsg 0): an input in a coordinate reference system cannot be placed in it");
        }
      const auto& G = d.geological_map;
      if (G.present) {
        if (!d.implicit.present) c.err("geological_map", "F47", "a model from a map is built by the implicit builder: add an implicit block");
        for (auto* s : {&G.units, &G.faults, &G.measurements})
          if (!s->input.empty() && !has_input(s->input)) c.err("geological_map", "I1", "input '" + s->input + "' is not listed under inputs");
        if (!G.dem_input.empty() && !has_input(G.dem_input)) c.err("geological_map.elevation_model", "I1", "input '" + G.dem_input + "' is not listed under inputs");
        if (G.units.field.empty()) c.err("geological_map.units.field", "F47", "name the attribute field that carries the unit");
        if (G.units.format == "geopackage" && G.units.layer.empty()) c.err("geological_map.units.layer", "F47", "a GeoPackage map names its layer");
        for (std::size_t k = 0; k < G.unit_order.size(); ++k) {
          const auto& [mu, base] = G.unit_order[k];
          if (base.empty()) {  // the oldest units: below the model's last horizon
            for (std::size_t q = k + 1; q < G.unit_order.size(); ++q)
              if (!G.unit_order[q].second.empty()) c.err("geological_map.unit_order", "F47", "units without a base are the oldest: list them last");
            continue;
          }
          if (d.horizon_index(base) < 0) c.err("geological_map.unit_order", "F47", "unknown horizon '" + base + "' (base of " + mu + ")");
        }
        if (!G.measurements.input.empty() && (G.dip_field.empty() || G.direction_field.empty()))
          c.err("geological_map.measurements", "F47", "measurements name their dip and direction fields and convention");
      }
    }
    // implicit (data-conditioned) builder
    if (d.implicit.present && d.implicit.data.empty())
      c.err("implicit.data", "I1", "the implicit builder needs data (implicit.data or a geological_map)");
    if (d.implicit.present) {
      for (std::size_t k = 0; k < d.implicit.data.size(); ++k) {
        auto& x = d.implicit.data[k];
        std::string b = "implicit.data[" + std::to_string(k) + "]";
        if (x.input.rfind("geological_map:", 0) != 0 &&
            std::none_of(d.inputs.begin(), d.inputs.end(), [&](const InputDesc& i) { return i.id == x.input; }))
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
    // STX-TN-11 (F29): solver and uncertainties
    if (d.implicit.present) {
      const auto& im = d.implicit;
      if (im.solver_type != "auto" && im.solver_type != "dense" && im.solver_type != "partition_of_unity" && im.solver_type != "discrete")
        c.err("implicit.solver.type", "F29", "solver type must be auto, dense, partition_of_unity or discrete");
      // STX-TN-14 (F35): free horizon values are solved on the discrete grid
      if (im.horizon_values == "free" && im.solver_type != "auto" && im.solver_type != "discrete")
        c.err("implicit.horizon_values", "F35", "free horizon values need the discrete solver (implicit.solver.type: discrete or auto)");
      if (im.horizon_values == "free" && !im.isochore_horizons.empty())
        c.err("implicit.isochore_horizons", "F35", "isochore horizons are not needed with free horizon values: remove them or use fixed values");
      for (auto* p : {&im.discrete_spacing, &im.discrete_vertical_spacing})
        if (p->present && !(p->si() > 0)) c.err("implicit.discrete", "F35", "a discrete grid spacing must be positive");
      if (im.discrete_regularisation.present && !(im.discrete_regularisation.si() > 0))
        c.err("implicit.discrete.regularisation", "F35", "the regularisation weight must be positive");
      if (im.discrete_orientation_weight.present && !(im.discrete_orientation_weight.si() > 0))
        c.err("implicit.discrete.orientation_weight", "F35", "the orientation weight must be positive");
      if (im.discrete_max_iterations.present && !(im.discrete_max_iterations.si() >= 10))
        c.err("implicit.discrete.max_iterations", "F35", "at least 10 iterations per level");
      // STX-TN-14 sections 2 and 3: fold frames and fault relations
      for (std::size_t k = 0; k < im.folds.size(); ++k) {
        auto& ff = im.folds[k];
        std::string bb = "implicit.folds[" + std::to_string(k) + "]";
        if (!ff.series.empty() && std::none_of(d.series.begin(), d.series.end(), [&](const SeriesDesc& sr) { return sr.name == ff.series; }))
          c.err(bb + ".series", "R1", "unknown series '" + ff.series + "'");
        if (!ff.foliation.empty() && std::none_of(d.inputs.begin(), d.inputs.end(), [&](const InputDesc& x) { return x.id == ff.foliation; }))
          c.err(bb + ".foliation", "R1", "unknown input '" + ff.foliation + "'");
        if (ff.order.present && !(ff.order.si() >= 1 && ff.order.si() <= 3)) c.err(bb + ".order", "I3", "the Fourier order must be 1, 2 or 3");
        if (ff.wavelength.present && !(ff.wavelength.si() > 0)) c.err(bb + ".wavelength", "I3", "the wavelength must be positive");
      }
      if (im.relation_distance.present && !(im.relation_distance.si() > 0))
        c.err("implicit.relation_distance", "I3", "the relation distance must be positive");
      if (im.points_per_patch.present && !(im.points_per_patch.si() >= 30))
        c.err("implicit.solver.points_per_patch", "F29", "a patch needs at least 30 interface data");
      if (im.overlap.present && !(im.overlap.si() > 1 && im.overlap.si() <= 3))
        c.err("implicit.solver.overlap", "F29", "overlap must be in (1, 3]: at 1 or below the patches leave gaps");
      for (auto& hn : im.isochore_horizons) {
        int h = d.horizon_index(hn);
        if (h < 0) { c.err("implicit.isochore_horizons", "F29", "unknown horizon '" + hn + "'"); continue; }
        int sr = d.series.empty() ? -1 : d.series_of(h);
        bool above = false;  // a shallower horizon of the same series, itself from the field
        for (int a = h - 1; a >= 0 && !above; --a)
          above = (d.series.empty() ? -1 : d.series_of(a)) == sr &&
                  std::find(im.isochore_horizons.begin(), im.isochore_horizons.end(), d.horizons[static_cast<std::size_t>(a)].name) == im.isochore_horizons.end();
        if (!above) c.err("implicit.isochore_horizons", "F29", "horizon '" + hn + "' needs a shallower horizon of its series that comes from the field");
      }
      for (std::size_t k = 0; k < im.data.size(); ++k)
        if (im.data[k].uncertainty.present && im.data[k].uncertainty.si() < 0)
          c.err("implicit.data[" + std::to_string(k) + "].uncertainty", "F29", "an uncertainty must be zero or positive");
    }
    for (std::size_t k = 0; k < d.wells.size(); ++k)
      if (d.wells[k].top_uncertainty.present && d.wells[k].top_uncertainty.si() < 0)
        c.err("wells[" + std::to_string(k) + "].top_uncertainty", "F29", "an uncertainty must be zero or positive");
    // STX-TN-11 (F30): throws, length and tips estimated from data; tips stated otherwise
    for (std::size_t k = 0; k < d.faults.size(); ++k) {
      auto& f = d.faults[k];
      std::string b = "faults[" + std::to_string(k) + "]";
      if (!f.throw_estimate) {
        if (f.upper.rule.empty()) c.err(b + ".extent.upper_tip", "F10", "upper_tip is required (or estimate the throw from data)");
        if (f.lower.rule.empty()) c.err(b + ".extent.lower_tip", "F10", "lower_tip is required (or estimate the throw from data)");
        continue;
      }
      if (!d.implicit.present) c.err(b + ".throw.estimate", "F30", "throws are estimated from data by the implicit builder: add an implicit section");
      if (f.geom_type != "sticks" && f.geom_type != "planar") c.err(b + ".throw.estimate", "F30", "throws are estimated on planar faults");
      if (f.geom_type == "sticks" && f.sticks_fit != "plane") c.err(b + ".throw.estimate", "F30", "throws are estimated on planar faults (fit: plane)");
      if (f.throw_max.present) c.err(b + ".throw.max", "F30", "give either throw.max or throw.estimate, not both");
      if (f.estimate_band.present && !(f.estimate_band.si() > 0)) c.err(b + ".throw.estimate.band", "F30", "band must be positive");
      if (f.estimate_spacing.present && !(f.estimate_spacing.si() > 0))
        c.err(b + ".throw.estimate.station_spacing", "F30", "station spacing must be positive");
      const auto& pt = f.profile_strike.type;
      if (pt != "elliptical" && pt != "flat_topped" && pt != "sine" && pt != "linear")
        c.err(b + ".throw.profile_strike", "F30", "the strike profile must be elliptical, flat_topped, sine or linear");
      if (f.strike_component()) c.err(b + ".slip_mode", "F30", "throws are estimated for dip-slip faults");
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
        if (d.implicit.present && f.geom_type != "sticks")
          c.err(b + ".geometry.type", "F14", "the implicit builder takes faults from sticks (geometry.type: sticks, fit: profile for a curved fault)");
        if (!shear_hw)
          c.err(b + ".hanging_wall.type", "F14", "a curved fault needs a shear or flow hanging wall (vertical_shear, inclined_shear or "
                                                 "fault_parallel_flow): "
                                                 "a rigid hanging wall would open a void above the curved fault");
        if (f.slip_mode != "dip_slip")
          c.err(b + ".slip_mode", "F14", "a curved fault moves by dip_slip");
        if ((f.upper.rule != "free" || f.lower.rule != "free") && f.geom_type != "sticks")
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
            c.err(b + ".relay.breached", "F15", "the connecting fault of this breached relay could not be placed: declare it "
                                                "explicitly, abutting both faults");
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

// ---- Stage 4B (STX-TN-15 sections 2.1 to 2.4): fault systems and breached relays --------------
// A fault system (polygonal network, flower structure, duplex) and a breached relay expand, before
// typing, into ordinary faults appended to `faults`: every rule and check then applies to each
// generated fault exactly as to a fault written by hand, and the resolved description lists them.
namespace {
std::string num(double v, int dec = 3) { return fmt_fixed(v, dec); }
struct SysParam {
  Param p;
  bool ok = false;
};
/// A scalar or tuple parameter of a fault system, in SI.
SysParam sys_param(const Node& sys, const char* key, const std::string& b, std::vector<Issue>& out, bool required = true) {
  SysParam r;
  const Node* n = sys.get(key);
  if (!n) {
    if (required) out.push_back({b + "." + key, "F43", std::string("a fault system needs ") + key, true, sys.line});
    return r;
  }
  std::string why;
  if (!parse_param(*n, b + "." + key, r.p, why)) {
    out.push_back({b + "." + key, "E1", why, true, n->line});
    return r;
  }
  r.ok = true;
  return r;
}
std::string sys_str(const Node& sys, const char* key) {
  const Node* n = sys.get(key);
  return n && n->is_string() ? n->str : std::string();
}
/// Evidence and source of a generated value: the evidence of the parameter it derives from.
std::string prov(const Param& from, const std::string& note) { return std::string(" | ") + evidence_name(from.evidence) + " | " + note; }
const Node* find_fault(const Node& doc, const std::string& name) {
  if (const Node* fs = doc.get("faults"))
    for (auto& f : fs->seq)
      if (const Node* n = f.get("name"); n && n->str == name) return &f;
  return nullptr;
}
double param_si(const Node* n, std::size_t i = 0) {
  if (!n) return kNaN;
  Param p;
  std::string why;
  if (!parse_param(*n, "", p, why) || p.values.size() <= i) return kNaN;
  return p.si(i);
}
/// One planar dip-slip fault as YAML (flow style), parsed into a node.
Node planar_fault(const std::string& name, double age, const std::string& pv, double x, double y, double depth, double strike_deg,
                  double dip_deg, double length, double throw_m, const std::string& sense, const std::string& upper,
                  const std::string& lower, const std::string& extra, double dip_taper = 0.25) {
  std::ostringstream o;
  o << "{name: " << name << ", age: \"" << num(age, 4) << " ordinal" << pv << "\",\n"
    << " geometry: {type: planar, reference_point: \"[" << num(x, 2) << ", " << num(y, 2) << "] m" << pv << "\", reference_depth: \""
    << num(depth, 2) << " m" << pv << "\", strike: \"" << num(std::fmod(strike_deg + 720.0, 360.0), 4) << " deg" << pv << "\", dip: \""
    << num(dip_deg, 4) << " deg" << pv << "\"},\n"
    << " extent: {upper_tip: " << upper << ", lower_tip: " << lower << ", length: \"" << num(length, 2) << " m" << pv << "\"},\n"
    << " throw: {max: \"" << num(throw_m, 3) << " m" << pv << "\", sense: " << sense
    << ", profile_dip: {type: flat_topped, taper: \"" << num(dip_taper, 2) << " fraction" << pv << "\"}, profile_strike: {type: flat_topped, taper: \"0.25 fraction"
    << pv << "\"}},\n slip_mode: dip_slip, hanging_wall: {type: rigid}" << extra << "}";
  return parse_document(o.str(), "fault system");
}
std::string tip_at(const std::string& horizon, const std::string& pv) {
  return "{rule: horizon, horizon: " + horizon + ", offset: \"0 m" + pv + "\"}";
}
}  // namespace

static void expand_fault_systems(Node& doc, std::vector<Issue>& out) {
  Node* faults = doc.get("faults");
  auto append = [&](Node f) {
    if (!faults) { doc.set("faults", Node::make_seq()); faults = doc.get("faults"); }
    faults->seq.push_back(std::move(f));
  };
  if (const Node* systems = doc.get("fault_systems")) {
    const Node sys_copy = *systems;  // appending to faults must not invalidate the systems being read
    std::set<std::string> names;
    for (std::size_t k = 0; k < sys_copy.seq.size(); ++k) {
      const Node& sy = sys_copy.seq[k];
      const std::string b = "fault_systems[" + std::to_string(k) + "]", name = sys_str(sy, "name"), type = sys_str(sy, "type");
      if (name.empty() || !names.insert(name).second) { out.push_back({b + ".name", "F43", "every fault system needs a unique name", true, sy.line}); continue; }
      auto age = sys_param(sy, "age", b, out);
      if (!age.ok) continue;
      if (type == "polygonal") {
        auto area = sys_param(sy, "area", b, out), dep = sys_param(sy, "reference_depth", b, out), sp = sys_param(sy, "spacing", b, out),
             dm = sys_param(sy, "dip_mean", b, out), ds = sys_param(sy, "dip_sd", b, out), tr = sys_param(sy, "throw_ratio", b, out),
             seed = sys_param(sy, "seed", b, out);
        const std::string top = sys_str(sy, "unit_top"), base = sys_str(sy, "unit_base");
        if (top.empty() || base.empty()) out.push_back({b, "F43", "a polygonal system needs unit_top and unit_base (the horizons that confine it)", true, sy.line});
        if (!(area.ok && dep.ok && sp.ok && dm.ok && ds.ok && tr.ok && seed.ok) || top.empty() || base.empty()) continue;
        if (area.p.values.size() != 4) { out.push_back({b + ".area", "F43", "area is [[x0, y0], [x1, y1]]", true, sy.line}); continue; }
        PolygonalSpec ps;
        ps.x0 = area.p.si(0); ps.y0 = area.p.si(1); ps.x1 = area.p.si(2); ps.y1 = area.p.si(3);
        ps.unit_top = 0; ps.unit_base = -1;  // the tips follow the horizons; the planner needs only the order
        ps.spacing = sp.p.si(); ps.dip_mean = rad2deg(dm.p.si()); ps.dip_sd = rad2deg(ds.p.si()); ps.throw_ratio = tr.p.si();
        if (ps.spacing < 100 || ps.spacing > 2000) { out.push_back({b + ".spacing", "F43", "polygonal fault spacing from 100 to 2000 m", true, sy.line}); continue; }
        if (ps.dip_mean < 30 || ps.dip_mean > 70) { out.push_back({b + ".dip_mean", "F43", "polygonal fault dips from 30 to 70 degrees", true, sy.line}); continue; }
        if (ps.throw_ratio <= 0 || ps.throw_ratio > 0.1) { out.push_back({b + ".throw_ratio", "F43", "throw ratio above 0 and at most 0.1 (C10 scaling)", true, sy.line}); continue; }
        PolygonalNetwork net;
        try {
          net = polygonal_network(ps, static_cast<std::uint64_t>(seed.p.values[0]));
        } catch (const Error& e) {
          out.push_back({b, "F43", e.what(), true, sy.line});
          continue;
        }
        if (net.faults.size() > 48) {
          out.push_back({b, "F43", "the network has " + std::to_string(net.faults.size()) + " faults; a model takes at most 48 from one system "
                                   "(reduce the area or widen the spacing)", true, sy.line});
          continue;
        }
        const std::string pv = prov(sp.p, "fault system " + name + " (seed " + num(seed.p.values[0], 0) + ")");
        for (std::size_t i = 0; i < net.faults.size(); ++i) {
          const auto& f = net.faults[i];
          if (f.throw_m > 100) { out.push_back({b, "F43", "polygonal fault throws are at most 100 m", true, sy.line}); break; }
          const Vec3 m = (f.a + f.b) * 0.5, d = f.b - f.a;
          // conjugate dips: alternate dip directions across the network
          const double strike = rad2deg(std::atan2(d.x, d.y)) + (i % 2 ? 180.0 : 0.0);
          char nm[16];
          std::snprintf(nm, sizeof nm, "_%02zu", i + 1);
          append(planar_fault(name + nm, age.p.si() + 0.001 * static_cast<double>(i), pv, m.x, m.y, dep.p.si(), strike, f.dip_deg, f.length,
                              f.throw_m, "normal", tip_at(top, pv), tip_at(base, pv), ""));
        }
      } else if (type == "flower") {
        auto hs = sys_param(sy, "vertical_component", b, out), spread = sys_param(sy, "spread", b, out), ns = sys_param(sy, "splays", b, out),
             bd = sys_param(sy, "branch_depth", b, out);
        const std::string master = sys_str(sy, "master"), sense = sys_str(sy, "sense");
        if (!(hs.ok && spread.ok && ns.ok && bd.ok)) continue;
        if (sense != "positive" && sense != "negative") { out.push_back({b + ".sense", "F43", "a flower is positive (transpression) or negative (transtension)", true, sy.line}); continue; }
        const Node* mf = find_fault(doc, master);
        if (!mf) { out.push_back({b + ".master", "F43", "unknown master fault '" + master + "'", true, sy.line}); continue; }
        const Node* mg = mf->get("geometry");
        const double mx = param_si(mg ? mg->get("reference_point") : nullptr, 0), my = param_si(mg ? mg->get("reference_point") : nullptr, 1),
                     mstr = param_si(mg ? mg->get("strike") : nullptr), mdip = param_si(mg ? mg->get("dip") : nullptr),
                     mlen = param_si(mf->get("extent") ? mf->get("extent")->get("length") : nullptr);
        if (!std::isfinite(mx) || !std::isfinite(mstr) || !std::isfinite(mdip) || !std::isfinite(mlen) || rad2deg(mdip) < 80) {
          out.push_back({b + ".master", "F43", "the master of a flower is a planar fault placed by reference point and strike, dipping at least 80 degrees", true, sy.line});
          continue;
        }
        FlowerSpec fs;
        fs.splays = static_cast<int>(ns.p.values[0]);
        fs.spread_deg = rad2deg(spread.p.si());
        fs.vertical_component = hs.p.si();
        fs.branch_depth = bd.p.si();
        FlowerResult fr;
        try {
          fr = flower_splays(fs);
        } catch (const Error& e) {
          out.push_back({b, "F43", e.what(), true, sy.line});
          continue;
        }
        // splays that are not vertical: the master continues the vertical ones to the surface
        std::vector<int> gen;
        for (int i = 0; i < fs.splays; ++i)
          if (fr.dips_deg[static_cast<std::size_t>(i)] < 89.999) gen.push_back(i);
        const double mage = param_si(mf->get("age"));
        if (!(age.p.si() > mage)) { out.push_back({b + ".age", "F43", "the splays abut the master: the system must be younger than '" + master + "'", true, sy.line}); continue; }
        const std::string pv = prov(hs.p, "fault system " + name + " (flower on " + master + ")");
        const double share = fs.vertical_component / static_cast<double>(gen.size());
        int j = 0;
        for (int i : gen) {
          const double f = static_cast<double>(i) / (fs.splays - 1) - 0.5;
          // a splay on the master's right (its dip side) dips back towards it: strike reversed
          const bool right = f > 0;
          const double strike = rad2deg(mstr) + (right ? 180.0 : 0.0);
          char nm[16];
          std::snprintf(nm, sizeof nm, "_%02d", ++j);
          append(planar_fault(name + nm, age.p.si() + 0.001 * j, pv, mx, my, fs.branch_depth, strike, fr.dips_deg[static_cast<std::size_t>(i)], mlen, share,
                              sense == "positive" ? "reverse" : "normal", "{rule: free}",
                              "{rule: depth, depth: \"" + num(fs.branch_depth + 50, 2) + " m" + pv + "\"}",
                              ", abuts: {fault: " + master + ", side: " + (right ? "hanging_wall" : "footwall") + "}"));
        }
      } else if (type == "duplex") {
        auto start = sys_param(sy, "start", b, out), az = sys_param(sy, "strike", b, out), nh = sys_param(sy, "horses", b, out),
             hl = sys_param(sy, "horse_length", b, out), ra = sys_param(sy, "ramp_angle", b, out), sl = sys_param(sy, "slip", b, out),
             fdp = sys_param(sy, "floor_depth", b, out), rdp = sys_param(sy, "roof_depth", b, out), len = sys_param(sy, "length", b, out);
        if (!(start.ok && az.ok && nh.ok && hl.ok && ra.ok && sl.ok && fdp.ok && rdp.ok && len.ok)) continue;
        const double zf = fdp.p.si(), zr = rdp.p.si();
        if (!(zf > zr + 1)) { out.push_back({b + ".floor_depth", "F43", "the floor of a duplex lies below its roof", true, sy.line}); continue; }
        DuplexSpec ds;
        ds.horses = static_cast<int>(nh.p.values[0]);
        ds.horse_length = hl.p.si();
        ds.ramp_angle_deg = rad2deg(ra.p.si());
        ds.slip = sl.p.si();
        ds.thickness = zf - zr;
        DuplexResult dr;
        try {
          dr = duplex(ds);
        } catch (const Error& e) {
          out.push_back({b, "F43", e.what(), true, sy.line});
          continue;
        }
        // Each horse: a thrust with a floor flat (hinterland side), a ramp and a roof flat (foreland
        // side), the flats in bedding so that their slip offsets no horizon, on a fault-parallel-flow
        // hanging wall (fault-bend folding over the ramp). The ramps dip to the right of the strike
        // (towards the hinterland); the system breaks forward, so each thrust is younger than the one
        // behind it and carries the older horses over its ramp (Boyer and Elliott 1982).
        const double st = az.p.si(), fx = -std::cos(st), fy = std::sin(st);
        const double h = 0.5 * (zf - zr) / std::tan(ra.p.si()), U = 4 * len.p.si();
        const std::string pv = prov(sl.p, "fault system " + name + " (" + dr.style + " duplex, shortening " + num(dr.shortening, 0) + " m)");
        for (int i = 0; i < ds.horses; ++i) {
          char nm[16];
          std::snprintf(nm, sizeof nm, "_%02d", i + 1);
          const double x = start.p.si(0) + fx * ds.horse_length * i, y = start.p.si(1) + fy * ds.horse_length * i;
          std::ostringstream o;
          o << "{name: " << name << nm << ", age: \"" << num(age.p.si() + 0.001 * i, 4) << " ordinal" << pv << "\",\n"
            << " geometry: {type: ramp_flat, reference_point: \"[" << num(x, 2) << ", " << num(y, 2) << "] m" << pv << "\", strike: \""
            << num(rad2deg(st), 4) << " deg" << pv << "\", profile: \"[[" << num(-U, 2) << ", " << num(zr, 2) << "], [" << num(-h, 2) << ", "
            << num(zr, 2) << "], [" << num(h, 2) << ", " << num(zf, 2) << "], [" << num(U, 2) << ", " << num(zf, 2) << "]] m" << pv << "\"},\n"
            << " extent: {upper_tip: {rule: free}, lower_tip: {rule: free}, length: \"" << num(len.p.si(), 2) << " m" << pv << "\"},\n"
            << " throw: {slip: \"" << num(ds.slip, 3) << " m" << pv << "\", sense: reverse, profile_strike: {type: flat_topped, taper: \"0.01 fraction"
            << pv << "\"}},\n slip_mode: dip_slip, hanging_wall: {type: fault_parallel_flow}}";
          append(parse_document(o.str(), "fault system"));
        }
      } else {
        out.push_back({b + ".type", "F43", "a fault system is polygonal, flower or duplex (got '" + type + "')", true, sy.line});
      }
    }
  }
  // Breached relays (STX-TN-15 section 2.4, I-023): a connecting fault across the ramp, at the middle
  // of the overlap. Both relay faults terminate against it, each on its own side, so the throw passes
  // from one to the other through the breach: the connecting fault is the surface where the
  // abandoned tips stop, and its offset is the throw of the faults it links (no slip of its own;
  // rigid hanging walls cannot carry the ramp's throw twice).
  if (!faults) return;
  const std::size_t n0 = faults->seq.size();
  for (std::size_t k = 0; k < n0; ++k) {
    Node* f = &faults->seq[k];
    Node* rl = f->get("relay");
    const Node* br = rl ? rl->get("breached") : nullptr;
    if (!br || br->type != Node::Type::Bool || !br->boolean) continue;
    const std::string b = "faults[" + std::to_string(k) + "].relay";
    const std::string fname = f->get("name") ? f->get("name")->str : "", gname = rl->get("fault") ? rl->get("fault")->str : "";
    const Node* g = find_fault(doc, gname);
    if (!g) continue;  // the typed rules report the unknown fault
    auto geo = [&](const Node* x, const char* key, std::size_t i = 0) {
      const Node* gg = x->get("geometry");
      return param_si(gg ? gg->get(key) : nullptr, i);
    };
    auto ext = [&](const Node* x) { return param_si(x->get("extent") ? x->get("extent")->get("length") : nullptr); };
    const double fx = geo(f, "reference_point", 0), fy = geo(f, "reference_point", 1), fs = geo(f, "strike"), fd = geo(f, "dip"),
                 fz = geo(f, "reference_depth"), fl = ext(f);
    const double gx = geo(g, "reference_point", 0), gy = geo(g, "reference_point", 1), gl = ext(g);
    const double fage = param_si(f->get("age")), gage = param_si(g->get("age"));
    if (!std::isfinite(fx + fy + fs + fd + fz + fl + gx + gy + gl + fage + gage)) {
      out.push_back({b + ".breached", "F15", "a breached relay needs two planar faults placed by reference point, strike, dip and depth, "
                                             "with their lengths", true, rl->line});
      continue;
    }
    const double sx = std::sin(fs), sy = std::cos(fs);
    const double c2 = (gx - fx) * sx + (gy - fy) * sy, sep = (gx - fx) * sy - (gy - fy) * sx;
    const double lo = std::max(-0.5 * fl, c2 - 0.5 * gl), hi = std::min(0.5 * fl, c2 + 0.5 * gl), mid = 0.5 * (lo + hi);
    if (!(hi > lo) || std::fabs(sep) < 1.0) continue;  // the typed relay rules report it
    const double ax = fx + sx * mid, ay = fy + sy * mid;     // F's trace at the middle of the overlap
    const double bx = ax + sy * sep, by = ay - sx * sep;     // G's trace there (offset normal to F's strike)
    const double cx = 0.5 * (ax + bx), cy = 0.5 * (ay + by);
    // the connecting fault strikes from F to G and contains F's slip vector, so that neither relay
    // fault carries its hanging wall across it; it reaches well past both traces at every depth (it
    // has no slip, so its own extent changes nothing)
    double cstrike = rad2deg(std::atan2(bx - ax, by - ay));
    const Vec3 tv{std::sin(deg2rad(cstrike)), std::cos(deg2rad(cstrike)), 0.0};
    const Vec3 dh{std::sin(fs + kPi / 2), std::cos(fs + kPi / 2), 0.0};
    const Vec3 slip = dh * std::cos(fd) + Vec3{0, 0, -std::sin(fd)};
    Vec3 nv = cross(tv, slip);
    if (nv.z < 0) nv = nv * -1.0;  // upward normal: towards the hanging wall
    const Vec3 right{tv.y, -tv.x, 0.0};
    if (dot(nv, right) < 0) cstrike += 180;  // dips to the right of its strike
    const double cdip = rad2deg(std::acos(std::clamp(nv.z / norm(nv), -1.0, 1.0)));
    const double cs = deg2rad(cstrike);
    auto right_of = [&](double x, double y) { return (x - cx) * std::cos(cs) - (y - cy) * std::sin(cs) > 0; };
    // F lies on the side of the breach towards its own reference point (behind the overlap); G opposite
    const bool f_right = right_of(fx, fy), g_right = right_of(gx, gy);
    if (f_right == g_right) {
      out.push_back({b + ".breached", "F15", "the relay faults do not lie on opposite sides of the breach", true, rl->line});
      continue;
    }
    const std::string cname = fname + "_" + gname + "_breach";
    Param ev;
    ev.evidence = Evidence::Assumed;
    if (const Node* tm = f->get("throw") ? f->get("throw")->get("max") : nullptr) {
      std::string why;
      parse_param(*tm, "", ev, why);
    }
    const std::string pv = prov(ev, "breached relay of " + fname + " and " + gname);
    // long enough to reach both traces over the model's depth range (the planes dip, so their
    // traces move with depth), with 250 m to spare at each end
    const Node* dom = doc.get("domain");
    double zmax = param_si(dom ? dom->get("depth_range") : nullptr, 1);
    if (!std::isfinite(zmax)) zmax = 2 * fz;
    const double clen = std::fabs(sep) + 2 * std::fabs(zmax) / std::max(std::tan(fd), 0.05) + 500;
    // older than both: they terminate against it
    Node cf = planar_fault(cname, std::min(fage, gage) - 0.5, pv, cx, cy, fz, cstrike, std::max(cdip, 1.0), clen, 0.0, "normal", "{rule: free}",
                           "{rule: free}", ", linking: true");
    append(std::move(cf));
    faults = doc.get("faults");
    auto set_abut = [&](const std::string& nm, bool right) {
      for (auto& x : faults->seq)
        if (x.get("name") && x.get("name")->str == nm) {
          if (x.get("abuts")) {
            out.push_back({b + ".breached", "F15", "fault '" + nm + "' of a breached relay already abuts another fault; it terminates "
                                                   "against the connecting fault instead", true, rl->line});
            continue;
          }
          Node ab = Node::make_map();
          ab.set("fault", Node::make_string(cname));
          ab.set("side", Node::make_string(right ? "hanging_wall" : "footwall"));
          x.set("abuts", ab);
          auto& mp = x.map;
          mp.erase(std::remove_if(mp.begin(), mp.end(), [](const auto& kv) { return kv.first == "relay"; }), mp.end());
        }
    };
    // the hanging wall of the connecting fault is on the right of its strike (right-hand rule)
    set_abut(fname, f_right);
    set_abut(gname, g_right);
  }
}

namespace {
/// A survey CSV: a header naming md, inc and azi with their units in brackets (md [m], inc [deg],
/// azi [deg]); units are taken from the file and never assumed.
bool read_survey_csv(const std::filesystem::path& file, std::vector<std::array<double, 3>>& out, std::string& why) {
  std::string text;
  try { text = read_file(file); } catch (const Error& e) { why = e.what(); return false; }
  std::vector<std::string> lines;
  {
    std::string cur;
    for (char ch : text) {
      if (ch == '\n') { if (!cur.empty() && cur.back() == '\r') cur.pop_back(); lines.push_back(cur); cur.clear(); }
      else cur += ch;
    }
    if (!cur.empty()) lines.push_back(cur);
  }
  auto split = [](const std::string& l) {
    std::vector<std::string> v;
    std::string c;
    for (char ch : l) { if (ch == ',') { v.push_back(c); c.clear(); } else c += ch; }
    v.push_back(c);
    for (auto& x : v) {
      while (!x.empty() && std::isspace(static_cast<unsigned char>(x.front()))) x.erase(x.begin());
      while (!x.empty() && std::isspace(static_cast<unsigned char>(x.back()))) x.pop_back();
    }
    return v;
  };
  std::size_t h = 0;
  while (h < lines.size() && (lines[h].empty() || lines[h][0] == '#')) ++h;
  if (h >= lines.size()) { why = file.generic_string() + ": no header line"; return false; }
  auto head = split(lines[h]);
  int col[3] = {-1, -1, -1};
  double scale[3] = {1, 1, 1};
  const char* names[3] = {"md", "inc", "azi"};
  for (std::size_t i = 0; i < head.size(); ++i) {
    std::string nm = head[i], unit;
    const auto lb = nm.find('['), rb = nm.find(']');
    if (lb != std::string::npos && rb != std::string::npos && rb > lb) { unit = nm.substr(lb + 1, rb - lb - 1); nm = nm.substr(0, lb); }
    while (!nm.empty() && std::isspace(static_cast<unsigned char>(nm.back()))) nm.pop_back();
    for (auto& ch : nm) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    for (int q = 0; q < 3; ++q)
      if (nm == names[q]) {
        col[q] = static_cast<int>(i);
        if (unit.empty()) { why = file.generic_string() + ": column '" + head[i] + "' states no unit (write e.g. 'md [m]', 'inc [deg]')"; return false; }
        const auto u = find_unit(unit);
        if (!u || (q == 0 && u->dim != Dimension::Length) || (q > 0 && u->dim != Dimension::Angle)) {
          why = file.generic_string() + ": unit '" + unit + "' of column '" + nm + "' is not a " + (q == 0 ? "length" : "angle") + " unit";
          return false;
        }
        scale[q] = u->to_si;
      }
  }
  for (int q = 0; q < 3; ++q)
    if (col[q] < 0) { why = file.generic_string() + ": no column '" + names[q] + "'"; return false; }
  for (std::size_t r = h + 1; r < lines.size(); ++r) {
    if (lines[r].empty() || lines[r][0] == '#') continue;
    auto v = split(lines[r]);
    std::array<double, 3> st{};
    for (int q = 0; q < 3; ++q) {
      const std::size_t c = static_cast<std::size_t>(col[q]);
      char* end = nullptr;
      const double x = c < v.size() ? std::strtod(v[c].c_str(), &end) : 0.0;
      if (c >= v.size() || !end || *end) { why = file.generic_string() + ": line " + std::to_string(r + 1) + " is not three numbers"; return false; }
      st[static_cast<std::size_t>(q)] = x * scale[q];
    }
    out.push_back(st);
  }
  if (out.size() < 2) { why = file.generic_string() + ": a survey needs at least two stations"; return false; }
  return true;
}
}  // namespace

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
  expand_fault_systems(doc, rep.issues);
  if (!rep.ok()) return rep;
  Description d;
  Typer t{Ctx{&rep.issues}, doc, d};
  t.run();
  t.semantic();
  for (auto& in : d.inputs) {
    std::filesystem::path p(in.path);
    in.resolved = p.is_absolute() ? p : base_dir / p;
  }
  // stage 4C: well files relative to the description; survey files read here (CSV with units)
  for (std::size_t k = 0; k < d.wells.size(); ++k) {
    auto& w = d.wells[k];
    for (auto& l : w.logs) {
      std::filesystem::path p(l.file);
      l.resolved = p.is_absolute() ? p : base_dir / p;
    }
    if (!w.survey_file.empty()) {
      std::filesystem::path p(w.survey_file);
      if (!p.is_absolute()) p = base_dir / p;
      std::string why;
      if (!read_survey_csv(p, w.survey, why))
        rep.issues.push_back({"wells[" + std::to_string(k) + "].trajectory.file", "F44", why, true, 0});
    }
  }
  if (!rep.ok()) return rep;
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
