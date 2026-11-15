// Implicit (data-conditioned) structure builder: Hermite RBF stratigraphic field, surfaces in the
// restored frame, faults from the v1 fault network, iterated restoration of the data.
#include "stratonix/engines/implicit.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include "stratonix/core/error.hpp"
#include "stratonix/core/linalg.hpp"
#include "stratonix/core/roots.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/wells.hpp"

namespace stratonix {

// ---- RBF field ---------------------------------------------------------------------------------
namespace {
inline double phi(const Vec3& d) {
  double r = norm(d);
  return r * r * r;
}
inline Vec3 dphi(const Vec3& d) {  // gradient of |d|^3 = 3 |d| d
  double r = norm(d);
  return d * (3 * r);
}
inline void hphi(const Vec3& d, double h[9]) {  // Hessian of |d|^3 = 3 (|d| I + d d^T / |d|)
  double r = norm(d);
  double dv[3] = {d.x, d.y, d.z};
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) h[3 * a + b] = 3 * ((a == b ? r : 0.0) + (r > 0 ? dv[a] * dv[b] / r : 0.0));
}
inline double comp(const Vec3& v, int k) { return k == 0 ? v.x : (k == 1 ? v.y : v.z); }
}  // namespace

RbfField::RbfField(std::vector<Value> values, std::vector<Gradient> grads, double nugget, Vec3 centre, double scale, double vertical_stretch)
    : centre_(centre), scale_(scale), stretch_(vertical_stretch) {
  // Orientations constrain the direction of the gradient only: two tangent constraints per
  // orientation (t . grad s = 0), so the magnitude of the gradient is free and the interface
  // values set the spacing of the horizons (unit-gradient constraints conflict with values
  // wherever layer thickness varies).
  const std::size_t ni = values.size();
  if (ni < 4) throw Error(ErrorCode::Build, "implicit field: at least four interface points are needed (have " + std::to_string(ni) + ")");
  for (auto& v : values) xi_.push_back(local(v.p));
  for (auto& g : grads) {
    Vec3 n = normalized(g.g);
    Vec3 t1 = normalized(cross(n, std::fabs(n.z) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0}));
    Vec3 t2 = cross(n, t1);
    for (const Vec3& t : {t1, t2}) {
      yj_.push_back(local(g.p));
      ej_.push_back(normalized(Vec3{t.x, t.y, t.z * stretch_}));  // tangent in stretched coordinates
    }
  }
  const std::size_t nd = yj_.size();
  const std::size_t n = ni + nd + 4;
  std::vector<double> A(n * n, 0.0), rhs(n, 0.0);
  auto at = [&](std::size_t r, std::size_t c) -> double& { return A[r * n + c]; };
  const std::size_t od = ni, op = ni + nd;
  const double nug = nugget / scale_;  // values are scaled with the coordinates
  for (std::size_t r = 0; r < ni; ++r) {
    for (std::size_t c = 0; c < ni; ++c) at(r, c) = phi(xi_[r] - xi_[c]);
    at(r, r) += nug;
    for (std::size_t j = 0; j < nd; ++j) at(r, od + j) = at(od + j, r) = -dot(ej_[j], dphi(xi_[r] - yj_[j]));
    at(r, op) = at(op, r) = 1;
    for (int k = 0; k < 3; ++k) at(r, op + 1 + k) = at(op + 1 + k, r) = comp(xi_[r], k);
    rhs[r] = values[r].v / scale_;
  }
  for (std::size_t l = 0; l < nd; ++l) {
    for (std::size_t j = 0; j < nd; ++j) {
      double h[9];
      hphi(yj_[l] - yj_[j], h);
      const Vec3 &el = ej_[l], &ej = ej_[j];
      double v = 0;
      for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y) v += comp(el, x) * h[3 * x + y] * comp(ej, y);
      at(od + l, od + j) = -v;
    }
    for (int k = 0; k < 3; ++k) at(od + l, op + 1 + k) = at(op + 1 + k, od + l) = comp(ej_[l], k);
    rhs[od + l] = 0;
  }
  solver_ = optimised_kernels() && n >= 256 ? "blocked LU (optimised)" : "Gaussian elimination (reference)";
  solve_dense(A, rhs, n);
  a_.assign(rhs.begin(), rhs.begin() + static_cast<long>(ni));
  beta_.assign(rhs.begin() + static_cast<long>(od), rhs.begin() + static_cast<long>(op));
  c0_ = rhs[op];
  c_ = {rhs[op + 1], rhs[op + 2], rhs[op + 3]};
}

double RbfField::value(const Vec3& xw) const {
  Vec3 x = local(xw);
  double s = c0_ + dot(c_, x);
  for (std::size_t i = 0; i < xi_.size(); ++i) s += a_[i] * phi(x - xi_[i]);
  for (std::size_t j = 0; j < yj_.size(); ++j) s -= beta_[j] * dot(ej_[j], dphi(x - yj_[j]));
  return s * scale_;
}

Vec3 RbfField::gradient(const Vec3& xw) const {
  Vec3 x = local(xw);
  Vec3 g = c_;
  for (std::size_t i = 0; i < xi_.size(); ++i) g = g + dphi(x - xi_[i]) * a_[i];
  for (std::size_t j = 0; j < yj_.size(); ++j) {
    double h[9];
    hphi(x - yj_[j], h);
    const Vec3& e = ej_[j];
    g = g - Vec3{h[0] * e.x + h[1] * e.y + h[2] * e.z, h[3] * e.x + h[4] * e.y + h[5] * e.z, h[6] * e.x + h[7] * e.y + h[8] * e.z} * beta_[j];
  }
  g.z *= stretch_;  // chain rule for the stretched vertical coordinate
  return g;
}

// ---- partition of unity (stage 3A, STX-TN-11 section 1): specified, built by Prompt 35 ----------
PuField::PuField(std::vector<RbfField::Value>, std::vector<RbfField::Gradient>, double, Vec3, double, double, double, double, double, double,
                 const SolverOptions&) {
  throw Error(ErrorCode::Build, "partition-of-unity solver: specified in STX-TN-11, not built yet (Prompt 35)");
}
double PuField::value(const Vec3&) const { return kNaN; }
Vec3 PuField::gradient(const Vec3&) const { return {kNaN, kNaN, kNaN}; }
std::size_t PuField::largest_patch() const { return 0; }
std::vector<double> PuField::leave_one_out() const { return {}; }
std::vector<double> rippa_loo(const std::vector<RbfField::Value>&, const std::vector<RbfField::Gradient>&, double, Vec3, double, double) {
  throw Error(ErrorCode::Build, "leave-one-out errors: specified in STX-TN-11, not built yet (Prompt 35)");
}

// ---- plane fit -------------------------------------------------------------------------------------
FittedPlane fit_plane(const std::vector<Polyline>& sticks) {
  FittedPlane fp;
  Vec3 c{0, 0, 0};
  for (auto& s : sticks)
    for (auto& p : s) { c = c + p; ++fp.points; }
  if (fp.points < 3) throw Error(ErrorCode::Build, "fault sticks: at least three points are needed to fit a plane");
  c = c * (1.0 / static_cast<double>(fp.points));
  double cov[9] = {0};
  for (auto& s : sticks)
    for (auto& p : s) {
      Vec3 d = p - c;
      double v[3] = {d.x, d.y, d.z};
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) cov[3 * a + b] += v[a] * v[b];
    }
  double val[3], vec[9];
  symmetric_eigen3(cov, val, vec);
  Vec3 n{vec[0], vec[3], vec[6]};  // smallest eigenvalue
  if (n.z < 0) n = n * -1.0;
  fp.centroid = c;
  fp.normal = normalized(n);
  fp.dip = std::acos(std::clamp(fp.normal.z, -1.0, 1.0));
  double dipdir = std::atan2(fp.normal.x, fp.normal.y);  // up-normal leans towards the dip direction
  fp.strike = std::fmod(dipdir - kPi / 2 + 4 * kPi, 2 * kPi);
  double ss = 0;
  for (auto& s : sticks)
    for (auto& p : s) ss += std::pow(dot(p - c, fp.normal), 2);
  fp.rms = std::sqrt(ss / static_cast<double>(fp.points));
  if (fp.dip < 1e-3) throw Error(ErrorCode::Build, "fault sticks define a horizontal plane; a fault must dip");
  return fp;
}

// ---- model -------------------------------------------------------------------------------------
ImplicitModel::ImplicitModel(std::shared_ptr<ParametricModel> inner, ModelInfo info, std::vector<DataMisfit> misfit)
    : inner_(std::move(inner)), info_(std::move(info)), misfit_(std::move(misfit)) {}

const void* ImplicitModel::capability(std::string_view name) const {
  if (name == "data_misfit") return &misfit_;
  return inner_->capability(name);
}

// ---- builder -------------------------------------------------------------------------------------
namespace {
desc::Param make_param(const std::string& path, std::vector<double> values, const std::string& unit, Dimension dim, double to_si,
                       desc::Evidence ev, const std::string& source) {
  desc::Param p;
  p.present = true;
  p.path = path;
  p.values = std::move(values);
  if (p.values.size() > 1) p.shape = {p.values.size()};
  p.unit = unit;
  p.dim = dim;
  p.to_si = to_si;
  p.evidence = ev;
  p.source = source;
  p.text = "";
  for (std::size_t i = 0; i < p.values.size(); ++i) p.text += (i ? ", " : "") + fmt_shortest(p.values[i]);
  if (p.values.size() > 1) p.text = "[" + p.text + "]";
  p.text += " " + unit + " | " + desc::evidence_name(ev) + " | " + source;
  return p;
}

desc::Evidence evidence_of(const std::string& s) {
  static const std::map<std::string, desc::Evidence> m = {{"measured", desc::Evidence::Measured}, {"interpreted", desc::Evidence::Interpreted},
                                                         {"calibrated", desc::Evidence::Calibrated}, {"upscaled", desc::Evidence::Upscaled},
                                                         {"literature_derived", desc::Evidence::LiteratureDerived}, {"assumed", desc::Evidence::Assumed},
                                                         {"user_instruction", desc::Evidence::UserInstruction}};
  auto it = m.find(s);
  return it == m.end() ? desc::Evidence::Assumed : it->second;
}

struct Datum {
  Vec3 present, restored;
  Vec3 normal{0, 0, 0};           ///< present-day normal (data)
  Vec3 normal_restored{0, 0, 0};  ///< normal carried through the restoration
  bool oriented = false;
  int horizon = -1;
  std::string input, evidence, kind;
};

/// Surfaces of the field (s = value_h) on a grid in the restored frame, by scanning each column
/// from the top and refining the first crossing (s decreases downwards).
std::vector<Grid2D> extract_surfaces(const RbfField& f, const std::vector<double>& vals, const GridSpec& g, double ztop, double zbot,
                                     std::size_t& clamped) {
  const std::size_t nh = vals.size(), K = 96;
  std::vector<Grid2D> out(nh, Grid2D(g));
  std::vector<std::size_t> clamp_rows(g.ny, 0);
  parallel_for(g.ny, 1, [&](std::size_t j0, std::size_t j1) {
    std::vector<double> sz(K + 1);
    for (std::size_t j = j0; j < j1; ++j)
      for (std::size_t i = 0; i < g.nx; ++i) {
        double x = g.x(i), y = g.y(j);
        for (std::size_t k = 0; k <= K; ++k) sz[k] = f.value({x, y, ztop + (zbot - ztop) * static_cast<double>(k) / K});
        for (std::size_t h = 0; h < nh; ++h) {
          double v = vals[h], z = kNaN;
          for (std::size_t k = 0; k < K && std::isnan(z); ++k) {
            double fa = sz[k] - v, fb = sz[k + 1] - v;
            if (fa >= 0 && fb < 0) {
              double za = ztop + (zbot - ztop) * static_cast<double>(k) / K, zb = ztop + (zbot - ztop) * static_cast<double>(k + 1) / K;
              z = illinois([&](double zz) { return f.value({x, y, zz}) - v; }, za, zb, fa, fb, 1e-6);
            }
          }
          if (std::isnan(z)) {  // no crossing inside the domain: hold the surface at the nearer boundary
            z = sz[0] < v ? ztop : zbot;
            ++clamp_rows[j];
          }
          out[h].at(i, j) = z;
        }
      }
  });
  clamped = 0;
  for (auto c : clamp_rows) clamped += c;
  // Horizons must not cross: enforce order (a deeper horizon never above a shallower one).
  for (std::size_t h = 1; h < nh; ++h)
    for (std::size_t k = 0; k < out[h].v.size(); ++k) out[h].v[k] = std::min(out[h].v[k], out[h - 1].v[k]);
  return out;
}
}  // namespace

ModelPtr build_implicit_model(const desc::Description& d, const BuildInputs& in, BuildLog& log) {
  const std::size_t nh = d.horizons.size();
  std::map<std::string, const desc::InputDesc*> inputs;
  for (auto& x : d.inputs) inputs[x.id] = &x;
  // 1. Data.
  std::vector<Datum> data;
  std::set<std::string> ignored_labels;
  std::vector<std::string> provenance;
  for (auto& e : d.implicit.data) {
    auto it = in.data.find(e.input);
    if (it == in.data.end()) throw Error(ErrorCode::Build, "implicit data input '" + e.input + "' was not loaded");
    const auto* inp = inputs.at(e.input);
    std::size_t used = 0;
    for (auto& dp : it->second) {
      int h = d.horizon_index(dp.label);
      if (h < 0) { ignored_labels.insert(dp.label); continue; }
      Datum dt;
      dt.present = dt.restored = dp.p;
      dt.horizon = h;
      dt.oriented = dp.oriented;
      dt.normal = dt.normal_restored = dp.n;
      dt.input = e.input;
      dt.evidence = inp->evidence.empty() ? "assumed" : inp->evidence;
      dt.kind = dp.oriented ? "orientation" : "point";
      data.push_back(dt);
      ++used;
    }
    provenance.push_back("Data " + e.input + " (" + inp->path + "): " + std::to_string(used) + " " + (e.kind == "orientations" ? "orientations" : "points") +
                         ", evidence " + (inp->evidence.empty() ? std::string("not stated (treated as assumed)") : inp->evidence) + ", source: " + inp->source);
  }
  if (!ignored_labels.empty()) {
    std::string s;
    for (auto& l : ignored_labels) s += (s.empty() ? "" : ", ") + (l.empty() ? std::string("(no label)") : l);
    log.note("implicit: data labelled " + s + " are not horizons and were not used for the stratigraphic field");
  }
  if (d.implicit.use_well_tops) {
    auto wells = make_well_paths(d);
    std::size_t nt = 0;
    for (std::size_t w = 0; w < d.wells.size(); ++w) {
      auto& wd = d.wells[w];
      if (wd.role != "control") continue;
      for (auto& [hn, p] : wd.tops) {
        int h = d.horizon_index(hn);
        if (h < 0) continue;
        double z = -p.si();  // tops are depths below z = 0 (as in the parametric engine)
        Datum dt;
        dt.present = dt.restored = wells[w].vertical() ? Vec3{wells[w].x(), wells[w].y(), z} : wells[w].at_elevation(z);
        dt.horizon = h;
        dt.input = "well:" + wd.name;
        dt.evidence = desc::evidence_name(p.evidence);
        dt.kind = "well_top";
        data.push_back(dt);
        ++nt;
      }
    }
    if (nt) provenance.push_back("Data: " + std::to_string(nt) + " control-well tops (evidence as labelled in the description)");
  }
  // 2. Horizon values: elevation at the reference well when it has the top, else the mean of the data.
  std::vector<double> vals(nh, kNaN), sum(nh, 0.0);
  std::vector<std::size_t> cnt(nh, 0);
  for (auto& dt : data)
    if (!dt.oriented) { sum[static_cast<std::size_t>(dt.horizon)] += dt.present.z; ++cnt[static_cast<std::size_t>(dt.horizon)]; }
  const desc::WellDesc* rw = d.ref_well.empty() ? nullptr : d.well(d.ref_well);
  for (std::size_t h = 0; h < nh; ++h) {
    if (rw)
      for (auto& [hn, p] : rw->tops)
        if (hn == d.horizons[h].name) vals[h] = -p.si();
    if (std::isnan(vals[h]) && cnt[h]) vals[h] = sum[h] / static_cast<double>(cnt[h]);
    if (std::isnan(vals[h])) throw Error(ErrorCode::Build, "implicit: horizon " + d.horizons[h].name + " has no data point and no reference-well top");
    if (h && !(vals[h] < vals[h - 1]))
      throw Error(ErrorCode::Build, "implicit: the field values of " + d.horizons[h - 1].name + " and " + d.horizons[h].name +
                                        " are not in stratigraphic order (check the data labels)");
  }
  // 3. Faults given by sticks: plane fitted to the sticks.
  desc::Description d2 = d;
  for (auto& f : d2.faults) {
    if (f.geom_type != "sticks") continue;
    auto it = in.lines.find(f.geom_input);
    if (it == in.lines.end()) throw Error(ErrorCode::Build, "fault " + f.name + ": sticks input '" + f.geom_input + "' was not loaded");
    FittedPlane fp = fit_plane(it->second);
    const auto* inp = inputs.at(f.geom_input);
    auto ev = evidence_of(inp->evidence);
    std::string src = "plane fitted to " + std::to_string(fp.points) + " stick points of input " + f.geom_input + " (RMS " + fmt_fixed(fp.rms, 3) + " m)";
    std::string b = "faults." + f.name + ".geometry";
    f.geom_type = "planar";
    f.reference_point = make_param(b + ".reference_point", {fp.centroid.x, fp.centroid.y}, "m", Dimension::Length, 1, ev, src);
    f.reference_depth = make_param(b + ".reference_depth", {-fp.centroid.z}, "m", Dimension::Length, 1, ev, src);
    f.strike = make_param(b + ".strike", {rad2deg(fp.strike)}, "deg", Dimension::Angle, kPi / 180, ev, src);
    f.dip = make_param(b + ".dip", {rad2deg(fp.dip)}, "deg", Dimension::Angle, kPi / 180, ev, src);
    log.value("implicit.fault." + f.name + ".fitted_plane", "strike " + fmt_fixed(rad2deg(fp.strike), 2) + " deg, dip " + fmt_fixed(rad2deg(fp.dip), 2) +
                                                            " deg, centre (" + fmt_fixed(fp.centroid.x, 1) + ", " + fmt_fixed(fp.centroid.y, 1) + ", " +
                                                            fmt_fixed(fp.centroid.z, 1) + "), RMS " + fmt_fixed(fp.rms, 3) + " m");
    provenance.push_back("Fault " + f.name + ": " + src + ", evidence " + desc::evidence_name(ev));
  }
  // 4. Iterate: field in the restored frame -> surfaces -> parametric model with faults -> restore data.
  const Domain dom = d.domain();
  double sp = d.implicit.surface_spacing.present ? d.implicit.surface_spacing.si() : 2 * d.grid_spacing.si();  // Stratonix default: 2 x grid
  GridSpec g = GridSpec::covering(dom, sp);
  // Maximum number of restoration passes; stop earlier when no datum moves by more than 1 cm.
  int iters = d.implicit.iterations.present ? static_cast<int>(d.implicit.iterations.si()) : 6;  // Stratonix default
  const double converged = 0.01;  // m, Stratonix default
  if (d.faults.empty()) iters = 1;
  double nugget = d.implicit.nugget.present ? d.implicit.nugget.si() : 0.0;
  double stretch = d.implicit.vertical_stretch.present ? d.implicit.vertical_stretch.si() : 5.0;  // Stratonix default
  Vec3 centre{dom.x0 + dom.sx / 2, dom.y0 + dom.sy / 2, 0.5 * (dom.z_top + dom.z_base)};
  double scale = 0.5 * std::max({dom.sx, dom.sy, dom.z_top - dom.z_base});
  for (std::size_t h = 0; h < nh; ++h) {
    d2.horizons[h].surface_input = "__implicit_" + d.horizons[h].name;
    d2.horizons[h].surface_format = "implicit";
    d2.horizons[h].surface_z_is = "elevation";
    d2.horizons[h].iso_from.clear();
    d2.horizons[h].iso_input.clear();
  }
  d2.primitives.clear();
  d2.implicit = desc::ImplicitDesc{};
  // The implicit field already interpolates the control-well tops; the parametric engine's
  // correction field must not adjust the surfaces further (misfits are reported, not hidden).
  for (auto& w : d2.wells)
    if (w.role == "control") w.role = "check";
  std::shared_ptr<ParametricModel> inner;
  double moved = 0;
  std::string solver;
  std::size_t unknowns = 0, clamped = 0;
  int passes = 1;
  for (int it = 0; it < iters; ++it) {
    std::vector<RbfField::Value> vs;
    std::vector<RbfField::Gradient> gs;
    for (auto& dt : data) {
      if (dt.oriented) gs.push_back({dt.restored, dt.normal_restored});
      else vs.push_back({dt.restored, vals[static_cast<std::size_t>(dt.horizon)]});
    }
    RbfField field(vs, gs, nugget, centre, scale, stretch);
    solver = field.solver();
    unknowns = field.unknowns();
    auto surfs = extract_surfaces(field, vals, g, dom.z_top, dom.z_base, clamped);
    BuildInputs in2 = in;
    for (std::size_t h = 0; h < nh; ++h) in2.grids[d2.horizons[h].surface_input] = std::move(surfs[h]);
    BuildLog l2;
    inner = std::make_shared<ParametricModel>(d2, in2, l2);
    if (it + 1 == iters || (it > 0 && moved < converged)) {
      for (auto& kv : l2.values) log.value(kv.first, kv.second);
      for (auto& n : l2.notes) log.note(n);
      passes = it + 1;
      break;
    }
    // Restore the data; an orientation is carried by restoring two points 5 m away in its plane.
    std::vector<Vec3> pres;
    for (auto& dt : data) {
      pres.push_back(dt.present);
      if (dt.oriented) {
        Vec3 t1 = normalized(cross(dt.normal, std::fabs(dt.normal.z) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0}));
        Vec3 t2 = cross(dt.normal, t1);
        pres.push_back(dt.present + t1 * 5.0);
        pres.push_back(dt.present + t2 * 5.0);
      }
    }
    std::vector<Vec3> rest(pres.size());
    inner->restore(pres, rest);
    moved = 0;
    std::size_t k = 0;
    for (auto& dt : data) {
      moved = std::max(moved, norm(rest[k] - dt.restored));
      dt.restored = rest[k];
      if (dt.oriented) {
        Vec3 n = normalized(cross(rest[k + 1] - rest[k], rest[k + 2] - rest[k]));
        if (n.z < 0) n = n * -1.0;
        dt.normal_restored = n;
        k += 3;
      } else {
        ++k;
      }
    }
    log.value("implicit.iteration_" + std::to_string(it + 1) + ".max_data_shift_m", fmt_fixed(moved, 4));
  }
  // 5. Misfit at every data point, through the final model (present-day positions).
  std::vector<DataMisfit> mis;
  std::vector<Vec3> pres;
  for (auto& dt : data) pres.push_back(dt.present);
  EvalResult er;
  inner->evaluate(pres, kEvalScalars, er);
  double maxp = 0, sump = 0, maxo = 0;
  std::size_t np = 0, no = 0;
  for (std::size_t i = 0; i < data.size(); ++i) {
    auto& dt = data[i];
    DataMisfit m{dt.input, dt.evidence, d.horizons[static_cast<std::size_t>(dt.horizon)].name, dt.kind, dt.present, 0};
    if (!dt.oriented) {
      m.misfit = er.s(static_cast<std::size_t>(dt.horizon), i);
      maxp = std::max(maxp, std::fabs(m.misfit));
      sump += std::fabs(m.misfit);
      ++np;
    } else {
      // horizon normal of the model at the point, from the scalar's gradient (finite differences)
      const double e = 2.0;
      std::vector<Vec3> q{dt.present + Vec3{e, 0, 0}, dt.present - Vec3{e, 0, 0}, dt.present + Vec3{0, e, 0}, dt.present - Vec3{0, e, 0},
                          dt.present + Vec3{0, 0, e}, dt.present - Vec3{0, 0, e}};
      EvalResult eq;
      inner->evaluate(q, kEvalScalars, eq);
      auto s = [&](std::size_t k) { return eq.s(static_cast<std::size_t>(dt.horizon), k); };
      Vec3 gr = normalized(Vec3{s(0) - s(1), s(2) - s(3), s(4) - s(5)});
      m.misfit = rad2deg(std::acos(std::clamp(dot(gr, dt.normal), -1.0, 1.0)));
      maxo = std::max(maxo, m.misfit);
      ++no;
    }
    mis.push_back(m);
  }
  log.value("implicit.data", std::to_string(np) + " points and tops, " + std::to_string(no) + " orientations");
  log.value("implicit.solver", solver + ", " + std::to_string(unknowns) + " unknowns, " + std::to_string(passes) + " pass(es)");
  if (passes == iters && iters > 1 && moved >= converged)
    log.note("implicit: data restoration not converged after " + std::to_string(iters) + " passes (last shift " + fmt_fixed(moved, 3) + " m)");
  log.value("implicit.surface_spacing_m", fmt_shortest(sp));
  log.value("implicit.vertical_stretch", fmt_shortest(stretch));
  log.value("implicit.misfit_points_max_m", fmt_fixed(maxp, 4));
  log.value("implicit.misfit_points_mean_m", fmt_fixed(np ? sump / static_cast<double>(np) : 0.0, 4));
  log.value("implicit.misfit_orientations_max_deg", fmt_fixed(maxo, 3));
  if (clamped) log.note("implicit: " + std::to_string(clamped) + " surface nodes had no crossing inside the domain and were held at its top or base");
  ModelInfo info = inner->info();
  info.engine = "implicit-rbf";
  info.provenance = provenance;
  return std::make_shared<ImplicitModel>(inner, info, std::move(mis));
}

namespace {
class ImplicitBuilder final : public StructureBuilder {
 public:
  bool accepts(const desc::Description& d) const override { return d.implicit.present; }
  ModelPtr build(const desc::Description& d, const BuildInputs& in, BuildLog& log) const override { return build_implicit_model(d, in, log); }
};
}  // namespace

STRATONIX_REGISTER(implicit_builder, registry::Kind::StructureBuilder, "implicit-rbf", "0.1.0",
                   "Data-conditioned: Hermite RBF stratigraphic field from points, orientations and tops; faults from the fault network",
                   [] { return std::make_unique<ImplicitBuilder>(); })

}  // namespace stratonix
