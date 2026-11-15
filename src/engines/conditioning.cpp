// Stage 3A (STX-TN-11): faults from data and hybrid conditioning. Specified; built by Prompts 35
// and 36.
#include "stratonix/engines/conditioning.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <set>

#include "stratonix/core/error.hpp"
#include "stratonix/core/linalg.hpp"
#include "stratonix/engines/parametric.hpp"
#include "stratonix/engines/wells.hpp"

namespace stratonix {

namespace {

}  // namespace

namespace {
/// Along-strike shapes of stage 2A (faults.cpp profile_value), on t in (0, 1).
double shape_value(const std::string& p, double taper, double t) {
  if (t <= 0 || t >= 1) return 0;
  if (p == "elliptical") return std::sqrt(std::max(0.0, 1 - (2 * t - 1) * (2 * t - 1)));
  if (p == "sine") return std::sin(kPi * t);
  if (p == "linear") return 1 - std::fabs(2 * t - 1);
  if (taper <= 0) return 1;  // flat_topped
  if (t < taper) return 0.5 * (1 - std::cos(kPi * t / taper));
  if (t > 1 - taper) return 0.5 * (1 - std::cos(kPi * (1 - t) / taper));
  return 1;
}
/// Least squares z = a + b x + c y (3 x 3 normal equations); returns false if singular.
bool fit_plane3(const std::vector<std::array<double, 3>>& p, double& a, double& b, double& c) {
  double S[3][3] = {{0}}, r[3] = {0};
  for (auto& q : p) {
    double v[3] = {1, q[0], q[1]};
    for (int i = 0; i < 3; ++i) {
      r[i] += v[i] * q[2];
      for (int j = 0; j < 3; ++j) S[i][j] += v[i] * v[j];
    }
  }
  std::vector<double> A(9), x(3);
  for (int i = 0; i < 3; ++i) {
    x[static_cast<std::size_t>(i)] = r[i];
    for (int j = 0; j < 3; ++j) A[static_cast<std::size_t>(3 * i + j)] = S[i][j];
  }
  // ill-conditioned (data in a line): reject
  double tr = S[0][0] + S[1][1] + S[2][2];
  double det = S[0][0] * (S[1][1] * S[2][2] - S[1][2] * S[2][1]) - S[0][1] * (S[1][0] * S[2][2] - S[1][2] * S[2][0]) +
               S[0][2] * (S[1][0] * S[2][1] - S[1][1] * S[2][0]);
  if (!(std::fabs(det) > 1e-12 * tr * tr * tr)) return false;
  try {
    solve_dense_reference(A, x, 3);
  } catch (...) {
    return false;
  }
  a = x[0];
  b = x[1];
  c = x[2];
  return true;
}
}  // namespace

std::vector<ThrowEstimate> estimate_throws(const FittedPlane& pl, const std::vector<HorizonPoint>& pts, std::size_t nh, const ThrowOptions& o) {
  const Vec3 sdir{std::sin(pl.strike), std::cos(pl.strike), 0}, ddir{std::cos(pl.strike), -std::sin(pl.strike), 0};
  const Vec3 n = pl.normal;
  const double sind = std::max(1e-6, std::sin(pl.dip)), cotd = std::cos(pl.dip) / sind, cz = pl.centroid.z;
  std::vector<ThrowEstimate> out;
  for (std::size_t h = 0; h < nh; ++h) {
    ThrowEstimate e;
    e.horizon = static_cast<int>(h);
    // (s, u, z, side) of the data of this horizon within the band
    struct Q { double s, u, z; bool hw; };
    std::vector<Q> q;
    double smin = 1e300, smax = -1e300;
    for (auto& hp : pts) {
      if (hp.horizon != static_cast<int>(h)) continue;
      Vec3 d = hp.p - pl.centroid;
      double dist = dot(d, n), dh = dist / sind;  // horizontal distance to the surface at the point's level
      if (std::fabs(dh) > o.band) continue;
      q.push_back({dot(d, sdir), dot(d, ddir), hp.p.z, dist > 0});
      smin = std::min(smin, q.back().s);
      smax = std::max(smax, q.back().s);
    }
    if (q.size() < 6) {
      e.reason = "fewer than 6 data within the band";
      out.push_back(e);
      continue;
    }
    // stations and cut-offs
    for (double sk = smin; sk <= smax + 1e-9; sk += o.station_spacing) {
      std::vector<std::array<double, 3>> fw, hw;
      for (auto& x : q)
        if (std::fabs(x.s - sk) <= o.station_spacing) (x.hw ? hw : fw).push_back({x.s - sk, x.u, x.z});
      if (fw.size() < 3 || hw.size() < 3) continue;
      double z[2];
      bool ok = true;
      for (int side = 0; side < 2 && ok; ++side) {
        double a, b, c;
        if (!fit_plane3(side ? hw : fw, a, b, c)) { ok = false; break; }
        z[side] = (a + c * cz * cotd) / (1 + c * cotd);  // where the side's plane meets the fault surface at s_k
      }
      if (!ok) continue;
      e.stations.push_back({sk, z[0] - z[1], fw.size(), hw.size()});
    }
    if (e.stations.size() < 3) {
      e.reason = "fewer than 3 stations with data on both sides";
      out.push_back(e);
      continue;
    }
    // profile fit: T(s) = T_max * shape((s - s0) / L + 1/2); T_max linear for given (s0, L)
    auto eval = [&](double s0, double L, double& tmax) {
      double sg = 0, sgt = 0;
      std::vector<double> g(e.stations.size());
      for (std::size_t k = 0; k < g.size(); ++k) {
        g[k] = shape_value(o.profile, o.taper, (e.stations[k].s - s0) / L + 0.5);
        sg += g[k] * g[k];
        sgt += g[k] * e.stations[k].throw_m;
      }
      tmax = sg > 0 ? sgt / sg : 0.0;
      double ss = 0;
      for (std::size_t k = 0; k < g.size(); ++k) ss += std::pow(e.stations[k].throw_m - tmax * g[k], 2);
      return std::sqrt(ss / static_cast<double>(g.size()));
    };
    const double span = std::max(e.stations.back().s - e.stations.front().s, o.station_spacing);
    double best = 1e300, bs0 = 0, bL = span, bt = 0;
    for (int i = 0; i <= 60; ++i)
      for (int j = 1; j <= 60; ++j) {
        double s0 = e.stations.front().s + span * i / 60.0, L = 2 * o.station_spacing + (2 * span) * j / 60.0, t;
        double r = eval(s0, L, t);
        if (r < best) { best = r; bs0 = s0; bL = L; bt = t; }
      }
    // refine: coordinate golden-section searches
    double ds = span / 60.0, dL = 2 * span / 60.0;
    for (int round = 0; round < 6; ++round) {
      auto golden = [&](double lo, double hi, auto f) {
        const double g = 0.5 * (std::sqrt(5.0) - 1);
        double a1 = hi - g * (hi - lo), b1 = lo + g * (hi - lo), fa = f(a1), fb = f(b1);
        for (int it = 0; it < 40; ++it) {
          if (fa < fb) { hi = b1; b1 = a1; fb = fa; a1 = hi - g * (hi - lo); fa = f(a1); }
          else { lo = a1; a1 = b1; fa = fb; b1 = lo + g * (hi - lo); fb = f(b1); }
        }
        return 0.5 * (lo + hi);
      };
      double t;
      bs0 = golden(bs0 - ds, bs0 + ds, [&](double x) { return eval(x, bL, t); });
      bL = golden(std::max(o.station_spacing, bL - dL), bL + dL, [&](double x) { return eval(bs0, x, t); });
      best = eval(bs0, bL, bt);
      ds *= 0.5;
      dL *= 0.5;
    }
    e.max_throw = bt;
    e.centre_s = bs0;
    e.length = bL;
    e.rms = best;
    e.ok = true;
    out.push_back(e);
  }
  return out;
}

StickProfile profile_from_sticks(const std::vector<Polyline>& sticks) {
  FittedPlane fp = fit_plane(sticks);
  StickProfile pr;
  pr.strike = fp.strike;
  const Vec3 ddir{std::cos(fp.strike), -std::sin(fp.strike), 0};
  // each stick as u(depth), sorted by depth
  std::vector<std::vector<std::pair<double, double>>> st;
  double dtop = -1e300, dbot = 1e300;  // depth range common to every stick
  for (auto& s : sticks) {
    std::vector<std::pair<double, double>> v;
    for (auto& p : s) v.push_back({-p.z, dot(p - fp.centroid, ddir)});
    if (v.size() < 2) continue;
    std::sort(v.begin(), v.end());
    dtop = std::max(dtop, v.front().first);
    dbot = std::min(dbot, v.back().first);
    st.push_back(std::move(v));
  }
  if (st.empty() || !(dbot > dtop + 100)) throw Error(ErrorCode::Build, "fault profile from sticks: the sticks share less than 100 m of depth");
  auto u_at = [](const std::vector<std::pair<double, double>>& v, double d) {
    auto it = std::lower_bound(v.begin(), v.end(), std::make_pair(d, -1e300));
    if (it == v.begin()) return v.front().second;
    if (it == v.end()) return v.back().second;
    auto pv = it - 1;
    double t = (d - pv->first) / std::max(1e-9, it->first - pv->first);
    return pv->second + t * (it->second - pv->second);
  };
  const double bin = 50.0;  // m: depth bins of the mean profile
  std::vector<double> dz, uz;
  for (double d = dtop; d <= dbot + 1e-9; d += bin) {
    double s = 0;
    for (auto& v : st) s += u_at(v, d);
    dz.push_back(d);
    uz.push_back(s / static_cast<double>(st.size()));
  }
  for (std::size_t k = 0; k < dz.size(); ++k) {
    std::size_t a = k ? k - 1 : 0, b = std::min(k + 1, dz.size() - 1);
    double du = uz[b] - uz[a], dd = dz[b] - dz[a];
    pr.depth.push_back(dz[k]);
    pr.u.push_back(uz[k]);
    pr.dip.push_back(std::atan2(dd, std::max(1e-9, du)));
  }
  pr.centroid = fp.centroid;
  pr.top_point = fp.centroid + ddir * uz.front();
  pr.top_point.z = -dz.front();
  // RMS distance of the stick points to the mean surface (normal to the local dip)
  double ss = 0;
  std::size_t np = 0;
  for (auto& v : st)
    for (auto& [d, u] : v) {
      if (d < dtop || d > dbot) continue;
      std::size_t k = std::min(dz.size() - 1, static_cast<std::size_t>(std::lround((d - dtop) / bin)));
      double um = 0;
      for (auto& w : st) um += u_at(w, d);
      um /= static_cast<double>(st.size());
      ss += std::pow((u - um) * std::sin(pr.dip[k]), 2);
      ++np;
    }
  pr.rms = np ? std::sqrt(ss / static_cast<double>(np)) : 0.0;
  return pr;
}

namespace {
struct CondDatum {
  Vec3 p, n;
  int horizon = -1;
  bool oriented = false;
  double sigma = 0;  ///< m; degrees for orientations
  std::string input, evidence, kind;
};
}  // namespace

ConditioningResult condition_to_data(const desc::Description& d, const BuildInputs& in, BuildLog& log) {
  ConditioningResult R;
  const auto& cd = d.conditioning;
  // 1. data, each with its uncertainty
  std::vector<CondDatum> data;
  std::set<std::string> used_inputs;
  for (auto& e : cd.data) {
    auto it = in.data.find(e.input);
    if (it == in.data.end()) { R.reason = "conditioning data input '" + e.input + "' was not loaded"; return R; }
    std::string ev = "assumed";
    for (auto& x : d.inputs)
      if (x.id == e.input && !x.evidence.empty()) ev = x.evidence;
    const bool orient = e.kind == "orientations";
    for (auto& dp : it->second) {
      int h = d.horizon_index(dp.label);
      if (h < 0) continue;
      CondDatum c;
      c.p = dp.p;
      c.n = dp.n;
      c.horizon = h;
      c.oriented = dp.oriented;
      c.input = e.input;
      c.evidence = ev;
      c.kind = dp.oriented ? "orientation" : "point";
      if (dp.sigma > 0) c.sigma = dp.sigma;
      else if (e.uncertainty.present) c.sigma = orient ? rad2deg(e.uncertainty.si()) : e.uncertainty.si();
      else if (orient && cd.orientation_uncertainty.present) c.sigma = rad2deg(cd.orientation_uncertainty.si());
      else if (!orient && cd.uncertainty.present) c.sigma = cd.uncertainty.si();
      if (!(c.sigma > 0)) { R.reason = "datum of input " + e.input + " has no positive uncertainty (give conditioning.uncertainty or a sigma column)"; return R; }
      data.push_back(c);
      used_inputs.insert(e.input);
    }
  }
  if (cd.use_well_tops) {
    auto wells = make_well_paths(d);
    for (std::size_t w = 0; w < d.wells.size(); ++w) {
      auto& wd = d.wells[w];
      if (wd.role != "control") continue;
      for (auto& [hn, p] : wd.tops) {
        int h = d.horizon_index(hn);
        if (h < 0) continue;
        double z = -p.si();
        CondDatum c;
        c.p = wells[w].vertical() ? Vec3{wells[w].x(), wells[w].y(), z} : wells[w].at_elevation(z);
        c.horizon = h;
        c.input = "well:" + wd.name;
        c.evidence = desc::evidence_name(p.evidence);
        c.kind = "well_top";
        c.sigma = wd.top_uncertainty.present ? wd.top_uncertainty.si() : (cd.uncertainty.present ? cd.uncertainty.si() : 0.0);
        if (!(c.sigma > 0)) { R.reason = "top " + hn + " of well " + wd.name + " has no positive uncertainty (give top_uncertainty or conditioning.uncertainty)"; return R; }
        data.push_back(c);
        used_inputs.insert(c.input);
      }
    }
  }
  const std::size_t m = data.size(), np = cd.parameters.size();
  if (m < np) { R.reason = "fewer data (" + std::to_string(m) + ") than parameters (" + std::to_string(np) + ")"; return R; }
  // 2. parameters (SI internally)
  desc::Description dd = d;
  dd.conditioning = desc::ConditioningDesc{};
  // control wells stay data, not a correction of the surfaces: the residuals must show the fit
  for (auto& w : dd.wells)
    if (w.role == "control") w.role = "check";
  std::vector<desc::Param*> P;
  std::vector<double> x(np), lo(np), hi(np);
  for (std::size_t k = 0; k < np; ++k) {
    auto& fp = cd.parameters[k];
    desc::Param* p = dd.param_by_name(fp.path);
    if (!p) { R.reason = "'" + fp.path + "' names no parameter of this description"; return R; }
    P.push_back(p);
    x[k] = p->si();
    lo[k] = fp.min.si();
    hi[k] = fp.max.si();
  }
  std::vector<Vec3> pts;
  for (auto& c : data) pts.push_back(c.p);
  // forward model: normalised residuals at x; false if the model cannot be built there
  auto forward = [&](const std::vector<double>& xv, std::vector<double>& r, std::vector<double>* raw) {
    for (std::size_t k = 0; k < np; ++k) P[k]->values[0] = xv[k] / P[k]->to_si;
    BuildLog lq;
    std::unique_ptr<ParametricModel> mdl;
    try {
      mdl = std::make_unique<ParametricModel>(dd, in, lq);
    } catch (const std::exception&) {
      return false;
    }
    EvalResult er;
    mdl->evaluate(pts, kEvalScalars, er);
    r.assign(m, 0.0);
    if (raw) raw->assign(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
      auto& c = data[i];
      double v;
      if (!c.oriented) {
        v = er.s(static_cast<std::size_t>(c.horizon), i);
      } else {
        const double e = 2.0;
        std::vector<Vec3> q{c.p + Vec3{e, 0, 0}, c.p - Vec3{e, 0, 0}, c.p + Vec3{0, e, 0}, c.p - Vec3{0, e, 0}, c.p + Vec3{0, 0, e}, c.p - Vec3{0, 0, e}};
        EvalResult eq;
        mdl->evaluate(q, kEvalScalars, eq);
        auto s = [&](std::size_t k) { return eq.s(static_cast<std::size_t>(c.horizon), k); };
        Vec3 g = normalized(Vec3{s(0) - s(1), s(2) - s(3), s(4) - s(5)});
        v = rad2deg(std::acos(std::clamp(dot(g, c.n), -1.0, 1.0)));
      }
      if (!std::isfinite(v)) return false;
      r[i] = v / c.sigma;
      if (raw) (*raw)[i] = v;
    }
    return true;
  };
  auto chi2_of = [](const std::vector<double>& r) {
    double s = 0;
    for (double v : r) s += v * v;
    return s;
  };
  auto jacobian = [&](const std::vector<double>& xv, const std::vector<double>& r, std::vector<double>& J) {
    J.assign(m * np, 0.0);
    for (std::size_t k = 0; k < np; ++k) {
      double h = 1e-3 * (hi[k] - lo[k]);
      std::vector<double> xs = xv;
      if (xs[k] + h > hi[k]) h = -h;
      xs[k] += h;
      std::vector<double> rs;
      if (!forward(xs, rs, nullptr)) return false;
      for (std::size_t i = 0; i < m; ++i) J[i * np + k] = (rs[i] - r[i]) / h;
    }
    return true;
  };
  const int maxit = cd.max_iterations.present ? static_cast<int>(cd.max_iterations.si()) : 50;  // Stratonix default
  std::vector<double> r, J;
  if (!forward(x, r, nullptr)) { R.reason = "the model cannot be built at the start values"; return R; }
  double chi2 = chi2_of(r), lambda = 1e-3;
  bool converged = false;
  int it = 0;
  for (; it < maxit && !converged; ++it) {
    if (!jacobian(x, r, J)) { R.reason = "the model cannot be built near the current values"; return R; }
    std::vector<double> A(np * np, 0.0), g(np, 0.0);
    for (std::size_t i = 0; i < m; ++i)
      for (std::size_t a = 0; a < np; ++a) {
        g[a] += J[i * np + a] * r[i];
        for (std::size_t b = 0; b < np; ++b) A[a * np + b] += J[i * np + a] * J[i * np + b];
      }
    // a parameter the data do not control (its column of J vanishes)
    double amax = 0;
    for (std::size_t a = 0; a < np; ++a) amax = std::max(amax, A[a * np + a]);
    for (std::size_t a = 0; a < np; ++a)
      if (!(A[a * np + a] > 1e-12 * amax) || A[a * np + a] == 0) {
        R.reason = "parameter " + cd.parameters[a].path + " is not controlled by the data (no datum is sensitive to it)";
        R.iterations = it;
        return R;
      }
    for (int tries = 0; tries < 20; ++tries) {
      std::vector<double> M = A, dx(np);
      for (std::size_t a = 0; a < np; ++a) {
        M[a * np + a] *= 1 + lambda;
        dx[a] = -g[a];
      }
      try {
        solve_dense_reference(M, dx, np);
      } catch (...) {
        lambda *= 10;
        continue;
      }
      std::vector<double> xn(np);
      double maxstep = 0;
      for (std::size_t a = 0; a < np; ++a) {
        xn[a] = std::clamp(x[a] + dx[a], lo[a], hi[a]);
        maxstep = std::max(maxstep, std::fabs(xn[a] - x[a]) / (hi[a] - lo[a]));
      }
      std::vector<double> rn;
      if (maxstep < 1e-6) { converged = true; break; }  // no further step (at a minimum or held at a bound)
      if (forward(xn, rn, nullptr) && chi2_of(rn) <= chi2) {
        double c2 = chi2_of(rn);
        converged = std::fabs(chi2 - c2) <= 1e-8 * std::max(1.0, chi2) || maxstep < 1e-6;
        x = xn;
        r = rn;
        chi2 = c2;
        lambda = std::max(lambda / 10, 1e-9);
        break;
      }
      lambda *= 10;
      if (lambda > 1e12) { converged = true; break; }  // no descent direction left: a minimum
    }
  }
  R.iterations = it;
  if (!converged) {
    R.reason = "the fit did not converge within " + std::to_string(maxit) + " iterations (reduced chi2 " +
               fmt_fixed(chi2 / std::max<double>(1.0, static_cast<double>(m - np)), 3) + ")";
    return R;
  }
  for (std::size_t k = 0; k < np; ++k) {
    double rng = hi[k] - lo[k];
    if (x[k] - lo[k] < 1e-6 * rng || hi[k] - x[k] < 1e-6 * rng) {
      R.reason = "parameter " + cd.parameters[k].path + " converged on the bound of its admissible range (" + fmt_shortest(x[k] / P[k]->to_si) + " " +
                 P[k]->unit + "): the best fit lies outside the range";
      return R;
    }
  }
  // 3. uncertainty, controlling data, leverage at the solution
  std::vector<double> raw;
  forward(x, r, &raw);
  if (!jacobian(x, r, J)) { R.reason = "the model cannot be built near the solution"; return R; }
  std::vector<double> A(np * np, 0.0), Cov(np * np, 0.0);
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t a = 0; a < np; ++a)
      for (std::size_t b = 0; b < np; ++b) A[a * np + b] += J[i * np + a] * J[i * np + b];
  for (std::size_t a = 0; a < np; ++a) Cov[a * np + a] = 1.0;
  try {
    std::vector<double> Ac = A;
    solve_dense_multi(Ac, Cov, np, np);
  } catch (...) {
    R.reason = "the parameters are not independently controlled by the data (J^T J is singular)";
    return R;
  }
  R.chi2 = chi2;
  R.chi2_reduced = m > np ? chi2 / static_cast<double>(m - np) : kNaN;
  std::string inputs_txt;
  for (auto& s2 : used_inputs) inputs_txt += (inputs_txt.empty() ? "" : ", ") + s2;
  const std::string src = "fitted to " + inputs_txt + " (" + std::to_string(m) + " data, Levenberg-Marquardt, reduced chi2 " +
                          fmt_fixed(R.chi2_reduced, 3) + ")";
  for (std::size_t k = 0; k < np; ++k) {
    FittedParameter fp;
    fp.path = cd.parameters[k].path;
    fp.unit = P[k]->unit;
    fp.start = d.param_by_name(fp.path)->values[0];
    fp.value = x[k] / P[k]->to_si;
    fp.sigma = std::sqrt(std::max(0.0, Cov[k * np + k])) / P[k]->to_si;
    fp.min = lo[k] / P[k]->to_si;
    fp.max = hi[k] / P[k]->to_si;
    // sensitivity of the parameter to each normalised datum: G = Cov J^T
    std::vector<std::pair<double, std::size_t>> gk;
    double tot = 0;
    for (std::size_t i = 0; i < m; ++i) {
      double gv = 0;
      for (std::size_t b = 0; b < np; ++b) gv += Cov[k * np + b] * J[i * np + b];
      gk.push_back({gv * gv, i});
      tot += gv * gv;
    }
    std::sort(gk.begin(), gk.end(), [](auto& a2, auto& b2) { return a2.first > b2.first; });
    for (std::size_t q = 0; q < std::min<std::size_t>(5, gk.size()); ++q) fp.controlling.push_back({gk[q].second, tot > 0 ? gk[q].first / tot : 0.0});
    // label the fitted value
    P[k]->values[0] = fp.value;
    P[k]->evidence = desc::Evidence::Calibrated;
    P[k]->source = src;
    P[k]->text = fmt_fixed(fp.value, 3) + " " + fp.unit + " | calibrated | " + src;
    log.value("conditioning." + fp.path, fmt_fixed(fp.value, 3) + " +- " + fmt_fixed(fp.sigma, 3) + " " + fp.unit + " (start " + fmt_shortest(fp.start) +
                                             ", range " + fmt_shortest(fp.min) + " to " + fmt_shortest(fp.max) + ")");
    R.parameters.push_back(fp);
  }
  for (std::size_t i = 0; i < m; ++i) {
    double h = 0;
    for (std::size_t a = 0; a < np; ++a)
      for (std::size_t b = 0; b < np; ++b) h += J[i * np + a] * Cov[a * np + b] * J[i * np + b];
    R.leverage.push_back(h);
    auto& c = data[i];
    DataMisfit dm{c.input, c.evidence, d.horizons[static_cast<std::size_t>(c.horizon)].name, c.kind, c.p, raw[i]};
    dm.sigma = c.sigma;
    dm.normalised = raw[i] / c.sigma;
    R.residuals.push_back(dm);
  }
  log.value("conditioning.data", std::to_string(m) + " data from " + inputs_txt);
  log.value("conditioning.iterations", std::to_string(R.iterations));
  log.value("conditioning.chi2_reduced", fmt_fixed(R.chi2_reduced, 4));
  R.fitted = dd;
  R.fitted.conditioning = d.conditioning;  // keep the block, so checks and reports know the fit
  R.converged = true;
  return R;
}

ConditionedModel::ConditionedModel(std::shared_ptr<ParametricModel> inner, ConditioningResult r)
    : inner_(std::move(inner)), r_(std::move(r)), info_(inner_->info()) {}
const ModelInfo& ConditionedModel::info() const { return info_; }
void ConditionedModel::evaluate(std::span<const Vec3> pts, unsigned fields, EvalResult& out) const { inner_->evaluate(pts, fields, out); }
void ConditionedModel::scalar_in_block(std::span<const Vec3> pts, BlockCode block, int horizon, std::span<double> out) const {
  inner_->scalar_in_block(pts, block, horizon, out);
}
void ConditionedModel::evaluate_fault(int fault, std::span<const Vec3> pts, FaultEval& out) const { inner_->evaluate_fault(fault, pts, out); }
void ConditionedModel::restore(std::span<const Vec3> pts, std::span<Vec3> out) const { inner_->restore(pts, out); }
const void* ConditionedModel::capability(std::string_view name) const {
  if (name == "conditioning") return &r_;
  if (name == "data_misfit") return &r_.residuals;
  return inner_->capability(name);
}

ModelPtr build_conditioned_model(const desc::Description& d, const BuildInputs& in, BuildLog& log) {
  auto r = condition_to_data(d, in, log);
  if (!r.converged) throw Error(ErrorCode::Build, "conditioning refused: " + r.reason);
  desc::Description fd = r.fitted;
  fd.conditioning = desc::ConditioningDesc{};
  auto inner = std::make_shared<ParametricModel>(fd, in, log);
  return std::make_shared<ConditionedModel>(inner, std::move(r));
}

}  // namespace stratonix
