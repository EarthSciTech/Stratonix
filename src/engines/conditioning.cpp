// Stage 3A (STX-TN-11): faults from data and hybrid conditioning. Specified; built by Prompts 35
// and 36.
#include "stratonix/engines/conditioning.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include "stratonix/core/error.hpp"
#include "stratonix/core/linalg.hpp"

namespace stratonix {

namespace {
[[noreturn]] void not_built(const char* what, const char* prompt) {
  throw Error(ErrorCode::Build, std::string(what) + ": specified in STX-TN-11, not built yet (" + prompt + ")");
}
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

ConditioningResult condition_to_data(const desc::Description&, const BuildInputs&, BuildLog&) {
  not_built("hybrid conditioning", "Prompt 36");
}

}  // namespace stratonix
