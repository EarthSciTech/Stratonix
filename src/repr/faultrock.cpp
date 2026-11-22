// Faults as rock and the stress field (stage 2F, STX-TN-10).
#include "stratonix/repr/faultrock.hpp"

#include <algorithm>
#include <cmath>

#include "stratonix/core/linalg.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/repr/repr.hpp"

namespace stratonix::repr {

namespace {
void add_missing(FaultRockPoint& p, const std::string& why) {
  if (p.missing.find(why) == std::string::npos) p.missing += (p.missing.empty() ? "" : "; ") + why;
}
/// Unit boundaries of a column as depths below its top horizon (b[h] for horizon h).
std::vector<double> boundaries(const ColumnQuery& cq, std::size_t nh, double x, double y) {
  std::vector<double> z(nh), b(nh);
  cq.column(x, y, z.data());
  for (std::size_t h = 0; h < nh; ++h) b[h] = z[0] - z[h];
  return b;
}
/// Model unit holding stratigraphic depth d (unit k between horizons k-1 and k; 0 above, nh below).
int unit_at(const std::vector<double>& b, double d) {
  int u = 0;
  for (std::size_t h = 0; h < b.size(); ++h)
    if (d > b[h]) u = static_cast<int>(h) + 1;
  return u;
}
}  // namespace

void stress_at(const StressInput& s, int unit, double z, double S6[6], double& pp) {
  for (int k = 0; k < 6; ++k) S6[k] = s.g[k] * z;
  pp = s.pp * z;
  if (!s.lateral || unit < 0 || static_cast<std::size_t>(unit) >= s.unit_mech.size()) return;
  const auto& m = s.unit_mech[static_cast<std::size_t>(unit)];
  const double E = m[0], nu = m[1], al = std::isfinite(m[2]) ? m[2] : 1.0;
  if (!std::isfinite(E) || !std::isfinite(nu)) return;
  // uniaxial strain plus tectonic strain (Thiercelin and Plumb 1994), sigma_v from the regime
  const double sv = S6[2], base = nu / (1 - nu) * (sv - al * pp) + al * pp, c = E / (1 - nu * nu);
  const double sH = base + c * (s.e_H + nu * s.e_h), sh = base + c * (s.e_h + nu * s.e_H);
  const double sn = std::sin(s.azimuth), cs = std::cos(s.azimuth);  // SHmax along (sin, cos)
  S6[0] = sH * sn * sn + sh * cs * cs;
  S6[1] = sH * cs * cs + sh * sn * sn;
  S6[3] = (sH - sh) * sn * cs;
  S6[4] = S6[5] = 0;
}

FaultRockPoint fault_rock_at(const Model& m, int fault, const Vec3& p, const FaultRockInput& in, const StressInput& s) {
  FaultRockPoint out;
  const auto& inf = m.info();
  const std::size_t nh = inf.horizons.size();
  // normal (towards the hanging wall) from the signed distance to the fault
  const double e = 0.5;
  std::vector<Vec3> q{p + Vec3{e, 0, 0}, p - Vec3{e, 0, 0}, p + Vec3{0, e, 0}, p - Vec3{0, e, 0}, p + Vec3{0, 0, e}, p - Vec3{0, 0, e}, p};
  FaultEval fe;
  m.evaluate_fault(fault, q, fe);
  Vec3 n{(fe.distance[0] - fe.distance[1]) / (2 * e), (fe.distance[2] - fe.distance[3]) / (2 * e), (fe.distance[4] - fe.distance[5]) / (2 * e)};
  if (norm(n) < 1e-9) { add_missing(out, "no fault normal"); return out; }
  n = normalized(n);
  const double displacement = fe.displacement[6];
  const double sin_dip = std::max(1e-6, std::hypot(n.x, n.y));
  // units and restored stratigraphic positions on each side
  std::vector<Vec3> sides{p - n * 0.25, p + n * 0.25}, rs(2);
  EvalResult er;
  m.evaluate(sides, kEvalUnit, er);
  out.fw_unit = er.unit[0];
  out.hw_unit = er.unit[1];
  m.restore(sides, rs);
  auto* cq = static_cast<const ColumnQuery*>(m.capability("column"));
  if (!cq || nh == 0) {
    add_missing(out, "no stratigraphic column (parametric engine only)");
  } else {
    auto bfw = boundaries(*cq, nh, rs[0].x, rs[0].y), bhw = boundaries(*cq, nh, rs[1].x, rs[1].y);
    std::vector<double> zf(nh), zh(nh);
    cq->column(rs[0].x, rs[0].y, zf.data());
    cq->column(rs[1].x, rs[1].y, zh.data());
    out.d_fw = zf[0] - rs[0].z;
    out.d_hw = zh[0] - rs[1].z;
    out.separation = out.d_fw - out.d_hw;
    out.fw_unit_col = unit_at(bfw, out.d_fw);
    out.hw_unit_col = unit_at(bhw, out.d_hw);
    for (std::size_t h = 0; h < nh; ++h)
      if (std::fabs(out.d_fw - bfw[h]) < 0.5 || std::fabs(out.d_hw - bhw[h]) < 0.5) out.near_horizon = true;
    auto clay = [&](int u) { return u >= 0 && static_cast<std::size_t>(u) < in.unit_clay.size() ? in.unit_clay[static_cast<std::size_t>(u)] : kNaN; };
    // shale gouge ratio over the window between the juxtaposed positions (footwall column)
    if (std::fabs(out.separation) < 1.0) {
      add_missing(out, "SGR: separation below 1 m");
    } else {
      const double lo = std::min(out.d_fw, out.d_hw), hi = std::max(out.d_fw, out.d_hw);
      std::vector<double> cuts{lo};
      for (double v : bfw)
        if (v > lo && v < hi) cuts.push_back(v);
      cuts.push_back(hi);
      std::sort(cuts.begin(), cuts.end());
      double sum = 0;
      bool ok = true;
      for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
        const double v = clay(unit_at(bfw, 0.5 * (cuts[k] + cuts[k + 1])));
        if (std::isnan(v)) { ok = false; break; }
        sum += v * (cuts[k + 1] - cuts[k]);
      }
      if (ok) out.sgr = sum / (hi - lo);
      else add_missing(out, "SGR: a unit in the window has no clay fraction");
    }
    // clay smear potential from both walls
    if (std::isnan(in.shale_threshold)) {
      add_missing(out, "CSP: no shale threshold");
    } else {
      double csp = 0;
      for (int side = 0; side < 2; ++side) {
        const auto& b = side == 0 ? bfw : bhw;
        const double d = side == 0 ? out.d_fw : out.d_hw;
        for (std::size_t u = 1; u < nh; ++u) {
          const double v = clay(static_cast<int>(u));
          if (std::isnan(v) || v < in.shale_threshold) continue;
          const double top = b[u - 1], base = b[u], dz = base - top;
          if (dz <= 0) continue;
          const double sv = d < top ? top - d : (d > base ? d - base : 0.0);
          csp += dz * dz / std::max(1.0, sv / sin_dip);
        }
      }
      out.csp = csp;
    }
  }
  // fault-zone thickness and seal capacity
  if (std::isnan(in.thickness_ratio)) add_missing(out, "thickness: no thickness-displacement ratio");
  else out.thickness = in.thickness_ratio * displacement;
  if (!in.seal) add_missing(out, "seal capacity: no calibration");
  else if (!std::isnan(out.sgr)) out.seal_bar = std::pow(10.0, 100.0 * out.sgr / in.seal_a - in.seal_c);
  // stress and fault response
  if (s.present) {
    const double z = std::max(0.0, -p.z);
    double S6[6], pp;
    stress_at(s, out.fw_unit, z, S6, pp);  // the footwall's unit carries the stress at the fault (stage 4B)
    const double S[9] = {S6[0], S6[3], S6[4], S6[3], S6[1], S6[5], S6[4], S6[5], S6[2]};
    const Vec3 t{S[0] * n.x + S[1] * n.y + S[2] * n.z, S[3] * n.x + S[4] * n.y + S[5] * n.z, S[6] * n.x + S[7] * n.y + S[8] * n.z};
    out.sigma_n = dot(t, n);
    out.tau = norm(t - n * out.sigma_n);
    const double sne = out.sigma_n - pp;
    double ev[3], evec[9];
    symmetric_eigen3(S, ev, evec);  // ascending
    const double s1 = ev[2], s3 = ev[0];
    out.ts = sne > 0 ? out.tau / sne : kNaN;
    out.td = s1 > s3 ? (s1 - out.sigma_n) / (s1 - s3) : kNaN;
    const double s1e = s1 - pp, s3e = s3 - pp;
    out.ts_max = s3e > 0 ? (s1e - s3e) / (2 * std::sqrt(s1e * s3e)) : kNaN;
    if (std::isnan(s.friction)) add_missing(out, "reactivation: no friction");
    else out.dp = sne - (out.tau - s.cohesion) / s.friction;
  } else {
    add_missing(out, "tendencies: no stress");
  }
  return out;
}

std::vector<FaultRockPoint> fault_rock_surface(const Model& m, const FaultSurface& fs, const FaultRockInput& in, const StressInput& s,
                                               double max_extent) {
  std::vector<FaultRockPoint> out(fs.nodes.size());
  const Domain& dm = m.info().domain;
  for (std::size_t k = 0; k < fs.nodes.size(); ++k) {
    const Vec3& p = fs.nodes[k];
    if (!(fs.extent[k] <= max_extent) || !dm.inside_xy(p.x, p.y) || p.z > dm.z_top || p.z < dm.z_base || std::isnan(p.z)) {
      out[k].valid = false;
      continue;
    }
    out[k] = fault_rock_at(m, fs.fault, p, in, s);
  }
  return out;
}

CurvatureGrids curvature(const Model& m, int horizon, const GridSpec& g, double thr) {
  CurvatureGrids c{Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g)};
  GridSpec gg{g.x0 - g.dx, g.y0 - g.dy, g.dx, g.dy, g.nx + 2, g.ny + 2};
  auto fs = horizon_grid_faulted(m, horizon, gg);
  const double dx = g.dx, dy = g.dy;
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i) {
      double z[3][3];
      bool ok = true;
      BlockCode b0 = fs.block[gg.index(i + 1, j + 1)];
      for (int a = -1; a <= 1 && ok; ++a)
        for (int bb = -1; bb <= 1 && ok; ++bb) {
          const std::size_t k = gg.index(static_cast<std::size_t>(static_cast<long>(i) + 1 + a), static_cast<std::size_t>(static_cast<long>(j) + 1 + bb));
          z[a + 1][bb + 1] = fs.z.v[k];
          ok = !std::isnan(fs.z.v[k]) && fs.block[k] == b0;
        }
      if (!ok) continue;
      // derivatives of elevation (central differences on the 3 x 3 window)
      const double zx = (z[2][1] - z[0][1]) / (2 * dx), zy = (z[1][2] - z[1][0]) / (2 * dy);
      const double zxx = (z[2][1] - 2 * z[1][1] + z[0][1]) / (dx * dx), zyy = (z[1][2] - 2 * z[1][1] + z[1][0]) / (dy * dy);
      const double zxy = (z[2][2] - z[0][2] - z[2][0] + z[0][0]) / (4 * dx * dy);
      // Roberts (2001) with convex upwards positive: a = -zxx/2, b = -zyy/2, c = -zxy, d = zx, e = zy
      const double A = -0.5 * zxx, B = -0.5 * zyy, C = -zxy, D = zx, E = zy;
      const double q = 1 + D * D + E * E;
      const double km = (A * (1 + E * E) + B * (1 + D * D) - C * D * E) / std::pow(q, 1.5);
      const double kg = (4 * A * B - C * C) / (q * q);
      const double disc = std::sqrt(std::max(0.0, km * km - kg));
      const std::size_t k = g.index(i, j);
      c.mean.v[k] = km;
      c.gaussian.v[k] = kg;
      c.kmax.v[k] = km + disc;
      c.kmin.v[k] = km - disc;
      c.most_pos.v[k] = (A + B) + std::sqrt((A - B) * (A - B) + C * C);
      c.most_neg.v[k] = (A + B) - std::sqrt((A - B) * (A - B) + C * C);
      const double s2 = D * D + E * E;
      if (s2 > 1e-14) {
        c.dip.v[k] = 2 * (A * D * D + B * E * E + C * D * E) / (s2 * std::pow(q, 1.5));
        c.strike.v[k] = 2 * (A * E * E + B * D * D - C * D * E) / (s2 * std::sqrt(q));
      }
      // fracture trend: the principal direction of the smaller curvature (slopes small)
      if (std::max(std::fabs(c.kmax.v[k]), std::fabs(c.kmin.v[k])) >= thr) {
        const double hxx = 2 * A, hyy = 2 * B, hxy = C;
        const double tr = 0.5 * (hxx + hyy), dd = std::sqrt(0.25 * (hxx - hyy) * (hxx - hyy) + hxy * hxy);
        const double l1 = tr + dd, l2 = tr - dd, lsmall = std::fabs(l1) < std::fabs(l2) ? l1 : l2;
        double vx = hxy, vy = lsmall - hxx;
        if (std::hypot(vx, vy) < 1e-18) { vx = lsmall - hyy; vy = hxy; }
        if (std::hypot(vx, vy) < 1e-18) { vx = 0; vy = 1; }
        double az = std::atan2(vx, vy) * 180.0 / kPi;
        if (az < 0) az += 180.0;
        if (az >= 180.0) az -= 180.0;
        c.trend.v[k] = az;
      }
    }
  return c;
}

}  // namespace stratonix::repr
