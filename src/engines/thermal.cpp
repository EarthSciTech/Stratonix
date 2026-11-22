// Stage 4B (STX-TN-15 section 6, decision 20): temperature and overpressure along 1D columns.
#include "stratonix/engines/thermal.hpp"

#include <algorithm>
#include <cmath>

#include "stratonix/core/error.hpp"
#include "stratonix/engines/structure.hpp"

namespace stratonix {

std::vector<double> steady_temperature(const std::vector<ThermalLayer>& layers, double t_surface, double q_surface) {
  // downward through each layer: T(z) = T0 + q z / k - A z^2 / (2 k), q decreasing by A h
  std::vector<double> T{t_surface};
  double t = t_surface, q = q_surface;
  for (auto& L : layers) {
    if (!(L.conductivity > 0) || L.thickness < 0) throw Error(ErrorCode::InvalidArgument, "steady conduction: positive conductivity and thickness");
    t += q * L.thickness / L.conductivity - L.heat_production * L.thickness * L.thickness / (2 * L.conductivity);
    q -= L.heat_production * L.thickness;
    T.push_back(t);
  }
  return T;
}

namespace {
/// Tridiagonal solve (Thomas), in place on d.
void thomas(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double>& d) {
  const std::size_t n = d.size();
  for (std::size_t i = 1; i < n; ++i) {
    const double m = a[i] / b[i - 1];
    b[i] -= m * c[i - 1];
    d[i] -= m * d[i - 1];
  }
  d[n - 1] /= b[n - 1];
  for (std::size_t i = n - 1; i-- > 0;) d[i] = (d[i] - c[i] * d[i + 1]) / b[i];
}
}  // namespace

std::vector<double> transient_temperature(const TransientSpec& s) {
  if (s.nodes < 3 || !(s.length > 0) || !(s.diffusivity > 0) || s.time < 0) throw Error(ErrorCode::InvalidArgument, "transient conduction: a column, diffusivity and time");
  const std::size_t n = static_cast<std::size_t>(s.nodes);
  const double dz = s.length / (n - 1);
  std::vector<double> T(n, s.t_initial);
  T.front() = s.t_top;
  T.back() = s.t_base;
  // Crank-Nicolson, the first steps implicit Euler (Rannacher start) against the step's oscillation
  const int steps = 2000;
  const double dt = s.time / steps, r = s.diffusivity * dt / (dz * dz);
  for (int it = 0; it < steps; ++it) {
    const double th = it < 4 ? 1.0 : 0.5;
    std::vector<double> a(n, 0), b(n, 1), c(n, 0), d(n);
    for (std::size_t i = 1; i + 1 < n; ++i) {
      a[i] = -th * r;
      b[i] = 1 + 2 * th * r;
      c[i] = -th * r;
      d[i] = T[i] + (1 - th) * r * (T[i - 1] - 2 * T[i] + T[i + 1]);
    }
    d[0] = s.t_top;
    d[n - 1] = s.t_base;
    thomas(a, b, c, d);
    T = d;
  }
  return T;
}

double gibson_overpressure_base(const GibsonSpec& s) {
  if (!(s.rate > 0) || !(s.duration > 0) || !(s.permeability > 0) || !(s.viscosity > 0) || !(s.compressibility > 0))
    throw Error(ErrorCode::InvalidArgument, "Gibson consolidation: positive rate, duration, permeability, viscosity and compressibility");
  // Excess pressure u in a layer growing at rate m on an impermeable base (Gibson 1958):
  //   du/dt = c_v d2u/dz2 + gamma' m,  c_v = k / (mu m_v),  u = 0 at the surface, du/dz = 0 at the base.
  // In zeta = z / h(t) (0 at the base, 1 at the surface): du/dt = (zeta m / h) du/dzeta + c_v / h^2 d2u/dzeta2 + gamma' m.
  const double cv = s.permeability / (s.viscosity * s.compressibility), gp = (s.rho_s - s.rho_w) * s.g;
  const std::size_t n = static_cast<std::size_t>(std::max(20, s.nodes));
  const double dzeta = 1.0 / (n - 1);
  std::vector<double> u(n, 0.0);
  // logarithmic time steps from a thin initial layer (its pressure is negligible)
  const double t0 = s.duration * 1e-6;
  const int steps = 4000;
  const double f = std::pow(s.duration / t0, 1.0 / steps);
  double t = t0;
  for (int it = 0; it < steps; ++it) {
    const double t1 = t * f, dt = t1 - t, h = s.rate * t1;
    const double D = cv / (h * h), v = s.rate / h;  // diffusion and the zeta-advection coefficient
    std::vector<double> a(n, 0), b(n, 1), c(n, 0), d(n);
    for (std::size_t i = 0; i + 1 < n; ++i) {
      const double z = static_cast<double>(i) * dzeta;
      // implicit: (u_i' - u_i)/dt = v z (u_{i+1}' - u_i')/dzeta + D (u_{i-1}' - 2 u_i' + u_{i+1}')/dzeta^2 + gamma' m
      const double adv = v * z / dzeta, dif = D / (dzeta * dzeta);
      if (i == 0) {
        b[i] = 1 + dt * 2 * dif;  // mirror node at the impermeable base
        c[i] = -dt * 2 * dif;
      } else {
        a[i] = -dt * dif;
        b[i] = 1 + dt * (2 * dif + adv);
        c[i] = -dt * (dif + adv);
      }
      d[i] = u[i] + dt * gp * s.rate;
    }
    d[n - 1] = 0;
    b[n - 1] = 1;
    a[n - 1] = 0;
    thomas(a, b, c, d);
    u = d;
    t = t1;
  }
  return u[0];
}

// ---- columns through the burial history ---------------------------------------------------------
namespace {
double prop(const std::vector<double>& v, int u) {
  if (v.empty()) return kNaN;
  if (u < 0) u = 0;
  for (int k = std::min<int>(u, static_cast<int>(v.size()) - 1); k >= 0; --k)
    if (std::isfinite(v[static_cast<std::size_t>(k)])) return v[static_cast<std::size_t>(k)];
  for (double x : v)
    if (std::isfinite(x)) return x;
  return kNaN;
}
/// Unit at a depth in one step of the burial history (depths per horizon, NaN before deposition).
int unit_at_depth(const std::vector<double>& depth, double z) {
  int top = -1;
  for (std::size_t h = 0; h < depth.size(); ++h) {
    if (!std::isfinite(depth[h])) continue;
    if (top < 0) top = static_cast<int>(h);
    if (z < depth[h]) return static_cast<int>(h);  // above horizon h: unit h (between h - 1 and h)
  }
  return static_cast<int>(depth.size());  // below the deepest horizon
}
/// Linear interpolation of two burial steps (fraction w of the way from a to b).
std::vector<double> lerp_depths(const std::vector<double>& a, const std::vector<double>& b, double w) {
  std::vector<double> o(a.size(), kNaN);
  for (std::size_t h = 0; h < a.size(); ++h) {
    if (std::isfinite(a[h]) && std::isfinite(b[h])) o[h] = a[h] + w * (b[h] - a[h]);
    else if (std::isfinite(b[h]) && w > 0.999999) o[h] = b[h];
    else if (std::isfinite(b[h]) && b[h] == 0) o[h] = 0;  // a horizon being deposited stays at the surface
  }
  return o;
}
constexpr double kMa = 1e6 * 365.25 * 86400;
}  // namespace

double easy_ro(const std::vector<std::pair<double, double>>& path) {
  // 20 parallel first-order reactions, A = 1e13 /s, E = 34..72 kcal/mol (Sweeney and Burnham 1990)
  static const double f[20] = {0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.06, 0.04, 0.04, 0.07, 0.06, 0.06, 0.06, 0.05, 0.05, 0.04, 0.03, 0.02, 0.02, 0.01};
  const double A = 1e13, R = 1.987e-3;  // kcal/mol/K
  double I[20] = {};
  for (std::size_t k = 1; k < path.size(); ++k) {
    const double dt = (path[k - 1].first - path[k].first) * kMa;  // ages decrease
    if (!(dt > 0)) continue;
    const int ns = 20;
    for (int s = 0; s < ns; ++s) {
      const double w = (s + 0.5) / ns, T = path[k - 1].second + w * (path[k].second - path[k - 1].second) + 273.15;
      for (int i = 0; i < 20; ++i) I[i] += A * std::exp(-(34 + 2 * i) / (R * T)) * dt / ns;
    }
  }
  double F = 0;
  for (int i = 0; i < 20; ++i) F += f[i] * (1 - std::exp(-I[i]));
  return std::exp(-1.6 + 3.7 * F);
}

ThermalHistory thermal_history(const Structure& st, const ThermalColumnSpec& s, double x, double y) {
  ThermalHistory H;
  auto hist = st.burial_history(x, y);
  if (hist.size() < 2) return H;
  const std::size_t nh = hist.back().depth.size();
  double zmax = 0;
  for (double z : hist.back().depth)
    if (std::isfinite(z)) zmax = std::max(zmax, z);
  const double L = s.column_depth > 0 ? s.column_depth : std::max(5000.0, 2 * zmax);
  const std::size_t n = 401;
  const double dz = L / (n - 1);
  auto qb = [&](double age) {
    if (s.q_history.empty() || !std::isfinite(age)) return s.q_base;
    auto q = s.q_history;
    std::sort(q.begin(), q.end());
    if (age <= q.front().first) return q.front().second;
    if (age >= q.back().first) return q.back().second;
    for (std::size_t k = 1; k < q.size(); ++k)
      if (age <= q[k].first) return q[k - 1].second + (age - q[k - 1].first) / (q[k].first - q[k - 1].first) * (q[k].second - q[k - 1].second);
    return s.q_base;
  };
  // node properties for a step: conductivity (geometric mean of grain and water by porosity) and
  // radiogenic heat (grains only)
  auto props = [&](const std::vector<double>& depth, std::vector<double>& k, std::vector<double>& A) {
    k.assign(n, 0);
    A.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
      const double z = i * dz;
      const int u = unit_at_depth(depth, z);
      double phi = st.mean_porosity(u, z, z);
      if (!std::isfinite(phi)) phi = 0;
      const double kg = prop(s.props.conductivity, u), ag = prop(s.props.heat_production, u);
      k[i] = std::pow(std::isfinite(kg) ? kg : 2.5, 1 - phi) * std::pow(s.k_water, phi);
      A[i] = (std::isfinite(ag) ? ag : 0.0) * (1 - phi);
    }
  };
  // steady column for given properties: from the surface down, q(z) = q_s - int A (finite volumes)
  auto steady = [&](const std::vector<double>& k, const std::vector<double>& A, double q_base, std::vector<double>& T) {
    double qs = q_base;
    for (std::size_t i = 0; i < n; ++i) qs += A[i] * dz * (i == 0 || i + 1 == n ? 0.5 : 1.0);
    T.assign(n, s.t_surface);
    double q = qs;
    for (std::size_t i = 1; i < n; ++i) {
      const double kk = 2 / (1 / k[i - 1] + 1 / k[i]);
      q -= 0.5 * (A[i - 1] + A[i]) * dz;
      T[i] = T[i - 1] + (q + 0.5 * (A[i - 1] + A[i]) * dz * 0.5) * dz / kk;  // midpoint flux of the cell
    }
    return qs;
  };
  std::vector<double> k, A, T;
  props(hist.front().depth, k, A);
  steady(k, A, qb(hist.front().age), T);
  const double kappa_c = 1e-6;  // thermal diffusivity scale for the step size (m2/s)
  std::vector<std::vector<std::pair<double, double>>> paths(nh);
  auto record = [&](const std::vector<double>& depth, double age, std::size_t step) {
    std::vector<double> th(nh, kNaN);
    for (std::size_t h = 0; h < nh; ++h) {
      if (!std::isfinite(depth[h])) continue;
      const double fz = std::clamp(depth[h] / dz, 0.0, static_cast<double>(n - 1));
      const std::size_t i = std::min(n - 2, static_cast<std::size_t>(fz));
      th[h] = T[i] + (fz - i) * (T[i + 1] - T[i]);
      if (std::isfinite(age)) paths[h].push_back({age, th[h]});
    }
    (void)step;
    return th;
  };
  for (std::size_t sidx = 0; sidx < hist.size(); ++sidx) {
    const auto& cur = hist[sidx];
    if (sidx > 0) {
      const auto& prev = hist[sidx - 1];
      const double a0 = prev.age, a1 = cur.age;
      const double dur = std::isfinite(a0) && std::isfinite(a1) ? (a0 - a1) * kMa : 0.0;
      if (dur > 0) {
        // implicit conduction, substeps of at most the diffusion time of 10 nodes, burial moving
        const int ns = std::clamp(static_cast<int>(dur / (100 * dz * dz / kappa_c)) + 1, 10, 400);
        const double dt = dur / ns;
        for (int it = 1; it <= ns; ++it) {
          const double w = static_cast<double>(it) / ns, age = a0 + w * (a1 - a0);
          const auto depth = lerp_depths(prev.depth, cur.depth, w);
          props(depth, k, A);
          const double rc = 2.4e6;  // volumetric heat capacity (J/m3/K), sediment
          std::vector<double> a(n, 0), b(n, 1), c(n, 0), d(n);
          for (std::size_t i = 1; i + 1 < n; ++i) {
            const double km = 2 / (1 / k[i - 1] + 1 / k[i]), kp = 2 / (1 / k[i] + 1 / k[i + 1]);
            a[i] = -dt * km / (rc * dz * dz);
            c[i] = -dt * kp / (rc * dz * dz);
            b[i] = 1 - a[i] - c[i];
            d[i] = T[i] + dt * A[i] / rc;
          }
          d[0] = s.t_surface;
          // base: prescribed heat flow (ghost node)
          const double kb = k[n - 1];
          a[n - 1] = -dt * 2 * kb / (rc * dz * dz);
          b[n - 1] = 1 - a[n - 1];
          d[n - 1] = T[n - 1] + dt * (A[n - 1] / rc + 2 * qb(age) / (rc * dz));
          thomas(a, b, c, d);
          T = d;
          if (s.vitrinite)
            for (std::size_t h = 0; h < nh; ++h) {
              if (!std::isfinite(depth[h])) continue;
              const double fz = std::clamp(depth[h] / dz, 0.0, static_cast<double>(n - 1));
              const std::size_t i = std::min(n - 2, static_cast<std::size_t>(fz));
              paths[h].push_back({age, T[i] + (fz - i) * (T[i + 1] - T[i])});
            }
        }
      }
    }
    H.labels.push_back(cur.label);
    H.ages.push_back(cur.age);
    H.temperature.push_back(record(cur.depth, cur.age, sidx));
    if (s.vitrinite) {
      std::vector<double> ro(nh, kNaN);
      for (std::size_t h = 0; h < nh; ++h)
        if (std::isfinite(cur.depth[h]) && paths[h].size() >= 2) ro[h] = easy_ro(paths[h]);
      H.ro.push_back(ro);
    }
  }
  for (std::size_t i = 0; i < n; ++i) { H.present_depth.push_back(i * dz); H.present_temperature.push_back(T[i]); }
  H.surface_flow_present = k[0] * (T[1] - T[0]) / dz;
  // the present column in steady state: the surface flow equals basal flow plus the radiogenic heat
  {
    std::vector<double> Ts;
    const double qs = steady(k, A, qb(0.0), Ts);
    double sumA = 0;
    for (std::size_t i = 0; i < n; ++i) sumA += A[i] * dz * (i == 0 || i + 1 == n ? 0.5 : 1.0);
    // flux back from the steady temperatures at the top cell, corrected by its own heat
    const double q_num = 2 / (1 / k[0] + 1 / k[1]) * (Ts[1] - Ts[0]) / dz + 0.5 * (A[0] + A[1]) * dz * 0.5;
    H.balance_error = std::fabs(q_num - (qb(0.0) + sumA)) / std::max(1e-12, qs);
  }
  return H;
}

PressureHistory pressure_history(const Structure& st, const PressureColumnSpec& s, double x, double y) {
  PressureHistory P;
  auto hist = st.burial_history(x, y);
  if (hist.size() < 2) return P;
  const std::size_t nh = hist.back().depth.size();
  const double gp = (s.rho_g - s.rho_w) * s.g;
  // column from the surface to the deepest deposited horizon of each step, in zeta = z / h
  const std::size_t n = 121;
  std::vector<double> u(n, 0.0);  // excess pressure at zeta_i = i / (n - 1) (0 surface, 1 base)
  auto base_depth = [&](const std::vector<double>& d) {
    double z = 0;
    for (double v : d)
      if (std::isfinite(v)) z = std::max(z, v);
    return z;
  };
  auto perm = [&](int unit, double phi) {
    const int law = unit >= 0 && static_cast<std::size_t>(unit) < s.props.perm_law.size() ? s.props.perm_law[static_cast<std::size_t>(unit)] : 0;
    const double k0 = prop(s.props.k0, unit), b = prop(s.props.b, unit), p0 = prop(s.props.phi0, unit);
    if (!std::isfinite(k0) || !std::isfinite(p0)) return 1e-18;
    if (law == 2) return k0 * std::pow(phi / p0, 3) * std::pow((1 - p0) / std::max(1e-6, 1 - phi), 2);
    return k0 * std::exp((std::isfinite(b) ? b : 10.0) * (phi - p0));
  };
  auto record = [&](const std::vector<double>& depth) {
    std::vector<double> o(nh, kNaN);
    const double h = base_depth(depth);
    for (std::size_t q = 0; q < nh; ++q) {
      if (!std::isfinite(depth[q]) || !(h > 0)) continue;
      const double fz = std::clamp(depth[q] / h * (n - 1), 0.0, static_cast<double>(n - 1));
      const std::size_t i = std::min(n - 2, static_cast<std::size_t>(fz));
      o[q] = u[i] + (fz - i) * (u[i + 1] - u[i]);
    }
    return o;
  };
  for (std::size_t sidx = 0; sidx < hist.size(); ++sidx) {
    const auto& cur = hist[sidx];
    if (sidx > 0) {
      const auto& prev = hist[sidx - 1];
      const double dur = std::isfinite(prev.age) && std::isfinite(cur.age) ? (prev.age - cur.age) * kMa : 0.0;
      if (dur > 0) {
        const int ns = 400;
        const double dt = dur / ns;
        for (int it = 1; it <= ns; ++it) {
          const double w0 = static_cast<double>(it - 1) / ns, w1 = static_cast<double>(it) / ns;
          const auto d0 = lerp_depths(prev.depth, cur.depth, w0), d1 = lerp_depths(prev.depth, cur.depth, w1);
          const double h0 = std::max(1.0, base_depth(d0)), h1 = std::max(1.0, base_depth(d1)), hd = (h1 - h0) / dt;
          // node coefficients: m_v from the Athy law in effective stress (phi = phi0 exp(-c sigma'/(gamma' )))
          // as c phi / ((1 - phi) gamma'), k from the unit's porosity law
          std::vector<double> mv(n), K(n);
          for (std::size_t i = 0; i < n; ++i) {
            const double z = h1 * i / (n - 1);
            const int unit = unit_at_depth(d1, z);
            double phi = st.mean_porosity(unit, z, z);
            if (!std::isfinite(phi)) phi = 0.05;
            const double phi2 = st.mean_porosity(unit, z + 1, z + 1);
            const double c = std::isfinite(phi2) && phi > 0 && phi2 > 0 ? std::max(1e-6, -std::log(phi2 / phi)) : 1e-4;  // 1/m
            mv[i] = std::max(1e-11, c * phi / ((1 - phi) * gp));
            K[i] = perm(unit, phi) / s.viscosity;
          }
          // m_v Du/Dt = 1/h^2 d/dzeta (K du/dzeta) + m_v gamma' hdot in the material frame (each node sees the load
          // of the sediment added above it); zeta from the surface
          const double dzeta = 1.0 / (n - 1);
          std::vector<double> a(n, 0), b(n, 1), c(n, 0), d(n);
          for (std::size_t i = 1; i < n; ++i) {
            const double zeta = static_cast<double>(i) * dzeta;
            // material velocity in zeta: burial carries every point down, d(zeta)/dt = hdot (1 - zeta) / h
            // (zero at the base, which is material); upwind: the material comes from smaller zeta
            const double adv = (1 - zeta) * hd / h1 / dzeta;
            const double km = 2 / (1 / K[i - 1] + 1 / K[i]);
            const double kp = i + 1 < n ? 2 / (1 / K[i] + 1 / K[i + 1]) : km;  // mirror at the base
            const double dif_m = km / (h1 * h1 * dzeta * dzeta), dif_p = kp / (h1 * h1 * dzeta * dzeta);
            if (i + 1 < n) {
              a[i] = dt * (-dif_m / mv[i] - adv);
              c[i] = dt * (-dif_p / mv[i]);
              b[i] = 1 + dt * (dif_m + dif_p) / mv[i] + dt * adv;
            } else {
              a[i] = dt * (-2 * dif_m / mv[i] - adv);
              b[i] = 1 + dt * 2 * dif_m / mv[i] + dt * adv;
            }
            d[i] = u[i] + dt * gp * std::max(0.0, hd);
          }
          a[0] = 0;
          b[0] = 1;
          c[0] = 0;
          d[0] = 0;
          thomas(a, b, c, d);
          u = d;
        }
      }
    }
    P.labels.push_back(cur.label);
    P.ages.push_back(cur.age);
    P.overpressure.push_back(record(cur.depth));
  }
  const auto& now = hist.back().depth;
  const double h = base_depth(now);
  double lith = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double z = h * i / (n - 1);
    const int unit = unit_at_depth(now, z);
    double phi = st.mean_porosity(unit, z, z);
    if (!std::isfinite(phi)) phi = 0.05;
    if (i > 0) lith += (s.rho_g * (1 - phi) + s.rho_w * phi) * s.g * h / (n - 1);
    const double hyd = s.rho_w * s.g * z;
    P.present_depth.push_back(z);
    P.present_hydrostatic.push_back(hyd);
    P.present_lithostatic.push_back(lith);
    P.present_pore.push_back(hyd + u[i]);
  }
  return P;
}

}  // namespace stratonix
