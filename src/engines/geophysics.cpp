// Stage 4A (STX-TN-14 section 7, decision D-048): gravity and magnetic forward modelling on right
// rectangular prisms. Vertical gravity of a prism after Nagy et al. (2000); the field of a
// uniformly magnetised prism from the second derivatives of the prism's Newtonian potential
// (Bhattacharyya 1964; Blakely 1995, chapter 9), induced magnetisation only.
//
// Discretisation: the domain is divided into vertical columns of `cell` x `cell`; in each column the
// unit boundaries at the column centre are located by bisection, and each run of one unit is one
// prism. Interfaces are therefore exact at the column centres, and the anomaly changes smoothly
// with the structural parameters (a cell-centre discretisation would change in steps of a cell).
#include "stratonix/engines/geophysics.hpp"

#include <algorithm>
#include <cmath>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"

namespace stratonix {

namespace {
constexpr double kGravity = 6.67430e-11;  // m^3 kg^-1 s^-2 (CODATA 2018)

template <class F>
double corner_sum(const Vec3& o, double x1, double x2, double y1, double y2, double z1, double z2, F&& f) {
  double s = 0;
  const double xs[2] = {x1 - o.x, x2 - o.x}, ys[2] = {y1 - o.y, y2 - o.y}, zs[2] = {z1 - o.z, z2 - o.z};
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k) {
        const double mu = ((i + j + k) % 2 == 0) ? -1.0 : 1.0;  // (-1)^(i+j+k+1)
        s += mu * f(xs[i], ys[j], zs[k]);
      }
  return s;
}
/// ln(a + r) with r = sqrt(a^2 + b^2 + c^2), stable for a < 0: ln(b^2 + c^2) - ln(r - a). Where
/// b = c = 0 the ln(b^2 + c^2) part is dropped: it cancels between the two corners that share b and c.
double log_ar(double a, double b, double c, double r) {
  if (a >= 0) return std::log(a + r);
  const double q = b * b + c * c;
  return q > 0 ? std::log(q) - std::log(r - a) : -std::log(r - a);
}

struct Prism { double x1, x2, y1, y2, z1, z2; int unit; };

/// Columns of prisms: per column, runs of one unit between the domain's base and top.
std::vector<Prism> discretise(const Model& m, double cell) {
  const auto& d = m.info().domain;
  if (!(cell > 0)) throw Error(ErrorCode::InvalidArgument, "geophysics: the cell size must be positive");
  const std::size_t nx = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(d.sx / cell - 1e-9)));
  const std::size_t ny = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(d.sy / cell - 1e-9)));
  const double dx = d.sx / static_cast<double>(nx), dy = d.sy / static_cast<double>(ny);
  const double H = d.z_top - d.z_base;
  const double step = std::min(cell / 2, 25.0);  // Stratonix default sampling along a column
  const std::size_t nz = std::max<std::size_t>(2, static_cast<std::size_t>(std::ceil(H / step)) + 1);
  std::vector<std::vector<Prism>> cols(nx * ny);
  parallel_for(nx * ny, 8, [&](std::size_t b, std::size_t e) {
    std::vector<Vec3> pts(nz);
    EvalResult er;
    for (std::size_t c = b; c < e; ++c) {
      const std::size_t i = c % nx, j = c / nx;
      const double x1 = d.x0 + static_cast<double>(i) * dx, y1 = d.y0 + static_cast<double>(j) * dy, xc = x1 + dx / 2, yc = y1 + dy / 2;
      for (std::size_t k = 0; k < nz; ++k) pts[k] = {xc, yc, d.z_top - H * static_cast<double>(k) / static_cast<double>(nz - 1)};
      m.evaluate(pts, kEvalUnit, er);
      double ztop = d.z_top;
      int u = er.unit[0];
      for (std::size_t k = 1; k < nz; ++k) {
        if (er.unit[k] == u) continue;
        // every boundary between pts[k-1] and pts[k], by bisection to 1 mm (a unit thinner than the
        // sampling step between them is found too)
        double zfrom = pts[k - 1].z;
        while (u != er.unit[k]) {
          double a = zfrom, bz = pts[k].z;  // unit at a is u, at bz is not
          while (a - bz > 1e-3) {
            const double zm = 0.5 * (a + bz);
            (m.unit_at({xc, yc, zm}) == u ? a : bz) = zm;
          }
          const double zbnd = 0.5 * (a + bz);
          cols[c].push_back({x1, x1 + dx, y1, y1 + dy, zbnd, ztop, u});
          ztop = zbnd;
          u = m.unit_at({xc, yc, bz});
          zfrom = bz;
        }
      }
      cols[c].push_back({x1, x1 + dx, y1, y1 + dy, d.z_base, ztop, u});
    }
  });
  std::vector<Prism> out;
  for (auto& c : cols)
    for (auto& p : c)
      if (p.z2 > p.z1) out.push_back(p);
  return out;
}
}  // namespace

double prism_gz(const Vec3& o, double x1, double x2, double y1, double y2, double z1, double z2, double rho) {
  const double s = corner_sum(o, x1, x2, y1, y2, z1, z2, [](double x, double y, double z) {
    const double r = std::sqrt(x * x + y * y + z * z);
    double t = 0;
    if (x != 0) t += x * log_ar(y, x, z, r);
    if (y != 0) t += y * log_ar(x, y, z, r);
    if (z != 0) t -= z * std::atan(x * y / (z * r));
    return t;
  });
  return kGravity * rho * s * 1e5;  // m/s^2 -> mGal; positive down
}

double prism_dt(const Vec3& o, double x1, double x2, double y1, double y2, double z1, double z2, const Vec3& M, const Vec3& F) {
  // K_ij: second derivatives of the prism's potential integral of 1/r; the field tensor
  // H = [[-Kxx, Kxy, Kxz], [Kxy, -Kyy, Kyz], [Kxz, Kyz, -Kzz]] gives B = mu0/(4 pi) H M
  double Kxx = 0, Kyy = 0, Kzz = 0, Kxy = 0, Kxz = 0, Kyz = 0;
  const double xs[2] = {x1 - o.x, x2 - o.x}, ys[2] = {y1 - o.y, y2 - o.y}, zs[2] = {z1 - o.z, z2 - o.z};
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k) {
        const double mu = ((i + j + k) % 2 == 0) ? -1.0 : 1.0;
        const double x = xs[i], y = ys[j], z = zs[k], r = std::sqrt(x * x + y * y + z * z);
        if (x != 0) Kxx += mu * std::atan(y * z / (x * r));
        if (y != 0) Kyy += mu * std::atan(x * z / (y * r));
        if (z != 0) Kzz += mu * std::atan(x * y / (z * r));
        Kxy += mu * log_ar(z, x, y, r);
        Kxz += mu * log_ar(y, x, z, r);
        Kyz += mu * log_ar(x, y, z, r);
      }
  const Vec3 B{-Kxx * M.x + Kxy * M.y + Kxz * M.z, Kxy * M.x - Kyy * M.y + Kyz * M.z, Kxz * M.x + Kyz * M.y - Kzz * M.z};
  return dot(B, F) * 1e-7 * 1e9;  // mu0 / (4 pi) = 1e-7 T m / A; T -> nT
}

Vec3 field_direction(double inc, double dec) {
  return {std::cos(inc) * std::sin(dec), std::cos(inc) * std::cos(dec), -std::sin(inc)};
}

std::vector<double> gravity_anomaly(const Model& m, const PotentialFieldSpec& s, const std::vector<Vec3>& stations) {
  auto prisms = discretise(m, s.cell);
  std::vector<double> contrast(m.info().units.size(), 0.0);
  for (std::size_t u = 0; u < contrast.size(); ++u)
    if (u < s.density.size() && std::isfinite(s.density[u])) contrast[u] = s.density[u] - s.background_density;
  std::vector<double> out(stations.size(), 0.0);
  parallel_for(stations.size(), 4, [&](std::size_t b, std::size_t e) {
    for (std::size_t k = b; k < e; ++k) {
      double g = 0;
      for (auto& p : prisms) {
        const double c = (p.unit >= 0 && static_cast<std::size_t>(p.unit) < contrast.size()) ? contrast[static_cast<std::size_t>(p.unit)] : 0.0;
        if (c != 0) g += prism_gz(stations[k], p.x1, p.x2, p.y1, p.y2, p.z1, p.z2, c);
      }
      out[k] = g;
    }
  });
  return out;
}

std::vector<double> magnetic_anomaly(const Model& m, const PotentialFieldSpec& s, const std::vector<Vec3>& stations) {
  auto prisms = discretise(m, s.cell);
  const Vec3 F = field_direction(s.inclination, s.declination);
  const double H0 = s.intensity * 1e-9 / (4e-7 * kPi);  // inducing field B0 / mu0 (A/m)
  std::vector<double> contrast(m.info().units.size(), 0.0);
  for (std::size_t u = 0; u < contrast.size(); ++u)
    if (u < s.susceptibility.size() && std::isfinite(s.susceptibility[u])) contrast[u] = s.susceptibility[u] - s.background_susceptibility;
  std::vector<double> out(stations.size(), 0.0);
  parallel_for(stations.size(), 4, [&](std::size_t b, std::size_t e) {
    for (std::size_t k = b; k < e; ++k) {
      double t = 0;
      for (auto& p : prisms) {
        const double c = (p.unit >= 0 && static_cast<std::size_t>(p.unit) < contrast.size()) ? contrast[static_cast<std::size_t>(p.unit)] : 0.0;
        if (c != 0) t += prism_dt(stations[k], p.x1, p.x2, p.y1, p.y2, p.z1, p.z2, F * (c * H0), F);
      }
      out[k] = t;
    }
  });
  return out;
}

PotentialFieldSpec potential_field_spec(const desc::Description& d, const ModelInfo& info) {
  const auto& g = d.geophysics;
  PotentialFieldSpec s;
  if (g.cell.present) {
    s.cell = g.cell.si();
  } else if (d.voxel_spacing.present && !d.voxel_spacing.si_list().empty()) {
    s.cell = d.voxel_spacing.si_list()[0];
  } else {
    s.cell = std::max(info.domain.sx, info.domain.sy) / 64.0;  // Stratonix default
  }
  if (g.background_density.present) s.background_density = g.background_density.si();
  if (g.background_susceptibility.present) s.background_susceptibility = g.background_susceptibility.si();
  if (g.inclination.present) s.inclination = g.inclination.si();
  if (g.declination.present) s.declination = g.declination.si();
  if (g.intensity.present) s.intensity = g.intensity.si() * 1e9;  // T -> nT
  s.density.assign(info.units.size(), kNaN);
  s.susceptibility.assign(info.units.size(), kNaN);
  for (std::size_t u = 0; u < info.units.size(); ++u)
    for (auto& gu : g.units)
      if (gu.unit == info.units[u].name) {
        if (gu.density.present) s.density[u] = gu.density.si();
        if (gu.susceptibility.present) s.susceptibility[u] = gu.susceptibility.si();
      }
  return s;
}

PotentialFieldResult potential_field_misfit(const Model& m, const desc::Description& d, const BuildInputs& in) {
  PotentialFieldResult R;
  if (!d.geophysics.present) return R;
  const auto spec = potential_field_spec(d, m.info());
  for (int kind = 0; kind < 2; ++kind) {
    const auto& set = kind == 0 ? d.geophysics.gravity : d.geophysics.magnetics;
    if (set.input.empty()) continue;
    auto it = in.data.find(set.input);
    if (it == in.data.end()) throw Error(ErrorCode::Build, "geophysics: input '" + set.input + "' was not loaded");
    std::vector<Vec3> st;
    for (auto& p : it->second) st.push_back(p.p);
    auto mod = kind == 0 ? gravity_anomaly(m, spec, st) : magnetic_anomaly(m, spec, st);
    const double sig_si = set.uncertainty.present ? set.uncertainty.si() / (kind == 0 ? 1e-5 : 1e-9) : kNaN;  // mGal or nT
    double chi = 0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < st.size(); ++i) {
      StationMisfit sm{kind == 0 ? "gravity" : "magnetics", set.input, st[i], it->second[i].value, mod[i], it->second[i].sigma > 0 ? it->second[i].sigma : sig_si};
      if (std::isfinite(sm.sigma) && sm.sigma > 0) {
        chi += std::pow((sm.modelled - sm.observed) / sm.sigma, 2);
        ++n;
      }
      R.stations.push_back(sm);
    }
    (kind == 0 ? R.n_gravity : R.n_magnetics) = st.size();
    if (n > 0) (kind == 0 ? R.chi2_gravity : R.chi2_magnetics) = chi / static_cast<double>(n);
  }
  return R;
}

}  // namespace stratonix
