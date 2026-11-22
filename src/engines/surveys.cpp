// Stage 4C (STX-TN-16 sections 1.3 and 1.4): deviation surveys by a stated method, and ages in wells.
#include "stratonix/engines/ages.hpp"
#include "stratonix/engines/surveys.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "stratonix/core/error.hpp"
#include "stratonix/core/geometry.hpp"

namespace stratonix {

namespace {
const double kPiS = 3.14159265358979323846;
double wrap_pi(double a) {
  while (a > kPiS) a -= 2 * kPiS;
  while (a < -kPiS) a += 2 * kPiS;
  return a;
}
}  // namespace

std::string survey_method_name(SurveyMethod m) {
  switch (m) {
    case SurveyMethod::MinimumCurvature: return "minimum curvature";
    case SurveyMethod::RadiusOfCurvature: return "radius of curvature";
    case SurveyMethod::BalancedTangential: return "balanced tangential";
    case SurveyMethod::AverageAngle: return "average angle";
    case SurveyMethod::Tangential: return "tangential";
  }
  return "";
}

std::vector<SurveyPoint> survey_positions(const std::vector<SurveyStation>& st, const SurveySpec& spec) {
  if (spec.reference == AzimuthReference::Unstated)
    throw Error(ErrorCode::Validation, "survey: the azimuth reference (true, grid or magnetic north) is not stated; it is never assumed");
  // azimuths to grid north: true = magnetic + declination; grid = true - convergence
  auto grid_azi = [&](double a) {
    if (spec.reference == AzimuthReference::MagneticNorth) a += spec.declination;
    if (spec.reference != AzimuthReference::GridNorth) a -= spec.convergence;
    return a;
  };
  std::vector<SurveyPoint> out;
  out.reserve(st.size());
  if (st.empty()) return out;
  SurveyPoint p;
  p.md = st[0].md;
  out.push_back(p);
  for (std::size_t i = 1; i < st.size(); ++i) {
    const double dmd = st[i].md - st[i - 1].md;
    if (!(dmd > 0)) throw Error(ErrorCode::Validation, "survey: measured depth must increase (station " + std::to_string(i + 1) + ")");
    const double i1 = st[i - 1].inc, i2 = st[i].inc, a1 = grid_azi(st[i - 1].azi), a2 = grid_azi(st[i].azi);
    // dogleg from the unit tangents (stable for small angles)
    const double tx = std::sin(i2) * std::cos(a2) - std::sin(i1) * std::cos(a1), ty = std::sin(i2) * std::sin(a2) - std::sin(i1) * std::sin(a1),
                 tz = std::cos(i2) - std::cos(i1);
    const double beta = 2 * std::asin(std::min(1.0, 0.5 * std::sqrt(tx * tx + ty * ty + tz * tz)));
    double dn = 0, de = 0, dv = 0;
    switch (spec.method) {
      case SurveyMethod::MinimumCurvature:
      case SurveyMethod::BalancedTangential: {
        const double rf = (spec.method == SurveyMethod::BalancedTangential || beta < 1e-9) ? 1.0 : 2.0 / beta * std::tan(beta / 2);
        dn = 0.5 * dmd * (std::sin(i1) * std::cos(a1) + std::sin(i2) * std::cos(a2)) * rf;
        de = 0.5 * dmd * (std::sin(i1) * std::sin(a1) + std::sin(i2) * std::sin(a2)) * rf;
        dv = 0.5 * dmd * (std::cos(i1) + std::cos(i2)) * rf;
        break;
      }
      case SurveyMethod::AverageAngle: {
        const double im = 0.5 * (i1 + i2), am = a1 + 0.5 * wrap_pi(a2 - a1);
        dn = dmd * std::sin(im) * std::cos(am);
        de = dmd * std::sin(im) * std::sin(am);
        dv = dmd * std::cos(im);
        break;
      }
      case SurveyMethod::Tangential:
        dn = dmd * std::sin(i2) * std::cos(a2);
        de = dmd * std::sin(i2) * std::sin(a2);
        dv = dmd * std::cos(i2);
        break;
      case SurveyMethod::RadiusOfCurvature: {
        const double di = i2 - i1, da = wrap_pi(a2 - a1);
        double h;
        if (std::fabs(di) < 1e-9) { dv = dmd * std::cos(i1); h = dmd * std::sin(i1); }
        else { dv = dmd * (std::sin(i2) - std::sin(i1)) / di; h = dmd * (std::cos(i1) - std::cos(i2)) / di; }
        if (std::fabs(da) < 1e-9) { dn = h * std::cos(a1); de = h * std::sin(a1); }
        else { dn = h * (std::sin(a1 + da) - std::sin(a1)) / da; de = h * (std::cos(a1) - std::cos(a1 + da)) / da; }
        break;
      }
    }
    p.md = st[i].md;
    p.north += dn;
    p.east += de;
    p.tvd += dv;
    p.dls = beta / dmd;
    out.push_back(p);
  }
  return out;
}

// ---- ages ------------------------------------------------------------------------------------
namespace {
std::vector<std::pair<double, double>> depth_age(const std::vector<AgeDatum>& data) {
  std::vector<std::pair<double, double>> v;
  for (auto& a : data) v.push_back({a.depth_base >= a.depth ? 0.5 * (a.depth + a.depth_base) : a.depth, a.age});
  std::sort(v.begin(), v.end());
  return v;
}
}  // namespace

double age_at_depth(const std::vector<AgeDatum>& data, double depth, AgeRule) {
  // constant rate per unit equals linear in depth between dated samples when no unit boundary lies
  // between them; the unit form is applied by the well correlation, which knows the units
  const auto v = depth_age(data);
  if (v.size() < 2 || depth < v.front().first || depth > v.back().first) return kNaN;
  for (std::size_t i = 1; i < v.size(); ++i)
    if (depth <= v[i].first) {
      const double d0 = v[i - 1].first, d1 = v[i].first;
      if (d1 == d0) return v[i].second;
      return v[i - 1].second + (depth - d0) / (d1 - d0) * (v[i].second - v[i - 1].second);
    }
  return kNaN;
}

double depth_of_age(const std::vector<AgeDatum>& data, double age, AgeRule) {
  const auto v = depth_age(data);
  for (std::size_t i = 1; i < v.size(); ++i) {
    const double a0 = v[i - 1].second, a1 = v[i].second;
    if ((age - a0) * (age - a1) <= 0 && a1 != a0) return v[i - 1].first + (age - a0) / (a1 - a0) * (v[i].first - v[i - 1].first);
    if (age == a0) return v[i - 1].first;
  }
  if (!v.empty() && age == v.back().second) return v.back().first;
  return kNaN;
}

Correlation correlate_by_age(const std::vector<WellAges>& wells, const std::vector<std::pair<std::string, double>>& horizon_ages, AgeRule rule) {
  Correlation c;
  for (auto& [h, a] : horizon_ages) c.horizons.push_back(h);
  for (auto& w : wells) c.wells.push_back(w.well);
  c.depth.assign(horizon_ages.size(), std::vector<double>(wells.size(), kNaN));
  for (std::size_t j = 0; j < wells.size(); ++j) {
    std::vector<AgeDatum> d = wells[j].data;
    std::sort(d.begin(), d.end(), [](const AgeDatum& a, const AgeDatum& b) { return a.depth < b.depth; });
    for (std::size_t i = 1; i < d.size(); ++i)
      if (d[i].age < d[i - 1].age - (d[i].sigma + d[i - 1].sigma))
        c.conflicts.push_back({wells[j].well, "age inversion: " + d[i].name + " (" + std::to_string(d[i].age) + " Ma) below " + d[i - 1].name + " (" +
                                                  std::to_string(d[i - 1].age) + " Ma)", d[i].depth});
    for (std::size_t h = 0; h < horizon_ages.size(); ++h) {
      c.depth[h][j] = depth_of_age(d, horizon_ages[h].second, rule);
      if (!std::isfinite(c.depth[h][j]))
        c.conflicts.push_back({wells[j].well, "horizon " + horizon_ages[h].first + " (" + std::to_string(horizon_ages[h].second) +
                                                  " Ma) lies outside the well's dated interval", kNaN});
    }
  }
  return c;
}

}  // namespace stratonix
