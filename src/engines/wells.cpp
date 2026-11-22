#include "stratonix/engines/wells.hpp"

#include "stratonix/engines/surveys.hpp"

namespace stratonix {

WellPath make_well_path(const desc::WellDesc& w, const desc::Description& d) {
  double x = w.head.si(0), y = w.head.si(1);
  double datum = w.datum_elevation.present ? w.datum_elevation.si() : 0.0;
  if (!w.survey_method.empty() && w.survey.size() >= 2) {
    // stage 4C (STX-TN-16 section 1.3): the stated method and azimuth reference
    SurveySpec sp;
    const std::string& m = w.survey_method;
    sp.method = m == "radius_of_curvature" ? SurveyMethod::RadiusOfCurvature : m == "balanced_tangential" ? SurveyMethod::BalancedTangential
              : m == "average_angle" ? SurveyMethod::AverageAngle : m == "tangential" ? SurveyMethod::Tangential : SurveyMethod::MinimumCurvature;
    const std::string& r = w.azimuth_reference;
    sp.reference = r == "true_north" ? AzimuthReference::TrueNorth : r == "grid_north" ? AzimuthReference::GridNorth
                 : r == "magnetic_north" ? AzimuthReference::MagneticNorth : AzimuthReference::Unstated;
    sp.declination = w.declination.present ? w.declination.si() : 0.0;
    sp.convergence = w.convergence.present ? w.convergence.si() : 0.0;
    std::vector<SurveyStation> st;
    for (auto& s : w.survey) st.push_back({s[0], s[1], s[2]});
    auto pos = survey_positions(st, sp);
    std::vector<double> md;
    std::vector<Vec3> net;
    for (auto& p : pos) { md.push_back(p.md); net.push_back({p.north, p.east, p.tvd}); }
    return WellPath::from_positions(w.name, x, y, datum, md, net, survey_method_name(sp.method));
  }
  if (w.stations.present) {
    std::vector<WellStation> st;
    for (std::size_t r = 0; r < w.stations.rows(); ++r)
      st.push_back({w.stations.values[3 * r] * w.stations.to_si, deg2rad(w.stations.values[3 * r + 1]),
                    deg2rad(w.stations.values[3 * r + 2])});
    return WellPath(w.name, x, y, datum, std::move(st));
  }
  double td = w.total_depth.present ? w.total_depth.si() : d.depth_range.si(1);
  return WellPath(w.name, x, y, datum, td);
}

std::vector<WellPath> make_well_paths(const desc::Description& d) {
  std::vector<WellPath> r;
  for (auto& w : d.wells) r.push_back(make_well_path(w, d));
  return r;
}

}  // namespace stratonix
