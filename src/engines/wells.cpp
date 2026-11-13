#include "stratonix/engines/wells.hpp"

namespace stratonix {

WellPath make_well_path(const desc::WellDesc& w, const desc::Description& d) {
  double x = w.head.si(0), y = w.head.si(1);
  double datum = w.datum_elevation.present ? w.datum_elevation.si() : 0.0;
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
