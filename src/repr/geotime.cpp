// Stage 4A (STX-TN-14 section 5): geological-time coordinates. Interface only until Prompt 46.
#include "stratonix/repr/geotime.hpp"

#include <utility>

#include "stratonix/core/error.hpp"

namespace stratonix {

GeologicalTime::GeologicalTime(const Model& m, std::vector<double> ages) : m_(m), ages_(std::move(ages)) { throw Error(ErrorCode::NotImplemented, "stage 4A: geological time (Prompt 46)"); }
GeoTimePoint GeologicalTime::forward(const Vec3&) const { throw Error(ErrorCode::NotImplemented, "stage 4A: geological time (Prompt 46)"); }
std::optional<Vec3> GeologicalTime::inverse(double, double, double, double) const { throw Error(ErrorCode::NotImplemented, "stage 4A: geological time (Prompt 46)"); }
std::vector<Hiatus> GeologicalTime::hiatuses(double, double) const { throw Error(ErrorCode::NotImplemented, "stage 4A: geological time (Prompt 46)"); }

}  // namespace stratonix
