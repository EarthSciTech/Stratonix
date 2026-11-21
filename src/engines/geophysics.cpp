// Stage 4A (STX-TN-14 section 7): gravity and magnetics. Interface only until Prompt 47.
#include "stratonix/engines/geophysics.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

double prism_gz(const Vec3&, double, double, double, double, double, double, double) { throw Error(ErrorCode::NotImplemented, "stage 4A: gravity and magnetics (Prompt 47)"); }
double prism_dt(const Vec3&, double, double, double, double, double, double, const Vec3&, const Vec3&) { throw Error(ErrorCode::NotImplemented, "stage 4A: gravity and magnetics (Prompt 47)"); }
Vec3 field_direction(double, double) { throw Error(ErrorCode::NotImplemented, "stage 4A: gravity and magnetics (Prompt 47)"); }
std::vector<double> gravity_anomaly(const Model&, const PotentialFieldSpec&, const std::vector<Vec3>&) { throw Error(ErrorCode::NotImplemented, "stage 4A: gravity and magnetics (Prompt 47)"); }
std::vector<double> magnetic_anomaly(const Model&, const PotentialFieldSpec&, const std::vector<Vec3>&) { throw Error(ErrorCode::NotImplemented, "stage 4A: gravity and magnetics (Prompt 47)"); }

}  // namespace stratonix
