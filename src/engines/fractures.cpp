// Stage 4B stubs (STX-TN-15): implemented by Prompt 53.
#include "stratonix/engines/fractures.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

double fracture_intensity(double, double, const FractureLaw&) {
  throw Error(ErrorCode::NotImplemented, "fracture_intensity: not implemented yet (Prompt 53)");
}

DfnResult sample_dfn(const BBox&, double, const std::vector<FractureSet>&, std::uint64_t, std::size_t) {
  throw Error(ErrorCode::NotImplemented, "sample_dfn: not implemented yet (Prompt 53)");
}

double realised_p32(const std::vector<Disk>&, const BBox&) {
  throw Error(ErrorCode::NotImplemented, "realised_p32: not implemented yet (Prompt 53)");
}

std::pair<double, double> horizontal_stress(double, double, double, double, double, double, double) {
  throw Error(ErrorCode::NotImplemented, "horizontal_stress: not implemented yet (Prompt 53)");
}

}  // namespace stratonix
