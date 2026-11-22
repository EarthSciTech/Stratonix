// Stage 4B stubs (STX-TN-15): implemented by Prompt 52.
#include "stratonix/engines/restoration3d.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

UnfoldResult flexural_slip_unfold(const std::vector<Vec3>&, double, double) {
  throw Error(ErrorCode::NotImplemented, "flexural_slip_unfold: not implemented yet (Prompt 52)");
}

VolumetricResult volumetric_restore(const Model&, const VolumetricSpec&) {
  throw Error(ErrorCode::NotImplemented, "volumetric_restore: not implemented yet (Prompt 52)");
}

double block_porosity(double, double, double) {
  throw Error(ErrorCode::NotImplemented, "block_porosity: not implemented yet (Prompt 52)");
}

double airy_tectonic_subsidence(double, double, double, double, double, double) {
  throw Error(ErrorCode::NotImplemented, "airy_tectonic_subsidence: not implemented yet (Prompt 52)");
}

double flexural_response(double, double, double, double, double) {
  throw Error(ErrorCode::NotImplemented, "flexural_response: not implemented yet (Prompt 52)");
}

double eroded_section(const std::vector<std::pair<double, double>>&, double, double) {
  throw Error(ErrorCode::NotImplemented, "eroded_section: not implemented yet (Prompt 52)");
}

}  // namespace stratonix
