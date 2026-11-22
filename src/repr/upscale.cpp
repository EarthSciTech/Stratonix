// Stage 4C stubs (STX-TN-16 section 1.5): implemented by Prompt 55.
#include "stratonix/repr/upscale.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix::repr {

std::string upscale_method_name(UpscaleMethod m) {
  switch (m) {
    case UpscaleMethod::Arithmetic: return "arithmetic mean";
    case UpscaleMethod::Harmonic: return "harmonic mean";
    case UpscaleMethod::Geometric: return "geometric mean";
    case UpscaleMethod::Rms: return "root mean square";
    case UpscaleMethod::Minimum: return "minimum";
    case UpscaleMethod::Maximum: return "maximum";
    case UpscaleMethod::MostFrequent: return "most frequent";
  }
  return "";
}
UpscaledCells upscale_log(const CornerPointGrid&, const PathLog&, UpscaleMethod, double) {
  throw Error(ErrorCode::NotImplemented, "upscale_log: not implemented yet (Prompt 55)");
}

}  // namespace stratonix::repr
