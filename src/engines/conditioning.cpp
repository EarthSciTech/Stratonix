// Stage 3A (STX-TN-11): faults from data and hybrid conditioning. Specified; built by Prompts 35
// and 36.
#include "stratonix/engines/conditioning.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

namespace {
[[noreturn]] void not_built(const char* what, const char* prompt) {
  throw Error(ErrorCode::Build, std::string(what) + ": specified in STX-TN-11, not built yet (" + prompt + ")");
}
}  // namespace

std::vector<ThrowEstimate> estimate_throws(const FittedPlane&, const std::vector<HorizonPoint>&, std::size_t, const ThrowOptions&) {
  not_built("throw estimation from horizon cut-offs", "Prompt 35");
}
StickProfile profile_from_sticks(const std::vector<Polyline>&) { not_built("fault profile from sticks", "Prompt 35"); }
ConditioningResult condition_to_data(const desc::Description&, const BuildInputs&, BuildLog&) {
  not_built("hybrid conditioning", "Prompt 36");
}

}  // namespace stratonix
