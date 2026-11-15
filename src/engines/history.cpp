// Deformation history (stage 2E, STX-TN-09). Implemented in Prompt 29.
#include "stratonix/engines/history.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

std::vector<double> event_ages(const desc::Description&) { return {}; }

desc::Description state_description(const desc::Description& d, const Model&, double) { return d; }

ModelPtr build_state(const desc::Description&, const Model&, double, BuildLog*) {
  throw Error(ErrorCode::Build, "model states at an event come with stage 2E (Prompt 29)");
}

}  // namespace stratonix
