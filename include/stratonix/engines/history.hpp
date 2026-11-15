// Deformation history (stage 2E, STX-TN-09): the event timeline and the model at any event.
#pragma once
#include <vector>

#include "stratonix/desc/description.hpp"
#include "stratonix/engines/builder.hpp"
#include "stratonix/model/model.hpp"

namespace stratonix {

/// Ages (Ma) of the timeline's events, oldest first; empty for a description without events.
std::vector<double> event_ages(const desc::Description& d);

/// The description of the model at time `age` (Ma): horizons deposited later removed, later
/// events removed, events under way scaled to their completed part; reference depths taken from
/// the full model `full` (STX-TN-09 section 2).
desc::Description state_description(const desc::Description& d, const Model& full, double age);

/// The model at time `age`, built by the same builder as the full model.
ModelPtr build_state(const desc::Description& d, const Model& full, double age, BuildLog* log = nullptr);

}  // namespace stratonix
