// Stage 4B stubs (STX-TN-15): implemented by Prompt 51.
#include "stratonix/engines/depositional.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

std::vector<Vec3> onlap_line(const std::function<double(double, double)>&, const std::function<double(double, double)>&, const GridSpec&) {
  throw Error(ErrorCode::NotImplemented, "onlap_line: not implemented yet (Prompt 51)");
}

std::vector<Clinoform> clinoform_set(const ClinoformSpec&) {
  throw Error(ErrorCode::NotImplemented, "clinoform_set: not implemented yet (Prompt 51)");
}

bool Geobody::contains(const Vec3&) const {
  throw Error(ErrorCode::NotImplemented, "Geobody::contains: not implemented yet (Prompt 51)");
}

std::vector<Geobody> place_channels(const ChannelSpec&, const UnitBox&, std::uint64_t) {
  throw Error(ErrorCode::NotImplemented, "place_channels: not implemented yet (Prompt 51)");
}

double net_to_gross(const std::vector<Geobody>&, const UnitBox&, double) {
  throw Error(ErrorCode::NotImplemented, "net_to_gross: not implemented yet (Prompt 51)");
}

}  // namespace stratonix
