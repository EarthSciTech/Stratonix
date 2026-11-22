// Stage 4B stubs (STX-TN-15): implemented by Prompt 53.
#include "stratonix/engines/thermal.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

std::vector<double> steady_temperature(const std::vector<ThermalLayer>&, double, double) {
  throw Error(ErrorCode::NotImplemented, "steady_temperature: not implemented yet (Prompt 53)");
}

std::vector<double> transient_temperature(const TransientSpec&) {
  throw Error(ErrorCode::NotImplemented, "transient_temperature: not implemented yet (Prompt 53)");
}

double gibson_overpressure_base(const GibsonSpec&) {
  throw Error(ErrorCode::NotImplemented, "gibson_overpressure_base: not implemented yet (Prompt 53)");
}

}  // namespace stratonix
