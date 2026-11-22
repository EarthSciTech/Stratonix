// Stage 4B stubs (STX-TN-15): implemented by Prompt 49.
#include "stratonix/engines/bodies.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

double SheetBody::thickness_at(double, double) const {
  throw Error(ErrorCode::NotImplemented, "SheetBody::thickness_at: not implemented yet (Prompt 49)");
}

double sheet_volume(const SheetBody&, double) {
  throw Error(ErrorCode::NotImplemented, "sheet_volume: not implemented yet (Prompt 49)");
}

double SillBody::thickness_at(double) const {
  throw Error(ErrorCode::NotImplemented, "SillBody::thickness_at: not implemented yet (Prompt 49)");
}

double sill_volume(const SillBody&) {
  throw Error(ErrorCode::NotImplemented, "sill_volume: not implemented yet (Prompt 49)");
}

double forced_fold_uplift(const SillBody&, double) {
  throw Error(ErrorCode::NotImplemented, "forced_fold_uplift: not implemented yet (Prompt 49)");
}

double forced_fold_volume(const SillBody&, double) {
  throw Error(ErrorCode::NotImplemented, "forced_fold_volume: not implemented yet (Prompt 49)");
}

InjectiteComplex build_injectites(const InjectiteSpec&, std::uint64_t) {
  throw Error(ErrorCode::NotImplemented, "build_injectites: not implemented yet (Prompt 49)");
}

}  // namespace stratonix
