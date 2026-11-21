// Stage 4A (STX-TN-14 section 1): discrete implicit field with free horizon values. Interface only
// until Prompt 43.
#include "stratonix/engines/fields.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

DiscreteField::DiscreteField(const std::vector<HorizonDatum>&, const std::vector<RbfField::Gradient>&, int, double, const BBox&,
                             const DiscreteOptions&, const std::vector<double>&) {
  throw Error(ErrorCode::NotImplemented, "stage 4A: discrete field (Prompt 43)");
}
double DiscreteField::value(const Vec3&) const { throw Error(ErrorCode::NotImplemented, "stage 4A: discrete field (Prompt 43)"); }
Vec3 DiscreteField::gradient(const Vec3&) const { throw Error(ErrorCode::NotImplemented, "stage 4A: discrete field (Prompt 43)"); }
std::string DiscreteField::summary() const { return {}; }

std::string horizon_order_problem(const std::vector<double>&, const std::vector<std::string>&) { throw Error(ErrorCode::NotImplemented, "stage 4A: discrete field (Prompt 43)"); }

}  // namespace stratonix
