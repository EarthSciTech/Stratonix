// Stage 4D (STX-TN-17 section 1, decision 24): geostatistics. Interfaces first (Prompt 58); the
// methods are built in Prompt 59.
#include "stratonix/engines/geostat.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix::geo {
namespace {
[[noreturn]] void todo(const char* what) {
  throw Error(ErrorCode::NotImplemented, std::string(what) + ": not implemented (stage 4D, Prompt 59)");
}
}  // namespace

double Variogram::sill() const { todo("Variogram::sill"); }
double Variogram::gamma(const Vec3&) const { todo("Variogram::gamma"); }
double Variogram::covariance(const Vec3&) const { todo("Variogram::covariance"); }
ExperimentalVariogram experimental_variogram(const std::vector<Vec3>&, const std::vector<double>&, const VariogramDirection&) {
  todo("experimental_variogram");
}
Variogram fit_variogram(const ExperimentalVariogram&, const Variogram&, bool, std::size_t) { todo("fit_variogram"); }
NormalScore::NormalScore(std::vector<double>, std::vector<double>, double, double) { todo("NormalScore"); }
double NormalScore::forward(double) const { todo("NormalScore::forward"); }
double NormalScore::back(double) const { todo("NormalScore::back"); }
std::vector<double> decluster_weights(const std::vector<Vec3>&, double) { todo("decluster_weights"); }
std::vector<double> sgs(const std::vector<Vec3>&, const std::vector<ConditioningPoint>&, const Variogram&, const SearchOptions&, std::uint64_t,
                        const std::vector<double>&) {
  todo("sgs");
}
std::vector<int> sis(const std::vector<Vec3>&, const std::vector<std::pair<std::size_t, int>>&, const std::vector<double>&, const std::vector<Variogram>&,
                     const SearchOptions&, std::uint64_t) {
  todo("sis");
}
ObjectResult object_facies(const std::vector<Vec3>&, const Vec3&, const Vec3&, const ObjectSpec&, const std::vector<std::pair<std::size_t, int>>&,
                           std::uint64_t, std::size_t) {
  todo("object_facies");
}

}  // namespace stratonix::geo
