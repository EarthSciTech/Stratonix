// Uncertainty (stage 3C, STX-TN-13). Interfaces opened by Prompt 40; implemented by Prompt 41.
#include "stratonix/package/uncertainty.hpp"

#include <stdexcept>

namespace stratonix::unc {
namespace {
[[noreturn]] void todo(const char* what) {
  throw std::runtime_error(std::string("not implemented yet: ") + what + " (STX-TN-13, Prompt 41)");
}
}  // namespace

std::uint64_t member_seed(std::uint64_t, std::uint64_t) { todo("member seeds"); }
Rng::Rng(std::uint64_t) : s_{0, 0, 0, 0} { todo("generator"); }
std::uint64_t Rng::next() { todo("generator"); }
double Rng::uniform() { todo("generator"); }
double normal_cdf(double) { todo("normal CDF"); }
double normal_quantile(double) { todo("normal quantile"); }
double Distribution::quantile(double) const { todo("inverse CDF"); }
double Distribution::mean() const { todo("distribution mean"); }
double Distribution::variance() const { todo("distribution variance"); }
EnsembleSpec ensemble_spec(const desc::Description&) { todo("ensemble specification"); }
std::vector<MemberSample> sample_members(const EnsembleSpec&) { todo("sampling"); }
std::vector<double> correlated_uniforms(const EnsembleSpec&, Rng&) { todo("correlated draws"); }
std::vector<double> cholesky(const std::vector<double>&, int) { todo("Cholesky factor"); }
desc::Node member_document(const desc::Node&, const EnsembleSpec&, const MemberSample&) { todo("member descriptions"); }
void set_parameter(desc::Node&, const std::string&, double, const std::string&, const std::string&, const std::string&) { todo("setting a parameter"); }
bool get_parameter(const desc::Node&, const std::string&, double*, std::string*) { todo("reading a parameter"); }
Stats statistics(std::vector<double>) { todo("statistics"); }
double first_order_index(const std::vector<double>&, const std::vector<double>&) { todo("first-order indices"); }
std::size_t EnsembleResult::ok() const { todo("ensemble results"); }
EnsembleResult run_ensemble(const std::filesystem::path&, const std::filesystem::path&, const EnsembleOptions&) { todo("the ensemble runner"); }
bool parse_key_value(const std::string&, double*, std::string*) { todo("key-number values"); }

}  // namespace stratonix::unc
