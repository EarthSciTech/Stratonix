// Stage 4A (STX-TN-14 section 4): Bayesian conditioning. Interface only until Prompt 45.
#include "stratonix/engines/posterior.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

double PriorSpec::log_density(double) const { throw Error(ErrorCode::NotImplemented, "stage 4A: posterior sampling (Prompt 45)"); }
PosteriorResult sample_ensemble(const LogDensity&, const std::vector<std::vector<double>>&, const PosteriorOptions&) { throw Error(ErrorCode::NotImplemented, "stage 4A: posterior sampling (Prompt 45)"); }
double split_rhat(const std::vector<std::vector<double>>&) { throw Error(ErrorCode::NotImplemented, "stage 4A: posterior sampling (Prompt 45)"); }
double effective_sample_size(const std::vector<std::vector<double>>&) { throw Error(ErrorCode::NotImplemented, "stage 4A: posterior sampling (Prompt 45)"); }
PosteriorResult sample_posterior(const desc::Description&, const BuildInputs&, const PosteriorOptions&, BuildLog&) { throw Error(ErrorCode::NotImplemented, "stage 4A: posterior sampling (Prompt 45)"); }

}  // namespace stratonix
