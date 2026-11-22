// Stage 4C stubs (STX-TN-16 section 2): implemented by Prompt 56.
#include "stratonix/io/seismic.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix::io {

double SegyVolume::sample(double, double, double) const { throw Error(ErrorCode::NotImplemented, "SegyVolume::sample: not implemented yet (Prompt 56)"); }
SegyVolume read_segy(const std::filesystem::path&, const SegyLayout&) { throw Error(ErrorCode::NotImplemented, "read_segy: not implemented yet (Prompt 56)"); }
std::vector<PickSet> read_picks(const std::filesystem::path&, const std::string&, const std::string&, const std::string&) {
  throw Error(ErrorCode::NotImplemented, "read_picks: not implemented yet (Prompt 56)");
}
StickSet read_fault_sticks(const std::filesystem::path&, const std::string&, const std::string&, const std::string&) {
  throw Error(ErrorCode::NotImplemented, "read_fault_sticks: not implemented yet (Prompt 56)");
}

}  // namespace stratonix::io
