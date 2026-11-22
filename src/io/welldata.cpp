// Stage 4C stubs (STX-TN-16 section 1): implemented by Prompt 55.
#include "stratonix/io/welldata.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix::io {

const LogCurve* WellLogSet::find(const std::string& m) const {
  for (auto& c : curves)
    if (c.mnemonic == m) return &c;
  return nullptr;
}
WellLogSet read_las(const std::filesystem::path&) { throw Error(ErrorCode::NotImplemented, "read_las: not implemented yet (Prompt 55)"); }
WellLogSet parse_las(const std::string&, const std::string&) { throw Error(ErrorCode::NotImplemented, "parse_las: not implemented yet (Prompt 55)"); }
std::vector<WellLogSet> read_dlis(const std::filesystem::path&) { throw Error(ErrorCode::NotImplemented, "read_dlis: not implemented yet (Prompt 55)"); }

}  // namespace stratonix::io
