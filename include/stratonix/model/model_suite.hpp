/// @file model_suite.hpp
/// @brief Interface conformance suite. Every engine must pass it.
#pragma once
#include <string>
#include <vector>

#include "stratonix/model/model.hpp"

namespace stratonix {

struct SuiteCase {
  std::string name;
  bool passed = false;
  std::string detail;
};
struct SuiteReport {
  std::vector<SuiteCase> cases;
  bool passed() const;
  std::string text() const;
};
struct SuiteOptions {
  std::size_t samples = 4000;   ///< random points in the domain
  unsigned seed = 12345;
  double tol = 1e-6;            ///< numerical tolerance (m)
};

/// Run the interface suite against a model:
///  S1 info is consistent (names, units tile horizons, domain valid, faults <= 64)
///  S2 evaluate fills exactly the requested fields with the right sizes
///  S3 unit index agrees with the horizon scalars (unit = number of horizons the point is below)
///  S4 scalars are ordered: s_0 <= s_1 <= ... (a point is never further above a shallower horizon than a deeper one: no crossing)
///  S5 scalar_in_block equals evaluate's scalar for points in that block
///  S6 block code bit f agrees with the sign of evaluate_fault distance
///  S7 results are identical in single and multi-threaded evaluation and in any batch split
///  S8 restore leaves points of the fixed block (if any) unmoved and keeps array sizes
///  S9 scalar decreases downward along verticals inside a block (monotone in z)
SuiteReport run_model_suite(const Model& m, const SuiteOptions& opt = {});

}  // namespace stratonix
