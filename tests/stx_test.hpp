// Minimal test harness (Catch2-style macros) so the test suite builds without network access.
// Each test executable registers TEST_CASEs and runs them all; CTest runs the executables.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace stxtest {
struct Case { const char* name; std::function<void()> fn; };
inline std::vector<Case>& cases() { static std::vector<Case> c; return c; }
inline int& failures() { static int f = 0; return f; }
inline int& checks() { static int c = 0; return c; }
struct Reg { Reg(const char* n, std::function<void()> f) { cases().push_back({n, std::move(f)}); } };
struct Abort {};
inline std::string data_dir() {
#ifdef STX_TEST_DATA_DIR
  return STX_TEST_DATA_DIR;
#else
  return ".";
#endif
}
inline std::string source_dir() {
#ifdef STX_SOURCE_DIR
  return STX_SOURCE_DIR;
#else
  return ".";
#endif
}
inline std::string tmp_dir() {
#ifdef STX_TMP_DIR
  return STX_TMP_DIR;
#else
  return "./test_tmp";
#endif
}
}  // namespace stxtest

#define STX_CAT2(a, b) a##b
#define STX_CAT(a, b) STX_CAT2(a, b)
#define TEST_CASE(name)                                                     \
  static void STX_CAT(stx_test_fn_, __LINE__)();                            \
  static stxtest::Reg STX_CAT(stx_test_reg_, __LINE__)(name, STX_CAT(stx_test_fn_, __LINE__)); \
  static void STX_CAT(stx_test_fn_, __LINE__)()

#define STX_FAIL_MSG(expr)                                                                   \
  do {                                                                                       \
    ++stxtest::failures();                                                                   \
    std::cerr << "  FAILED: " << __FILE__ << ":" << __LINE__ << ": " << expr << "\n";        \
  } while (0)

#define CHECK(cond)                                    \
  do {                                                 \
    ++stxtest::checks();                               \
    if (!(cond)) STX_FAIL_MSG(#cond);                  \
  } while (0)
#define REQUIRE(cond)                                  \
  do {                                                 \
    ++stxtest::checks();                               \
    if (!(cond)) { STX_FAIL_MSG(#cond); throw stxtest::Abort{}; } \
  } while (0)
#define CHECK_NEAR(a, b, tol)                                                                    \
  do {                                                                                           \
    ++stxtest::checks();                                                                         \
    double stx_a_ = (a), stx_b_ = (b);                                                           \
    if (!(std::fabs(stx_a_ - stx_b_) <= (tol))) {                                                \
      std::ostringstream stx_s_;                                                                 \
      stx_s_ << #a << " = " << stx_a_ << " vs " << #b << " = " << stx_b_ << " (tol " << (tol) << ")"; \
      STX_FAIL_MSG(stx_s_.str());                                                                \
    }                                                                                            \
  } while (0)
#define CHECK_THROWS(stmt)                                       \
  do {                                                           \
    ++stxtest::checks();                                         \
    bool stx_thrown_ = false;                                    \
    try { stmt; } catch (...) { stx_thrown_ = true; }            \
    if (!stx_thrown_) STX_FAIL_MSG("expected exception: " #stmt); \
  } while (0)

namespace stxtest {
/// A case runs when the filter is empty or one of its '|'-separated parts occurs in its name.
inline bool selected(const std::string& name, const std::string& filter) {
  if (filter.empty()) return true;
  std::size_t a = 0;
  while (a <= filter.size()) {
    const std::size_t b = std::min(filter.find('|', a), filter.size());
    const std::string part = filter.substr(a, b - a);
    if (!part.empty() && name.find(part) != std::string::npos) return true;
    a = b + 1;
  }
  return false;
}
}  // namespace stxtest
#define STX_TEST_MAIN()                                                                 \
  int main(int argc, char** argv) {                                                     \
    std::string filter = argc > 1 ? argv[1] : "";                                       \
    int run = 0;                                                                        \
    for (auto& c : stxtest::cases()) {                                                  \
      if (!stxtest::selected(c.name, filter)) continue;                                  \
      ++run;                                                                            \
      int before = stxtest::failures();                                                 \
      try { c.fn(); } catch (stxtest::Abort&) {                                         \
      } catch (std::exception& e) { STX_FAIL_MSG(std::string("exception: ") + e.what()); } \
      std::cout << (stxtest::failures() == before ? "PASS " : "FAIL ") << c.name << "\n"; \
    }                                                                                   \
    std::cout << run << " test case(s), " << stxtest::checks() << " check(s), "        \
              << stxtest::failures() << " failure(s)\n";                                \
    return stxtest::failures() ? 1 : 0;                                                 \
  }
