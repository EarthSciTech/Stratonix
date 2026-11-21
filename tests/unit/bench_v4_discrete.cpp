// Stage 4A benchmark (STX-TN-14 section 1): the discrete field with free horizon values on
// 1 000 000 interface points and 50 000 orientations over 20 km x 20 km, three horizons, no
// faults, on a 161 x 161 x 41 grid, within 600 s and 4 GB on the CI Linux release machine.
// Debug and sanitiser builds skip it (exit code 77).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <memory>
#include <random>

#include "../stx_test.hpp"
#include "stratonix/engines/fields.hpp"
#if defined(__linux__)
#include <sys/resource.h>
#endif

using namespace stratonix;

namespace {
double surf(int h, double x, double y) {
  double r2 = (x - 10000) * (x - 10000) + (y - 10000) * (y - 10000);
  return -(1000.0 + 150.0 * h + 0.004 * h * (x - 10000)) + 200.0 * std::exp(-r2 / (2 * 5000.0 * 5000.0));
}
double peak_rss_gb() {
#if defined(__linux__)
  rusage u{};
  getrusage(RUSAGE_SELF, &u);
  return static_cast<double>(u.ru_maxrss) / (1024.0 * 1024.0);
#else
  return 0;
#endif
}
}  // namespace

int main() {
#if !defined(NDEBUG) || defined(STX_SANITIZE)
  std::printf("SKIP bench_v4_discrete: optimised (release) builds only\n");
  return 77;
#else
  std::mt19937_64 rng(1000000);
  std::uniform_real_distribution<double> U(50, 19950);
  std::vector<HorizonDatum> pts;
  std::vector<RbfField::Gradient> grads;
  pts.reserve(1000000);
  for (int i = 0; i < 1000000; ++i) {
    int h = i % 3;
    double x = U(rng), y = U(rng), z = surf(h, x, y);
    pts.push_back({{x, y, z}, h});
    if (i % 20 == 0) {
      const double e = 1;
      double zx = (surf(h, x + e, y) - surf(h, x - e, y)) / (2 * e), zy = (surf(h, x, y + e) - surf(h, x, y - e)) / (2 * e);
      grads.push_back({{x, y, z}, normalized(Vec3{-zx, -zy, 1})});
    }
  }
  BBox box;
  box.add({0, 0, -1600});
  box.add({20000, 20000, -700});
  DiscreteOptions o;
  o.spacing = 125;
  o.vertical_spacing = 22.5;
  auto t0 = std::chrono::steady_clock::now();
  std::unique_ptr<DiscreteField> fp;
  try {
    fp = std::make_unique<DiscreteField>(pts, grads, 0, surf(0, 10000, 10000), box, o);
  } catch (const std::exception& e) {
    std::printf("FAILED: %s\n", e.what());
    return 1;
  }
  const DiscreteField& f = *fp;
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  auto& v = f.horizon_values();
  double err = 0;
  for (int h = 1; h < 3; ++h) err = std::max(err, std::fabs(f.value({10000, 10000, surf(h, 10000, 10000)}) - v[static_cast<std::size_t>(h)]));
  std::printf("bench_v4_discrete: %zu x %zu x %zu nodes, %zu unknowns, %d iterations, residual %.2e, %.1f s, peak %.2f GB, "
              "value error at the centre %.2f m\n",
              f.nx(), f.ny(), f.nz(), f.unknowns(), f.iterations(), f.relative_residual(), secs, peak_rss_gb(), err);
  int fails = 0;
  if (secs > 600) { std::printf("FAILED: %.1f s > 600 s\n", secs); ++fails; }
  if (peak_rss_gb() > 4.0) { std::printf("FAILED: %.2f GB > 4 GB\n", peak_rss_gb()); ++fails; }
  if (err > 5.0) { std::printf("FAILED: value error %.2f m > 5 m\n", err); ++fails; }
  return fails ? 1 : 0;
#endif
}
