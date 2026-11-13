// Reference-versus-optimised parity (Prompt 6). Every optimised kernel is compared with its plain
// reference implementation on the synthetic models and on randomised inputs, including points on
// fault planes, at fault tips and outside the domain; and results must not depend on the number
// of threads.
#include <random>

#include "../stx_test.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/parametric.hpp"
#include "stratonix/package/bench.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

using namespace stratonix;

struct Mode {
  explicit Mode(bool opt) { set_optimised_kernels(opt); }
  ~Mode() { set_optimised_kernels(true); }
};

static std::vector<std::pair<std::string, ModelPtr>> models() {
  register_all_plugins();
  std::vector<std::pair<std::string, ModelPtr>> r;
  auto d = desc::load_description(std::string(stxtest::source_dir()) + "/examples/synthetic_fold_two_faults.yaml");
  BuildOptions o;
  r.push_back({"synthetic", build_only(d, o).model});
  desc::Description p;
  auto rep = desc::validate_description_text(production_model_yaml(200), ".", &p);
  REQUIRE(rep.ok());
  r.push_back({"production", build_only(p, o).model});
  return r;
}

/// Random points in and around the domain plus edge cases: on each fault plane (sigma = 0), at
/// the upper and lower tips, and outside the domain.
static std::vector<Vec3> test_points(const Model& m, std::size_t n, unsigned seed) {
  const Domain& d = m.info().domain;
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> ux(d.x0 - 0.1 * d.sx, d.x0 + 1.1 * d.sx), uy(d.y0 - 0.1 * d.sy, d.y0 + 1.1 * d.sy),
      uz(d.z_base - 200, d.z_top + 200), u01(0, 1);
  std::vector<Vec3> p;
  for (std::size_t i = 0; i < n; ++i) p.push_back({ux(rng), uy(rng), uz(rng)});
  auto* net = static_cast<const FaultNetwork*>(m.capability("faults"));
  for (auto& f : net->faults()) {
    for (int k = 0; k < 200; ++k) {
      double a = (u01(rng) - 0.5) * f.length;
      double bu, bl;
      f.tips(a, bu, bl);
      p.push_back(f.at(a, bu + (bl - bu) * u01(rng)));  // on the plane
      p.push_back(f.at(a, bu));                          // at the upper tip
      p.push_back(f.at(a, bl));                          // at the lower tip
      p.push_back(f.at(0.5 * f.length * (k % 2 ? 1 : -1), bu + (bl - bu) * u01(rng)));  // strike tip
    }
  }
  p.push_back({d.x0 - 5000, d.y0 - 5000, d.z_top + 1000});  // far outside
  return p;
}

TEST_CASE("K1 point evaluation: optimised equals reference bitwise") {
  for (auto& [name, m] : models()) {
    auto p = test_points(*m, 20000, 3);
    EvalResult a, b;
    { Mode r(false); m->evaluate(p, kEvalAll, a); }
    { Mode o(true); m->evaluate(p, kEvalAll, b); }
    CHECK(a.unit == b.unit);
    CHECK(a.block == b.block);
    double worst = 0;
    for (std::size_t i = 0; i < a.scalar.size(); ++i) worst = std::max(worst, std::fabs(a.scalar[i] - b.scalar[i]));
    std::cout << "  " << name << ": " << p.size() << " points, max scalar difference " << worst << " m\n";
    CHECK(worst == 0);
    std::vector<Vec3> ra(p.size()), rb(p.size());
    { Mode r(false); m->restore(p, ra); }
    { Mode o(true); m->restore(p, rb); }
    CHECK(ra == rb);
    for (int f = 0; f < static_cast<int>(m->info().faults.size()); ++f) {
      FaultEval fa, fb;
      { Mode r(false); m->evaluate_fault(f, p, fa); }
      { Mode o(true); m->evaluate_fault(f, p, fb); }
      CHECK(fa.distance == fb.distance);
      CHECK(fa.extent == fb.extent);
    }
  }
}

TEST_CASE("K1 results do not depend on the number of threads") {
  for (auto& [name, m] : models()) {
    auto p = test_points(*m, 50000, 5);
    EvalResult a, b, c;
    set_thread_count(1);
    m->evaluate(p, kEvalAll, a);
    set_thread_count(3);
    m->evaluate(p, kEvalAll, b);
    set_thread_count(8);
    m->evaluate(p, kEvalAll, c);
    set_thread_count(0);
    CHECK(a.scalar == b.scalar);
    CHECK(a.scalar == c.scalar);
    CHECK(a.block == c.block);
  }
}

TEST_CASE("K2 per-block horizon grids: optimised within 1e-6 m of reference") {
  for (auto& [name, m] : models()) {
    GridSpec g = GridSpec::covering(m->info().domain, m->info().domain.sx / 120);
    auto blocks = repr::list_blocks(*m, 16, 16);
    double worst = 0;
    std::size_t defined_mismatch = 0;
    for (std::size_t k = 0; k < std::min<std::size_t>(3, blocks.codes.size()); ++k)
      for (int h = 0; h < static_cast<int>(m->info().horizons.size()); ++h) {
        Grid2D a, b;
        { Mode r(false); a = repr::horizon_grid_block(*m, h, blocks.codes[k], g); }
        { Mode o(true); b = repr::horizon_grid_block(*m, h, blocks.codes[k], g); }
        for (std::size_t i = 0; i < g.size(); ++i) {
          if (std::isnan(a.v[i]) != std::isnan(b.v[i])) { ++defined_mismatch; continue; }
          if (!std::isnan(a.v[i])) worst = std::max(worst, std::fabs(a.v[i] - b.v[i]));
        }
      }
    std::cout << "  " << name << ": max difference " << worst << " m, " << defined_mismatch << " definition mismatches\n";
    CHECK(worst < 1e-6);
    CHECK(defined_mismatch == 0);
  }
}

TEST_CASE("K3 faulted horizon grids: optimised within 1e-6 m where both are defined") {
  for (auto& [name, m] : models()) {
    GridSpec g = GridSpec::covering(m->info().domain, m->info().domain.sx / 160);
    double worst = 0;
    std::size_t mism = 0, total = 0;
    for (int h = 0; h < static_cast<int>(m->info().horizons.size()); ++h) {
      repr::FaultedSurface a, b;
      { Mode r(false); a = repr::horizon_grid_faulted(*m, h, g); }
      { Mode o(true); b = repr::horizon_grid_faulted(*m, h, g); }
      for (std::size_t i = 0; i < g.size(); ++i) {
        ++total;
        if (std::isnan(a.z.v[i]) != std::isnan(b.z.v[i])) { ++mism; continue; }
        if (!std::isnan(a.z.v[i])) worst = std::max(worst, std::fabs(a.z.v[i] - b.z.v[i]));
      }
    }
    double frac = static_cast<double>(mism) / total;
    std::cout << "  " << name << ": max difference " << worst << " m; gap-mask mismatch " << mism << " of " << total << " nodes\n";
    CHECK(worst < 1e-6);
    CHECK(frac < 0.002);
  }
}

TEST_CASE("K4 cut-off lines: optimised within 1e-4 m of reference") {
  for (auto& [name, m] : models()) {
    double worst = 0;
    std::size_t n = 0;
    for (int f = 0; f < std::min(3, static_cast<int>(m->info().faults.size())); ++f) {
      auto fs = repr::fault_surface(*m, f, 50);
      for (int h = 0; h < static_cast<int>(m->info().horizons.size()); ++h) {
        std::vector<repr::CutoffLine> a, b;
        { Mode r(false); a = repr::cutoff_lines(*m, fs, h); }
        { Mode o(true); b = repr::cutoff_lines(*m, fs, h); }
        for (int s = 0; s < 2; ++s)
          for (std::size_t i = 0; i < fs.na; ++i) {
            const Vec3 &p = a[s].by_column[i], &q = b[s].by_column[i];
            if (std::isnan(p.x) || std::isnan(q.x)) { CHECK(std::isnan(p.x) == std::isnan(q.x)); continue; }
            worst = std::max(worst, norm(p - q));
            ++n;
          }
      }
    }
    std::cout << "  " << name << ": " << n << " cut-off points, max difference " << worst << " m\n";
    CHECK(worst < 1e-4);
  }
}

TEST_CASE("K5 voxel grids: optimised equals reference") {
  for (auto& [name, m] : models()) {
    const Domain& d = m->info().domain;
    repr::VoxelGrid a, b;
    { Mode r(false); a = repr::voxel_grid(*m, d.sx / 40, d.sy / 40, (d.z_top - d.z_base) / 40, true); }
    { Mode o(true); b = repr::voxel_grid(*m, d.sx / 40, d.sy / 40, (d.z_top - d.z_base) / 40, true); }
    CHECK(a.unit == b.unit);
    CHECK(a.block == b.block);
  }
}

STX_TEST_MAIN()
