#include "stratonix/model/model.hpp"

#include <cmath>
#include <random>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/model/model_suite.hpp"

namespace stratonix {

int Model::unit_at(const Vec3& p) const {
  EvalResult r;
  evaluate(std::span<const Vec3>(&p, 1), kEvalUnit, r);
  return r.unit[0];
}
BlockCode Model::block_at(const Vec3& p) const {
  EvalResult r;
  evaluate(std::span<const Vec3>(&p, 1), kEvalBlock, r);
  return r.block[0];
}

// ---- flat layers ---------------------------------------------------------------------------
FlatLayersModel::FlatLayersModel(ModelInfo info, std::vector<double> horizon_z)
    : info_(std::move(info)), z_(std::move(horizon_z)) {
  if (z_.size() != info_.horizons.size()) throw Error(ErrorCode::InvalidArgument, "one elevation per horizon is required");
  for (std::size_t i = 1; i < z_.size(); ++i)
    if (!(z_[i] < z_[i - 1])) throw Error(ErrorCode::InvalidArgument, "horizon elevations must decrease downwards");
  if (!info_.faults.empty()) throw Error(ErrorCode::InvalidArgument, "the flat-layers model has no faults");
  info_.engine = "flat-layers";
}
void FlatLayersModel::evaluate(std::span<const Vec3> pts, unsigned fields, EvalResult& out) const {
  std::size_t n = pts.size(), nh = z_.size();
  out.n = n;
  out.unit.assign(fields & kEvalUnit ? n : 0, 0);
  out.block.assign(fields & kEvalBlock ? n : 0, 0);
  out.scalar.assign(fields & kEvalScalars ? n * nh : 0, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    int u = 0;
    for (std::size_t h = 0; h < nh; ++h) {
      double s = pts[i].z - z_[h];
      if (fields & kEvalScalars) out.scalar[h * n + i] = s;
      if (s < 0) u = static_cast<int>(h) + 1;
    }
    if (fields & kEvalUnit) out.unit[i] = u;
  }
}
void FlatLayersModel::scalar_in_block(std::span<const Vec3> pts, BlockCode, int horizon, std::span<double> out) const {
  for (std::size_t i = 0; i < pts.size(); ++i) out[i] = pts[i].z - z_.at(static_cast<std::size_t>(horizon));
}
void FlatLayersModel::evaluate_fault(int, std::span<const Vec3>, FaultEval&) const {
  throw Error(ErrorCode::InvalidArgument, "the flat-layers model has no faults");
}
void FlatLayersModel::restore(std::span<const Vec3> pts, std::span<Vec3> out) const {
  for (std::size_t i = 0; i < pts.size(); ++i) out[i] = pts[i];
}

// ---- interface suite -----------------------------------------------------------------------
bool SuiteReport::passed() const {
  for (auto& c : cases)
    if (!c.passed) return false;
  return !cases.empty();
}
std::string SuiteReport::text() const {
  std::ostringstream s;
  for (auto& c : cases) s << (c.passed ? "PASS " : "FAIL ") << c.name << (c.detail.empty() ? "" : " - " + c.detail) << "\n";
  return s.str();
}

SuiteReport run_model_suite(const Model& m, const SuiteOptions& opt) {
  SuiteReport rep;
  auto add = [&](std::string name, bool ok, std::string detail = {}) { rep.cases.push_back({std::move(name), ok, std::move(detail)}); };
  const auto& inf = m.info();
  const std::size_t nh = inf.horizons.size(), nf = inf.faults.size();
  // S1
  {
    bool ok = nh > 0 && inf.units.size() == nh + 1 && inf.domain.sx > 0 && inf.domain.sy > 0 &&
              inf.domain.z_top > inf.domain.z_base && nf <= 64;
    for (std::size_t k = 1; ok && k < inf.units.size(); ++k)
      ok = inf.units[k].top_horizon == static_cast<int>(k) - 1 &&
           (k == nh ? inf.units[k].base_horizon == -1 : inf.units[k].base_horizon == static_cast<int>(k));
    add("S1 info consistent", ok, ok ? "" : "units must tile horizons: unit k lies between horizons k-1 and k");
  }
  // sample points
  std::mt19937_64 rng(opt.seed);
  const Domain& d = inf.domain;
  std::uniform_real_distribution<double> ux(d.x0, d.x0 + d.sx), uy(d.y0, d.y0 + d.sy), uz(d.z_base, d.z_top);
  std::vector<Vec3> pts(opt.samples);
  for (auto& p : pts) p = {ux(rng), uy(rng), uz(rng)};
  EvalResult all, one;
  m.evaluate(pts, kEvalAll, all);
  // S2
  {
    EvalResult r;
    m.evaluate(pts, kEvalUnit, r);
    bool ok = all.n == pts.size() && all.unit.size() == pts.size() && all.block.size() == pts.size() &&
              all.scalar.size() == pts.size() * nh && r.unit.size() == pts.size() && r.scalar.empty() && r.block.empty();
    add("S2 evaluate fills requested fields", ok);
  }
  // S3, S4
  {
    std::size_t bad3 = 0, bad4 = 0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
      if (unit_from_scalars(&all.scalar[i], all.n, nh) != all.unit[i]) ++bad3;
      for (std::size_t h = 1; h < nh; ++h)
        if (all.s(h, i) < all.s(h - 1, i) - opt.tol) { ++bad4; break; }
    }
    add("S3 unit agrees with scalars", bad3 == 0, std::to_string(bad3) + " mismatches");
    add("S4 no horizon crossing at sampled points", bad4 == 0, std::to_string(bad4) + " points");
  }
  // S5
  {
    std::size_t bad = 0;
    double worst = 0;
    std::vector<double> s(1);
    for (std::size_t i = 0; i < pts.size(); i += 7)
      for (std::size_t h = 0; h < nh; ++h) {
        m.scalar_in_block(std::span<const Vec3>(&pts[i], 1), all.block[i], static_cast<int>(h), s);
        double e = std::fabs(s[0] - all.s(h, i));
        worst = std::max(worst, e);
        if (e > opt.tol) ++bad;
      }
    add("S5 scalar_in_block matches evaluate in own block", bad == 0, "max diff " + fmt_shortest(worst) + " m");
  }
  // S6
  {
    std::size_t bad = 0;
    for (std::size_t f = 0; f < nf; ++f) {
      FaultEval fe;
      m.evaluate_fault(static_cast<int>(f), pts, fe);
      if (fe.distance.size() != pts.size() || fe.extent.size() != pts.size() || fe.displacement.size() != pts.size()) { ++bad; continue; }
      for (std::size_t i = 0; i < pts.size(); ++i) {
        if (std::fabs(fe.distance[i]) < 1e-6) continue;
        bool hw = fe.distance[i] > 0;
        bool bit = (all.block[i] >> f) & 1u;
        if (hw != bit) ++bad;
      }
    }
    add("S6 block bits agree with fault distance sign", bad == 0, std::to_string(bad) + " mismatches");
  }
  // S7
  {
    unsigned saved = thread_count();
    set_thread_count(1);
    EvalResult a;
    m.evaluate(pts, kEvalAll, a);
    set_thread_count(4);
    EvalResult b;
    m.evaluate(pts, kEvalAll, b);
    set_thread_count(saved);
    bool ok = a.unit == b.unit && a.block == b.block && a.scalar == b.scalar && a.scalar == all.scalar;
    // batch split: evaluate the first 100 points separately
    EvalResult c;
    std::size_t k = std::min<std::size_t>(100, pts.size());
    m.evaluate(std::span<const Vec3>(pts.data(), k), kEvalAll, c);
    for (std::size_t i = 0; ok && i < k; ++i)
      for (std::size_t h = 0; h < nh; ++h) ok = ok && c.s(h, i) == all.s(h, i);
    add("S7 deterministic across threads and batch splits", ok);
  }
  // S8
  {
    std::vector<Vec3> r(pts.size());
    m.restore(pts, r);
    bool ok = true;
    for (auto& p : r) ok = ok && std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
    add("S8 restore returns finite positions", ok);
  }
  // S9 monotone along verticals (sample a few columns inside one block)
  {
    std::size_t bad = 0;
    for (std::size_t i = 0; i < 50 && i < pts.size(); ++i) {
      std::vector<Vec3> col;
      for (int k = 0; k <= 40; ++k) col.push_back({pts[i].x, pts[i].y, d.z_top - (d.z_top - d.z_base) * k / 40.0});
      for (std::size_t h = 0; h < nh; ++h) {
        std::vector<double> s(col.size());
        m.scalar_in_block(col, all.block[i], static_cast<int>(h), s);
        for (std::size_t k = 1; k < s.size(); ++k)
          if (s[k] > s[k - 1] + opt.tol) { ++bad; break; }
      }
    }
    add("S9 scalar monotone with depth inside a block", bad == 0, std::to_string(bad) + " columns");
  }
  return rep;
}

}  // namespace stratonix
