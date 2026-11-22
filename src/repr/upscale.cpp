// Stage 4C (STX-TN-16 section 1.5): logs upscaled to the cells a well passes through, by a stated
// method, weighted by the measured-depth length of each sample in the cell.
#include "stratonix/repr/upscale.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

#include "stratonix/core/error.hpp"

namespace stratonix::repr {

std::string upscale_method_name(UpscaleMethod m) {
  switch (m) {
    case UpscaleMethod::Arithmetic: return "arithmetic mean";
    case UpscaleMethod::Harmonic: return "harmonic mean";
    case UpscaleMethod::Geometric: return "geometric mean";
    case UpscaleMethod::Rms: return "root mean square";
    case UpscaleMethod::Minimum: return "minimum";
    case UpscaleMethod::Maximum: return "maximum";
    case UpscaleMethod::MostFrequent: return "most frequent";
  }
  return "";
}

namespace {
/// Cell containing p: the column from the grid's areal lattice, the layer by bilinear
/// interpolation of each cell's top and base corner elevations at p. -1 outside.
long cell_of(const CornerPointGrid& g, const Vec3& p) {
  const double fi = (p.x - g.x0) / g.dx, fj = (p.y - g.y0) / g.dy;
  if (fi < 0 || fj < 0 || fi >= g.nx || fj >= g.ny) return -1;
  const int i = static_cast<int>(fi), j = static_cast<int>(fj);
  const double u = fi - i, v = fj - j;
  auto bil = [&](int k, int bottom) {
    const double z00 = g.zcorn[g.zcorn_index(i, j, k, 0 | bottom)], z10 = g.zcorn[g.zcorn_index(i, j, k, 1 | bottom)];
    const double z01 = g.zcorn[g.zcorn_index(i, j, k, 2 | bottom)], z11 = g.zcorn[g.zcorn_index(i, j, k, 3 | bottom)];
    return (1 - u) * (1 - v) * z00 + u * (1 - v) * z10 + (1 - u) * v * z01 + u * v * z11;
  };
  for (int k = 0; k < g.nz; ++k) {
    const double zt = bil(k, 0), zb = bil(k, 4);
    if (p.z <= std::max(zt, zb) && p.z >= std::min(zt, zb)) return static_cast<long>(g.cell(i, j, k));
  }
  return -1;
}
}  // namespace

UpscaledCells upscale_log(const CornerPointGrid& g, const PathLog& log, UpscaleMethod method, double min_valid) {
  if (log.p.size() != log.v.size() || log.md.size() != log.v.size())
    throw Error(ErrorCode::InvalidArgument, "upscale_log: positions, measured depths and values differ in length");
  UpscaledCells out;
  out.method = upscale_method_name(method);
  const std::size_t n = log.v.size();
  // measured-depth length of each sample: half the distance to its neighbours
  std::vector<double> w(n, 0.0);
  for (std::size_t s = 0; s < n; ++s) {
    const double a = s > 0 ? 0.5 * (log.md[s] - log.md[s - 1]) : (n > 1 ? 0.5 * (log.md[1] - log.md[0]) : 1.0);
    const double b = s + 1 < n ? 0.5 * (log.md[s + 1] - log.md[s]) : a;
    w[s] = std::max(0.0, a + b);
  }
  struct Acc { std::vector<std::size_t> samples; };
  std::map<long, Acc> cells;
  std::vector<long> order;
  std::size_t outside = 0;
  for (std::size_t s = 0; s < n; ++s) {
    const long c = cell_of(g, log.p[s]);
    if (c < 0) { ++outside; continue; }
    if (!cells.count(c)) order.push_back(c);
    cells[c].samples.push_back(s);
  }
  if (outside) out.notes.push_back(std::to_string(outside) + " samples outside the grid");
  for (long c : order) {
    const auto& smp = cells[c].samples;
    const std::size_t cc = static_cast<std::size_t>(c);
    if (!g.actnum.empty() && !g.actnum[cc]) { out.notes.push_back("cell " + std::to_string(c) + " is inactive: not filled"); continue; }
    double wt = 0, wv = 0;
    std::vector<std::pair<double, double>> vals;  // value, weight
    for (auto s : smp) {
      wt += w[s];
      if (std::isfinite(log.v[s])) { wv += w[s]; vals.push_back({log.v[s], w[s]}); }
    }
    double mdt = log.md[smp.front()], mdb = log.md[smp.front()];
    for (auto s : smp) { mdt = std::min(mdt, log.md[s]); mdb = std::max(mdb, log.md[s]); }
    if (wt <= 0 || wv / wt < min_valid || vals.empty()) {
      out.missing.push_back(cc);
      out.notes.push_back("cell " + std::to_string(c) + ": " + std::to_string(vals.size()) + " valid of " + std::to_string(smp.size()) +
                          " samples, below the valid fraction " + std::to_string(min_valid) + ": missing");
      continue;
    }
    double r = std::numeric_limits<double>::quiet_NaN(), sw = 0, acc = 0;
    switch (method) {
      case UpscaleMethod::Arithmetic: for (auto& [v, ww] : vals) { acc += v * ww; sw += ww; } r = acc / sw; break;
      case UpscaleMethod::Rms: for (auto& [v, ww] : vals) { acc += v * v * ww; sw += ww; } r = std::sqrt(acc / sw); break;
      case UpscaleMethod::Harmonic:
      case UpscaleMethod::Geometric: {
        std::size_t bad = 0;
        for (auto& [v, ww] : vals) {
          if (!(v > 0)) { ++bad; continue; }
          acc += (method == UpscaleMethod::Harmonic ? 1.0 / v : std::log(v)) * ww;
          sw += ww;
        }
        if (bad) out.notes.push_back("cell " + std::to_string(c) + ": " + std::to_string(bad) + " values not positive, left out of the " + out.method);
        if (sw > 0) r = method == UpscaleMethod::Harmonic ? sw / acc : std::exp(acc / sw);
        break;
      }
      case UpscaleMethod::Minimum: r = vals.front().first; for (auto& [v, ww] : vals) r = std::min(r, v); break;
      case UpscaleMethod::Maximum: r = vals.front().first; for (auto& [v, ww] : vals) r = std::max(r, v); break;
      case UpscaleMethod::MostFrequent: {
        std::map<long long, double> f;
        for (auto& [v, ww] : vals) f[std::llround(v)] += ww;
        double best = -1;
        std::size_t ties = 0;
        for (auto& [k, ww] : f) {
          if (ww > best + 1e-12) { best = ww; r = static_cast<double>(k); ties = 0; }
          else if (std::fabs(ww - best) <= 1e-12) ++ties;
        }
        if (ties) out.notes.push_back("cell " + std::to_string(c) + ": " + std::to_string(ties + 1) + " values equally frequent; the smallest kept");
        break;
      }
    }
    if (!std::isfinite(r)) { out.missing.push_back(cc); continue; }
    out.cells.push_back(cc);
    out.value.push_back(r);
    out.md_top.push_back(mdt);
    out.md_base.push_back(mdb);
  }
  return out;
}

}  // namespace stratonix::repr
