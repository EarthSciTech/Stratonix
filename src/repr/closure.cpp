// Closure analysis (priority flood from the crest) and contour lines (marching squares).
#include <cmath>
#include <cstdint>
#include <map>
#include <queue>

#include "stratonix/repr/repr.hpp"

namespace stratonix::repr {

Closure closure_analysis(const Grid2D& g, bool seal) {
  Closure c;
  const auto& sp = g.spec;
  const std::size_t n = sp.size();
  c.mask = Grid2D(sp, 0.0);
  std::size_t crest = n;
  for (std::size_t k = 0; k < n; ++k) {
    if (std::isnan(g.v[k])) { c.mask.v[k] = kNaN; continue; }
    if (crest == n || g.v[k] > g.v[crest]) crest = k;
  }
  if (crest == n) return c;
  auto ij = [&](std::size_t k) { return std::make_pair(k % sp.nx, k / sp.nx); };
  auto is_exit = [&](std::size_t k) {
    auto [i, j] = ij(k);
    if (i == 0 || j == 0 || i + 1 == sp.nx || j + 1 == sp.ny) return true;
    if (!seal) {
      const std::size_t nb[4] = {k - 1, k + 1, k - sp.nx, k + sp.nx};
      for (auto q : nb)
        if (std::isnan(g.v[q])) return true;
    }
    return false;
  };
  // Priority flood: always expand the highest unvisited neighbour; the running minimum when the
  // first exit node is reached is the spill elevation.
  std::vector<char> seen(n, 0);
  std::priority_queue<std::pair<double, std::size_t>> pq;
  pq.push({g.v[crest], crest});
  seen[crest] = 1;
  double level = g.v[crest];
  std::size_t level_at = crest, last = crest;
  bool exit_found = false;
  while (!pq.empty()) {
    auto [z, k] = pq.top();
    pq.pop();
    if (z < level) { level = z; level_at = k; }
    last = k;
    if (is_exit(k)) { exit_found = true; break; }
    auto [i, j] = ij(k);
    std::size_t nb[4];
    int m = 0;
    if (i > 0) nb[m++] = k - 1;
    if (i + 1 < sp.nx) nb[m++] = k + 1;
    if (j > 0) nb[m++] = k - sp.nx;
    if (j + 1 < sp.ny) nb[m++] = k + sp.nx;
    for (int q = 0; q < m; ++q)
      if (!seen[nb[q]] && !std::isnan(g.v[nb[q]])) { seen[nb[q]] = 1; pq.push({g.v[nb[q]], nb[q]}); }
  }
  (void)last;
  c.found = true;
  c.crest_z = g.v[crest];
  c.crest = {sp.x(crest % sp.nx), sp.y(crest / sp.nx), c.crest_z};
  c.spill_z = level;
  c.spill_point = {sp.x(level_at % sp.nx), sp.y(level_at / sp.nx), level};
  c.spill_at_boundary = exit_found && [&] {
    auto [i, j] = ij(level_at);
    return i == 0 || j == 0 || i + 1 == sp.nx || j + 1 == sp.ny;
  }();
  c.closure = c.crest_z - c.spill_z;
  // Connected region above the spill elevation containing the crest.
  std::vector<std::size_t> stack{crest};
  std::vector<char> in(n, 0);
  in[crest] = 1;
  std::size_t count = 0;
  while (!stack.empty()) {
    std::size_t k = stack.back();
    stack.pop_back();
    ++count;
    c.mask.v[k] = 1;
    auto [i, j] = ij(k);
    std::size_t nb[4];
    int m = 0;
    if (i > 0) nb[m++] = k - 1;
    if (i + 1 < sp.nx) nb[m++] = k + 1;
    if (j > 0) nb[m++] = k - sp.nx;
    if (j + 1 < sp.ny) nb[m++] = k + sp.nx;
    for (int q = 0; q < m; ++q) {
      std::size_t kk = nb[q];
      if (!in[kk] && !std::isnan(g.v[kk]) && g.v[kk] > c.spill_z) { in[kk] = 1; stack.push_back(kk); }
    }
  }
  c.area = static_cast<double>(count) * sp.dx * sp.dy;
  return c;
}

TriMesh grid_to_mesh(const Grid2D& g) {
  const auto& sp = g.spec;
  TriMesh m;
  std::vector<std::uint32_t> id(sp.size(), UINT32_MAX);
  for (std::size_t k = 0; k < sp.size(); ++k)
    if (!std::isnan(g.v[k])) {
      id[k] = static_cast<std::uint32_t>(m.v.size());
      m.v.push_back({sp.x(k % sp.nx), sp.y(k / sp.nx), g.v[k]});
    }
  for (std::size_t j = 0; j + 1 < sp.ny; ++j)
    for (std::size_t i = 0; i + 1 < sp.nx; ++i) {
      std::uint32_t a = id[sp.index(i, j)], b = id[sp.index(i + 1, j)], c = id[sp.index(i + 1, j + 1)], d = id[sp.index(i, j + 1)];
      int n = (a != UINT32_MAX) + (b != UINT32_MAX) + (c != UINT32_MAX) + (d != UINT32_MAX);
      if (n == 4) { m.t.push_back({a, b, c}); m.t.push_back({a, c, d}); }
      else if (n == 3) {
        std::uint32_t v[3];
        int k = 0;
        for (auto x : {a, b, c, d})
          if (x != UINT32_MAX) v[k++] = x;
        m.t.push_back({v[0], v[1], v[2]});
      }
    }
  return m;
}

std::vector<Polyline> contour_lines(const Grid2D& g, double level) {
  const auto& sp = g.spec;
  using Key = std::pair<long long, long long>;  // edge identifier
  std::vector<std::pair<Key, Key>> segs;
  std::map<Key, Vec3> pos;
  auto ekey = [&](std::size_t i1, std::size_t j1, std::size_t i2, std::size_t j2) {
    long long a = static_cast<long long>(i1 + sp.nx * j1), b = static_cast<long long>(i2 + sp.nx * j2);
    return Key{std::min(a, b), std::max(a, b)};
  };
  auto interp = [&](std::size_t i1, std::size_t j1, std::size_t i2, std::size_t j2) {
    double v1 = g.at(i1, j1) - level, v2 = g.at(i2, j2) - level;
    double t = v1 / (v1 - v2);
    Key k = ekey(i1, j1, i2, j2);
    pos[k] = {sp.x(i1) + t * (sp.x(i2) - sp.x(i1)), sp.y(j1) + t * (sp.y(j2) - sp.y(j1)), level};
    return k;
  };
  for (std::size_t j = 0; j + 1 < sp.ny; ++j)
    for (std::size_t i = 0; i + 1 < sp.nx; ++i) {
      double v[4] = {g.at(i, j), g.at(i + 1, j), g.at(i + 1, j + 1), g.at(i, j + 1)};
      if (std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2]) || std::isnan(v[3])) continue;
      int code = 0;
      for (int c = 0; c < 4; ++c)
        if (v[c] > level) code |= 1 << c;
      if (code == 0 || code == 15) continue;
      std::size_t ci[4] = {i, i + 1, i + 1, i}, cj[4] = {j, j, j + 1, j + 1};
      std::vector<Key> ks;
      for (int e = 0; e < 4; ++e) {
        int e2 = (e + 1) % 4;
        if (((code >> e) & 1) != ((code >> e2) & 1)) ks.push_back(interp(ci[e], cj[e], ci[e2], cj[e2]));
      }
      if (ks.size() == 2) segs.push_back({ks[0], ks[1]});
      else if (ks.size() == 4) {  // saddle: resolve with the centre value
        double centre = 0.25 * (v[0] + v[1] + v[2] + v[3]);
        if ((centre > level) == ((code & 1) != 0)) { segs.push_back({ks[0], ks[3]}); segs.push_back({ks[1], ks[2]}); }
        else { segs.push_back({ks[0], ks[1]}); segs.push_back({ks[2], ks[3]}); }
      }
    }
  // Join segments into polylines.
  std::multimap<Key, std::size_t> at;
  for (std::size_t s = 0; s < segs.size(); ++s) { at.insert({segs[s].first, s}); at.insert({segs[s].second, s}); }
  std::vector<char> used(segs.size(), 0);
  std::vector<Polyline> out;
  for (std::size_t s = 0; s < segs.size(); ++s) {
    if (used[s]) continue;
    used[s] = 1;
    std::vector<Key> chain{segs[s].first, segs[s].second};
    for (int dir = 0; dir < 2; ++dir) {
      for (;;) {
        Key end = dir == 0 ? chain.back() : chain.front();
        bool ext = false;
        auto rng = at.equal_range(end);
        for (auto it = rng.first; it != rng.second; ++it) {
          if (used[it->second]) continue;
          used[it->second] = 1;
          Key nxt = segs[it->second].first == end ? segs[it->second].second : segs[it->second].first;
          if (dir == 0) chain.push_back(nxt); else chain.insert(chain.begin(), nxt);
          ext = true;
          break;
        }
        if (!ext) break;
      }
    }
    Polyline pl;
    for (auto& k : chain) pl.push_back(pos[k]);
    out.push_back(std::move(pl));
  }
  return out;
}

}  // namespace stratonix::repr
