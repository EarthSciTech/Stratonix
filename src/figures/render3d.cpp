// Software 3D renderer for figure views: orthographic camera, z-buffer, Lambert shading,
// 2x supersampling. Deterministic (no threads, fixed arithmetic order).
#include <algorithm>
#include <cmath>
#include <limits>

#include "stratonix/figures/figures.hpp"

namespace stratonix::fig {

Image render_meshes(const std::vector<Mesh3>& meshes, double az_deg, double el_deg, double ve, int W, int H, const BBox& box) {
  const int S = 2;  // supersampling
  const int w = W * S, h = H * S;
  const double az = deg2rad(az_deg), el = deg2rad(el_deg);
  Vec3 ctr = (box.lo + box.hi) * 0.5;
  // Camera looks from (azimuth, elevation) towards the centre.
  Vec3 fwd{-std::sin(az) * std::cos(el), -std::cos(az) * std::cos(el), -std::sin(el)};
  Vec3 right = normalized(cross(fwd, {0, 0, 1}));
  Vec3 up = cross(right, fwd);
  Vec3 light = normalized(fwd * -1.0 + right * -0.5 + up * 0.8);
  auto view = [&](const Vec3& p) {
    Vec3 q{p.x - ctr.x, p.y - ctr.y, (p.z - ctr.z) * ve};
    return Vec3{dot(q, right), dot(q, up), dot(q, fwd)};
  };
  // Fit the meshes (clipped to the box) to the image.
  double umin = 1e300, umax = -1e300, vmin = 1e300, vmax = -1e300;
  for (auto& m : meshes)
    for (auto& p : m.mesh.v) {
      if (!std::isfinite(p.z)) continue;
      Vec3 v = view({std::clamp(p.x, box.lo.x, box.hi.x), std::clamp(p.y, box.lo.y, box.hi.y), std::clamp(p.z, box.lo.z, box.hi.z)});
      umin = std::min(umin, v.x); umax = std::max(umax, v.x); vmin = std::min(vmin, v.y); vmax = std::max(vmax, v.y);
    }
  if (!(umax > umin) || !(vmax > vmin)) { umin = vmin = -1; umax = vmax = 1; }
  double sc = 0.94 * std::min(w / (umax - umin), h / (vmax - vmin));
  double ox = 0.5 * (w - (umax - umin) * sc), oy = 0.5 * (h - (vmax - vmin) * sc);
  std::vector<float> zb(static_cast<std::size_t>(w) * h, std::numeric_limits<float>::infinity());
  std::vector<std::uint8_t> col(static_cast<std::size_t>(w) * h * 4, 0);
  // Rasterise one mesh into (zbuf, cbuf), nearest fragment wins.
  auto draw = [&](const Mesh3& m, std::vector<float>& zbuf, std::vector<std::uint8_t>& cbuf, const std::vector<float>* occl) {
    std::vector<Vec3> sv(m.mesh.v.size());
    for (std::size_t i = 0; i < sv.size(); ++i) {
      Vec3 v = view(m.mesh.v[i]);
      sv[i] = {ox + (v.x - umin) * sc, oy + (vmax - v.y) * sc, v.z};
    }
    for (auto& t : m.mesh.t) {
      const Vec3 &a = sv[t[0]], &b = sv[t[1]], &c = sv[t[2]];
      Vec3 wa{m.mesh.v[t[0]].x, m.mesh.v[t[0]].y, m.mesh.v[t[0]].z * ve};
      Vec3 wb{m.mesh.v[t[1]].x, m.mesh.v[t[1]].y, m.mesh.v[t[1]].z * ve};
      Vec3 wc{m.mesh.v[t[2]].x, m.mesh.v[t[2]].y, m.mesh.v[t[2]].z * ve};
      Vec3 n = normalized(cross(wb - wa, wc - wa));
      double lam = 0.35 + 0.65 * std::fabs(dot(n, light));
      std::uint8_t r = static_cast<std::uint8_t>(std::min(255.0, m.color.r * lam));
      std::uint8_t g = static_cast<std::uint8_t>(std::min(255.0, m.color.g * lam));
      std::uint8_t bl = static_cast<std::uint8_t>(std::min(255.0, m.color.b * lam));
      double area = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
      if (!(std::fabs(area) > 1e-12)) continue;
      int x0 = std::max(0, static_cast<int>(std::floor(std::min({a.x, b.x, c.x}))));
      int x1 = std::min(w - 1, static_cast<int>(std::ceil(std::max({a.x, b.x, c.x}))));
      int y0 = std::max(0, static_cast<int>(std::floor(std::min({a.y, b.y, c.y}))));
      int y1 = std::min(h - 1, static_cast<int>(std::ceil(std::max({a.y, b.y, c.y}))));
      for (int y = y0; y <= y1; ++y)
        for (int x = x0; x <= x1; ++x) {
          double px = x + 0.5, py = y + 0.5;
          double w0 = ((b.x - px) * (c.y - py) - (b.y - py) * (c.x - px)) / area;
          double w1 = ((c.x - px) * (a.y - py) - (c.y - py) * (a.x - px)) / area;
          double w2 = 1 - w0 - w1;
          if (w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9) continue;
          float depth = static_cast<float>(w0 * a.z + w1 * b.z + w2 * c.z);
          std::size_t k = static_cast<std::size_t>(y) * w + x;
          if (occl && depth >= (*occl)[k]) continue;
          if (depth < zbuf[k]) {
            zbuf[k] = depth;
            cbuf[4 * k] = r; cbuf[4 * k + 1] = g; cbuf[4 * k + 2] = bl; cbuf[4 * k + 3] = 255;
          }
        }
    }
  };
  for (auto& m : meshes)
    if (m.alpha >= 1) draw(m, zb, col, nullptr);
  // See-through meshes: each in its own layer (hidden by opaque meshes), blended in list order.
  for (auto& m : meshes) {
    if (m.alpha >= 1) continue;
    std::vector<float> lz(zb.size(), std::numeric_limits<float>::infinity());
    std::vector<std::uint8_t> lc(col.size(), 0);
    draw(m, lz, lc, &zb);
    for (std::size_t k = 0; k < lz.size(); ++k) {
      if (!lc[4 * k + 3]) continue;
      double al = m.alpha;
      if (!col[4 * k + 3]) {  // nothing behind: draw against white
        for (int c = 0; c < 3; ++c) col[4 * k + c] = static_cast<std::uint8_t>(std::lround(al * lc[4 * k + c] + (1 - al) * 255));
        col[4 * k + 3] = 255;
      } else {
        for (int c = 0; c < 3; ++c) col[4 * k + c] = static_cast<std::uint8_t>(std::lround(al * lc[4 * k + c] + (1 - al) * col[4 * k + c]));
      }
    }
  }
  // Edge darkening where depth jumps (outlines), then downsample.
  Image im;
  im.w = W;
  im.h = H;
  im.rgba.assign(static_cast<std::size_t>(W) * H * 4, 0);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      unsigned acc[4] = {0, 0, 0, 0};
      for (int dy = 0; dy < S; ++dy)
        for (int dx = 0; dx < S; ++dx) {
          std::size_t k = static_cast<std::size_t>(y * S + dy) * w + (x * S + dx);
          for (int c = 0; c < 4; ++c) acc[c] += col[4 * k + c];
        }
      std::size_t o = (static_cast<std::size_t>(y) * W + x) * 4;
      unsigned a = acc[3] / (S * S);
      im.rgba[o + 3] = static_cast<std::uint8_t>(a);
      if (acc[3] > 0)
        for (int c = 0; c < 3; ++c) im.rgba[o + c] = static_cast<std::uint8_t>(acc[c] * 255 / acc[3]);
    }
  return im;
}

}  // namespace stratonix::fig
