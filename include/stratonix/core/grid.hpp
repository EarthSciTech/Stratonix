/// @file grid.hpp
/// @brief Regular 2D grids (unrotated in v1) with contiguous storage, x fastest.
#pragma once
#include <cmath>
#include <cstddef>
#include <vector>

#include "stratonix/core/geometry.hpp"

namespace stratonix {

/// Geometry of a regular grid: node (i, j) is at (x0 + i*dx, y0 + j*dy).
struct GridSpec {
  double x0 = 0, y0 = 0, dx = 1, dy = 1;
  std::size_t nx = 0, ny = 0;
  std::size_t size() const { return nx * ny; }
  double x(std::size_t i) const { return x0 + static_cast<double>(i) * dx; }
  double y(std::size_t j) const { return y0 + static_cast<double>(j) * dy; }
  std::size_t index(std::size_t i, std::size_t j) const { return i + nx * j; }
  /// Grid covering a domain at a spacing (node-centred, includes both edges).
  static GridSpec covering(const Domain& d, double spacing) {
    GridSpec g; g.x0 = d.x0; g.y0 = d.y0; g.dx = spacing; g.dy = spacing;
    g.nx = static_cast<std::size_t>(std::llround(d.sx / spacing)) + 1;
    g.ny = static_cast<std::size_t>(std::llround(d.sy / spacing)) + 1;
    return g;
  }
  bool operator==(const GridSpec& o) const {
    return x0 == o.x0 && y0 == o.y0 && dx == o.dx && dy == o.dy && nx == o.nx && ny == o.ny;
  }
};

/// Scalar field on a regular grid. Undefined nodes hold NaN.
struct Grid2D {
  GridSpec spec;
  std::vector<double> v;
  Grid2D() = default;
  explicit Grid2D(const GridSpec& s, double fill = kNaN) : spec(s), v(s.size(), fill) {}
  double& at(std::size_t i, std::size_t j) { return v[spec.index(i, j)]; }
  double at(std::size_t i, std::size_t j) const { return v[spec.index(i, j)]; }
  /// Bilinear interpolation; NaN if outside or if any corner is undefined.
  double sample(double x, double y) const {
    double fx = (x - spec.x0) / spec.dx, fy = (y - spec.y0) / spec.dy;
    if (fx < -1e-9 || fy < -1e-9 || fx > spec.nx - 1 + 1e-9 || fy > spec.ny - 1 + 1e-9) return kNaN;
    fx = std::clamp(fx, 0.0, double(spec.nx - 1)); fy = std::clamp(fy, 0.0, double(spec.ny - 1));
    std::size_t i = std::min<std::size_t>(std::size_t(fx), spec.nx > 1 ? spec.nx - 2 : 0);
    std::size_t j = std::min<std::size_t>(std::size_t(fy), spec.ny > 1 ? spec.ny - 2 : 0);
    if (spec.nx == 1 || spec.ny == 1) return at(std::min(i, spec.nx - 1), std::min(j, spec.ny - 1));
    double tx = fx - i, ty = fy - j;
    double a = at(i, j), b = at(i + 1, j), c = at(i, j + 1), d = at(i + 1, j + 1);
    return (a * (1 - tx) + b * tx) * (1 - ty) + (c * (1 - tx) + d * tx) * ty;
  }
};

}  // namespace stratonix
