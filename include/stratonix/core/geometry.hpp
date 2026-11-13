/// @file geometry.hpp
/// @brief Points, vectors, boxes and polylines. SI units; z is elevation (positive up).
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace stratonix {

constexpr double kPi = 3.14159265358979323846;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
inline double deg2rad(double d) { return d * kPi / 180.0; }
inline double rad2deg(double r) { return r * 180.0 / kPi; }

/// 3D point or vector (metres). x east, y north, z elevation.
struct Vec3 {
  double x = 0, y = 0, z = 0;
  constexpr Vec3() = default;
  constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
  Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
  Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
  Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
  bool operator==(const Vec3& o) const { return x == o.x && y == o.y && z == o.z; }
};
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalized(const Vec3& a) { double n = norm(a); return n > 0 ? a / n : a; }

/// Horizontal unit vector for an azimuth in radians, clockwise from north.
inline Vec3 azimuth_vector(double az_rad) { return {std::sin(az_rad), std::cos(az_rad), 0.0}; }

/// Axis-aligned box.
struct BBox {
  Vec3 lo{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
          std::numeric_limits<double>::max()};
  Vec3 hi{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
          -std::numeric_limits<double>::max()};
  void add(const Vec3& p) {
    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
  }
  bool contains(const Vec3& p, double eps = 0) const {
    return p.x >= lo.x - eps && p.x <= hi.x + eps && p.y >= lo.y - eps && p.y <= hi.y + eps &&
           p.z >= lo.z - eps && p.z <= hi.z + eps;
  }
  bool valid() const { return lo.x <= hi.x; }
};

/// Model domain: a box in x, y and an elevation range [z_min, z_max].
struct Domain {
  double x0 = 0, y0 = 0, sx = 0, sy = 0;  ///< origin and size in x, y
  double z_top = 0, z_base = 0;            ///< elevations (z_top > z_base)
  double spacing = 50;                     ///< default grid spacing
  BBox box() const {
    BBox b; b.add({x0, y0, z_base}); b.add({x0 + sx, y0 + sy, z_top}); return b;
  }
  bool inside_xy(double x, double y, double eps = 1e-9) const {
    return x >= x0 - eps && x <= x0 + sx + eps && y >= y0 - eps && y <= y0 + sy + eps;
  }
};

using Polyline = std::vector<Vec3>;

}  // namespace stratonix
