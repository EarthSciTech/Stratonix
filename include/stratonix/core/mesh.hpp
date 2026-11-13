/// @file mesh.hpp
/// @brief Triangulated surface meshes with quality checks.
#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"

namespace stratonix {

/// Triangulated mesh; vertex array and triangle index array (0-based).
struct TriMesh {
  std::vector<Vec3> v;
  std::vector<std::array<std::uint32_t, 3>> t;
  bool empty() const { return t.empty(); }
  double area() const;
};

/// Mesh quality summary.
struct MeshQuality {
  std::size_t triangles = 0;
  std::size_t degenerate = 0;        ///< area below eps or repeated vertex
  std::size_t boundary_edges = 0;    ///< edges used by one triangle
  std::size_t nonmanifold_edges = 0; ///< edges used by more than two triangles
  double min_area = 0, max_area = 0;
  bool watertight() const { return boundary_edges == 0 && nonmanifold_edges == 0; }
};

MeshQuality mesh_quality(const TriMesh& m, double area_eps = 1e-6);

}  // namespace stratonix
