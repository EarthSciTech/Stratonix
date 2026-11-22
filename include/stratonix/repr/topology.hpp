/// @file topology.hpp
/// @brief Stage 4A (STX-TN-14 section 6): the block-and-contact graph of a model and the
/// watertightness of its surfaces.
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"
#include "stratonix/model/model.hpp"
#include "stratonix/repr/repr.hpp"

namespace stratonix {

/// Node: a unit in a fault block.
struct ContactNode { int unit = -1; int block = -1; std::size_t cells = 0; };
/// Edge: two nodes in contact, the surface between them and its sampled area (m^2).
struct ContactEdge {
  int a = -1, b = -1;        ///< node indices, a < b
  std::string surface;       ///< "horizon:<name>", "fault:<name>" or "unconformity:<name>"; non-adjacent units in one
                             ///< block are "unconformity" when the series change or a unit between them terminates
  double area = 0;
};
struct ContactGraph {
  std::vector<ContactNode> nodes;
  std::vector<ContactEdge> edges;
  double spacing = 0;
  std::vector<std::string> unit_names;  ///< the model's unit names (for the JSON)
  std::string json() const;
};

/// Sample the model on a lattice of the given spacing and record face-neighbour contacts.
ContactGraph contact_graph(const Model& m, double spacing);
/// Same nodes and edges (areas ignored); `diff` lists edges present in only one graph.
bool same_topology(const ContactGraph& a, const ContactGraph& b, std::vector<std::string>* diff = nullptr);

struct WatertightGap { std::string surface; Vec3 at; double gap = 0; };
struct WatertightReport {
  double max_gap = 0;               ///< m
  double tolerance = 0;             ///< m
  std::vector<WatertightGap> gaps;  ///< those above tolerance (largest first, at most 20)
  bool ok() const { return max_gap <= tolerance; }
};
/// Horizon edges: on a surface grid of the given spacing, every horizon node next to a gap or another
/// block is followed along the horizon (bisection) to where it leaves its block; the gap is the
/// distance from there to the nearest fault surface, other horizon or domain boundary.
WatertightReport watertightness(const Model& m, double spacing, double tolerance);
/// The same on a given surface grid, with the faulted horizon surfaces supplied (checks reuse theirs).
WatertightReport watertightness(const Model& m, const GridSpec& g, const std::function<const repr::FaultedSurface&(int)>& surface, double tolerance);

}  // namespace stratonix
