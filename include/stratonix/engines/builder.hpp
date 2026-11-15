/// @file builder.hpp
/// @brief Structure-builder plugins turn a validated description into a Model.
#pragma once
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "stratonix/core/grid.hpp"
#include "stratonix/desc/description.hpp"
#include "stratonix/model/model.hpp"
#include "stratonix/registry/registry.hpp"

namespace stratonix {

/// Messages and numbers a builder reports (corrections applied, resolved tips, ...).
struct BuildLog {
  std::vector<std::string> notes;
  std::vector<std::pair<std::string, std::string>> values;  ///< key/value for the manifest
  void note(std::string s) { notes.push_back(std::move(s)); }
  void value(std::string k, std::string v) { values.emplace_back(std::move(k), std::move(v)); }
};

/// A data point for data-conditioned builders: position, optional unit normal (orientation
/// data, pointing up), and a label (horizon or fault name; empty if the data entry names it).
struct DataPoint {
  Vec3 p;
  Vec3 n{0, 0, 0};
  bool oriented = false;
  std::string label;
  double sigma = -1;  ///< stated uncertainty (m; degrees for orientations); < 0 = none (stage 3A)
};

/// Data read from description inputs, keyed by input id (filled by the pipeline's importers).
struct BuildInputs {
  std::map<std::string, Grid2D> grids;
  std::map<std::string, std::vector<DataPoint>> data;    ///< points and orientations
  std::map<std::string, std::vector<Polyline>> lines;    ///< fault sticks
};

/// Structure builder plugin interface.
class StructureBuilder : public registry::PluginBase {
 public:
  /// True if this builder can build the description.
  virtual bool accepts(const desc::Description& d) const = 0;
  virtual ModelPtr build(const desc::Description& d, const BuildInputs& in, BuildLog& log) const = 0;
};

/// Build a model with the first registered builder that accepts the description, preferring
/// `preferred` if given. Throws Error(Build) if none accepts it.
ModelPtr build_model(const desc::Description& d, const BuildInputs& in, BuildLog& log,
                     const std::string& preferred = {});

/// Common ModelInfo filled from a description (horizons, units, domain, CRS, basis).
ModelInfo model_info_from(const desc::Description& d);

/// Register this layer's plugins (called by stratonix::register_all_plugins()).
void register_engine_plugins();

}  // namespace stratonix
