// Builder selection, ModelInfo from descriptions, and the flat-layers builder plugin.
#include "stratonix/engines/builder.hpp"

#include <cmath>

#include "stratonix/core/error.hpp"

namespace stratonix {

ModelInfo model_info_from(const desc::Description& d) {
  ModelInfo mi;
  mi.name = d.model.name;
  mi.version = d.model.version;
  mi.title = d.model.title;
  mi.basis = d.is_design() ? EvidenceBasis::Design : EvidenceBasis::Interpreted;
  mi.crs = d.crs;
  mi.domain = d.domain();
  for (std::size_t i = 0; i < d.horizons.size(); ++i)
    mi.horizons.push_back({d.horizons[i].name, static_cast<int>(i), d.horizons[i].intended_offset,
                           d.series_of(static_cast<int>(i)) >= 0 ? d.series[static_cast<std::size_t>(d.series_of(static_cast<int>(i)))].name : ""});
  const std::size_t nh = d.horizons.size();
  mi.units.resize(nh + 1);
  mi.units[0] = {d.above_top_name.empty() ? "above_" + d.horizons[0].name : d.above_top_name, -1, 0, "#E8EEF4"};
  mi.units[nh] = {d.below_base_name.empty() ? "below_" + d.horizons[nh - 1].name : d.below_base_name,
                  static_cast<int>(nh) - 1, -1, "#B9B4AA"};
  for (auto& u : d.units) {  // colours of the outer units when the units list names them
    if (u.base == "model_base" && u.name == mi.units[nh].name && !u.color.empty()) mi.units[nh].color = u.color;
    if (u.top == "model_top" && u.name == mi.units[0].name && !u.color.empty()) mi.units[0].color = u.color;
  }
  for (std::size_t k = 1; k < nh; ++k) {
    mi.units[k].top_horizon = static_cast<int>(k) - 1;
    mi.units[k].base_horizon = static_cast<int>(k);
    for (auto& u : d.units)
      if (d.horizon_index(u.top) == static_cast<int>(k) - 1) { mi.units[k].name = u.name; mi.units[k].color = u.color; }
  }
  return mi;
}

ModelPtr build_model(const desc::Description& d, const BuildInputs& in, BuildLog& log, const std::string& preferred) {
  auto& reg = registry::Registry::global();
  std::vector<const registry::Entry*> cands;
  if (!preferred.empty()) {
    if (auto* e = reg.find(registry::Kind::StructureBuilder, preferred)) cands.push_back(e);
    else throw Error(ErrorCode::InvalidArgument, "no structure builder named '" + preferred + "'");
  }
  for (auto* e : reg.list(registry::Kind::StructureBuilder)) cands.push_back(e);
  for (auto* e : cands) {
    auto b = reg.create<StructureBuilder>(registry::Kind::StructureBuilder, e->name);
    if (b->accepts(d)) {
      log.value("structure_builder", e->name + "@" + e->version);
      return b->build(d, in, log);
    }
  }
  throw Error(ErrorCode::Build, "no registered structure builder accepts this description");
}

namespace {
/// Flat layers: descriptions with no structure primitives, no imported surfaces and no faults.
class FlatLayersBuilder final : public StructureBuilder {
 public:
  bool accepts(const desc::Description& d) const override {
    if (d.implicit.present || !d.primitives.empty() || !d.faults.empty() || d.regional_type != "flat") return false;
    if (!d.series.empty() || d.salt.present() || d.compaction.present()) return false;  // stages 2C, 2D: parametric builder
    for (auto& h : d.horizons)
      if (!h.surface_input.empty() || !h.iso_input.empty() || !h.iso_trend.empty()) return false;
    return true;
  }
  ModelPtr build(const desc::Description& d, const BuildInputs&, BuildLog& log) const override {
    ModelInfo mi = model_info_from(d);
    std::vector<double> z;
    const desc::WellDesc* rw = d.ref_well.empty() ? nullptr : d.well(d.ref_well);
    double prev = kNaN;
    for (auto& h : d.horizons) {
      double depth = kNaN;
      if (rw)
        for (auto& [n, p] : rw->tops)
          if (n == h.name) depth = p.si();
      if (std::isnan(depth) && h.depth.present) depth = h.depth.si();
      if (std::isnan(depth) && h.thickness_above.present) depth = prev + h.thickness_above.si();
      if (std::isnan(depth) && h.iso_thickness.present) {
        int j = d.horizon_index(h.iso_from);
        depth = -z.at(static_cast<std::size_t>(j)) + h.iso_thickness.si();
      }
      z.push_back(-depth);
      prev = depth;
    }
    log.note("flat-layers builder: horizons are horizontal planes");
    return std::make_shared<FlatLayersModel>(std::move(mi), std::move(z));
  }
};
}  // namespace

STRATONIX_REGISTER(flat_layers_builder, registry::Kind::StructureBuilder, "flat-layers", "1.0.0",
                   "Horizontal unfaulted layers (reference implementation of the model interface)",
                   [] { return std::make_unique<FlatLayersBuilder>(); })

}  // namespace stratonix
