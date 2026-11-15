// Parametric engine: Model interface over Structure + FaultNetwork.
#include "stratonix/engines/parametric.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"

namespace stratonix {

ParametricModel::ParametricModel(const desc::Description& d, const BuildInputs& in, BuildLog& log)
    : info_(model_info_from(d)) {
  info_.engine = "parametric";
  structure_ = std::make_shared<Structure>(d, in, log);
  for (std::size_t u = 1; u < info_.units.size() && u < info_.horizons.size(); ++u) info_.units[u].termination = structure_->unit_termination(u);
  faults_ = std::make_shared<FaultNetwork>(d, *structure_, log);
  for (std::size_t u : structure_->erosional_surfaces()) {
    bool undisplaced = !faults_->empty();
    for (auto& f : faults_->faults())
      undisplaced = undisplaced && (f.tmax == 0 || (!f.growth_f.empty() && u < f.growth_f.size() && f.growth_f[u] == 0));
    if (undisplaced) {
      post_erosion_.push_back(u);
      for (std::size_t k = u + 1; k < info_.units.size() && k < structure_->horizons(); ++k)  // may be truncated
        if (info_.units[k].termination.find("truncation") == std::string::npos)
          info_.units[k].termination = info_.units[k].termination.empty() ? "truncation" : "truncation, " + info_.units[k].termination;
      log.note("horizon " + d.horizons[u].name + ": an erosional surface no fault displaces; it truncates the faulted geometry");
    }
  }
  for (auto& f : faults_->faults()) {
    FaultInfo fi;
    fi.name = f.name;
    fi.index = f.index;
    fi.age = static_cast<int>(f.age);
    fi.ref_point = f.p0;
    fi.strike_dir = f.s;
    fi.dip_dir = f.dd;
    fi.normal = f.n;
    fi.half_length = 0.5 * f.length;
    fi.b_min = std::numeric_limits<double>::infinity();
    fi.b_max = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < f.a_tab.size(); ++i) {
      double bu, bl;
      f.tips(f.a_tab[i], bu, bl);
      fi.b_min = std::min(fi.b_min, bu);
      fi.b_max = std::max(fi.b_max, bl);
    }
    fi.slip_mode = f.kin->name();
    if (f.curved) {
      // Approximate plane: the chord from where the fault meets the model top to its far end.
      fi.geometry = d.faults[static_cast<std::size_t>(f.index)].geom_type;
      double bu, bl;
      f.tips(0, bu, bl);
      Vec3 top = f.at(0, bu), far = f.at(0, bl);
      Vec3 dd = normalized(far - top);
      double dip = std::asin(std::clamp(-dd.z, -1.0, 1.0));
      fi.ref_point = top;
      fi.dip_dir = dd;
      fi.normal = f.dh * std::sin(dip) + Vec3{0, 0, std::cos(dip)};
      fi.b_min = 0;
      fi.b_max = norm(far - top);
      for (int i = 0; i <= 400; ++i) {
        double u = bu + (bl - bu) * i / 400.0;
        fi.profile.emplace_back(u, f.curved->zf(u));
      }
    }
    fi.overlapping = f.curved ? f.curved->fpf && f.curved->heave < 0 : (f.tmax < 0 ? -f.slip.z : f.slip.z) > 1e-12;
    for (double g : f.growth_f) fi.overlapping = fi.overlapping || g < 0;  // inverted layers move up the fault
    const auto& all = faults_->faults();
    if (f.abuts >= 0) fi.relations.push_back((f.transfer ? "splay of " : "abuts ") + all[static_cast<std::size_t>(f.abuts)].name);
    if (f.relay_with >= 0) fi.relations.push_back("relays with " + all[static_cast<std::size_t>(f.relay_with)].name);
    info_.faults.push_back(fi);
  }
}

const void* ParametricModel::capability(std::string_view name) const {
  if (name == "structure") return structure_.get();
  if (name == "stratigraphy") return structure_->has_stratigraphy() ? static_cast<const StratigraphyQuery*>(structure_.get()) : nullptr;
  if (name == "faults") return faults_.get();
  return nullptr;
}

void ParametricModel::eval_point(const Vec3& p, BlockCode* block, Vec3* restored, double* scalars) const {
  Vec3 r = faults_->empty() ? p : faults_->restore(p, block);
  if (!post_erosion_.empty() && above_post_erosion(p)) r = p;
  if (faults_->empty() && block) *block = 0;
  if (restored) *restored = r;
  if (scalars) {
    structure_->z_all(r.x, r.y, scalars);
    for (std::size_t h = 0; h < structure_->horizons(); ++h) scalars[h] = r.z - scalars[h];
    if (!post_erosion_.empty()) post_truncate(scalars);
  }
}

void ParametricModel::evaluate(std::span<const Vec3> pts, unsigned fields, EvalResult& out) const {
  const std::size_t n = pts.size(), nh = structure_->horizons();
  out.n = n;
  out.unit.assign(fields & kEvalUnit ? n : 0, 0);
  out.block.assign(fields & kEvalBlock ? n : 0, 0);
  out.scalar.assign(fields & kEvalScalars ? n * nh : 0, 0.0);
  auto body = [&](std::size_t a, std::size_t b) {
    std::vector<double> s(nh);
    for (std::size_t i = a; i < b; ++i) {
      BlockCode bc = 0;
      eval_point(pts[i], &bc, nullptr, s.data());
      if (fields & kEvalBlock) out.block[i] = bc;
      if (fields & kEvalScalars)
        for (std::size_t h = 0; h < nh; ++h) out.scalar[h * n + i] = s[h];
      if (fields & kEvalUnit) out.unit[i] = unit_from_scalars(s.data(), 1, nh);
    }
  };
  if (optimised_kernels() && n >= 2048) parallel_for(n, 1024, body);
  else body(0, n);
}

void ParametricModel::scalar_in_block(std::span<const Vec3> pts, BlockCode block, int horizon, std::span<double> out) const {
  if (horizon < 0 || static_cast<std::size_t>(horizon) >= structure_->horizons())
    throw Error(ErrorCode::InvalidArgument, "horizon index out of range");
  auto body = [&](std::size_t a, std::size_t b) {
    for (std::size_t i = a; i < b; ++i) {
      Vec3 r = faults_->empty() ? pts[i] : faults_->restore(pts[i], nullptr, &block);
      if (!post_erosion_.empty() && above_post_erosion(pts[i])) r = pts[i];
      out[i] = r.z - structure_->z(static_cast<std::size_t>(horizon), r.x, r.y);
      for (std::size_t u : post_erosion_)
        if (static_cast<std::size_t>(horizon) > u) out[i] = std::max(out[i], r.z - structure_->z(u, r.x, r.y));
    }
  };
  if (optimised_kernels() && pts.size() >= 2048) parallel_for(pts.size(), 1024, body);
  else body(0, pts.size());
}

void ParametricModel::evaluate_fault(int f, std::span<const Vec3> pts, FaultEval& out) const {
  if (f < 0 || static_cast<std::size_t>(f) >= faults_->size()) throw Error(ErrorCode::InvalidArgument, "fault index out of range");
  const Fault& F = faults_->faults()[static_cast<std::size_t>(f)];
  out.distance.resize(pts.size());
  out.extent.resize(pts.size());
  out.displacement.resize(pts.size());
  auto body = [&](std::size_t a0, std::size_t b0) {
    for (std::size_t i = a0; i < b0; ++i) {
      Vec3 q = faults_->to_fault_frame(pts[i], f);
      double a, b;
      F.profile_coords(q, a, b);
      out.distance[i] = F.sigma(q);
      out.extent[i] = F.extent_coord(a, b);
      bool active = true;
      if (faults_->layered()) {
        // the material's level comes from the full restoration of the present point
        out.displacement[i] = faults_->layer_displacement(f, pts[i], &active);
      } else if (F.curved) {
        out.displacement[i] = F.sigma(q) > 0 ? norm(q - faults_->undo(f, q, 1)) : 0.0;
      } else {
        Vec3 r = faults_->relative_displacement(f, q, &active);
        out.displacement[i] = active ? norm(r) : 0.0;
      }
    }
  };
  if (optimised_kernels() && pts.size() >= 2048) parallel_for(pts.size(), 1024, body);
  else body(0, pts.size());
}

void ParametricModel::restore(std::span<const Vec3> pts, std::span<Vec3> out) const {
  for (std::size_t i = 0; i < pts.size(); ++i) {
    out[i] = faults_->empty() ? pts[i] : faults_->restore(pts[i], nullptr);
    if (!post_erosion_.empty() && above_post_erosion(pts[i])) out[i] = pts[i];
  }
}

namespace {
class ParametricBuilder final : public StructureBuilder {
 public:
  bool accepts(const desc::Description& d) const override { return !d.implicit.present; }
  ModelPtr build(const desc::Description& d, const BuildInputs& in, BuildLog& log) const override {
    return std::make_shared<ParametricModel>(d, in, log);
  }
};
}  // namespace

STRATONIX_REGISTER(parametric_builder, registry::Kind::StructureBuilder, "parametric", "1.0.0",
                   "Fold primitives, imported surfaces, isochores, faults with kinematics and frame rule",
                   [] { return std::make_unique<ParametricBuilder>(); })

}  // namespace stratonix
