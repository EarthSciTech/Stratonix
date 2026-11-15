#include "stratonix/registry/registry.hpp"

#include <algorithm>

#include "stratonix/core/error.hpp"

namespace stratonix::registry {

const char* kind_name(Kind k) {
  switch (k) {
    case Kind::StructureBuilder: return "structure-builder";
    case Kind::Kinematics: return "kinematics";
    case Kind::HangingWall: return "hanging-wall";
    case Kind::Exporter: return "exporter";
    case Kind::Importer: return "importer";
    case Kind::Extractor: return "extractor";
    case Kind::Upscaler: return "upscaler";
    default: return "check";
  }
}

Registry& Registry::global() {
  static Registry r;
  return r;
}
void Registry::add(Entry e) {
  std::lock_guard<std::mutex> lk(m_);
  for (auto& x : entries_)
    if (x.kind == e.kind && x.name == e.name) return;
  entries_.push_back(std::move(e));
}
const Entry* Registry::find(Kind k, const std::string& name) const {
  std::lock_guard<std::mutex> lk(m_);
  for (auto& x : entries_)
    if (x.kind == k && x.name == name) return &x;
  return nullptr;
}
std::vector<const Entry*> Registry::list(Kind k) const {
  std::vector<const Entry*> r;
  for (auto* e : all())
    if (e->kind == k) r.push_back(e);
  return r;
}
std::vector<const Entry*> Registry::all() const {
  std::lock_guard<std::mutex> lk(m_);
  std::vector<const Entry*> r;
  for (auto& x : entries_) r.push_back(&x);
  std::sort(r.begin(), r.end(), [](const Entry* a, const Entry* b) {
    if (a->kind != b->kind) return static_cast<int>(a->kind) < static_cast<int>(b->kind);
    return a->name < b->name;
  });
  return r;
}
void Registry::missing(Kind k, const std::string& name) {
  throw Error(ErrorCode::InvalidArgument, std::string("no ") + kind_name(k) + " plugin named '" + name + "'");
}
void Registry::wrong_type(Kind k, const std::string& name) {
  throw Error(ErrorCode::Internal, std::string(kind_name(k)) + " plugin '" + name + "' has the wrong type");
}

}  // namespace stratonix::registry
