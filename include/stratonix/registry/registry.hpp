/// @file registry.hpp
/// @brief Plugin registry. v1: compile-time self-registration. Later: run-time loading through
/// the C ABI in plugin_abi.h.
#pragma once
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace stratonix::registry {

/// Plugin kinds known to the registry.
enum class Kind { StructureBuilder, Kinematics, HangingWall, Exporter, Importer, Extractor, Check };
const char* kind_name(Kind k);

/// Base class of every plugin object.
struct PluginBase {
  virtual ~PluginBase() = default;
};

/// A registered plugin: kind, name, semantic version, one-line description and a factory.
struct Entry {
  Kind kind;
  std::string name;
  std::string version;
  std::string description;
  std::string origin = "built-in";  ///< "built-in" or the path of a plugin library
  std::function<std::unique_ptr<PluginBase>()> make;
};

/// Process-wide registry. Thread-safe.
class Registry {
 public:
  static Registry& global();
  /// Add an entry; a second registration of the same kind and name is ignored.
  void add(Entry e);
  const Entry* find(Kind k, const std::string& name) const;
  std::vector<const Entry*> list(Kind k) const;
  std::vector<const Entry*> all() const;
  /// Create a plugin of type T (derived from PluginBase); throws Error if missing or wrong type.
  template <class T>
  std::unique_ptr<T> create(Kind k, const std::string& name) const {
    auto* e = find(k, name);
    if (!e) missing(k, name);
    auto p = e->make();
    auto* t = dynamic_cast<T*>(p.get());
    if (!t) wrong_type(k, name);
    p.release();
    return std::unique_ptr<T>(t);
  }

 private:
  [[noreturn]] static void missing(Kind k, const std::string& name);
  [[noreturn]] static void wrong_type(Kind k, const std::string& name);
  mutable std::mutex m_;
  std::vector<Entry> entries_;
};

/// Helper object whose constructor registers an entry (used by STRATONIX_REGISTER).
struct Registrar {
  explicit Registrar(Entry e) { Registry::global().add(std::move(e)); }
};

// Built-in plugins are registered by stratonix::register_all_plugins() (stratonix/package/plugins.hpp),
// which lives in the top library layer because it references every layer's plugins.

}  // namespace stratonix::registry

/// Self-registration. Use once per plugin in its .cpp file. `id` must be a unique C identifier;
/// the macro also defines the C-linkage function stratonix_plugin_anchor_<id>(), which the
/// layer's register function calls so the linker keeps the object. C linkage makes the anchor
/// the same symbol whichever namespace declares or defines it (MSVC and GCC agree).
#define STRATONIX_REGISTER(id, kind, name, version, description, factory)                       \
  namespace {                                                                                    \
  const ::stratonix::registry::Registrar stratonix_registrar_##id{                              \
      ::stratonix::registry::Entry{kind, name, version, description, "built-in", factory}};      \
  }                                                                                              \
  extern "C" void stratonix_plugin_anchor_##id() {}

/// Declares an anchor at file scope (outside any function) in the file that calls it.
#define STRATONIX_DECLARE_PLUGIN_ANCHOR(id) extern "C" void stratonix_plugin_anchor_##id()

/// Calls an anchor declared with STRATONIX_DECLARE_PLUGIN_ANCHOR.
#define STRATONIX_PLUGIN_ANCHOR(id) ::stratonix_plugin_anchor_##id()
