#include "stratonix/engines/builder.hpp"
#include "stratonix/engines/faults.hpp"

namespace stratonix {
void register_engine_plugins() {
  STRATONIX_PLUGIN_ANCHOR(flat_layers_builder);
  STRATONIX_PLUGIN_ANCHOR(parametric_builder);
  STRATONIX_PLUGIN_ANCHOR(implicit_builder);
  register_fault_plugins();
}
}  // namespace stratonix
