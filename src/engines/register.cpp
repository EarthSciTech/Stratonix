#include "stratonix/engines/builder.hpp"

namespace stratonix {
void register_engine_plugins() {
  STRATONIX_PLUGIN_ANCHOR(flat_layers_builder);
}
}  // namespace stratonix
