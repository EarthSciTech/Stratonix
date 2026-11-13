#include "stratonix/io/export.hpp"
#include "stratonix/io/io.hpp"

namespace stratonix::io {
void register_grid_importers();
void register_io_plugins() {
  register_grid_importers();
  register_exporters();
  register_importers();
}
}  // namespace stratonix::io
