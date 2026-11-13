#include "stratonix/io/io.hpp"

namespace stratonix::io {
void register_grid_importers();
void register_io_plugins() { register_grid_importers(); }
}  // namespace stratonix::io
