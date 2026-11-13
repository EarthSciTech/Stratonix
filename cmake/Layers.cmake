# Layers above the model interface. Levels: engines 3, repr 4, checks/io/figures 5, package 6.

# ---- Layer 3: engines ----------------------------------------------------------------------
add_library(stratonix_engines STATIC
  src/engines/builder.cpp
  src/engines/structure.cpp
  src/engines/faults.cpp
  src/engines/parametric.cpp
  src/engines/wells.cpp
  src/engines/register.cpp)
target_link_libraries(stratonix_engines PUBLIC stratonix_model stratonix_desc stratonix_registry)
stx_layer(stratonix_engines 3)

# ---- Layer 4: representations --------------------------------------------------------------
add_library(stratonix_repr STATIC
  src/repr/surfaces.cpp
  src/repr/faults_wells.cpp
  src/repr/closure.cpp)
target_link_libraries(stratonix_repr PUBLIC stratonix_model)
stx_layer(stratonix_repr 4)

# ---- Layer 5: checks and provenance --------------------------------------------------------
add_library(stratonix_checks STATIC
  src/checks/manifest.cpp
  src/checks/checks.cpp)
target_link_libraries(stratonix_checks PUBLIC stratonix_engines stratonix_repr)
stx_layer(stratonix_checks 5)

# ---- Layer 5: input/output plugins ---------------------------------------------------------
add_library(stratonix_io STATIC
  src/io/grids.cpp
  src/io/register.cpp)
target_link_libraries(stratonix_io PUBLIC stratonix_engines stratonix_repr)
stx_layer(stratonix_io 5)

# ---- Layer 6: delivery packaging and pipeline ---------------------------------------------
add_library(stratonix_package STATIC
  src/package/pipeline.cpp
  src/package/bench.cpp)
target_link_libraries(stratonix_package PUBLIC stratonix_checks stratonix_io)
stx_layer(stratonix_package 6)
