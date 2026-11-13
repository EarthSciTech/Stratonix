# Layers above the model interface. Levels: engines 3, repr 4, checks/io/figures 5, package 6.

# ---- Layer 3: engines ----------------------------------------------------------------------
add_library(stratonix_engines STATIC
  src/engines/builder.cpp
  src/engines/register.cpp)
target_link_libraries(stratonix_engines PUBLIC stratonix_model stratonix_desc stratonix_registry)
stx_layer(stratonix_engines 3)

# ---- Layer 5: checks and provenance --------------------------------------------------------
add_library(stratonix_checks STATIC
  src/checks/manifest.cpp)
target_link_libraries(stratonix_checks PUBLIC stratonix_engines)
stx_layer(stratonix_checks 5)

# ---- Layer 6: delivery packaging and pipeline ---------------------------------------------
add_library(stratonix_package STATIC
  src/package/pipeline.cpp)
target_link_libraries(stratonix_package PUBLIC stratonix_checks)
stx_layer(stratonix_package 6)
