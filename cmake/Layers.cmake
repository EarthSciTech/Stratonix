# Layers above the model interface. Levels: engines 3, repr 4, checks/io/figures 5, package 6.

# ---- Layer 3: engines ----------------------------------------------------------------------
add_library(stratonix_engines STATIC
  src/engines/builder.cpp
  src/engines/structure.cpp
  src/engines/history.cpp
  src/engines/faults.cpp
  src/engines/parametric.cpp
  src/engines/implicit.cpp
  src/engines/conditioning.cpp
  src/engines/wells.cpp
  src/engines/register.cpp)
target_link_libraries(stratonix_engines PUBLIC stratonix_model stratonix_desc stratonix_registry)
stx_layer(stratonix_engines 3)

# ---- Layer 4: representations --------------------------------------------------------------
add_library(stratonix_repr STATIC
  src/repr/surfaces.cpp
  src/repr/faultrock.cpp
  src/repr/faults_wells.cpp
  src/repr/closure.cpp
  src/repr/grid.cpp)
target_link_libraries(stratonix_repr PUBLIC stratonix_model stratonix_registry)
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
  src/io/formats.cpp
  src/io/exporters.cpp
  src/io/neutral.cpp
  src/io/importers.cpp
  src/io/grid_export.cpp
  src/io/register.cpp)
target_link_libraries(stratonix_io PUBLIC stratonix_engines stratonix_repr)
if(STRATONIX_WITH_FESAPI)
  target_sources(stratonix_io PRIVATE src/io/resqml_fesapi.cpp)
  target_link_libraries(stratonix_io PRIVATE stx::fesapi)
  target_compile_definitions(stratonix_io PRIVATE STRATONIX_HAVE_FESAPI=1)
endif()
stx_layer(stratonix_io 5)

# ---- Layer 5: figures ----------------------------------------------------------------------
add_library(stratonix_figures STATIC
  src/figures/canvas.cpp
  src/figures/raster.cpp
  src/figures/pptx.cpp
  src/figures/render3d.cpp
  src/figures/plots.cpp
  src/figures/sheet.cpp)
target_link_libraries(stratonix_figures PUBLIC stratonix_repr PRIVATE stx::zlib)
if(STRATONIX_WITH_CAIRO)
  target_link_libraries(stratonix_figures PRIVATE stx::cairo)
  target_compile_definitions(stratonix_figures PRIVATE STRATONIX_HAVE_CAIRO=1)
endif()
stx_layer(stratonix_figures 5)

# ---- Layer 6: delivery packaging and pipeline ---------------------------------------------
add_library(stratonix_package STATIC
  src/package/pipeline.cpp
  src/package/bench.cpp
  src/package/delivery.cpp
  src/package/uncertainty.cpp)
target_link_libraries(stratonix_package PUBLIC stratonix_checks stratonix_io stratonix_figures)
stx_layer(stratonix_package 6)
