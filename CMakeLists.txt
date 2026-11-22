cmake_minimum_required(VERSION 3.24)
project(stratonix VERSION 2026.3 DESCRIPTION "Stratonix geomodelling toolkit" LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(STRATONIX_WITH_CAIRO "Rasterise figures (PNG, TIFF) with Cairo" ON)
option(STRATONIX_WITH_PROJ "Enable coordinate reprojection with PROJ" OFF)
option(STRATONIX_WITH_SQLITE "Read GeoPackage maps through SQLite when it is found" ON)
option(STRATONIX_WITH_FESAPI "RESQML 2.0.1 grid export through FESAPI (FESAPI_ROOT: its install prefix)" OFF)
option(STRATONIX_SANITIZE "Build with address and undefined-behaviour sanitisers" OFF)
option(STRATONIX_BUILD_TESTS "Build tests" ON)
option(STRATONIX_BUILD_BENCHMARKS "Build benchmarks" ON)
option(STRATONIX_STATIC_RUNTIME "Link the C++ runtime statically where possible" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(StratonixOptions)
include(Dependencies)
include(LayerCheck)
include(EmbedFile)

# Build information compiled into the library.
execute_process(COMMAND git rev-parse --short=12 HEAD WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_VARIABLE STX_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT STX_GIT_COMMIT)
  set(STX_GIT_COMMIT "unknown")
endif()

# ---- Layer 0: core -------------------------------------------------------------------------
add_library(stratonix_core STATIC src/core/core.cpp src/core/well.cpp src/core/linalg.cpp src/core/faultsystems.cpp)
target_include_directories(stratonix_core PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_definitions(stratonix_core PRIVATE
  STRATONIX_VERSION="${PROJECT_VERSION}" STRATONIX_GIT_COMMIT="${STX_GIT_COMMIT}"
  STRATONIX_BUILD_TYPE="$<CONFIG>")
find_package(Threads REQUIRED)
target_link_libraries(stratonix_core PUBLIC Threads::Threads)
if(WIN32)
  target_link_libraries(stratonix_core PUBLIC psapi)
endif()
if(STRATONIX_WITH_PROJ)
  target_link_libraries(stratonix_core PRIVATE PROJ::proj)
  target_compile_definitions(stratonix_core PRIVATE STRATONIX_WITH_PROJ=1)
endif()
stx_layer(stratonix_core 0)

# ---- Layer 1: description and registry -----------------------------------------------------
stx_embed_file(${CMAKE_CURRENT_SOURCE_DIR}/schema/stratonix_schema_v1.json
               ${CMAKE_CURRENT_BINARY_DIR}/generated/schema_embed.cpp
               "stratonix::desc" "schema_v1_text")
add_library(stratonix_desc STATIC src/desc/document.cpp src/desc/description.cpp
            ${CMAKE_CURRENT_BINARY_DIR}/generated/schema_embed.cpp)
target_link_libraries(stratonix_desc PUBLIC stratonix_core PRIVATE stx::yaml)
stx_layer(stratonix_desc 1)

add_library(stratonix_registry STATIC src/registry/registry.cpp)
target_link_libraries(stratonix_registry PUBLIC stratonix_core)
stx_layer(stratonix_registry 1)

# ---- Layer 2: model interface ---------------------------------------------------------------
add_library(stratonix_model STATIC src/model/model.cpp)
target_link_libraries(stratonix_model PUBLIC stratonix_core)
stx_layer(stratonix_model 2)

# ---- Higher layers (added stage by stage) ---------------------------------------------------
include(cmake/Layers.cmake OPTIONAL)

# ---- C API and application ------------------------------------------------------------------
add_subdirectory(capi)
add_subdirectory(apps/stratonix)

stx_check_layers()

if(STRATONIX_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
if(STRATONIX_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

include(Packaging)
