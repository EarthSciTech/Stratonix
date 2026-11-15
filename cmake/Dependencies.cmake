# Third-party dependencies. Each is found through vcpkg/CMake config first, then pkg-config,
# then a plain library search, and exposed as an stx:: alias target.
find_package(PkgConfig QUIET)
include(CheckCSourceCompiles)

# A static archive can go into the shared C API library only if it was built as position-
# independent code (Ubuntu 22.04's libz.a and libyaml.a are not). Test by linking a shared object.
function(stx_static_is_pic out lib)
  set(CMAKE_REQUIRED_FLAGS "-fPIC")
  # whole archive: every object must be position-independent, not only the one a probe touches
  set(CMAKE_REQUIRED_LINK_OPTIONS "-shared;-Wl,--whole-archive;${lib};-Wl,--no-whole-archive")
  set(CMAKE_REQUIRED_QUIET ON)
  string(MAKE_C_IDENTIFIER "STX_PIC_ALL_${lib}" var)
  check_c_source_compiles("void stx_probe(void) {}" ${var})
  set(${out} ${${var}} PARENT_SCOPE)
endfunction()

# libyaml (MIT): YAML and JSON reader. Prefer a static library so the program has no
# runtime dependency to install.
find_package(yaml CONFIG QUIET)
if(NOT TARGET yaml AND NOT WIN32)
  find_path(STX_YAML_INC_S yaml.h)
  find_library(STX_YAML_LIB_S NAMES libyaml.a)
  if(STX_YAML_LIB_S)
    stx_static_is_pic(STX_YAML_STATIC_PIC ${STX_YAML_LIB_S})
  endif()
  if(STX_YAML_INC_S AND STX_YAML_LIB_S AND STX_YAML_STATIC_PIC)
    add_library(stx_yaml_static STATIC IMPORTED)
    set_target_properties(stx_yaml_static PROPERTIES IMPORTED_LOCATION ${STX_YAML_LIB_S}
                          INTERFACE_INCLUDE_DIRECTORIES ${STX_YAML_INC_S})
    add_library(stx::yaml ALIAS stx_yaml_static)
  endif()
endif()
if(TARGET stx::yaml)
elseif(TARGET yaml)
  add_library(stx::yaml ALIAS yaml)
else()
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(STX_YAML QUIET IMPORTED_TARGET yaml-0.1)
  endif()
  if(TARGET PkgConfig::STX_YAML)
    add_library(stx::yaml ALIAS PkgConfig::STX_YAML)
  else()
    find_path(STX_YAML_INC yaml.h REQUIRED)
    find_library(STX_YAML_LIB NAMES yaml libyaml REQUIRED)
    add_library(stx_yaml_imported UNKNOWN IMPORTED)
    set_target_properties(stx_yaml_imported PROPERTIES IMPORTED_LOCATION ${STX_YAML_LIB}
                          INTERFACE_INCLUDE_DIRECTORIES ${STX_YAML_INC})
    add_library(stx::yaml ALIAS stx_yaml_imported)
  endif()
endif()

# zlib (zlib licence): PNG, zip and PPTX writers. On Linux and macOS the static archive is used
# if it is position-independent (checked before any target is created); otherwise, and on
# Windows (vcpkg), the library CMake's FindZLIB finds.
if(NOT WIN32)
  find_path(STX_ZLIB_INC_S zlib.h)
  find_library(STX_ZLIB_LIB_S NAMES libz.a)
  if(STX_ZLIB_INC_S AND STX_ZLIB_LIB_S)
    stx_static_is_pic(STX_ZLIB_STATIC_PIC ${STX_ZLIB_LIB_S})
  endif()
endif()
if(STX_ZLIB_STATIC_PIC)
  add_library(stx_zlib_static STATIC IMPORTED)
  set_target_properties(stx_zlib_static PROPERTIES IMPORTED_LOCATION ${STX_ZLIB_LIB_S}
                        INTERFACE_INCLUDE_DIRECTORIES ${STX_ZLIB_INC_S})
  add_library(stx::zlib ALIAS stx_zlib_static)
  message(STATUS "zlib: static ${STX_ZLIB_LIB_S}")
else()
  find_package(ZLIB REQUIRED)
  add_library(stx::zlib ALIAS ZLIB::ZLIB)
  message(STATUS "zlib: ${ZLIB_LIBRARIES}")
endif()

# Cairo (LGPL-2.1 / MPL-1.1): optional raster backend, linked dynamically.
if(STRATONIX_WITH_CAIRO)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(STX_CAIRO QUIET IMPORTED_TARGET cairo)
  endif()
  if(TARGET PkgConfig::STX_CAIRO)
    add_library(stx::cairo ALIAS PkgConfig::STX_CAIRO)
  else()
    find_path(STX_CAIRO_INC cairo.h PATH_SUFFIXES cairo)
    find_library(STX_CAIRO_LIB NAMES cairo)
    if(STX_CAIRO_INC AND STX_CAIRO_LIB)
      add_library(stx_cairo_imported UNKNOWN IMPORTED)
      set_target_properties(stx_cairo_imported PROPERTIES IMPORTED_LOCATION ${STX_CAIRO_LIB}
                            INTERFACE_INCLUDE_DIRECTORIES ${STX_CAIRO_INC})
      add_library(stx::cairo ALIAS stx_cairo_imported)
    else()
      message(WARNING "Cairo not found: PNG and TIFF figures disabled (SVG and PPTX still written)")
      set(STRATONIX_WITH_CAIRO OFF)
    endif()
  endif()
endif()

# FESAPI (Apache 2.0): optional RESQML writer (stage 3B), linked dynamically. Built from source by
# the CI job 'resqml'; FESAPI_ROOT points to its install prefix.
if(STRATONIX_WITH_FESAPI)
  find_path(STX_FESAPI_INC fesapi/common/DataObjectRepository.h HINTS ${FESAPI_ROOT}/include ENV FESAPI_ROOT PATH_SUFFIXES include)
  find_library(STX_FESAPI_LIB NAMES FesapiCpp FesapiCppd HINTS ${FESAPI_ROOT}/lib ${FESAPI_ROOT}/lib64 ENV FESAPI_ROOT PATH_SUFFIXES lib lib64)
  if(NOT STX_FESAPI_INC OR NOT STX_FESAPI_LIB)
    message(FATAL_ERROR "STRATONIX_WITH_FESAPI: FESAPI not found (set FESAPI_ROOT to its install prefix)")
  endif()
  add_library(stx_fesapi_imported UNKNOWN IMPORTED)
  set_target_properties(stx_fesapi_imported PROPERTIES IMPORTED_LOCATION ${STX_FESAPI_LIB} INTERFACE_INCLUDE_DIRECTORIES ${STX_FESAPI_INC})
  add_library(stx::fesapi ALIAS stx_fesapi_imported)
  message(STATUS "FESAPI: ${STX_FESAPI_LIB}")
endif()

if(STRATONIX_WITH_PROJ)
  find_package(PROJ CONFIG REQUIRED)
endif()
