add_executable(stratonix main.cpp cli_commands.cpp)
target_link_libraries(stratonix PRIVATE stratonix_package)
stx_layer(stratonix 7)
install(TARGETS stratonix RUNTIME DESTINATION bin)
if(WIN32)
  # Ship the DLLs of dynamically linked dependencies (Cairo is LGPL and linked dynamically).
  install(FILES $<TARGET_RUNTIME_DLLS:stratonix> DESTINATION bin OPTIONAL)
endif()
# Stage 4C (STX-TN-16 section 4): release packages carry the RESQML exporter. FESAPI's library is
# shipped beside the executable (lib on Linux, found through the RPATH; bin on Windows), with its
# own dependencies on Windows (HDF5, minizip, zlib from vcpkg: STRATONIX_RUNTIME_DIRS). On Linux,
# HDF5 and minizip come from the system (libhdf5, libminizip), as stated in the README.
if(STRATONIX_WITH_FESAPI)
  # FESAPI's licence text (Apache 2.0) goes with its binary (copied into its prefix by CI)
  if(EXISTS "${FESAPI_ROOT}/LICENSE")
    install(FILES "${FESAPI_ROOT}/LICENSE" DESTINATION share/doc/stratonix RENAME FESAPI_LICENSE.txt)
  endif()
  if(WIN32)
    set(STX_RT_DIRS "${STRATONIX_RUNTIME_DIRS}")
    get_filename_component(STX_FESAPI_DLL_DIR "${STX_FESAPI_DLL}" DIRECTORY)
    install(CODE "
      file(GET_RUNTIME_DEPENDENCIES LIBRARIES \"${STX_FESAPI_DLL}\" RESOLVED_DEPENDENCIES_VAR stx_deps UNRESOLVED_DEPENDENCIES_VAR stx_unres
           DIRECTORIES \"${STX_FESAPI_DLL_DIR}\" ${STX_RT_DIRS}
           PRE_EXCLUDE_REGEXES \"^api-ms-.*\" \"^ext-ms-.*\" POST_EXCLUDE_REGEXES \".*[Ss][Yy][Ss][Tt][Ee][Mm]32.*\")
      file(INSTALL \${stx_deps} DESTINATION \"\${CMAKE_INSTALL_PREFIX}/bin\")
      if(stx_unres)
        message(STATUS \"FESAPI dependencies not found (not shipped): \${stx_unres}\")
      endif()")
  else()
    set_target_properties(stratonix PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
    get_filename_component(STX_FESAPI_REAL "${STX_FESAPI_LIB}" REALPATH)
    file(GLOB STX_FESAPI_SO "${STX_FESAPI_LIB}*")
    install(FILES ${STX_FESAPI_SO} ${STX_FESAPI_REAL} DESTINATION lib)
  endif()
endif()
