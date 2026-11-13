add_library(stratonix_capi STATIC src/stratonix_c.cpp src/capi_extensions.cpp)
target_include_directories(stratonix_capi PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                          PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(stratonix_capi PUBLIC stratonix_package)
stx_layer(stratonix_capi 7)

# Shared library with the C interface (for other languages later).
add_library(stratonix_c SHARED src/stratonix_c.cpp src/capi_extensions.cpp)
target_include_directories(stratonix_c PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(stratonix_c PRIVATE stratonix_package)
set_target_properties(stratonix_c PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1
                      WINDOWS_EXPORT_ALL_SYMBOLS ON POSITION_INDEPENDENT_CODE ON)
stx_layer(stratonix_c 7)
install(TARGETS stratonix_c RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES include/stratonix_c.h DESTINATION include)

if(STRATONIX_BUILD_TESTS)
  add_executable(capi_test tests/capi_test.c)
  target_link_libraries(capi_test PRIVATE stratonix_capi)
  set_target_properties(capi_test PROPERTIES LINKER_LANGUAGE CXX)
endif()
