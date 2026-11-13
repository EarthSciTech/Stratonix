add_executable(stratonix main.cpp cli_commands.cpp)
target_link_libraries(stratonix PRIVATE stratonix_package)
stx_layer(stratonix 7)
install(TARGETS stratonix RUNTIME DESTINATION bin)
if(WIN32)
  # Ship the DLLs of dynamically linked dependencies (Cairo is LGPL and linked dynamically).
  install(FILES $<TARGET_RUNTIME_DLLS:stratonix> DESTINATION bin OPTIONAL)
endif()
