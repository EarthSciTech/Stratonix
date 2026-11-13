# Tests added stage by stage.
stx_test(test_engines stratonix_package)
add_executable(test_regression unit/test_regression.cpp)
target_link_libraries(test_regression PRIVATE stratonix_package)
add_test(NAME regression_prototype COMMAND test_regression)
set_tests_properties(regression_prototype PROPERTIES SKIP_RETURN_CODE 77)
stx_test(test_repr stratonix_package)
stx_test(test_parity stratonix_package)
stx_test(test_figures stratonix_package)
stx_test(test_export stratonix_package)
# Independent check of the figure files (Python: Pillow, tifffile, python-pptx); skipped if absent.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  add_test(NAME figures_independent_check
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_figures.py ${STX_TEST_TMP}/figures_synthetic/figures)
  set_tests_properties(figures_independent_check PROPERTIES DEPENDS test_figures SKIP_RETURN_CODE 77)
endif()
if(Python3_Interpreter_FOUND)
  add_test(NAME petrel_set_independent_check
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_petrel_set.py ${STX_TEST_TMP}/export_petrel_a)
  set_tests_properties(petrel_set_independent_check PROPERTIES DEPENDS test_export SKIP_RETURN_CODE 77)
endif()
# The example builds with every check passing (exit code 0).
add_test(NAME cli_build_example COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix>
         -DDESC=${PROJECT_SOURCE_DIR}/examples/synthetic_fold_two_faults.yaml -DOUT=${STX_TEST_TMP}/cli_build_example
         -P ${CMAKE_CURRENT_SOURCE_DIR}/RunBuild.cmake)
# LeakSanitizer: ignore process-lifetime caches of fontconfig/Cairo (see lsan.supp).
get_property(_stx_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${_stx_tests} PROPERTIES ENVIRONMENT "LSAN_OPTIONS=suppressions=${CMAKE_CURRENT_SOURCE_DIR}/lsan.supp")
