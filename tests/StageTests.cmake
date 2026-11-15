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
stx_test(test_interop stratonix_package)
stx_test(test_package stratonix_package stx::zlib)
stx_test(test_implicit stratonix_package)
stx_test(test_sheet stratonix_package)
# Stage 2A analytical tests (STX-TN-05), written before the features (Prompt 15) and made to pass
# by Prompts 16 to 18.
stx_test(test_v2_curved stratonix_package)
stx_test(test_v2_interaction stratonix_package)
stx_test(test_v2_movement stratonix_package)
stx_test(test_v2_folding stratonix_package)
stx_test(test_v2_strat stratonix_package)
stx_test(test_v2_salt stratonix_package)
stx_test(test_v2_history stratonix_package)
stx_test(test_v2_faultrock stratonix_package)
set_tests_properties(test_v2_faultrock PROPERTIES WILL_FAIL TRUE)  # until Prompt 33
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
if(Python3_Interpreter_FOUND)
  add_test(NAME neutral_set_independent_check
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_neutral_set.py ${STX_TEST_TMP}/export_interop)
  set_tests_properties(neutral_set_independent_check PROPERTIES DEPENDS test_interop SKIP_RETURN_CODE 77)
endif()
# The example builds with every check passing (exit code 0).
add_test(NAME cli_build_example COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix>
         -DDESC=${PROJECT_SOURCE_DIR}/examples/synthetic_fold_two_faults.yaml -DOUT=${STX_TEST_TMP}/cli_build_example
         -P ${CMAKE_CURRENT_SOURCE_DIR}/RunBuild.cmake)
# The worked examples of the user guide build with every check passing.
foreach(ex one_fault dome_blind_fault implicit_from_samples layered_faults_sheet listric_rollover fault_system_relay_splay thrust_repeated_section strike_slip_offset oblique_fault trishear_fold fault_bend_fold detachment_folds unconformity_onlap_pinchout growth_fault prograding_wedge salt_pillow_withdrawal salt_diapir_crestal_faults compaction_drape rift_inversion_history)
  add_test(NAME cli_build_${ex} COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix>
           -DDESC=${PROJECT_SOURCE_DIR}/examples/${ex}.yaml -DOUT=${STX_TEST_TMP}/cli_build_${ex}
           -P ${CMAKE_CURRENT_SOURCE_DIR}/RunBuild.cmake)
endforeach()
# Every v1 command works unchanged on a data-conditioned (implicit) model.
set(_imp ${PROJECT_SOURCE_DIR}/examples/implicit_from_samples.yaml)
add_test(NAME cli_implicit_section COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix> -DOUT=${STX_TEST_TMP}/cli_implicit_section
         "-DARGS=section;${_imp};--line;500,500\;7500,5500" -P ${CMAKE_CURRENT_SOURCE_DIR}/RunCmd.cmake)
add_test(NAME cli_implicit_export COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix> -DOUT=${STX_TEST_TMP}/cli_implicit_export
         "-DARGS=export;${_imp};--formats;petrel,vtk,csv,gempy" -P ${CMAKE_CURRENT_SOURCE_DIR}/RunCmd.cmake)
add_test(NAME cli_implicit_package COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix> -DOUT=${STX_TEST_TMP}/cli_implicit_package
         "-DARGS=package;${_imp};--formats;petrel;--no-figures;--zip" -P ${CMAKE_CURRENT_SOURCE_DIR}/RunCmd.cmake)
add_test(NAME cli_implicit_suite COMMAND stratonix suite ${_imp})
# LeakSanitizer: ignore process-lifetime caches of fontconfig/Cairo (see lsan.supp). Last, so it
# reaches every test defined above.
get_property(_stx_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${_stx_tests} PROPERTIES ENVIRONMENT "LSAN_OPTIONS=suppressions=${CMAKE_CURRENT_SOURCE_DIR}/lsan.supp")
