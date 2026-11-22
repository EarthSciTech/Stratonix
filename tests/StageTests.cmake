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
# Stage 3A analytical tests (STX-TN-11), written before the features (Prompt 34), passing since
# Prompt 36.
stx_test(test_v3_data stratonix_package)
set_tests_properties(test_v3_data PROPERTIES TIMEOUT 5400)
# Stage 3B analytical tests (STX-TN-12), written before the features (Prompt 37) and made to pass
# by Prompts 38 and 39 (all passing since Prompt 39).
stx_test(test_v3_grids stratonix_package)
set_tests_properties(test_v3_grids PROPERTIES TIMEOUT 3600)
if(STRATONIX_WITH_FESAPI)
  target_compile_definitions(test_v3_grids PRIVATE STRATONIX_HAVE_FESAPI=1)
endif()
# Stage 3C analytical tests (STX-TN-13), written before the features (Prompt 40) and made to pass
# by Prompt 41 (all passing since Prompt 41).
stx_test(test_v3_uncertainty stratonix_package)
set_tests_properties(test_v3_uncertainty PROPERTIES TIMEOUT 3600)
# Stage 4A analytical tests (STX-TN-14), written before the features (Prompt 42) and made to pass
# by Prompts 43 to 47.
stx_test(test_v4_engine stratonix_package)
set_tests_properties(test_v4_engine PROPERTIES TIMEOUT 3600)
# Stage 4A: the fold-frame command on a sine fold (axis 000, wavelength 4 km)
add_test(NAME cli_fold_frame COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix> -DOUT=${STX_TEST_TMP}/cli_fold_frame.json
         "-DARGS=fold-frame;${CMAKE_CURRENT_SOURCE_DIR}/data/fold_bedding.csv" -P ${CMAKE_CURRENT_SOURCE_DIR}/RunCmd.cmake)
# Stage 4A benchmark: 1 000 000 points within 600 s and 4 GB (release builds; others skip),
# passing since Prompt 43.
stx_test(bench_v4_discrete stratonix_package)
set_tests_properties(bench_v4_discrete PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 1200)
if(STRATONIX_SANITIZE)
  target_compile_definitions(bench_v4_discrete PRIVATE STX_SANITIZE=1)
endif()
# Stage 3A benchmark: 100 000 points within 300 s (release builds; others skip with code 77).
stx_test(bench_implicit_100k stratonix_package)
set_tests_properties(bench_implicit_100k PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900)
if(STRATONIX_SANITIZE)
  target_compile_definitions(bench_implicit_100k PRIVATE STX_SANITIZE=1)
endif()
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
# Stage 3B: the GRDECL of the gridded examples read back by an independent reader (test tool).
if(Python3_Interpreter_FOUND)
  foreach(ex fault_system_relay_splay rift_inversion_history)
    add_test(NAME grdecl_independent_${ex}
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/read_grdecl.py ${STX_TEST_TMP}/cli_build_${ex})
    set_tests_properties(grdecl_independent_${ex} PROPERTIES DEPENDS cli_build_${ex} SKIP_RETURN_CODE 77)
  endforeach()
endif()
# Stage 3C: the 200-member ensemble of the fault-system example, every member's checks recorded,
# and the sensitivity ranking of the analytical case (optimised builds only: 208 model builds).
if(Python3_Interpreter_FOUND AND NOT STRATONIX_SANITIZE AND (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo"))
  add_test(NAME cli_ensemble_fault_system COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix> -DOUT=${STX_TEST_TMP}/cli_ensemble_fault_system
           "-DARGS=ensemble;${PROJECT_SOURCE_DIR}/examples/fault_system_uncertainty.yaml;--no-figures" -P ${CMAKE_CURRENT_SOURCE_DIR}/RunCmd.cmake)
  set_tests_properties(cli_ensemble_fault_system PROPERTIES TIMEOUT 3600)
  add_test(NAME ensemble_ranking_fault_system
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_ensemble.py ${STX_TEST_TMP}/cli_ensemble_fault_system 200 20261115)
  set_tests_properties(ensemble_ranking_fault_system PROPERTIES DEPENDS cli_ensemble_fault_system)
endif()
# The example builds with every check passing (exit code 0).
add_test(NAME cli_build_example COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix>
         -DDESC=${PROJECT_SOURCE_DIR}/examples/synthetic_fold_two_faults.yaml -DOUT=${STX_TEST_TMP}/cli_build_example
         -P ${CMAKE_CURRENT_SOURCE_DIR}/RunBuild.cmake)
# The worked examples of the user guide build with every check passing.
foreach(ex one_fault dome_blind_fault implicit_from_samples layered_faults_sheet listric_rollover fault_system_relay_splay thrust_repeated_section strike_slip_offset oblique_fault trishear_fold fault_bend_fold detachment_folds unconformity_onlap_pinchout growth_fault prograding_wedge salt_pillow_withdrawal salt_diapir_crestal_faults compaction_drape rift_inversion_history implicit_unconformity_from_samples implicit_relay_splay_from_samples implicit_detachment_folds_from_samples implicit_faults_series_from_data hybrid_fault_fold_conditioning fault_system_uncertainty)
  add_test(NAME cli_build_${ex} COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix>
           -DDESC=${PROJECT_SOURCE_DIR}/examples/${ex}.yaml -DOUT=${STX_TEST_TMP}/cli_build_${ex}
           -P ${CMAKE_CURRENT_SOURCE_DIR}/RunBuild.cmake)
endforeach()
foreach(ex implicit_unconformity_from_samples implicit_relay_splay_from_samples implicit_detachment_folds_from_samples implicit_faults_series_from_data)
  set_tests_properties(cli_build_${ex} PROPERTIES TIMEOUT 3600)
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
# Stage 3A: the same commands on an implicit model with series and a fault whose throw comes from the data.
set(_imp3 ${PROJECT_SOURCE_DIR}/examples/implicit_faults_series_from_data.yaml)
add_test(NAME cli_implicit3_section COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix> -DOUT=${STX_TEST_TMP}/cli_implicit3_section
         "-DARGS=section;${_imp3};--line;200,5000\;9800,5000" -P ${CMAKE_CURRENT_SOURCE_DIR}/RunCmd.cmake)
add_test(NAME cli_implicit3_export COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:stratonix> -DOUT=${STX_TEST_TMP}/cli_implicit3_export
         "-DARGS=export;${_imp3};--formats;petrel,vtk,csv,gempy,loopstructural" -P ${CMAKE_CURRENT_SOURCE_DIR}/RunCmd.cmake)
set_tests_properties(cli_implicit3_section cli_implicit3_export PROPERTIES TIMEOUT 3600)
# LeakSanitizer: ignore process-lifetime caches of fontconfig/Cairo (see lsan.supp). Last, so it
# reaches every test defined above.
get_property(_stx_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${_stx_tests} PROPERTIES ENVIRONMENT "LSAN_OPTIONS=suppressions=${CMAKE_CURRENT_SOURCE_DIR}/lsan.supp")
