# Unit, integration and regression tests. Each test file is one executable run by CTest.
set(STX_TEST_TMP ${CMAKE_CURRENT_BINARY_DIR}/tmp)
file(MAKE_DIRECTORY ${STX_TEST_TMP})

function(stx_test name)
  add_executable(${name} unit/${name}.cpp)
  target_link_libraries(${name} PRIVATE ${ARGN})
  target_compile_definitions(${name} PRIVATE
    STX_SOURCE_DIR="${PROJECT_SOURCE_DIR}" STX_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
    STX_TMP_DIR="${STX_TEST_TMP}")
  add_test(NAME ${name} COMMAND ${name})
endfunction()

stx_test(test_core stratonix_core)
stx_test(test_desc stratonix_desc)
stx_test(test_model stratonix_package)
include(${CMAKE_CURRENT_SOURCE_DIR}/StageTests.cmake OPTIONAL)

# Architecture: lower layers never link higher ones.
add_test(NAME layer_check COMMAND ${CMAKE_COMMAND} -DTABLE=${CMAKE_BINARY_DIR}/layer_dependencies.txt
         -P ${PROJECT_SOURCE_DIR}/cmake/LayerCheckScript.cmake)

# C API
add_test(NAME capi_test COMMAND ${CMAKE_COMMAND} -E rm -rf ${STX_TEST_TMP}/capi)
set_tests_properties(capi_test PROPERTIES FIXTURES_SETUP capi_tmp)
add_test(NAME capi_test_run COMMAND capi_test ${PROJECT_SOURCE_DIR}/examples/flat_three_layers.yaml
         ${PROJECT_SOURCE_DIR}/examples/broken/broken_missing_evidence.yaml
         ${PROJECT_SOURCE_DIR}/examples/synthetic_fold_two_faults.yaml ${STX_TEST_TMP}/capi)
set_tests_properties(capi_test_run PROPERTIES FIXTURES_REQUIRED capi_tmp)
add_test(NAME capi_example COMMAND stratonix_capi_example ${PROJECT_SOURCE_DIR}/examples/synthetic_fold_two_faults.yaml)

# Command line: validate passes on the example and fails clearly on the three broken files.
set(STX_EX ${PROJECT_SOURCE_DIR}/examples)
add_test(NAME cli_validate_example COMMAND stratonix validate ${STX_EX}/synthetic_fold_two_faults.yaml)
foreach(b missing_evidence overlapping_units tip_below_base)
  add_test(NAME cli_validate_broken_${b} COMMAND stratonix validate ${STX_EX}/broken/broken_${b}.yaml)
  set_tests_properties(cli_validate_broken_${b} PROPERTIES WILL_FAIL TRUE)
endforeach()
add_test(NAME cli_version COMMAND stratonix version)
add_test(NAME cli_plugins COMMAND stratonix plugins)
