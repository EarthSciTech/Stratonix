# Tests added stage by stage.
stx_test(test_engines stratonix_package)
add_executable(test_regression unit/test_regression.cpp)
target_link_libraries(test_regression PRIVATE stratonix_package)
add_test(NAME regression_prototype COMMAND test_regression)
set_tests_properties(regression_prototype PROPERTIES SKIP_RETURN_CODE 77)
