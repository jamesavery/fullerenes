add_executable(test_csr test_csr.cc)
target_link_libraries(test_csr fullerenes GTest::gtest_main)
gtest_discover_tests(test_csr)
