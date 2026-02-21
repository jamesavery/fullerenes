add_executable(gc-transform-test gc-transform-test.cc)
target_link_libraries(gc-transform-test fullerenes GTest::GTest GTest::Main)
add_test(NAME gc-transform-test COMMAND gc-transform-test)
