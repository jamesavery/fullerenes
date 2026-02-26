add_executable(gc-transform-test gc-transform-test.cc)
target_link_libraries(gc-transform-test fullerenes GTest::GTest GTest::Main)
add_test(NAME gc-transform-test COMMAND gc-transform-test)

add_executable(spiral-roundtrip-test spiral-roundtrip-test.cc)
target_link_libraries(spiral-roundtrip-test fullerenes GTest::GTest GTest::Main)
add_test(NAME spiral-roundtrip-test COMMAND spiral-roundtrip-test)

add_executable(symmetry-test symmetry-test.cc)
target_link_libraries(symmetry-test fullerenes GTest::GTest)
add_test(NAME symmetry-test COMMAND symmetry-test)

add_executable(deltahedron-test deltahedron-test.cc)
target_link_libraries(deltahedron-test fullerenes GTest::GTest GTest::Main)
add_test(NAME deltahedron-test COMMAND deltahedron-test)
