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
target_include_directories(deltahedron-test PRIVATE ${CMAKE_SOURCE_DIR}/data)
target_compile_definitions(deltahedron-test PRIVATE C44_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
add_test(NAME deltahedron-test COMMAND deltahedron-test)

add_executable(buckinverse-phase1-test buckinverse-phase1-test.cc)
target_link_libraries(buckinverse-phase1-test fullerenes GTest::GTest GTest::Main)
add_test(NAME buckinverse-phase1-test COMMAND buckinverse-phase1-test)

add_executable(buckinverse-phase2-test buckinverse-phase2-test.cc)
target_link_libraries(buckinverse-phase2-test fullerenes GTest::GTest GTest::Main)
add_test(NAME buckinverse-phase2-test COMMAND buckinverse-phase2-test)

add_executable(buckinverse-phase3-test buckinverse-phase3-test.cc)
target_link_libraries(buckinverse-phase3-test fullerenes GTest::GTest GTest::Main)
add_test(NAME buckinverse-phase3-test COMMAND buckinverse-phase3-test)

add_executable(buckinverse-roundtrip-test buckinverse-roundtrip-test.cc)
target_link_libraries(buckinverse-roundtrip-test fullerenes GTest::GTest GTest::Main)
add_test(NAME buckinverse-roundtrip-test COMMAND buckinverse-roundtrip-test)

add_executable(deltahedron-geometry-test deltahedron-geometry-test.cc)
target_link_libraries(deltahedron-geometry-test fullerenes GTest::GTest GTest::Main)
add_test(NAME deltahedron-geometry-test COMMAND deltahedron-geometry-test)

add_executable(convex-hull-test convex-hull-test.cc)
target_link_libraries(convex-hull-test fullerenes GTest::GTest GTest::Main)
add_test(NAME convex-hull-test COMMAND convex-hull-test)
