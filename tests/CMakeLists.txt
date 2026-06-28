add_executable(gc-transform-test gc-transform-test.cc)
target_link_libraries(gc-transform-test fullerenes GTest::GTest GTest::Main)
add_test(NAME gc-transform-test COMMAND gc-transform-test)

add_executable(surface-geodesics-test surface-geodesics-test.cc)
target_link_libraries(surface-geodesics-test fullerenes buckygen GTest::GTest GTest::Main)
add_test(NAME surface-geodesics-test COMMAND surface-geodesics-test)

add_executable(unfold-straighten-test unfold-straighten-test.cc)
target_link_libraries(unfold-straighten-test fullerenes buckygen GTest::GTest GTest::Main)
add_test(NAME unfold-straighten-test COMMAND unfold-straighten-test)

# Standalone (no GTest): stresses BuckyGen child-process teardown; exit 0 = pass,
# 1 = assertion, 2 = watchdog timeout. CLI-tunable for heavy manual runs.
add_executable(buckygen-stress-test buckygen-stress-test.cc)
target_link_libraries(buckygen-stress-test fullerenes buckygen)
add_test(NAME buckygen-stress-test COMMAND buckygen-stress-test)

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
# Disabled from ctest: takes ~10 minutes. Run manually when needed.
# add_test(NAME deltahedron-geometry-test COMMAND deltahedron-geometry-test)

add_executable(convex-hull-test convex-hull-test.cc)
target_link_libraries(convex-hull-test fullerenes GTest::GTest GTest::Main)
add_test(NAME convex-hull-test COMMAND convex-hull-test)

add_executable(orientation-test orientation-test.cc)
target_link_libraries(orientation-test fullerenes GTest::GTest GTest::Main)
add_test(NAME orientation-test COMMAND orientation-test)

add_executable(ply-io-test ply-io-test.cc)
target_link_libraries(ply-io-test fullerenes GTest::GTest GTest::Main)
add_test(NAME ply-io-test COMMAND ply-io-test)

add_executable(cone-finder-test cone-finder-test.cc)
target_link_libraries(cone-finder-test fullerenes GTest::GTest GTest::Main)
add_test(NAME cone-finder-test COMMAND cone-finder-test)

add_executable(csr-test csr-test.cc)
target_link_libraries(csr-test fullerenes GTest::GTest GTest::Main)
add_test(NAME csr-test COMMAND csr-test)

add_executable(dense-graph-test dense-graph-test.cc)
target_link_libraries(dense-graph-test fullerenes GTest::GTest GTest::Main)
add_test(NAME dense-graph-test COMMAND dense-graph-test)

add_executable(view-test view-test.cc)
# The view/batch layer is header-only (std::vector fallback when SYCL is off),
# so view-test only needs the SYCL runtime lib in a SYCL build.
if(ENABLE_SYCL)
  target_link_libraries(view-test fullerenes sycl_fullerene_lib GTest::GTest GTest::Main)
else()
  target_link_libraries(view-test fullerenes GTest::GTest GTest::Main)
endif()
add_test(NAME view-test COMMAND view-test)

add_executable(delaunay-test delaunay-test.cc)
target_link_libraries(delaunay-test fullerenes GTest::GTest GTest::Main)
target_compile_definitions(delaunay-test PRIVATE FULLERENE_ROOT="${CMAKE_SOURCE_DIR}")
add_test(NAME delaunay-test COMMAND delaunay-test)
add_executable(curvature-flow-test curvature-flow-test.cc)
target_link_libraries(curvature-flow-test fullerenes GTest::GTest GTest::Main)
add_test(NAME curvature-flow-test COMMAND curvature-flow-test)
add_executable(extract-graph extract_graph.cc)
target_link_libraries(extract-graph fullerenes)
target_compile_definitions(extract-graph PRIVATE FULLERENE_ROOT="${CMAKE_SOURCE_DIR}")

add_executable(planargraph-ascii-test planargraph-ascii-test.cc)
target_link_libraries(planargraph-ascii-test fullerenes GTest::GTest GTest::Main)
add_test(NAME planargraph-ascii-test COMMAND planargraph-ascii-test)
