cmake_minimum_required(VERSION 3.5)
project(buckinverse CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the fullerenes library from the main project build
set(FULLERENES_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)
set(FULLERENES_INCLUDE ${FULLERENES_ROOT}/include)
set(FULLERENES_LIBDIR ${FULLERENES_ROOT}/build2/src/c++)

add_library(buckinverse buckinverse.cc)
target_include_directories(buckinverse PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FULLERENES_INCLUDE})
target_link_libraries(buckinverse
    ${FULLERENES_LIBDIR}/libfullerenes.so
    gfortran)
target_link_directories(buckinverse PUBLIC ${FULLERENES_LIBDIR})

# Phase 1 tests
add_executable(test_phase1 test_phase1.cc)
target_link_libraries(test_phase1 buckinverse)
set_target_properties(test_phase1 PROPERTIES
    BUILD_RPATH ${FULLERENES_LIBDIR})

# Phase 2 tests
add_executable(test_phase2 test_phase2.cc)
target_link_libraries(test_phase2 buckinverse)
set_target_properties(test_phase2 PROPERTIES
    BUILD_RPATH ${FULLERENES_LIBDIR})

enable_testing()
add_test(NAME phase1 COMMAND test_phase1)
add_test(NAME phase2 COMMAND test_phase2)
