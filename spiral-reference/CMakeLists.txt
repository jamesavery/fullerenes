cmake_minimum_required(VERSION 3.2.0)
project(spiral_reference)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)   # See below (1)

include_directories(AFTER SYSTEM .)

add_subdirectory(libgraph)
add_subdirectory(programs)


