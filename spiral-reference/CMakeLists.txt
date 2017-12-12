cmake_minimum_required(VERSION 3.2.0)
project(spiral_reference)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)   # See below (1)

set(libs)

include_directories(AFTER SYSTEM .)

find_package(GSL)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  list(APPEND libs ${GSL_LIBRARIES})
  add_definitions(-DHAS_GSL)
endif(GSL_FOUND)

add_subdirectory(libgraph)
add_subdirectory(programs)

