# Suppress -soname: pgacclnk (nvc++ CUDA linker wrapper) misinterprets the
# soname value as an input file, causing "same as output file" linker errors.
set(CMAKE_SHARED_LIBRARY_SONAME_CXX_FLAG "")

set(SYCL_FULLERENE_SOURCES
  dualize.cc
  forcefield-optimize.cc
  tutte.cc
  spherical-projection.cc
  hessian.cc
  eigen.cc
  sycl-util-impl.cc
  geometry-properties.cc
  nop_kernel.cc
)

add_library(sycl_fullerene_lib SHARED ${SYCL_FULLERENE_SOURCES})
target_compile_options(sycl_fullerene_lib PRIVATE ${SYCL_FLAGS})
target_link_options(sycl_fullerene_lib PRIVATE ${SYCL_FLAGS})
target_link_libraries(sycl_fullerene_lib PRIVATE ${SYCL_RUNTIME_LIBS})
set_target_properties(sycl_fullerene_lib PROPERTIES NO_SONAME 1)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(forcefield-optimize.cc hessian.cc
    PROPERTIES COMPILE_OPTIONS "-fpermissive")
endif()

if(APPLE)
  target_link_libraries(sycl_fullerene_lib PRIVATE fullerenes)
endif()
