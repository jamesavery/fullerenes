# The nvc++ CUDA linker wrapper (pgacclnk) misinterprets the -soname value as an
# input file ("same as output file" errors), so the SONAME must be suppressed for
# the cuda-nvcxx backend ONLY. For every other backend (X86/omp, AMD, generic) a
# normal SONAME is required: without it the executables record a RELATIVE
# DT_NEEDED path to libsycl_fullerene_lib.so and fail to load depending on the
# working directory (e.g. mdspan-test "cannot open ../../src/sycl/...").
set(SYCL_SUPPRESS_SONAME OFF)
if(SYCL_TARGETS STREQUAL "NVIDIA" AND SYCL_CUDA_BACKEND STREQUAL "cuda-nvcxx")
  set(SYCL_SUPPRESS_SONAME ON)
  set(CMAKE_SHARED_LIBRARY_SONAME_CXX_FLAG "")
endif()

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
if(SYCL_SUPPRESS_SONAME)
  set_target_properties(sycl_fullerene_lib PROPERTIES NO_SONAME 1)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(forcefield-optimize.cc hessian.cc
    PROPERTIES COMPILE_OPTIONS "-fpermissive")
endif()

if(APPLE)
  target_link_libraries(sycl_fullerene_lib PRIVATE fullerenes)
endif()
