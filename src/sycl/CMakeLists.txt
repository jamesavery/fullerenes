# The nvc++ CUDA linker wrapper (pgacclnk) misinterprets the -soname value as an
# input file ("same as output file" errors), so the SONAME must be suppressed for
# the cuda-nvcxx backend ONLY. For every other backend (X86/omp, AMD, generic) a
# normal SONAME is required: without it the executables record a RELATIVE
# DT_NEEDED path to libsycl_fullerene_lib.so and fail to load depending on the
# working directory (e.g. mdspan-test "cannot open ../../src/sycl/...").
set(SYCL_SUPPRESS_SONAME OFF)
if(SYCL_HAS_NVCXX)
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
# SYCL_TARGETS_FLAG must stay quoted: it may carry an escaped ';' separating
# AdaptiveCpp backends, which an unquoted expansion would split.
target_compile_options(sycl_fullerene_lib PRIVATE ${SYCL_FLAGS} "${SYCL_TARGETS_FLAG}")
target_link_options(sycl_fullerene_lib PRIVATE ${SYCL_FLAGS} "${SYCL_TARGETS_FLAG}")
target_link_libraries(sycl_fullerene_lib PRIVATE ${SYCL_RUNTIME_LIBS})
if(SYCL_SUPPRESS_SONAME)
  set_target_properties(sycl_fullerene_lib PROPERTIES NO_SONAME 1)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(forcefield-optimize.cc hessian.cc
    PROPERTIES COMPILE_OPTIONS "-fpermissive")
endif()

# Per-translation-unit optimization override, for bisecting codegen bugs:
#
#   cmake -DSYCL_TU_OPT_FILE=forcefield-optimize.cc -DSYCL_TU_OPT_LEVEL=-O1 ...
#
# This is how the open forcefield_optimize NaN bug was localized: under the
# cuda-nvcxx backend that kernel returns all-NaN coordinates when nvc++
# compiles it at -O2 or above, while -O0/-O1 give correct geometry.  Building
# ONLY this file at -O1 with every other SYCL TU at -O3 makes
# eigen-functor-test pass, which pins the fault inside this one file rather
# than in the hessian/eigensolve kernels downstream of it.
if(DEFINED SYCL_TU_OPT_FILE OR DEFINED SYCL_TU_OPT_LEVEL)
  # Half an override silently builds everything at -O3 and reads as "the
  # hypothesis was wrong", which is the one answer a bisect must not invent.
  if(NOT (DEFINED SYCL_TU_OPT_FILE AND DEFINED SYCL_TU_OPT_LEVEL))
    message(FATAL_ERROR
      "SYCL_TU_OPT_FILE and SYCL_TU_OPT_LEVEL must be given together "
      "(got FILE='${SYCL_TU_OPT_FILE}' LEVEL='${SYCL_TU_OPT_LEVEL}').")
  endif()
  # APPEND, not set_source_files_properties: that REPLACES COMPILE_OPTIONS, and
  # the -fpermissive set above is on forcefield-optimize.cc -- the file this
  # knob exists to bisect -- so replacing it breaks the build under GNU instead
  # of answering the question.
  set_property(SOURCE ${SYCL_TU_OPT_FILE}
    APPEND PROPERTY COMPILE_OPTIONS ${SYCL_TU_OPT_LEVEL})
  message(STATUS "SYCL per-TU optimization override: ${SYCL_TU_OPT_FILE} at ${SYCL_TU_OPT_LEVEL}")
endif()

if(APPLE)
  target_link_libraries(sycl_fullerene_lib PRIVATE fullerenes)
endif()
