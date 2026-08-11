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
# CMake's auto-generated `-Wl,-soname,...` confuses pgacclnk (the nvc++
# CUDA linker wrapper) — see commit history. We disable it via NO_SONAME
# and pass an explicit basename soname ourselves so dependents embed
# `libsycl_fullerene_lib.so` in their NEEDED entry rather than the
# relative build-tree path.
set_target_properties(sycl_fullerene_lib PROPERTIES NO_SONAME 1)
target_link_options(sycl_fullerene_lib PRIVATE
  "LINKER:-soname,libsycl_fullerene_lib.so")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(forcefield-optimize.cc hessian.cc
    PROPERTIES COMPILE_OPTIONS "-fpermissive")
endif()

# Workaround: nvc++ (cuda-nvcxx backend) miscompiles forcefield-optimize.cc
# at -O2/-O3 with AdaptiveCpp 20.x / nvc++ 26.x, corrupting local-memory
# coord3d arrays during energy() reductions inside GSS line search and
# producing NaN geometries on the GPU. -O1 produces correct code.
# Remove this override once a compiler upgrade resolves the issue.
if(ACPP_TARGETS MATCHES "cuda-nvcxx")
  set_source_files_properties(forcefield-optimize.cc
    PROPERTIES COMPILE_FLAGS "-O1")
endif()

if(APPLE)
  target_link_libraries(sycl_fullerene_lib PRIVATE fullerenes)
endif()
