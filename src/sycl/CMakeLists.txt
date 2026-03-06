list(APPEND SYCL_LIBRARY_SOURCES
  dualize
  forcefield-optimize
  tutte
  spherical-projection
  hessian
  eigen
  sycl-util-impl
  geometry-properties
)

list(APPEND SYCL_SOURCE_CUSTOM_FLAGS 
""    # dualize.cc
""    # forcefield-optimize.cc
""    # tutte.cc
""    # spherical-projection.cc
""    # hessian.cc
""    # eigen.cc
""    # sycl-util-impl.cc
""    # geometry-properties.cc
)


add_library(sycl_fullerene_lib SHARED nop_kernel.cc)
target_compile_options(sycl_fullerene_lib PRIVATE ${SYCL_FLAGS})
target_link_options(sycl_fullerene_lib PRIVATE ${SYCL_FLAGS})
target_link_libraries(sycl_fullerene_lib PRIVATE ${SYCL_RUNTIME_LIBS})
if(SYCL_TARGETS STREQUAL "NVIDIA" AND SYCL_CUDA_BACKEND STREQUAL "cuda-nvcxx")
  set_target_properties(sycl_fullerene_lib PROPERTIES NO_SONAME ON)
endif()

foreach(SYCL_SOURCE FLAGS IN ZIP_LISTS SYCL_LIBRARY_SOURCES SYCL_SOURCE_CUSTOM_FLAGS)
  add_library(${SYCL_SOURCE}_lib SHARED ${SYCL_SOURCE}.cc)
  target_compile_options(${SYCL_SOURCE}_lib PRIVATE ${SYCL_FLAGS} ${FLAGS})
  target_link_options(${SYCL_SOURCE}_lib PRIVATE ${SYCL_FLAGS})
  target_link_libraries(${SYCL_SOURCE}_lib PRIVATE ${SYCL_RUNTIME_LIBS})
  if(SYCL_TARGETS STREQUAL "NVIDIA" AND SYCL_CUDA_BACKEND STREQUAL "cuda-nvcxx")
    set_target_properties(${SYCL_SOURCE}_lib PROPERTIES NO_SONAME ON)
  endif()
  if(APPLE)
    target_link_libraries(${SYCL_SOURCE}_lib PRIVATE fullerenes)
  endif()
endforeach()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(forcefield-optimize_lib PRIVATE -fpermissive)
  target_compile_options(hessian_lib PRIVATE -fpermissive)
endif()

target_link_libraries(
  sycl_fullerene_lib 
  PRIVATE forcefield-optimize_lib 
  PRIVATE dualize_lib tutte_lib 
  PRIVATE spherical-projection_lib 
  PRIVATE hessian_lib eigen_lib 
  PRIVATE sycl-util-impl_lib 
  PRIVATE geometry-properties_lib 
)