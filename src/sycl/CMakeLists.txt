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

set(SYCL_ENABLE_SONAME TRUE)
if(CMAKE_CXX_COMPILER MATCHES "acpp" AND ACPP_TARGETS MATCHES "cuda-nvcxx")
  set(SYCL_ENABLE_SONAME FALSE)
endif()

function(configure_sycl_shared_library target_name)
  target_compile_options(${target_name} PRIVATE ${SYCL_FLAGS} ${ARGN})
  target_link_options(${target_name} PRIVATE ${SYCL_FLAGS})
  target_link_libraries(${target_name} PRIVATE ${SYCL_RUNTIME_LIBS})
  if(SYCL_ENABLE_SONAME)
    set_target_properties(${target_name} PROPERTIES VERSION 1.0.0 SOVERSION 1)
  endif()
  if(CMAKE_CXX_COMPILER MATCHES "acpp" AND ACPP_TARGETS MATCHES "cuda-nvcxx")
    add_custom_command(TARGET ${target_name} PRE_LINK
      COMMAND ${CMAKE_COMMAND} -E rm -f $<TARGET_FILE:${target_name}>
    )
  endif()
endfunction()


add_library(sycl_fullerene_lib SHARED nop_kernel.cc)
configure_sycl_shared_library(sycl_fullerene_lib)

foreach(SYCL_SOURCE FLAGS IN ZIP_LISTS SYCL_LIBRARY_SOURCES SYCL_SOURCE_CUSTOM_FLAGS)
  add_library(${SYCL_SOURCE}_lib SHARED ${SYCL_SOURCE}.cc)
  configure_sycl_shared_library(${SYCL_SOURCE}_lib ${FLAGS})
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