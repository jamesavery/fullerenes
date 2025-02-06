list(APPEND SYCL_LIBRARY_SOURCES
  dualize
  forcefield-optimize
  tutte
  spherical-projection
  hessian
  eigen
  sycl-util-impl
  primitives-instantiations
  geometry-properties
  lobpcg-v1
  lobpcg-v2
  lobpcg-v3
  lobpcg-v4
  lobpcg-v5
  lobpcg-v6
  lobpcg-v7
  lobpcg-v8
)

list(APPEND SYCL_SOURCE_CUSTOM_FLAGS 
"-Xcuda-ptxas --maxrregcount=40"    # dualize.cc
"-Xcuda-ptxas --maxrregcount=80"    # forcefield-optimize.cc
"-Xcuda-ptxas --maxrregcount=40"    # tutte.cc
"-Xcuda-ptxas --maxrregcount=40"    # spherical-projection.cc
"-Xcuda-ptxas --maxrregcount=80"     # hessian.cc
"-Xcuda-ptxas --maxrregcount=80"     # eigen.cc
"-Xcuda-ptxas --maxrregcount=40 -ltbb"   # sycl-util-impl.cc
"-Xcuda-ptxas --maxrregcount=40 -ltbb"    # primitives.cc
"-Xcuda-ptxas --maxrregcount=40"    # geometry-properties.cc
"-Xcuda-ptxas --maxrregcount=80"    # lobpcg-v1.cc
"-Xcuda-ptxas --maxrregcount=80"    # lobpcg-v2.cc
"-Xcuda-ptxas --maxrregcount=80"    # lobpcg-v3.cc
"-Xcuda-ptxas --maxrregcount=80"    # lobpcg-v4.cc
"-Xcuda-ptxas --maxrregcount=80"    # lobpcg-v5.cc
"-Xcuda-ptxas --maxrregcount=80"    # lobpcg-v6.cc
"-Xcuda-ptxas --maxrregcount=80"    # lobpcg-v7.cc
"-Xcuda-ptxas --maxrregcount=80"    # lobpcg-v8.cc
)


add_library(sycl_fullerene_lib SHARED nop_kernel.cc)
set_target_properties(sycl_fullerene_lib PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${DPCPP_FLAGS}")
target_link_libraries(sycl_fullerene_lib PRIVATE ${CMAKE_CXX_FLAGS} ${DPCPP_FLAGS})

foreach(SYCL_SOURCE FLAGS IN ZIP_LISTS SYCL_LIBRARY_SOURCES SYCL_SOURCE_CUSTOM_FLAGS)
    add_library(${SYCL_SOURCE}_lib SHARED ${SYCL_SOURCE}.cc)
    set_target_properties(${SYCL_SOURCE}_lib PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${DPCPP_FLAGS} ${FLAGS}")
    target_link_libraries(${SYCL_SOURCE}_lib PRIVATE ${CMAKE_CXX_FLAGS} ${DPCPP_FLAGS} ${FLAGS} cudart cuda cusolver cublas cusparse)
endforeach()

target_link_libraries(
  sycl_fullerene_lib 
  PRIVATE forcefield-optimize_lib 
  PRIVATE dualize_lib tutte_lib 
  PRIVATE spherical-projection_lib 
  PRIVATE hessian_lib eigen_lib 
  PRIVATE sycl-util-impl_lib 
  PRIVATE primitives-instantiations_lib 
  PRIVATE geometry-properties_lib 
  PRIVATE lobpcg-v1_lib 
  PRIVATE lobpcg-v2_lib
  PRIVATE lobpcg-v3_lib
  PRIVATE lobpcg-v4_lib
  PRIVATE lobpcg-v5_lib
  PRIVATE lobpcg-v6_lib
  PRIVATE lobpcg-v7_lib
  PRIVATE lobpcg-v8_lib
)

add_subdirectory(linalg/batched/QR)