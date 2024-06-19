list(APPEND SYCL_LIBRARY_SOURCES
  forcefield-optimize
  dualize
  tutte
  spherical-projection
  sycl-isomer-batch
  hessian
  eigen
  sycl-util-impl
)

list(APPEND SYCL_SOURCE_CUSTOM_FLAGS 
"-Xcuda-ptxas --maxrregcount=80"    # forcefield-optimize.cc
"-Xcuda-ptxas --maxrregcount=40"    # dualize.cc
"-Xcuda-ptxas --maxrregcount=40"    # tutte.cc
"-Xcuda-ptxas --maxrregcount=40"    # spherical-projection.cc
"-Xcuda-ptxas --maxrregcount=40"    # sycl-isomer-batch.cc
"-Xcuda-ptxas --maxrregcount=255"     # hessian.cc
"-Xcuda-ptxas --maxrregcount=255"     # eigen.cc
"-Xcuda-ptxas --maxrregcount=40"    # sycl-util-impl.cc
)


add_library(sycl_fullerene_lib SHARED nop_kernel.cc)
set_target_properties(sycl_fullerene_lib PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${DPCPP_FLAGS}")
target_link_libraries(sycl_fullerene_lib ${DPCPP_FLAGS})

foreach(SYCL_SOURCE FLAGS IN ZIP_LISTS SYCL_LIBRARY_SOURCES SYCL_SOURCE_CUSTOM_FLAGS)
    add_library(${SYCL_SOURCE}_lib SHARED ${SYCL_SOURCE}.cc)
    set_target_properties(${SYCL_SOURCE}_lib PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${DPCPP_FLAGS} ${FLAGS}")
    target_link_libraries(${SYCL_SOURCE}_lib ${DPCPP_FLAGS} ${FLAGS})
endforeach()

target_link_libraries(sycl_fullerene_lib forcefield-optimize_lib dualize_lib tutte_lib spherical-projection_lib sycl-isomer-batch_lib hessian_lib eigen_lib sycl-util-impl_lib)

#add_library(sycl_fullerene_lib hessian.cc forcefield-optimize.cc dualize.cc tutte.cc spherical-projection.cc sycl-isomer-batch.cc nop_kernel.cc eigen.cc)

