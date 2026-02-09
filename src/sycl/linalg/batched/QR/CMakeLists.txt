list(APPEND SYCL_BATCHED_QR_SOURCES
  chol_qr
)

foreach(SYCL_BATCHED_QR_SOURCE IN LISTS SYCL_BATCHED_QR_SOURCES)
    add_library(${SYCL_BATCHED_QR_SOURCE}_lib SHARED ${SYCL_BATCHED_QR_SOURCE}.cc)
    set_target_properties(${SYCL_BATCHED_QR_SOURCE}_lib PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${DPCPP_FLAGS} -I/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/math_libs/12.6/targets/x86_64-linux/include")
    target_link_libraries(${SYCL_BATCHED_QR_SOURCE}_lib PRIVATE ${CMAKE_CXX_FLAGS} ${DPCPP_FLAGS} nvidia_blas cudart cuda)
endforeach()
