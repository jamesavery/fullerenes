
add_library(nvidia_blas SHARED level_3.cc)
set_target_properties(nvidia_blas PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${DPCPP_FLAGS} -I/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/math_libs/12.6/targets/x86_64-linux/include")
target_link_libraries(nvidia_blas PRIVATE ${CMAKE_CXX_FLAGS} ${DPCPP_FLAGS} cublas cusparse cusolver cudart cuda)