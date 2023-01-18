link_libraries(fortran_opt fullerenes)

if(CMAKE_CUDA_COMPILER)
  add_executable(queue-test queue-test.cc)
  target_link_libraries(queue-test cuda_lib fullerenes)
endif()
