link_libraries(fortran_opt fullerenes)

if(CMAKE_CUDA_COMPILER)
  add_executable(tutte-kernel-test tutte-kernel-test.cc)
  add_executable(X0-kernel-test X0-kernel-test.cc)
  target_link_libraries(tutte-kernel-test cuda_lib)
  target_link_libraries(X0-kernel-test cuda_lib)
endif()
