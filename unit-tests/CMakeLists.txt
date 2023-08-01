if(ENABLE_CUDA)
  find_package(CUDAToolkit REQUIRED)
  link_libraries(fortran_opt fullerenes cuda_lib ${CMAKE_CUDA_LIBS})
  set(EXECUTABLES
  batch-clear-append-test
  async-execution-test
  validity-test
  race-test
  dual-test
  tutte-test
  tutte-test-random
  X0-test
  X0-test-random
  refill-race-test
  scan-test
  async-pipeline-test
  eigen-test
  produce-figure-data
  queue-test
  properties-test
  volume-test
  surface-area-test
  eccentricity-test
  minimum-problem
  eigensolve-special
)

foreach(EXECUTABLE ${EXECUTABLES})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.cc)
endforeach()

  if(CUSOLVER_LIBRARY)
    add_executable(hessian-test hessian-test.cc)
    target_link_libraries(hessian-test cuda_lib fullerenes  ${CMAKE_CUDA_LIBS})
  endif()  

  foreach(EXECUTABLE ${EXECUTABLES})
    target_link_libraries(${EXECUTABLE} cuda_lib fullerenes  ${CMAKE_CUDA_LIBS})  
  endforeach()
  
endif()
