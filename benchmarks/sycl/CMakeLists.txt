set(SYCL_EXECUTABLES
  baseline_pipeline
  full_pipeline
  full_pipeline_strong_scaling  
  sequential_dualize
  lockstep_dualize
  matmul
  lobpcg-proto
  segmented_reduce
  general_dualize
)

foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  set_source_files_properties(${EXECUTABLE}.cc PROPERTIES COMPILE_FLAGS "${DPCPP_FLAGS}")    
endforeach()
  
foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.cc)  
  target_include_directories(${EXECUTABLE} PRIVATE ${DPCPP_INCLUDE} ${DPCPP_INCLUDE}/sycl)
  target_link_libraries(${EXECUTABLE} PRIVATE ${SYCL_FULLERENE_LIBS} sycl fullerenes)
endforeach()
