set(SYCL_EXECUTABLES
  baseline_pipeline
  full_pipeline
  full_pipeline_strong_scaling  
  sequential_dualize
  lockstep_dualize
  general_dualize
)

foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.cc)  
  target_compile_options(${EXECUTABLE} PRIVATE ${SYCL_FLAGS})
  target_link_options(${EXECUTABLE} PRIVATE ${SYCL_FLAGS})
  target_link_libraries(${EXECUTABLE} PRIVATE ${SYCL_FULLERENE_LIBS} ${SYCL_RUNTIME_LIBS} fullerenes)
endforeach()