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
  # SYCL_TARGETS_FLAG stays quoted: it may carry an escaped ';' separating
  # AdaptiveCpp backends, which an unquoted expansion would split.  Without
  # it these binaries get no target at all and fall back to acpp's configured
  # default, i.e. a different backend than the library they link.
  target_compile_options(${EXECUTABLE} PRIVATE ${SYCL_FLAGS} "${SYCL_TARGETS_FLAG}")
  target_link_options(${EXECUTABLE} PRIVATE ${SYCL_FLAGS} "${SYCL_TARGETS_FLAG}")
  target_link_libraries(${EXECUTABLE} PRIVATE ${SYCL_FULLERENE_LIBS} ${SYCL_RUNTIME_LIBS} fullerenes)
endforeach()