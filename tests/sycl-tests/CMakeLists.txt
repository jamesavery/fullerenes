set(SYCL_EXECUTABLES
  hello-sycl
  sycl-pipeline-test
  sycl-pipeline-herd
  sycl-pipeline-multinode
  eigen-test-sycl
)

foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  set_source_files_properties(${EXECUTABLE}.cc PROPERTIES COMPILE_FLAGS "${DPCPP_FLAGS} ${DPCPP_TARGET_FLAGS}")    
endforeach()
  
foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.cc)  
  target_include_directories(${EXECUTABLE} PRIVATE ${DPCPP_INCLUDE} ${DPCPP_INCLUDE}/sycl)
  target_link_libraries(${EXECUTABLE} PRIVATE sycl_fullerene_lib sycl fullerenes m)
  target_link_options(${EXECUTABLE} PRIVATE     
    "SHELL:${DPCPP_LINK_FLAGS}"
    "SHELL:${DPCPP_FLAGS}"
    ${PTXAS_LINK_FLAGS}
  )
endforeach()

