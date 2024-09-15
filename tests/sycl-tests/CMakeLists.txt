set(SYCL_EXECUTABLES
  #hello-sycl
  #sycl-pipeline-test
  #sycl-pipeline-herd
  #sycl-pipeline-multinode
  #eigen-test-sycl
  #dual-test-sycl
  #dual-validation-sycl
  fullerene-queue-test
  condition-functor-test
  tutte-functor-test
  to-polyhedron
)

foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  set_source_files_properties(${EXECUTABLE}.cc PROPERTIES COMPILE_FLAGS "${DPCPP_FLAGS}")    
endforeach()
  
foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.cc)  
  target_include_directories(${EXECUTABLE} PRIVATE ${DPCPP_INCLUDE} ${DPCPP_INCLUDE}/sycl)
  target_link_libraries(${EXECUTABLE} PRIVATE sycl_fullerene_lib dualize_lib sycl fullerenes ${GTEST_LIBRARIES})
endforeach()



add_test(NAME fullerene-queue-test COMMAND fullerene-queue-test)
add_test(NAME condition-functor-test COMMAND condition-functor-test)
add_test(NAME tutte-functor-test COMMAND tutte-functor-test)