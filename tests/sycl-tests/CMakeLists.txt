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
  eigen-functor-test
  to-polyhedron
  forcefield-functor-test
  sycl-vector-test
  mdspan-test
)

foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  set_source_files_properties(${EXECUTABLE}.cc PROPERTIES COMPILE_FLAGS "${DPCPP_FLAGS}")    
endforeach()
  
foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.cc)  
  target_include_directories(${EXECUTABLE} PRIVATE ${DPCPP_INCLUDE} ${DPCPP_INCLUDE}/sycl)
  target_link_libraries(${EXECUTABLE} PRIVATE ${SYCL_FULLERENE_LIBS} sycl fullerenes ${GTEST_LIBRARIES})
endforeach()

message(STATUS "SYCL_EXECUTABLES: ${SYCL_EXECUTABLES}")
message(STATUS "SYCL_FULLERENE_LIBS: ${SYCL_FULLERENE_LIBS}")

add_test(NAME fullerene-queue-test COMMAND fullerene-queue-test)
add_test(NAME condition-functor-test COMMAND condition-functor-test)
add_test(NAME tutte-functor-test COMMAND tutte-functor-test)
add_test(NAME eigen-functor-test COMMAND eigen-functor-test)
add_test(NAME mdspan-test COMMAND mdspan-test)
