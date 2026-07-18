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
  sycl-smoke-test
  mdspan-test
  pipeline-view-test
)

if(GTest_FOUND)
  foreach(EXECUTABLE ${SYCL_EXECUTABLES})
    add_executable(${EXECUTABLE} ${EXECUTABLE}.cc)
    target_compile_options(${EXECUTABLE} PRIVATE ${SYCL_FLAGS})
    target_link_options(${EXECUTABLE} PRIVATE ${SYCL_FLAGS})
    target_link_libraries(${EXECUTABLE} PRIVATE ${SYCL_FULLERENE_LIBS} ${SYCL_RUNTIME_LIBS} fullerenes GTest::GTest GTest::Main)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${EXECUTABLE} PRIVATE -fpermissive)
    endif()
  endforeach()

  message(STATUS "SYCL_EXECUTABLES: ${SYCL_EXECUTABLES}")
  message(STATUS "SYCL_FULLERENE_LIBS: ${SYCL_FULLERENE_LIBS}")

  add_test(NAME fullerene-queue-test COMMAND fullerene-queue-test)
  add_test(NAME condition-functor-test COMMAND condition-functor-test)
  add_test(NAME tutte-functor-test COMMAND tutte-functor-test)
  add_test(NAME eigen-functor-test COMMAND eigen-functor-test)
  add_test(NAME sycl-smoke-test COMMAND sycl-smoke-test)
  add_test(NAME mdspan-test COMMAND mdspan-test)
else()
  message(STATUS "GTest not found -- skipping GTest-based tests in tests/sycl-tests/")
endif()
