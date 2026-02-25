add_executable(gc-benchmark gc-benchmark.cc)
target_link_libraries(gc-benchmark fullerenes)

add_executable(gc-spirals gc-spirals.cc)
target_link_libraries(gc-spirals fullerenes)

add_executable(gc-spiral-benchmark gc-spiral-benchmark.cc)
target_link_libraries(gc-spiral-benchmark fullerenes)

add_executable(symmetry-benchmark symmetry-benchmark.cc)
target_link_libraries(symmetry-benchmark fullerenes)

link_libraries(fortran_opt fullerenes)
add_executable(fullerene-optimize-seq fullerene-optimize-seq.cc)
add_executable(generate_benchmark generate_benchmark.cc)

target_link_libraries(generate_benchmark fullerenes)
target_link_libraries(fullerene-optimize-seq fullerenes)

if(ENABLE_SYCL)
	target_link_libraries(generate_benchmark ${SYCL_FULLERENE_LIBS} ${SYCL_RUNTIME_LIBS})
	target_link_libraries(fullerene-optimize-seq ${SYCL_FULLERENE_LIBS} ${SYCL_RUNTIME_LIBS})
endif()
