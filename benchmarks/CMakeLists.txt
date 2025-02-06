
link_libraries(fortran_opt fullerenes ${SYCL_FULLERENE_LIBS})
add_executable(fullerene-optimize-seq fullerene-optimize-seq.cc)
add_executable(generate_benchmark generate_benchmark.cc)

target_link_libraries(generate_benchmark fullerenes ${SYCL_FULLERENE_LIBS})
target_link_libraries(fullerene-optimize-seq  fullerenes ${SYCL_FULLERENE_LIBS})
