link_libraries(fortran_opt fullerenes)
add_executable(fullerene-optimize-seq fullerene-optimize-seq.cc)
add_executable(isomerspace-optimize-par isomerspace-optimize-par.cc)
target_link_libraries(isomerspace-optimize-par cuda_lib)