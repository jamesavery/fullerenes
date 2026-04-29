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

add_executable(patch_diag patch_diag.cc)
target_link_libraries(patch_diag fullerenes)


add_executable(step_mol2 step_mol2.cc)
target_link_libraries(step_mol2 fullerenes)

add_executable(bench_epopt bench_epopt.cc)
target_link_libraries(bench_epopt fullerenes)

add_executable(bench_alexandrov bench_alexandrov.cc)
target_link_libraries(bench_alexandrov fullerenes)
if(OpenMP_CXX_FOUND)
  target_link_libraries(bench_alexandrov OpenMP::OpenMP_CXX)
endif()

add_executable(diag_selfloop diag_selfloop.cc)
target_link_libraries(diag_selfloop fullerenes)

add_executable(scan_topo scan_topo.cc)
target_link_libraries(scan_topo fullerenes)

add_executable(diag_mol2 diag_mol2.cc)
target_link_libraries(diag_mol2 fullerenes)

add_executable(scan_nonsimp scan_nonsimp.cc)
target_link_libraries(scan_nonsimp fullerenes)

add_executable(probe_multi probe_multi.cc)
target_link_libraries(probe_multi fullerenes)

add_executable(trace_palc trace_palc.cc)
target_link_libraries(trace_palc fullerenes)

add_executable(scan_diag_max scan_diag_max.cc)
target_link_libraries(scan_diag_max fullerenes)
if(OpenMP_CXX_FOUND)
  target_link_libraries(scan_diag_max OpenMP::OpenMP_CXX)
endif()

add_executable(spiral_mol2 spiral_mol2.cc)
target_link_libraries(spiral_mol2 fullerenes)

add_executable(iso_mol2 iso_mol2.cc)
target_link_libraries(iso_mol2 fullerenes buckygen)

add_executable(idt_lengths idt_lengths.cc)
target_link_libraries(idt_lengths fullerenes)

add_executable(perturb_solve perturb_solve.cc)
target_link_libraries(perturb_solve fullerenes)

add_executable(margin_sweep margin_sweep.cc)
target_link_libraries(margin_sweep fullerenes)

add_executable(replay_perstep replay_perstep.cc)
target_link_libraries(replay_perstep fullerenes)

add_executable(bench_quality bench_quality.cc)
target_link_libraries(bench_quality fullerenes)

add_executable(bench_quality_pipeline bench_quality_pipeline.cc)
target_link_libraries(bench_quality_pipeline fullerenes)
