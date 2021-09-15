enable_language(CXX CUDA)

add_library(cuda_lib SHARED kernel_concurrent.cu)

add_executable(benchmark_cuda_forcefield benchmark_cuda_forcefield.cu)
target_link_libraries(benchmark_cuda_forcefield cuda_lib)



