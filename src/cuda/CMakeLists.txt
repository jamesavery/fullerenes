enable_language(CXX CUDA)

set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} "--restrict --fmad true --ftz true  --use_fast_math --maxrregcount=128 --expt-relaxed-constexpr -O3 -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}/../../../../../examples/OpenACC/SDK/include/"
)

add_library(cuda_lib SHARED kernel_shared.cu)

add_executable(benchmark_cuda_forcefield benchmark_cuda_forcefield.cu)
target_link_libraries(benchmark_cuda_forcefield cuda_lib)



