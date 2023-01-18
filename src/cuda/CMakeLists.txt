enable_language(CXX CUDA)


add_library(cuda_lib SHARED isomer_queue.cu cu_array.cu misc_cuda.cu forcefield.cu tutte.cu spherical_projection.cu dualize.cu benchmark_functions.cu)
set_property(TARGET cuda_lib PROPERTY CUDA_ARCHITECTURES "86")
#set_target_properties(cuda_lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#set_target_properties(cuda_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
#set_target_properties(cuda_lib PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

