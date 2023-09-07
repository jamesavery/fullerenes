
include_directories(/opt/intel/oneapi/compiler/2021.1.1/linux/include/sycl)
set_source_files_properties(forcefield-optimise.cpp PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fsycl -fsycl-targets=nvptx64-nvidia-cuda")
set_source_files_properties(dualise.cpp PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fsycl -fsycl-targets=nvptx64-nvidia-cuda")
add_library(sycl_fullerene_lib SHARED forcefield-optimise.cpp dualise.cpp)