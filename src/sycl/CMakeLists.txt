set_source_files_properties(forcefield-optimise.cc PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fsycl -fsycl-targets=nvptx64-nvidia-cuda")
set_source_files_properties(dualise.cc PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fsycl -fsycl-targets=nvptx64-nvidia-cuda")
set_source_files_properties(tutte.cc PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fsycl -fsycl-targets=nvptx64-nvidia-cuda")
add_library(sycl_fullerene_lib SHARED forcefield-optimise.cc dualise.cc tutte.cc)