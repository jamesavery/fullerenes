
cmake_minimum_required(VERSION 3.1)
project(Fullerenes)

set(CMAKE_CXX_STANDARD 17)
set(FULLERENE_DATABASE_PATH ${CMAKE_SOURCE_DIR}/database/ CACHE PATH "Path to fullerene database root")
set(FORTRAN_NMAX 5000 CACHE STRING "Maximum vertex count for Fortran routines using statically allocated memory")

set(BLA_VENDOR Intel10_64lp_seq)
find_package(BLAS)
find_package(LAPACK)
find_package(GSL)
include(CheckLanguage)
check_language(CUDA)

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# C++ configure file
configure_file(
  ${CMAKE_SOURCE_DIR}/include/fullerenes/config.hh.in
  ${CMAKE_BINARY_DIR}/include/fullerenes/config.hh
  )

# Fortran configure file
configure_file(
  ${CMAKE_SOURCE_DIR}/include/fullerenes/config.f.in
  ${CMAKE_SOURCE_DIR}/src/fortran/config.f
)


include_directories(AFTER SYSTEM include)
include_directories(${CMAKE_BINARY_DIR}/include/)


# Fortran compiler flags
set(CMAKE_Fortran_FLAGS "-fPIC -cpp -O2 -m64 -mcmodel=medium")
# C++ compiler flags
set(CMAKE_CXX_FLAGS "-std=c++1z -fPIC -O2 -m64 -pthread")

add_subdirectory(src/c++)
add_subdirectory(src/fortran)
add_subdirectory(src/contrib)
add_subdirectory(playground)
add_subdirectory(programs)
add_subdirectory(benchmarks)
add_subdirectory(tests)
add_subdirectory(unit-tests)

# TODO: Add OpenMP Isomerspace implementation, choose at runtime if both openmp and cuda present. OpenACC too?


if(CMAKE_CUDA_COMPILER)
  message("CUDA found: " ${CMAKE_CUDA_COMPILER} ".")
#set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} "--restrict --fmad true --ftz true  --use_fast_math --maxrregcount=64 --resource-usage -Xptxas --expt-relaxed-constexpr -I /opt/nvidia/hpc_sdk/Linux_x86_64/21.2/cuda/11.0/targets/x86_64-linux/include/ -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}/../../../../../examples/OpenACC/SDK/include/"
#NOTE: -gencode arch argument is EXTREMELY important to get right, the program will simply crash with a cudaError: invalid device function if the architecture does not match the GPU
set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} "--std c++17 --use_fast_math --fmad true --ftz true --maxrregcount=64 -lineinfo -Xptxas -v --resource-usage -gencode arch=compute_86,code=sm_86  -I /opt/nvidia/hpc_sdk/Linux_x86_64/22.2/cuda/11.6/targets/x86_64-linux/include/ -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}/../../../../../examples/OpenACC/SDK/include/"
)
  
  add_subdirectory(src/cuda)
endif()





