
cmake_minimum_required(VERSION 3.22)
project(Fullerenes)

set(CMAKE_CXX_STANDARD 17)
set(FULLERENE_DATABASE_PATH ${CMAKE_SOURCE_DIR}/database/ CACHE PATH "Path to fullerene database root")
set(FORTRAN_NMAX 5000 CACHE STRING "Maximum vertex count for Fortran routines using statically allocated memory")

set(BLA_VENDOR Intel10_64lp_seq)
find_package(BLAS)
find_package(LAPACK)
find_package(GSL)

find_program(ICPX_COMPILER icpx)

if(ICPX_COMPILER)
  message("Found icpx compiler: ${ICPX_COMPILER}")
else()
  message("Could not find icpx compiler")
endif()


include(CheckLanguage)
include(CMakeDependentOption)
check_language(CUDA)
cmake_dependent_option(ENABLE_CUDA "Build CUDA GPU code" ON "CMAKE_CUDA_COMPILER" OFF)

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

include_directories(AFTER SYSTEM include)
include_directories(${CMAKE_BINARY_DIR}/include/)


# Fortran compiler flags
set(CMAKE_Fortran_FLAGS "-fPIC -cpp -O3 -m64 -mcmodel=medium")
# C++ compiler flags


set(CMAKE_CXX_FLAGS "-std=c++1z -fPIC -O3 -m64 -pthread -fopenmp")



# TODO: Add OpenMP Isomerspace implementation, choose at runtime if both openmp and cuda present. OpenACC too?


if(ENABLE_CUDA)
  enable_language(CUDA)
  set(CUDA_INCLUDE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  
  message("CUDA include path at:" ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
  message("CUDA found: " ${CMAKE_CUDA_COMPILER} ".")
#set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} "--restrict --fmad true --ftz true  --use_fast_math --maxrregcount=64 --resource-usage -Xptxas --expt-relaxed-constexpr -I /opt/nvidia/hpc_sdk/Linux_x86_64/21.2/cuda/11.0/targets/x86_64-linux/include/ -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}/../../../../../examples/OpenACC/SDK/include/"
#NOTE: -gencode arch argument is EXTREMELY important to get right, the program will simply crash with a cudaError: invalid device function if the architecture does not match the GPU

set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} "--std c++17 --use_fast_math --expt-relaxed-constexpr -lineinfo -Xptxas -O3,-v --resource-usage -gencode arch=compute_86,code=sm_86  -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"

)
set(CUDA_MAXREGCOUNT "--maxrregcount=80")

  # Optional extra CUDA libraries
  find_library(CUSOLVER_LIBRARY cusolver HINTS ${CUDAToolkit_LIBRARaY_DIR})
  find_library(CUSPARSE_LIBRARY cusparse HINTS ${CUDAToolkit_LIBRARY_DIR})
  if( CUSOLVER_LIBRARY )
    set(CMAKE_CUDA_LIBS  ${CMAKE_CUDA_LIBS} ${CUSOLVER_LIBRARY})
    set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} "-DCUSOLVER=1")
    message("CUSOLVER library found: " ${CUSOLVER_LIBRARY} ".")
  endif()

  if( CUSPARSE_LIBRARY )
    set(CMAKE_CUDA_LIBS  ${CMAKE_CUDA_LIBS} ${CUSPARSE_LIBRARY})
    set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} "-DCUSPARSE=1")
    message("CUSPARSE library found: " ${CUSPARSE_LIBRARY} ".")
  endif()

  add_subdirectory(src/cuda)
endif()

configure_file(
  ${CMAKE_SOURCE_DIR}/include/fullerenes/config.hh.in
  ${CMAKE_BINARY_DIR}/include/fullerenes/config.hh
  )

# Fortran configure file
configure_file(
  ${CMAKE_SOURCE_DIR}/include/fullerenes/config.f.in
  ${CMAKE_SOURCE_DIR}/src/fortran/config.f
)


add_subdirectory(src/fortran)
add_subdirectory(src/contrib)
add_subdirectory(playground)
add_subdirectory(src/c++)

if(ICPX_COMPILER)
  set(PREVIOUS_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER ${ICPX_COMPILER})
  add_subdirectory(src/sycl)
  set(CMAKE_CXX_COMPILER ${PREVIOUS_CXX_COMPILER})
endif()

add_subdirectory(programs)
add_subdirectory(benchmarks)
add_subdirectory(unit-tests)




