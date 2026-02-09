
cmake_minimum_required(VERSION 3.5)
include(CheckLanguage)
include(CMakeDependentOption)
project(Fullerenes)
enable_language(CXX)
enable_language(Fortran)
set(CMAKE_GENERATOR "Unix Makefiles" CACHE INTERNAL "" FORCE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(FULLERENE_DATABASE_PATH ${CMAKE_SOURCE_DIR}/database/ CACHE PATH "Path to fullerene database root")
set(FORTRAN_NMAX 5000 CACHE STRING "Maximum vertex count for Fortran routines using statically allocated memory")

set(BLA_VENDOR Intel10_64lp_seq)
find_package(BLAS)
find_package(LAPACK)
find_package(GSL)
find_package(GTest REQUIRED)


enable_testing()
include_directories(${GTEST_INCLUDE_DIRS})

# Shared between CUDA and SYCL
set(GPU_MAXREGCOUNT 80)

# Do we have CUDA compilers?
check_language(CUDA)
cmake_dependent_option(ENABLE_CUDA "Build CUDA GPU code" ON "CMAKE_CUDA_COMPILER" OFF)
if(CMAKE_CUDA_COMPILER)
  set(CMAKE_CUDA_ARCHITECTURES 86 CACHE INTEGER "Architecture for your NVIDIA GPU: 50|60|70|75|80|86|89|90. See https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/")
endif()

# Can we do a SYCL?
find_program(SYCL_LS sycl-ls)
if(SYCL_LS)
  get_filename_component(SYCL_BIN_DIR ${SYCL_LS} DIRECTORY)
  set(SYCL_DIR ${SYCL_BIN_DIR}/../)
  find_program(DPCPP_COMPILER NAMES ${SYCL_BIN_DIR}/icpx ${SYCL_BIN_DIR}/clang++)

  cmake_dependent_option(ENABLE_SYCL "Build SYCL GPU code" ON "DPCPP_COMPILER" OFF)  
  
  if(DPCPP_COMPILER)
    message("Found dpc++ compiler: ${DPCPP_COMPILER}")
  else()
    message("Could not find dpc++ compiler")
  endif()
  
endif()

if(ENABLE_SYCL)
  if(CMAKE_CUDA_COMPILER)
    set(SYCL_DEFAULT_TARGET "NVIDIA")
  else()
    set(SYCL_DEFAULT_TARGET "AMD")
  endif()
  
  set(SYCL_TARGETS ${SYCL_DEFAULT_TARGET} CACHE STRING "AMD|NVIDIA|X86 (Currently only single target supported.") # TODO: Multi-target compiles.
  message("SYCL compile enabled with ${SYCL_TARGETS} target")

  set(DPCPP_FLAGS "-Wno-pass-failed -Wunused-variable -Wno-unused-command-line-argument -Wno-deprecated-declarations -Wno-recommended-option -O3 -fsycl -fsycl-unnamed-lambda") # Use regardless of target.
  
  if(SYCL_TARGETS STREQUAL "AMD") 
    add_definitions(-DSOLVER_BACKEND=ROCSOLVER_BACKEND)
    set(DPCPP_TARGET_FLAGS "-fsycl-targets=amdgcn-amd-amdhsa,x86_64 -Xsycl-target-backend=amdgcn-amd-amdhsa --offload-arch=gfx90a -march=native")
  elseif(SYCL_TARGETS STREQUAL "NVIDIA")
    add_definitions(-DSOLVER_BACKEND=CUSOLVER_BACKEND)
    #set(PTXAS_COMPILE_FLAGS    " --verbose -g --maxrregcount=${GPU_MAXREGCOUNT} --opt-level=3 --generate-line-info --warn-on-double-precision-use --warn-on-spills")
    set(PTXAS_LINK_FLAGS
    "SHELL:-Xcuda-ptxas --verbose"
    "SHELL:-Xcuda-ptxas --maxrregcount=255"
    "SHELL:-Xcuda-ptxas --opt-level=3"
    "SHELL:-Xcuda-ptxas --generate-line-info"
    "SHELL:-Xcuda-ptxas --warn-on-double-precision-use"
    "SHELL:-Xcuda-ptxas --warn-on-spills"
    )
    set(DPCPP_TARGET_FLAGS "-Xcuda-ptxas --generate-line-info -fsycl-targets=nvptx64-nvidia-cuda -Xsycl-target-backend=nvptx64-nvidia-cuda --cuda-gpu-arch=sm_${CMAKE_CUDA_ARCHITECTURES} --cuda-path=$ENV{CUDA_PATH}")
    set(DPCPP_LINK_FLAGS   "-O3 -fsycl -fsycl-unnamed-lambda ${DPCPP_TARGET_FLAGS}")    
  elseif(SYCL_TARGETS STREQUAL "X86")
    set(DPCPP_TARGET_FLAGS "-fsycl-targets=x86_64 -march=native -mcpu=native")
  endif()

  set(DPCPP_INCLUDE ${SYCL_DIR}include)
  message("Using SYCL include path: ${DPCPP_INCLUDE}")  
  set(DPCPP_FLAGS "${DPCPP_FLAGS} -fsycl -fsycl-unnamed-lambda ${DPCPP_TARGET_FLAGS} -I${DPCPP_INCLUDE} -I${DPCPP_INCLUDE}/sycl -I$ENV{CUDA_INCLUDE}")

endif()


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
if(DEFINED ENV{CRAYPE_VERSION})
  message("Cray compiler environment detected")
  set(CMAKE_Fortran_FLAGS "-fPIC -O3")
else()
 set(CMAKE_Fortran_FLAGS "-fPIC -cpp -O3 -m64 -mcmodel=medium")
endif()

# C++ compiler flags

set(CMAKE_CXX_FLAGS "-fPIC -O3 -m64 -pthread -fopenmp")



if(ENABLE_CUDA)
  enable_language(CUDA)
  find_package(CUDAToolkit)
  set(CUDA_INCLUDE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  
  # Include the directory containing curand.h
  include_directories(${CURAND_INCLUDE_DIR})
  message("CUDA include path at:" ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  message("")
  message("CUDA libraries at:" ${CUDAToolkit_LIBRARY_DIR})
  message("CUDA Binaries at:" ${CUDAToolkit_BIN_DIR})
  message("CUDA include path at:" ${CUDAToolkit_INCLUDE_DIRS})
  include_directories(${CUDAToolkit_INCLUDE_DIRS})

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
  message("CUDA found: " ${CMAKE_CUDA_COMPILER} ".")
#NOTE: -gencode arch argument is EXTREMELY important to get right, the program will simply crash with a cudaError: invalid device function if the architecture does not match the GPU

set(CUDA_MAXREGCOUNT "--maxrregcount=80")

  message("Include Path :))) " $ENV{CUDA_INCLUDE})
  # Optional extra CUDA libraries
  find_library(CUSOLVER_LIBRARY cusolver HINTS "/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/math_libs/12.6/targets/x86_64-linux/lib")
  find_library(CUSPARSE_LIBRARY cusparse HINTS "/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/math_libs/12.6/targets/x86_64-linux/lib")
  find_library(CUBLAS_LIBRARY cublas HINTS "/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/math_libs/12.6/targets/x86_64-linux/lib")
  if( CUSOLVER_LIBRARY )
    set(CMAKE_CUDA_LIBS  ${CMAKE_CUDA_LIBS} ${CUSOLVER_LIBRARY})
    message("CUSOLVER library found: " ${CUSOLVER_LIBRARY} ".")
  endif()

  if( CUSPARSE_LIBRARY )
    set(CMAKE_CUDA_LIBS  ${CMAKE_CUDA_LIBS} ${CUSPARSE_LIBRARY})
    message("CUSPARSE library found: " ${CUSPARSE_LIBRARY} ".")
  endif()

  if( CUBLAS_LIBRARY )
    set(CMAKE_CUDA_LIBS  ${CMAKE_CUDA_LIBS} ${CUBLAS_LIBRARY})
    message("CUBLAS library found: " ${CUBLAS_LIBRARY} ".")
  endif()
  
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

set(SYCL_FULLERENE_LIBS 
  sycl_fullerene_lib 
  forcefield-optimize_lib 
  dualize_lib 
  tutte_lib 
  spherical-projection_lib 
  hessian_lib 
  eigen_lib 
  sycl-util-impl_lib 
  primitives-instantiations_lib 
  geometry-properties_lib
  lobpcg-v1_lib
  lobpcg-v2_lib
  lobpcg-v3_lib
  lobpcg-v4_lib
  lobpcg-v5_lib
  lobpcg-v6_lib
  lobpcg-v7_lib
  lobpcg-v8_lib
  )

set(SYCL_FULLERENE_LIBS ${SYCL_FULLERENE_LIBS} CACHE STRING "SYCL_FULLERENE_LIBS")

if(ENABLE_SYCL)
  set(PREVIOUS_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER ${DPCPP_COMPILER})
  add_subdirectory(src/sycl)
  add_subdirectory(src/linalg/wrappers)
  add_subdirectory(tests/sycl-tests)
  add_subdirectory(benchmarks/sycl)
  set(CMAKE_CXX_COMPILER ${PREVIOUS_CXX_COMPILER})
endif()

message("C++ FLAGS: " ${CMAKE_CXX_FLAGS})
message("DPCPP_FLAGS: " ${DPCPP_FLAGS})
message("FLAGS: " ${FLAGS})


#set(DPCPP_FLAGS "${DPCPP_FLAGS} -tis")

#get_target_property(TESTVAR sycl_fullerene_lib LINK_FLAGS)
#message("sycl_fullerene_lib LINK_FLAGS: " ${TESTVAR})

add_subdirectory(programs)
add_subdirectory(benchmarks)
add_subdirectory(tests)




