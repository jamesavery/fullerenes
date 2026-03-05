
cmake_minimum_required(VERSION 3.5)
include(CheckLanguage)
include(CMakeDependentOption)
project(Fullerenes)
enable_language(CXX)
enable_language(Fortran)
set(CMAKE_GENERATOR "Unix Makefiles" CACHE INTERNAL "" FORCE)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(FULLERENE_DATABASE_PATH ${CMAKE_SOURCE_DIR}/database/ CACHE PATH "Path to fullerene database root")
set(FORTRAN_NMAX 5000 CACHE STRING "Maximum vertex count for Fortran routines using statically allocated memory")

set(BLA_VENDOR Intel10_64lp_seq)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
  unset(BLA_VENDOR)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()
find_package(GSL)
find_package(GTest REQUIRED)
find_library(TBB_LIB tbb)


enable_testing()
include_directories(${GTEST_INCLUDE_DIRS})

if(TBB_LIB)
  message(STATUS "Found TBB runtime: ${TBB_LIB}")
else()
  message(STATUS "Could not find TBB runtime library (tbb)")
endif()

# Shared between CUDA and SYCL
set(GPU_MAXREGCOUNT 80)

# Do we have CUDA compilers?
check_language(CUDA)
cmake_dependent_option(ENABLE_CUDA "Build CUDA GPU code" ON "CMAKE_CUDA_COMPILER" OFF)
if(CMAKE_CUDA_COMPILER)
  set(CMAKE_CUDA_ARCHITECTURES 86 CACHE INTEGER "Architecture for your NVIDIA GPU: 50|60|70|75|80|86|89|90. See https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/")
endif()

# Can we do SYCL with AdaptiveCpp compilers?
set(ADAPTIVECPP_ROOT "/opt/adaptivecpp" CACHE PATH "AdaptiveCpp installation prefix")

find_program(
  SYCL_COMPILER
  NAMES acpp syclcc
  HINTS ${ADAPTIVECPP_ROOT} ENV ADAPTIVECPP_ROOT ENV ACPP_ROOT
  PATH_SUFFIXES bin
)

if(NOT SYCL_COMPILER)
  # Fallback: search PATH in case wrappers are already exported.
  find_program(SYCL_COMPILER NAMES acpp syclcc)
endif()

if(SYCL_COMPILER)
  set(ENABLE_SYCL_DEFAULT ON)
else()
  set(ENABLE_SYCL_DEFAULT OFF)
endif()

option(ENABLE_SYCL "Build SYCL GPU code (requires AdaptiveCpp)" ${ENABLE_SYCL_DEFAULT})

if(SYCL_COMPILER)
  message(STATUS "Found AdaptiveCpp compiler: ${SYCL_COMPILER}")
else()
  message(STATUS "Could not find an AdaptiveCpp compiler (looked for: acpp, syclcc)")
endif()

if(ENABLE_SYCL)
  if(NOT SYCL_COMPILER)
    message(FATAL_ERROR "ENABLE_SYCL=ON requires AdaptiveCpp compiler (acpp or syclcc)")
  endif()

  if(CMAKE_CUDA_COMPILER)
    set(SYCL_DEFAULT_TARGET "NVIDIA")
  elseif(APPLE)
    set(SYCL_DEFAULT_TARGET "X86")
  else()
    set(SYCL_DEFAULT_TARGET "AMD")
  endif()
  
  set(SYCL_TARGETS ${SYCL_DEFAULT_TARGET} CACHE STRING "AMD|NVIDIA|X86|GENERIC (Target preset for AdaptiveCpp)")
  message("SYCL compile enabled with ${SYCL_TARGETS} target")

  set(SYCL_BASE_FLAGS -O3)
  set(SYCL_EXTRA_FLAGS "" CACHE STRING "Extra SYCL compiler flags")
  
  if(SYCL_TARGETS STREQUAL "AMD") 
    add_definitions(-DSOLVER_BACKEND=ROCSOLVER_BACKEND)
    set(SYCL_HIP_ARCH "gfx90a" CACHE STRING "HIP GPU architecture for AdaptiveCpp AMD builds")
    set(ACPP_TARGETS_DEFAULT "hip:${SYCL_HIP_ARCH}")
  elseif(SYCL_TARGETS STREQUAL "NVIDIA")
    add_definitions(-DSOLVER_BACKEND=CUSOLVER_BACKEND)
    set(SYCL_CUDA_ARCH "86" CACHE STRING "CUDA architecture for NVIDIA builds (e.g. 80, 86, 90)")
    set(SYCL_CUDA_BACKEND "cuda-nvcxx" CACHE STRING "AdaptiveCpp NVIDIA backend (cuda|cuda-nvcxx)")
    if(SYCL_CUDA_BACKEND STREQUAL "cuda")
      set(ACPP_TARGETS_DEFAULT "cuda:sm_${SYCL_CUDA_ARCH}")
    else()
      # nvc++ backend avoids clang CUDA header/toolchain incompatibilities.
      set(ACPP_TARGETS_DEFAULT "cuda-nvcxx:sm_${SYCL_CUDA_ARCH}")
    endif()
  elseif(SYCL_TARGETS STREQUAL "X86")
    # library-only backend sidesteps known accelerated CPU compiler crashes.
    set(ACPP_TARGETS_DEFAULT "omp.library-only")
  elseif(SYCL_TARGETS STREQUAL "GENERIC")
    set(ACPP_TARGETS_DEFAULT "generic")
  endif()

  # Keep backend target selection aligned with SYCL_TARGETS to avoid stale cache values.
  set(ACPP_TARGETS "${ACPP_TARGETS_DEFAULT}" CACHE STRING "AdaptiveCpp backend targets" FORCE)

  set(SYCL_FLAGS ${SYCL_BASE_FLAGS} --acpp-targets=${ACPP_TARGETS})
  if(SYCL_TARGETS STREQUAL "NVIDIA" AND SYCL_CUDA_BACKEND STREQUAL "cuda" AND CMAKE_CUDA_COMPILER)
    get_filename_component(CUDA_BIN_DIR ${CMAKE_CUDA_COMPILER} DIRECTORY)
    get_filename_component(CUDA_TOOLKIT_ROOT ${CUDA_BIN_DIR} DIRECTORY)
    list(APPEND SYCL_FLAGS --acpp-cuda-path=${CUDA_TOOLKIT_ROOT})
  endif()
  if(SYCL_EXTRA_FLAGS)
    separate_arguments(SYCL_EXTRA_FLAGS_LIST UNIX_COMMAND "${SYCL_EXTRA_FLAGS}")
    list(APPEND SYCL_FLAGS ${SYCL_EXTRA_FLAGS_LIST})
  endif()
  # AdaptiveCpp runtime: try acpp-rt first (brew/modern), fall back to sycl (DPC++)
  find_library(ACPP_RT_LIB NAMES acpp-rt hipSYCL-rt sycl
    HINTS ${ADAPTIVECPP_ROOT}/lib ENV ADAPTIVECPP_ROOT ENV ACPP_ROOT
  )
  if(ACPP_RT_LIB)
    set(SYCL_RUNTIME_LIBS ${ACPP_RT_LIB})
  else()
    set(SYCL_RUNTIME_LIBS sycl)
  endif()
  if(TBB_LIB)
    list(APPEND SYCL_RUNTIME_LIBS ${TBB_LIB})
  endif()
  if(SYCL_TARGETS STREQUAL "NVIDIA" AND SYCL_CUDA_BACKEND STREQUAL "cuda-nvcxx")
    find_library(OPENMP_RUNTIME_LIB NAMES omp)
    if(NOT OPENMP_RUNTIME_LIB)
      find_file(
        OPENMP_RUNTIME_LIB
        NAMES libomp.so.5
        PATHS /lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu
      )
    endif()
    if(OPENMP_RUNTIME_LIB)
      list(APPEND SYCL_RUNTIME_LIBS ${OPENMP_RUNTIME_LIB})
    endif()
  endif()

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
 set(CMAKE_Fortran_FLAGS "-fPIC -cpp -O3 -m64 -mcmodel=small")
endif()

# C++ compiler flags

set(CMAKE_CXX_FLAGS "-fPIC -O3 -m64 -pthread")

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
  geometry-properties_lib
  )

set(SYCL_FULLERENE_LIBS ${SYCL_FULLERENE_LIBS} CACHE STRING "SYCL_FULLERENE_LIBS")

if(ENABLE_SYCL)
  set(PREVIOUS_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  set(PREVIOUS_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  if(SYCL_TARGETS STREQUAL "NVIDIA" AND SYCL_CUDA_BACKEND STREQUAL "cuda-nvcxx")
    string(REPLACE "-fopenmp" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    string(REGEX REPLACE " +" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  endif()
  set(CMAKE_CXX_COMPILER ${SYCL_COMPILER})
  add_subdirectory(src/sycl)
  add_subdirectory(tests/sycl-tests)
  add_subdirectory(benchmarks/sycl)
  set(CMAKE_CXX_COMPILER ${PREVIOUS_CXX_COMPILER})
  set(CMAKE_CXX_FLAGS ${PREVIOUS_CXX_FLAGS})
endif()

message("C++ FLAGS: " ${CMAKE_CXX_FLAGS})
message("SYCL_FLAGS: " ${SYCL_FLAGS})
message("FLAGS: " ${FLAGS})


#set(SYCL_FLAGS "${SYCL_FLAGS} -tis")

#get_target_property(TESTVAR sycl_fullerene_lib LINK_FLAGS)
#message("sycl_fullerene_lib LINK_FLAGS: " ${TESTVAR})

add_subdirectory(programs)
add_subdirectory(benchmarks)
add_subdirectory(tests)




