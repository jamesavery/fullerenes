
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

  set(DPCPP_FLAGS "-Wno-pass-failed -Wunused-variable -Wno-unused-command-line-argument -Wno-deprecated-declarations -Wno-experimental-option -Wno-c++20-extensions -O3 -fsycl -fsycl-unnamed-lambda") # Use regardless of target.
  
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
    set(DPCPP_TARGET_FLAGS "--cuda-path=$ENV{CUDA_PATH}")
    set(DPCPP_LINK_FLAGS   "-O3 -fsycl -fsycl-unnamed-lambda ${DPCPP_TARGET_FLAGS}")    
  elseif(SYCL_TARGETS STREQUAL "X86")
    set(DPCPP_TARGET_FLAGS "-fsycl-targets=x86_64 -march=native -mcpu=native")
  endif()

  set(DPCPP_INCLUDE ${SYCL_DIR}include)
  message("Using SYCL include path: ${DPCPP_INCLUDE}")  
  set(DPCPP_FLAGS "${DPCPP_FLAGS} ${DPCPP_TARGET_FLAGS} -I${DPCPP_INCLUDE} -I${DPCPP_INCLUDE}/sycl -I$ENV{CUDA_INCLUDE}")

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

# oneDPL (DPC++ STL) headers are not installed in /opt/dpcpp-cuda by default.
# If a oneAPI oneDPL install exists, add it explicitly.
set(ONEDPL_ROOT "/opt/intel/oneapi/dpl/latest")
if(EXISTS "${ONEDPL_ROOT}/include/oneapi/dpl")
    message(STATUS "Adding oneDPL include dir: ${ONEDPL_ROOT}/include")
    include_directories("${ONEDPL_ROOT}/include")
endif()


# Fortran compiler flags
if(DEFINED ENV{CRAYPE_VERSION})
  message("Cray compiler environment detected")
  set(CMAKE_Fortran_FLAGS "-fPIC -O3")
else()
 set(CMAKE_Fortran_FLAGS "-fPIC -cpp -O3 -m64 -mcmodel=medium")
endif()

# C++ compiler flags

set(CMAKE_CXX_FLAGS "-fPIC -O3 -m64 -pthread -fopenmp")

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
  )

set(SYCL_FULLERENE_LIBS ${SYCL_FULLERENE_LIBS} CACHE STRING "SYCL_FULLERENE_LIBS")

if(ENABLE_SYCL)
  set(PREVIOUS_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER ${DPCPP_COMPILER})
  add_subdirectory(src/sycl)
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




