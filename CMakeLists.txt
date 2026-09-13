
cmake_minimum_required(VERSION 3.5)
include(CheckLanguage)
include(CMakeDependentOption)
project(Fullerenes)
enable_language(CXX)

# Core-only mode: build just libfullerenes + libbuckygen (no Fortran,
# playground, programs, benchmarks, tests, SYCL, pybind).  This is the mode a
# parent project uses when consuming this repo as a git submodule via
# add_subdirectory.  NOTE for parents: cmake_minimum_required(3.5) leaves
# CMP0077 at OLD in this directory, so set the gate as a cache variable:
#   set(FULLERENES_CORE_ONLY ON CACHE BOOL "" FORCE)
option(FULLERENES_CORE_ONLY
  "Build only libfullerenes + buckygen (skip Fortran, playground, programs, benchmarks, tests, SYCL, pybind)"
  OFF)
# A parent that links libfullerenes into its own executables, and would rather
# not ship the library beside them, can build it static -- again as a cache
# variable:  set(FULLERENES_SHARED OFF CACHE BOOL "" FORCE)
option(FULLERENES_SHARED
  "Build libfullerenes as a shared library (OFF: a static archive)"
  ON)

# The legacy Fortran program (src/fortran: the Schwerdtfeger-Wirz-Avery
# "Fullerene" v4.5) and graph_fortran.cc, the C++ side of its ABI.  Nothing
# else in the tree links or calls Fortran -- every algorithm it carried that
# the library needs has a C++ home (claude-projects/unfortran) -- so it is
# opt-in: the default build, tests included, needs no Fortran compiler.  Off
# means enable_language(Fortran) never runs.  Forced off under core-only.
cmake_dependent_option(FULLERENES_LEGACY_FORTRAN
  "Build the legacy Fortran program and its C ABI (needs a Fortran compiler)"
  OFF "NOT FULLERENES_CORE_ONLY" OFF)

if(FULLERENES_LEGACY_FORTRAN)
  enable_language(Fortran)
  set(FORTRAN_NMAX 5000 CACHE STRING "Maximum vertex count for Fortran routines using statically allocated memory")
endif()
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  # Only when this repo is the top-level project: never rewrite a parent
  # project's generator cache from a submodule.
  set(CMAKE_GENERATOR "Unix Makefiles" CACHE INTERNAL "" FORCE)
endif()
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(FULLERENE_DATABASE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/database/ CACHE PATH "Path to fullerene database root")

find_package(GSL)
find_package(GTest)
find_library(TBB_LIB tbb)


enable_testing()
if(GTest_FOUND)
  include_directories(${GTEST_INCLUDE_DIRS})
else()
  message(STATUS "GTest not found -- GTest-based unit tests will be skipped")
endif()

if(TBB_LIB)
  message(STATUS "Found TBB runtime: ${TBB_LIB}")
else()
  message(STATUS "Could not find TBB runtime library (tbb)")
endif()

# Shared between CUDA and SYCL
set(GPU_MAXREGCOUNT 80)

# --- Which GPU is actually in this machine? ---------------------------------
# The evidence is what the vendor tools REPORT, not whether they are installed:
# a box can carry rocminfo and /opt/rocm from a package with no AMD GPU in it
# (rocminfo then exits nonzero: "ROCk module is NOT loaded"), and NVIDIA's
# nvc++ backend does not put nvcc on PATH, so check_language(CUDA) says
# nothing about the two RTX cards it drives.  So: an NVIDIA GPU is one that
# nvidia-smi lists with a compute capability; an AMD GPU is a gfx agent
# rocminfo lists (gfx000 is a CPU agent).  Each yields the architecture too,
# which becomes the default for the *_ARCH knobs below; an explicit
# -D on any of them still wins, these are only the defaults.
find_program(NVIDIA_SMI nvidia-smi)
find_program(ROCMINFO rocminfo HINTS ENV ROCM_PATH /opt/rocm PATH_SUFFIXES bin)
mark_as_advanced(NVIDIA_SMI ROCMINFO)
set(NVIDIA_GPU_ARCH "")   # e.g. 89 (sm_89), from the first GPU nvidia-smi lists
set(AMD_GPU_ARCH "")      # e.g. gfx90a, from the first GPU agent rocminfo lists
if(NVIDIA_SMI)
  execute_process(COMMAND ${NVIDIA_SMI} --query-gpu=compute_cap --format=csv,noheader
                  RESULT_VARIABLE _rc OUTPUT_VARIABLE _out ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(_rc EQUAL 0 AND _out MATCHES "^([0-9]+)\\.([0-9]+)")
    set(NVIDIA_GPU_ARCH "${CMAKE_MATCH_1}${CMAKE_MATCH_2}")
  endif()
endif()
if(ROCMINFO)
  execute_process(COMMAND ${ROCMINFO} RESULT_VARIABLE _rc OUTPUT_VARIABLE _out ERROR_QUIET)
  if(_rc EQUAL 0 AND _out MATCHES "(gfx[1-9][0-9a-f]*)")
    set(AMD_GPU_ARCH "${CMAKE_MATCH_1}")
  endif()
endif()
if(NVIDIA_GPU_ARCH)
  message(STATUS "GPU: NVIDIA, compute capability sm_${NVIDIA_GPU_ARCH} (nvidia-smi)")
endif()
if(AMD_GPU_ARCH)
  message(STATUS "GPU: AMD, ${AMD_GPU_ARCH} (rocminfo)")
endif()
if(NOT NVIDIA_GPU_ARCH AND NOT AMD_GPU_ARCH)
  message(STATUS "GPU: none reported by nvidia-smi or rocminfo")
endif()

# Do we have CUDA compilers?
check_language(CUDA)
cmake_dependent_option(ENABLE_CUDA "Build CUDA GPU code" ON "CMAKE_CUDA_COMPILER" OFF)
if(CMAKE_CUDA_COMPILER)
  if(NVIDIA_GPU_ARCH)
    set(_cuda_arch_default ${NVIDIA_GPU_ARCH})
  else()
    set(_cuda_arch_default 86)
  endif()
  set(CMAKE_CUDA_ARCHITECTURES ${_cuda_arch_default} CACHE STRING "Architecture for your NVIDIA GPU: 50|60|70|75|80|86|89|90 (default: what nvidia-smi reports). See https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/")
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
if(FULLERENES_CORE_ONLY)
  # Core-only must not let an auto-detected acpp flip the library ABI via
  # FULLERENES_ENABLE_SYCL in config.hh below.
  set(ENABLE_SYCL OFF)
endif()

if(SYCL_COMPILER)
  message(STATUS "Found AdaptiveCpp compiler: ${SYCL_COMPILER}")
else()
  message(STATUS "Could not find an AdaptiveCpp compiler (looked for: acpp, syclcc)")
endif()

# --- Python bindings (optional) -------------------------------------------
# Build the pybind11 module under src/pybind. Built as a subdirectory, it
# INHERITS this build's compiler/standard/include dirs (it crosses the
# libfullerenes .so ABI boundary, so it must match the library). Opt-in. Extra
# deps when ON: pybind11 >= 2.13 (pass -Dpybind11_DIR=$(python -m pybind11
# --cmakedir)) and the Python development module -- find_package(... REQUIRED) in
# src/pybind, so configure FAILS FAST if either is missing (like ENABLE_SYCL).
# The binding-drift gate additionally wants clang + Python 3.11+ (tomllib) but
# self-skips if clang is absent. The python ctests register only when this is ON.
option(ENABLE_PYTHON "Build the Python (pybind11) bindings (src/pybind)" OFF)

if(ENABLE_SYCL)
  if(NOT SYCL_COMPILER)
    message(FATAL_ERROR "ENABLE_SYCL=ON requires AdaptiveCpp compiler (acpp or syclcc)")
  endif()
  # batch::Batch storage is SyclVector (USM) in every translation unit of a
  # SYCL build -- host-compiled ones (programs/, benchmarks/) included, so
  # host-populated batch buffers are accessible to SYCL kernels.  The switch
  # is FULLERENES_ENABLE_SYCL in the generated config.hh (configure_file
  # below), which out-of-tree consumers see too; it used to be a compile
  # definition here, invisible to them.
  set(FULLERENES_ENABLE_SYCL ON)

  # The default target is the GPU the machine reports (above); a CUDA
  # toolchain with no GPU listed still means NVIDIA (it compiles); and with
  # neither, the CPU backend, which always builds -- never a guess at a
  # vendor that would then demand its toolkit.  -DSYCL_TARGETS= overrides.
  if(NVIDIA_GPU_ARCH OR CMAKE_CUDA_COMPILER)
    set(SYCL_DEFAULT_TARGET "NVIDIA")
  elseif(AMD_GPU_ARCH)
    set(SYCL_DEFAULT_TARGET "AMD")
  else()
    set(SYCL_DEFAULT_TARGET "X86")
  endif()

  # SYCL_TARGETS takes either a convenience preset -- AMD, NVIDIA, X86, GENERIC,
  # expanded below using the architecture the machine reports -- or an
  # AdaptiveCpp target list passed through verbatim, so everything acpp can
  # express is reachable without adding an option per capability:
  #
  #   -DSYCL_TARGETS=NVIDIA                 the detected NVIDIA GPU
  #   -DSYCL_TARGETS=cuda:sm_80,sm_89       one library for A100 and RTX 4090
  #   -DSYCL_TARGETS=omp+cuda:sm_89         CPU and NVIDIA in one build
  #   -DSYCL_TARGETS=generic                portable IR, JIT at first run
  #
  # Backends are separated by '+' or ';' and architectures within a backend by
  # ','.  acpp's own separator is ';', but that is CMake's list separator and
  # the shell's command separator, so -DSYCL_TARGETS=omp;cuda:sm_89 needs
  # escaping in both.  '+' needs none; the two are accepted interchangeably.
  set(SYCL_TARGETS ${SYCL_DEFAULT_TARGET} CACHE STRING "AdaptiveCpp targets: a preset (AMD|NVIDIA|X86|GENERIC) or a target list such as cuda:sm_80,sm_89 or omp+cuda:sm_89 (default: the GPU nvidia-smi/rocminfo report, else X86)")

  set(SYCL_BASE_FLAGS -O3)
  set(SYCL_EXTRA_FLAGS "" CACHE STRING "Extra SYCL compiler flags")

  # A ';'-separated value arrives as a CMake list; rejoin it, then normalise
  # both accepted separators to '+' so the rest of this block sees one string.
  string(REPLACE ";" "+" _sycl_targets_in "${SYCL_TARGETS}")
  string(TOUPPER "${_sycl_targets_in}" _sycl_targets_preset)

  if(_sycl_targets_preset STREQUAL "AMD")
    if(AMD_GPU_ARCH)
      set(_hip_arch_default ${AMD_GPU_ARCH})
    else()
      set(_hip_arch_default "gfx90a")
    endif()
    set(SYCL_HIP_ARCH ${_hip_arch_default} CACHE STRING "HIP GPU architecture(s) for the AMD preset, comma-separated (default: what rocminfo reports)")
    set(ACPP_TARGETS_RESOLVED "hip:${SYCL_HIP_ARCH}")
  elseif(_sycl_targets_preset STREQUAL "NVIDIA")
    if(NVIDIA_GPU_ARCH)
      set(_cuda_arch_default ${NVIDIA_GPU_ARCH})
    else()
      set(_cuda_arch_default 86)
    endif()
    set(SYCL_CUDA_ARCH ${_cuda_arch_default} CACHE STRING "CUDA architecture(s) for the NVIDIA preset, comma-separated, e.g. 89 or 80,90 (default: what nvidia-smi reports)")
    set(SYCL_CUDA_BACKEND "cuda-nvcxx" CACHE STRING "AdaptiveCpp NVIDIA backend (cuda|cuda-nvcxx)")
    # "80,90" -> "sm_80,sm_90": acpp wants the prefix on every architecture.
    string(REPLACE "," ";" _cuda_arch_list "${SYCL_CUDA_ARCH}")
    set(_cuda_arch_spec "")
    foreach(_arch IN LISTS _cuda_arch_list)
      string(STRIP "${_arch}" _arch)
      if(_cuda_arch_spec)
        string(APPEND _cuda_arch_spec ",sm_${_arch}")
      else()
        set(_cuda_arch_spec "sm_${_arch}")
      endif()
    endforeach()
    set(ACPP_TARGETS_RESOLVED "${SYCL_CUDA_BACKEND}:${_cuda_arch_spec}")
  elseif(_sycl_targets_preset STREQUAL "X86")
    # library-only backend sidesteps known accelerated CPU compiler crashes.
    set(ACPP_TARGETS_RESOLVED "omp.library-only")
  elseif(_sycl_targets_preset STREQUAL "GENERIC")
    set(ACPP_TARGETS_RESOLVED "generic")
  else()
    # Not a preset: an acpp target list, passed through as given.
    string(REPLACE "+" ";" ACPP_TARGETS_RESOLVED "${_sycl_targets_in}")
  endif()

  # Every backend must be one acpp knows, so a mistyped preset is caught here
  # rather than passed through as a target list: -DSYCL_TARGETS=NVIDA would
  # otherwise configure cleanly, leave SOLVER_BACKEND undefined and
  # --acpp-cuda-path underived, and surface much later as an opaque acpp error.
  foreach(_backend IN LISTS ACPP_TARGETS_RESOLVED)
    if(NOT _backend MATCHES "^(omp|cuda|cuda-nvcxx|hip|generic)([.:]|$)")
      message(FATAL_ERROR
        "SYCL_TARGETS: '${_backend}' (from '${SYCL_TARGETS}') is neither a "
        "preset (AMD|NVIDIA|X86|GENERIC) nor an AdaptiveCpp backend "
        "(omp|cuda|cuda-nvcxx|hip|generic).  Check the spelling.")
    endif()
  endforeach()

  # generic is the SSCP flow, not another backend to add alongside the rest:
  # combining it with an ahead-of-time target fails at link time, so reject it
  # here where the message can say why.
  list(LENGTH ACPP_TARGETS_RESOLVED _n_backends)
  if(_n_backends GREATER 1)
    foreach(_backend IN LISTS ACPP_TARGETS_RESOLVED)
      if(_backend MATCHES "^[Gg]eneric")
        message(FATAL_ERROR
          "SYCL_TARGETS: 'generic' cannot be combined with other targets "
          "(got '${SYCL_TARGETS}').  It selects the portable-IR SSCP flow that "
          "JITs at runtime, and linking it beside an ahead-of-time target "
          "fails.  Use generic alone, or list only ahead-of-time targets.")
      endif()
    endforeach()
  endif()

  # What the rest of the build needs to know about the resolved targets.  These
  # replace tests on the preset name, which cannot answer the question once the
  # value is a list.
  set(SYCL_HAS_NVCXX OFF)
  set(SYCL_HAS_CUDA  OFF)
  set(SYCL_HAS_HIP   OFF)
  foreach(_backend IN LISTS ACPP_TARGETS_RESOLVED)
    if(_backend MATCHES "^cuda-nvcxx")
      set(SYCL_HAS_NVCXX ON)
      set(SYCL_HAS_CUDA ON)
    elseif(_backend MATCHES "^cuda")
      set(SYCL_HAS_CUDA ON)
    elseif(_backend MATCHES "^hip")
      set(SYCL_HAS_HIP ON)
    endif()
  endforeach()

  # One dense-solver backend per build; CUDA wins if a list names both.
  if(SYCL_HAS_CUDA)
    add_definitions(-DSOLVER_BACKEND=CUSOLVER_BACKEND)
  elseif(SYCL_HAS_HIP)
    add_definitions(-DSOLVER_BACKEND=ROCSOLVER_BACKEND)
  endif()

  # Keep backend target selection aligned with SYCL_TARGETS to avoid stale cache values.
  set(ACPP_TARGETS "${ACPP_TARGETS_RESOLVED}" CACHE STRING "AdaptiveCpp backend targets" FORCE)
  message(STATUS "SYCL compile enabled with targets: ${ACPP_TARGETS}")

  # Several architectures in one backend (cuda:sm_80,sm_89) are comma-separated
  # and need no care.  Several BACKENDS are separated by ';', which is also
  # CMake's list separator, so the targets flag is kept OUT of SYCL_FLAGS and
  # handed to consumers on its own, with the separator backslash-escaped.
  # Consumers must expand it QUOTED -- "${SYCL_TARGETS_FLAG}".  That is the
  # whole trick: an escaped ';' survives being stored in COMPILE_OPTIONS and
  # written to the build rule, but an unquoted ${...} expansion on the way
  # there consumes the escape and splits the element, after which every backend
  # past the first reaches the compiler as a stray input file
  # ("--acpp-targets=omp cuda:sm_89").
  string(REPLACE ";" "\\;" _acpp_targets_arg "${ACPP_TARGETS}")
  set(SYCL_TARGETS_FLAG "--acpp-targets=${_acpp_targets_arg}")
  set(SYCL_FLAGS ${SYCL_BASE_FLAGS})
  if(SYCL_HAS_CUDA AND NOT SYCL_HAS_NVCXX AND CMAKE_CUDA_COMPILER)
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
  if(SYCL_HAS_NVCXX)
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


# Get the current working branch (CURRENT_SOURCE_DIR: report this repo's SHA
# even when consumed as a submodule)
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# C++ configure file
include_directories(AFTER SYSTEM include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/)

# Fortran compiler flags
if(FULLERENES_LEGACY_FORTRAN)
  if(DEFINED ENV{CRAYPE_VERSION})
    message("Cray compiler environment detected")
    set(CMAKE_Fortran_FLAGS "-fPIC -O3")
  else()
    set(CMAKE_Fortran_FLAGS "-fPIC -cpp -O3 -m64 -mcmodel=small")
  endif()
endif()

# C++ compiler flags

# `-m64` means x86_64 to Clang's Darwin driver. It was harmless on native
# Intel/Linux builds but silently produces x86_64 objects when a parent uses
# osxcross for arm64. CMAKE_OSX_ARCHITECTURES already carries the target.
if(CF_OSXCROSS_BUILD OR
   (APPLE AND CMAKE_OSX_ARCHITECTURES MATCHES "(^|;)arm64(;|$)") OR
   CMAKE_CXX_COMPILER_TARGET MATCHES "^arm64-apple-")
  set(_FULLERENES_M64 "")
else()
  set(_FULLERENES_M64 "-m64")
endif()
set(CMAKE_CXX_FLAGS "-fPIC -O3 ${_FULLERENES_M64} -pthread")

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # Clang (through at least 18.x) still reports __cpp_concepts=201907L (the
  # C++20 value) even in -std=gnu++23 mode -- it hasn't implemented the core
  # defect-report resolutions libstdc++ checks for. libstdc++'s <expected>
  # header gates its entire contents on __cpp_concepts >= 202002L, so
  # std::expected silently vanishes under Clang+libstdc++ even though
  # Clang's concepts support is otherwise sufficient. Force the macro so
  # <expected> becomes visible; this is a known, widely used workaround and
  # does not paper over any real missing capability.
  string(APPEND CMAKE_CXX_FLAGS " -Wno-builtin-macro-redefined -D__cpp_concepts=202002L")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/fullerenes/config.hh.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/fullerenes/config.hh
  )

if(FULLERENES_LEGACY_FORTRAN)
  # Fortran configure file.  NOTE: writes into the SOURCE tree -- gated so a
  # checkout that never builds the Fortran stays clean.
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fullerenes/config.f.in
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fortran/config.f
  )
endif()


if(FULLERENES_LEGACY_FORTRAN)
  add_subdirectory(src/fortran)
endif()
add_subdirectory(src/contrib)
if(NOT FULLERENES_CORE_ONLY)
  add_subdirectory(playground)
endif()
add_subdirectory(src/c++)

set(SYCL_FULLERENE_LIBS 
  sycl_fullerene_lib 
  )

set(SYCL_FULLERENE_LIBS ${SYCL_FULLERENE_LIBS} CACHE STRING "SYCL_FULLERENE_LIBS")

if(ENABLE_SYCL)
  set(PREVIOUS_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  set(PREVIOUS_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  set(PREVIOUS_CXX_STANDARD ${CMAKE_CXX_STANDARD})
  if(SYCL_HAS_NVCXX)
    # This backend drives NVIDIA's nvc++ (not clang) for the SYCL translation
    # units, but CMAKE_CXX_FLAGS was computed against the HOST compiler.  When
    # the host compiler is clang -- or acpp, which CMake identifies as Clang --
    # the block above appended clang-only switches that nvc++ rejects outright
    # ("nvc++-Error-Unknown switch: -Wno-builtin-macro-redefined"), failing
    # every SYCL TU.  Strip them: the __cpp_concepts override exists solely to
    # unhide libstdc++'s <expected> under clang, and nvc++ neither needs nor
    # understands it.
    string(REPLACE "-Wno-builtin-macro-redefined" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    string(REPLACE "-D__cpp_concepts=202002L" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    string(REPLACE "-fopenmp" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    string(REGEX REPLACE " +" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    # acpp's bundled compiler for the clang-based backends (generic/SSCP, AMD hip,
    # cuda clang) is clang-14, which does not accept CMake's gcc-derived
    # -std=gnu++23 (clang-14 only knows the pre-standard spelling gnu++2b for
    # C++23). CMake computed the standard flag against the host compiler and does
    # not re-derive it when CMAKE_CXX_COMPILER is swapped to acpp below. So
    # disable the automatic standard flag for the acpp subdirs and request
    # gnu++2b explicitly; the X86/omp.library-only backend's host gcc also
    # accepts gnu++2b, so this is safe across all non-nvcxx backends.
    set(CMAKE_CXX_STANDARD)
    string(APPEND CMAKE_CXX_FLAGS " -std=gnu++2b")
  endif()
  set(CMAKE_CXX_COMPILER ${SYCL_COMPILER})
  add_subdirectory(src/sycl)
  add_subdirectory(tests/sycl-tests)
  add_subdirectory(benchmarks/sycl)
  set(CMAKE_CXX_COMPILER ${PREVIOUS_CXX_COMPILER})
  set(CMAKE_CXX_FLAGS ${PREVIOUS_CXX_FLAGS})
  set(CMAKE_CXX_STANDARD ${PREVIOUS_CXX_STANDARD})
endif()

message("C++ FLAGS: " ${CMAKE_CXX_FLAGS})
message("SYCL_FLAGS: " ${SYCL_FLAGS})
message("FLAGS: " ${FLAGS})


#set(SYCL_FLAGS "${SYCL_FLAGS} -tis")

#get_target_property(TESTVAR sycl_fullerene_lib LINK_FLAGS)
#message("sycl_fullerene_lib LINK_FLAGS: " ${TESTVAR})

if(NOT FULLERENES_CORE_ONLY)
  add_subdirectory(programs)
  add_subdirectory(benchmarks)
  add_subdirectory(tests)
endif()

# Python bindings last: src/pybind links the `fullerenes` target (defined in
# add_subdirectory(src/c++) above) and registers its own ctests.
if(ENABLE_PYTHON AND NOT FULLERENES_CORE_ONLY)
  add_subdirectory(src/pybind)
endif()


