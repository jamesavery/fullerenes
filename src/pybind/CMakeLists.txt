cmake_minimum_required(VERSION 3.18)

# ---------------------------------------------------------------------------
# Dual-mode build. This file builds the _fullerenes extension two ways:
#   * STANDALONE   -- the dev loop / `pip install` (scikit-build): we are the
#                     top-level project, so we own the toolchain (default clang++
#                     to match the prebuilt libfullerenes.so) and link a PREBUILT
#                     .so found by path.
#   * SUBDIRECTORY -- add_subdirectory(src/pybind) from the parent fullerenes
#                     CMake (ENABLE_PYTHON=ON): we INHERIT the parent's compiler,
#                     standard and global include dirs, and link the in-tree
#                     `fullerenes` target (which propagates its transitive deps).
# Detected from whether we are the top-level project. CMAKE_SOURCE_DIR and
# CMAKE_CURRENT_SOURCE_DIR are both set from -S before project(), so the test is
# valid even above the standalone project() call below.
# ---------------------------------------------------------------------------
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(PYF_STANDALONE TRUE)
else()
  set(PYF_STANDALONE FALSE)
endif()

if(PYF_STANDALONE)
  # ABI: the module exchanges std::vector/std::span/std::string/matrix<T> with
  # libfullerenes.so across the .so boundary, so it MUST be built with the same
  # compiler/standard the library was. Standalone we own the toolchain and
  # default to clang++ (the compiler the prebuilt .so uses); as a subdirectory we
  # never touch it -- we inherit the parent's, which built the library, so the
  # module always matches. This MUST precede project(), which caches the compiler.
  if(NOT CMAKE_CXX_COMPILER AND NOT DEFINED CACHE{CMAKE_CXX_COMPILER})
    set(CMAKE_CXX_COMPILER clang++ CACHE FILEPATH "C++ compiler (match libfullerenes.so)")
  endif()
  project(pyfullerenes LANGUAGES CXX)
  set(CMAKE_CXX_STANDARD 23)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS ON)            # gnu++23
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
  endif()
endif()

# ---------------------------------------------------------------------------
# The library to link, its headers, and its runtime location differ by mode.
#   standalone -> a PREBUILT libfullerenes.so found by path. Parameterised on
#                 FULLERENE_ROOT / FULLERENE_BUILD_DIR for out-of-tree / sdist
#                 builds; the default ../.. is correct at both the old
#                 claude-projects/pybind and the promoted src/pybind (same depth
#                 below the repo root).
#   subdir     -> the in-tree `fullerenes` target. The parent's global
#                 include_directories(include; build/include) already cover the
#                 headers + generated config.hh, and linking the target pulls in
#                 its transitive deps (buckygen, fortran, BLAS/LAPACK, OpenMP).
# ---------------------------------------------------------------------------
if(PYF_STANDALONE)
  get_filename_component(_default_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
  set(FULLERENE_ROOT      "${_default_root}"          CACHE PATH "fullerenes repo root")
  set(FULLERENE_BUILD_DIR "${FULLERENE_ROOT}/build"   CACHE PATH "fullerenes build tree")

  set(FULL_INC        "${FULLERENE_ROOT}/include")
  set(FULL_BUILD_INC  "${FULLERENE_BUILD_DIR}/include")          # generated config.hh
  set(FULL_LIBDIR     "${FULLERENE_BUILD_DIR}/src/c++")
  set(LIBFULLERENES   "${FULL_LIBDIR}/libfullerenes.so")

  if(NOT EXISTS "${LIBFULLERENES}")
    message(FATAL_ERROR
      "libfullerenes.so not found at ${LIBFULLERENES}.\n"
      "FULLERENE_ROOT defaulted to ${FULLERENE_ROOT} (the source dir's ../..). "
      "This works for an in-tree or editable build, but a non-editable / sdist "
      "build (plain `pip install .`) configures from a copied temp dir where ../.. "
      "is wrong. Build the parent library first, and for an out-of-tree build pass "
      "-DFULLERENE_ROOT=/path/to/fullerenes (or -DFULLERENE_BUILD_DIR=.../build).")
  endif()

  set(PYF_LIB          "${LIBFULLERENES}")
  set(PYF_PRIVATE_INC  "${FULL_INC}" "${FULL_BUILD_INC}")
  set(PYF_RPATH_DIR    "${FULL_LIBDIR}")
else()
  set(PYF_LIB          fullerenes)                  # in-tree target
  # The `fullerenes` target carries no usage-requirement include dirs, so the
  # binding's headers are found via the parent's directory-scope
  # include_directories(include; build/include) -- which apply because this
  # subdir is added at root scope after them. If the parent ever moves to
  # target-based includes, add them here explicitly.
  set(PYF_PRIVATE_INC  "")
  set(PYF_RPATH_DIR    "$<TARGET_FILE_DIR:fullerenes>")
endif()

# ---------------------------------------------------------------------------
# Python + pybind11 (both modes; the parent CMake does not look these up). Use
# the pybind11 CMake package -- pass -Dpybind11_DIR=$(python -m pybind11
# --cmakedir); the system pybind11 2.9.1 is too old for Python 3.13.
# ---------------------------------------------------------------------------
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_fullerenes
    src/module.cc
    src/fullerene_graph.cc
    src/fullerene_dual.cc
    src/isomerdb.cc
    src/enumeration.cc
    src/polyhedron.cc
    src/deltahedron.cc
    src/symmetry.cc
    src/bindings_generated.cc)

# Pin gnu++23 on the target itself, so the module matches the library's ABI even
# if the (inherited) parent default ever differs.
set_target_properties(_fullerenes PROPERTIES
    CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS ON)

target_include_directories(_fullerenes PRIVATE
    ${PYF_PRIVATE_INC} ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(_fullerenes PRIVATE ${PYF_LIB} ${CMAKE_DL_LIBS})

set_target_properties(_fullerenes PROPERTIES
    BUILD_RPATH   "${PYF_RPATH_DIR}"
    INSTALL_RPATH "${PYF_RPATH_DIR}")

# ---------------------------------------------------------------------------
# Tests. The full pytest suite (zero-copy, IsomerDB, naming, enumeration,
# geometry, symmetry, the fail-closed binding-drift gate, the .pyi<->module
# gate) registered as one ctest. It registers only when this directory is
# configured -- standalone, or in the parent build under ENABLE_PYTHON=ON -- so
# a parent build with the bindings OFF carries no python tests.
#
# ctest does NOT build: run it after `cmake --build ... --target _fullerenes`
# (or a full build), else the drift gate reads a stale/absent .so. tests/conftest.py
# puts the source tree -- where the NOT-SKBUILD branch below drops the freshly-built
# .so -- on sys.path ahead of any editable install, so the suite always imports the
# just-built module.
#
# The drift gate inside the suite shells out to `gen_bindings.py check`, which reads
# the library and the generated config.hh from the canonical <root>/build (its
# defaults). That is correct for the standard build tree; a library built into a
# different directory would need the gate pointed at it (gen_bindings.py check
# --so <path>; the build-include is likewise assumed at <root>/build/include).
# ---------------------------------------------------------------------------
add_test(NAME pybind_pytest
         COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/tests -q)

# Install into the importable package dir (scikit-build-core places this inside
# the wheel's fullerenes/ package; relative to site-packages).
install(TARGETS _fullerenes DESTINATION fullerenes)

# Ship the type stubs explicitly (don't rely on wheel.packages copying the source
# dir) so PEP 561 type-checkers see them regardless of build/packaging path.
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/fullerenes/py.typed
              ${CMAKE_CURRENT_SOURCE_DIR}/fullerenes/_fullerenes.pyi
        DESTINATION fullerenes)

# For a plain in-tree cmake build (not under scikit-build), also drop the .so
# directly into the source package dir so `import fullerenes` works without an
# install step. Guarded so wheel builds don't pollute the source tree.
if(NOT DEFINED SKBUILD)
  set_target_properties(_fullerenes PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/fullerenes")
endif()
