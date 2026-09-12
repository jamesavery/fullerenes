#!/usr/bin/env bash
#
# Configure the fullerenes build on this machine (Ubuntu 22.04, 2x RTX 4090).
#
#   ./cmake_JRHNBR.sh [build-dir] [extra -D args...]   # default: ./build, AOT
#   SYCL_MODE=jit ./cmake_JRHNBR.sh                    # JIT instead of AOT
#
# then:
#
#   cmake --build build -j 16
#
# ---------------------------------------------------------------------------
# Why these settings and not CMake's defaults
# ---------------------------------------------------------------------------
#
# CXX = acpp (AdaptiveCpp).  CMake's default /usr/bin/c++ is g++ 11.4, which
#   has no <expected> -- that landed in GCC 12.  dense_linalg.hh includes it,
#   so a default configure dies on the first library TU with
#   "fatal error: expected: No such file or directory".  acpp wraps clang-20,
#   which has <expected>, and CMake identifies it as Clang, so the
#   -D__cpp_concepts=202002L workaround in CMakeLists.txt applies and
#   libstdc++'s <expected> becomes visible.  Using acpp for the WHOLE build
#   (not just src/sycl) keeps one ABI across the libfullerenes boundary.
#
# SYCL_TARGETS = NVIDIA + cuda + sm_89 -- ahead-of-time (AOT) codegen, with
#   clang-20 as the device compiler.  Kernels are compiled for the RTX 4090 at
#   build time, so nothing is JIT-compiled on first run.  SYCL_MODE=jit gets
#   the JIT build instead.
#
#   AdaptiveCpp can also drive NVIDIA's nvc++ as the device compiler
#   (-DSYCL_CUDA_BACKEND=cuda-nvcxx, still supported by CMakeLists).  This script
#   does not offer it: it has not compiled since HPC SDK 26.5 (2026-09-07) broke
#   its bundled stdexec against libstdc++-12, and before that it silently
#   produced all-NaN hessians.  Pass the -D yourself if you want to retry it.
#
#   Do NOT expect CMakeLists' defaults to do this for you: it picks AMD
#   (hip:gfx90a) on Linux whenever no nvcc is on PATH -- wrong for this box --
#   and its NVIDIA default arch is sm_86, not the 4090's sm_89.
#
#   Why the clang "cuda:" backend needs CUDA_TOOLKIT pinned: clang-20 supports
#   CUDA only up to 12.8, and /usr/local/cuda points at 13.3.  Left to its own
#   devices clang fails inside NVIDIA's headers ("error: expected function body
#   after function declarator" in crt/math_functions.hpp), which is why this
#   script names /usr/local/cuda-12.8 explicitly rather than trusting the
#   default.  A newer clang would lift the ceiling; until then, keep the pin.
#
# CC = clang-20.  buckygen and mgmres are C; matching clang-20 keeps the
#   toolchain uniform.  gcc would also work (C ABI), it is just tidier.
#
# Fortran = gfortran (11.4).  There is no clang Fortran here, and gfortran
#   links fine against clang-built C++.
#
set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${1:-$SRC_DIR/build}"

ACPP=/opt/adaptivecpp/bin/acpp
CLANG_C=/usr/lib/llvm-20/bin/clang
CUDA_TOOLKIT=/usr/local/cuda-12.8  # clang-20 supports CUDA only up to 12.8
CUDA_ARCH=89                       # RTX 4090.  3090 = 86, A100 = 80, H100 = 90.
SYCL_MODE="${SYCL_MODE:-aot}"      # aot | jit

if [ ! -x "$ACPP" ]; then
    echo "error: no AdaptiveCpp compiler at $ACPP" >&2
    echo "       install it, or point ACPP at your acpp/syclcc." >&2
    exit 1
fi
if [ ! -x "$CLANG_C" ]; then
    echo "error: no clang-20 C compiler at $CLANG_C" >&2
    exit 1
fi

case "$SYCL_MODE" in
    aot)
        # Kernels compiled to sm_$CUDA_ARCH at build time, with clang-20 as the
        # device compiler.  It needs a CUDA toolkit clang understands: clang-20
        # supports CUDA only up to 12.8, while every HPC SDK toolkit here is
        # 13.x.  CMakeLists derives --acpp-cuda-path from CMAKE_CUDA_COMPILER in
        # this branch, so pointing that at 12.8's nvcc is what selects the right
        # headers and libdevice.
        if [ ! -x "$CUDA_TOOLKIT/bin/nvcc" ]; then
            echo "error: AOT mode needs a CUDA toolkit <= 12.8, not found at" >&2
            echo "       $CUDA_TOOLKIT" >&2
            echo "       Install one, or re-run with:  SYCL_MODE=jit $0" >&2
            exit 1
        fi
        SYCL_ARGS=(-DSYCL_TARGETS=NVIDIA
                   -DSYCL_CUDA_BACKEND=cuda
                   -DSYCL_CUDA_ARCH="$CUDA_ARCH"
                   -DCMAKE_CUDA_COMPILER="$CUDA_TOOLKIT/bin/nvcc")
        ;;
    jit)
        # SSCP: one portable IR in the binary, JIT-compiled to the device on
        # first run.  Builds much faster, costs a JIT pass at startup.
        #
        # Named jit here for the axis this switch actually selects -- when the
        # kernels are compiled.  AdaptiveCpp calls the same thing "generic"
        # (acpp-core.json default-targets, --acpp-targets=generic), after the
        # other property it has: one binary that runs on any backend.  That is
        # the name to search their docs for, and it is what SYCL_TARGETS below
        # still takes.
        SYCL_ARGS=(-DSYCL_TARGETS=GENERIC)
        ;;
    *)
        echo "error: SYCL_MODE must be 'aot' or 'jit', got '$SYCL_MODE'" >&2
        exit 1
        ;;
esac

# CMake refuses to change CMAKE_CXX_COMPILER in an existing cache, so drop the
# cache and its compiler probes.  Deliberately NOT 'rm -rf $BUILD_DIR': any
# benchmark output or checkpoint you left in there survives.  Object files are
# stale after a compiler switch but CMake rebuilds them.
rm -f  "$BUILD_DIR/CMakeCache.txt"
rm -rf "$BUILD_DIR/CMakeFiles"

cmake -S "$SRC_DIR" -B "$BUILD_DIR" \
      -DCMAKE_CXX_COMPILER="$ACPP" \
      -DCMAKE_C_COMPILER="$CLANG_C" \
      -DCMAKE_Fortran_COMPILER=gfortran \
      -DENABLE_SYCL=ON \
      "${SYCL_ARGS[@]}" \
      "${@:2}"

# --- post-configure warnings ------------------------------------------------
#
# clang-20 needs its own omp.h (package libomp-20-dev).  The libomp runtime is
# already installed, so only the header is missing -- and find_package(OpenMP)
# then fails SILENTLY: the if(OpenMP_CXX_FOUND) guards in src/c++, benchmarks
# and programs just skip the pragmas, and the ~8 files using #pragma omp run
# single-threaded with no error anywhere.  Check loudly instead.
if grep -q '^OpenMP_CXX_FLAGS:STRING=NOTFOUND' "$BUILD_DIR/CMakeCache.txt" 2>/dev/null; then
    echo
    echo "WARNING: OpenMP was NOT found -- every #pragma omp region will run"
    echo "         serially, with no further warning at build time."
    echo "         Fix:  sudo apt install libomp-20-dev   (then re-run this script)"
fi

# The isomer database is gitignored external data (.gitignore:17).  Without it
# a large fraction of the test suite fails with IsomerDB::readPDB errors.
if [ ! -d "$SRC_DIR/database/All" ]; then
    echo
    echo "WARNING: $SRC_DIR/database/All is missing.  IsomerDB-dependent tests"
    echo "         and tools will throw 'cannot open database file'."
fi

echo
if [ "$SYCL_MODE" = aot ]; then
    echo "Configured $BUILD_DIR  (SYCL: AOT, $(grep '^ACPP_TARGETS' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2))."
    echo "Kernels are compiled to sm_$CUDA_ARCH at build time by clang-20 against"
    echo "$CUDA_TOOLKIT -- src/sycl is slower to build, and nothing is"
    echo "JIT-compiled at first run."
else
    echo "Configured $BUILD_DIR  (SYCL: JIT, AdaptiveCpp generic/SSCP -- kernels"
    echo "compiled to the device at first run)."
fi
echo "Build with:"
echo "    cmake --build $BUILD_DIR -j $(nproc)"
