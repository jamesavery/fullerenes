#pragma once

// Compatibility shim for C++17 parallel algorithms.
// Apple libc++ does not implement <execution> policies.

#include <algorithm>
#include <numeric>

#if __has_include(<execution>)
  #include <execution>
#endif

#if defined(__cpp_lib_parallel_algorithm) && __cpp_lib_parallel_algorithm >= 201603L
  #define FULLERENE_PAR_UNSEQ std::execution::par_unseq,
#else
  #define FULLERENE_PAR_UNSEQ
#endif
