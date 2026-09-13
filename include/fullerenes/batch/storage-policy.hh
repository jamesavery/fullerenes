#pragma once

// Storage backend policy for Batch<View>.
//
// Selects between USM-managed allocation (SYCL builds) and
// std::vector (CPU-only builds).
//
// Use BatchAlloc<T> as the owning container type in Batch<View> and
// Owned<View>. It provides:
//   - resize(n) / resize(n, value)
//   - data()  -> T*
//   - size()  -> size_t
//   - std::span<T> conversion
//
// Keep all higher-level code backend-agnostic by going through this alias.

// The choice is the BUILD's, recorded in the generated config.hh by
// ENABLE_SYCL, never inferred from which compiler is reading this header: a
// translation unit compiled by acpp in an ENABLE_SYCL=OFF tree, or by clang
// against an ENABLE_SYCL=ON library, must still agree with the library on
// Batch's layout.  (Before 2026-09-06 the acpp/SYCL predefined macros also
// selected SyclVector, so acpp as host compiler with SYCL off produced
// objects that could not link.)
#include "fullerenes/config.hh"

#if defined(FULLERENES_ENABLE_SYCL)
#  define BATCH_STORAGE_USE_SYCL 1
#else
#  define BATCH_STORAGE_USE_SYCL 0
#endif

#if BATCH_STORAGE_USE_SYCL
#  include <fullerenes/sycl-headers/sycl-vector.hh>
   template<typename T>
   using BatchAlloc = SyclVector<T>;
#else
#  include <vector>
   template<typename T>
   using BatchAlloc = std::vector<T>;
#endif

#include <span>
