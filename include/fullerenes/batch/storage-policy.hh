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

#if defined(FULLERENES_ENABLE_SYCL)
#  define BATCH_STORAGE_USE_SYCL 1
#elif defined(SYCL_LANGUAGE_VERSION)
#  define BATCH_STORAGE_USE_SYCL 1
#elif defined(__ACPP__) || defined(__ADAPTIVECPP__) || defined(__HIPSYCL__) || defined(__OPENSYCL__) || defined(SYCL_IMPLEMENTATION_ONEAPI)
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
