#pragma once
// Caller-owned scratch workspace for SYCL kernels.
//
// A Workspace is a non-owning view over a contiguous byte buffer plus the
// device the buffer lives on. Storage is typically a SyclVector<std::byte>
// in caller code, sized via the corresponding kernel's *_buffer_size(...)
// query.  Each kernel call wraps the workspace in a fresh BumpAllocator
// (cheap; just a pointer + size) and bump-allocates per-call scratch from
// it.
//
// Design notes:
//  - Pass Workspace by value: it's two pointers and a Device, cheap to copy,
//    no aliasing concerns.
//  - Default-constructed Workspace is empty(); kernels that need zero
//    scratch (e.g. when the impl uses local_accessor for everything) accept
//    a default-constructed Workspace.
//  - Workspace does NOT own its storage. Lifetime of the underlying
//    SyclVector<std::byte> must outlive the kernel call. For pipelines, one
//    persistent SyclVector<std::byte> per pipeline thread, sized to
//    max(buffer_size for every kernel touched), gives the same
//    "persistent across calls" property the old FunctorArrays had.

#include <fullerenes/mempool.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

#include <algorithm>
#include <cstddef>
#include <initializer_list>

struct Workspace {
    std::span<std::byte> bytes{};
    Device          device{};

    constexpr Workspace() = default;
    constexpr Workspace(std::span<std::byte> b, Device d) : bytes(b), device(d) {}

    // Convenience: build a Workspace view over a SyclVector<std::byte>.
    Workspace(SyclVector<std::byte>& storage, Device d)
        : bytes(std::span<std::byte>(storage.data(), storage.size())), device(d) {}

    BumpAllocator bump() const { return BumpAllocator(bytes); }
    size_t        capacity() const { return bytes.size(); }
    bool          empty()    const { return bytes.size() == 0; }
};

// Convenience: pick the maximum buffer-size requirement across multiple
// kernels, the typical pattern for a pipeline (one workspace, reused
// across all kernel calls in a thread).
inline constexpr size_t max_buffer_size(std::initializer_list<size_t> sizes) {
    size_t n = 0;
    for (size_t s : sizes) n = std::max(n, s);
    return n;
}
