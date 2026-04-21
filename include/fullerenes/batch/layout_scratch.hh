#pragma once

// LayoutScratch<T>: stage-local scratch buffer for 2D layout coordinates
// (e.g. Tutte embedding output). Sized as `capacity * Nf` where Nf is the
// per-entry face count (dual vertex count).
//
// Owned separately from Batch<TriangulationView> so that the batch payload
// carries only topology. The Tutte functor consumes a TriangulationView
// batch plus a LayoutScratchView<T>; the spherical-projection functor
// consumes both plus a destination PolyhedronView<T> batch.

#include "fullerenes/batch/storage-policy.hh"
#include "fullerenes/geometry.hh"

#include <cassert>
#include <cstddef>
#include <span>

namespace batch {

template <typename T>
struct LayoutScratchView {
    std::span<coord2<T>> points;   // capacity_ * Nf contiguous 2D coords
    int Nf    = 0;
    int size_ = 0;                 // active entries
    int capacity_ = 0;

    int size()     const { return size_; }
    int capacity() const { return capacity_; }
    int per_entry() const { return Nf; }
    bool empty()   const { return size_ == 0; }

    // Coordinates for entry `i`.
    std::span<coord2<T>> operator[](std::size_t i) const {
        assert(int(i) < size_);
        return points.subspan(std::size_t(Nf) * i, std::size_t(Nf));
    }

    LayoutScratchView slice(std::size_t offset, int count) const {
        assert(int(offset) + count <= size_);
        return LayoutScratchView{
            points.subspan(std::size_t(Nf) * offset, std::size_t(Nf) * count),
            Nf, count, count,
        };
    }
};

template <typename T = double>
class LayoutScratch {
public:
    LayoutScratch() = default;

    LayoutScratch(int Nf, int capacity)
        : Nf_(Nf), size_(0), capacity_(capacity) {
        allocate_();
    }

    void reserve(int new_capacity) {
        if (new_capacity <= capacity_) return;
        capacity_ = new_capacity;
        allocate_();
    }
    void resize(int new_size) {
        assert(new_size <= capacity_);
        size_ = new_size;
    }
    void clear() { size_ = 0; }

    LayoutScratchView<T> view() const {
        return LayoutScratchView<T>{
            std::span<coord2<T>>(const_cast<coord2<T>*>(points_.data()),
                                 std::size_t(Nf_) * capacity_),
            Nf_, size_, capacity_,
        };
    }
    LayoutScratchView<T> slice(std::size_t offset, int count) const {
        return view().slice(offset, count);
    }
    std::span<coord2<T>> operator[](std::size_t i) const { return view()[i]; }

    int  size()      const { return size_; }
    int  capacity()  const { return capacity_; }
    int  per_entry() const { return Nf_; }
    bool empty()     const { return size_ == 0; }

private:
    void allocate_() {
        points_.resize(std::size_t(Nf_) * capacity_);
    }

    int Nf_       = 0;
    int size_     = 0;
    int capacity_ = 0;
    BatchAlloc<coord2<T>> points_{};
};

} // namespace batch
