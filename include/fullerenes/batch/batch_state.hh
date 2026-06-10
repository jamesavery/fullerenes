#pragma once

// BatchState / BatchStateView: persistent per-entry runtime metadata that is
// kept parallel to a Batch<View> but has no view-level presence. This includes
// buckygen IDs, status flags, iteration counters, and the valid-index mapping
// used by queue filtering.
//
// The split is intentional: Batch<View> owns only the tuple fields of View.
// Anything that is not a topology/geometry field of the view -- e.g. status
// or IDs -- lives here, so view wrappers stay pure data owners and the
// metadata can be queried/filtered independently of the batch payload.

#include "fullerenes/batch/storage-policy.hh"
#include "fullerenes/sycl-headers/sycl-status-enum.hh"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>

namespace batch {

// Non-owning view of BatchState. Trivially copyable.
struct BatchStateView {
    std::span<uint64_t>   id;          // isomer id (e.g. buckygen index)
    std::span<StatusFlag> status;      // status bitset
    std::span<int32_t>    iteration;   // iteration counter (e.g. FF CG iters)
    std::span<int32_t>    valid_index; // compact index mapping for filtering

    int size() const { return int(id.size()); }
    bool empty() const { return id.empty(); }

    BatchStateView slice(std::size_t offset, int count) const {
        assert(int(offset) + count <= size());
        return BatchStateView{
            id.subspan(offset, count),
            status.subspan(offset, count),
            iteration.subspan(offset, count),
            valid_index.subspan(offset, count),
        };
    }
};

// Owning BatchState. Allocates four parallel arrays of length `capacity`.
class BatchState {
public:
    BatchState() = default;

    explicit BatchState(int capacity) : size_(0), capacity_(capacity) {
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

    // Push one entry.  Returns the index of the newly inserted entry.
    int push_back(uint64_t isomer_id,
                  StatusFlag flag      = StatusFlag{},
                  int32_t    iter      = 0,
                  int32_t    valid_idx = -1) {
        assert(size_ < capacity_);
        int i = size_++;
        id_[i]          = isomer_id;
        status_[i]      = flag;
        iteration_[i]   = iter;
        valid_index_[i] = valid_idx;
        return i;
    }

    BatchStateView view() const {
        auto n = std::size_t(size_);
        return BatchStateView{
            std::span<uint64_t>  (const_cast<uint64_t*>  (id_.data()),          n),
            std::span<StatusFlag>(const_cast<StatusFlag*>(status_.data()),      n),
            std::span<int32_t>   (const_cast<int32_t*>   (iteration_.data()),   n),
            std::span<int32_t>   (const_cast<int32_t*>   (valid_index_.data()), n),
        };
    }
    BatchStateView slice(std::size_t offset, int count) const {
        return view().slice(offset, count);
    }

    int  size()     const { return size_; }
    int  capacity() const { return capacity_; }
    bool empty()    const { return size_ == 0; }

    // Direct field accessors (sized to `size_`).
    std::span<uint64_t>   id()          const { return view().id; }
    std::span<StatusFlag> status()      const { return view().status; }
    std::span<int32_t>    iteration()   const { return view().iteration; }
    std::span<int32_t>    valid_index() const { return view().valid_index; }

    // Raw slot access (capacity-sized, bypasses logical size).  Used by
    // BatchQueue<V> which addresses slots directly.
    void write_slot(int slot, uint64_t id_v, StatusFlag flag, int32_t iter, int32_t vi = -1) {
        assert(slot >= 0 && slot < capacity_);
        id_[slot]          = id_v;
        status_[slot]      = flag;
        iteration_[slot]   = iter;
        valid_index_[slot] = vi;
    }
    void read_slot(int slot, uint64_t& id_v, StatusFlag& flag, int32_t& iter) const {
        assert(slot >= 0 && slot < capacity_);
        id_v = id_[slot];
        flag = status_[slot];
        iter = iteration_[slot];
    }

private:
    void allocate_() {
        id_         .resize(capacity_);
        status_     .resize(capacity_);
        iteration_  .resize(capacity_);
        valid_index_.resize(capacity_);
    }

    int size_     = 0;
    int capacity_ = 0;
    BatchAlloc<uint64_t>   id_{};
    BatchAlloc<StatusFlag> status_{};
    BatchAlloc<int32_t>    iteration_{};
    BatchAlloc<int32_t>    valid_index_{};
};

} // namespace batch
