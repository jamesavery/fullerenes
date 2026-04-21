#pragma once

// BatchQueue<V>: generic circular queue over any batchable view.
//
// Composition (not inheritance): a BatchQueue<V> owns a Batch<V> plus a
// BatchState (the persistent per-entry metadata from Phase 5) aligned at
// identical capacity. The queue adds two integer indices (front_, back_)
// to express circular-storage semantics: non-empty entries form a
// contiguous range in modular index space.
//
// Transfer helpers (batch::queue_push) move entries between a batch
// and a queue, selecting with a predicate over (StatusFlag, iteration)
// and optionally branding the source with a "consumed" status after
// the transfer. Replaces FullereneQueue and QueueUtil in the new layer.
//
// This phase is pure host code: the loops run serially over host-backed
// BatchAlloc storage. SYCL-kernelized transfer is a later phase; the
// interface is preserved so those kernels can back the same entry
// points without rewiring callers.

#include "fullerenes/batch/batch.hh"
#include "fullerenes/batch/batch_state.hh"
#include "fullerenes/batch/batchable.hh"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <type_traits>
#include <utility>

namespace batch {

// ---------------------------------------------------------------------------
// Default predicate: match every non-EMPTY entry.
// ---------------------------------------------------------------------------
struct MatchAll {
    bool operator()(const StatusFlag&, int32_t /*iteration*/) const { return true; }
};

// Match entries whose status has all of `required_on` set and none of
// `required_off`. Optional iteration range [iter_min, iter_max].
struct StatusPredicate {
    int  required_on   = 0;
    int  required_off  = 0;
    int32_t iter_min   = std::numeric_limits<int32_t>::min();
    int32_t iter_max   = std::numeric_limits<int32_t>::max();

    StatusPredicate() = default;
    StatusPredicate(StatusEnum on, StatusEnum off = StatusEnum(0))
        : required_on(int(on)), required_off(int(off)) {}

    bool operator()(const StatusFlag& f, int32_t iter) const {
        int v = static_cast<int>(f);
        if ((v & required_on)  != required_on)  return false;
        if ((v & required_off) != 0)            return false;
        if (iter < iter_min || iter > iter_max) return false;
        return true;
    }
};

// ---------------------------------------------------------------------------
// BatchQueue<V>
// ---------------------------------------------------------------------------
template <class V>
    requires batchable_view<V>
class BatchQueue {
public:
    using view_type = V;

    BatchQueue() = default;

    BatchQueue(int N, int capacity, int dmax = V::default_dmax)
        : batch_(N, capacity, dmax), state_(capacity),
          front_(-1), back_(-1) {}

    // --- Accessors ---
    int  N()        const { return batch_.N(); }
    int  dmax()     const { return batch_.dmax(); }
    int  capacity() const { return batch_.capacity(); }
    int  size()     const { return size_; }
    bool empty()    const { return size_ == 0; }
    int  front_index() const { return front_; }
    int  back_index()  const { return back_; }

    // Raw storage access (slot-indexed, not logical-indexed).
    Batch<V>&          storage()       { return batch_; }
    const Batch<V>&    storage() const { return batch_; }
    BatchState&        state()         { return state_; }
    const BatchState&  state()   const { return state_; }

    // Logical index -> physical slot in the circular buffer.
    int slot_of(int logical) const {
        assert(logical >= 0 && logical < size_);
        return (front_ + logical) % capacity();
    }

    V at(int logical) const {
        return batch_.view_capacity()[std::size_t(slot_of(logical))];
    }

    void clear() {
        size_  = 0;
        front_ = -1;
        back_  = -1;
        batch_.clear();
        state_.clear();
        // Make all slots EMPTY so predicates see a clean queue.
        auto S = state_.view();
        // state_.clear() already reset size_; but the underlying storage still
        // holds previous values. Callers iterating `capacity()` should only
        // trust entries between front_ and back_, so clearing size is enough.
        (void)S;
    }

    // Append one entry, optionally with user-provided metadata.
    // Grows capacity (doubling) if the queue is full.
    int push_back(const V& v,
                  uint64_t id           = 0,
                  StatusFlag flag       = StatusFlag{},
                  int32_t iteration     = 0) {
        if (size_ == capacity()) grow_(capacity() == 0 ? 1 : capacity() * 2);
        int slot;
        if (size_ == 0) {
            slot = 0;
            front_ = 0;
            back_  = 0;
        } else {
            slot = (back_ + 1) % capacity();
            back_ = slot;
        }
        write_entry_(slot, v, id, flag, iteration);
        ++size_;
        return slot;
    }

    // Remove and return the logical-front entry into `out_v`, `out_id`,
    // `out_flag`, `out_iter`.  Returns true on success.
    bool pop_front(V& out_v,
                   uint64_t& out_id,
                   StatusFlag& out_flag,
                   int32_t& out_iter) {
        if (empty()) return false;
        int slot = front_;
        read_entry_(slot, out_v, out_id, out_flag, out_iter);
        --size_;
        if (size_ == 0) {
            front_ = -1;
            back_  = -1;
        } else {
            front_ = (front_ + 1) % capacity();
        }
        return true;
    }

    // Discard the logical-front entry without copying it out.
    bool discard_front() {
        if (empty()) return false;
        --size_;
        if (size_ == 0) {
            front_ = -1;
            back_  = -1;
        } else {
            front_ = (front_ + 1) % capacity();
        }
        return true;
    }

    // Grow to at least `new_capacity`, preserving logical ordering so that
    // after resize front_ == 0 and entries occupy [0, size_).
    void resize(int new_capacity) {
        if (new_capacity <= capacity()) return;
        grow_(new_capacity);
    }

private:
    // Copy one entry's fields from `src` into slot `slot` of batch_,
    // plus metadata into state_.
    void write_entry_(int slot, const V& v,
                      uint64_t id, StatusFlag flag, int32_t iter) {
        auto bv = batch_.view_capacity();  // full-capacity span view
        V dst = bv[std::size_t(slot)];
        copy_fields_(dst, v);
        // Metadata: state_ tracks logical size; but to keep slot-indexed
        // writes we use the capacity-length underlying arrays directly.
        state_.write_slot(slot, id, flag, iter);
    }

    void read_entry_(int slot, V& out_v,
                     uint64_t& out_id, StatusFlag& out_flag, int32_t& out_iter) {
        auto bv = batch_.view_capacity();
        V src = bv[std::size_t(slot)];
        copy_fields_(out_v, src);
        state_.read_slot(slot, out_id, out_flag, out_iter);
    }

    static void copy_fields_(V& dst, const V& src) {
        auto dt = dst.to_tuple();
        auto st = src.to_tuple();
        detail::for_each_field<V>([&](auto Ic) {
            constexpr std::size_t k = Ic;
            auto& d = std::get<k>(dt);
            const auto& s = std::get<k>(st);
            const std::size_t n = std::min<std::size_t>(d.size(), s.size());
            for (std::size_t i = 0; i < n; ++i) d[i] = s[i];
        });
    }

    void grow_(int new_capacity) {
        // Linearize into new storage starting at slot 0.
        Batch<V>   nb(batch_.N(), new_capacity, batch_.dmax());
        BatchState ns(new_capacity);
        int old_size = size_;
        for (int i = 0; i < old_size; ++i) {
            int src_slot = (front_ + i) % capacity();
            V src = batch_.view_capacity()[std::size_t(src_slot)];
            nb.push_back(src);
            uint64_t id; StatusFlag f; int32_t it;
            state_.read_slot(src_slot, id, f, it);
            ns.push_back(id, f, it, -1);
        }
        batch_ = std::move(nb);
        state_ = std::move(ns);
        front_ = old_size == 0 ? -1 : 0;
        back_  = old_size == 0 ? -1 : old_size - 1;
    }

    Batch<V>   batch_{};
    BatchState state_{};
    int        size_  = 0;
    int        front_ = -1;
    int        back_  = -1;
};

// ---------------------------------------------------------------------------
// Transfer helpers.  Predicate signature: bool(const StatusFlag&, int32_t).
// ---------------------------------------------------------------------------

// Batch -> Queue: copy every src entry satisfying `predicate` into the
// queue. If `consumed_status != 0`, OR that flag into the source entry.
template <class V, class Pred = MatchAll>
    requires batchable_view<V>
int queue_push(BatchQueue<V>& dst_queue,
               Batch<V>& src_batch,
               BatchState& src_state,
               Pred predicate = Pred{},
               StatusEnum consumed_status = StatusEnum(0)) {
    assert(src_batch.size() == src_state.size());
    auto sv = src_batch.view();
    auto st = src_state.view();
    int pushed = 0;
    for (int i = 0; i < sv.size(); ++i) {
        if (!predicate(st.status[i], st.iteration[i])) continue;
        dst_queue.push_back(sv[std::size_t(i)],
                            st.id[i],
                            st.status[i],
                            st.iteration[i]);
        if (int(consumed_status) != 0)
            st.status[i] = StatusFlag(int(st.status[i]) | int(consumed_status));
        ++pushed;
    }
    return pushed;
}

// Queue -> Batch: pop-from-front every queue entry satisfying `predicate`
// and append it to `dst_batch`/`dst_state`, until dst is full or the
// queue front no longer matches. If `consumed_status != 0`, OR that
// flag into the newly-written dst entry.
template <class V, class Pred = MatchAll>
    requires batchable_view<V>
int queue_push(Batch<V>& dst_batch,
               BatchState& dst_state,
               BatchQueue<V>& src_queue,
               Pred predicate = Pred{},
               StatusEnum consumed_status = StatusEnum(0)) {
    assert(dst_batch.size() == dst_state.size());
    int transferred = 0;
    while (!src_queue.empty() && dst_batch.size() < dst_batch.capacity()) {
        // Peek at logical-front.
        int slot = src_queue.front_index();
        StatusFlag f; int32_t it; uint64_t id;
        src_queue.state().read_slot(slot, id, f, it);
        if (!predicate(f, it)) break;
        // Pop and append.
        V src = src_queue.storage().view_capacity()[std::size_t(slot)];
        dst_batch.push_back(src);
        if (int(consumed_status) != 0)
            f = StatusFlag(int(f) | int(consumed_status));
        dst_state.push_back(id, f, it, -1);
        src_queue.discard_front();
        ++transferred;
    }
    return transferred;
}

} // namespace batch
