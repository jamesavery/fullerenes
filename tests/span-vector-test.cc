// Test suite for Spanify::SpanStack<T>, Spanify::BitSpan, and
// Spanify::SpanVector<T> (span_vector.hh).
//
// The load-bearing property is vector equivalence: a SpanStack driven by any
// op sequence within capacity must be indistinguishable from a
// std::vector<int> driven by the same sequence, and a BitSpan from a
// std::vector<bool>.  Overflow must be loud (sticky flag, rejected push,
// storage untouched -- checked via a guard slot just past the span), never
// silent truncation or growth.  BitSpan's storage layout (word i/64, bit
// i%64) is pinned against the raw word image, not just self-consistency.

// Included FIRST: this test is the header's self-containment witness.
#include <fullerenes/span_vector.hh>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

using Spanify::SpanStack;
using Spanify::BitSpan;
using Spanify::SpanVector;

namespace {

constexpr int kGuard = 0x5EED;

// Backing storage with one guard slot PAST the span: a write outside the
// live region -- e.g. a push that stores before checking capacity -- lands
// in the guard slot and is caught by guard_intact().
struct StackRig {
    std::vector<int> storage;
    SpanStack<int>   s;
    explicit StackRig(long cap)
        : storage((std::size_t)cap + 1, kGuard),
          s(std::span<int>(storage).first((std::size_t)cap)) {}
    bool guard_intact() const { return storage.back() == kGuard; }
};

long popcount(const BitSpan& b) {
    long c = 0;
    for (long i = 0; i < b.capacity(); i++) c += b.test(i);
    return c;
}

}  // namespace

// ---------------------------------------------------------------------------
// SpanStack
// ---------------------------------------------------------------------------

TEST(SpanStack, FreshStateIsEmpty) {
    StackRig r(8);   // storage pre-filled with garbage: must be ignored
    EXPECT_TRUE(r.s.empty());
    EXPECT_EQ(r.s.size(), 0);
    EXPECT_EQ(r.s.capacity(), 8);
    EXPECT_FALSE(r.s.overflowed());
    EXPECT_EQ(r.s.peak(), 0);
    EXPECT_TRUE(r.s.live().empty());
}

TEST(SpanStack, LifoOrderAndStorageImage) {
    StackRig r(4);
    ASSERT_TRUE(r.s.push_back(10));
    ASSERT_TRUE(r.s.push_back(20));
    ASSERT_TRUE(r.s.push_back(30));
    // Layout pin: the live prefix occupies storage[0..n) bottom-to-top.
    EXPECT_EQ(r.storage[0], 10);
    EXPECT_EQ(r.storage[1], 20);
    EXPECT_EQ(r.storage[2], 30);
    EXPECT_EQ(r.s.back(), 30); r.s.pop_back();
    EXPECT_EQ(r.s.back(), 20); r.s.pop_back();
    ASSERT_TRUE(r.s.push_back(40));
    EXPECT_EQ(r.s.back(), 40); r.s.pop_back();
    EXPECT_EQ(r.s.back(), 10); r.s.pop_back();
    EXPECT_TRUE(r.s.empty());
    EXPECT_TRUE(r.guard_intact());
}

namespace {

// Mirrored random op sequence at one capacity: SpanStack vs std::vector<int>.
// Pushes are attempted UNCONDITIONALLY so the at-capacity rejection, the
// latch, and the mirrored peak are part of the asserted contract on every
// path (small caps make the boundary a common event; cap 0 and 1 cover the
// degenerate stacks).
void check_lifo_equivalence(long cap, uint32_t seed, long ops) {
    StackRig r(cap);
    std::vector<int> ref;
    long ref_peak = 0;
    bool expect_latch = false;

    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> op(0, 9), val(-1000, 1000);

    for (long i = 0; i < ops; i++) {
        const int o = op(rng);
        if (o < 6) {                                   // 60%: attempt a push
            const int  x    = val(rng);
            const bool room = (long)ref.size() < cap;
            ASSERT_EQ(r.s.push_back(x), room);         // the contract itself
            if (room) {
                ref.push_back(x);
                ref_peak = std::max(ref_peak, (long)ref.size());
            } else {
                expect_latch = true;
            }
        } else if (o < 9) {                            // 30%: pop when non-empty
            if (!ref.empty()) {
                ASSERT_EQ(r.s.back(), ref.back());
                r.s.pop_back();
                ref.pop_back();
            }
        } else {                                       // 10%: clear (fresh state)
            r.s.clear();
            ref.clear();
            expect_latch = false;                      // clear resets the latch
        }
        ASSERT_EQ(r.s.size(), (long)ref.size());
        ASSERT_EQ(r.s.overflowed(), expect_latch);
        ASSERT_EQ(r.s.peak(), ref_peak);               // lifetime: never reset
        const std::span<const int> lv = r.s.live();
        ASSERT_TRUE(std::equal(lv.begin(), lv.end(), ref.begin(), ref.end()));
    }
    EXPECT_TRUE(r.guard_intact());
}

}  // namespace

TEST(SpanStack, LifoEquivalence) {
    for (long cap : {0L, 1L, 4L, 8L, 64L})
        check_lifo_equivalence(cap, 0xf00d, 20000);
}

TEST(SpanStack, ExactlyFullIsNotOverflowed) {
    StackRig r(4);
    for (int i = 0; i < 4; i++) ASSERT_TRUE(r.s.push_back(i));
    EXPECT_EQ(r.s.size(), r.s.capacity());
    EXPECT_FALSE(r.s.overflowed());                    // full is NOT overflowed
    EXPECT_EQ(r.s.peak(), 4);

    EXPECT_FALSE(r.s.push_back(99));                   // rejected, not truncated
    EXPECT_TRUE(r.s.overflowed());                     // now latched
    EXPECT_EQ(r.s.peak(), 4);                          // rejection never bumps peak
    EXPECT_EQ(r.s.size(), 4);
    for (int i = 0; i < 4; i++) EXPECT_EQ(r.storage[i], i);  // contents untouched
    EXPECT_TRUE(r.guard_intact());                     // a rejected push writes NOTHING
}

TEST(SpanStack, OverflowLatchIsSticky) {
    StackRig r(4);
    for (int i = 0; i < 4; i++) ASSERT_TRUE(r.s.push_back(i));
    EXPECT_FALSE(r.s.push_back(99));
    EXPECT_TRUE(r.s.overflowed());

    r.s.pop_back();
    EXPECT_TRUE(r.s.overflowed());                     // still latched after pop
    ASSERT_TRUE(r.s.push_back(5));                     // room again: push succeeds
    EXPECT_TRUE(r.s.overflowed());                     // latch survives success

    r.s.clear();
    EXPECT_FALSE(r.s.overflowed());                    // clear resets the latch
    EXPECT_TRUE(r.s.empty());
    EXPECT_TRUE(r.guard_intact());
}

TEST(SpanStack, PeakIsLifetimeHighWater) {
    StackRig r(8);
    for (int i = 0; i < 5; i++) r.s.push_back(i);
    EXPECT_EQ(r.s.peak(), 5);
    r.s.pop_back(); r.s.pop_back(); r.s.pop_back();
    EXPECT_EQ(r.s.peak(), 5);                          // not reduced by pops
    r.s.push_back(9);
    EXPECT_EQ(r.s.peak(), 5);                          // below old peak: unchanged
    r.s.clear();
    EXPECT_EQ(r.s.peak(), 5);                          // survives clear (lifetime)
    for (int i = 0; i < 7; i++) r.s.push_back(i);
    EXPECT_EQ(r.s.peak(), 7);                          // new high water
}

TEST(SpanStack, PopAndClearLeaveStorageUntouched) {
    StackRig r(4);
    r.s.push_back(7);
    r.s.push_back(8);
    r.s.pop_back();
    EXPECT_EQ(r.storage[1], 8);                        // pop leaves the slot alone
    r.s.clear();
    EXPECT_EQ(r.storage[0], 7);                        // clear leaves storage alone
    EXPECT_EQ(r.storage[1], 8);
}

TEST(SpanStack, DefaultConstructedIsBoundedAndLoud) {
    SpanStack<int> s;                                  // no storage bound at all
    EXPECT_EQ(s.capacity(), 0);
    EXPECT_TRUE(s.empty());
    EXPECT_FALSE(s.push_back(1));                      // must reject, not write
    EXPECT_TRUE(s.overflowed());
    EXPECT_EQ(s.peak(), 0);
}

TEST(SpanStack, CopiesShareStorageAndForkCounters) {
    StackRig r(4);
    r.s.push_back(7);
    SpanStack<int> c = r.s;                            // TC copy: same storage
    ASSERT_TRUE(c.push_back(8));
    EXPECT_EQ(r.storage[1], 8);                        // wrote the SHARED storage
    EXPECT_EQ(c.size(), 2);
    EXPECT_EQ(r.s.size(), 1);                          // ...but counters are per-copy
    EXPECT_EQ(r.s.back(), 7);
}

// ---------------------------------------------------------------------------
// BitSpan
// ---------------------------------------------------------------------------

TEST(BitSpan, ClearAllEstablishesFreshState) {
    std::vector<std::uint64_t> words(BitSpan::words_for(200), ~0ull);  // garbage
    BitSpan b{std::span<std::uint64_t>(words)};
    EXPECT_EQ(b.capacity(), 64 * BitSpan::words_for(200));  // 200 rounds up to 256
    b.clear_all();
    EXPECT_EQ(popcount(b), 0);
}

// Pins the bit <-> (word, bit) map against the raw word image, so a
// self-consistent relayout (bit- or word-reversal) cannot pass: the layout
// is part of the contract (device transfers compare raw words).
TEST(BitSpan, StorageLayoutIsWordThenBit) {
    std::vector<std::uint64_t> words(3, 0);
    BitSpan b{std::span<std::uint64_t>(words)};
    b.set(0);    EXPECT_EQ(words[0], 1ull);
    b.set(63);   EXPECT_EQ(words[0], 1ull | (1ull << 63));
    b.set(64);   EXPECT_EQ(words[1], 1ull);
    b.set(191);  EXPECT_EQ(words[2], 1ull << 63);
    b.clear(63); EXPECT_EQ(words[0], 1ull);
    EXPECT_EQ(words[1], 1ull);                         // clear touches ONE word
}

TEST(BitSpan, SetTestClearAcrossWordBoundaries) {
    std::vector<std::uint64_t> words(3, 0);
    BitSpan b{std::span<std::uint64_t>(words)};
    const long probes[] = {0, 1, 62, 63, 64, 65, 127, 128, 191};

    for (long i : probes) {
        b.set(i);
        EXPECT_TRUE(b.test(i));
    }
    EXPECT_EQ(popcount(b), (long)std::size(probes));   // ONLY the probed bits

    b.clear(64);
    EXPECT_FALSE(b.test(64));
    EXPECT_TRUE(b.test(63));                           // neighbours intact
    EXPECT_TRUE(b.test(65));
    b.clear(0);
    EXPECT_FALSE(b.test(0));
    EXPECT_TRUE(b.test(1));

    b.clear_all();
    EXPECT_EQ(popcount(b), 0);
}

TEST(BitSpan, TestAndSetReturnsPreviousThenSets) {
    std::vector<std::uint64_t> words(2, 0);
    BitSpan b{std::span<std::uint64_t>(words)};
    EXPECT_FALSE(b.test_and_set(5));                   // first claim: was clear
    EXPECT_TRUE(b.test(5));                            // ...and is now set
    EXPECT_TRUE(b.test_and_set(5));                    // second claim: was set
    EXPECT_TRUE(b.test(5));                           // ...and stays set
    EXPECT_FALSE(b.test_and_set(64));                  // cross-word claim
    EXPECT_EQ(words[1], 1ull);                         // raw image agrees
    EXPECT_EQ(words[0], 1ull << 5);
}

// Mirrored random op sequence vs std::vector<bool>, comparing EVERY
// addressable bit each iteration -- bits in [kBits, capacity()) are never
// indexed by an op and must stay false throughout.
TEST(BitSpan, VectorBoolEquivalence) {
    constexpr long kBits  = 200;                       // non-multiple of 64
    const     long kWords = BitSpan::words_for(kBits);
    std::vector<std::uint64_t> words(kWords, 0);
    BitSpan b{std::span<std::uint64_t>(words)};
    std::vector<bool> ref(kBits, false);

    std::mt19937 rng(0xbeef);
    std::uniform_int_distribution<int>  op(0, 9);
    std::uniform_int_distribution<long> idx(0, kBits - 1);

    for (long i = 0; i < 20000; i++) {
        const int  o = op(rng);
        const long k = idx(rng);
        if (o < 4)      { b.set(k);   ref[k] = true;  }
        else if (o < 7) { b.clear(k); ref[k] = false; }
        else if (o < 9) {                              // claim: previous value mirrored
            ASSERT_EQ(b.test_and_set(k), (bool)ref[k]);
            ref[k] = true;
        }
        else            { b.clear_all(); ref.assign(kBits, false); }

        for (long j = 0; j < kBits; j++) ASSERT_EQ(b.test(j), ref[j]);
        for (long j = kBits; j < b.capacity(); j++) ASSERT_FALSE(b.test(j));
    }
}

TEST(BitSpan, DefaultConstructedHasNoBits) {
    BitSpan b;
    EXPECT_EQ(b.capacity(), 0);
    b.clear_all();                                     // no-op, no write
}

// ---------------------------------------------------------------------------
// SpanVector (the pre-existing type: its hand-written Rule of 5 must repoint
// the view on copy/move -- the one subtle behaviour in this header)
// ---------------------------------------------------------------------------

TEST(SpanVector, CopyOfOwningRepointsToOwnStorage) {
    SpanVector<int> a(std::vector<int>{1, 2, 3});
    SpanVector<int> b(a);
    ASSERT_TRUE(b.owns_memory());
    ASSERT_NE(a.data(), b.data());                     // deep copy, repointed
    b[0] = 9;
    EXPECT_EQ(a[0], 1);                                // originals untouched
    EXPECT_EQ(b[0], 9);
}

TEST(SpanVector, MoveStealsStorageAndEmptiesSource) {
    SpanVector<int> a(std::vector<int>{1, 2, 3});
    const int* p = a.data();
    SpanVector<int> b(std::move(a));
    EXPECT_EQ(b.data(), p);                            // storage stolen, view repointed
    EXPECT_EQ(b[1], 2);
    EXPECT_TRUE(a.empty());                            // moved-from view cleared
}

TEST(SpanVector, ViewModeSharesExternalStorage) {
    std::vector<int> ext{4, 5};
    SpanVector<int> v{std::span<int>(ext)};
    SpanVector<int> w(v);                              // copy of a view: still a view
    ASSERT_FALSE(w.owns_memory());
    w[0] = 6;
    EXPECT_EQ(ext[0], 6);                              // shared external storage
}
