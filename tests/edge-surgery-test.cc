// Embedded-edge surgery: insert_edge_at / remove_edge_at / the slot cores,
// the GraphView vertex-pair wrappers, twin maintenance, and the row gap
// words (dense_graph.hh).
//
// Oracles, in order of strength:
//  1. The MAINTAINED twin must equal a from-scratch recompute -- maintenance
//     and re-derivation must be indistinguishable.  Two independent
//     recomputes are used: Owned::compute_twin (which now throws on
//     asymmetry in every build configuration) and a test-local reference.
//  2. twin_is_valid() -- itself proven to bite by one planted fault per
//     clause (each plant is detectable by exactly one clause, so deleting
//     any clause is caught), plus a last-slot plant against a truncated
//     scan.
//  3. A twin-free surgery sequence compared row-for-row against a naive
//     vector<vector> model.
// Every contract violation is asserted by its CODE, and through
// expect_surgery_error, which also verifies the throwing call left the
// graph byte-identical (the atomicity @post).

#include <gtest/gtest.h>
#include <random>

#include "fullerenes/graph.hh"
#include "fullerenes/owned.hh"

using G_t     = Owned<GraphView>;
using arcix_t = GraphView::arcix_t;
using Spanify::graph_surgery_error;
using Code    = Spanify::graph_surgery_error::Code;
static constexpr uint8_t no_slot = GraphView::no_slot;

// ---------------------------------------------------------------------------
// Builders and oracles
// ---------------------------------------------------------------------------

// Cyclic n-gon with room to grow: row of i = {i+1, i-1} (mod n).
static G_t make_cycle(int n, uint8_t dmax) {
    G_t G(n, dmax);
    for (node_t i = 0; i < n; ++i) {
        G.push_back(i, (i + 1) % n);
        G.push_back(i, (i + n - 1) % n);
    }
    return G;
}

// Maintained twin == Owned::compute_twin's recompute, entrywise over live arcs.
static bool twin_matches_recompute(const G_t& G) {
    G_t H(G);
    H.compute_twin();
    for (node_t u = 0; u < G.N; ++u)
        for (int i = 0; i < G.degree(u); ++i)
            if (G.twin[u * G.dmax + i] != H.twin[u * H.dmax + i]) return false;
    return true;
}

// The atomicity + code oracle: f() must throw graph_surgery_error with
// exactly `code`, leaving adjacency, degrees and twin byte-identical.
template<typename F>
static void expect_surgery_error(G_t& G, Code code, F&& f) {
    const G_t orig(G);
    try {
        f();
        FAIL() << "expected graph_surgery_error, got no throw";
    } catch (const graph_surgery_error& e) {
        EXPECT_EQ(int(e.code), int(code)) << e.what();
    }
    EXPECT_TRUE(G == orig) << "throwing call changed the adjacency";
    for (node_t u = 0; u < G.N; ++u) EXPECT_EQ(G.degree(u), orig.degree(u));
    if (G.has_twin())
        for (size_t a = 0; a < G.twin.size(); ++a)
            EXPECT_EQ(G.twin[a], orig.twin[a]) << "throwing call changed twin[" << a << "]";
}

// ---------------------------------------------------------------------------
// Typed suite: the arithmetic at both vertex-id widths
// ---------------------------------------------------------------------------

// Minimal owner for a bare RSRAdjacencyView<K> (Owned<> is node_t-only).
template<typename K>
struct TestRSR {
    int dmax;
    std::vector<K>       nbr;
    std::vector<uint8_t> dg, tw;
    Spanify::RSRAdjacencyView<K> G;
    TestRSR(int n, int dmax_, bool with_twin)
        : dmax(dmax_),
          nbr(size_t(n) * dmax_, K(-1)), dg(n, 0),
          tw(with_twin ? size_t(n) * dmax_ : 0, 0),
          G(K(n), dmax_, std::span<K>(nbr), std::span<uint8_t>(dg),
            with_twin ? std::span<uint8_t>(tw) : std::span<uint8_t>()) {}

    void cycle() {
        const int n = G.N;
        for (int i = 0; i < n; ++i) {
            G.push_back(K(i), K((i + 1) % n));
            G.push_back(K(i), K((i + n - 1) % n));
        }
    }
    // Test-local reference recompute -- independent of Owned::compute_twin.
    void recompute_twin() {
        for (int u = 0; u < int(G.N); ++u)
            for (int i = 0; i < G.degree(K(u)); ++i) {
                const int j = G.find(G.nbrs(K(u))[i], K(u));
                ASSERT_GE(j, 0) << "asymmetric adjacency in test builder";
                tw[size_t(u) * dmax + i] = uint8_t(j);
            }
    }
    bool twin_matches_fresh() const {
        for (int u = 0; u < int(G.N); ++u)
            for (int i = 0; i < G.degree(K(u)); ++i) {
                const int j = G.find(G.nbrs(K(u))[i], K(u));
                if (j < 0 || G.twin[size_t(u) * dmax + i] != uint8_t(j)) return false;
            }
        return true;
    }
};

template<typename K> class SurgeryWidth : public ::testing::Test {};
using Widths = ::testing::Types<int32_t, uint16_t>;
TYPED_TEST_SUITE(SurgeryWidth, Widths);

TYPED_TEST(SurgeryWidth, ComputeAndValidate) {
    TestRSR<TypeParam> T(6, 4, true);
    T.cycle();
    T.recompute_twin();
    EXPECT_TRUE(T.G.twin_is_valid());
    EXPECT_TRUE(T.twin_matches_fresh());
    TestRSR<TypeParam> F(6, 4, false);
    F.cycle();
    EXPECT_TRUE(F.G.twin_is_valid());          // vacuous without twin
}

TYPED_TEST(SurgeryWidth, InsertRemoveRoundTrip) {
    TestRSR<TypeParam> T(6, 4, true);
    T.cycle();
    T.recompute_twin();

    auto [a_uv, a_vu] = T.G.insert_edge_at({TypeParam(0), 0}, {TypeParam(3), 0});
    EXPECT_EQ(int(T.G.target(a_uv)), 3);
    EXPECT_TRUE(T.G.twin_is_valid());
    EXPECT_TRUE(T.twin_matches_fresh());

    const std::vector<TypeParam> nbr_before = T.nbr;
    auto [c_u, c_v] = T.G.remove_edge_at(a_uv);
    EXPECT_TRUE(T.G.twin_is_valid());
    EXPECT_TRUE(T.twin_matches_fresh());
    T.G.insert_edge_at(c_u, c_v);
    EXPECT_EQ(T.nbr, nbr_before);              // exact restoration (s > 0)
    EXPECT_TRUE(T.G.twin_is_valid());
    EXPECT_TRUE(T.twin_matches_fresh());
}

TYPED_TEST(SurgeryWidth, RandomSurgeryStress) {
    constexpr int N = 10, DMAX = 8, STEPS = 300;
    TestRSR<TypeParam> T(N, DMAX, true);
    T.cycle();
    T.recompute_twin();

    std::mt19937 rng(2026);
    std::uniform_int_distribution<int> vtx(0, N - 1);
    for (int step = 0; step < STEPS; ++step) {
        const TypeParam u = TypeParam(vtx(rng)), v = TypeParam(vtx(rng));
        if (u == v) continue;
        const int s = T.G.find(u, v);
        if (s >= 0) {
            T.G.remove_edge_at({u, uint8_t(s)});
        } else if (T.G.degree(u) > 0 && T.G.degree(u) < DMAX &&
                   T.G.degree(v) > 0 && T.G.degree(v) < DMAX) {
            std::uniform_int_distribution<int> su(0, T.G.degree(u) - 1),
                                               sv(0, T.G.degree(v) - 1);
            T.G.insert_edge_at({u, uint8_t(su(rng))}, {v, uint8_t(sv(rng))});
        }
        ASSERT_TRUE(T.G.twin_is_valid()) << "step " << step;
        if (step % 25 == 0)
            ASSERT_TRUE(T.twin_matches_fresh()) << "step " << step;
    }
    EXPECT_TRUE(T.twin_matches_fresh());
}

// ---------------------------------------------------------------------------
// Corner semantics
// ---------------------------------------------------------------------------

TEST(EdgeSurgery, InsertEdgeAtMaintainsTwin) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();

    // Chord 0-3, landing after arc 0->1 (slot 0) and after arc 3->4 (slot 0).
    auto [a_uv, a_vu] = G.insert_edge_at({0, 0}, {3, 0});
    EXPECT_EQ(G.target(a_uv), 3);
    EXPECT_EQ(G.target(a_vu), 0);
    EXPECT_EQ(int(a_uv.second), 1);            // after slot 0
    EXPECT_EQ(int(a_vu.second), 1);

    EXPECT_EQ(G.degree(0), 3);
    EXPECT_EQ(G.nbrs(0)[1], 3);                // row 0: {1, 3, 5}
    EXPECT_EQ(G.nbrs(3)[1], 0);                // row 3: {4, 0, 2}

    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));
}

TEST(EdgeSurgery, RemoveReturnsVacatedCornersAndRoundTrips) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    G.insert_edge_at({0, 0}, {3, 0});
    const G_t orig(G);

    // Remove the chord: arc 0->3 sits at slot 1.
    auto [c_u, c_v] = G.remove_edge_at({0, 1});
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));
    EXPECT_EQ(c_u, (arcix_t{0, 0}));           // arc preceding the gap
    EXPECT_EQ(c_v, (arcix_t{3, 0}));

    // Re-inserting at the vacated corners restores the graph exactly.
    G.insert_edge_at(c_u, c_v);
    EXPECT_TRUE(G == orig);
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));
}

TEST(EdgeSurgery, RemoveSlotZeroWrapsCorner) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    G.insert_edge_at({0, 0}, {3, 0});          // row 0: {1, 3, 5}

    // Remove arc 0->1 at slot 0: predecessor in cyclic order was the last
    // slot, which shifts down -- corner (0, d_new-1) = (0, 1).
    auto [c_u, c_v] = G.remove_edge_at({0, 0});
    EXPECT_EQ(c_u, (arcix_t{0, 1}));
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));

    // Re-insert at the corners: same cyclic order (possibly rotated rows).
    G.insert_edge_at(c_u, c_v);
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));
    // Cyclic check: in 0's rotation the successor of 1 is 3 again.
    arcix_t a01 = G.find_arc(0, 1);
    EXPECT_EQ(G.target(G.next(a01)), 3);
}

TEST(EdgeSurgery, EmptiedRowHasNoCorner) {
    G_t G(2, 2);
    G.push_back(0, 1);
    G.push_back(1, 0);
    G.compute_twin();

    auto [c_u, c_v] = G.remove_edge_at({0, 0});
    EXPECT_EQ(c_u.second, no_slot);
    EXPECT_EQ(c_v.second, no_slot);
    EXPECT_EQ(G.degree(0), 0);
    EXPECT_EQ(G.degree(1), 0);
    EXPECT_TRUE(G.twin_is_valid());            // vacuously: no live arcs
}

// One endpoint emptied, the other keeps arcs: the returned pair is mixed.
TEST(EdgeSurgery, MixedCornerOnDegreeOneEndpoint) {
    G_t G(3, 4);                               // path 0-1-2
    G.push_back(0, 1);
    G.push_back(1, 0); G.push_back(1, 2);
    G.push_back(2, 1);
    G.compute_twin();

    auto [c_u, c_v] = G.remove_edge_at({0, 0});
    EXPECT_EQ(c_u, (arcix_t{0, no_slot}));     // 0 is now isolated
    EXPECT_EQ(c_v, (arcix_t{1, 0}));           // 1 keeps arc 1->2 at slot 0
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));
}

TEST(EdgeSurgery, DmaxBoundaryInsertSucceeds) {
    G_t G = make_cycle(6, 3);                  // deg 2 everywhere, one slot spare
    G.compute_twin();
    G.insert_edge_at({0, 0}, {3, 0});          // fills both rows to dmax
    EXPECT_EQ(G.degree(0), 3);
    EXPECT_EQ(G.degree(3), 3);
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));
}

// ---------------------------------------------------------------------------
// The validator's teeth: one planted fault per clause
// ---------------------------------------------------------------------------

// A raw row edit (surgery minus maintenance) leaves a stale twin the
// validator must detect -- the @post of the construction-phase primitives.
TEST(EdgeSurgery, ValidatorBitesOnRawRowEdit) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    ASSERT_TRUE(G.twin_is_valid());
    G.insert_at(0, 3, 1);                      // one-sided, twin-oblivious
    EXPECT_FALSE(G.twin_is_valid());
}

// Clause 1: the entry names a dead slot of the target row.
TEST(EdgeSurgery, ValidatorBitesOnDeadSlot) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    G.twin[0 * G.dmax + 0] = uint8_t(G.degree(1));   // dead but in storage
    EXPECT_FALSE(G.twin_is_valid());
}

// Clause 2: the entry names a live slot of the target row that does not
// point back.
TEST(EdgeSurgery, ValidatorBitesOnWrongReverse) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    // Arc 0->1 (slot 0); row of 1 is {2, 0}: slot 0 points at 2, not 0.
    G.twin[0 * G.dmax + 0] = 0;
    EXPECT_FALSE(G.twin_is_valid());
}

// Clause 3 in isolation needs a state where every entry points back at the
// right vertex but the pairing is not an involution -- possible only with
// parallel arcs, built through the construction primitives.
TEST(EdgeSurgery, ValidatorBitesOnBrokenInvolution) {
    G_t G(2, 4);
    G.push_back(0, 1); G.push_back(0, 1);      // two parallel arcs 0->1
    G.push_back(1, 0); G.push_back(1, 0);
    G.compute_twin();
    // Cross the pairing: 0/0 <-> 1/0 and 0/1 <-> 1/1 become a 4-cycle.
    G.twin[0 * G.dmax + 0] = 0;  G.twin[1 * G.dmax + 0] = 1;
    G.twin[0 * G.dmax + 1] = 1;  G.twin[1 * G.dmax + 1] = 0;
    // Clauses 1-2 hold at every arc (both slots of the other row name the
    // right vertex); only the involution clause can refuse this.
    EXPECT_FALSE(G.twin_is_valid());
}

// The scan must reach every slot, not just slot 0.
TEST(EdgeSurgery, ValidatorBitesAtLastSlot) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    const int last = G.degree(0) - 1;
    G.twin[0 * G.dmax + last] = uint8_t(G.degree(G.nbrs(0)[last]));
    EXPECT_FALSE(G.twin_is_valid());
}

// The validator must survive what it validates: garbage rows report
// invalid, never crash (a live prefix over the K(-1) padding).
TEST(EdgeSurgery, ValidatorSurvivesGarbageRows) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    G.resize_row(0, 4);                        // live prefix now covers padding
    EXPECT_FALSE(G.twin_is_valid());
}

// ---------------------------------------------------------------------------
// Contract violations: every guard, by code, leaving the graph unchanged
// ---------------------------------------------------------------------------

TEST(EdgeSurgery, GuardSelfLoopInsert) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    expect_surgery_error(G, Code::SelfLoop, [&]{ G.insert_edge_slots(2, 0, 2, 0); });
    expect_surgery_error(G, Code::SelfLoop, [&]{ G.insert_edge({2, 2}); });
}

// A pre-existing self-loop (constructible through the row primitives) must
// be REFUSED by removal, not silently double-closed into a corrupted row.
TEST(EdgeSurgery, GuardSelfLoopRemove) {
    G_t G(3, 4);
    G.push_back(0, 0);                         // a one-slot loop
    G.push_back(0, 1); G.push_back(1, 0);
    expect_surgery_error(G, Code::SelfLoop, [&]{ G.remove_edge({0, 0}); });
    EXPECT_TRUE(G.edge_exists({0, 1}));        // the unrelated edge survived
}

TEST(EdgeSurgery, GuardRowFullEitherSide) {
    // v full, u not: the v half of the room check must fire on its own.
    G_t H(4, 2);                               // path 0-1-2 at dmax=2: 1 is full
    H.push_back(0, 1); H.push_back(1, 0);
    H.push_back(1, 2); H.push_back(2, 1);
    expect_surgery_error(H, Code::RowFull, [&]{ H.insert_edge({3, 1}); });
    // both full:
    G_t F = make_cycle(6, 2);
    expect_surgery_error(F, Code::RowFull, [&]{ F.insert_edge({0, 3}); });
    expect_surgery_error(F, Code::RowFull, [&]{ F.insert_edge_at({0, 0}, {3, 0}); });
}

TEST(EdgeSurgery, GuardSlotOutOfRange) {
    G_t G = make_cycle(6, 4);
    expect_surgery_error(G, Code::SlotOutOfRange,
                         [&]{ G.insert_edge_slots(0, G.degree(0) + 1, 3, 0); });
}

TEST(EdgeSurgery, GuardVertexOutOfRange) {
    G_t G = make_cycle(6, 4);
    expect_surgery_error(G, Code::VertexOutOfRange, [&]{ G.insert_edge({0, 6}); });
    expect_surgery_error(G, Code::VertexOutOfRange, [&]{ G.remove_edge({0, 6}); });
    expect_surgery_error(G, Code::VertexOutOfRange, [&]{ G.insert_edge_hint({6, 0}); });
}

TEST(EdgeSurgery, GuardInsertOntoHalfEdge) {
    G_t G = make_cycle(6, 4);
    G.push_back(0, 3);                         // half an edge: corruption
    expect_surgery_error(G, Code::AsymmetricAdjacency,
                         [&]{ G.insert_edge_slots(0, 0, 3, 0); });
}

TEST(EdgeSurgery, GuardDeadArcEitherArgument) {
    G_t G = make_cycle(6, 4);
    const uint8_t dead_u = uint8_t(G.degree(0)), dead_v = uint8_t(G.degree(3));
    expect_surgery_error(G, Code::NotLiveArc, [&]{ G.insert_edge_at({0, dead_u}, {3, 0}); });
    expect_surgery_error(G, Code::NotLiveArc, [&]{ G.insert_edge_at({0, 0}, {3, dead_v}); });
    expect_surgery_error(G, Code::NotLiveArc, [&]{ G.remove_edge_at({0, dead_u}); });
}

TEST(EdgeSurgery, GuardWrongSlotsInRange) {
    G_t G = make_cycle(6, 4);
    // Slots in range, but slot 0 of row 0 targets 1, not 3.
    expect_surgery_error(G, Code::NotMutualPair,
                         [&]{ G.remove_edge_slots(0, 0, 3, 0); });
}

TEST(EdgeSurgery, GuardEdgeAlreadyPresent) {
    G_t G = make_cycle(6, 4);
    expect_surgery_error(G, Code::EdgePresent, [&]{ G.insert_edge_at({0, 1}, {1, 1}); });
}

TEST(EdgeSurgery, GuardAbsentSuccessorThrows) {
    G_t G = make_cycle(6, 4);
    // 4 is not a neighbour of 0 -- the caller's picture of the rotation is
    // wrong, and the old silent-append behaviour is exactly what must NOT
    // happen.
    expect_surgery_error(G, Code::SuccessorNotNeighbour,
                         [&]{ G.insert_edge({0, 3}, /*suc_uv=*/4, /*suc_vu=*/-1); });
}

TEST(EdgeSurgery, GuardReverseArcNeedsTwin) {
    G_t G = make_cycle(6, 4);
    expect_surgery_error(G, Code::TwinAbsent, [&]{ (void)G.reverse_arc({0, 0}); });
}

TEST(EdgeSurgery, ComputeTwinThrowsOnAsymmetry) {
    G_t G(2, 2);
    G.push_back(0, 1);                         // no reverse
    EXPECT_THROW(G.compute_twin(), graph_surgery_error);
}

// ---------------------------------------------------------------------------
// The vertex-pair wrappers
// ---------------------------------------------------------------------------

TEST(EdgeSurgery, WrapperEnsureIsIdempotent) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    const G_t orig(G);
    EXPECT_TRUE(G.insert_edge({0, 1}));        // present -> true, no-op
    EXPECT_TRUE(G == orig);
    EXPECT_TRUE(G.insert_edge_hint({0, 1}, 5, 2));
    EXPECT_TRUE(G == orig);
    EXPECT_TRUE(G.twin_is_valid());
}

// The successful removal path, twin-free and twinned, with the rows checked.
TEST(EdgeSurgery, WrapperRemoveSucceeds) {
    G_t G = make_cycle(6, 4);
    EXPECT_TRUE(G.remove_edge({0, 1}));
    EXPECT_EQ(G.degree(0), 1);
    EXPECT_EQ(G.nbrs(0)[0], 5);
    EXPECT_EQ(G.degree(1), 1);
    EXPECT_EQ(G.nbrs(1)[0], 2);

    G_t H = make_cycle(6, 4);
    H.compute_twin();
    EXPECT_TRUE(H.remove_edge({0, 1}));
    EXPECT_TRUE(H.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(H));
}

// Positioning through named successors, BOTH sides asserted.
TEST(EdgeSurgery, WrapperPositionsBothSides) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    // Insert 3 before 5 in 0's row, and 0 before 2 in 3's row.
    EXPECT_FALSE(G.insert_edge({0, 3}, /*suc_uv=*/5, /*suc_vu=*/2));
    EXPECT_EQ(G.nbrs(0)[1], 3);                // row 0: {1, 3, 5}
    EXPECT_EQ(G.nbrs(3)[1], 0);                // row 3: {4, 0, 2}
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));
}

// Append (-1) means the LAST slot of each row, not the first.
TEST(EdgeSurgery, WrapperAppendLandsLast) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();
    EXPECT_FALSE(G.insert_edge({0, 3}));
    EXPECT_EQ(G.nbrs(0).back(), 3);
    EXPECT_EQ(G.nbrs(3).back(), 0);
    EXPECT_TRUE(twin_matches_recompute(G));
}

TEST(EdgeSurgery, HintPositionsWhenPresentAppendsWhenNot) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();

    // u-side hint present (5 at slot 1 of row 0), v-side hint present (2 at
    // slot 1 of row 3): both position.
    EXPECT_FALSE(G.insert_edge_hint({0, 3}, /*hint_uv=*/5, /*hint_vu=*/2));
    EXPECT_EQ(G.nbrs(0)[1], 3);
    EXPECT_EQ(G.nbrs(3)[1], 0);
    EXPECT_TRUE(twin_matches_recompute(G));

    // Hint 4 is NOT a neighbour of 1: hint semantics append (the windup
    // idiom) where insert_edge's assertion semantics would throw.
    EXPECT_FALSE(G.insert_edge_hint({1, 4}, /*hint_uv=*/4, /*hint_vu=*/-1));
    EXPECT_EQ(G.nbrs(1).back(), 4);
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));
}

TEST(EdgeSurgery, WrapperAppendsAndStartsIsolatedRows) {
    G_t G(3, 4);
    G.compute_twin();                          // empty but present
    EXPECT_FALSE(G.insert_edge({0, 1}));       // both rows empty: slot 0
    EXPECT_FALSE(G.insert_edge({1, 2}));
    EXPECT_FALSE(G.insert_edge({2, 0}));
    EXPECT_EQ(G.degree(0), 2);
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));
}

TEST(EdgeSurgery, WrapperRemoveAbsentAndAsymmetric) {
    G_t G = make_cycle(6, 4);
    EXPECT_FALSE(G.remove_edge({0, 3}));       // absent: false, no throw

    G.push_back(0, 3);                         // half an edge: corruption
    expect_surgery_error(G, Code::NotMutualPair, [&]{ G.remove_edge({0, 3}); });
}

// ---------------------------------------------------------------------------
// Attribute mirroring and twin-free surgery
// ---------------------------------------------------------------------------

TEST(EdgeSurgery, AttributeShiftMirrorsSurgery) {
    G_t G = make_cycle(6, 4);
    G.compute_twin();

    // Tag every arc of vertex 0 with its target; the attribute row must
    // keep matching the neighbour row across surgery.
    std::vector<int> attr(G.dmax, -1);
    for (int i = 0; i < G.degree(0); ++i) attr[i] = G.nbrs(0)[i];

    int old_deg = G.degree(0);
    auto [a_uv, a_vu] = G.insert_edge_at({0, 0}, {3, 0});
    Spanify::row_open_gap(std::span<int>(attr), a_uv.second, old_deg);
    attr[a_uv.second] = G.target(a_uv);
    for (int i = 0; i < G.degree(0); ++i) EXPECT_EQ(attr[i], G.nbrs(0)[i]);

    old_deg = G.degree(0);
    G.remove_edge_at(a_uv);
    Spanify::row_close_gap(std::span<int>(attr), a_uv.second, old_deg);
    for (int i = 0; i < G.degree(0); ++i) EXPECT_EQ(attr[i], G.nbrs(0)[i]);
}

// Twin-free surgery (remove_edge_at's find fallback, both cores' twin-free
// branches) against a naive vector<vector> reference model.
TEST(EdgeSurgery, TwinFreeSurgeryMatchesModel) {
    constexpr int N = 8, DMAX = 6, STEPS = 200;
    G_t G = make_cycle(N, DMAX);
    std::vector<std::vector<node_t>> model(N);
    for (node_t u = 0; u < N; ++u)
        model[u] = {node_t((u + 1) % N), node_t((u + N - 1) % N)};

    std::mt19937 rng(7);
    std::uniform_int_distribution<int> vtx(0, N - 1);
    for (int step = 0; step < STEPS; ++step) {
        const node_t u = vtx(rng), v = vtx(rng);
        if (u == v) continue;
        const int s = G.find(u, v);
        if (s >= 0) {
            G.remove_edge_at({u, uint8_t(s)});
            model[u].erase(model[u].begin() + s);
            const auto rv = std::find(model[v].begin(), model[v].end(), u);
            model[v].erase(rv);
        } else if (G.degree(u) > 0 && G.degree(u) < DMAX &&
                   G.degree(v) > 0 && G.degree(v) < DMAX) {
            std::uniform_int_distribution<int> su(0, G.degree(u) - 1),
                                               sv(0, G.degree(v) - 1);
            const int cu = su(rng), cv = sv(rng);
            G.insert_edge_at({u, uint8_t(cu)}, {v, uint8_t(cv)});
            model[u].insert(model[u].begin() + cu + 1, v);
            model[v].insert(model[v].begin() + cv + 1, u);
        }
        for (node_t w = 0; w < N; ++w) {
            ASSERT_EQ(size_t(G.degree(w)), model[w].size()) << "step " << step;
            for (int i = 0; i < G.degree(w); ++i)
                ASSERT_EQ(G.nbrs(w)[i], model[w][i])
                    << "step " << step << " row " << int(w) << " slot " << i;
        }
    }
}

// ---------------------------------------------------------------------------
// Owner hygiene: a computed twin survives owner reallocation
// ---------------------------------------------------------------------------

TEST(EdgeSurgery, OwnedResizeKeepsTwinShape) {
    G_t G = make_cycle(4, 4);
    G.compute_twin();
    G.resize(8);                               // twin must grow alongside
    EXPECT_EQ(G.twin.size(), G.neighbours.size());
    EXPECT_TRUE(G.twin_is_valid());            // new rows have no live arcs
    EXPECT_FALSE(G.insert_edge({6, 7}));       // surgery into the new region
    EXPECT_FALSE(G.insert_edge({3, 6}));       // and across old/new
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_TRUE(twin_matches_recompute(G));
}

TEST(EdgeSurgery, OwnedPushPopKeepTwinShape) {
    G_t G = make_cycle(4, 4);
    G.compute_twin();
    G.push_back(std::vector<node_t>{});        // a new isolated vertex row
    EXPECT_EQ(G.twin.size(), G.neighbours.size());
    EXPECT_TRUE(G.twin_is_valid());
    EXPECT_FALSE(G.insert_edge({4, 0}));
    EXPECT_TRUE(twin_matches_recompute(G));
    G.remove_edge({4, 0});
    G.pop_back();
    EXPECT_EQ(G.twin.size(), G.neighbours.size());
    EXPECT_TRUE(G.twin_is_valid());
}

// ---------------------------------------------------------------------------
// Random surgery through every public entry point
// ---------------------------------------------------------------------------

TEST(EdgeSurgery, RandomSurgeryStress) {
    constexpr int N = 10, DMAX = 8, STEPS = 500;
    G_t G = make_cycle(N, DMAX);
    G.compute_twin();

    std::mt19937 rng(42);
    std::uniform_int_distribution<int> vtx(0, N - 1), pick(0, 3);

    for (int step = 0; step < STEPS; ++step) {
        const node_t u = vtx(rng), v = vtx(rng);
        if (u == v) continue;

        if (G.edge_exists({u, v})) {
            if (pick(rng) == 0) EXPECT_TRUE(G.remove_edge({u, v}));
            else                G.remove_edge_at(G.find_arc(u, v));
        } else if (G.degree(u) < DMAX && G.degree(v) < DMAX) {
            switch (G.degree(u) == 0 || G.degree(v) == 0 ? 1 : pick(rng)) {
              case 0: {                        // corner-addressed
                std::uniform_int_distribution<int> su(0, G.degree(u) - 1),
                                                   sv(0, G.degree(v) - 1);
                G.insert_edge_at({u, uint8_t(su(rng))}, {v, uint8_t(sv(rng))});
                break;
              }
              case 1: G.insert_edge({u, v}); break;              // append
              case 2:                          // position before the first
                G.insert_edge({u, v}, G.nbrs(u).empty() ? -1 : G.nbrs(u)[0], -1);
                break;
              default:                         // hint, possibly absent
                G.insert_edge_hint({u, v}, node_t(vtx(rng)), node_t(vtx(rng)));
            }
        }

        ASSERT_TRUE(G.twin_is_valid()) << "step " << step
            << " (u=" << int(u) << ", v=" << int(v) << ")";
        if (step % 25 == 0)
            ASSERT_TRUE(twin_matches_recompute(G)) << "step " << step;
    }
    EXPECT_TRUE(twin_matches_recompute(G));
}
