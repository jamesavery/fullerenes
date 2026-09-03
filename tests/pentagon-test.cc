// Tests for the FullereneDualView pentagon list.
//
// The list is a TYPE INVARIANT: every non-empty FullereneDualView carries the
// 12 degree-5 vertices, strictly ascending, in a 12-slot span.  Producers
// establish it (derive_pentagons -- all pentagon logic lives on the VIEW; the
// FullereneDual owner is a thin allocation wrapper), surgery maintains it
// (substitute_pentagons, atomic), and pentagons_valid() is the boundary gate.
// Covered here: derivation (both failure directions), the gate (each leg
// isolated), every guard of the substitution with byte-unchanged state after
// each throw, the k = 0 and k = 12 extremes, ownership (copy, move, and the
// N == 0 escape), both enumeration queues, staleness after relabelling, the
// batch contract's first constant-size field, and one coupled
// adjacency-surgery + substitution round trip.

#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

#include <fullerenes/triangulation.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/isomerdb.hh>

using Code = pentagon_error::Code;
using pentagon_buf = pentagon_storage_t<FullereneDualView>;
static constexpr int n_pent = FullereneDualView::n_pentagons;

// dual_slot hands out the STRONGEST view the destination can back
// (dual_slot_t names the law; these pin the instances).
static_assert(std::is_same_v<BuckyGen::dual_slot_t<FullereneDual>,  FullereneDualView>);
static_assert(std::is_same_v<BuckyGen::dual_slot_t<Triangulation>,  TriangulationView>);
static_assert(std::is_same_v<BuckyGen::dual_slot_t<Graph>,          TriangulationView>);

// C20's dual: 12 vertices, all pentagons.
static FullereneDual C20_dual() {
    return FullereneDual(20, std::vector<int>{0,1,2,3,4,5,6,7,8,9,10,11});
}

// One buckygen dual of C_n, filled through the strong enumeration path.
static FullereneDual buckygen_dual(int n) {
    auto Q = BuckyGen::start(n, false);
    FullereneDual g;
    EXPECT_TRUE(BuckyGen::next_fullerene(Q, g));
    BuckyGen::stop(Q);
    return g;
}

static std::vector<node_t> pentagon_list(const FullereneDualView& D) {
    return {D.pentagons.begin(), D.pentagons.end()};
}

static std::vector<node_t> vertices_of_degree(const FullereneDualView& D, int d) {
    std::vector<node_t> vs;
    for (node_t u = 0; u < D.N; u++)
        if (D.degree(u) == d) vs.push_back(u);
    return vs;
}

// The atomicity + code oracle: f() must throw pentagon_error with `code`,
// and leave the pentagon list byte-identical.
template<class F>
static void expect_pentagon_error(const FullereneDualView& D, Code code, F&& f) {
    const auto before = pentagon_list(D);
    try {
        f();
        FAIL() << "expected pentagon_error, got no throw";
    } catch (const pentagon_error& e) {
        EXPECT_EQ(e.code, code) << e.what();
    }
    EXPECT_EQ(before, pentagon_list(D));
}

// ---------------------------------------------------------------------------
// Derivation and the validity gate
// ---------------------------------------------------------------------------

TEST(Pentagons, DeriveEstablishesAscendingList) {
    FullereneDual g = C20_dual();
    ASSERT_EQ(g.N, 12);
    EXPECT_TRUE(g.pentagons_valid());
    for (int i = 0; i < n_pent; i++)
        EXPECT_EQ(g.pentagons[i], i);   // C20: every dual vertex is a pentagon
}

TEST(Pentagons, ValidityGateRejectsTampering) {
    FullereneDual g = C20_dual();
    std::swap(g.pentagons[3], g.pentagons[7]);   // right set, wrong order
    EXPECT_FALSE(g.pentagons_valid());
    std::swap(g.pentagons[3], g.pentagons[7]);
    EXPECT_TRUE(g.pentagons_valid());
    g.pentagons[0] = 11;                          // wrong multiset
    EXPECT_FALSE(g.pentagons_valid());
    g.derive_pentagons();
    EXPECT_TRUE(g.pentagons_valid());
}

TEST(Pentagons, ValidityGateSizeCheckDecidesAlone) {
    // A correctly-derived buffer re-viewed through a short span: only the
    // size leg of the gate can produce the false.
    FullereneDual g = C20_dual();
    pentagon_buf buf;
    FullereneDualView full(TriangulationView(g), buf);
    full.derive_pentagons();
    ASSERT_TRUE(full.pentagons_valid());
    FullereneDualView cut(TriangulationView(g), std::span<node_t>(buf).first(7));
    EXPECT_FALSE(cut.pentagons_valid());
}

TEST(Pentagons, DeriveThrowsOnTooFewPentagons) {
    // A tetrahedron: a perfectly good triangulation, zero degree-5 vertices.
    std::vector<node_t> nbrs = {1,2,3,-1,-1,-1, 0,3,2,-1,-1,-1,
                                0,1,3,-1,-1,-1, 0,2,1,-1,-1,-1};
    std::vector<uint8_t> deg = {3,3,3,3};
    pentagon_buf buf; buf.fill(-1);
    FullereneDualView D(TriangulationView(4, 6, nbrs, deg), buf);
    expect_pentagon_error(D, Code::NotFullereneDual,
                          [&]{ D.derive_pentagons(); });

    // The Owned establishment path refuses the same graph at construction.
    auto make_bad = [&]{ return FullereneDual(TriangulationView(4, 6, nbrs, deg)); };
    EXPECT_THROW(make_bad(), pentagon_error);
}

TEST(Pentagons, DeriveThrowsOnTooManyPentagons) {
    // Thirteen degree-5 rows (degrees are all derive consults): the fill must
    // refuse BEFORE writing, not overrun the 12-slot buffer.
    std::vector<node_t> nbrs(13*6, -1);
    std::vector<uint8_t> deg(13, 5);
    pentagon_buf buf; buf.fill(-1);
    FullereneDualView D(TriangulationView(13, 6, nbrs, deg), buf);
    expect_pentagon_error(D, Code::NotFullereneDual,
                          [&]{ D.derive_pentagons(); });
}

TEST(Pentagons, DeriveThrowsOnMisshapenSpan) {
    FullereneDual g = C20_dual();
    pentagon_buf buf; buf.fill(-1);
    FullereneDualView D(TriangulationView(g), std::span<node_t>(buf).first(7));
    expect_pentagon_error(D, Code::PentagonSpanWrongSize,
                          [&]{ D.derive_pentagons(); });
    EXPECT_FALSE(D.pentagons_valid());
}

// ---------------------------------------------------------------------------
// Atomic substitution -- on a dual that has hexagons (C24: 12 + 2)
// ---------------------------------------------------------------------------

struct C24 : public ::testing::Test {
    FullereneDual g;
    std::vector<node_t> hexagons;
    void SetUp() override {
        g = buckygen_dual(24);
        ASSERT_EQ(g.N, 14);
        ASSERT_TRUE(g.pentagons_valid());
        hexagons = vertices_of_degree(g, 6);
        ASSERT_EQ(hexagons.size(), 2u);
    }
};

TEST_F(C24, SubstituteIsSortedSetAlgebra) {
    const auto original = pentagon_list(g);
    const std::array<node_t,2> dep{original[0], original[5]};
    const std::array<node_t,2> arr{hexagons[1], hexagons[0]};  // unsorted input

    g.substitute_pentagons(dep, arr);
    auto expect = original;
    expect.erase(std::find(expect.begin(), expect.end(), dep[0]));
    expect.erase(std::find(expect.begin(), expect.end(), dep[1]));
    expect.insert(expect.end(), {hexagons[0], hexagons[1]});
    std::sort(expect.begin(), expect.end());
    EXPECT_EQ(pentagon_list(g), expect);
    // List-level only: the degrees now disagree, and the boundary gate says so.
    EXPECT_FALSE(g.pentagons_valid());

    // The inverse substitution restores the invariant exactly.
    g.substitute_pentagons(arr, dep);
    EXPECT_EQ(pentagon_list(g), original);
    EXPECT_TRUE(g.pentagons_valid());
}

TEST_F(C24, IdentitySubstitution) {
    const auto original = pentagon_list(g);
    g.substitute_pentagons({}, {});
    EXPECT_EQ(pentagon_list(g), original);
    EXPECT_TRUE(g.pentagons_valid());
}

TEST_F(C24, SubstituteGuards) {
    const node_t p0 = g.pentagons[0], p1 = g.pentagons[1];
    const node_t h0 = hexagons[0],    h1 = hexagons[1];

    expect_pentagon_error(g, Code::SubstitutionUnbalanced, [&]{
        const std::array<node_t,1> d{p0};
        g.substitute_pentagons(d, {}); });
    expect_pentagon_error(g, Code::PentagonAbsent, [&]{
        const std::array<node_t,1> d{h0}, a{h1};       // departure not a pentagon
        g.substitute_pentagons(d, a); });
    expect_pentagon_error(g, Code::DeparturesNotDistinct, [&]{
        const std::array<node_t,2> d{p0, p0}, a{h0, h1};
        g.substitute_pentagons(d, a); });
    expect_pentagon_error(g, Code::PentagonPresent, [&]{
        const std::array<node_t,1> d{p0}, a{p1};       // arrival already listed
        g.substitute_pentagons(d, a); });
    expect_pentagon_error(g, Code::ArrivalsNotDistinct, [&]{
        const std::array<node_t,2> d{p0, p1}, a{h0, h0};
        g.substitute_pentagons(d, a); });
    expect_pentagon_error(g, Code::PentagonVertexOutOfRange, [&]{
        const std::array<node_t,1> d{p0}, a{node_t(g.N)};
        g.substitute_pentagons(d, a); });
    EXPECT_TRUE(g.pentagons_valid());
}

TEST_F(C24, SubstituteRejectsOversizedRequest) {
    // k > n_pentagons must refuse before the fixed scratch arrays are touched.
    std::array<node_t, n_pent + 1> d{}, a{};
    expect_pentagon_error(g, Code::SubstitutionTooLong,
                          [&]{ g.substitute_pentagons(d, a); });
}

TEST_F(C24, InvariantViolationFailsLoudly) {
    // Sortedness is the producers' documented @inv, deliberately NOT
    // re-scanned per call (the orientation discipline).  What the op promises
    // on a violated list is no silent memory corruption: an argument guard
    // misfires honestly, or the O(1) conservation tripwire refuses before
    // the bounded merge.  Pinned on one concrete tampered shape.
    std::reverse(g.pentagons.begin(), g.pentagons.end());
    const std::array<node_t,1> d{g.pentagons[0]}, a{hexagons[0]};
    const auto tampered = pentagon_list(g);
    EXPECT_THROW(g.substitute_pentagons(d, a), pentagon_error);
    EXPECT_EQ(tampered, pentagon_list(g));   // pre-write refusal: unchanged

    // Restored, the very same substitution goes through.
    std::reverse(g.pentagons.begin(), g.pentagons.end());
    ASSERT_TRUE(g.pentagons_valid());
    g.substitute_pentagons(d, a);
    EXPECT_FALSE(g.pentagons_valid());       // list-level only, as designed
}

TEST(Pentagons, ReplaceAllPentagons) {
    // k = n_pentagons, the maximal case: every slot of the scratch arrays in
    // use.  C60's dual (32 vertices) has 20 hexagons to draw arrivals from.
    FullereneDual g = buckygen_dual(60);
    ASSERT_EQ(g.N, 32);
    const auto original = pentagon_list(g);
    const auto hexes    = vertices_of_degree(g, 6);
    ASSERT_GE((int)hexes.size(), n_pent);

    const std::span<const node_t> arr(hexes.data(), n_pent);
    g.substitute_pentagons(original, arr);
    std::vector<node_t> expect(arr.begin(), arr.end());
    std::sort(expect.begin(), expect.end());
    EXPECT_EQ(pentagon_list(g), expect);

    g.substitute_pentagons(arr, original);
    EXPECT_EQ(pentagon_list(g), original);
    EXPECT_TRUE(g.pentagons_valid());
}

// ---------------------------------------------------------------------------
// Ownership: the thin wrapper's whole contract
// ---------------------------------------------------------------------------

TEST(Pentagons, OwnedCopyEstablishesOverOwnStorage) {
    FullereneDual a = C20_dual();
    FullereneDual b(a);                          // FullereneDual copy
    EXPECT_TRUE(b.pentagons_valid());
    EXPECT_NE(a.pentagons.data(), b.pentagons.data());
    EXPECT_EQ(b.pentagons.data(), b.owned_pentagons.data());

    FullereneDual c(static_cast<const GraphView&>(a));   // foreign-view copy
    EXPECT_TRUE(c.pentagons_valid());
}

TEST(Pentagons, MoveCarriesListAndClearsSource) {
    FullereneDual a = C20_dual();
    FullereneDual b(std::move(a));
    EXPECT_TRUE(b.pentagons_valid());
    EXPECT_EQ(b.pentagons.data(), b.owned_pentagons.data());
    EXPECT_TRUE(a.pentagons.empty());

    FullereneDual c;
    c = std::move(b);
    EXPECT_TRUE(c.pentagons_valid());
    EXPECT_EQ(c.pentagons.data(), c.owned_pentagons.data());
    EXPECT_TRUE(b.pentagons.empty());
}

TEST(Pentagons, EmptyAndAllocatingConstruction) {
    // N == 0: nothing to establish, nothing thrown.
    FullereneDual empty;
    EXPECT_EQ(empty.N, 0);
    FullereneDual copied{GraphView{}};
    EXPECT_EQ(copied.N, 0);

    // Allocate-then-fill: the list is construction-stale, and the gate is the
    // thing that says so.
    FullereneDual scratch(12);
    EXPECT_EQ(scratch.N, 12);
    EXPECT_FALSE(scratch.pentagons_valid());
}

TEST(Pentagons, RelabellingStalesTheList) {
    // Owned's graph-restructuring helpers know nothing of pentagons: a
    // relabelling stales the list, and the producer re-derives at its
    // boundary.  (C24, so the reversal genuinely moves the degree-5 set.)
    FullereneDual g = buckygen_dual(24);
    Permutation pi(g.N);
    for (int u = 0; u < g.N; u++) pi[u] = g.N - 1 - u;
    g.apply_permutation(pi);
    EXPECT_FALSE(g.pentagons_valid());
    g.derive_pentagons();
    EXPECT_TRUE(g.pentagons_valid());
}

// ---------------------------------------------------------------------------
// The enumeration paths: strong and weak destinations, both queue kinds
// ---------------------------------------------------------------------------

TEST(Pentagons, EnumerationEstablishesInvariant) {
    for (int N : {20, 24, 26, 28, 30}) {
        auto Q = BuckyGen::start(N, false);
        FullereneDual g;
        int64_t count = 0;
        while (BuckyGen::next_fullerene(Q, g)) {
            ASSERT_TRUE(g.pentagons_valid())
                << "C" << N << " isomer " << count << " has an invalid pentagon list";
            count++;
        }
        EXPECT_EQ(count, IsomerDB::number_isomers(N, "Any", false)) << "C" << N;
    }
}

TEST(Pentagons, HerdEnumerationEstablishesInvariant) {
    // The parallel production path must establish the same invariant.
    BuckyGen::buckyherd_queue HQ(28, /*Nchunks=*/1, /*Nworkers=*/1,
                                 /*IPR=*/false, /*only_nontrivial=*/false);
    FullereneDual g;
    int64_t count = 0;
    while (HQ.next_fullerene(g)) {
        ASSERT_TRUE(g.pentagons_valid()) << "herd isomer " << count;
        count++;
    }
    EXPECT_EQ(count, IsomerDB::number_isomers(28, "Any", false));
}

TEST(Pentagons, WeakDestinationStillFills) {
    // A pentagon-less owner is filled as a plain triangulation -- same graphs,
    // no faked invariant.
    auto Q = BuckyGen::start(20, false);
    Triangulation T;
    ASSERT_TRUE(BuckyGen::next_fullerene(Q, T));
    BuckyGen::stop(Q);
    ASSERT_EQ(T.N, 12);
    for (node_t u = 0; u < T.N; u++) EXPECT_EQ(T.degree(u), 5);
}

// ---------------------------------------------------------------------------
// Batch: the contract's first constant-size field
// ---------------------------------------------------------------------------

TEST(Pentagons, BatchCarriesPentagons) {
    FullereneDual g = C20_dual();

    batch::Batch<FullereneDualView> B(g.N, 2, g.dmax);
    B.push_back(g);
    B.push_back(g);
    ASSERT_EQ(B.size(), 2);

    for (int i = 0; i < B.size(); i++) {
        FullereneDualView v = B[i];
        ASSERT_EQ((int)v.pentagons.size(), n_pent);
        EXPECT_TRUE(v.pentagons_valid());
        EXPECT_EQ(pentagon_list(v), pentagon_list(g));
    }
    // Entries are backed by disjoint slab slices, not shared storage.
    EXPECT_NE(B[0].pentagons.data(), B[1].pentagons.data());

    // A slice re-views the same slab entries, pentagons included.
    auto S = B.slice(1, 1);
    EXPECT_EQ(S[0].pentagons.data(), B[1].pentagons.data());
    EXPECT_EQ(pentagon_list(S[0]), pentagon_list(g));
}

// ---------------------------------------------------------------------------
// The loop closed: adjacency surgery and list surgery meet at the gate
// ---------------------------------------------------------------------------

TEST(Pentagons, SurgeryAndSubstitutionAgree) {
    // Degree-level round trip (this is a unit test of the degree <-> list
    // coupling, not of a planar-embedding pipeline): demote two adjacent
    // hexagons by removing their shared edge, promote two non-adjacent
    // pentagons by inserting one, then perform the matching list
    // substitution -- the boundary gate must close.
    FullereneDual g = buckygen_dual(60);
    const auto pents = pentagon_list(g);
    const auto hexes = vertices_of_degree(g, 6);

    node_t h1 = -1, h2 = -1;
    for (node_t h : hexes) {
        for (node_t v : g[h])
            if (g.degree(v) == 6) { h1 = h; h2 = v; break; }
        if (h1 >= 0) break;
    }
    ASSERT_GE(h1, 0);

    node_t p1 = -1, p2 = -1;
    for (size_t i = 0; i < pents.size() && p1 < 0; i++)
        for (size_t j = i+1; j < pents.size(); j++)
            if (g.find(pents[i], pents[j]) < 0) { p1 = pents[i]; p2 = pents[j]; break; }
    ASSERT_GE(p1, 0);

    g.remove_edge({h1, h2});        // 6,6 -> 5,5
    g.insert_edge({p1, p2});        // 5,5 -> 6,6
    EXPECT_FALSE(g.pentagons_valid());

    const std::array<node_t,2> dep{p1, p2}, arr{h1, h2};
    g.substitute_pentagons(dep, arr);
    EXPECT_TRUE(g.pentagons_valid());
}
