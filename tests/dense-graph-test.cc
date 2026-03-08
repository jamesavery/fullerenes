// Test suite for Spanify::DenseGraph and freeze/thaw conversions.

#include <gtest/gtest.h>
#include <fullerenes/dense_graph.hh>

#include <fullerenes/fullerenegraph.hh>
#include <fullerenes/isomerdb.hh>

namespace S = Spanify;

// Helper: get a fullerene graph for testing.
FullereneGraph get_fg(int N) {
    if (N == 20) return FullereneGraph::C20();
    bool IPR = (N >= 60);
    IsomerDB db = IsomerDB::readPDB(N, IPR);
    return IsomerDB::makeIsomer(N, db.entries[0]);
}

// ---------------------------------------------------------------------------
// DenseGraph: verify Graph::neighbours (which IS a DenseGraph with dmax=10) is correct
// ---------------------------------------------------------------------------

class DenseFromGraph : public ::testing::TestWithParam<int> {};

TEST_P(DenseFromGraph, CubicGraphCorrect) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);

    // FG.neighbours is now DenseGraph with dmax=3 (cubic restride)
    EXPECT_EQ(FG.size(), N);
    for (int v = 0; v < N; ++v) {
        EXPECT_EQ(FG.degree(v), 3);
        auto nbs = FG.nbrs(v);
        EXPECT_EQ(int(nbs.size()), 3);
    }
}

TEST_P(DenseFromGraph, DualGraphCorrect) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);
    PlanarGraph dual = FG.dual_graph();

    EXPECT_EQ(dual.size(), N / 2 + 2);
    for (int v = 0; v < dual.N; ++v) {
        int d = dual.degree(v);
        EXPECT_GE(d, 5);
        EXPECT_LE(d, 6);
    }
}

INSTANTIATE_TEST_SUITE_P(Sizes, DenseFromGraph, ::testing::Values(20, 60, 80));

// ---------------------------------------------------------------------------
// DenseGraph mutation operations
// ---------------------------------------------------------------------------

TEST(DenseMutation, PushBack) {
    S::OwnedDenseGraph<> g(4, 3);

    g.push_back(0, 1);
    g.push_back(0, 2);
    g.push_back(0, 3);

    EXPECT_EQ(g.degree(0), 3);
    auto nbs = g.nbrs(0);
    EXPECT_EQ(nbs[0], 1);
    EXPECT_EQ(nbs[1], 2);
    EXPECT_EQ(nbs[2], 3);
}

TEST(DenseMutation, InsertAt) {
    S::OwnedDenseGraph<> g(1, 4);

    g.push_back(0, 1);
    g.push_back(0, 3);
    g.insert_at(0, 2, 1);  // Insert 2 between 1 and 3

    EXPECT_EQ(g.degree(0), 3);
    auto nbs = g.nbrs(0);
    EXPECT_EQ(nbs[0], 1);
    EXPECT_EQ(nbs[1], 2);
    EXPECT_EQ(nbs[2], 3);
}

TEST(DenseMutation, EraseAt) {
    S::OwnedDenseGraph<> g(1, 4);

    g.push_back(0, 1);
    g.push_back(0, 2);
    g.push_back(0, 3);
    g.erase_at(0, 1);  // Remove middle element (2)

    EXPECT_EQ(g.degree(0), 2);
    auto nbs = g.nbrs(0);
    EXPECT_EQ(nbs[0], 1);
    EXPECT_EQ(nbs[1], 3);
}

TEST(DenseMutation, Find) {
    S::OwnedDenseGraph<> g(1, 3);

    g.push_back(0, 5);
    g.push_back(0, 10);
    g.push_back(0, 15);

    EXPECT_EQ(g.find(0, 5), 0);
    EXPECT_EQ(g.find(0, 10), 1);
    EXPECT_EQ(g.find(0, 15), 2);
    EXPECT_EQ(g.find(0, 99), -1);
}

TEST(DenseMutation, SpanAndAssign) {
    S::OwnedDenseGraph<> g(2, 3);

    // push_back through graph method
    g.push_back(0, 1);
    g.push_back(0, 2);
    g.push_back(0, 3);
    EXPECT_EQ(g[0].size(), 3u);
    EXPECT_EQ(g[0][0], 1);
    EXPECT_EQ(g[0][1], 2);

    // assign_row with initializer_list
    g.assign_row(1, {10, 20, 30});
    EXPECT_EQ(g[1].size(), 3u);
    EXPECT_EQ(g[1][2], 30);

    // span to vector conversion
    std::vector<int> v(g[0].begin(), g[0].end());
    EXPECT_EQ(v.size(), 3u);
    EXPECT_EQ(v[0], 1);

    // clear_row
    g.clear_row(0);
    EXPECT_EQ(g[0].size(), 0u);
    EXPECT_EQ(g.degree(0), 0);

    // assign_row with span
    g.assign_row(0, g[1]);
    EXPECT_EQ(g[0].size(), 3u);
    EXPECT_EQ(g[0][0], 10);
}

// ---------------------------------------------------------------------------
// freeze: Dense -> CSR
// ---------------------------------------------------------------------------

class FreezeTest : public ::testing::TestWithParam<int> {};

TEST_P(FreezeTest, CubicGraphRoundTrip) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);

    // FG.neighbours is DenseGraph with dmax=3; freeze it to CSR
    auto csr = S::freeze(static_cast<const neighbours_t&>(FG));

    EXPECT_EQ(csr.N, N);
    EXPECT_EQ(csr.n_arcs, 3 * N);
    EXPECT_TRUE(S::validate(S::as_view(csr)));

    for (int v = 0; v < N; ++v) {
        auto nbs = S::neighbours(S::as_view(csr), v);
        auto old_nbs = FG.nbrs(v);
        ASSERT_EQ(int(nbs.size()), int(old_nbs.size()));
        for (int j = 0; j < int(nbs.size()); ++j)
            EXPECT_EQ(int(nbs[j]), old_nbs[j]);
    }
}

TEST_P(FreezeTest, DualGraphRoundTrip) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);
    PlanarGraph dual = FG.dual_graph();

    auto csr = S::freeze(static_cast<const neighbours_t&>(dual));

    EXPECT_EQ(csr.N, dual.N);
    EXPECT_EQ(csr.n_arcs, 3 * N);
    EXPECT_TRUE(S::validate(S::as_view(csr)));
}

INSTANTIATE_TEST_SUITE_P(Sizes, FreezeTest, ::testing::Values(20, 60, 80));

// ---------------------------------------------------------------------------
// thaw: CSR -> Dense
// ---------------------------------------------------------------------------

class ThawTest : public ::testing::TestWithParam<int> {};

TEST_P(ThawTest, CubicFreezeThawRoundTrip) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);

    auto csr = S::freeze(static_cast<const neighbours_t&>(FG));
    auto dense2 = S::thaw(csr, 3);

    EXPECT_EQ(int(dense2.N), N);
    for (int v = 0; v < N; ++v) {
        EXPECT_EQ(dense2.degree(v), FG.degree(v));
        auto nbs1 = FG.nbrs(v);
        auto nbs2 = dense2.nbrs(v);
        for (int j = 0; j < FG.degree(v); ++j)
            EXPECT_EQ(nbs2[j], nbs1[j]);
    }
}

TEST_P(ThawTest, DualFreezeThawRoundTrip) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);
    PlanarGraph dual = FG.dual_graph();

    auto csr = S::freeze(static_cast<const neighbours_t&>(dual));
    auto dense2 = S::thaw(csr, 6);

    EXPECT_EQ(int(dense2.N), dual.N);
    for (int v = 0; v < dual.N; ++v) {
        EXPECT_EQ(dense2.degree(v), dual.degree(v));
        auto nbs1 = dual.nbrs(v);
        auto nbs2 = dense2.nbrs(v);
        for (int j = 0; j < dual.degree(v); ++j)
            EXPECT_EQ(nbs2[j], nbs1[j]);
    }
}

INSTANTIATE_TEST_SUITE_P(Sizes, ThawTest, ::testing::Values(20, 60, 80));

// ---------------------------------------------------------------------------
// to_vectors round-trip
// ---------------------------------------------------------------------------

class ToVectorsRoundTrip : public ::testing::TestWithParam<int> {};

TEST_P(ToVectorsRoundTrip, CubicGraph) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);

    auto adj = FG.to_vectors();

    ASSERT_EQ(int(adj.size()), N);
    for (int v = 0; v < N; ++v) {
        auto nbs = FG.nbrs(v);
        ASSERT_EQ(int(adj[v].size()), int(nbs.size()));
        for (int j = 0; j < int(adj[v].size()); ++j)
            EXPECT_EQ(adj[v][j], nbs[j]);
    }
}

INSTANTIATE_TEST_SUITE_P(Sizes, ToVectorsRoundTrip, ::testing::Values(20, 60, 80));
