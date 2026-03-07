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
// DenseGraph from neighbours_t
// ---------------------------------------------------------------------------

class DenseFromNeighbours : public ::testing::TestWithParam<int> {};

TEST_P(DenseFromNeighbours, CubicGraphMatchesOriginal) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);

    auto dense = S::dense_from_neighbours<3>(FG.N, FG.neighbours);

    EXPECT_EQ(dense.Nv, N);
    for (int v = 0; v < N; ++v) {
        EXPECT_EQ(dense.degree(v), 3);
        auto nbs = dense.nbrs(v);
        auto old_nbs = FG.nbrs(v);
        ASSERT_EQ(int(nbs.size()), int(old_nbs.size()));
        for (int j = 0; j < int(nbs.size()); ++j)
            EXPECT_EQ(int(nbs[j]), old_nbs[j]);
    }
}

TEST_P(DenseFromNeighbours, DualGraphMatchesOriginal) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);
    PlanarGraph dual = FG.dual_graph();

    auto dense = S::dense_from_neighbours<6>(dual.N, dual.neighbours);

    EXPECT_EQ(dense.Nv, N / 2 + 2);
    for (int v = 0; v < dual.N; ++v) {
        auto nbs = dense.nbrs(v);
        auto old_nbs = dual.nbrs(v);
        ASSERT_EQ(int(nbs.size()), int(old_nbs.size()));
        for (int j = 0; j < int(nbs.size()); ++j)
            EXPECT_EQ(int(nbs[j]), old_nbs[j]);
    }
}

INSTANTIATE_TEST_SUITE_P(Sizes, DenseFromNeighbours, ::testing::Values(20, 60, 80));

// ---------------------------------------------------------------------------
// DenseGraph mutation operations
// ---------------------------------------------------------------------------

TEST(DenseMutation, PushBack) {
    S::DenseGraph<3> g(4);

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
    S::DenseGraph<4> g(1);

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
    S::DenseGraph<4> g(1);

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
    S::DenseGraph<3> g(1);

    g.push_back(0, 5);
    g.push_back(0, 10);
    g.push_back(0, 15);

    EXPECT_EQ(g.find(0, 5), 0);
    EXPECT_EQ(g.find(0, 10), 1);
    EXPECT_EQ(g.find(0, 15), 2);
    EXPECT_EQ(g.find(0, 99), -1);
}

// ---------------------------------------------------------------------------
// freeze: Dense -> CSR
// ---------------------------------------------------------------------------

class FreezeTest : public ::testing::TestWithParam<int> {};

TEST_P(FreezeTest, CubicGraphRoundTrip) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);

    auto dense = S::dense_from_neighbours<3>(FG.N, FG.neighbours);
    auto csr = S::freeze(dense);

    EXPECT_EQ(csr.N, N);
    EXPECT_EQ(csr.n_arcs, 3 * N);
    EXPECT_TRUE(S::validate(S::as_view(csr)));

    // Check CSR neighbours match original
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

    auto dense = S::dense_from_neighbours<6>(dual.N, dual.neighbours);
    auto csr = S::freeze(dense);

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

    // neighbours -> Dense -> CSR -> Dense
    auto dense1 = S::dense_from_neighbours<3>(FG.N, FG.neighbours);
    auto csr = S::freeze(dense1);
    auto dense2 = S::thaw<3>(csr);

    EXPECT_EQ(dense2.Nv, dense1.Nv);
    for (int v = 0; v < N; ++v) {
        EXPECT_EQ(dense2.degree(v), dense1.degree(v));
        auto nbs1 = dense1.nbrs(v);
        auto nbs2 = dense2.nbrs(v);
        for (int j = 0; j < dense1.degree(v); ++j)
            EXPECT_EQ(nbs2[j], nbs1[j]);
    }
}

TEST_P(ThawTest, DualFreezeThawRoundTrip) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);
    PlanarGraph dual = FG.dual_graph();

    auto dense1 = S::dense_from_neighbours<6>(dual.N, dual.neighbours);
    auto csr = S::freeze(dense1);
    auto dense2 = S::thaw<6>(csr);

    EXPECT_EQ(dense2.Nv, dense1.Nv);
    for (int v = 0; v < dual.N; ++v) {
        EXPECT_EQ(dense2.degree(v), dense1.degree(v));
        auto nbs1 = dense1.nbrs(v);
        auto nbs2 = dense2.nbrs(v);
        for (int j = 0; j < dense1.degree(v); ++j)
            EXPECT_EQ(nbs2[j], nbs1[j]);
    }
}

INSTANTIATE_TEST_SUITE_P(Sizes, ThawTest, ::testing::Values(20, 60, 80));

// ---------------------------------------------------------------------------
// neighbours_t round-trip: neighbours -> Dense -> neighbours
// ---------------------------------------------------------------------------

class NeighboursRoundTrip : public ::testing::TestWithParam<int> {};

TEST_P(NeighboursRoundTrip, CubicGraph) {
    int N = GetParam();
    FullereneGraph FG = get_fg(N);

    auto dense = S::dense_from_neighbours<3>(FG.N, FG.neighbours);
    auto adj = S::to_neighbours(dense);

    ASSERT_EQ(int(adj.size()), N);
    for (int v = 0; v < N; ++v) {
        ASSERT_EQ(adj[v].size(), FG.neighbours[v].size());
        for (int j = 0; j < int(adj[v].size()); ++j)
            EXPECT_EQ(adj[v][j], FG.neighbours[v][j]);
    }
}

INSTANTIATE_TEST_SUITE_P(Sizes, NeighboursRoundTrip, ::testing::Values(20, 60, 80));
