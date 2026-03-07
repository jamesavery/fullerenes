// Test suite for Spanify::PlanarCSR — validates against existing Graph classes.
//
// For each test graph (C20, C60 Ih, C80 IPR#1), we:
//   1. Build the old Graph via the existing API
//   2. Convert to PlanarCSR via from_neighbours()
//   3. Verify structural invariants (validate())
//   4. Check degree/neighbours match the old graph
//   5. Check next/prev/next_on_face match the old graph
//   6. Check face computation matches old compute_faces_oriented()
//   7. Also test the dual graph (variable degree 5-6) via dual_graph()

#include <gtest/gtest.h>
#include <fullerenes/planar_csr.hh>

#include <fullerenes/fullerenegraph.hh>
#include <fullerenes/isomerdb.hh>

namespace S = Spanify;  // Short alias — no 'using namespace' to avoid collision with ::Graph

// Helper: build a CSRGraph from an existing old-style Graph object
S::CSRGraph csr_from_graph(const ::Graph& G) {
    return S::from_neighbours(G.N, G.neighbours);
}

// Helper: get a fullerene graph for testing. Uses readPDB (never readBinary/getIsomer).
FullereneGraph get_fullerene(int N) {
    if (N == 20) return FullereneGraph::C20();
    bool IPR = (N >= 60);
    IsomerDB db = IsomerDB::readPDB(N, IPR);
    return IsomerDB::makeIsomer(N, db.entries[0]);
}

// ---------------------------------------------------------------------------
// Structural invariants
// ---------------------------------------------------------------------------

class CSRStructure : public ::testing::TestWithParam<int> {};

TEST_P(CSRStructure, ValidateFromCubicGraph) {
    int N = GetParam();

    FullereneGraph FG = get_fullerene(N);

    auto csr = csr_from_graph(FG);

    EXPECT_EQ(csr.N, N);
    EXPECT_EQ(csr.n_arcs, 3 * N);
    EXPECT_EQ(int(csr.offsets.size()), N + 1);
    EXPECT_EQ(int(csr.values.size()), 3 * N);
    EXPECT_EQ(int(csr.twin.size()), 3 * N);

    EXPECT_TRUE(S::validate(S::as_view(csr)));

    for (int v = 0; v < N; ++v)
        EXPECT_EQ(S::degree(S::as_view(csr), v), 3);
}

TEST_P(CSRStructure, ValidateFromDualGraph) {
    int N = GetParam();

    FullereneGraph FG = get_fullerene(N);
    PlanarGraph dual = FG.dual_graph();
    auto csr = csr_from_graph(dual);

    EXPECT_EQ(csr.N, N / 2 + 2);
    EXPECT_EQ(csr.n_arcs, 3 * N);
    EXPECT_TRUE(S::validate(S::as_view(csr)));

    auto view = S::as_view(csr);
    for (int v = 0; v < csr.N; ++v) {
        int d = S::degree(view, v);
        EXPECT_GE(d, 5);
        EXPECT_LE(d, 6);
    }
}

INSTANTIATE_TEST_SUITE_P(Sizes, CSRStructure, ::testing::Values(20, 60, 80));

// ---------------------------------------------------------------------------
// Navigation: compare CSR next/prev with old Graph
// ---------------------------------------------------------------------------

class CSRNavigation : public ::testing::TestWithParam<int> {};

TEST_P(CSRNavigation, NeighboursMatch) {
    int N = GetParam();
    FullereneGraph FG = get_fullerene(N);

    auto csr = csr_from_graph(FG);
    auto view = S::as_view(csr);

    for (int v = 0; v < N; ++v) {
        auto nbs = S::neighbours(view, v);
        auto old_nbs = FG.nbrs(v);
        ASSERT_EQ(int(nbs.size()), int(old_nbs.size()));
        for (int j = 0; j < int(nbs.size()); ++j)
            EXPECT_EQ(int(nbs[j]), old_nbs[j]);
    }
}

TEST_P(CSRNavigation, NextPrevMatch) {
    int N = GetParam();
    FullereneGraph FG = get_fullerene(N);

    auto csr = csr_from_graph(FG);
    auto view = S::as_view(csr);

    for (int u = 0; u < N; ++u) {
        for (int v : FG.nbrs(u)) {
            EXPECT_EQ(int(S::next(view, u, v)), FG.next(u, v))
                << "next(" << u << ", " << v << ") mismatch";
            EXPECT_EQ(int(S::prev(view, u, v)), FG.prev(u, v))
                << "prev(" << u << ", " << v << ") mismatch";
        }
    }
}

TEST_P(CSRNavigation, NextOnFaceMatch) {
    int N = GetParam();
    FullereneGraph FG = get_fullerene(N);

    auto csr = csr_from_graph(FG);
    auto view = S::as_view(csr);

    for (int u = 0; u < N; ++u) {
        for (int v : FG.nbrs(u)) {
            EXPECT_EQ(int(S::next_on_face(view, u, v)), FG.next_on_face(u, v))
                << "next_on_face(" << u << ", " << v << ") mismatch";
        }
    }
}

TEST_P(CSRNavigation, DualNeighboursMatch) {
    int N = GetParam();
    FullereneGraph FG = get_fullerene(N);

    PlanarGraph dual = FG.dual_graph();
    auto csr = csr_from_graph(dual);
    auto view = S::as_view(csr);

    for (int v = 0; v < dual.N; ++v) {
        auto nbs = S::neighbours(view, v);
        auto old_nbs = dual.nbrs(v);
        ASSERT_EQ(int(nbs.size()), int(old_nbs.size()));
        for (int j = 0; j < int(nbs.size()); ++j)
            EXPECT_EQ(int(nbs[j]), old_nbs[j]);
    }
}

TEST_P(CSRNavigation, DualNextPrevMatch) {
    int N = GetParam();
    FullereneGraph FG = get_fullerene(N);

    PlanarGraph dual = FG.dual_graph();
    auto csr = csr_from_graph(dual);
    auto view = S::as_view(csr);

    for (int u = 0; u < dual.N; ++u) {
        for (int v : dual.nbrs(u)) {
            EXPECT_EQ(int(S::next(view, u, v)), dual.next(u, v))
                << "dual next(" << u << ", " << v << ") mismatch";
            EXPECT_EQ(int(S::prev(view, u, v)), dual.prev(u, v))
                << "dual prev(" << u << ", " << v << ") mismatch";
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Sizes, CSRNavigation, ::testing::Values(20, 60, 80));

// ---------------------------------------------------------------------------
// Face computation
// ---------------------------------------------------------------------------

class CSRFaces : public ::testing::TestWithParam<int> {};

TEST_P(CSRFaces, CubicFaceCount) {
    int N = GetParam();
    FullereneGraph FG = get_fullerene(N);

    auto csr = csr_from_graph(FG);
    auto view = S::as_view(csr);
    auto fd = S::compute_faces(view);

    EXPECT_EQ(int(fd.n_faces), N / 2 + 2);
}

TEST_P(CSRFaces, CubicFaceSizes) {
    int N = GetParam();
    FullereneGraph FG = get_fullerene(N);

    auto csr = csr_from_graph(FG);
    auto view = S::as_view(csr);
    auto fd = S::compute_faces(view);

    int n5 = 0, n6 = 0;
    for (int f = 0; f < fd.n_faces; ++f) {
        int sz = S::face_size(view, fd, f);
        if (sz == 5) ++n5;
        else if (sz == 6) ++n6;
        else FAIL() << "Unexpected face size " << sz << " for face " << f;
    }
    EXPECT_EQ(n5, 12);
    EXPECT_EQ(n6, N / 2 + 2 - 12);
}

TEST_P(CSRFaces, CubicFaceBisimulation) {
    int N = GetParam();
    FullereneGraph FG = get_fullerene(N);

    auto csr = csr_from_graph(FG);
    auto view = S::as_view(csr);
    auto fd = S::compute_faces(view);

    // For each CSR face, trace the same face in the old graph starting from
    // the same directed edge and verify identical vertex sequences.
    for (int f = 0; f < fd.n_faces; ++f) {
        auto verts = S::face_vertices(view, fd, f);

        // Get the starting arc's (u,v) from CSR
        int p0 = fd.face_start[f];
        int u = S::source(view, p0);
        int v = csr.values[p0];

        // Trace the same face in the old graph
        std::vector<int> old_verts;
        int cu = u, cv = v;
        do {
            old_verts.push_back(cu);
            int w = FG.next_on_face(cu, cv);
            cu = cv;
            cv = w;
        } while (cu != u || cv != v);

        ASSERT_EQ(verts.size(), old_verts.size())
            << "Face " << f << " size mismatch (arc " << u << "->" << v << ")";
        for (int i = 0; i < int(verts.size()); ++i)
            EXPECT_EQ(int(verts[i]), old_verts[i])
                << "Face " << f << " vertex " << i << " mismatch";
    }
}

TEST_P(CSRFaces, DualFaceCount) {
    int N = GetParam();
    FullereneGraph FG = get_fullerene(N);

    PlanarGraph dual = FG.dual_graph();
    auto csr = csr_from_graph(dual);
    auto view = S::as_view(csr);
    auto fd = S::compute_faces(view);

    // Dual of cubic graph has N faces (one per cubic vertex, all triangles)
    EXPECT_EQ(int(fd.n_faces), N);
}

TEST_P(CSRFaces, DualAllFacesAreTriangles) {
    int N = GetParam();
    FullereneGraph FG = get_fullerene(N);

    PlanarGraph dual = FG.dual_graph();
    auto csr = csr_from_graph(dual);
    auto view = S::as_view(csr);
    auto fd = S::compute_faces(view);

    for (int f = 0; f < fd.n_faces; ++f)
        EXPECT_EQ(S::face_size(view, fd, f), 3) << "Face " << f << " is not a triangle";
}

INSTANTIATE_TEST_SUITE_P(Sizes, CSRFaces, ::testing::Values(20, 60, 80));

// ---------------------------------------------------------------------------
// Twin array correctness
// ---------------------------------------------------------------------------

TEST(CSRTwin, TwinIsInvolution) {
    FullereneGraph FG = FullereneGraph::C20();
    auto csr = csr_from_graph(FG);

    for (int p = 0; p < csr.n_arcs; ++p) {
        EXPECT_NE(csr.twin[p], -1) << "Arc " << p << " has no twin";
        EXPECT_EQ(csr.twin[csr.twin[p]], p)
            << "twin is not an involution at arc " << p;
    }
}

TEST(CSRTwin, TwinReversesArc) {
    FullereneGraph FG = FullereneGraph::C20();
    auto csr = csr_from_graph(FG);
    auto view = S::as_view(csr);

    for (int u = 0; u < csr.N; ++u) {
        for (int p = csr.offsets[u]; p < csr.offsets[u+1]; ++p) {
            int v = csr.values[p];
            int tp = csr.twin[p];
            EXPECT_EQ(csr.values[tp], u)
                << "twin of arc " << u << "->" << v << " doesn't point back to " << u;
            EXPECT_EQ(int(S::source(view, tp)), v)
                << "source of twin arc should be " << v;
        }
    }
}

// ---------------------------------------------------------------------------
// Arc-based navigation
// ---------------------------------------------------------------------------

TEST(CSRArcNav, NextArcCyclesBackCubic) {
    FullereneGraph FG = FullereneGraph::C20();
    auto csr = csr_from_graph(FG);
    auto view = S::as_view(csr);

    for (int u = 0; u < csr.N; ++u) {
        int p = csr.offsets[u];
        int q = p;
        for (int i = 0; i < 3; ++i)
            q = S::next_arc(view, q, u);
        EXPECT_EQ(q, p) << "3 next_arc steps didn't cycle back at vertex " << u;
    }
}

TEST(CSRArcNav, FaceTraversalCyclesBack) {
    FullereneGraph FG = FullereneGraph::C20();
    auto csr = csr_from_graph(FG);
    auto view = S::as_view(csr);

    for (int p = 0; p < csr.n_arcs; ++p) {
        int q = p;
        int steps = 0;
        do {
            q = S::next_on_face_arc(view, q);
            ++steps;
        } while (q != p && steps <= 10);
        EXPECT_EQ(q, p) << "Face traversal from arc " << p
                         << " didn't cycle back within 10 steps";
        EXPECT_GE(steps, 5);
        EXPECT_LE(steps, 6);
    }
}

// ---------------------------------------------------------------------------
// CSRBuilder
// ---------------------------------------------------------------------------

TEST(CSRBuilderTest, BuildC20FromEdges) {
    FullereneGraph FG = FullereneGraph::C20();

    S::CSRBuilder<int32_t> builder(20);
    for (int u = 0; u < 20; ++u)
        for (int v : FG.nbrs(u))
            builder.add_arc(u, v);

    auto csr = builder.freeze();
    EXPECT_TRUE(S::validate(S::as_view(csr)));
    EXPECT_EQ(csr.N, 20);
    EXPECT_EQ(csr.n_arcs, 60);

    auto view = S::as_view(csr);
    for (int v = 0; v < 20; ++v) {
        auto nbs = S::neighbours(view, v);
        auto old_nbs = FG.nbrs(v);
        ASSERT_EQ(int(nbs.size()), int(old_nbs.size()));
        for (int j = 0; j < int(nbs.size()); ++j)
            EXPECT_EQ(int(nbs[j]), old_nbs[j]);
    }
}
