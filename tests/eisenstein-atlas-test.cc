// eisenstein-atlas-test -- the flat atlas's own gates (stage 5b of the
// span promotion):
//
//   Determinism  -- build_atlas output is a pure function of its input:
//                   two builds compare equal on every flat table (no
//                   hash-iteration order anywhere), and the anchors are
//                   in scan order (nondecreasing host cell) -- the
//                   canonical-selection property itself, not just
//                   run-to-run agreement.
//   Resolution   -- resolve_all completes on every fixture.
//   Location     -- locate_sample at ALL THREE corners of every face
//                   lands on a lattice point claimed by the host cell
//                   with that corner's own vertex id.  Corner slots 1-2
//                   force genuine cross-cell traces (slot 0 alone is a
//                   near-no-op on some fixtures), so a single corrupted
//                   half-edge transition fails here.
//
// Fixtures: C60-IPR#0 / C80-IPR#0 / C28#0 have depth-0 routing (every
// face edge directly anchored); C60 general #0 (C1 symmetry) and
// C100#30 add depth-1 and depth-2 chains, so the midpoint-hop loop in
// resolve_face executes.  Charted complexes come from realize_dual
// (the Alexandrov-realized iDT): the raw dual_idt can hit the
// documented folded-development refusal on some isomers, which is the
// pipeline's business, not the atlas's.

#include "gtest/gtest.h"

#include "fullerenes/eisenstein_atlas.hh"
#include "fullerenes/eisenstein_paint.hh"
#include "fullerenes/eisenstein_paint_geometry.hh"   // realize_dual, DualPolytope
#include "fullerenes/buckygen-wrapper.hh"

#include <cstring>

using namespace eisenstein_paint;

namespace {

struct Charted {
    SortedDual   S;
    DualPolytope P;
    SurfaceParametrization param;
};

Triangulation nth_dual(int N, bool IPR, int idx) {
    auto Q = BuckyGen::start(N, IPR, false);
    Graph G;
    Triangulation T;
    bool found = false;
    for (int i = 0; BuckyGen::next_fullerene(Q, G); i++)
        if (i == idx) { T = Triangulation(G); found = true; break; }
    BuckyGen::stop(Q);
    EXPECT_TRUE(found) << "isomer C" << N << " #" << idx << " not found";
    return T;
}

Charted chart(int N, int idx, bool ipr) {
    Charted c;
    Triangulation T = nth_dual(N, ipr, idx);
    c.S = sorted_dual(T);
    c.P = realize_dual(c.S);
    c.param = parametrize(c.P.D, c.S);
    return c;
}

// Padding-free tables (ints and arrays of ints) byte-compare; the
// Eisenstein-carrying tables (he_trans, anchors) compare element-wise
// (LatticeIsometry has operator==; AnchorEdge's is defaulted) so struct
// padding can never enter the verdict.
template <class Vec>
bool bytes_equal(const Vec& a, const Vec& b) {
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(),
                                     a.size() * sizeof(a[0])) == 0);
}

void expect_deterministic(const SurfaceParametrization& P, const char* tag) {
    CellAtlas A1 = build_atlas(P);
    CellAtlas A2 = build_atlas(P);
    EXPECT_TRUE(A1.he_trans == A2.he_trans) << tag;
    EXPECT_TRUE(A1.anchors == A2.anchors) << tag;
    EXPECT_TRUE(bytes_equal(A1.tface, A2.tface)) << tag;
    EXPECT_TRUE(bytes_equal(A1.arc_face, A2.arc_face)) << tag;
    EXPECT_TRUE(bytes_equal(A1.arc_edge, A2.arc_edge)) << tag;
    EXPECT_TRUE(bytes_equal(A1.tedge, A2.tedge)) << tag;
    EXPECT_TRUE(bytes_equal(A1.edge_faces, A2.edge_faces)) << tag;
    EXPECT_TRUE(bytes_equal(A1.anchor_of_edge, A2.anchor_of_edge)) << tag;
    EXPECT_TRUE(bytes_equal(A1.bfs_parent_edge, A2.bfs_parent_edge)) << tag;
    EXPECT_TRUE(bytes_equal(A1.bfs_via_face, A2.bfs_via_face)) << tag;
    EXPECT_TRUE(bytes_equal(A1.bfs_depth, A2.bfs_depth)) << tag;
    // The canonical-selection property itself: anchors were chosen by
    // scanning cells in id order, so their host cells are nondecreasing.
    for (size_t i = 1; i < A1.anchors.size(); i++)
        EXPECT_LE(A1.anchors[i - 1].cell, A1.anchors[i].cell) << tag << " anchor " << i;
}

void expect_resolves_and_locates(const SurfaceParametrization& P, const char* tag) {
    CellAtlas A = build_atlas(P);
    ASSERT_NO_THROW(resolve_all(A)) << tag;
    for (int fi = 0; fi < (int)A.tface.size(); fi++) {
        // Every corner sample (den = 1, weight on slot k) must land on a
        // lattice point claimed by its host cell with the corner's own
        // vertex id.  Slots 1-2 route through trace_segment's cell
        // crossings, exercising he_trans on every fixture.
        for (int k = 0; k < 3; k++) {
            CellPoint cp = locate_sample(A, fi, k == 0, k == 1, k == 2, 1);
            ASSERT_EQ(cp.pos.den, 1) << tag << " face " << fi << " corner " << k;
            const LatticePoint* q = A.V.claim(cp.cell, cp.pos.num);
            ASSERT_NE(q, nullptr) << tag << " face " << fi << " corner " << k;
            EXPECT_EQ(q->vid, A.tface[fi][k]) << tag << " face " << fi << " corner " << k;
        }
    }
}

}  // namespace

TEST(EisensteinAtlas, DeterministicBuild) {
    const Charted a = chart(60, 0, true);
    const Charted b = chart(80, 0, true);
    const Charted c = chart(28, 0, false);
    const Charted d = chart(60, 0, false);
    const Charted e = chart(100, 30, false);
    expect_deterministic(a.param, "C60-IPR#0");
    expect_deterministic(b.param, "C80-IPR#0");
    expect_deterministic(c.param, "C28#0");
    expect_deterministic(d.param, "C60#0");
    expect_deterministic(e.param, "C100#30");
}

TEST(EisensteinAtlas, ResolveAllAndLocateCorners) {
    const Charted a = chart(60, 0, true);
    const Charted b = chart(80, 0, true);
    const Charted c = chart(28, 0, false);
    const Charted d = chart(60, 0, false);    // C1; depth-1 routing chains
    const Charted e = chart(100, 30, false);  // depth-2 routing chains
    expect_resolves_and_locates(a.param, "C60-IPR#0");
    expect_resolves_and_locates(b.param, "C80-IPR#0");
    expect_resolves_and_locates(c.param, "C28#0");
    expect_resolves_and_locates(d.param, "C60#0");
    expect_resolves_and_locates(e.param, "C100#30");
}
