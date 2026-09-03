// delaunay-unfold-test -- executing gates for unfold_iDT (delaunay_unfold)
// and the cell_developments builder:
//
//   Placement -- the best-of-seeds DFS unfolding places EVERY live charted
//                cell, with the GOLDEN tear count (10 on every fixture; an
//                unglued placement shows ~3x as many tears).
//   Gluing    -- every placed cell has >= 2 corners at the first-sighting
//                anchor positions the DFS glues to -- the cross-cell gate
//                an identity-placement mutant fails (a per-cell isometry
//                check alone cannot see gluing).
//   Isometry  -- per placed cell, the local->global map recomputed from the
//                corner triple (align on the 0->1 edge, anchored at corner
//                0) is a ROTATION and reproduces corners 1 and 2 AND every
//                lattice entry -- pinning slot pairing, chirality
//                (mirrored placements refused), and the apply arithmetic.
//   Builder   -- cell_developments agrees with the charts on every fixture
//                (corners; the accepted chart frame is among the cell's
//                enumerated developments) -- the builder's only executing
//                gate in this repo.
//   C980      -- the development-count refutation witness: GC(7,0) of the
//                icosahedron has 20 live iDT cells of side norm 49, EACH
//                admitting exactly 4 lattice developments, enumerated in
//                the canonical order (7,0), (5,3), (3,5), (0,7) (the
//                retired fixed-slot builder threw EMBED here).
//   Refusal   -- unfold_iDT throws on a mismatched (D, V) pair.
//
// Fixtures run on the RAW dual_idt (multi-edge delta-complexes included:
// C128#0's raw iDT is non-simplicial), so the slot-keyed pairing is
// exercised where label-keyed lookups would be ambiguous.

#include "gtest/gtest.h"

#include "fullerenes/delaunay_unfold.hh"
#include "fullerenes/eisenstein_paint.hh"
#include "fullerenes/eisenstein_paint_tables.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

using namespace eisenstein_paint;

namespace {

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

struct Charted {
    SortedDual            S;
    DelaunayTriangulation D;
    SurfaceParametrization param;
};

Charted chart(Triangulation T) {
    Charted c;
    c.S     = sorted_dual(T);
    c.D     = dual_idt(c.S);
    c.param = parametrize(c.D, c.S);
    return c;
}

// The placement + gluing + isometry gates on one charted complex.
// expect_tears is a GOLDEN (10 on every current fixture): an unfolding
// that stops gluing to the first-sighting anchors roughly triples it.
void expect_unfolds(const Charted& c, const char* tag, int expect_tears = 10) {
    const ParamTablesView V = c.param.view();

    int live = 0;
    for (int f = 0; f < V.nf; ++f) live += V.cell_live(f);

    const LatticeUnfolding U = unfold_iDT(c.D, V, -1, 20);
    EXPECT_EQ(U.n_cells, live) << tag << ": unfolding did not place all cells";
    EXPECT_EQ((int)U.cells.size(), U.n_cells) << tag;
    EXPECT_EQ(U.n_tears, expect_tears)
        << tag << ": tear count drifted from the golden";
    EXPECT_EQ(U.n_cones, V.n_cones) << tag;

    for (const LatticeUnfolding::UnfoldedCell& F : U.cells) {
        const ChartView ch = V.chart(F.cell_id);
        ASSERT_TRUE(ch.live()) << tag << " cell " << F.cell_id;
        EXPECT_EQ(F.c0, ch.corners[0]) << tag << " cell " << F.cell_id;
        EXPECT_EQ(F.c1, ch.corners[1]) << tag << " cell " << F.cell_id;
        EXPECT_EQ(F.c2, ch.corners[2]) << tag << " cell " << F.cell_id;

        // GLUING: the DFS anchors every non-seed cell's entering edge at
        // the first-sighting positions, so at least two corners of every
        // placed cell (all three of the seed) must sit exactly there.
        // This is the cross-cell gate -- an unfolding placing each cell
        // in isolation fails it immediately.
        int anchored = 0;
        for (auto [cone, p] : { std::pair<int, Eisenstein>{F.c0, F.P0},
                                {F.c1, F.P1}, {F.c2, F.P2} }) {
            const auto it = U.cone_positions.find(cone);
            ASSERT_TRUE(it != U.cone_positions.end())
                << tag << " cell " << F.cell_id << ": corner cone untracked";
            anchored += (int)(it->second.front() == p);
        }
        EXPECT_GE(anchored, 2)
            << tag << " cell " << F.cell_id
            << ": placed cell not glued to the first-sighting anchors";

        // ISOMETRY: recompute the local->global map INDEPENDENTLY from
        // the corner triple -- align the local 0->1 edge onto the global
        // one, anchor corner 0.  A development is orientation-preserving,
        // so the alignment must be a ROTATION, and it must reproduce
        // corners 1 and 2 (pinning slot order) and every entry.
        const auto L = ch.frame_points();
        ASSERT_EQ((L[1] - L[0]).norm2(), (F.P1 - F.P0).norm2())
            << tag << " cell " << F.cell_id << ": corner edge norm mismatch";
        const LatticeIsometry T = align(L[1] - L[0], F.P1 - F.P0);
        ASSERT_FALSE(T.reflect)
            << tag << " cell " << F.cell_id << ": mirrored placement";
        auto g = [&](Eisenstein p) { return (p - L[0]) * T.u + F.P0; };
        ASSERT_EQ(g(L[1]), F.P1) << tag << " cell " << F.cell_id;
        ASSERT_EQ(g(L[2]), F.P2)
            << tag << " cell " << F.cell_id
            << ": corner 2 not the rotation image -- slot pairing broken";

        const auto local = ch.entries;
        ASSERT_EQ(local.size(), F.entries.size()) << tag << " cell " << F.cell_id;
        for (size_t i = 0; i < local.size(); ++i) {
            EXPECT_EQ(g(local[i].pos()), F.entries[i].first)
                << tag << " cell " << F.cell_id << " entry " << i;
            EXPECT_EQ(local[i].vid, F.entries[i].second)
                << tag << " cell " << F.cell_id << " entry " << i;
        }
    }
}

// The builder gate: the pre-selection developments agree with the
// post-selection charts -- same live set and corners, and every accepted
// chart frame appears among that cell's enumerated developments.
void expect_developments(const Charted& c, const char* tag) {
    const ParamTablesView V = c.param.view();
    const CellDevelopments cd = cell_developments(c.D, c.S);
    ASSERT_EQ(cd.nf, V.nf) << tag;
    EXPECT_EQ(cd.n_cones, V.n_cones) << tag;
    for (int f = 0; f < V.nf; ++f) {
        const CellCorners cc = cd.corners(f);
        ASSERT_EQ(cc.c0 < 0, !V.cell_live(f))
            << tag << " cell " << f << ": dead-slot disagreement";
        if (!V.cell_live(f)) continue;
        const ChartView ch = V.chart(f);
        EXPECT_EQ(cc.c0, ch.corners[0]) << tag << " cell " << f;
        EXPECT_EQ(cc.c1, ch.corners[1]) << tag << " cell " << f;
        EXPECT_EQ(cc.c2, ch.corners[2]) << tag << " cell " << f;
        bool found = false;
        for (int k = 0; k < cd.n_developments(f) && !found; ++k) {
            const CellFrame Fk = cd.development(f, k).frame;
            found = Fk.p1a == ch.frame.p1a && Fk.p1b == ch.frame.p1b &&
                    Fk.p2a == ch.frame.p2a && Fk.p2b == ch.frame.p2b;
        }
        EXPECT_TRUE(found)
            << tag << " cell " << f << ": accepted frame not among its "
            << cd.n_developments(f) << " developments";
    }
}

// The RETIRED float development construction, frozen as a census
// reference (2026-08-07, when enumerate_developments went exact via
// place_third_eis_total): per sector-0 base P1, the apex from atan2 of
// the base + the interior angle + L20*cos/sin + lattice rounding,
// accepted by the exact norm/orientation identities.  Never used by the
// library; its whole value is being the unchanging pre-exactness
// formula the census below compares against.
std::vector<std::array<Eisenstein, 3>>
frozen_float_developments(const DelaunayTriangulation& D, int f) {
    const int h0 = D.f_he[f], h1 = D.he_next[h0], h2 = D.he_next[h1];
    const double L01 = D.he_length[h0], L20 = D.he_length[h2];
    const long   N01 = (long)std::lround(L01 * L01);
    const long   N12 = (long)std::lround(D.he_length[h1] * D.he_length[h1]);
    const long   N20 = (long)std::lround(L20 * L20);
    const double alpha_0 = D.he_angle[h0];
    std::vector<std::array<Eisenstein, 3>> out;
    const Eisenstein P0(0, 0);
    for (Eisenstein P1 : sector0_reps_of_norm((int)N01)) {
        auto [P1x, P1y] = P1.coord();
        const double theta_02 = std::atan2(P1y, P1x) + alpha_0;
        const Eisenstein P2(std::pair<double, double>{L20 * std::cos(theta_02),
                                                      L20 * std::sin(theta_02)});
        if ((long)P2.norm2() == N20 && (long)(P2 - P1).norm2() == N12 &&
            wedge(P1, P2) > 0)
            out.push_back({P0, P1, P2});
    }
    return out;
}

// Census: the exact enumeration must reproduce the retired float
// construction development-for-development (count, frames, order) on
// every live cell — a float-vs-exact disagreement would mean the float
// formula was mis-placing an apex on real inputs.
void expect_exact_matches_frozen_float(const Charted& c, const char* tag) {
    const CellDevelopments cd = cell_developments(c.D, c.S);
    for (int f = 0; f < cd.nf; ++f) {
        if (cd.corners(f).c0 < 0) continue;
        const auto frozen = frozen_float_developments(c.D, f);
        ASSERT_EQ(cd.n_developments(f), (int)frozen.size())
            << tag << " cell " << f << ": exact/frozen-float count differs";
        for (int k = 0; k < (int)frozen.size(); ++k) {
            const CellFrame Fk = cd.development(f, k).frame;
            EXPECT_EQ(Eisenstein(Fk.p1a, Fk.p1b), frozen[(size_t)k][1])
                << tag << " cell " << f << " dev " << k << " (P1)";
            EXPECT_EQ(Eisenstein(Fk.p2a, Fk.p2b), frozen[(size_t)k][2])
                << tag << " cell " << f << " dev " << k << " (P2)";
        }
    }
}

}  // namespace

TEST(DelaunayUnfold, ExactDevelopmentsMatchFrozenFloat) {
    expect_exact_matches_frozen_float(chart(nth_dual(20, false, 0)), "C20#0");
    expect_exact_matches_frozen_float(chart(nth_dual(60, false, 0)), "C60#0");
    expect_exact_matches_frozen_float(chart(nth_dual(128, false, 0)), "C128#0");
    expect_exact_matches_frozen_float(chart(nth_dual(20, false, 0).GCtransform(7, 0)),
                                      "C980");
}

TEST(DelaunayUnfold, C20) {
    Charted c = chart(nth_dual(20, false, 0));
    expect_unfolds(c, "C20#0");
    expect_developments(c, "C20#0");
}
TEST(DelaunayUnfold, C60) {
    Charted c = chart(nth_dual(60, false, 0));
    expect_unfolds(c, "C60#0");
    expect_developments(c, "C60#0");
}
TEST(DelaunayUnfold, C80IPR) {
    Charted c = chart(nth_dual(80, true, 0));
    expect_unfolds(c, "C80-IPR#0");
    expect_developments(c, "C80-IPR#0");
}
TEST(DelaunayUnfold, C128) {
    Charted c = chart(nth_dual(128, false, 0));
    expect_unfolds(c, "C128#0");
    expect_developments(c, "C128#0");
}

TEST(DelaunayUnfold, RefusesMismatchedTables) {
    Charted a = chart(nth_dual(20, false, 0));
    Charted b = chart(nth_dual(60, false, 0));
    EXPECT_THROW(unfold_iDT(a.D, b.param.view(), -1, 1), std::runtime_error);
    EXPECT_THROW(unfold_iDT(b.D, a.param.view(), -1, 1), std::runtime_error);
}

TEST(DelaunayUnfold, C980RefutationWitness) {
    Triangulation T20  = nth_dual(20, false, 0);
    Triangulation T980 = T20.GCtransform(7, 0);

    Charted c = chart(std::move(T980));

    // The lib-level refutation gate, as a CENSUS: every one of the 20
    // live side-norm-49 cells admits exactly 4 lattice developments
    // (80 total), enumerated in the canonical sector-0 order.
    const CellDevelopments cd = cell_developments(c.D, c.S);
    int live = 0;
    for (int f = 0; f < cd.nf; ++f)
        if (cd.corners(f).c0 >= 0) {
            ++live;
            EXPECT_EQ(cd.n_developments(f), 4)
                << "C980 cell " << f << ": expected exactly 4 developments";
        }
    EXPECT_EQ(live, 20) << "C980: expected 20 live cells";
    EXPECT_EQ(cd.dev_first[cd.nf], 80) << "C980: development census drifted";

    // Enumeration-ORDER golden on the first live cell: the P1 sequence
    // is sector0_reps_of_norm(49)'s order.  Selection depends on this
    // order (the first gluing development wins), so it is pinned.
    const Eisenstein want[4] = { {7, 0}, {5, 3}, {3, 5}, {0, 7} };
    for (int f = 0; f < cd.nf; ++f) {
        if (cd.corners(f).c0 < 0) continue;
        for (int k = 0; k < 4; ++k) {
            const CellFrame Fk = cd.development(f, k).frame;
            EXPECT_EQ(Eisenstein(Fk.p1a, Fk.p1b), want[k])
                << "C980 cell " << f << " development " << k
                << ": enumeration order drifted";
        }
        break;
    }

    expect_unfolds(c, "C980=GC(7,0)");
    expect_developments(c, "C980=GC(7,0)");
}
