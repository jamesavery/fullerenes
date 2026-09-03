// Law tests for the rotation-system words on RSRAdjacencyView
// (dense_graph.hh: target, next/prev = sigma/sigma^-1, reverse_arc =
// alpha, next_on_face, arcid/arc_of, find_arc).
//
// The laws are pinned here rather than left to incidental use: sigma's
// orbit is exactly the rotation, prev inverts next, alpha is an
// involution that flips the arc and agrees with an independently derived
// twin, next_on_face walks the library's ONE face convention (the
// vertex-pair form's prev(v,u)) with an anchored orbit so the mirror
// walk cannot pass, face orbits partition the arcs and satisfy Euler's
// formula, and the flat arc id round-trips.  Typed over both vertex-id
// widths, on two oriented planar fixtures of distinct degree profiles.
//
// The fixture is deliberately test-local (raw arrays; the twin derived
// from the ROTATION LISTS, not through the view's own find, so the
// alpha-vs-scan comparison is differential; padding slots poisoned with
// a plausible id so a degree-overrun scan cannot hide behind -1).
// Consolidating it with edge-surgery-test's TestRSR and the tetrahedron
// fixture oriented-surface-test also carries is queued as refactor debt
// rather than done here.

#include <gtest/gtest.h>
#include <algorithm>
#include <set>
#include <vector>

#include "fullerenes/dense_graph.hh"

using Spanify::RSRAdjacencyView;

template <typename K>
struct ArcNav : public ::testing::Test {
    using View = RSRAdjacencyView<K>;
    using arcix = typename View::arcix_t;

    struct Fixture {
        const char* name;
        std::vector<std::vector<int>> rot;
        int faces;
    };
    // Tetrahedron: all degree 3.  5-wheel (centre 0, rim 1..5 CCW, rim
    // rotations (next-rim, centre, prev-rim)): degrees 3 and 5, and rows
    // narrower than dmax, so padding slots are genuinely exercised.
    static std::vector<Fixture> fixtures() {
        return {{"tetrahedron",
                 {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}}, 4},
                {"wheel5",
                 {{1, 2, 3, 4, 5},
                  {2, 0, 5}, {3, 0, 1}, {4, 0, 2}, {5, 0, 3}, {1, 0, 4}}, 6}};
    }

    std::vector<K>       nbrs;
    std::vector<uint8_t> deg;
    std::vector<uint8_t> twin;
    int N = 0, dmax = 0;

    // Build from oriented rotations.  The twin is derived from the
    // ROTATION LISTS (std::find over rot), independent of the view's own
    // scan, so comparing reverse_arc against find_arc is a differential
    // check and not a tautology.  Padding slots hold vertex 0 -- a
    // plausible id, so a scan that overruns the degree returns a WRONG
    // slot instead of quietly missing.
    View build(const std::vector<std::vector<int>>& rot) {
        N    = (int)rot.size();
        dmax = 0;
        for (const auto& r : rot) dmax = std::max(dmax, (int)r.size());
        nbrs.assign((size_t)N * dmax, K(0));
        deg.assign((size_t)N, 0);
        twin.assign((size_t)N * dmax, 0);
        for (int u = 0; u < N; u++) {
            deg[(size_t)u] = (uint8_t)rot[(size_t)u].size();
            for (size_t i = 0; i < rot[(size_t)u].size(); i++)
                nbrs[(size_t)u * dmax + i] = (K)rot[(size_t)u][i];
        }
        for (int u = 0; u < N; u++)
            for (size_t i = 0; i < rot[(size_t)u].size(); i++) {
                const auto& rv = rot[(size_t)rot[(size_t)u][i]];
                const auto  it = std::find(rv.begin(), rv.end(), u);
                EXPECT_NE(it, rv.end()) << "asymmetric rotation in fixture";
                twin[(size_t)u * dmax + i] = (uint8_t)(it - rv.begin());
            }
        return View(K(N), dmax, nbrs, deg, twin);
    }

    // Every live arc of the graph, read off the graph itself.
    std::vector<arcix> arcs(const View& G) const {
        std::vector<arcix> as;
        for (K u = 0; u < G.N; ++u)
            for (int i = 0; i < G.degree(u); ++i)
                as.push_back({u, (uint8_t)i});
        return as;
    }
};

using Widths = ::testing::Types<int32_t, uint16_t>;
TYPED_TEST_SUITE(ArcNav, Widths);

TYPED_TEST(ArcNav, SigmaOrbitIsTheRotation) {
    for (const auto& f : this->fixtures()) {
        SCOPED_TRACE(f.name);
        auto G = this->build(f.rot);
        for (int u = 0; u < this->N; u++) {
            const int d = G.degree(TypeParam(u));
            // next's orbit from slot 0 visits every slot once and closes.
            auto a = typename decltype(G)::arcix_t{TypeParam(u), 0};
            std::set<int> seen;
            for (int k = 0; k < d; k++) { seen.insert(a.second); a = G.next(a); }
            EXPECT_EQ((int)seen.size(), d);
            EXPECT_EQ(a.second, 0);       // sigma^degree = id
        }
        for (auto a : this->arcs(G)) {
            EXPECT_EQ(G.prev(G.next(a)), a);   // prev . next = id
            EXPECT_EQ(G.next(G.prev(a)), a);
        }
    }
}

TYPED_TEST(ArcNav, AlphaIsAnInvolutionThatFlipsTheArc) {
    for (const auto& f : this->fixtures()) {
        SCOPED_TRACE(f.name);
        auto G = this->build(f.rot);
        ASSERT_TRUE(G.twin_is_valid());
        for (auto a : this->arcs(G)) {
            const auto r = G.reverse_arc(a);
            EXPECT_EQ(G.target(a), r.first);       // head of a = tail of alpha(a)
            EXPECT_EQ(G.target(r), a.first);       // and vice versa
            EXPECT_EQ(G.reverse_arc(r), a);        // involution
        }
    }
}

TYPED_TEST(ArcNav, AlphaAgreesWithTheFindScan) {
    // The twin-table read against the view's own scan.  Differential: the
    // fixture derived its twin from the rotation lists, not through find,
    // so the two sides share no code path below the arrays.
    for (const auto& f : this->fixtures()) {
        SCOPED_TRACE(f.name);
        auto G = this->build(f.rot);
        for (auto a : this->arcs(G))
            EXPECT_EQ(G.reverse_arc(a), G.find_arc(G.target(a), a.first));
    }
}

TYPED_TEST(ArcNav, PhiWalksTheSameFaceAsTheVertexPairForm) {
    // The library has ONE face convention -- the vertex-pair form's
    // next_on_face(u,v) = prev(v,u) (graphview.hh, planar_csr.hh) -- and
    // the arc form must walk it, not its mirror.  Orbit counts cannot see
    // the difference (a face and its mirror have the same size), so the
    // agreement is asserted arc by arc.
    for (const auto& f : this->fixtures()) {
        SCOPED_TRACE(f.name);
        auto G = this->build(f.rot);
        for (auto a : this->arcs(G)) {
            const auto v = G.target(a);
            const auto expect = G.prev(G.find_arc(v, a.first));   // prev(v,u)
            EXPECT_EQ(G.next_on_face(a), expect);
        }
    }
}

TYPED_TEST(ArcNav, PhiOrbitIsAnchored) {
    // One face spelled out, so a mirror walk fails on CONTENT: on the
    // 5-wheel the face left of 0->1 is the triangle 0,1,2, i.e. successive
    // targets 1, 2, 0.
    auto G = this->build(this->fixtures()[1].rot);
    auto a = typename decltype(G)::arcix_t{TypeParam(0), 0};      // 0 -> 1
    const std::vector<int> heads = {1, 2, 0};
    for (int h : heads) {
        EXPECT_EQ((int)G.target(a), h);
        a = G.next_on_face(a);
    }
    EXPECT_EQ(a, (typename decltype(G)::arcix_t{TypeParam(0), 0}));  // closed
}

TYPED_TEST(ArcNav, PhiOrbitsAreTheFaces) {
    // next_on_face's orbits partition the arcs into faces, and a
    // consistently oriented sphere embedding satisfies Euler's formula --
    // so an inconsistent fixture rotation fails HERE, not silently.
    for (const auto& f : this->fixtures()) {
        SCOPED_TRACE(f.name);
        auto G = this->build(f.rot);
        const auto all = this->arcs(G);
        std::set<size_t> unvisited;
        for (auto a : all) unvisited.insert(G.arcid(a));
        int orbits = 0;
        while (!unvisited.empty()) {
            auto a = G.arc_of(*unvisited.begin());
            const auto start = a;
            do {
                // Fatal: a double visit means the orbits do not partition,
                // and continuing would spin the outer loop forever.
                ASSERT_EQ(unvisited.erase(G.arcid(a)), 1u);
                a = G.next_on_face(a);
            } while (!(a == start));
            ++orbits;
            ASSERT_LE(orbits, (int)all.size());   // in-binary termination guard
        }
        EXPECT_EQ(orbits, f.faces);
        // Euler: N - E + F = 2, with E = arcs/2.
        EXPECT_EQ(this->N - (int)all.size() / 2 + orbits, 2);
    }
}

TYPED_TEST(ArcNav, ArcIdRoundTrips) {
    for (const auto& f : this->fixtures()) {
        SCOPED_TRACE(f.name);
        auto G = this->build(f.rot);
        for (auto a : this->arcs(G)) {
            EXPECT_EQ(G.arc_of(G.arcid(a)), a);
            EXPECT_EQ(G.arcid(a), G.arcid(a.first, a.second));
        }
    }
}

TYPED_TEST(ArcNav, FindArcOnAbsentEdgeCarriesNoSlot) {
    // The 5-wheel has absent edges (the tetrahedron is complete): rim
    // vertices 1 and 3 are not adjacent.
    auto G = this->build(this->fixtures()[1].rot);
    EXPECT_EQ(G.find_arc(TypeParam(1), TypeParam(3)).second,
              RSRAdjacencyView<TypeParam>::no_slot);
}
