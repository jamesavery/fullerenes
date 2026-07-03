// star-unfold-test: TriangulationView::star_unfold (@anchor tri-star-unfold).
//
// The deep invariants (exact closure, area identity) are asserted inside the
// method; these tests exercise the outcome codes, re-check the promised
// post-conditions externally, and verify polygon simplicity for the
// maximally tie-degenerate icosahedral source (outside the hypotheses of
// the classical Aronov-O'Rourke nonoverlap theorem).  Full-population
// evidence: the C20-C100 all-source sweep in
// claude-projects/visualization/validation/validate_star_unfold.cc.

#include "fullerenes/triangulation.hh"
#include "fullerenes/spiral.hh"

#include <gtest/gtest.h>

using SU = Triangulation::star_unfolding;

static Triangulation from_name(const std::string& spiral) {
    spiral_nomenclature sn(spiral);
    return Triangulation(sn);
}

// The (a,b)-coordinate orientation test carries over from cartesian exactly:
// the Eisenstein basis map is linear and orientation-preserving.
static long long orient(const Eisenstein& o, const Eisenstein& p, const Eisenstein& q) {
    return (long long)(p.first - o.first) * (q.second - o.second) -
           (long long)(p.second - o.second) * (q.first - o.first);
}

static bool properly_crossing(const Eisenstein& a, const Eisenstein& b,
                              const Eisenstein& c, const Eisenstein& d) {
    const long long d1 = orient(a, b, c), d2 = orient(a, b, d);
    const long long d3 = orient(c, d, a), d4 = orient(c, d, b);
    return ((d1 > 0) != (d2 > 0)) && d1 != 0 && d2 != 0 &&
           ((d3 > 0) != (d4 > 0)) && d3 != 0 && d4 != 0;
}

static void expect_star_postconditions(const Triangulation& T, node_t source) {
    const SU st = T.star_unfold(source);
    ASSERT_EQ(st.code, SU::Code::Ok);

    // the outline is the standard CW labelled 2(k-1)-gon: even entries are
    // source copies, odd entries the k-1 other cones, each exactly once
    size_t k = 0;
    for (node_t u = 0; u < (node_t)T.N; ++u) k += T.degree(u) != 6;
    ASSERT_EQ(st.outline.size(), 2 * (k - 1));

    std::vector<node_t> leaves;
    for (size_t i = 0; i < st.outline.size(); ++i) {
        if (i % 2 == 0) EXPECT_EQ(st.outline[i].second, source);
        else leaves.push_back(st.outline[i].second);
    }
    std::sort(leaves.begin(), leaves.end());
    EXPECT_TRUE(std::adjacent_find(leaves.begin(), leaves.end()) == leaves.end());

    // the walk starts at the origin; cut lengths derive from the outline
    EXPECT_EQ(st.outline.front().first, Eisenstein(0, 0));
    for (size_t i = 0; i + 1 < st.outline.size(); i += 2)
        EXPECT_GT((st.outline[i + 1].first - st.outline[i].first).norm2(), 0);
}

TEST(StarUnfold, C20AllSources) {
    const Triangulation T =
        from_name("[1,2,3,4,5,6,7,8,9,10,11,12]-fullerene");
    for (node_t u = 0; u < (node_t)T.N; ++u)
        expect_star_postconditions(T, u);     // every C20 vertex is a cone
}

TEST(StarUnfold, C60C70ByName) {
    expect_star_postconditions(
        from_name("[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene"), 0);
    expect_star_postconditions(
        from_name("[1,7,9,11,13,15,27,29,31,33,35,37]-fullerene"), 0);
}

// The star unfolding needs no iDT and no Alexandrov solve: it must succeed
// on the isomer whose iDT is multi-edged AND whose polytope is
// intrinsically degenerate.
TEST(StarUnfold, DegenerateMultiEdgeC96) {
    const Triangulation T = from_name(
        "[GS:1,2,12,17,23,29,35,40,41,45,49,50]-fullerene");
    expect_star_postconditions(T, 0);
    EXPECT_TRUE(T.star_unfold(0).cuts_globally_shortest);
}

// Icosahedral C60: maximal length ties at the source; the boundary must
// still be a simple polygon (0 proper crossings).
TEST(StarUnfold, SimplicityUnderIhTies) {
    const Triangulation T =
        from_name("[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene");
    const SU st = T.star_unfold(0);
    ASSERT_EQ(st.code, SU::Code::Ok);
    const auto& b = st.outline;
    const size_t n = b.size();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j) {
            if (j == i + 1 || (i == 0 && j == n - 1)) continue;
            EXPECT_FALSE(properly_crossing(b[i].first, b[(i + 1) % n].first,
                                           b[j].first, b[(j + 1) % n].first))
                << "outline edges " << i << " and " << j << " cross";
        }
}

TEST(StarUnfold, SourceNotConeCode) {
    const Triangulation T =
        from_name("[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene");
    node_t flat = 0;
    while (T.degree(flat) != 6) ++flat;       // C60 has 20 flat vertices
    const SU st = T.star_unfold(flat);
    EXPECT_EQ(st.code, SU::Code::SourceNotCone);
    EXPECT_EQ(st.metadata, std::to_string(flat));
}
