// Phase 1 tests: navigation primitives, seed construction, path computation.

#include "buckinverse.hh"
#include <iostream>
#include <cassert>

using namespace buckinverse;
using namespace std;

// =====================================================================
// Test helpers
// =====================================================================

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (cond) { tests_passed++; } \
    else { tests_failed++; cerr << "FAIL: " << msg << " (line " << __LINE__ << ")\n"; } \
} while(0)

static bool isValidTriangulation(const Triangulation& t) {
    // All faces should be triangles
    // Euler: V - E + F = 2, E = 3(V-2)/1 for triangulation... no:
    // For a triangulation of the sphere: F = 2V - 4, E = 3V - 6
    int V = t.N;
    int E = 0;
    for (int v = 0; v < V; ++v) E += t.degree(v);
    E /= 2;  // each edge counted twice
    int F = E - V + 2;  // Euler formula
    bool euler_ok = (F == 2 * V - 4) && (E == 3 * V - 6);
    return euler_ok;
}

static bool allDegreesValid(const Triangulation& t) {
    for (int v = 0; v < t.N; ++v) {
        int d = t.degree(v);
        if (d != 5 && d != 6) return false;
    }
    return true;
}

static int countDeg5(const Triangulation& t) {
    int count = 0;
    for (int v = 0; v < t.N; ++v)
        if (t.degree(v) == 5) count++;
    return count;
}

// =====================================================================
// Test: Seed construction
// =====================================================================

void test_seeds() {
    cout << "--- Seed construction ---\n";

    // C20
    {
        auto c20 = makeSeedC20();
        CHECK(c20.N == 12, "C20 has 12 vertices");
        CHECK(allDegreesValid(c20), "C20 all degrees 5-6");
        CHECK(countDeg5(c20) == 12, "C20 has 12 degree-5 vertices (all)");
        CHECK(isValidTriangulation(c20), "C20 is a valid triangulation");
        CHECK(c20.is_oriented, "C20 is oriented");
        CHECK(identifySeed(c20) == SeedType::C20, "C20 identified as C20");
    }

    // C28
    {
        auto c28 = makeSeedC28();
        CHECK(c28.N == 16, "C28 has 16 vertices");
        CHECK(allDegreesValid(c28), "C28 all degrees 5-6");
        CHECK(countDeg5(c28) == 12, "C28 has 12 degree-5 vertices");
        CHECK(isValidTriangulation(c28), "C28 is a valid triangulation");
        CHECK(c28.is_oriented, "C28 is oriented");
        CHECK(identifySeed(c28) == SeedType::C28, "C28 identified as C28");
    }

    // C30
    {
        auto c30 = makeSeedC30();
        CHECK(c30.N == 17, "C30 has 17 vertices");
        CHECK(allDegreesValid(c30), "C30 all degrees 5-6");
        CHECK(countDeg5(c30) == 12, "C30 has 12 degree-5 vertices");
        CHECK(isValidTriangulation(c30), "C30 is a valid triangulation");
        CHECK(c30.is_oriented, "C30 is oriented");
        CHECK(identifySeed(c30) == SeedType::C30, "C30 identified as C30");
    }
}

// =====================================================================
// Test: Navigation primitives
// =====================================================================

void test_navigation() {
    cout << "--- Navigation primitives ---\n";

    auto c20 = makeSeedC20();

    // In C20, all vertices are degree-5.
    // Verify advanceCW wraps correctly.
    for (int v = 0; v < c20.N; ++v) {
        CHECK(c20.degree(v) == 5, "C20 vertex " + to_string(v) + " degree 5");

        // advanceCW by 0 should give the same neighbor
        node_t nb0 = c20.neighbours[v][0];
        CHECK(advanceCW(c20, v, nb0, 0) == nb0,
              "advanceCW(v," + to_string(nb0) + ",0) == " + to_string(nb0));

        // advanceCW by degree should wrap back
        CHECK(advanceCW(c20, v, nb0, 5) == nb0,
              "advanceCW wraps at degree");

        // next == advanceCW by 1
        CHECK(c20.next(v, nb0) == advanceCW(c20, v, nb0, 1),
              "next == advanceCW(1)");

        // prev == advanceCW by degree-1
        CHECK(c20.prev(v, nb0) == advanceCW(c20, v, nb0, 4),
              "prev == advanceCW(4)");
    }

    // Test straightAhead on C20 (all degree-5).
    // For degree-5, DRight: advance 3, DLeft: advance 2.
    {
        node_t u = 0;
        node_t v = c20.neighbours[u][0];
        node_t sa_r = straightAhead(c20, Dir::DRight, v, u);
        node_t sa_l = straightAhead(c20, Dir::DLeft, v, u);
        CHECK(sa_r == advanceCW(c20, v, u, 3),
              "straightAhead DRight at deg-5 == advance 3");
        CHECK(sa_l == advanceCW(c20, v, u, 2),
              "straightAhead DLeft at deg-5 == advance 2");
    }

    // Test sideNbr
    {
        node_t u = 0;
        node_t v = c20.neighbours[u][0];
        CHECK(sideNbr(c20, Dir::DRight, u, v) == c20.prev(u, v),
              "sideNbr DRight == prev");
        CHECK(sideNbr(c20, Dir::DLeft, u, v) == c20.next(u, v),
              "sideNbr DLeft == next");
    }

    // Test on a graph with degree-6 vertices (C28)
    auto c28 = makeSeedC28();
    {
        // Find a degree-6 vertex
        node_t v6 = -1;
        for (int v = 0; v < c28.N; ++v)
            if (c28.degree(v) == 6) { v6 = v; break; }
        CHECK(v6 >= 0, "C28 has a degree-6 vertex");

        node_t nb = c28.neighbours[v6][0];
        // At degree-6, straightAhead should always advance 3
        CHECK(straightAhead(c28, Dir::DRight, v6, nb) == advanceCW(c28, v6, nb, 3),
              "straightAhead DRight at deg-6 == advance 3");
        CHECK(straightAhead(c28, Dir::DLeft, v6, nb) == advanceCW(c28, v6, nb, 3),
              "straightAhead DLeft at deg-6 == advance 3");

        // turnAhead: DRight advance 2, DLeft advance 4
        CHECK(turnAhead(c28, Dir::DRight, v6, nb) == advanceCW(c28, v6, nb, 2),
              "turnAhead DRight at deg-6 == advance 2");
        CHECK(turnAhead(c28, Dir::DLeft, v6, nb) == advanceCW(c28, v6, nb, 4),
              "turnAhead DLeft at deg-6 == advance 4");
    }
}

// =====================================================================
// Test: Path computation
// =====================================================================

void test_paths() {
    cout << "--- Path computation ---\n";

    // Use C28 seed (16 dual vertices, has both deg-5 and deg-6) for path tests.
    auto c28 = makeSeedC28();

    // Find all degree-5 vertices
    auto d5 = deg5vertices(c28);
    CHECK((int)d5.size() == 12, "C28 deg5vertices returns 12");

    // Test straight path computation from each degree-5 vertex
    int valid_L0_count = 0;
    for (node_t u : d5) {
        for (int ni = 0; ni < c28.degree(u); ++ni) {
            node_t v = c28.neighbours[u][ni];
            for (Dir d : {Dir::DRight, Dir::DLeft}) {
                auto pi = computeStraightPath(c28, u, v, d, 3);  // L0: 3 entries
                if (pi.valid) {
                    CHECK((int)pi.path.size() == 3, "L0 path has 3 entries");
                    CHECK((int)pi.parallel.size() == 2, "L0 parallel has 2 entries");
                    CHECK(c28.arc_ix(pi.path[0], pi.path[1]) >= 0,
                          "path[0]-path[1] are adjacent");
                    CHECK(c28.arc_ix(pi.path[1], pi.path[2]) >= 0,
                          "path[1]-path[2] are adjacent");
                    valid_L0_count++;
                }
            }
        }
    }
    cout << "  Found " << valid_L0_count << " valid L0 paths in C28\n";

    // Test longer straight paths (L1: 4 entries)
    int valid_L1_count = 0;
    for (node_t u : d5) {
        for (int ni = 0; ni < c28.degree(u); ++ni) {
            node_t v = c28.neighbours[u][ni];
            for (Dir d : {Dir::DRight, Dir::DLeft}) {
                auto pi = computeStraightPath(c28, u, v, d, 4);  // L1: 4 entries
                if (pi.valid) {
                    CHECK((int)pi.path.size() == 4, "L1 path has 4 entries");
                    CHECK((int)pi.parallel.size() == 3, "L1 parallel has 3 entries");
                    valid_L1_count++;
                }
            }
        }
    }
    cout << "  Found " << valid_L1_count << " valid L1 paths in C28\n";

    // Test bent zero path
    int valid_B00_count = 0;
    for (node_t u : d5) {
        for (int ni = 0; ni < c28.degree(u); ++ni) {
            node_t v = c28.neighbours[u][ni];
            for (Dir d : {Dir::DRight, Dir::DLeft}) {
                auto pi = computeBentZeroPath(c28, u, v, d);
                if (pi.valid) {
                    CHECK((int)pi.path.size() == 5, "B00 path has 5 entries");
                    CHECK((int)pi.parallel.size() == 3, "B00 parallel has 3 entries");
                    valid_B00_count++;
                }
            }
        }
    }
    cout << "  Found " << valid_B00_count << " valid B(0,0) paths in C28\n";
}

// =====================================================================
// Test: Degree-5 enumeration helper
// =====================================================================

void test_deg5() {
    cout << "--- Degree-5 enumeration ---\n";

    auto c20 = makeSeedC20();
    auto d5 = deg5vertices(c20);
    CHECK((int)d5.size() == 12, "C20 has 12 degree-5 vertices");

    auto c28 = makeSeedC28();
    d5 = deg5vertices(c28);
    CHECK((int)d5.size() == 12, "C28 has 12 degree-5 vertices");

    auto c30 = makeSeedC30();
    d5 = deg5vertices(c30);
    CHECK((int)d5.size() == 12, "C30 has 12 degree-5 vertices");
}

// =====================================================================
// Test: Symmetry of navigation
// =====================================================================

void test_navigation_consistency() {
    cout << "--- Navigation consistency ---\n";

    auto c28 = makeSeedC28();

    // For every edge (u,v), verify that straightAhead and turnAhead
    // produce valid vertices (in range, actually neighbors of u).
    for (int u = 0; u < c28.N; ++u) {
        for (node_t v : c28.neighbours[u]) {
            for (Dir d : {Dir::DRight, Dir::DLeft}) {
                node_t sa = straightAhead(c28, d, u, v);
                CHECK(sa >= 0 && sa < c28.N,
                      "straightAhead in range");
                CHECK(c28.arc_ix(u, sa) >= 0,
                      "straightAhead result is a neighbor of u");

                node_t ta = turnAhead(c28, d, u, v);
                CHECK(ta >= 0 && ta < c28.N,
                      "turnAhead in range");
                CHECK(c28.arc_ix(u, ta) >= 0,
                      "turnAhead result is a neighbor of u");

                node_t sn = sideNbr(c28, d, u, v);
                CHECK(sn >= 0 && sn < c28.N,
                      "sideNbr in range");
                CHECK(c28.arc_ix(u, sn) >= 0,
                      "sideNbr result is a neighbor of u");
            }
        }
    }
}

// =====================================================================
// Main
// =====================================================================

int main() {
    test_seeds();
    test_navigation();
    test_paths();
    test_deg5();
    test_navigation_consistency();

    cout << "\n=== Results: " << tests_passed << " passed, "
         << tests_failed << " failed ===\n";
    return tests_failed > 0 ? 1 : 0;
}
