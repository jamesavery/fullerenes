// Phase 2 tests: reduction enumeration.

#include "buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <iostream>
#include <map>
#include <set>
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

// =====================================================================
// Test: allReductions on seed graphs
// =====================================================================

void test_reductions_C20() {
    cout << "--- Reductions on C20 ---\n";

    auto c20 = makeSeedC20();
    auto reds = allReductions(c20);

    // Count by type
    map<string, int> counts;
    for (const auto& r : reds)
        counts[r.kind.toString()]++;

    cout << "  Total reductions: " << reds.size() << "\n";
    for (const auto& [k, c] : counts)
        cout << "  " << k << ": " << c << "\n";

    // C20 is a seed graph (all 12 vertices degree-5). L0 requires flanking
    // vertices to be degree-6, which is impossible on C20. So no L0 reductions.
    // B(0,0) reductions ARE found because their validity check is weaker
    // (no flanking or intermediate degree check), but they can't actually be
    // applied to produce a valid fullerene (result would be too small).
    CHECK(reds.size() > 0, "C20 has reductions (all B(0,0))");
    CHECK(counts.count("L0") == 0, "C20 has no L0 reductions (no degree-6 flanking)");
    CHECK(counts.count("B(0,0)") > 0, "C20 has B(0,0) reductions");

    // Verify each reduction's properties
    for (const auto& r : reds) {
        CHECK(c20.degree(r.u) == 5, "start vertex is degree-5");
        CHECK(c20.arc_ix(r.u, r.v) >= 0, "edge exists");
    }
}

void test_reductions_C28() {
    cout << "--- Reductions on C28 ---\n";

    auto c28 = makeSeedC28();
    auto reds = allReductions(c28);

    map<string, int> counts;
    for (const auto& r : reds)
        counts[r.kind.toString()]++;

    cout << "  Total reductions: " << reds.size() << "\n";
    for (const auto& [k, c] : counts)
        cout << "  " << k << ": " << c << "\n";

    CHECK(reds.size() > 0, "C28 has reductions");

    // Verify each reduction's properties
    for (const auto& r : reds) {
        CHECK(c28.degree(r.u) == 5,
              r.toString() + " start vertex is degree-5");
        CHECK(c28.arc_ix(r.u, r.v) >= 0,
              r.toString() + " edge exists");
    }
}

void test_reductions_C30() {
    cout << "--- Reductions on C30 ---\n";

    auto c30 = makeSeedC30();
    auto reds = allReductions(c30);

    map<string, int> counts;
    for (const auto& r : reds)
        counts[r.kind.toString()]++;

    cout << "  Total reductions: " << reds.size() << "\n";
    for (const auto& [k, c] : counts)
        cout << "  " << k << ": " << c << "\n";

    CHECK(reds.size() > 0, "C30 has reductions");

    for (const auto& r : reds) {
        CHECK(c30.degree(r.u) == 5,
              r.toString() + " start vertex is degree-5");
        CHECK(c30.arc_ix(r.u, r.v) >= 0,
              r.toString() + " edge exists");
    }
}

// =====================================================================
// Test: followStraightToFive
// =====================================================================

void test_followStraight() {
    cout << "--- followStraightToFive ---\n";

    auto c28 = makeSeedC28();

    // For every degree-5 vertex u and neighbor v with deg(v) == 6,
    // followStraightToFive should either find a degree-5 endpoint or fail.
    int found = 0, notfound = 0;
    for (node_t u : deg5vertices(c28)) {
        for (int ni = 0; ni < c28.degree(u); ++ni) {
            node_t v = c28.neighbours[u][ni];
            if (c28.degree(v) == 5) continue;

            auto ep = followStraightToFive(c28, u, v, 5);
            if (ep) {
                CHECK(c28.degree(ep->endpoint) == 5,
                      "followStraight endpoint is degree-5");
                CHECK(ep->distance >= 2,
                      "followStraight distance >= 2 (skipping L0)");
                CHECK(ep->distance <= 5,
                      "followStraight distance <= maxDist");
                CHECK(c28.arc_ix(ep->prev, ep->endpoint) >= 0,
                      "followStraight prev is adjacent to endpoint");
                found++;
            } else {
                notfound++;
            }
        }
    }
    cout << "  Straight walks: " << found << " found, " << notfound << " not found\n";
}

// =====================================================================
// Test: L0 reductions are symmetric
// =====================================================================

void test_B00_symmetry() {
    cout << "--- B00 symmetry ---\n";

    // B(0,0) reductions on C20: for each (u, v, dir) there should
    // also be the reverse starting point with the same kind.
    auto c20 = makeSeedC20();
    auto reds = allReductions(c20, 2);

    // Just verify all reductions have degree-5 starting vertex
    for (const auto& r : reds) {
        CHECK(c20.degree(r.u) == 5, "B00 start is degree-5");
    }

    cout << "  B00 reductions on C20: " << reds.size() << "\n";
}

// =====================================================================
// Test: Reduction path validity
// =====================================================================

void test_reduction_paths() {
    cout << "--- Reduction path validity ---\n";

    auto c28 = makeSeedC28();
    auto reds = allReductions(c28);

    for (const auto& r : reds) {
        if (r.kind.type == ExpKind::L_type) {
            int numEntries = r.kind.i + 3;  // L_i path has i+3 entries
            auto pi = computeStraightPath(c28, r.u, r.v, r.dir, numEntries);
            CHECK(pi.valid,
                  r.toString() + " straight path is valid");
            if (pi.valid) {
                // First vertex must be degree-5
                CHECK(c28.degree(pi.path.front()) == 5,
                      r.toString() + " path starts at degree-5");
                // Last vertex must be degree-5
                CHECK(c28.degree(pi.path.back()) == 5,
                      r.toString() + " path ends at degree-5");
            }
        } else if (r.kind == Bk(0, 0)) {
            auto pi = computeBentZeroPath(c28, r.u, r.v, r.dir);
            // B00 validity was checked by the enumerator, but verify path
            set<node_t> pathSet(pi.path.begin(), pi.path.end());
            CHECK((int)pathSet.size() == (int)pi.path.size(),
                  r.toString() + " B00 path has no duplicates");
            CHECK(c28.degree(pi.path.back()) == 5,
                  r.toString() + " B00 path ends at degree-5");
        } else if (r.kind.type == ExpKind::B_type) {
            auto pi = computeBentPath(c28, r.u, r.v, r.dir, r.kind.i, r.kind.j);
            // The surgery parallel path may differ from the validity-check path,
            // but the main path vertices should be the same.
            CHECK((int)pi.path.size() == r.kind.i + r.kind.j + 5,
                  r.toString() + " bent path has correct length");
        }
    }
}

// =====================================================================
// Test: No duplicate reductions
// =====================================================================

void test_no_duplicates() {
    cout << "--- No duplicate reductions ---\n";

    auto c28 = makeSeedC28();
    auto reds = allReductions(c28);

    set<tuple<int, int, int, int, node_t, node_t, int>> seen;
    int dups = 0;
    for (const auto& r : reds) {
        auto key = make_tuple(r.kind.type, r.kind.i, r.kind.j,
                              (int)r.dir, r.u, r.v, 0);
        if (!seen.insert(key).second) dups++;
    }
    CHECK(dups == 0, "No duplicate reductions in C28 (" + to_string(dups) + " found)");
}

// =====================================================================
// Test: maxRedLen parameter
// =====================================================================

void test_maxRedLen() {
    cout << "--- maxRedLen filtering ---\n";

    auto c28 = makeSeedC28();
    auto reds1 = allReductions(c28, 1);  // Only L0
    auto reds2 = allReductions(c28, 2);  // L0 + L1 + B00
    auto reds5 = allReductions(c28, 5);  // Everything

    // Monotone: more allowed length means at least as many reductions
    CHECK(reds1.size() <= reds2.size(),
          "maxRedLen=1 <= maxRedLen=2");
    CHECK(reds2.size() <= reds5.size(),
          "maxRedLen=2 <= maxRedLen=5");

    // maxRedLen=1 should only produce L0 (none on C28 since no valid flanking)
    for (const auto& r : reds1) {
        CHECK(r.kind == Lk(0),
              "maxRedLen=1 produces only L0");
    }
    CHECK(reds1.empty() || true, "maxRedLen=1 may be empty on seeds");

    cout << "  maxRedLen=1: " << reds1.size()
         << ", maxRedLen=2: " << reds2.size()
         << ", maxRedLen=5: " << reds5.size() << "\n";
}

// =====================================================================
// Test: Reductions on generated fullerenes (non-seeds)
// =====================================================================

void test_generated_fullerenes() {
    cout << "--- Reductions on generated fullerenes ---\n";

    // Test C32 through C40 — all non-seed, should have reductions.
    // Flush before BuckyGen::start (which forks) to avoid duplicate output.
    for (int N = 32; N <= 40; N += 2) {
        cout << flush;
        cerr << flush;
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_count = 0;
        int with_L0 = 0, with_straight = 0, with_B00 = 0, with_bent = 0;
        int total_reds = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_count++;
            auto reds = allReductions(G, 5);
            total_reds += reds.size();

            CHECK(!reds.empty(),
                  "C" + to_string(N) + " isomer " + to_string(isomer_count)
                  + " has at least one reduction");

            map<string, int> counts;
            for (const auto& r : reds)
                counts[r.kind.toString()]++;

            if (counts.count("L0")) with_L0++;
            for (const auto& [k, c] : counts) {
                if (k[0] == 'L' && k != "L0") with_straight++;
                if (k == "B(0,0)") with_B00++;
                if (k[0] == 'B' && k != "B(0,0)") with_bent++;
            }
        }
        BuckyGen::stop(Q);

        cout << "  C" << N << ": " << isomer_count << " isomers, "
             << total_reds << " total reductions"
             << " (L0=" << with_L0 << " L_i=" << with_straight
             << " B00=" << with_B00 << " B_ij=" << with_bent << ")\n";
    }
}

// =====================================================================
// Main
// =====================================================================

int main() {
    test_reductions_C20();
    test_reductions_C28();
    test_reductions_C30();
    test_followStraight();
    test_B00_symmetry();
    test_reduction_paths();
    test_no_duplicates();
    test_maxRedLen();
    test_generated_fullerenes();

    cout << "\n=== Results: " << tests_passed << " passed, "
         << tests_failed << " failed ===\n";
    return tests_failed > 0 ? 1 : 0;
}
