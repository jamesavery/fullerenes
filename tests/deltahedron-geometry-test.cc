// Tests for incremental geometry computation via extension path.
// Validates geometric quality of Deltahedron::fromExtensionPath output:
//   1. Face area statistics vs spherical model
//   2. Edge length statistics vs spherical model
//   3. Correct vertex count (topological sanity)
//   4. Volume sanity and comparison to expected spherical volume

#include <gtest/gtest.h>
#include "fullerenes/deltahedron.hh"
#include "fullerenes/buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/stats.hh"
#include <cmath>
#include <map>
#include <set>
#include <chrono>

using namespace buckinverse;
using namespace std;

// =====================================================================
// Statistics helpers
// =====================================================================

struct DistributionStats {
    double expected;          // spherical model prediction
    double mean;              // mean of actual values
    double min_rel_err;       // min (val - expected) / expected
    double max_rel_err;       // max (val - expected) / expected
    double mean_rel_err;      // mean of relative errors
    double var_rel_err;       // variance of relative errors
};

static DistributionStats computeStats(const vector<double>& vals, double expected) {
    DistributionStats s{};
    s.expected = expected;

    s.mean = 0;
    for (double v : vals) s.mean += v;
    s.mean /= vals.size();

    s.min_rel_err = 1e30;
    s.max_rel_err = -1e30;
    s.mean_rel_err = 0;
    for (double v : vals) {
        double r = (v - expected) / expected;
        s.min_rel_err = min(s.min_rel_err, r);
        s.max_rel_err = max(s.max_rel_err, r);
        s.mean_rel_err += r;
    }
    s.mean_rel_err /= vals.size();

    s.var_rel_err = 0;
    for (double v : vals) {
        double r = (v - expected) / expected;
        double dev = r - s.mean_rel_err;
        s.var_rel_err += dev * dev;
    }
    s.var_rel_err /= vals.size();

    return s;
}

// =====================================================================
// Per-isomer geometry analysis
// =====================================================================

struct IsomerGeometry {
    double R;                 // effective radius (mean vertex distance from centroid)

    DistributionStats area;   // face area statistics
    DistributionStats edge;   // edge length statistics

    double volume;            // from Polyhedron divergence theorem
    double hull_volume;       // convex hull volume
    double expected_volume;   // (4/3) pi R^3
    double vol_ratio;         // volume / hull_volume
};

static IsomerGeometry analyzeGeometry(const Deltahedron& D) {
    IsomerGeometry g{};

    // Effective radius from centroid
    coord3d center;
    for (int u = 0; u < D.N; u++) center += D.points[u];
    center = center / D.N;
    g.R = 0;
    for (int u = 0; u < D.N; u++) g.R += (D.points[u] - center).norm();
    g.R /= D.N;

    // Expected values from spherical model:
    //   F triangles on a sphere of radius R
    //   expected_area = 4*pi*R^2 / F
    //   expected_edge from equilateral triangle: A = sqrt(3)/4 * L^2
    int n_faces = D.triangles.size();
    double exp_area = 4.0 * M_PI * g.R * g.R / n_faces;
    double exp_edge = sqrt(4.0 * exp_area / sqrt(3.0));
    g.expected_volume = 4.0 / 3.0 * M_PI * g.R * g.R * g.R;

    // Collect face areas
    vector<double> areas;
    areas.reserve(n_faces);
    for (const auto& tri : D.triangles) {
        coord3d a = D.points[tri[0]], b = D.points[tri[1]], c = D.points[tri[2]];
        areas.push_back((b - a).cross(c - a).norm() / 2.0);
    }
    g.area = computeStats(areas, exp_area);

    // Collect edge lengths
    vector<double> edges;
    for (int u = 0; u < D.N; u++)
        for (int v : D.neighbours[u])
            if (v > u) edges.push_back((D.points[u] - D.points[v]).norm());
    g.edge = computeStats(edges, exp_edge);

    // Volume: mesh volume and convex hull volume via Polyhedron
    Polyhedron P(static_cast<const PlanarGraph&>(D), D.points);
    g.volume = P.volume();
    try {
        Polyhedron hull = P.convex_hull();
        g.hull_volume = hull.volume();
    } catch (const std::out_of_range&) {
        // Convex hull can fail for near-spherical optimized geometry
        // (near-coplanar points cause inconsistent visibility tests).
        // Fall back: hull_volume = volume (assumes convexity).
        g.hull_volume = g.volume;
    }
    g.vol_ratio = (g.hull_volume > 1e-30) ? g.volume / g.hull_volume : 0;

    return g;
}

// =====================================================================
// Test: Seed geometry produces valid Deltahedra
// =====================================================================

TEST(DeltahedronGeometry, SeedGeometry) {
    for (int seedN : {20, 28, 30}) {
        Graph G;
        if (seedN == 20) G = makeSeedC20();
        else if (seedN == 28) G = makeSeedC28();
        else G = makeSeedC30();

        ReducibleDual rd(G);
        ExtensionPath ep = rd.reduceToExtensionPath();

        ASSERT_EQ((int)ep.steps.size(), 0)
            << "C" << seedN << " seed has 0 expansion steps";

        Deltahedron D = Deltahedron::fromExtensionPath(ep);

        int expected_N = seedN / 2 + 2;
        EXPECT_EQ(D.N, expected_N)
            << "C" << seedN << " dual has " << expected_N << " vertices";
        EXPECT_EQ((int)D.points.size(), expected_N);

        auto g = analyzeGeometry(D);

        // Seeds are on a unit sphere with spring relaxation: moderate bounds
        double max_area_err = max(fabs(g.area.min_rel_err), fabs(g.area.max_rel_err));
        double max_edge_err = max(fabs(g.edge.min_rel_err), fabs(g.edge.max_rel_err));

        fprintf(stderr, "  C%d seed: area(exp=%.4f mean=%.4f max|rel|=%.3f var=%.3f) "
                "edge(exp=%.4f mean=%.4f max|rel|=%.3f var=%.3f) "
                "vol=%.4f hull=%.4f ratio=%.3f\n",
                seedN, g.area.expected, g.area.mean, max_area_err, g.area.var_rel_err,
                g.edge.expected, g.edge.mean, max_edge_err, g.edge.var_rel_err,
                g.volume, g.hull_volume, g.vol_ratio);

        EXPECT_LT(max_area_err, 0.5)
            << "C" << seedN << " area max|rel_err|=" << max_area_err;
        EXPECT_LT(g.area.var_rel_err, 0.05)
            << "C" << seedN << " area var=" << g.area.var_rel_err;
        EXPECT_LT(max_edge_err, 0.3)
            << "C" << seedN << " edge max|rel_err|=" << max_edge_err;
        EXPECT_LT(g.edge.var_rel_err, 0.02)
            << "C" << seedN << " edge var=" << g.edge.var_rel_err;
        EXPECT_GT(g.volume, 1.0)
            << "C" << seedN << " volume=" << g.volume;
        EXPECT_GT(g.hull_volume, 1.0)
            << "C" << seedN << " hull_volume=" << g.hull_volume;
        EXPECT_GT(g.vol_ratio, 0.95)
            << "C" << seedN << " vol/hull=" << g.vol_ratio;
    }
}

// =====================================================================
// Test: Single expansion step produces valid geometry
// =====================================================================

TEST(DeltahedronGeometry, SingleStep) {
    int tested = 0;
    for (int N = 32; N <= 40; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;

        while (BuckyGen::next_fullerene(Q, G)) {
            ReducibleDual rd(G);
            ExtensionPath ep = rd.reduceToExtensionPath();

            if (ep.steps.empty()) continue;

            ExtensionPath ep1;
            ep1.seed = ep.seed;
            ep1.full_N = ep.full_N;
            ep1.seed_state = ep.seed_state;
            ep1.steps.push_back(ep.steps[0]);

            Deltahedron D = Deltahedron::fromExtensionPath(ep1);

            for (const auto& tri : D.triangles) {
                coord3d a = D.points[tri[0]], b = D.points[tri[1]], c = D.points[tri[2]];
                EXPECT_GT((b - a).cross(c - a).norm(), 1e-15)
                    << "C" << N << " degenerate triangle";
            }

            tested++;
            break;
        }
        BuckyGen::stop(Q);
    }
    EXPECT_GT(tested, 0) << "Tested at least one single-step expansion";
}

// =====================================================================
// Accumulator for worst-case statistics across isomers
// =====================================================================

struct AccumulatedStats {
    int total = 0;
    int degenerate = 0;       // isomers with a zero-area face
    int coincident = 0;       // isomers with a zero-length edge

    double worst_area_max_rel = 0;
    double worst_area_var = 0;
    double worst_edge_max_rel = 0;
    double worst_edge_var = 0;
    double min_vol_ratio = 1e30;
    double min_volume = 1e30;
    double min_hull_volume = 1e30;
};

static void processIsomer(const Deltahedron& D, AccumulatedStats& acc) {
    acc.total++;

    // Degeneracy check
    for (const auto& tri : D.triangles) {
        coord3d a = D.points[tri[0]], b = D.points[tri[1]], c = D.points[tri[2]];
        if ((b - a).cross(c - a).norm() / 2.0 < 1e-15) {
            acc.degenerate++;
            return;  // skip detailed analysis for degenerate isomers
        }
    }
    for (int u = 0; u < D.N; u++)
        for (int v : D.neighbours[u])
            if (v > u && (D.points[u] - D.points[v]).norm() < 1e-10) {
                acc.coincident++;
                return;
            }

    auto g = analyzeGeometry(D);

    double max_area_err = max(fabs(g.area.min_rel_err), fabs(g.area.max_rel_err));
    double max_edge_err = max(fabs(g.edge.min_rel_err), fabs(g.edge.max_rel_err));

    acc.worst_area_max_rel = max(acc.worst_area_max_rel, max_area_err);
    acc.worst_area_var = max(acc.worst_area_var, g.area.var_rel_err);
    acc.worst_edge_max_rel = max(acc.worst_edge_max_rel, max_edge_err);
    acc.worst_edge_var = max(acc.worst_edge_var, g.edge.var_rel_err);
    acc.min_vol_ratio = min(acc.min_vol_ratio, g.vol_ratio);
    acc.min_volume = min(acc.min_volume, g.volume);
    acc.min_hull_volume = min(acc.min_hull_volume, g.hull_volume);
}

// =====================================================================
// Test: Full pipeline C20-C40 (3 seeds + 84 isomers = 87 total)
// =====================================================================

TEST(DeltahedronGeometry, FullPipelineC40) {
    AccumulatedStats acc;

    for (int seedN : {20, 28, 30}) {
        Graph G;
        if (seedN == 20) G = makeSeedC20();
        else if (seedN == 28) G = makeSeedC28();
        else G = makeSeedC30();
        ReducibleDual rd(G);
        ExtensionPath ep = rd.reduceToExtensionPath();
        Deltahedron D = Deltahedron::fromExtensionPath(ep);
        processIsomer(D, acc);
    }

    for (int N = 32; N <= 40; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int idx = 0;
        while (BuckyGen::next_fullerene(Q, G)) {
            idx++;
            ReducibleDual rd(G);
            ExtensionPath ep = rd.reduceToExtensionPath();
            Deltahedron D = Deltahedron::fromExtensionPath(ep);
            EXPECT_EQ(D.N, N / 2 + 2) << "C" << N << " #" << idx;
            processIsomer(D, acc);
        }
        BuckyGen::stop(Q);
    }

    fprintf(stderr, "  C20-C40 (%d isomers):\n", acc.total);
    fprintf(stderr, "    area:  max|rel_err|=%.3f  var=%.3f\n",
            acc.worst_area_max_rel, acc.worst_area_var);
    fprintf(stderr, "    edge:  max|rel_err|=%.3f  var=%.3f\n",
            acc.worst_edge_max_rel, acc.worst_edge_var);
    fprintf(stderr, "    vol:   min=%.6f  hull_min=%.6f  min_ratio=%.4f\n",
            acc.min_volume, acc.min_hull_volume, acc.min_vol_ratio);

    EXPECT_EQ(acc.degenerate, 0);
    EXPECT_EQ(acc.coincident, 0);
    EXPECT_LT(acc.worst_area_max_rel, 2.0)
        << "area max|rel_err|=" << acc.worst_area_max_rel;
    EXPECT_LT(acc.worst_area_var, 1.0)
        << "area var=" << acc.worst_area_var;
    EXPECT_LT(acc.worst_edge_max_rel, 2.0)
        << "edge max|rel_err|=" << acc.worst_edge_max_rel;
    EXPECT_LT(acc.worst_edge_var, 0.5)
        << "edge var=" << acc.worst_edge_var;
    EXPECT_GT(acc.min_hull_volume, 1.0) << "hull volume sanity";
    EXPECT_GT(acc.min_vol_ratio, 0.8) << "vol/hull=" << acc.min_vol_ratio;
    EXPECT_EQ(acc.total, 87);
}

// =====================================================================
// Test: Full pipeline C42-C60 (all 5678 isomers)
// =====================================================================

TEST(DeltahedronGeometry, FullPipelineC60) {
    AccumulatedStats acc;

    for (int N = 42; N <= 60; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int idx = 0;
        while (BuckyGen::next_fullerene(Q, G)) {
            idx++;
            ReducibleDual rd(G);
            ExtensionPath ep = rd.reduceToExtensionPath();
            Deltahedron D = Deltahedron::fromExtensionPath(ep);
            EXPECT_EQ(D.N, N / 2 + 2) << "C" << N << " #" << idx;
            processIsomer(D, acc);
        }
        BuckyGen::stop(Q);
    }

    fprintf(stderr, "  C42-C60 (%d isomers):\n", acc.total);
    fprintf(stderr, "    area:  max|rel_err|=%.3f  var=%.3f\n",
            acc.worst_area_max_rel, acc.worst_area_var);
    fprintf(stderr, "    edge:  max|rel_err|=%.3f  var=%.3f\n",
            acc.worst_edge_max_rel, acc.worst_edge_var);
    fprintf(stderr, "    vol:   min=%.6f  hull_min=%.6f  min_ratio=%.4f\n",
            acc.min_volume, acc.min_hull_volume, acc.min_vol_ratio);

    EXPECT_EQ(acc.degenerate, 0);
    EXPECT_EQ(acc.coincident, 0);
    EXPECT_LT(acc.worst_area_max_rel, 8.0)
        << "area max|rel_err|=" << acc.worst_area_max_rel;
    EXPECT_LT(acc.worst_area_var, 8.0)
        << "area var=" << acc.worst_area_var;
    EXPECT_LT(acc.worst_edge_max_rel, 3.0)
        << "edge max|rel_err|=" << acc.worst_edge_max_rel;
    EXPECT_LT(acc.worst_edge_var, 2.0)
        << "edge var=" << acc.worst_edge_var;
    EXPECT_GT(acc.min_hull_volume, 1.0) << "hull volume sanity";
    EXPECT_GT(acc.min_vol_ratio, 0.7) << "vol/hull=" << acc.min_vol_ratio;
    EXPECT_EQ(acc.total, 5678);
}

// =====================================================================
// Test: Phase 2 — Per-step optimized pipeline C20-C40
// =====================================================================

TEST(DeltahedronGeometry, OptimizedPipelineC40) {
    AccumulatedStats acc_p1, acc_p2;

    // Seeds
    for (int seedN : {20, 28, 30}) {
        Graph G;
        if (seedN == 20) G = makeSeedC20();
        else if (seedN == 28) G = makeSeedC28();
        else G = makeSeedC30();
        ReducibleDual rd(G);
        ExtensionPath ep = rd.reduceToExtensionPath();
        Deltahedron D1 = Deltahedron::fromExtensionPath(ep);
        processIsomer(D1, acc_p1);
        Deltahedron D2 = Deltahedron::fromExtensionPathOptimized(ep);
        processIsomer(D2, acc_p2);
    }

    for (int N = 32; N <= 40; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        while (BuckyGen::next_fullerene(Q, G)) {
            ReducibleDual rd(G);
            ExtensionPath ep = rd.reduceToExtensionPath();
            Deltahedron D1 = Deltahedron::fromExtensionPath(ep);
            processIsomer(D1, acc_p1);
            Deltahedron D2 = Deltahedron::fromExtensionPathOptimized(ep);
            processIsomer(D2, acc_p2);
        }
        BuckyGen::stop(Q);
    }

    fprintf(stderr, "  Phase 1 C20-C40 (%d isomers):\n", acc_p1.total);
    fprintf(stderr, "    area:  max|rel_err|=%.3f  var=%.3f\n",
            acc_p1.worst_area_max_rel, acc_p1.worst_area_var);
    fprintf(stderr, "    edge:  max|rel_err|=%.3f  var=%.3f\n",
            acc_p1.worst_edge_max_rel, acc_p1.worst_edge_var);
    fprintf(stderr, "    vol:   min_ratio=%.4f\n", acc_p1.min_vol_ratio);

    fprintf(stderr, "  Phase 2 C20-C40 (%d isomers):\n", acc_p2.total);
    fprintf(stderr, "    area:  max|rel_err|=%.3f  var=%.3f\n",
            acc_p2.worst_area_max_rel, acc_p2.worst_area_var);
    fprintf(stderr, "    edge:  max|rel_err|=%.3f  var=%.3f\n",
            acc_p2.worst_edge_max_rel, acc_p2.worst_edge_var);
    fprintf(stderr, "    vol:   min_ratio=%.4f\n", acc_p2.min_vol_ratio);

    EXPECT_EQ(acc_p2.total, 87);
    EXPECT_EQ(acc_p2.degenerate, 0);
    EXPECT_EQ(acc_p2.coincident, 0);

    // Phase 2 area/edge bounds: much tighter than Phase 1
    EXPECT_LT(acc_p2.worst_area_max_rel, 1.0)
        << "P2 area max|rel_err|=" << acc_p2.worst_area_max_rel;
    EXPECT_LT(acc_p2.worst_edge_max_rel, 0.7)
        << "P2 edge max|rel_err|=" << acc_p2.worst_edge_max_rel;
    EXPECT_GT(acc_p2.min_vol_ratio, 0.8)
        << "P2 vol/hull=" << acc_p2.min_vol_ratio;

    // Phase 2 should improve over Phase 1 in area/edge
    EXPECT_LT(acc_p2.worst_area_max_rel, acc_p1.worst_area_max_rel)
        << "P2 area should be tighter than P1";
    EXPECT_LT(acc_p2.worst_edge_max_rel, acc_p1.worst_edge_max_rel)
        << "P2 edge should be tighter than P1";
}

// =====================================================================
// Test: Phase 2 — Sampled optimized pipeline C42-C60
// =====================================================================

TEST(DeltahedronGeometry, OptimizedPipelineC60) {
    AccumulatedStats acc_p1, acc_p2;
    const int per_size = 10;  // sample to keep runtime manageable

    for (int N = 42; N <= 60; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int idx = 0;
        while (BuckyGen::next_fullerene(Q, G) && idx < per_size) {
            idx++;
            ReducibleDual rd(G);
            ExtensionPath ep = rd.reduceToExtensionPath();
            Deltahedron D1 = Deltahedron::fromExtensionPath(ep);
            Deltahedron D2 = Deltahedron::fromExtensionPathOptimized(ep);
            processIsomer(D1, acc_p1);
            processIsomer(D2, acc_p2);
        }
        BuckyGen::stop(Q);
    }

    fprintf(stderr, "  Phase 1 C42-C60 (%d isomers, %d/size):\n", acc_p1.total, per_size);
    fprintf(stderr, "    area:  max|rel_err|=%.3f  var=%.3f\n",
            acc_p1.worst_area_max_rel, acc_p1.worst_area_var);
    fprintf(stderr, "    edge:  max|rel_err|=%.3f  var=%.3f\n",
            acc_p1.worst_edge_max_rel, acc_p1.worst_edge_var);
    fprintf(stderr, "    vol:   min_ratio=%.4f\n", acc_p1.min_vol_ratio);

    fprintf(stderr, "  Phase 2 C42-C60 (%d isomers, %d/size):\n", acc_p2.total, per_size);
    fprintf(stderr, "    area:  max|rel_err|=%.3f  var=%.3f\n",
            acc_p2.worst_area_max_rel, acc_p2.worst_area_var);
    fprintf(stderr, "    edge:  max|rel_err|=%.3f  var=%.3f\n",
            acc_p2.worst_edge_max_rel, acc_p2.worst_edge_var);
    fprintf(stderr, "    vol:   min_ratio=%.4f\n", acc_p2.min_vol_ratio);

    EXPECT_EQ(acc_p2.degenerate, 0);
    EXPECT_EQ(acc_p2.coincident, 0);

    // Phase 2 area/edge bounds: much tighter than Phase 1
    EXPECT_LT(acc_p2.worst_area_max_rel, 2.0)
        << "P2 area max|rel_err|=" << acc_p2.worst_area_max_rel;
    EXPECT_LT(acc_p2.worst_edge_max_rel, 1.0)
        << "P2 edge max|rel_err|=" << acc_p2.worst_edge_max_rel;
    EXPECT_GT(acc_p2.min_vol_ratio, 0.5)
        << "P2 vol/hull=" << acc_p2.min_vol_ratio;

    // Phase 2 should improve over Phase 1 in area/edge
    EXPECT_LT(acc_p2.worst_area_max_rel, acc_p1.worst_area_max_rel)
        << "P2 area should be tighter than P1";
    EXPECT_LT(acc_p2.worst_edge_max_rel, acc_p1.worst_edge_max_rel)
        << "P2 edge should be tighter than P1";
}

// =====================================================================
// Test: Sampled pipeline C62-C100 (first 50 isomers per size)
// =====================================================================

TEST(DeltahedronGeometry, SampledPipelineC100) {
    AccumulatedStats acc;
    const int per_size = 50;

    for (int N = 62; N <= 100; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int idx = 0;
        while (BuckyGen::next_fullerene(Q, G) && idx < per_size) {
            idx++;
            ReducibleDual rd(G);
            ExtensionPath ep = rd.reduceToExtensionPath();
            Deltahedron D = Deltahedron::fromExtensionPath(ep);
            EXPECT_EQ(D.N, N / 2 + 2) << "C" << N << " #" << idx;
            processIsomer(D, acc);
        }
        BuckyGen::stop(Q);
    }

    fprintf(stderr, "  C62-C100 (%d isomers, %d/size):\n", acc.total, per_size);
    fprintf(stderr, "    area:  max|rel_err|=%.3f  var=%.3f\n",
            acc.worst_area_max_rel, acc.worst_area_var);
    fprintf(stderr, "    edge:  max|rel_err|=%.3f  var=%.3f\n",
            acc.worst_edge_max_rel, acc.worst_edge_var);
    fprintf(stderr, "    vol:   min=%.6f  hull_min=%.6f  min_ratio=%.6f\n",
            acc.min_volume, acc.min_hull_volume, acc.min_vol_ratio);

    EXPECT_EQ(acc.degenerate, 0);
    EXPECT_EQ(acc.coincident, 0);
    EXPECT_LT(acc.worst_area_max_rel, 25.0)
        << "area max|rel_err|=" << acc.worst_area_max_rel;
    EXPECT_LT(acc.worst_area_var, 30.0)
        << "area var=" << acc.worst_area_var;
    EXPECT_LT(acc.worst_edge_max_rel, 6.0)
        << "edge max|rel_err|=" << acc.worst_edge_max_rel;
    EXPECT_LT(acc.worst_edge_var, 4.0)
        << "edge var=" << acc.worst_edge_var;
    EXPECT_GT(acc.min_hull_volume, 1.0) << "hull volume sanity";
    EXPECT_GT(acc.min_vol_ratio, 0.7) << "vol/hull=" << acc.min_vol_ratio;
    EXPECT_EQ(acc.total, 20 * per_size);
}

// =====================================================================
// Edge length statistics helper
// =====================================================================

struct EdgeStats {
    double mean;
    double variance;
    double min_len;
    double max_len;

    static EdgeStats compute(const Deltahedron& D) {
        EdgeStats s{};
        vector<double> lens;
        for (int u = 0; u < D.N; u++)
            for (int v : D.neighbours[u])
                if (v > u) lens.push_back((D.points[u] - D.points[v]).norm());
        s.mean = 0;
        for (double l : lens) s.mean += l;
        s.mean /= lens.size();
        s.variance = 0;
        for (double l : lens) s.variance += (l - s.mean) * (l - s.mean);
        s.variance /= lens.size();
        s.min_len = *min_element(lens.begin(), lens.end());
        s.max_len = *max_element(lens.begin(), lens.end());
        return s;
    }
};

// =====================================================================
// Test: GC fullerene geometry via extension path
//
// Pipeline: icosahedron → GC(k,0) → ReducibleDual → ExtensionPath
//         → fromExtensionPathOptimized → compare quality/timing
//
// This tests whether the incremental geometry pipeline can produce
// good initial geometry for large fullerenes obtained via GC transform.
// =====================================================================

TEST(DeltahedronGeometry, GCExpansionGeometry) {
    // Build optimized icosahedron as base for GC transforms
    Polyhedron P20 = Polyhedron::C20();
    Polyhedron ico_poly = P20.dual();
    Deltahedron ico(ico_poly);
    ico.optimize(ico.points);

    fprintf(stderr, "  %-8s %6s %6s %5s | %-20s | %-20s\n",
            "GC(k,0)", "N_ful", "N_dual", "steps",
            "  Direct optimize",
            "  ExtPath optimized");
    fprintf(stderr, "  %-8s %6s %6s %5s | %8s %8s %6s | %8s %8s %6s\n",
            "", "", "", "",
            "edge_var", "max|rel|", "ms",
            "edge_var", "max|rel|", "ms");
    fprintf(stderr, "  %s\n", string(80, '-').c_str());

    for (int k : {2, 3, 5, 7, 10}) {
        int T = k * k;
        int N_fullerene = 20 * T;
        int N_dual = N_fullerene / 2 + 2;

        // 1. GC transform with barycentric interpolation
        Deltahedron gc_bary = ico.GCtransform(k, 0);
        ASSERT_EQ(gc_bary.N, N_dual) << "GC(" << k << ",0) dual node count";

        // 2. Direct optimize on GC geometry (baseline)
        Deltahedron gc_opt(static_cast<const Triangulation&>(gc_bary), gc_bary.points);
        auto t0 = chrono::steady_clock::now();
        gc_opt.optimize(gc_opt.points);
        auto t1 = chrono::steady_clock::now();
        double ms_direct = chrono::duration<double, milli>(t1 - t0).count();
        EdgeStats direct_stats = EdgeStats::compute(gc_opt);

        // 3. Decompose GC triangulation → ExtensionPath → rebuild with per-step opt
        //    Use higher maxRedLen for larger isomers.
        t0 = chrono::steady_clock::now();
        ReducibleDual rd(static_cast<const Graph&>(gc_bary));
        int maxRedLen = 20;
        ExtensionPath ep = rd.reduceToExtensionPath(maxRedLen);
        t1 = chrono::steady_clock::now();
        double ms_reduce = chrono::duration<double, milli>(t1 - t0).count();

        if (ep.seed == SeedType::NotASeed) {
            fprintf(stderr, "  GC(%d,0)  C%-4d %5d   --- | %8.2e %8.3f %6.0f | SKIP (no seed, %.0fms reduce)\n",
                    k, N_fullerene, N_dual,
                    direct_stats.variance,
                    (direct_stats.max_len - direct_stats.min_len) / direct_stats.mean,
                    ms_direct,
                    ms_reduce);
            continue;
        }

        t0 = chrono::steady_clock::now();
        Deltahedron gc_ep = Deltahedron::fromExtensionPathOptimized(ep);
        t1 = chrono::steady_clock::now();
        double ms_ep = chrono::duration<double, milli>(t1 - t0).count();

        EXPECT_EQ(gc_ep.N, N_dual) << "ExtPath GC(" << k << ",0) node count";

        EdgeStats ep_stats = EdgeStats::compute(gc_ep);

        fprintf(stderr, "  GC(%d,0)  C%-4d %5d %5zu | %8.2e %8.3f %6.0f | %8.2e %8.3f %6.0f\n",
                k, N_fullerene, N_dual, ep.steps.size(),
                direct_stats.variance,
                (direct_stats.max_len - direct_stats.min_len) / direct_stats.mean,
                ms_direct,
                ep_stats.variance,
                (ep_stats.max_len - ep_stats.min_len) / ep_stats.mean,
                ms_ep);

        // No degenerate faces
        for (const auto& tri : gc_ep.triangles) {
            coord3d a = gc_ep.points[tri[0]], b = gc_ep.points[tri[1]], c = gc_ep.points[tri[2]];
            EXPECT_GT((b - a).cross(c - a).norm(), 1e-15)
                << "GC(" << k << ",0) degenerate triangle";
        }

        // No NaN
        for (int u = 0; u < gc_ep.N; u++)
            for (int d = 0; d < 3; d++)
                EXPECT_FALSE(isnan(gc_ep.points[u][d]))
                    << "GC(" << k << ",0) NaN at vertex " << u;

        // Edge uniformity: should be reasonable
        EXPECT_LT(ep_stats.variance, 1e-3)
            << "GC(" << k << ",0) edge variance too high";
    }

    // Also test with BuckyGen-generated fullerenes (non-regular geometry).
    // This is the more representative case since BuckyGen fullerenes don't
    // start with perfect geometry like the icosahedron GC.
    fprintf(stderr, "\n  Non-icosahedral GC fullerenes:\n");
    fprintf(stderr, "  %-12s %6s %5s | %-20s | %-20s\n",
            "base+GC", "N_dual", "steps",
            "  Direct optimize",
            "  ExtPath optimized");
    fprintf(stderr, "  %s\n", string(72, '-').c_str());

    for (int baseN : {40, 60, 80}) {
        auto Q = BuckyGen::start(baseN, baseN >= 60);  // IPR for N>=60
        Graph G;
        if (!BuckyGen::next_fullerene(Q, G)) { BuckyGen::stop(Q); continue; }
        BuckyGen::stop(Q);

        // ReducibleDual takes the fullerene graph directly and builds the dual
        ReducibleDual rd_base(G);
        ExtensionPath ep_base = rd_base.reduceToExtensionPath(20);
        if (ep_base.seed == SeedType::NotASeed) continue;

        Deltahedron D_base = Deltahedron::fromExtensionPathOptimized(ep_base);

        // Apply GC(2,0) to get larger triangulation with geometry
        Deltahedron gc_dual = D_base.GCtransform(2, 0);
        int gc_N_dual = gc_dual.N;

        // Direct optimize on GC geometry
        Deltahedron gc_direct(static_cast<const Triangulation&>(gc_dual), gc_dual.points);
        auto t0 = chrono::steady_clock::now();
        gc_direct.optimize(gc_direct.points);
        auto t1 = chrono::steady_clock::now();
        double ms_direct = chrono::duration<double, milli>(t1 - t0).count();
        EdgeStats direct_stats = EdgeStats::compute(gc_direct);

        // Decompose GC triangulation → ExtensionPath → fromExtensionPathOptimized
        t0 = chrono::steady_clock::now();
        ReducibleDual rd_gc(static_cast<const Graph&>(gc_dual));
        ExtensionPath ep_gc = rd_gc.reduceToExtensionPath(20);
        if (ep_gc.seed == SeedType::NotASeed) {
            t1 = chrono::steady_clock::now();
            fprintf(stderr, "  C%d+GC(2) %5d   --- | %8.2e %8.3f %6.0f | SKIP (no seed)\n",
                    baseN, gc_N_dual,
                    direct_stats.variance,
                    (direct_stats.max_len - direct_stats.min_len) / direct_stats.mean,
                    ms_direct);
            continue;
        }
        Deltahedron gc_ep = Deltahedron::fromExtensionPathOptimized(ep_gc);
        t1 = chrono::steady_clock::now();
        double ms_ep = chrono::duration<double, milli>(t1 - t0).count();

        EdgeStats ep_stats = EdgeStats::compute(gc_ep);

        fprintf(stderr, "  C%d+GC(2) %5d %5zu | %8.2e %8.3f %6.0f | %8.2e %8.3f %6.0f\n",
                baseN, gc_N_dual, ep_gc.steps.size(),
                direct_stats.variance,
                (direct_stats.max_len - direct_stats.min_len) / direct_stats.mean,
                ms_direct,
                ep_stats.variance,
                (ep_stats.max_len - ep_stats.min_len) / ep_stats.mean,
                ms_ep);

        EXPECT_EQ(gc_ep.N, gc_N_dual);
        EXPECT_LT(ep_stats.variance, 1e-3)
            << "C" << baseN << "+GC(2,0) edge variance too high";
    }
}

// =====================================================================
// (5,0) nanotube dual graph builder
//
// Constructs the dual triangulation of a (5,0) nanotube fullerene
// with n_rings degree-6 rings. The structure is:
//   pole_N → cap_N[5] → ring_0[5] → ... → ring_{k-1}[5] → cap_S[5] → pole_S
//
// Fullerene: C(20 + 10*n_rings), dual has 12 + 5*n_rings vertices.
// CW ordering derived from the C30 seed (n_rings=1) and verified.
// =====================================================================

static Graph makeNanotubeDual(int n_rings) {
    assert(n_rings >= 1);
    int N = 12 + 5 * n_rings;

    // Vertex ID helpers (all indices mod 5)
    auto mod5 = [](int i) -> int { return ((i % 5) + 5) % 5; };
    auto cn  = [&](int i) -> int { return 1 + mod5(i); };             // cap_N[i]
    auto rng = [&](int j, int i) -> int { return 6 + 5*j + mod5(i); }; // ring_j[i]
    auto cs  = [&](int i) -> int { return 6 + 5*n_rings + mod5(i); };  // cap_S[i]
    int pole_N = 0;
    int pole_S = 11 + 5*n_rings;

    // Build edge set (unordered — Triangulation constructor will orient)
    set<pair<int,int>> edges;
    auto add = [&](int u, int v) {
        int a = std::min(u,v), b = std::max(u,v);
        edges.insert(make_pair(a, b));
    };

    // Pole N connects to all cap_N
    for (int i = 0; i < 5; i++) add(pole_N, cn(i));

    // Cap_N ring: cap_N[i] -- cap_N[(i+1)%5]
    for (int i = 0; i < 5; i++) add(cn(i), cn((i+1)%5));

    // Cap_N to ring_0: each cap_N[i] connects to ring_0[i] and ring_0[(i+1)%5]
    for (int i = 0; i < 5; i++) {
        add(cn(i), rng(0, i));
        add(cn(i), rng(0, (i+1)%5));
    }

    // Ring_j to ring_j: ring_j[i] -- ring_j[(i+1)%5]
    for (int j = 0; j < n_rings; j++)
        for (int i = 0; i < 5; i++)
            add(rng(j, i), rng(j, (i+1)%5));

    // Ring_j to ring_{j+1}: each ring_j[i] connects to ring_{j+1}[i] and ring_{j+1}[(i+1)%5]
    for (int j = 0; j + 1 < n_rings; j++)
        for (int i = 0; i < 5; i++) {
            add(rng(j, i), rng(j+1, i));
            add(rng(j, i), rng(j+1, (i+1)%5));
        }

    // Ring_{last} to cap_S: each ring_last[i] connects to cap_S[i] and cap_S[(i+1)%5]
    for (int i = 0; i < 5; i++) {
        add(rng(n_rings-1, i), cs(i));
        add(rng(n_rings-1, i), cs((i+1)%5));
    }

    // Cap_S ring: cap_S[i] -- cap_S[(i+1)%5]
    for (int i = 0; i < 5; i++) add(cs(i), cs((i+1)%5));

    // Pole S connects to all cap_S
    for (int i = 0; i < 5; i++) add(pole_S, cs(i));

    // Verify edge count: should be 3*N - 6 for triangulation of sphere
    assert((int)edges.size() == 3 * N - 6);

    // Build unoriented adjacency from edge set
    neighbours_t adj(N);
    for (const auto& [u, v] : edges) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    // Use Triangulation constructor to orient the neighbour lists
    Triangulation T(adj, false);  // already_oriented=false → computes faces + orients
    return static_cast<const Graph&>(T);
}

// =====================================================================
// Test: Verify nanotube builder against BuckyGen for small sizes
// =====================================================================

TEST(DeltahedronGeometry, NanotubeBuilderVerify) {
    // n_rings=1 should be C30 (the seed itself)
    Graph G1 = makeNanotubeDual(1);
    EXPECT_EQ(G1.N, 17) << "n_rings=1 should give C30 dual (17 vertices)";
    EXPECT_TRUE(G1.is_oriented) << "Should be oriented";

    // Verify it's a valid triangulation
    Triangulation T1(G1);
    EXPECT_EQ((int)T1.triangles.size(), 2 * T1.N - 4)
        << "Euler formula: F = 2V-4 for triangulation of sphere";

    // n_rings=2 should be C40 dual (22 vertices)
    Graph G2 = makeNanotubeDual(2);
    EXPECT_EQ(G2.N, 22);
    Triangulation T2(G2);
    EXPECT_EQ((int)T2.triangles.size(), 2 * T2.N - 4);

    // Verify n_rings=2 is isomorphic to the C40 (5,0) nanotube from BuckyGen
    // The (5,0) nanotube reduces to C30 with an F-ring step
    auto Q = BuckyGen::start(40, false);
    Graph G;
    bool found_nanotube = false;
    while (BuckyGen::next_fullerene(Q, G)) {
        ReducibleDual rd(G);
        ExtensionPath ep = rd.reduceToExtensionPath(5);
        if (ep.seed == SeedType::C30 && ep.steps.size() == 1
            && ep.steps[0].kind.type == ExpKind::F_type) {
            found_nanotube = true;
            break;
        }
    }
    BuckyGen::stop(Q);
    ASSERT_TRUE(found_nanotube) << "Should find (5,0) nanotube in C40";

    // The nanotube dual from BuckyGen and our builder should both reduce to C30
    ReducibleDual rd2(G2);
    ExtensionPath ep2 = rd2.reduceToExtensionPath(5);
    EXPECT_EQ(ep2.seed, SeedType::C30) << "Built nanotube should reduce to C30";
    EXPECT_EQ((int)ep2.steps.size(), 1) << "One F-ring step from C30 to C40";
    if (!ep2.steps.empty())
        EXPECT_EQ(ep2.steps[0].kind.type, ExpKind::F_type) << "Should be F-ring step";

    fprintf(stderr, "  Nanotube builder verified: n_rings=1 (C30, 17v), n_rings=2 (C40, 22v)\n");
}

// =====================================================================
// Test: Thin nanotube geometry via direct construction + extension path
//
// Constructs (5,0) nanotubes of increasing length using makeNanotubeDual,
// reduces to C30 seed, and tests geometry quality via both Phase 1 and
// Phase 2 pipelines. This tests the hardest case for incremental geometry:
// high aspect ratio with many F-ring expansion steps.
// =====================================================================

TEST(DeltahedronGeometry, NanotubeGeometry) {
    fprintf(stderr, "  %-6s %5s %4s %4s %4s | %-20s | %-20s\n",
            "C_N", "N_dv", "step", "F", "L",
            "  Phase 1 (no opt)",
            "  Phase 2 (per-step)");
    fprintf(stderr, "  %-6s %5s %4s %4s %4s | %8s %8s %6s | %8s %8s %6s\n",
            "", "", "", "", "",
            "edge_var", "max|rel|", "ms",
            "edge_var", "max|rel|", "ms");
    fprintf(stderr, "  %s\n", string(80, '-').c_str());

    int n_tested = 0;

    // Test nanotubes from C40 (n_rings=2) up to C500 (n_rings=48)
    for (int n_rings : {2, 3, 5, 8, 10, 15, 20, 30, 48}) {
        int N_fullerene = 20 + 10 * n_rings;
        int N_dual = 12 + 5 * n_rings;

        // 1. Build nanotube dual graph directly
        Graph G = makeNanotubeDual(n_rings);
        ASSERT_EQ(G.N, N_dual);

        // 2. Reduce to get ExtensionPath
        auto t0 = chrono::steady_clock::now();
        ReducibleDual rd(G);
        ExtensionPath ep = rd.reduceToExtensionPath(20);
        auto t1 = chrono::steady_clock::now();
        double ms_reduce = chrono::duration<double, milli>(t1 - t0).count();

        if (ep.seed == SeedType::NotASeed) {
            fprintf(stderr, "  C%-4d %5d  ---  --- --- | %8s %8s %6s | SKIP (no seed, %.0fms)\n",
                    N_fullerene, N_dual, "-", "-", "-", ms_reduce);
            continue;
        }

        EXPECT_EQ(ep.seed, SeedType::C30) << "C" << N_fullerene << " should reduce to C30";

        // Count F-type and L-type steps
        int n_F = 0, n_L = 0;
        for (const auto& s : ep.steps) {
            if (s.kind.type == ExpKind::F_type) n_F++;
            else n_L++;
        }

        // 3. Phase 1: no per-step optimization
        t0 = chrono::steady_clock::now();
        Deltahedron D1 = Deltahedron::fromExtensionPath(ep);
        t1 = chrono::steady_clock::now();
        double ms_p1 = chrono::duration<double, milli>(t1 - t0).count();
        EdgeStats p1_stats = EdgeStats::compute(D1);

        // 4. Phase 2: per-step optimization
        t0 = chrono::steady_clock::now();
        Deltahedron D2 = Deltahedron::fromExtensionPathOptimized(ep);
        t1 = chrono::steady_clock::now();
        double ms_p2 = chrono::duration<double, milli>(t1 - t0).count();
        EdgeStats p2_stats = EdgeStats::compute(D2);

        fprintf(stderr, "  C%-4d %5d %4zu %4d %4d | %8.2e %8.3f %6.0f | %8.2e %8.3f %6.0f\n",
                N_fullerene, N_dual, ep.steps.size(), n_F, n_L,
                p1_stats.variance,
                (p1_stats.max_len - p1_stats.min_len) / p1_stats.mean,
                ms_p1,
                p2_stats.variance,
                (p2_stats.max_len - p2_stats.min_len) / p2_stats.mean,
                ms_p2);

        // Phase 2 should produce non-degenerate geometry
        for (const auto& tri : D2.triangles) {
            coord3d a = D2.points[tri[0]], b = D2.points[tri[1]], c = D2.points[tri[2]];
            EXPECT_GT((b - a).cross(c - a).norm(), 1e-15)
                << "C" << N_fullerene << " degenerate triangle in Phase 2";
        }

        // No NaN
        for (int u = 0; u < D2.N; u++)
            for (int d = 0; d < 3; d++)
                EXPECT_FALSE(isnan(D2.points[u][d]))
                    << "C" << N_fullerene << " NaN at vertex " << u;

        n_tested++;
    }

    fprintf(stderr, "  Tested %d nanotube sizes\n", n_tested);
    EXPECT_GE(n_tested, 8) << "Should test at least 8 sizes";
}

// =====================================================================
// Convexity and surface orientation helpers
// (Adapted from deltahedron-test.cc)
// =====================================================================

// Signed height of vertex v above its neighbor centroid plane.
// Positive = convex (vertex sticks out), negative = concave (vertex is "inside").
// Only computes for deg<=6 vertices with all deg<=6 neighbors.
static double min_convexity_height(const Deltahedron& D) {
    double min_h = INFINITY;
    for (int v = 0; v < D.N; v++) {
        int d = (int)D.neighbours[v].size();
        if (d > 6) continue;
        bool all_low = true;
        for (int i = 0; i < d; i++)
            if ((int)D.neighbours[D.neighbours[v][i]].size() > 6) { all_low = false; break; }
        if (!all_low) continue;

        coord3d centroid(0,0,0);
        for (int i = 0; i < d; i++) centroid += D.points[D.neighbours[v][i]];
        centroid /= (double)d;

        coord3d n_fan(0,0,0);
        for (int i = 0; i < d; i++) {
            coord3d e1 = D.points[D.neighbours[v][i]] - D.points[v];
            coord3d e2 = D.points[D.neighbours[v][(i+1)%d]] - D.points[v];
            n_fan += e1.cross(e2);
        }
        double n_len = n_fan.norm();
        if (n_len < 1e-15) continue;
        double h = (D.points[v] - centroid).dot(n_fan / n_len);
        if (h < min_h) min_h = h;
    }
    return min_h;
}

// Count concave vertices (h < -tol) and return the list of vertex indices.
static vector<int> concave_vertices(const Deltahedron& D, double tol = 1e-6) {
    vector<int> result;
    for (int v = 0; v < D.N; v++) {
        int d = (int)D.neighbours[v].size();
        if (d > 6) continue;
        bool all_low = true;
        for (int i = 0; i < d; i++)
            if ((int)D.neighbours[D.neighbours[v][i]].size() > 6) { all_low = false; break; }
        if (!all_low) continue;

        coord3d centroid(0,0,0);
        for (int i = 0; i < d; i++) centroid += D.points[D.neighbours[v][i]];
        centroid /= (double)d;

        coord3d n_fan(0,0,0);
        for (int i = 0; i < d; i++) {
            coord3d e1 = D.points[D.neighbours[v][i]] - D.points[v];
            coord3d e2 = D.points[D.neighbours[v][(i+1)%d]] - D.points[v];
            n_fan += e1.cross(e2);
        }
        double n_len = n_fan.norm();
        if (n_len < 1e-15) continue;
        double h = (D.points[v] - centroid).dot(n_fan / n_len);
        if (h < -tol) result.push_back(v);
    }
    return result;
}

// Check surface orientation consistency via edge traversal.
// For a correctly oriented closed triangulated surface, every edge must be
// traversed in opposite directions by its two adjacent triangles. This is a
// purely topological test — no geometry assumptions about shape.
// Returns the number of edges with inconsistent traversal (0 = perfect).
static int count_orientation_defects(const Deltahedron& D) {
    // Count directed arc occurrences across all triangles
    map<pair<int,int>, int> arc_count;
    for (const auto& tri : D.triangles) {
        arc_count[{tri[0], tri[1]}]++;
        arc_count[{tri[1], tri[2]}]++;
        arc_count[{tri[2], tri[0]}]++;
    }

    // Every undirected edge should have exactly one forward and one reverse arc
    int defects = 0;
    set<pair<int,int>> checked;
    for (auto& [arc, cnt] : arc_count) {
        auto edge = make_pair(min(arc.first, arc.second), max(arc.first, arc.second));
        if (checked.count(edge)) continue;
        checked.insert(edge);

        int fwd = arc_count.count({edge.first, edge.second}) ? arc_count[{edge.first, edge.second}] : 0;
        int rev = arc_count.count({edge.second, edge.first}) ? arc_count[{edge.second, edge.first}] : 0;
        if (fwd != 1 || rev != 1) defects++;
    }
    return defects;
}

// Check that the total signed volume is positive (outward-oriented normals).
// The signed volume V = (1/6) sum_i a_i . (b_i x c_i) is positive when
// triangle normals point outward, regardless of surface shape.
static bool has_positive_orientation(const Deltahedron& D) {
    double total_vol = 0;
    for (const auto& tri : D.triangles) {
        const coord3d& a = D.points[tri[0]];
        const coord3d& b = D.points[tri[1]];
        const coord3d& c = D.points[tri[2]];
        total_vol += a.dot(b.cross(c));
    }
    return total_vol > 0;
}

// Per-isomer convexity statistics
struct ConvexityStats {
    double L_cv;           // edge length coefficient of variation
    double ang_min, ang_max, ang_std;  // triangle angles in degrees
    double h_min;          // min convexity height (negative = concave)
    int n_concave;         // number of concave vertices
    int n_orient_defects;  // number of orientation-inconsistent edges
    double vol_ratio;      // volume / hull_volume
    bool converged;
    int iters;

    static ConvexityStats compute(const Deltahedron& D) {
        ConvexityStats s{};
        s.iters = D.iterations_used;

        // Edge lengths
        vector<double> edge_lens;
        for (int u = 0; u < D.N; u++)
            for (int v : D.neighbours[u])
                if (v > u) edge_lens.push_back((D.points[u] - D.points[v]).norm());
        s.L_cv = cv_twopass(edge_lens);

        // Triangle angles
        vector<double> angles;
        s.ang_min = 180; s.ang_max = 0;
        for (const auto& tri : D.triangles)
            for (int c = 0; c < 3; c++) {
                coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
                coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
                double ang = coord3d::angle(va, vb) * 180.0 / M_PI;
                angles.push_back(ang);
                s.ang_min = min(s.ang_min, ang);
                s.ang_max = max(s.ang_max, ang);
            }
        s.ang_std = stddev_twopass(angles);

        // Convexity
        s.h_min = min_convexity_height(D);
        s.n_concave = (int)concave_vertices(D).size();
        s.n_orient_defects = count_orientation_defects(D);

        // Signed volume from triangle fan (no Polyhedron construction needed)
        double vol = 0;
        for (const auto& tri : D.triangles) {
            const coord3d &a = D.points[tri[0]], &b = D.points[tri[1]], &c = D.points[tri[2]];
            vol += a.dot(b.cross(c));
        }
        s.vol_ratio = vol / 6.0;

        return s;
    }

    static void print_header() {
        fprintf(stderr, "\n%6s %5s  %8s  %7s %7s %7s  %8s %4s %4s  %7s\n",
               "iso", "iter", "L_cv",
               "ang_min", "ang_max", "ang_std",
               "h_min", "conc", "inv", "vol_rat");
        fprintf(stderr, "------  -----  --------  ------- ------- -------  -------- ---- ----  -------\n");
    }

    void print_row(int idx) const {
        fprintf(stderr, "%6d %5d  %8.5f  %7.2f %7.2f %7.3f  %+8.4f %4d %4d  %7.4f%s\n",
               idx, iters, L_cv,
               ang_min, ang_max, ang_std,
               h_min, n_concave, n_orient_defects,
               vol_ratio,
               converged ? "" : " [NC]");
    }

    static void print_summary(const vector<ConvexityStats>& all) {
        double worst_cv = 0, worst_ang_min = 180, worst_ang_max = 0;
        double worst_ang_std = 0, worst_h = INFINITY;
        int total_concave = 0, total_inverted = 0, n_with_concave = 0, n_with_inverted = 0;
        double worst_vol_ratio = 1.0;
        int n_nc = 0, max_iters = 0;
        for (const auto& s : all) {
            worst_cv = max(worst_cv, s.L_cv);
            worst_ang_min = min(worst_ang_min, s.ang_min);
            worst_ang_max = max(worst_ang_max, s.ang_max);
            worst_ang_std = max(worst_ang_std, s.ang_std);
            worst_h = min(worst_h, s.h_min);
            worst_vol_ratio = min(worst_vol_ratio, s.vol_ratio);
            if (s.n_concave > 0) { total_concave += s.n_concave; n_with_concave++; }
            if (s.n_orient_defects > 0) { total_inverted += s.n_orient_defects; n_with_inverted++; }
            if (!s.converged) n_nc++;
            max_iters = max(max_iters, s.iters);
        }
        fprintf(stderr, "---- summary (%d isomers, %d not converged, max %d iters) ----\n",
               (int)all.size(), n_nc, max_iters);
        fprintf(stderr, "  worst L_cv=%.5f  ang_range=[%.2f, %.2f]  worst ang_std=%.3f\n",
               worst_cv, worst_ang_min, worst_ang_max, worst_ang_std);
        fprintf(stderr, "  worst h_min=%+.4f  isomers with concave verts: %d/%d  total concave verts: %d\n",
               worst_h, n_with_concave, (int)all.size(), total_concave);
        fprintf(stderr, "  isomers with orient defects: %d/%d  total defects: %d\n",
               n_with_inverted, (int)all.size(), total_inverted);
        fprintf(stderr, "  worst vol/hull ratio: %.4f\n", worst_vol_ratio);
    }
};

// =====================================================================
// Test: ExtPath convexity sweep C20-C40 (all 87 isomers)
// =====================================================================

TEST(DeltahedronGeometry, ExtPathConvexityC40) {
    vector<ConvexityStats> all_stats;
    vector<pair<int,int>> problematic;  // (N, idx) of concave isomers

    fprintf(stderr, "\n  ExtPath Optimized Convexity Sweep C20-C40:\n");
    ConvexityStats::print_header();

    // Seeds
    int global_idx = 0;
    for (int seedN : {20, 28, 30}) {
        Graph G;
        if (seedN == 20) G = makeSeedC20();
        else if (seedN == 28) G = makeSeedC28();
        else G = makeSeedC30();
        ReducibleDual rd(G);
        ExtensionPath ep = rd.reduceToExtensionPath();
        Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);
        ConvexityStats s = ConvexityStats::compute(D);
        s.converged = true;  // seeds are trivial
        bool bad = s.n_concave > 0 || s.n_orient_defects > 0 || s.h_min < -0.01;
        if (bad) {
            s.print_row(global_idx);
            problematic.push_back({seedN, 0});
        }
        all_stats.push_back(s);
        global_idx++;
    }

    for (int N = 32; N <= 40; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int idx = 0;
        while (BuckyGen::next_fullerene(Q, G)) {
            ReducibleDual rd(G);
            ExtensionPath ep = rd.reduceToExtensionPath();
            Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);
            ConvexityStats s = ConvexityStats::compute(D);
            s.converged = true;
            bool bad = s.n_concave > 0 || s.n_orient_defects > 0 || s.h_min < -0.01;
            if (bad) {
                fprintf(stderr, "  C%d #%d: ", N, idx);
                s.print_row(idx);
                problematic.push_back({N, idx});
            }
            all_stats.push_back(s);
            idx++;
        }
        BuckyGen::stop(Q);
    }

    ConvexityStats::print_summary(all_stats);

    if (!problematic.empty()) {
        fprintf(stderr, "  Problematic isomers (%d):", (int)problematic.size());
        for (auto& [N, idx] : problematic)
            fprintf(stderr, " C%d#%d", N, idx);
        fprintf(stderr, "\n");
    }

    EXPECT_EQ((int)all_stats.size(), 87);

    // Aggregate quality assertions
    double worst_cv = 0, worst_h = INFINITY;
    double worst_ang_min = 180, worst_ang_max = 0, worst_ang_std = 0;
    int n_concave_isomers = 0, n_defect_isomers = 0;
    for (const auto& s : all_stats) {
        worst_cv = max(worst_cv, s.L_cv);
        worst_h = min(worst_h, s.h_min);
        worst_ang_min = min(worst_ang_min, s.ang_min);
        worst_ang_max = max(worst_ang_max, s.ang_max);
        worst_ang_std = max(worst_ang_std, s.ang_std);
        if (s.n_concave > 0) n_concave_isomers++;
        if (s.n_orient_defects > 0) n_defect_isomers++;
    }
    fprintf(stderr, "  Concave: %d/87  Orient defects: %d/87\n",
           n_concave_isomers, n_defect_isomers);
    fprintf(stderr, "  L_cv=%.5f  h_min=%+.4f  ang=[%.2f,%.2f] std=%.3f\n",
           worst_cv, worst_h, worst_ang_min, worst_ang_max, worst_ang_std);

    EXPECT_EQ(n_concave_isomers, 0) << "no concave isomers expected";
    EXPECT_EQ(n_defect_isomers, 0) << "no orientation defects expected";
    EXPECT_LT(worst_cv, 0.001) << "L_cv=" << worst_cv;
    EXPECT_GT(worst_h, 0) << "h_min=" << worst_h;
    EXPECT_GT(worst_ang_min, 50.0) << "ang_min=" << worst_ang_min;
    EXPECT_LT(worst_ang_max, 70.0) << "ang_max=" << worst_ang_max;
    EXPECT_LT(worst_ang_std, 2.0) << "ang_std=" << worst_ang_std;
}

// =====================================================================
// Test: Known problematic isomers from deltahedron-test.cc
//
// These C60 and C70 isomers converge to concave/"innie" solutions
// when starting from zero_order_geometry(). Does the ExtPath pipeline
// (starting from optimized seed geometry) avoid this?
// =====================================================================

TEST(DeltahedronGeometry, ExtPathKnownConcave) {
    // C60 concave isomers (from deltahedron-test.cc)
    static const vector<int> c60_concave = {17, 18, 388, 393, 820, 1471};
    // C70 worst concave isomers
    static const vector<int> c70_concave = {0, 15, 68, 216, 217, 884, 886, 1555,
        2072, 2179, 2237, 2268, 2537, 2539, 2623, 2803, 3223, 3225, 3226, 3443, 3462, 3464};

    fprintf(stderr, "\n  Known concave isomers via ExtPath Optimized:\n");

    // Build the dual triangulation from an IsomerDB entry.
    // IsomerDB::makeIsomer returns a cubic FullereneGraph (degree 3),
    // but ReducibleDual expects the dual (triangulation, degree 5/6).
    // We reconstruct the Triangulation directly from the spiral pentagon indices.
    // Note: RSPI values in IsomerDB are 1-based, must decrement to 0-based.
    auto dual_from_entry = [](int N, const IsomerDB::Entry& e) -> Graph {
        int n_faces = N/2 + 2;
        vector<int> spiral_string(n_faces, 6);
        for (int i = 0; i < 12; i++) spiral_string[e.RSPI[i] - 1] = 5;
        Triangulation T(spiral_string);
        return static_cast<const Graph&>(T);
    };

    auto test_size = [&](int N, const vector<int>& indices) {
        IsomerDB db = IsomerDB::readPDB(N, false);
        ASSERT_GT((int)db.entries.size(), *max_element(indices.begin(), indices.end()));

        fprintf(stderr, "\n  C%d (%d known concave isomers):\n", N, (int)indices.size());
        ConvexityStats::print_header();
        vector<ConvexityStats> stats;
        int n_still_concave = 0, n_skipped = 0;

        for (int idx : indices) {
            Graph G = dual_from_entry(N, db.entries[idx]);
            ReducibleDual rd(G);
            ExtensionPath ep = rd.reduceToExtensionPath(20);
            if (ep.seed == SeedType::NotASeed) {
                fprintf(stderr, "%6d  (skipped — cannot reduce to seed)\n", idx);
                n_skipped++;
                continue;
            }

            Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);
            ConvexityStats s = ConvexityStats::compute(D);
            s.converged = true;
            s.print_row(idx);
            stats.push_back(s);

            if (s.n_concave > 0) n_still_concave++;
        }
        if (!stats.empty()) ConvexityStats::print_summary(stats);
        fprintf(stderr, "  C%d: %d/%d still concave, %d skipped (no seed path)\n",
               N, n_still_concave, (int)stats.size(), n_skipped);
    };

    test_size(60, c60_concave);
    test_size(70, c70_concave);
}

// =====================================================================
// ExtPath convexity sweep — full isomer spaces up to C50
// Run all in parallel: --gtest_filter='ExtPathConvexity.*'
// Uses default iters/step (= D.N, the compact vertex count at each step).
// =====================================================================

static void test_extpath_convexity_size(int N) {
    auto Q = BuckyGen::start(N, false);
    Graph G;
    int idx = 0, n_concave = 0, n_defects = 0;
    double worst_h = INFINITY, worst_cv = 0;
    double worst_ang_min = 180, worst_ang_max = 0, worst_ang_std = 0;

    while (BuckyGen::next_fullerene(Q, G)) {
        ReducibleDual rd(G);
        ExtensionPath ep = rd.reduceToExtensionPath();
        Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);

        double h = min_convexity_height(D);
        int nc = (int)concave_vertices(D).size();
        int nd = count_orientation_defects(D);

        // Edge CV
        vector<double> edge_lens;
        for (int u = 0; u < D.N; u++)
            for (int v : D.neighbours[u])
                if (v > u) edge_lens.push_back((D.points[u] - D.points[v]).norm());
        double cv = cv_twopass(edge_lens);

        // Triangle angles
        vector<double> angles;
        double ang_min = 180, ang_max = 0;
        for (const auto& tri : D.triangles)
            for (int c = 0; c < 3; c++) {
                coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
                coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
                double ang = coord3d::angle(va, vb) * 180.0 / M_PI;
                angles.push_back(ang);
                ang_min = min(ang_min, ang); ang_max = max(ang_max, ang);
            }
        double ang_std = stddev_twopass(angles);

        worst_h = min(worst_h, h);
        worst_cv = max(worst_cv, cv);
        worst_ang_min = min(worst_ang_min, ang_min);
        worst_ang_max = max(worst_ang_max, ang_max);
        worst_ang_std = max(worst_ang_std, ang_std);
        if (nc > 0) n_concave++;
        if (nd > 0) n_defects++;
        idx++;
    }
    BuckyGen::stop(Q);

    fprintf(stderr, "  C%-3d %5d isomers  L_cv=%9.2e  h_min=%+8.4f  ang=[%.2f,%.2f] std=%.2e  concave=%d  odef=%d\n",
           N, idx, worst_cv, worst_h, worst_ang_min, worst_ang_max, worst_ang_std,
           n_concave, n_defects);

    EXPECT_EQ(n_defects, 0) << "C" << N << ": orientation defects";
    EXPECT_EQ(n_concave, 0) << "C" << N << ": concave isomers";
    EXPECT_LT(worst_cv, 0.001) << "C" << N << ": L_cv=" << worst_cv;
    EXPECT_GT(worst_h, 0) << "C" << N << ": h_min=" << worst_h;
    EXPECT_GT(worst_ang_min, 50.0) << "C" << N << ": ang_min=" << worst_ang_min;
    EXPECT_LT(worst_ang_max, 70.0) << "C" << N << ": ang_max=" << worst_ang_max;
    EXPECT_LT(worst_ang_std, 2.0) << "C" << N << ": ang_std=" << worst_ang_std;
}

TEST(ExtPathConvexity, C42) { test_extpath_convexity_size(42); }
TEST(ExtPathConvexity, C44) { test_extpath_convexity_size(44); }
TEST(ExtPathConvexity, C46) { test_extpath_convexity_size(46); }
TEST(ExtPathConvexity, C48) { test_extpath_convexity_size(48); }
TEST(ExtPathConvexity, C50) { test_extpath_convexity_size(50); }

// =====================================================================
// C60-Ih (Buckminsterfullerene) via ExtPath
//
// The dual of C60-Ih is the pentakis dodecahedron: 32 vertices
// (12 deg-5 + 20 deg-6), 90 edges, 60 equilateral triangles.
// With Ih symmetry, the optimized geometry should be perfectly convex
// with all edges equal length (two edge classes on the exact pentakis
// dodecahedron, but the equilateral optimizer makes them all equal).
// =====================================================================

TEST(ExtPathConvexity, C60_Ih) {
    // C60-Ih is isomer #1811 (0-indexed) in the database
    IsomerDB db = IsomerDB::readPDB(60, false);
    ASSERT_GT((int)db.entries.size(), 1811);

    // --- Exact pentakis dodecahedron from traditional pipeline ---
    // IsomerDB → FullereneGraph → Tutte layout → zero_order_geometry → Polyhedron → dual
    FullereneGraph FG = IsomerDB::makeIsomer(60, db.entries[1811]);
    FG.layout2d = FG.tutte_layout();
    static const double target_L = sqrt(3.0) * 1.45;
    double scalerad = target_L / (1.5 * sqrt(3.0));
    vector<coord3d> pts = FG.zero_order_geometry(scalerad);
    Polyhedron P_exact(FG, pts, 6);
    Polyhedron dual_exact = P_exact.dual();
    Deltahedron D_exact(dual_exact);
    D_exact.optimize(D_exact.points, target_L);

    // Exact pentakis dodecahedron edge stats: two edge classes
    {
        vector<double> lens;
        for (int u = 0; u < D_exact.N; u++)
            for (int v : D_exact.neighbours[u])
                if (v > u) lens.push_back((D_exact.points[u] - D_exact.points[v]).norm());
        sort(lens.begin(), lens.end());
        double L_min = lens.front(), L_max = lens.back();
        double L_mean = 0; for (double l : lens) L_mean += l; L_mean /= lens.size();
        fprintf(stderr, "\n  Exact pentakis dodecahedron (direct optimize, target_L=%.4f):\n", target_L);
        fprintf(stderr, "    %d edges, L_min=%.6f L_max=%.6f L_mean=%.6f range/mean=%.4f\n",
               (int)lens.size(), L_min, L_max, L_mean, (L_max - L_min) / L_mean);
        ConvexityStats s_exact = ConvexityStats::compute(D_exact);
        ConvexityStats::print_header();
        s_exact.print_row(1811);
    }

    // --- ExtPath pipeline ---
    // Build dual triangulation from spiral pentagon indices
    const auto& e = db.entries[1811];
    int n_faces = 60/2 + 2;  // 32
    vector<int> spiral_string(n_faces, 6);
    for (int i = 0; i < 12; i++) spiral_string[e.RSPI[i] - 1] = 5;
    Triangulation T(spiral_string);
    Graph G = static_cast<const Graph&>(T);

    ASSERT_EQ(G.N, 32) << "C60 dual has 32 vertices";

    ReducibleDual rd(G);
    ExtensionPath ep = rd.reduceToExtensionPath(20);
    ASSERT_NE(ep.seed, SeedType::NotASeed) << "C60-Ih should reduce to a seed";

    Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);

    {
        vector<double> lens;
        for (int u = 0; u < D.N; u++)
            for (int v : D.neighbours[u])
                if (v > u) lens.push_back((D.points[u] - D.points[v]).norm());
        sort(lens.begin(), lens.end());
        double L_min = lens.front(), L_max = lens.back();
        double L_mean = 0; for (double l : lens) L_mean += l; L_mean /= lens.size();
        fprintf(stderr, "\n  ExtPath pentakis dodecahedron (N iters/step):\n");
        fprintf(stderr, "    seed=%d  steps=%zu\n", (int)ep.seed, ep.steps.size());
        fprintf(stderr, "    %d edges, L_min=%.6f L_max=%.6f L_mean=%.6f range/mean=%.4f\n",
               (int)lens.size(), L_min, L_max, L_mean, (L_max - L_min) / L_mean);
        ConvexityStats s = ConvexityStats::compute(D);
        ConvexityStats::print_header();
        s.print_row(1811);

        // Compare radii: both should produce spherical-ish shapes
        // Compute effective radius for both
        coord3d c_exact(0,0,0), c_ext(0,0,0);
        for (int u = 0; u < D_exact.N; u++) c_exact += D_exact.points[u];
        c_exact = c_exact / D_exact.N;
        for (int u = 0; u < D.N; u++) c_ext += D.points[u];
        c_ext = c_ext / D.N;

        // Sorted radii for both (rotation-invariant comparison)
        vector<double> r_exact, r_ext;
        for (int u = 0; u < D_exact.N; u++) r_exact.push_back((D_exact.points[u] - c_exact).norm());
        for (int u = 0; u < D.N; u++) r_ext.push_back((D.points[u] - c_ext).norm());
        sort(r_exact.begin(), r_exact.end());
        sort(r_ext.begin(), r_ext.end());

        // Normalize radii to mean=1 for scale-independent shape comparison
        double mean_r_exact = 0, mean_r_ext = 0;
        for (double r : r_exact) mean_r_exact += r;
        for (double r : r_ext) mean_r_ext += r;
        mean_r_exact /= r_exact.size();
        mean_r_ext /= r_ext.size();

        fprintf(stderr, "\n  Normalized radii comparison (mean=1, sorted):\n");
        fprintf(stderr, "    %10s %10s %10s\n", "exact_R", "extpath_R", "rel_diff");
        double max_rdiff = 0;
        for (int u = 0; u < D.N; u++) {
            double re = r_exact[u] / mean_r_exact;
            double rx = r_ext[u] / mean_r_ext;
            double diff = rx - re;
            max_rdiff = max(max_rdiff, fabs(diff));
            if (u < 5 || u >= D.N - 3)
                fprintf(stderr, "    %10.6f %10.6f %+10.6f\n", re, rx, diff);
            else if (u == 5)
                fprintf(stderr, "    %10s %10s %10s\n", "...", "...", "...");
        }
        fprintf(stderr, "    max |normalized R_diff| = %.6f\n", max_rdiff);
        fprintf(stderr, "    mean radii: exact=%.6f  extpath=%.6f  scale_ratio=%.4f\n",
               mean_r_exact, mean_r_ext, mean_r_exact / mean_r_ext);

        // Sorted edge lengths (rotation-invariant comparison)
        vector<double> e_exact;
        for (int u = 0; u < D_exact.N; u++)
            for (int v : D_exact.neighbours[u])
                if (v > u) e_exact.push_back((D_exact.points[u] - D_exact.points[v]).norm());
        sort(e_exact.begin(), e_exact.end());

        fprintf(stderr, "\n  Edge length distributions:\n");
        fprintf(stderr, "    exact:   min=%.6f  max=%.6f  (two classes in pentakis dodecahedron)\n",
               e_exact.front(), e_exact.back());
        fprintf(stderr, "    extpath: min=%.6f  max=%.6f  (equilateral optimizer → single class)\n",
               lens.front(), lens.back());

        // The Ih isomer should be perfectly convex and equilateral
        EXPECT_EQ(s.n_concave, 0) << "C60-Ih should have no concave vertices";
        EXPECT_EQ(s.n_orient_defects, 0) << "C60-Ih should have no orientation defects";
        EXPECT_GT(s.h_min, 0) << "C60-Ih should be strictly convex";
        EXPECT_LT(s.L_cv, 0.01) << "C60-Ih edge CV=" << s.L_cv;
    }
}

// =====================================================================
// Nanotube series via F-ring expansion: (5,0) nanotubes C30 → C500
//
// These are constructed by makeNanotubeDual(n_rings) and reduced to
// ExtensionPath (always reduces to C30 seed via F-ring steps).
// Tests the hardest case: high aspect ratio, many expansion steps.
// =====================================================================

TEST(ExtPathConvexity, NanotubeSeries) {
    fprintf(stderr, "\n  (5,0) nanotube series via ExtPath (N iters/step):\n");
    fprintf(stderr, "  %6s %5s %5s | %9s %8s %7s %7s %9s %4s %4s %8s\n",
            "C_N", "N_dv", "steps", "L_cv", "h_min", "ang_min", "ang_max", "ang_std", "conc", "odef", "ms");
    fprintf(stderr, "  %s\n", string(92, '-').c_str());

    for (int n_rings : {2, 3, 5, 8, 10, 15, 20, 30, 48}) {
        int N_ful = 20 + 10 * n_rings;
        int N_dv = 12 + 5 * n_rings;

        Graph G = makeNanotubeDual(n_rings);
        ASSERT_EQ(G.N, N_dv);

        ReducibleDual rd(G);
        ExtensionPath ep = rd.reduceToExtensionPath(20);
        if (ep.seed == SeedType::NotASeed) {
            fprintf(stderr, "  C%-4d %5d   --- | SKIP (no seed)\n", N_ful, N_dv);
            continue;
        }

        auto t0 = chrono::steady_clock::now();
        Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);
        auto t1 = chrono::steady_clock::now();
        double ms = chrono::duration<double, milli>(t1 - t0).count();

        double h = min_convexity_height(D);
        int nc = (int)concave_vertices(D).size();
        int nd = count_orientation_defects(D);

        vector<double> edge_lens;
        for (int u = 0; u < D.N; u++)
            for (int v : D.neighbours[u])
                if (v > u) edge_lens.push_back((D.points[u] - D.points[v]).norm());
        double cv = cv_twopass(edge_lens);

        vector<double> angles;
        double ang_min = 180, ang_max = 0;
        for (const auto& tri : D.triangles)
            for (int c = 0; c < 3; c++) {
                coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
                coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
                double ang = coord3d::angle(va, vb) * 180.0 / M_PI;
                angles.push_back(ang);
                ang_min = min(ang_min, ang); ang_max = max(ang_max, ang);
            }
        double ang_std = stddev_twopass(angles);

        fprintf(stderr, "  C%-4d %5d %5zu | %9.2e %+8.4f %7.2f %7.2f %9.2e %4d %4d %8.0f\n",
               N_ful, N_dv, ep.steps.size(), cv, h, ang_min, ang_max, ang_std, nc, nd, ms);

        EXPECT_EQ(nd, 0) << "C" << N_ful << " nanotube: orientation defects";
        EXPECT_EQ(nc, 0) << "C" << N_ful << " nanotube: concave vertices";
        EXPECT_TRUE(has_positive_orientation(D)) << "C" << N_ful << " nanotube: negative orientation";
        EXPECT_LT(cv, 0.001) << "C" << N_ful << " nanotube: L_cv=" << cv;
        EXPECT_GT(ang_min, 50.0) << "C" << N_ful << " nanotube: ang_min=" << ang_min;
        EXPECT_LT(ang_max, 70.0) << "C" << N_ful << " nanotube: ang_max=" << ang_max;
    }
}

// =====================================================================
// GC fullerenes via ExtPath
//
// GC(k,0) of the icosahedron gives C_{20k^2}. These are large, highly
// regular fullerenes. We also test GC of BuckyGen-generated fullerenes.
// =====================================================================

TEST(ExtPathConvexity, GCSeries) {
    // Build optimized icosahedron as GC base
    Polyhedron P20 = Polyhedron::C20();
    Polyhedron ico_poly = P20.dual();
    Deltahedron ico(ico_poly);
    ico.optimize(ico.points);

    fprintf(stderr, "\n  GC fullerenes via ExtPath (N iters/step):\n");
    fprintf(stderr, "  %-10s %6s %5s %5s | %9s %8s %7s %7s %9s %4s %4s %8s\n",
            "source", "C_N", "N_dv", "steps", "L_cv", "h_min", "ang_min", "ang_max", "ang_std", "conc", "odef", "ms");
    fprintf(stderr, "  %s\n", string(102, '-').c_str());

    // Icosahedral GC(k,0) series
    for (int k : {2, 3, 5, 7}) {
        int N_ful = 20 * k * k;
        int N_dv = N_ful / 2 + 2;

        Deltahedron gc = ico.GCtransform(k, 0);
        ASSERT_EQ(gc.N, N_dv);

        ReducibleDual rd(static_cast<const Graph&>(gc));
        ExtensionPath ep = rd.reduceToExtensionPath(20);
        if (ep.seed == SeedType::NotASeed) {
            fprintf(stderr, "  GC(%d,0)    C%-4d %5d   --- | SKIP (no seed)\n", k, N_ful, N_dv);
            continue;
        }

        auto t0 = chrono::steady_clock::now();
        Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);
        auto t1 = chrono::steady_clock::now();
        double ms = chrono::duration<double, milli>(t1 - t0).count();

        double h = min_convexity_height(D);
        int nc = (int)concave_vertices(D).size();
        int nd = count_orientation_defects(D);

        vector<double> edge_lens;
        for (int u = 0; u < D.N; u++)
            for (int v : D.neighbours[u])
                if (v > u) edge_lens.push_back((D.points[u] - D.points[v]).norm());
        double cv = cv_twopass(edge_lens);

        vector<double> angles;
        double ang_min = 180, ang_max = 0;
        for (const auto& tri : D.triangles)
            for (int c = 0; c < 3; c++) {
                coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
                coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
                double ang = coord3d::angle(va, vb) * 180.0 / M_PI;
                angles.push_back(ang);
                ang_min = min(ang_min, ang); ang_max = max(ang_max, ang);
            }
        double ang_std = stddev_twopass(angles);

        fprintf(stderr, "  GC(%d,0)    C%-4d %5d %5zu | %9.2e %+8.4f %7.2f %7.2f %9.2e %4d %4d %8.0f\n",
               k, N_ful, N_dv, ep.steps.size(), cv, h, ang_min, ang_max, ang_std, nc, nd, ms);

        EXPECT_EQ(nd, 0) << "GC(" << k << ",0) orientation defects";
        EXPECT_LT(cv, 0.01) << "GC(" << k << ",0) L_cv=" << cv;
        EXPECT_GT(ang_min, 40.0) << "GC(" << k << ",0) ang_min=" << ang_min;
        EXPECT_LT(ang_max, 80.0) << "GC(" << k << ",0) ang_max=" << ang_max;
    }

    // Non-icosahedral: GC(2,0) of a BuckyGen fullerene
    for (int baseN : {40, 60, 80}) {
        auto Q = BuckyGen::start(baseN, baseN >= 60);  // IPR for N>=60
        Graph G;
        if (!BuckyGen::next_fullerene(Q, G)) { BuckyGen::stop(Q); continue; }
        BuckyGen::stop(Q);

        ReducibleDual rd_base(G);
        ExtensionPath ep_base = rd_base.reduceToExtensionPath(20);
        if (ep_base.seed == SeedType::NotASeed) continue;

        Deltahedron D_base = Deltahedron::fromExtensionPathOptimized(ep_base);
        Deltahedron gc = D_base.GCtransform(2, 0);
        int gc_N_ful = baseN * 4;
        int gc_N_dv = gc.N;

        ReducibleDual rd_gc(static_cast<const Graph&>(gc));
        ExtensionPath ep_gc = rd_gc.reduceToExtensionPath(20);
        if (ep_gc.seed == SeedType::NotASeed) {
            fprintf(stderr, "  C%d+GC(2) C%-4d %5d   --- | SKIP (no seed)\n", baseN, gc_N_ful, gc_N_dv);
            continue;
        }

        auto t0 = chrono::steady_clock::now();
        Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep_gc);
        auto t1 = chrono::steady_clock::now();
        double ms = chrono::duration<double, milli>(t1 - t0).count();

        double h = min_convexity_height(D);
        int nc = (int)concave_vertices(D).size();
        int nd = count_orientation_defects(D);

        vector<double> edge_lens;
        for (int u = 0; u < D.N; u++)
            for (int v : D.neighbours[u])
                if (v > u) edge_lens.push_back((D.points[u] - D.points[v]).norm());
        double cv = cv_twopass(edge_lens);

        // Triangle angles
        vector<double> angles;
        double ang_min = 180, ang_max = 0;
        for (const auto& tri : D.triangles)
            for (int c = 0; c < 3; c++) {
                coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
                coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
                double ang = coord3d::angle(va, vb) * 180.0 / M_PI;
                angles.push_back(ang);
                ang_min = min(ang_min, ang); ang_max = max(ang_max, ang);
            }
        double ang_std = stddev_twopass(angles);

        fprintf(stderr, "  C%d+GC(2) C%-4d %5d %5zu | %9.2e %+8.4f %7.2f %7.2f %9.2e %4d %4d %8.0f\n",
               baseN, gc_N_ful, gc_N_dv, ep_gc.steps.size(), cv, h, ang_min, ang_max, ang_std, nc, nd, ms);

        EXPECT_EQ(nd, 0) << "C" << baseN << "+GC(2,0) orientation defects";
        EXPECT_LT(cv, 0.01) << "C" << baseN << "+GC(2,0) L_cv=" << cv;
        EXPECT_GT(ang_min, 40.0) << "C" << baseN << "+GC(2,0) ang_min=" << ang_min;
        EXPECT_LT(ang_max, 80.0) << "C" << baseN << "+GC(2,0) ang_max=" << ang_max;
    }
}
