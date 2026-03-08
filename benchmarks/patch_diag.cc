// Per-step geometry diagnostic for specific fullerene isomers.
// Shows quality at each sub-phase within every expansion step,
// especially BEFORE full-graph CG (the "patched" phase).
//
// Usage:
//   patch_diag N idx            — buckygen isomer (0-based index)
//   patch_diag nanotube N       — (5,0) nanotube with N carbon atoms
//
// Build: cmake --build build --target patch_diag
// Run:   ./build/benchmarks/patch_diag 60 1217

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/buckinverse.hh"
#include <cstdio>
#include <cmath>
#include <vector>
#include <set>
#include <string>
#include <cassert>

using namespace std;
using namespace buckinverse;

// =====================================================================
// (5,0) nanotube dual graph builder
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
    int last = n_rings - 1;

    // Build oriented adjacency lists directly (CCW as seen from outside).
    neighbours_t adj(N, GRAPH_DMAX);

    // pole_N: CCW from outside = cn(0), cn(1), ..., cn(4)
    for (int i = 0; i < 5; i++) adj.push_back(pole_N, cn(i));

    // cn(i) (deg 5): CCW neighbors
    for (int i = 0; i < 5; i++)
        adj.assign_row(cn(i), {pole_N, cn(i+1), rng(0, i+1), rng(0, i), cn(i-1)});

    // rng(0, i): connects up to cn(i-1), cn(i) and down to next layer
    if (n_rings == 1) {
        for (int i = 0; i < 5; i++)
            adj.assign_row(rng(0, i), {cn(i-1), cn(i), rng(0, i+1), cs(i+1), cs(i), rng(0, i-1)});
    } else {
        for (int i = 0; i < 5; i++)
            adj.assign_row(rng(0, i), {cn(i-1), cn(i), rng(0, i+1), rng(1, i+1), rng(1, i), rng(0, i-1)});
    }

    // Interior rings: rng(j, i) for 1 <= j <= last-1
    for (int j = 1; j < last; j++)
        for (int i = 0; i < 5; i++)
            adj.assign_row(rng(j, i), {rng(j-1, i-1), rng(j-1, i), rng(j, i+1), rng(j+1, i+1), rng(j+1, i), rng(j, i-1)});

    // rng(last, i): connects down to cs
    if (n_rings >= 2)
        for (int i = 0; i < 5; i++)
            adj.assign_row(rng(last, i), {rng(last-1, i-1), rng(last-1, i), rng(last, i+1), cs(i+1), cs(i), rng(last, i-1)});

    // cs(i) (deg 5): connected up to rng(last,i-1) and rng(last,i)
    for (int i = 0; i < 5; i++)
        adj.assign_row(cs(i), {pole_S, cs(i-1), rng(last, i-1), rng(last, i), cs(i+1)});

    // pole_S (deg 5): CCW from outside (below) = cs(4), cs(3), ..., cs(0)
    for (int i = 4; i >= 0; i--) adj.push_back(pole_S, cs(i));

    Graph G(adj);
    assert(G.is_consistently_oriented());
    return G;
}

// =====================================================================
// Quality statistics for a Deltahedron
// =====================================================================
struct QStats {
    double edge_cv, edge_re;
    double ang_min, ang_max, ang_re;
    double h_min;
    int n_concave;
};

static QStats quality(const Deltahedron& D) {
    QStats q{};

    // Edge lengths
    double sum = 0, sum2 = 0; int ne = 0;
    double lmin = 1e30, lmax = 0;
    for (int u = 0; u < D.N; u++)
        for (int v : D.nbrs(u))
            if (v > u) {
                double l = (D.points[u] - D.points[v]).norm();
                sum += l; sum2 += l*l; ne++;
                lmin = min(lmin, l); lmax = max(lmax, l);
            }
    double L = sum / ne;
    q.edge_cv = sqrt(max(0.0, sum2/ne - L*L)) / L;
    q.edge_re = max(fabs(lmax - L), fabs(L - lmin)) / L;

    // Triangle angles
    q.ang_min = 180; q.ang_max = 0;
    for (const auto& tri : D.triangles)
        for (int c = 0; c < 3; c++) {
            coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
            coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
            double ang = coord3d::angle(va, vb) * 180.0 / M_PI;
            q.ang_min = min(q.ang_min, ang); q.ang_max = max(q.ang_max, ang);
        }
    q.ang_re = max(60.0 - q.ang_min, q.ang_max - 60.0) / 60.0;

    // Convexity
    q.h_min = 1e30; q.n_concave = 0;
    for (int v = 0; v < D.N; v++) {
        int d = (int)D.degree(v);
        if (d > 6) continue;
        coord3d cen(0,0,0);
        for (int i = 0; i < d; i++) cen += D.points[D.nbrs(v)[i]];
        cen /= (double)d;
        coord3d nf(0,0,0);
        for (int i = 0; i < d; i++) {
            coord3d e1 = D.points[D.nbrs(v)[i]] - D.points[v];
            coord3d e2 = D.points[D.nbrs(v)[(i+1)%d]] - D.points[v];
            nf += e1.cross(e2);
        }
        double nl = nf.norm();
        if (nl < 1e-15) continue;
        double h = (D.points[v] - cen).dot(nf / nl);
        if (h < q.h_min) q.h_min = h;
        if (h < -1e-6) q.n_concave++;
    }
    return q;
}

static void write_mol2(const Deltahedron& D, const string& path) {
    Polyhedron P(static_cast<const PlanarGraph&>(D), D.points);
    Polyhedron::to_file(P, path);
}

static void print_header() {
    fprintf(stderr, "  %4s %7s %4s | %9s %8s | %7s %7s %8s | %7s %4s\n",
            "step", "phase", "N",
            "edge_cv", "edge_re",
            "ang_min", "ang_max", "ang_re",
            "h_min", "conc");
    fprintf(stderr, "  %s\n", string(85, '-').c_str());
}

static void print_row(int step, const char* phase, const Deltahedron& D) {
    QStats q = quality(D);
    fprintf(stderr, "  %4d %-7s %4d | %9.2e %7.2f%% | %7.2f %7.2f %7.2f%% | %+7.4f %4d\n",
            step, phase, D.N,
            q.edge_cv, q.edge_re * 100,
            q.ang_min, q.ang_max, q.ang_re * 100,
            q.h_min, q.n_concave);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  %s N idx           — buckygen isomer (0-based index)\n", argv[0]);
        fprintf(stderr, "  %s nanotube N      — (5,0) nanotube with N carbon atoms\n", argv[0]);
        return 1;
    }

    Graph G;
    int N_fullerene;
    string label;
    bool is_nanotube = (string(argv[1]) == "nanotube");

    if (is_nanotube) {
        N_fullerene = atoi(argv[2]);
        int n_rings = (N_fullerene - 20) / 10;
        if (n_rings < 1 || 20 + 10 * n_rings != N_fullerene) {
            fprintf(stderr, "Nanotube N must be 30, 40, 50, ... (20 + 10k)\n");
            return 1;
        }
        G = makeNanotubeDual(n_rings);
        label = "C" + to_string(N_fullerene) + "_nanotube";
        fprintf(stderr, "C%d (5,0) nanotube: %d dual vertices, %d rings\n",
                N_fullerene, G.N, n_rings);
    } else {
        N_fullerene = atoi(argv[1]);
        int target_idx = atoi(argv[2]);
        auto Q = BuckyGen::start(N_fullerene, false);
        int idx = 0;
        bool found = false;
        while (BuckyGen::next_fullerene(Q, G)) {
            if (idx == target_idx) { found = true; break; }
            idx++;
        }
        BuckyGen::stop(Q);
        if (!found) { fprintf(stderr, "C%d #%d not found\n", N_fullerene, target_idx); return 1; }
        label = "C" + to_string(N_fullerene) + "_" + to_string(target_idx);
        fprintf(stderr, "C%d #%d: %d dual vertices\n", N_fullerene, target_idx, G.N);
    }

    // Reduce to extension path
    ReducibleDual rd(G);
    ExtensionPath ep = rd.reduceToExtensionPath(20);
    if (ep.seed == SeedType::NotASeed) {
        fprintf(stderr, "Failed to reduce to seed\n");
        return 1;
    }
    const char* seed_name = ep.seed == SeedType::C20 ? "C20" :
                            ep.seed == SeedType::C28 ? "C28" : "C30";
    fprintf(stderr, "Seed: %s, %zu expansion steps\n\n", seed_name, ep.steps.size());

    // Run with diagnostic callback
    fprintf(stderr, "=== Sub-phase quality at every expansion step ===\n");
    fprintf(stderr, "  Phases: placed = after strip insert + lift\n");
    fprintf(stderr, "          reflected = after concavity reflection\n");
    fprintf(stderr, "          patched = after patch Newton optimize (BEFORE full-graph CG)\n");
    fprintf(stderr, "          cg = after full-graph CG relaxation\n\n");
    print_header();

    Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep, 0, stderr,
        [&](int step, const char* phase, const Deltahedron& D_snap) {
            print_row(step, phase, D_snap);

            // Write mol2 for key phases: seed, patched (before CG), cg, final
            if (strcmp(phase, "seed") == 0 ||
                strcmp(phase, "patched") == 0 ||
                strcmp(phase, "cg") == 0 ||
                strcmp(phase, "final") == 0) {
                char path[256];
                snprintf(path, sizeof(path), "/tmp/%s_step%02d_%s.mol2",
                         label.c_str(), step, phase);
                write_mol2(D_snap, path);
            }
        });

    fprintf(stderr, "\nmol2 files in /tmp/%s_step*_{seed,patched,cg,final}.mol2\n", label.c_str());
    return 0;
}
