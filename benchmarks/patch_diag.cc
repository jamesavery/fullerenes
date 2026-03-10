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
#include <chrono>

using namespace std;
using namespace buckinverse;

// =====================================================================
// (5,0) nanotube dual graph builder
// =====================================================================
static Graph makeNanotubeDual(int n_rings) {
    assert(n_rings >= 1);
    int N = 12 + 5 * n_rings;

    auto mod5 = [](int i) -> int { return ((i % 5) + 5) % 5; };
    auto cn  = [&](int i) -> int { return 1 + mod5(i); };
    auto rng = [&](int j, int i) -> int { return 6 + 5*j + mod5(i); };
    auto cs  = [&](int i) -> int { return 6 + 5*n_rings + mod5(i); };
    int pole_N = 0;
    int pole_S = 11 + 5*n_rings;

    set<pair<int,int>> edges;
    auto add = [&](int u, int v) {
        int a = min(u,v), b = max(u,v);
        edges.insert(make_pair(a, b));
    };

    for (int i = 0; i < 5; i++) add(pole_N, cn(i));
    for (int i = 0; i < 5; i++) add(cn(i), cn((i+1)%5));
    for (int i = 0; i < 5; i++) {
        add(cn(i), rng(0, i));
        add(cn(i), rng(0, (i+1)%5));
    }
    for (int j = 0; j < n_rings; j++)
        for (int i = 0; i < 5; i++)
            add(rng(j, i), rng(j, (i+1)%5));
    for (int j = 0; j + 1 < n_rings; j++)
        for (int i = 0; i < 5; i++) {
            add(rng(j, i), rng(j+1, i));
            add(rng(j, i), rng(j+1, (i+1)%5));
        }
    for (int i = 0; i < 5; i++) {
        add(rng(n_rings-1, i), cs(i));
        add(rng(n_rings-1, i), cs((i+1)%5));
    }
    for (int i = 0; i < 5; i++) add(cs(i), cs((i+1)%5));
    for (int i = 0; i < 5; i++) add(pole_S, cs(i));

    assert((int)edges.size() == 3 * N - 6);

    neighbours_t adj(N);
    for (const auto& [u, v] : edges) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    Triangulation T(adj, false);
    return static_cast<const Graph&>(T);
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
        for (int v : D.neighbours[u])
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
        int d = (int)D.neighbours[v].size();
        if (d > 6) continue;
        coord3d cen(0,0,0);
        for (int i = 0; i < d; i++) cen += D.points[D.neighbours[v][i]];
        cen /= (double)d;
        coord3d nf(0,0,0);
        for (int i = 0; i < d; i++) {
            coord3d e1 = D.points[D.neighbours[v][i]] - D.points[v];
            coord3d e2 = D.points[D.neighbours[v][(i+1)%d]] - D.points[v];
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
        fprintf(stderr, "  %s N idx [--method CG|LBFGS|ST] [--step-tol T] [--final-tol T]\n", argv[0]);
        fprintf(stderr, "  %s nanotube N [--method CG|LBFGS|ST] [--step-tol T] [--final-tol T]\n", argv[0]);
        return 1;
    }

    // Parse optional flags from end of argv
    auto opt_method = OptMethod::CG;
    double step_tol = 1e-3, final_tol = 1e-5;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--method") && i+1 < argc) {
            i++;
            if (!strcmp(argv[i], "LBFGS")) opt_method = OptMethod::LBFGS;
            else if (!strcmp(argv[i], "ST")) opt_method = OptMethod::STEIHAUG;
        } else if (!strcmp(argv[i], "--step-tol") && i+1 < argc) {
            step_tol = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--final-tol") && i+1 < argc) {
            final_tol = atof(argv[++i]);
        }
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

    const char* method_name = opt_method == OptMethod::CG ? "CG" :
                              opt_method == OptMethod::LBFGS ? "LBFGS" : "ST";
    fprintf(stderr, "Method: %s, step_tol=%.0e, final_tol=%.0e\n\n", method_name, step_tol, final_tol);

    Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep, stderr,
        [&](int step, const char* phase, const Deltahedron& D_snap) {
            print_row(step, phase, D_snap);

            // Write mol2 for key phases: seed, patched (before relaxation), relaxed, final
            if (strcmp(phase, "seed") == 0 ||
                strcmp(phase, "patched") == 0 ||
                strcmp(phase, "relaxed") == 0 ||
                strcmp(phase, "final") == 0) {
                char path[256];
                snprintf(path, sizeof(path), "/tmp/%s_step%02d_%s.mol2",
                         label.c_str(), step, phase);
                write_mol2(D_snap, path);
            }
        },
        opt_method, step_tol, final_tol);

    fprintf(stderr, "\nResult: %s  ang=%.2e  conc=%d  gmax_L=%.2e\n",
            opt_result_name(D.final_opt_result), D.max_angle_relerr(),
            D.count_concave(), D.final_gmax_L);
    fprintf(stderr, "mol2 files in /tmp/%s_step*_{seed,patched,relaxed,final}.mol2\n", label.c_str());

    // Primitive cost measurement: run each method for a fixed budget, measure time and eval counts.
    fprintf(stderr, "\n=== Primitive cost ratios at N=%d ===\n", D.N);
    for (auto m : {OptMethod::CG, OptMethod::LBFGS, OptMethod::STEIHAUG}) {
        Deltahedron Dt = D;  // copy converged geometry
        Dt.opt_method = m;
        Dt.opt_k_flat = 0;
        long long budget = 20LL * Dt.N * Dt.N;
        auto t0 = chrono::high_resolution_clock::now();
        Dt.optimize(Dt.points, 0, 1e-15, {}, budget);  // unreachable tol, run to budget
        double ms = chrono::duration<double,milli>(chrono::high_resolution_clock::now() - t0).count();
        long long work_new = (long long)Dt.n_energy_evals + 2LL * Dt.n_grad_evals + 7LL * Dt.n_hv_evals;
        const char* mn = m == OptMethod::CG ? "CG" : m == OptMethod::LBFGS ? "LBFGS" : "ST";
        fprintf(stderr, "  %6s: %6.0f ms, %d iters, n_E=%d n_G=%d n_Hv=%d  work=%lld\n"
                        "          us/E=%.1f  us/G=%.1f  us/Hv=%.1f\n",
                mn, ms, Dt.iterations_used, Dt.n_energy_evals, Dt.n_grad_evals, Dt.n_hv_evals,
                work_new,
                Dt.n_energy_evals > 0 ? ms * 1000.0 / Dt.n_energy_evals : 0.0,
                Dt.n_grad_evals > 0 ? ms * 1000.0 / Dt.n_grad_evals : 0.0,
                Dt.n_hv_evals > 0 ? ms * 1000.0 / Dt.n_hv_evals : 0.0);
    }

    return 0;
}
