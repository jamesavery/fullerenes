// Dump per-step mol2 files for a specific C60 isomer.
// Strategy: build partial ExtensionPaths (first k steps) and call fromExtensionPathOptimized
// on each, so we see the actual optimized geometry at each intermediate size.
//
// Build: cd /home/avery/work/fullerenes/build2 && g++ -std=c++20 -O2 -I../include -I../claude-projects/buckinverse /tmp/step_mol2.cc -Lsrc/c++ -lfullerenes -o /tmp/step_mol2
// Run:   cd /home/avery/work/fullerenes/build2 && LD_LIBRARY_PATH=src/c++ /tmp/step_mol2

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/deltahedron.hh"
#include "fullerenes/buckinverse.hh"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace std;
using namespace buckinverse;

static void write_mol2(const Deltahedron& D, const char* path) {
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); return; }
    int ne = 0;
    for (int u = 0; u < D.N; u++)
        for (int v : D.nbrs(u))
            if (v > u) ne++;
    fprintf(f, "@<TRIPOS>MOLECULE\nDeltahedron\n%d %d 0 0 0\nSMALL\nNO_CHARGES\n\n", D.N, ne);
    fprintf(f, "@<TRIPOS>ATOM\n");
    for (int u = 0; u < D.N; u++) {
        const char* type = ((int)D.degree(u) == 5) ? "N" : "C";
        fprintf(f, "%d %s%d %f %f %f %s 1 UNK 0\n",
                u+1, type, u, D.points[u][0], D.points[u][1], D.points[u][2], type);
    }
    fprintf(f, "@<TRIPOS>BOND\n");
    int bid = 1;
    for (int u = 0; u < D.N; u++)
        for (int v : D.nbrs(u))
            if (v > u) fprintf(f, "%d %d %d 1\n", bid++, u+1, v+1);
    fclose(f);
}

static void print_stats(const Deltahedron& D, const char* label) {
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
    double edge_cv = sqrt(max(0.0, sum2/ne - L*L)) / L;
    double edge_re = max(fabs(lmax - L), fabs(L - lmin)) / L;

    double amin = 180, amax = 0;
    for (const auto& tri : D.triangles)
        for (int c = 0; c < 3; c++) {
            coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
            coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
            double ang = coord3d::angle(va, vb) * 180.0 / M_PI;
            amin = min(amin, ang); amax = max(amax, ang);
        }
    double ang_re = max(60.0 - amin, amax - 60.0) / 60.0;

    double hmin = 1e30; int nconc = 0;
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
        if (h < hmin) hmin = h;
        if (h < -1e-6) nconc++;
    }

    fprintf(stderr, "  %-20s N=%2d edge_cv=%.3e edge_re=%.2f%% ang=[%.2f,%.2f] ang_re=%.2f%% h_min=%+.4f conc=%d gmax*L=%.2e\n",
            label, D.N, edge_cv, edge_re*100, amin, amax, ang_re*100, hmin, nconc, D.final_gmax_L);
}

int main() {
    int N = 60;
    int target_idx = 1217;

    // Find target isomer via buckygen
    auto Q = BuckyGen::start(N, false);
    Graph G;
    int idx = 0;
    bool found = false;
    while (BuckyGen::next_fullerene(Q, G)) {
        if (idx == target_idx) { found = true; break; }
        idx++;
    }
    BuckyGen::stop(Q);
    if (!found) { fprintf(stderr, "Not found\n"); return 1; }

    // Reduce to extension path
    ReducibleDual rd(G);
    ExtensionPath ep = rd.reduceToExtensionPath(20);
    fprintf(stderr, "C%d idx=%d seed=%d steps=%d\n", N, target_idx, (int)ep.seed, (int)ep.steps.size());

    // Build full result with logging
    fprintf(stderr, "\n=== Full pipeline (with CG log) ===\n");
    Deltahedron D_full = Deltahedron::fromExtensionPathOptimized(ep, 0, stderr);
    write_mol2(D_full, "/tmp/c60_1217_final.mol2");
    print_stats(D_full, "final");

    // Build partial extension paths: after step 0, 1, 2, ..., n-1
    fprintf(stderr, "\n=== Step-by-step intermediates ===\n");
    for (int k = 0; k <= (int)ep.steps.size(); k++) {
        ExtensionPath partial = ep;
        partial.steps.resize(k);
        // Recompute full_N for partial path
        int partial_N = 0;
        switch (ep.seed) {
            case SeedType::C20: partial_N = 12; break;
            case SeedType::C28: partial_N = 16; break;
            case SeedType::C30: partial_N = 17; break;
            default: break;
        }
        for (int s = 0; s < k; s++)
            partial_N += (int)partial.steps[s].strip.size();
        // But full_N in the ExtensionPath is the allocation size (vertex IDs go up to full_N-1)
        // so we keep it as-is; fromExtensionPathOptimized only uses vertices that exist.
        // Actually the issue is that the vertex IDs in later steps reference vertices
        // that don't exist yet. But we've truncated steps, so that's fine.

        Deltahedron D = Deltahedron::fromExtensionPathOptimized(partial);
        char path[256];
        snprintf(path, sizeof(path), "/tmp/c60_%d_step%02d.mol2", target_idx, k);
        write_mol2(D, path);
        char label[64];
        snprintf(label, sizeof(label), "step %d (N=%d)", k, D.N);
        print_stats(D, label);
    }

    return 0;
}
