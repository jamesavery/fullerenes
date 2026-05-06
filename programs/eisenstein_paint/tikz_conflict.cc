// programs/eisenstein_paint/tikz_conflict
//
// Find an interior-ambiguity conflict (one hex vertex_id claimed by
// >= 2 cells as strictly interior in their lattice maps) and dump
// TikZ side-by-side LatticeMap pictures of the conflicting cells,
// with the disputed vertex highlighted.
//
// Usage:
//   eisenstein_paint_tikz_conflict N IPR idx [vertex_id|-1] [out_prefix]
//
// vertex_id = -1 -> pick the first interior-ambiguity conflict found.

#include "fullerenes/eisenstein_paint.hh"
#include "fullerenes/barycentric.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>

namespace ep = eisenstein_paint;

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr,
            "usage: eisenstein_paint_tikz_conflict N IPR idx "
            "[vertex_id|-1] [prefix]\n");
        return 1;
    }
    const int   N        = std::atoi(argv[1]);
    const bool  IPR      = std::atoi(argv[2]) != 0;
    const int   idx      = std::atoi(argv[3]);
    int         target_v = (argc > 4) ? std::atoi(argv[4]) : -1;
    const char* prefix   = (argc > 5) ? argv[5] : "conflict";

    BuckyGen::buckygen_queue Q = BuckyGen::start(N, IPR, false);
    Triangulation T;
    int i = 0;
    Graph G;
    while (BuckyGen::next_fullerene(Q, G)) {
        if (i == idx) { T = Triangulation(G); BuckyGen::stop(Q); break; }
        ++i;
    }
    if (T.N == 0) { std::fprintf(stderr, "isomer not found\n"); return 1; }

    auto [T_sorted, D] = ep::prepare_inputs(T);
    auto cells         = ep::embed_all_cells(D, T_sorted);

    struct Claim { int cell_id; Eisenstein pos; bool is_interior; };
    std::vector<std::vector<Claim>> claims(T_sorted.N);
    std::vector<ep::LatticeMap> lmaps(cells.size());
    for (size_t fi = 0; fi < cells.size(); ++fi) {
        if (!cells[fi].ok) continue;
        lmaps[fi] = ep::enumerate_cell_lattice(cells[fi], T_sorted);
        const ep::Cell& F = cells[fi];
        for (const auto& [p, v] : lmaps[fi].entries) {
            const IntBary bw = integer_barycentric(p, F.P0, F.P1, F.P2);
            const bool interior = (bw.n0 > 0 && bw.n1 > 0 && bw.n2 > 0);
            claims[v].push_back({(int)fi, p, interior});
        }
    }

    if (target_v < 0) {
        for (int v = 12; v < T_sorted.N; ++v) {
            int n_int = 0;
            for (const auto& c : claims[v]) if (c.is_interior) ++n_int;
            if (n_int >= 2) { target_v = v; break; }
        }
    }
    if (target_v < 0) {
        std::printf("no interior-ambiguity conflict found in C%d %s idx %d\n",
                    N, IPR ? "IPR" : "gen", idx);
        return 0;
    }

    std::vector<int> interior_cells;
    for (const auto& c : claims[target_v])
        if (c.is_interior) interior_cells.push_back(c.cell_id);

    std::printf("CONFLICT: hex %d claimed as interior by %zu cells:",
                target_v, interior_cells.size());
    for (int f : interior_cells) std::printf(" %d", f);
    std::printf("\n");
    for (int f : interior_cells) {
        const ep::Cell& F = cells[f];
        std::printf("  cell %d: c0=%d c1=%d c2=%d  P1=(%d,%d) P2=(%d,%d)\n",
                    f, F.c0, F.c1, F.c2,
                    F.P1.first, F.P1.second, F.P2.first, F.P2.second);
        for (const auto& c : claims[target_v]) {
            if (c.cell_id == f) {
                const IntBary bw = integer_barycentric(c.pos,
                                                       F.P0, F.P1, F.P2);
                std::printf("    hex %d at lattice (%d, %d) "
                            "bary=(%ld, %ld, %ld)/g=%ld\n",
                            target_v, c.pos.first, c.pos.second,
                            bw.n0, bw.n1, bw.n2, bw.denom);
            }
        }
    }

    for (int f : interior_cells) {
        char path[256];
        std::snprintf(path, sizeof path, "%s_v%d_f%d.tex",
                      prefix, target_v, f);
        std::ofstream out(path);
        ep::dump_lattice_map_tikz(cells[f], lmaps[f], T_sorted, out,
                                  0.7, target_v);
        std::printf("  wrote %s\n", path);
    }
    return 0;
}
