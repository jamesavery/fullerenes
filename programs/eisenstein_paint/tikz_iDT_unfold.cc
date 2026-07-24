// programs/eisenstein_paint/tikz_iDT_unfold
//
// Unfold the iDT triangulation into the Eisenstein plane via spanning-
// tree DFS and emit a TikZ picture + a debug report.
//
// Usage:
//   eisenstein_paint_tikz_iDT_unfold N IPR idx [start_cell|-1] [out_prefix] [seeds=K]
//
// seeds=K  try K different seed cells (default 20 = all cells) and keep
//          the unfolding with the fewest cone tears.
//
// Outputs:
//   <prefix>.tex    TikZ standalone of the spanning-tree unfolding
//   <prefix>.txt    debug report (placed cells, tears, missing hexes)

#include "fullerenes/eisenstein_paint_geometry.hh"
#include "fullerenes/delaunay_unfold.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

namespace ep = eisenstein_paint;

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr,
            "usage: eisenstein_paint_tikz_iDT_unfold N IPR idx "
            "[start_cell|-1] [prefix] [seeds=K]\n");
        return 1;
    }
    const int   N      = std::atoi(argv[1]);
    const bool  IPR    = std::atoi(argv[2]) != 0;
    const int   idx    = std::atoi(argv[3]);
    const int   start  = (argc > 4) ? std::atoi(argv[4]) : -1;
    const char* prefix = (argc > 5) ? argv[5] : "idt_unfold";
    int         seeds  = 20;
    for (int ai = 6; ai < argc; ++ai)
        if (std::strncmp(argv[ai], "seeds=", 6) == 0)
            seeds = std::atoi(argv[ai] + 6);

    BuckyGen::buckygen_queue Q = BuckyGen::start(N, IPR, false);
    Triangulation T;
    int i = 0;
    Graph G;
    while (BuckyGen::next_fullerene(Q, G)) {
        if (i == idx) { T = Triangulation(G); BuckyGen::stop(Q); break; }
        ++i;
    }
    if (T.N == 0) { std::fprintf(stderr, "isomer not found\n"); return 1; }

    ep::SortedDual S_d = ep::sorted_dual(T);
    ep::DualPolytope Pl = ep::realize_dual(S_d);
    const Triangulation& T_sorted = S_d.T;
    const DelaunayTriangulation& D = Pl.D;
    auto cells_em      = ep::embed_all_cells(D, T_sorted);
    std::vector<ep::LatticeMap> lmaps(cells_em.size());
    for (size_t fi = 0; fi < cells_em.size(); ++fi)
        if (cells_em[fi].ok)
            lmaps[fi] = ep::enumerate_cell_lattice(cells_em[fi], T_sorted);

    // Adapter: bundle (Cell, LatticeMap) into the lib-side CellPlacement
    // type unfold_iDT consumes.
    std::vector<CellPlacement> cells(cells_em.size());
    for (size_t fi = 0; fi < cells_em.size(); ++fi) {
        const ep::Cell& F = cells_em[fi];
        cells[fi] = { F.cell_id, F.corners[0], F.corners[1], F.corners[2],
                      F.P[0], F.P[1], F.P[2],
                      F.ok ? lmaps[fi].entries
                           : std::vector<std::pair<Eisenstein, int>>{},
                      F.ok };
    }

    LatticeUnfolding U = unfold_iDT(D, cells, start, seeds);

    char tikz_path[256], txt_path[256];
    std::snprintf(tikz_path, sizeof tikz_path, "%s.tex", prefix);
    std::snprintf(txt_path,  sizeof txt_path,  "%s.txt", prefix);
    std::ofstream tikz(tikz_path);
    std::ofstream txt(txt_path);
    dump_lattice_unfolding_tikz(U, T_sorted, tikz, txt);

    std::printf("C%d %s idx %d: placed %d iDT cells, %d cone tears, "
                "%zu cones covered\n",
                N, IPR ? "IPR" : "gen", idx, U.n_cells, U.n_tears,
                U.cone_positions.size());
    std::printf("  wrote %s\n  wrote %s\n", tikz_path, txt_path);
    return 0;
}
