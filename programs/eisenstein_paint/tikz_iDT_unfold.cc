// programs/eisenstein_paint/tikz_iDT_unfold
//
// Unfold the iDT triangulation into the Eisenstein plane via spanning-
// tree DFS and emit a TikZ picture + a debug report.
//
// Usage:
//   eisenstein_paint_tikz_iDT_unfold N IPR idx [start_cell|-1] [out_prefix] [seeds=K]
//
// seeds=K  try K different seed cells (default 20) and keep the
//          unfolding placing the most cells, ties by fewest cone tears.
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

    try {
        ep::SortedDual S_d = ep::sorted_dual(T);
        ep::DualPolytope Pl = ep::realize_dual(S_d);
        const Triangulation& T_sorted = S_d.T;
        const DelaunayTriangulation& D = Pl.D;

        // The per-cell charts unfold_iDT consumes, as the paint
        // pipeline's flat tables.
        ep::SurfaceParametrization param = ep::parametrize(D, S_d);

        LatticeUnfolding U = unfold_iDT(D, param.view(), start, seeds);

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
    } catch (const ep::PaintError& e) {
        std::fprintf(stderr, "tikz_iDT_unfold: paint stage %s: %s\n",
                     ep::code_name(e.code), e.what());
        return 1;
    }
    return 0;
}
