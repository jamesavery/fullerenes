// programs/eisenstein_paint/tikz_paint
//
// Emit TikZ pictures of one iDT cell's LatticeMap and its unfolded-
// with-3-neighbours view for inspection of the Eisenstein -> vertex_id
// mapping.
//
// Usage:
//   eisenstein_paint_tikz_paint N IPR idx cell_id [out_prefix]
//
// Outputs:
//   <prefix>_f<id>_lmap.tex
//   <prefix>_f<id>_unfolded.tex

#include "fullerenes/eisenstein_paint_geometry.hh"
#include "fullerenes/delaunay_unfold.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <cstdio>
#include <cstdlib>
#include <fstream>

namespace ep = eisenstein_paint;

int main(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr,
            "usage: eisenstein_paint_tikz_paint N IPR idx cell_id [prefix]\n");
        return 1;
    }
    const int   N       = std::atoi(argv[1]);
    const bool  IPR     = std::atoi(argv[2]) != 0;
    const int   idx     = std::atoi(argv[3]);
    int         cell_id = std::atoi(argv[4]);
    const char* prefix  = (argc > 5) ? argv[5] : "paint";

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

        // The flat chart tables both the prints and the unfolded view
        // consume.
        ep::SurfaceParametrization param = ep::parametrize(D, S_d);
        const ep::ParamTablesView V = param.view();

        if (cell_id < 0) {
            for (int f = 0; f < D.nf; ++f)
                if (V.cell_live(f)) { cell_id = f; break; }
        }
        if (cell_id < 0 || cell_id >= D.nf || !V.cell_live(cell_id)) {
            std::fprintf(stderr,
                "cell %d not live.  Live cell ids: ", cell_id);
            for (int f = 0; f < D.nf; ++f)
                if (V.cell_live(f)) std::fprintf(stderr, "%d ", f);
            std::fprintf(stderr, "\n");
            return 1;
        }

        const auto ids = V.corner_ids(cell_id);
        const auto P   = V.frame_points(cell_id);
        std::printf("cell %d: c0=%d c1=%d c2=%d  P0=(%d,%d) P1=(%d,%d) "
                    "P2=(%d,%d)  %zu lattice pts\n",
                    cell_id, ids[0], ids[1], ids[2],
                    P[0].first, P[0].second,
                    P[1].first, P[1].second,
                    P[2].first, P[2].second,
                    V.cell_entries(cell_id).size());
        int hh = D.f_he[cell_id];
        for (int k = 0; k < 3; ++k, hh = D.he_next[hh]) {
            const int f_opp = D.he_face[hh ^ 1];
            if (f_opp < 0 || !V.cell_live(f_opp)) {
                std::printf("  neighbour[%d]: (none)\n", k);
                continue;
            }
            const auto nids = V.corner_ids(f_opp);
            std::printf("  neighbour[%d] cell %d: c0=%d c1=%d c2=%d  "
                        "%zu lattice pts\n",
                        k, f_opp, nids[0], nids[1], nids[2],
                        V.cell_entries(f_opp).size());
        }

        // Strip-level view for the lmap picture only -- the tables
        // intentionally omit strip data, so this is the ONE
        // construction-vocabulary call.
        ep::Cell F = ep::embed_cell(D, T_sorted, cell_id);
        ep::LatticeMap F_lmap = ep::enumerate_cell_lattice(F, T_sorted);

        char path[256];
        std::snprintf(path, sizeof path, "%s_f%d_lmap.tex", prefix, cell_id);
        std::ofstream out1(path);
        ep::dump_lattice_map_tikz(F, F_lmap, T_sorted, out1);
        std::printf("wrote %s\n", path);

        std::snprintf(path, sizeof path, "%s_f%d_unfolded.tex", prefix, cell_id);
        std::ofstream out2(path);
        dump_neighbour_unfolding_tikz(D, V, cell_id, out2);
        std::printf("wrote %s\n", path);
    } catch (const ep::PaintError& e) {
        std::fprintf(stderr, "tikz_paint: paint stage %s: %s\n",
                     ep::code_name(e.code), e.what());
        return 1;
    }
    return 0;
}
