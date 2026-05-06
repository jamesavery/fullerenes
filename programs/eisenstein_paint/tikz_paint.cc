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

#include "fullerenes/eisenstein_paint.hh"
#include "fullerenes/delaunay_unfold.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <array>
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

    auto [T_sorted, D] = ep::prepare_inputs(T);
    auto cells         = ep::embed_all_cells(D, T_sorted);

    if (cell_id < 0) {
        for (int f = 0; f < D.nf; ++f)
            if (cells[f].ok) { cell_id = f; break; }
    }
    if (cell_id < 0 || cell_id >= D.nf || !cells[cell_id].ok) {
        std::fprintf(stderr,
            "cell %d not live.  Live cell ids: ", cell_id);
        for (int f = 0; f < D.nf; ++f)
            if (cells[f].ok) std::fprintf(stderr, "%d ", f);
        std::fprintf(stderr, "\n");
        return 1;
    }
    const ep::Cell& F = cells[cell_id];
    ep::LatticeMap F_lmap = ep::enumerate_cell_lattice(F, T_sorted);

    // Adapter: bundle (Cell, LatticeMap) into the lib-side CellPlacement
    // type the unfold tools consume.
    auto to_placement = [](const ep::Cell& F, const ep::LatticeMap& lmap) {
        return CellPlacement{
            F.cell_id, F.c0, F.c1, F.c2,
            F.P0, F.P1, F.P2,
            F.ok ? lmap.entries : std::vector<std::pair<Eisenstein, int>>{},
            F.ok };
    };
    CellPlacement F_place = to_placement(F, F_lmap);

    // 3 adjacent iDT cells across F's 3 edges.
    const int hh[3] = { D.f_he[cell_id],
                        D.he_next[D.f_he[cell_id]],
                        D.he_next[D.he_next[D.f_he[cell_id]]] };
    std::array<NeighbourCell, 3> neighbours;
    for (int k = 0; k < 3; ++k) {
        const int h_opp = hh[k] ^ 1;
        const int f_opp = D.he_face[h_opp];
        if (f_opp < 0 || f_opp >= D.nf || !cells[f_opp].ok) continue;
        const ep::Cell& Fn = cells[f_opp];
        ep::LatticeMap  Ln = ep::enumerate_cell_lattice(Fn, T_sorted);
        neighbours[k].placement = to_placement(Fn, Ln);
        neighbours[k].valid     = true;
    }

    std::printf("cell %d: c0=%d c1=%d c2=%d  P0=(%d,%d) P1=(%d,%d) P2=(%d,%d)  "
                "%zu lattice pts\n",
                F.cell_id, F.c0, F.c1, F.c2,
                F.P0.first, F.P0.second,
                F.P1.first, F.P1.second,
                F.P2.first, F.P2.second,
                F_lmap.entries.size());
    for (int k = 0; k < 3; ++k) {
        const auto& nd = neighbours[k];
        if (!nd.valid) { std::printf("  neighbour[%d]: (none)\n", k); continue; }
        std::printf("  neighbour[%d] cell %d: c0=%d c1=%d c2=%d  "
                    "%zu lattice pts\n",
                    k, nd.placement.cell_id,
                    nd.placement.c0, nd.placement.c1, nd.placement.c2,
                    nd.placement.lattice_points.size());
    }

    char path[256];
    std::snprintf(path, sizeof path, "%s_f%d_lmap.tex", prefix, cell_id);
    std::ofstream out1(path);
    ep::dump_lattice_map_tikz(F, F_lmap, T_sorted, out1);
    std::printf("wrote %s\n", path);

    std::snprintf(path, sizeof path, "%s_f%d_unfolded.tex", prefix, cell_id);
    std::ofstream out2(path);
    dump_neighbour_unfolding_tikz(F_place, neighbours, T_sorted, out2);
    std::printf("wrote %s\n", path);
    return 0;
}
