// programs/eisenstein_paint/tikz_cell
//
// Emit TikZ pictures of the two iDT cells sharing one edge, plus the
// canonical strip for each of those cells' three CCW arcs (so the
// reader can compare the canonical strip frame against F's local
// frame after the affine transform).
//
// Usage:
//   eisenstein_paint_tikz_cell [N] [IPR] [idx] [arc|-1] [out_prefix]
//
// Defaults: N=200, IPR=1, idx=0, arc=-1 (auto-pick longest strip),
//           out_prefix=cell.  Emits <prefix>_f<id>.tex per cell, plus
//           <prefix>_f<id>_canon_<edge>.tex per CCW arc.

#include "fullerenes/eisenstein_paint_geometry.hh"
#include "fullerenes/delaunay_strip.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

namespace ep = eisenstein_paint;

static Triangulation load_isomer(int N, int idx, bool IPR) {
    BuckyGen::buckygen_queue Q = BuckyGen::start(N, IPR, false);
    Triangulation T;
    int i = 0;
    Graph G;
    while (BuckyGen::next_fullerene(Q, G)) {
        if (i == idx) { T = Triangulation(G); BuckyGen::stop(Q); return T; }
        ++i;
    }
    BuckyGen::stop(Q);
    std::fprintf(stderr, "load_isomer: C%d %s idx %d not found\n",
                 N, IPR ? "IPR" : "gen", idx);
    std::abort();
}

int main(int argc, char** argv) {
    const int   N      = (argc > 1) ? std::atoi(argv[1]) : 200;
    const bool  IPR    = (argc > 2) ? std::atoi(argv[2]) != 0 : true;
    const int   idx    = (argc > 3) ? std::atoi(argv[3]) : 0;
    int         arc    = (argc > 4) ? std::atoi(argv[4]) : -1;
    const char* prefix = (argc > 5) ? argv[5] : "cell";

    Triangulation T     = load_isomer(N, idx, IPR);
    ep::SortedDual S_d  = ep::sorted_dual(T);
    ep::DualPolytope P  = ep::realize_dual(S_d);
    const Triangulation& T_sorted = S_d.T;
    const DelaunayTriangulation& D = P.D;
    auto strips        = unfold_all_arc_strips(D, T_sorted);

    if (arc < 0) {
        int best_h = -1, best_sum = -1;
        for (int h = 0; h < (int)strips.size(); h += 2) {
            const Strip& S = strips[h];
            if (!S.ok) continue;
            const int s = S.v_position.first + S.v_position.second;
            if (s > best_sum) { best_sum = s; best_h = h; }
        }
        if (best_h < 0) {
            std::fprintf(stderr, "no live strips found\n");
            return 1;
        }
        arc = best_h;
        std::printf("auto-picked arc %d (longest a+b=%d)\n", arc, best_sum);
    }
    if (arc < 0 || arc >= (int)strips.size() || !strips[arc].ok) {
        std::fprintf(stderr, "arc %d invalid or dead\n", arc);
        return 1;
    }

    const int h_a = arc;
    const int h_b = h_a ^ 1;
    const int f_a = D.he_face[h_a];
    const int f_b = D.he_face[h_b];
    const Strip& S = strips[h_a];
    std::printf("C%d %s idx %d arc %d: u=%d v=%d v_pos=(%d,%d) norm^2=%d  "
                "cells=[%d, %d]\n",
                N, IPR ? "IPR" : "gen", idx, h_a, S.u, S.v,
                S.v_position.first, S.v_position.second,
                S.v_position.norm2(), f_a, f_b);

    for (int f : {f_a, f_b}) {
        if (f < 0 || f >= D.nf || D.f_he[f] < 0) {
            std::fprintf(stderr, "cell %d invalid or dead\n", f);
            continue;
        }
        ep::Cell F = ep::embed_cell(D, T_sorted, f);
        if (!F.ok) {
            std::fprintf(stderr, "embed_cell(%d) failed\n", f);
            continue;
        }
        char path[256];
        std::snprintf(path, sizeof path, "%s_f%d.tex", prefix, f);
        std::ofstream out(path);
        if (!out) {
            std::fprintf(stderr, "cannot open %s for writing\n", path);
            continue;
        }
        ep::dump_cell_tikz(F, T_sorted, out);
        std::printf("  cell %d: c0=%d c1=%d c2=%d  P1=(%d,%d) P2=(%d,%d) -> %s\n",
                    f, F.c0, F.c1, F.c2,
                    F.P1.first, F.P1.second, F.P2.first, F.P2.second,
                    path);

        // Canonical strips for each CCW arc.
        const int hh0 = D.f_he[f];
        const int hh1 = D.he_next[hh0];
        const int hh2 = D.he_next[hh1];
        const int  arcs[3] = { hh0, hh1, hh2 };
        const char* tags[3] = { "01", "12", "20" };
        for (int k = 0; k < 3; ++k) {
            const int hk = arcs[k];
            const Strip& Sk = strips[hk];
            char spath[256];
            std::snprintf(spath, sizeof spath, "%s_f%d_canon_%s.tex",
                          prefix, f, tags[k]);
            std::ofstream sout(spath);
            if (!sout) continue;
            dump_strip_tikz(Sk, T_sorted, sout);
            std::printf("    canonical strip edge %s arc %d "
                        "(u=%d v=%d v_pos=(%d,%d)) -> %s\n",
                        tags[k], hk, Sk.u, Sk.v,
                        Sk.v_position.first, Sk.v_position.second, spath);
        }
    }
    return 0;
}
