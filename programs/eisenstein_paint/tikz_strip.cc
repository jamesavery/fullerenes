// programs/eisenstein_paint/tikz_strip
//
// Emit a TikZ picture of one unfolded iDT-arc strip.
//
// Usage:
//   eisenstein_paint_tikz_strip [N] [IPR] [idx] [arc|-1] [out.tex]
//
// All arguments optional; defaults: N=200, IPR=1, idx=0, arc=-1
// (auto-pick longest), out.tex=strip.tex.

#include "fullerenes/eisenstein_paint.hh"
#include "fullerenes/delaunay_strip.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <cstdio>
#include <cstdlib>
#include <fstream>

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
    const int   N        = (argc > 1) ? std::atoi(argv[1]) : 200;
    const bool  IPR      = (argc > 2) ? std::atoi(argv[2]) != 0 : true;
    const int   idx      = (argc > 3) ? std::atoi(argv[3]) : 0;
    int         arc      = (argc > 4) ? std::atoi(argv[4]) : -1;
    const char* out_path = (argc > 5) ? argv[5] : "strip.tex";

    Triangulation T    = load_isomer(N, idx, IPR);
    auto [T_sorted, D] = ep::prepare_inputs(T);
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

    const Strip& S = strips[arc];
    std::printf("C%d %s idx %d arc %d: u=%d v=%d v_pos=(%d,%d) "
                "norm^2=%d strip-vertices=%zu\n",
                N, IPR ? "IPR" : "gen", idx, arc, S.u, S.v,
                S.v_position.first, S.v_position.second, S.v_position.norm2(),
                S.left.size() + S.right.size());

    std::ofstream out(out_path);
    if (!out) {
        std::fprintf(stderr, "cannot open %s for writing\n", out_path);
        return 1;
    }
    dump_strip_tikz(S, T_sorted, out);
    std::printf("wrote %s\n", out_path);
    return 0;
}
