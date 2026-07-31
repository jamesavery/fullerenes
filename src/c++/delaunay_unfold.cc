#include "fullerenes/delaunay_unfold.hh"
#include "fullerenes/eisenstein_tikz.hh"
#include "fullerenes/eisenstein.hh"

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"

#include <algorithm>
#include <array>
#include <cstdarg>
#include <cstdio>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

using tikz::BBox;
using tikz::cart;
using eisenstein_paint::LatticePoint;
using eisenstein_paint::ParamTablesView;

namespace {

[[noreturn]] void unfold_throw(const char* fmt, ...) {
    char buf[512];
    va_list ap; va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    throw std::runtime_error(buf);
}

void require_parametrizes(const DelaunayTriangulation& D,
                          const ParamTablesView& V, const char* who) {
    if (V.nf != D.nf)
        unfold_throw("%s: tables/complex cell-count mismatch (V.nf=%d, "
                     "D.nf=%d) -- V must be the parametrization of this D",
                     who, V.nf, D.nf);
}

// An orientation-preserving lattice gluing of one directed segment onto
// another exists iff their edge vectors are ASSOCIATES (unit quotient):
// equal norm AND the aligning rotation z_to * conj(z_from) / N integral.
// Equal norm alone is NOT sufficient -- the smallest norm with several
// ideal classes is 49, where equal-norm vectors sit in unrelated orbits
// (rotation exists for neither branch) or in conjugate orbits (only a
// REFLECTED gluing exists, and a development never mirrors).  This is
// the non-throwing guard for isometry_from_segments, whose failures are
// deep-invariant logic_errors.
bool rotation_glues(Eisenstein z_from, Eisenstein z_to) {
    const long N = z_from.norm2();
    if (N == 0 || z_to.norm2() != N) return false;
    const Eisenstein num = z_to * z_from.complex_conj();
    return num.first % N == 0 && num.second % N == 0;
}

}  // namespace

// =====================================================================
// dump_neighbour_unfolding_tikz -- f + 3 unfolded neighbours.
// =====================================================================

void dump_neighbour_unfolding_tikz(const DelaunayTriangulation& D,
                                   const ParamTablesView& V,
                                   int f,
                                   std::ostream& os,
                                   double scale)
{
    // For each neighbour across f's three edges: compute the lattice
    // isometry from its frame to f's frame that glues the shared edge,
    // identified by cycle SLOT on both sides (multi-edge-sound).  Then
    // transform the neighbour's lattice points into f's frame and draw
    // them colour-coded.
    require_parametrizes(D, V, "dump_neighbour_unfolding_tikz");
    if (f < 0 || f >= D.nf || !V.cell_live(f))
        unfold_throw("dump_neighbour_unfolding_tikz: cell %d not charted", f);

    struct Transformed {
        int cell_id;
        Eisenstein P0, P1, P2;                          // corners in f's frame
        std::vector<std::pair<Eisenstein, int>> entries; // (pos in f's frame, vertex_id)
    };
    std::vector<Transformed> adj;

    const eisenstein_paint::ChartView chf = V.chart(f);
    const auto FP = chf.frame_points();

    int hh = D.f_he[f];
    for (int k = 0; k < 3; ++k, hh = D.he_next[hh]) {
        const int twin = hh ^ 1;
        const int g    = D.he_face[twin];
        if (g < 0 || g >= D.nf || !V.cell_live(g)) continue;
        // hh runs u -> v with u at f's corner slot k; the twin runs
        // v -> u with v at g's corner slot k_g.
        const int k_g = D.cycle_slot(twin);
        const eisenstein_paint::ChartView chg = V.chart(g);
        const auto GP = chg.frame_points();
        const Eisenstein Pa = FP[k],             Pb = FP[(k + 1) % 3];   // u, v in f
        const Eisenstein Qa = GP[(k_g + 1) % 3], Qb = GP[k_g];           // u, v in g
        // Both sides are genuine developments of the same edge, so the
        // rotation always exists here; the guard is defensive.
        if (!rotation_glues(Qb - Qa, Pb - Pa)) {
            std::fprintf(stderr,
                "dump_neighbour_unfolding_tikz: neighbour %d edge %d cannot "
                "be glued (norm or orbit mismatch)\n", g, k);
            continue;
        }
        const LatticeIsometry T = isometry_from_segments(Qa, Qb, Pa, Pb);
        Transformed t;
        t.cell_id = g;
        t.P0 = T.apply(GP[0]);
        t.P1 = T.apply(GP[1]);
        t.P2 = T.apply(GP[2]);
        for (const LatticePoint& e : chg.entries)
            t.entries.push_back({T.apply(e.pos()), e.vid});
        adj.push_back(std::move(t));
    }

    os << "\\documentclass[tikz,border=4pt]{standalone}\n"
          "\\usepackage{tikz}\n"
          "\\begin{document}\n"
          "\\begin{tikzpicture}[scale=" << scale
       << ", every node/.style={font=\\tiny}]\n";
    os << "% face " << f << " + 3 unfolded neighbours\n";

    BBox bb;
    bb.bump(FP[0]); bb.bump(FP[1]); bb.bump(FP[2]);
    for (const LatticePoint& e : chf.entries) bb.bump(e.pos());
    for (const auto& t : adj) {
        bb.bump(t.P0); bb.bump(t.P1); bb.bump(t.P2);
        for (const auto& [p, v] : t.entries) bb.bump(p);
    }
    const double margin = std::max(3.0, 0.2 * bb.span());
    bb.xmin -= margin; bb.xmax += margin;
    bb.ymin -= margin; bb.ymax += margin;
    tikz::emit_grid(os, bb);

    // f's lattice triangle in red; neighbours in lighter shades.
    const char* nb_tri_col[3] = { "blue!40", "green!40!black", "orange!60!black" };
    tikz::emit_lattice_triangle(os, FP[0], FP[1], FP[2]);
    for (size_t i = 0; i < adj.size(); ++i) {
        const auto& t = adj[i];
        auto [x0, y0] = cart(t.P0);
        auto [x1, y1] = cart(t.P1);
        auto [x2, y2] = cart(t.P2);
        os << "  \\draw[" << nb_tri_col[i] << ",line width=0.7pt] ("
           << x0 << "," << y0 << ") -- (" << x1 << "," << y1
           << ") -- (" << x2 << "," << y2 << ") -- cycle;\n";
    }

    auto face_col = [&](int tag) {
        switch (tag) {
            case 0: return "red!75!black";       // f itself
            case 1: return "blue!70!black";
            case 2: return "green!55!black";
            case 3: return "orange!75!black";
            default: return "black";
        }
    };

    // Draw EVERY lattice-point occurrence (not just first ownership) so
    // inconsistent placements show as MULTIPLE dots for the same id.
    std::vector<std::tuple<Eisenstein, int, int>> all_occ;
    auto cone_or_not = [&](int v) { return v < V.n_cones; };
    for (const LatticePoint& e : chf.entries)
        all_occ.push_back({e.pos(), e.vid, 0});
    for (size_t i = 0; i < adj.size(); ++i)
        for (const auto& [p, v] : adj[i].entries) all_occ.push_back({p, v, (int)i + 1});

    // Detect inconsistencies: same vertex_id at multiple distinct positions.
    std::unordered_map<int, std::set<std::pair<int,int>>> seen_positions;
    for (const auto& [p, v, tag] : all_occ)
        seen_positions[v].insert({p.first, p.second});

    os << "% lattice points (vertex_id shown; red ring = vertex mapped to "
          "≥ 2 distinct lattice positions across the 4 faces)\n";
    for (const auto& [p, v, tag] : all_occ) {
        auto [x, y] = cart(p);
        bool inconsistent = seen_positions[v].size() > 1;
        bool is_cone      = cone_or_not(v);
        const char* shape = is_cone ? "rectangle" : "circle";
        const char* fill  = face_col(tag);
        const char* draw  = inconsistent ? "draw=red,line width=1pt"
                                          : "draw=black";
        os << "  \\node[" << shape << "," << draw << ",fill=" << fill
           << ",inner sep=1.4pt,label=below right:{\\tiny " << v
           << "}] at (" << x << "," << y << ") {};\n";
    }

    os << "\\end{tikzpicture}\n\\end{document}\n";
}


// =====================================================================
// unfold_iDT -- spanning-tree DFS over cell adjacency.
// =====================================================================

namespace {

// One greedy DFS unfolding starting from a specific seed cell.  Each
// reachable charted cell is placed AT MOST once (via the half-edge
// that first reaches it and can be glued -- the header lists the skip
// conditions).  Cone deficits produce "tears": a cone placed by face A
// at position P_A may get placed by face B at a different P_B; we
// record all distinct positions per cone and count the tears.
LatticeUnfolding unfold_from_seed(const DelaunayTriangulation& D,
                                  const ParamTablesView& V,
                                  int seed_cell_id)
{
    LatticeUnfolding U;
    U.n_cones = V.n_cones;
    if (seed_cell_id < 0 || seed_cell_id >= D.nf || !V.cell_live(seed_cell_id))
        return U;

    // Per-cone FIRST position (the position recorded by the first cell
    // to place this cone in the DFS).  Used as the twin-arc anchor when
    // expanding into a new cell.
    std::map<int, Eisenstein> first_pos;
    auto record_cone = [&](int cone, Eisenstein p) {
        auto& all = U.cone_positions[cone];
        for (const auto& q : all) if (q == p) return;
        all.push_back(p);
    };

    // Place one cell: its chart mapped into the global frame by the
    // (orientation-preserving) lattice isometry g.
    auto place = [&](int cell_id, const eisenstein_paint::ChartView& ch,
                     int parent_id, const LatticeIsometry& g)
    {
        LatticeUnfolding::UnfoldedCell F;
        const auto FP = ch.frame_points();
        F.cell_id = cell_id;
        F.c0 = ch.corners[0]; F.c1 = ch.corners[1]; F.c2 = ch.corners[2];
        F.P0 = g.apply(FP[0]); F.P1 = g.apply(FP[1]); F.P2 = g.apply(FP[2]);
        F.parent_cell_id = parent_id;
        for (const LatticePoint& e : ch.entries)
            F.entries.push_back({g.apply(e.pos()), e.vid});

        for (auto [cone, p] : std::initializer_list<std::pair<int, Eisenstein>>{
                {F.c0, F.P0}, {F.c1, F.P1}, {F.c2, F.P2}})
        {
            auto it = first_pos.find(cone);
            if (it == first_pos.end()) first_pos[cone] = p;
            else if (it->second != p) ++U.n_tears;
            record_cone(cone, p);
        }
        U.cells.push_back(std::move(F));
        ++U.n_cells;
    };

    std::vector<bool> cell_placed(D.nf, false);

    // Seed: place at its local lattice positions (identity isometry).
    place(seed_cell_id, V.chart(seed_cell_id), /*parent=*/-1, LatticeIsometry{});
    cell_placed[seed_cell_id] = true;

    // DFS stack of half-edges pointing INTO unplaced cells from placed
    // ones (LIFO order).
    std::vector<int> stack;
    {
        int h0 = D.f_he[seed_cell_id];
        int h1 = D.he_next[h0];
        int h2 = D.he_next[h1];
        stack.push_back(h0 ^ 1);
        stack.push_back(h1 ^ 1);
        stack.push_back(h2 ^ 1);
    }

    while (!stack.empty()) {
        int h = stack.back();
        stack.pop_back();
        int fid = D.he_face[h];
        if (fid < 0 || fid >= D.nf) continue;
        if (cell_placed[fid]) continue;
        if (!V.cell_live(fid)) continue;

        // h runs u -> v inside fid, with u at fid's corner slot
        // D.cycle_slot(h) -- the slot pairing identifies the shared edge
        // by position in the face cycle, not by corner label.  The GLOBAL
        // anchors stay label-keyed on purpose: gluing to each cone's
        // first-sighting position is what defines the tear semantics.
        // Both anchors always exist -- the entering edge's endpoints are
        // corners of the already-placed twin cell, which recorded them.
        int u = D.he_origin[h];
        int v = D.he_origin[D.he_next[h]];
        auto itu = first_pos.find(u);
        auto itv = first_pos.find(v);
        if (itu == first_pos.end() || itv == first_pos.end())
            unfold_throw("unfold_iDT: entering edge %d -> %d of cell %d has "
                         "an unrecorded anchor -- placement invariant broken",
                         u, v, fid);
        Eisenstein Pu_g = itu->second, Pv_g = itv->second;
        const int  k_c  = D.cycle_slot(h);
        const eisenstein_paint::ChartView ch = V.chart(fid);
        const auto FP   = ch.frame_points();
        Eisenstein Pu_l = FP[k_c];
        Eisenstein Pv_l = FP[(k_c + 1) % 3];
        // Skip when no orientation-preserving gluing exists: after a tear
        // the label-keyed anchors can be ANY equal-norm vector -- in a
        // different ideal-class orbit (no lattice alignment at all) or in
        // the conjugate orbit (only a REFLECTED alignment, and a
        // development never mirrors; e.g. the other parallel copy of a
        // multi-edge, C128#0 cell 24).
        if (!rotation_glues(Pv_l - Pu_l, Pv_g - Pu_g)) continue;
        const LatticeIsometry g =
            isometry_from_segments(Pu_l, Pv_l, Pu_g, Pv_g);

        place(fid, ch, D.he_face[h ^ 1], g);
        cell_placed[fid] = true;

        int hh0 = D.f_he[fid];
        int hh1 = D.he_next[hh0];
        int hh2 = D.he_next[hh1];
        for (int hh : { hh0, hh1, hh2 }) {
            int other = hh ^ 1;
            int ofid  = D.he_face[other];
            if (ofid >= 0 && ofid < D.nf && !cell_placed[ofid]) {
                stack.push_back(other);
            }
        }
    }
    return U;
}

}  // namespace

LatticeUnfolding unfold_iDT(const DelaunayTriangulation& D,
                            const ParamTablesView& V,
                            int start_cell_id,
                            int dfs_max_seeds)
{
    require_parametrizes(D, V, "unfold_iDT");

    // Collect candidate seed cells.
    std::vector<int> seeds;
    if (start_cell_id >= 0 && start_cell_id < D.nf && V.cell_live(start_cell_id)) {
        seeds.push_back(start_cell_id);
    }
    if (dfs_max_seeds > (int)seeds.size() || seeds.empty()) {
        for (int f = 0; f < D.nf; ++f) {
            if (!V.cell_live(f)) continue;
            bool already = false;
            for (int s : seeds) if (s == f) { already = true; break; }
            if (already) continue;
            seeds.push_back(f);
            if ((int)seeds.size() >= std::max(1, dfs_max_seeds)) break;
        }
    }
    if (seeds.empty()) return LatticeUnfolding{};

    auto better = [](const LatticeUnfolding& a, const LatticeUnfolding& b) {
        if (a.n_cells != b.n_cells) return a.n_cells > b.n_cells;
        return a.n_tears < b.n_tears;
    };

    LatticeUnfolding best = unfold_from_seed(D, V, seeds[0]);
    for (size_t i = 1; i < seeds.size(); ++i) {
        LatticeUnfolding u = unfold_from_seed(D, V, seeds[i]);
        if (better(u, best)) best = std::move(u);
    }
    return best;
}

// =====================================================================
// dump_lattice_unfolding_tikz
// =====================================================================

void dump_lattice_unfolding_tikz(const LatticeUnfolding& U,
                                 const Triangulation& T_sorted,
                                 std::ostream& tikz,
                                 std::ostream& dbg,
                                 double scale)
{
    // Bounding box.
    BBox bb;
    for (const auto& F : U.cells) {
        bb.bump(F.P0); bb.bump(F.P1); bb.bump(F.P2);
        for (const auto& [p, v] : F.entries) bb.bump(p);
    }
    const double margin = std::max(3.0, 0.2 * bb.span());
    bb.xmin -= margin; bb.xmax += margin;
    bb.ymin -= margin; bb.ymax += margin;

    // Duplicates: vertex_id -> list of (global pos, cell_id).
    std::map<int, std::vector<std::pair<Eisenstein, int>>> occurrences;
    for (const auto& F : U.cells)
        for (const auto& [p, v] : F.entries) occurrences[v].push_back({p, F.cell_id});

    auto is_duplicate = [&](int v) -> bool {
        auto it = occurrences.find(v);
        if (it == occurrences.end()) return false;
        std::set<std::pair<int, int>> distinct;
        for (const auto& [p, f] : it->second) distinct.insert({p.first, p.second});
        return distinct.size() > 1;
    };

    tikz << "\\documentclass[tikz,border=4pt]{standalone}\n"
            "\\usepackage{tikz}\n"
            "\\begin{document}\n"
            "\\begin{tikzpicture}[scale=" << scale
         << ", every node/.style={font=\\tiny}]\n";
    tikz << "% iDT unfolding: " << U.n_cells << " faces placed, "
         << U.n_tears << " cone tears (spanning-tree DFS)\n";

    tikz::emit_grid(tikz, bb);

    // Draw cell triangles in muted blue.
    for (const auto& F : U.cells) {
        auto [x0, y0] = cart(F.P0);
        auto [x1, y1] = cart(F.P1);
        auto [x2, y2] = cart(F.P2);
        tikz << "  \\draw[blue!60!black,line width=0.5pt] ("
             << x0 << "," << y0 << ") -- (" << x1 << "," << y1
             << ") -- (" << x2 << "," << y2 << ") -- cycle;\n";
    }

    // Draw lattice points; duplicates circled in red.
    for (const auto& F : U.cells) {
        for (const auto& [p, v] : F.entries) {
            auto [x, y] = cart(p);
            bool cone = v < U.n_cones;
            bool dup  = is_duplicate(v);
            const char* shape = cone ? "rectangle" : "circle";
            const char* fill  = cone ? "red!75!black" : "blue!70!black";
            const char* drawspec = dup ? "draw=red,line width=1.1pt"
                                        : "draw=black!50";
            double sep = dup ? 2.2 : 1.3;
            tikz << "  \\node[" << shape << "," << drawspec << ",fill="
                 << fill << ",inner sep=" << sep
                 << "pt,label=below right:{\\tiny " << v
                 << "}] at (" << x << "," << y << ") {};\n";
        }
    }
    tikz << "\\end{tikzpicture}\n\\end{document}\n";

    // Debug report.
    dbg << "iDT unfolding report (spanning-tree DFS)\n";
    dbg << "  placed iDT faces: " << U.n_cells << "\n";
    dbg << "  cone tears (cone placed at >1 position): "
        << U.n_tears << "\n";
    dbg << "  cones covered: " << U.cone_positions.size() << "\n";
    dbg << "\n-- cones (positions may differ per face due to tears) --\n";
    for (const auto& [c, ps] : U.cone_positions) {
        dbg << "    cone " << c << ":";
        for (const auto& p : ps)
            dbg << " (" << p.first << "," << p.second << ")";
        dbg << "\n";
    }
    dbg << "\n-- placed faces (parent = -1 means seed) --\n";
    for (const auto& F : U.cells) {
        dbg << "    face " << F.cell_id << " parent=" << F.parent_cell_id
            << "  c0=" << F.c0 << " c1=" << F.c1 << " c2=" << F.c2
            << "  P0=(" << F.P0.first << "," << F.P0.second << ")"
            << "  P1=(" << F.P1.first << "," << F.P1.second << ")"
            << "  P2=(" << F.P2.first << "," << F.P2.second << ")  "
            << F.entries.size() << " lattice pts\n";
    }

    // Duplicate non-cone vertices (cones tear by design; anything ELSE
    // at two distinct positions is a placement inconsistency).
    dbg << "\n-- DUPLICATE non-cone vertices (at >=2 distinct global positions) --\n";
    int n_dup = 0;
    for (const auto& [v, occ] : occurrences) {
        if (v < U.n_cones) continue;
        std::set<std::pair<int, int>> distinct;
        for (const auto& [p, f] : occ) distinct.insert({p.first, p.second});
        if (distinct.size() > 1) {
            ++n_dup;
            dbg << "    vertex " << v << " at " << distinct.size()
                << " positions:\n";
            for (const auto& [p, f] : occ) {
                dbg << "       face " << f << " pos (" << p.first
                    << ", " << p.second << ")\n";
            }
        }
    }
    if (n_dup == 0) dbg << "    (none)\n";

    // Missing vertices (never placed by any placed cell -- e.g. the
    // entire subtree carrying them was stranded by a tear).
    dbg << "\n-- MISSING vertices (T_sorted vertex never placed) --\n";
    int n_miss = 0;
    for (int v = 0; v < T_sorted.N; ++v) {
        if (occurrences.find(v) == occurrences.end()) {
            ++n_miss;
            dbg << "    vertex " << v << "\n";
        }
    }
    if (n_miss == 0) dbg << "    (none)\n";

    dbg << "\nsummary: " << n_dup << " duplicates, " << n_miss
        << " missing; unfolding placed " << U.n_cells << " faces\n";
}
