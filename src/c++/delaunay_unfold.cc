#include "fullerenes/delaunay_unfold.hh"
#include "fullerenes/eisenstein_tikz.hh"
#include "fullerenes/eisenstein_ops.hh"

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

namespace {

[[noreturn]] void unfold_throw(const char* fmt, ...) {
    char buf[512];
    va_list ap; va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    throw std::runtime_error(buf);
}

}  // namespace

// =====================================================================
// dump_neighbour_unfolding_tikz -- F + 3 unfolded neighbours.
// =====================================================================

void dump_neighbour_unfolding_tikz(const CellPlacement& F,
                                   const std::array<NeighbourCell, 3>& neighbours,
                                   const Triangulation& T_sorted,
                                   std::ostream& os,
                                   double scale)
{
    // For each neighbour: compute an affine mapping from its frame to
    // F's frame such that the shared edge's two cones coincide.  Then
    // transform the neighbour's lattice points into F's frame and draw
    // them colour-coded.

    struct Transformed {
        int cell_id;
        Eisenstein P0, P1, P2;                          // corners in F's frame
        std::vector<std::pair<Eisenstein, int>> entries; // (pos in F's frame, vertex_id)
    };
    std::vector<Transformed> adj;

    // F's edges, by index: 0 = c0-c1, 1 = c1-c2, 2 = c2-c0.
    int Fcor[3][2] = { { F.c0, F.c1 }, { F.c1, F.c2 }, { F.c2, F.c0 } };
    Eisenstein FP[3][2] = { { F.P0, F.P1 }, { F.P1, F.P2 }, { F.P2, F.P0 } };

    auto find_corner_pos = [&](const CellPlacement& N, int cone) -> Eisenstein {
        if (cone == N.c0) return N.P0;
        if (cone == N.c1) return N.P1;
        if (cone == N.c2) return N.P2;
        return Eisenstein(0, 0);
    };

    for (int k = 0; k < 3; ++k) {
        const NeighbourCell& nd = neighbours[k];
        if (!nd.valid || !nd.placement.ok) continue;
        const CellPlacement& N = nd.placement;
        int c_a = Fcor[k][0], c_b = Fcor[k][1];
        Eisenstein Pa = FP[k][0], Pb = FP[k][1];
        Eisenstein Qa = find_corner_pos(N, c_a);
        Eisenstein Qb = find_corner_pos(N, c_b);
        Eisenstein z_src = Qb - Qa;
        Eisenstein z_dst = Pb - Pa;
        if (z_src.norm2() != z_dst.norm2()) {
            std::fprintf(stderr,
                "dump_neighbour_unfolding_tikz: neighbour %d edge %d norm "
                "mismatch (src %d vs dst %d)\n",
                N.cell_id, k, z_src.norm2(), z_dst.norm2());
            continue;
        }
        D6Affine T = align(z_src, z_dst);
        auto apply = [&](Eisenstein p) {
            Eisenstein q = p - Qa;
            if (T.reflect) q = complex_conj(q);
            q = q * T.unit;
            return q + Pa;
        };
        Transformed t;
        t.cell_id = N.cell_id;
        t.P0 = apply(N.P0);
        t.P1 = apply(N.P1);
        t.P2 = apply(N.P2);
        for (const auto& [p, v] : N.lattice_points)
            t.entries.push_back({apply(p), v});
        adj.push_back(std::move(t));
    }

    os << "\\documentclass[tikz,border=4pt]{standalone}\n"
          "\\usepackage{tikz}\n"
          "\\begin{document}\n"
          "\\begin{tikzpicture}[scale=" << scale
       << ", every node/.style={font=\\tiny}]\n";
    os << "% face " << F.cell_id << " + 3 unfolded neighbours\n";

    BBox bb;
    bb.bump(F.P0); bb.bump(F.P1); bb.bump(F.P2);
    for (const auto& [p, v] : F.lattice_points) bb.bump(p);
    for (const auto& t : adj) {
        bb.bump(t.P0); bb.bump(t.P1); bb.bump(t.P2);
        for (const auto& [p, v] : t.entries) bb.bump(p);
    }
    const double margin = std::max(3.0, 0.2 * bb.span());
    bb.xmin -= margin; bb.xmax += margin;
    bb.ymin -= margin; bb.ymax += margin;
    tikz::emit_grid(os, bb);

    // F's lattice triangle in red; neighbours in lighter shades.
    const char* nb_tri_col[3] = { "blue!40", "green!40!black", "orange!60!black" };
    tikz::emit_lattice_triangle(os, F.P0, F.P1, F.P2);
    for (size_t i = 0; i < adj.size(); ++i) {
        const auto& t = adj[i];
        auto [x0, y0] = cart(t.P0);
        auto [x1, y1] = cart(t.P1);
        auto [x2, y2] = cart(t.P2);
        os << "  \\draw[" << nb_tri_col[i] << ",line width=0.7pt] ("
           << x0 << "," << y0 << ") -- (" << x1 << "," << y1
           << ") -- (" << x2 << "," << y2 << ") -- cycle;\n";
    }

    // Walk F's lattice points first; track {vertex_id -> (pos, owner_face_tag)}.
    // owner_face_tag: 0 = F, 1..3 = neighbours by index in adj.
    std::unordered_map<int, std::pair<Eisenstein, int>> ownership;
    auto mark = [&](Eisenstein p, int v, int tag) {
        auto it = ownership.find(v);
        if (it == ownership.end()) ownership[v] = {p, tag};
        // else: keep first; disagreement will show up as two draw-calls
        // at different lattice positions.
    };
    for (const auto& [p, v] : F.lattice_points) mark(p, v, 0);
    for (size_t i = 0; i < adj.size(); ++i)
        for (const auto& [p, v] : adj[i].entries) mark(p, v, (int)i + 1);

    auto face_col = [&](int tag) {
        switch (tag) {
            case 0: return "red!75!black";       // F itself
            case 1: return "blue!70!black";
            case 2: return "green!55!black";
            case 3: return "orange!75!black";
            default: return "black";
        }
    };

    // Draw EVERY lattice-point occurrence (not just first ownership) so
    // inconsistent placements show as MULTIPLE dots for the same id.
    std::vector<std::tuple<Eisenstein, int, int>> all_occ;
    auto cone_or_not = [&](int v) { return v < 12; };
    for (const auto& [p, v] : F.lattice_points) all_occ.push_back({p, v, 0});
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

// One greedy DFS unfolding starting from a specific seed cell.  Every
// reachable ok cell is placed exactly once (via the arc that first
// reaches it in the DFS traversal).  Cone deficits produce "tears": a
// cone placed by face A at position P_A may get placed by face B at a
// different P_B; we record all distinct positions per cone and count
// the tears.
LatticeUnfolding unfold_from_seed(const DelaunayTriangulation& D,
                                  const std::vector<CellPlacement>& cells,
                                  int seed_cell_id)
{
    LatticeUnfolding U;
    if (seed_cell_id < 0 || seed_cell_id >= D.nf || !cells[seed_cell_id].ok)
        return U;

    auto corner_pos = [&](const CellPlacement& F, int cone) -> Eisenstein {
        if (cone == F.c0) return F.P0;
        if (cone == F.c1) return F.P1;
        if (cone == F.c2) return F.P2;
        unfold_throw("unfold_iDT: cone %d not a corner of cell %d",
                     cone, F.cell_id);
    };

    // Per-cone FIRST position (the position recorded by the first cell
    // to place this cone in the DFS).  Used as the twin-arc anchor when
    // expanding into a new cell.
    std::map<int, Eisenstein> first_pos;
    auto record_cone = [&](int cone, Eisenstein p) {
        auto& all = U.cone_positions[cone];
        for (const auto& q : all) if (q == p) return;
        all.push_back(p);
    };

    auto place = [&](int cell_id, Eisenstein P0g, Eisenstein P1g,
                      Eisenstein P2g, int parent_id,
                      Eisenstein anchor_local, Eisenstein anchor_global,
                      const D6Affine& T)
    {
        LatticeUnfolding::UnfoldedCell g;
        const CellPlacement& F = cells[cell_id];
        g.cell_id        = cell_id;
        g.c0 = F.c0; g.c1 = F.c1; g.c2 = F.c2;
        g.P0 = P0g;  g.P1 = P1g;  g.P2 = P2g;
        g.parent_cell_id = parent_id;
        for (const auto& [p, v] : F.lattice_points) {
            Eisenstein q = p - anchor_local;
            if (T.reflect) q = complex_conj(q);
            q = q * T.unit + anchor_global;
            g.entries.push_back({q, v});
        }
        U.cells.push_back(std::move(g));
        ++U.n_cells;

        for (auto [cone, p] : std::initializer_list<std::pair<int, Eisenstein>>{
                {F.c0, P0g}, {F.c1, P1g}, {F.c2, P2g}})
        {
            auto it = first_pos.find(cone);
            if (it == first_pos.end()) first_pos[cone] = p;
            else if (it->second != p) ++U.n_tears;
            record_cone(cone, p);
        }
    };

    std::vector<bool> cell_placed(D.nf, false);

    // Seed: place at its local lattice positions in the global frame.
    {
        const CellPlacement& F0 = cells[seed_cell_id];
        D6Affine T_id{ .unit = Eisenstein(1, 0), .reflect = false };
        place(seed_cell_id, F0.P0, F0.P1, F0.P2,
              /*parent=*/-1, F0.P0, F0.P0, T_id);
        cell_placed[seed_cell_id] = true;
    }

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
        if (!cells[fid].ok) continue;
        const CellPlacement& F = cells[fid];

        int h_next = D.he_next[h];
        int u = D.he_origin[h];
        int v = D.he_origin[h_next];
        auto itu = first_pos.find(u);
        auto itv = first_pos.find(v);
        if (itu == first_pos.end() || itv == first_pos.end()) continue;
        Eisenstein Pu_g = itu->second, Pv_g = itv->second;
        Eisenstein Pu_l = corner_pos(F, u);
        Eisenstein Pv_l = corner_pos(F, v);
        Eisenstein z_local  = Pv_l - Pu_l;
        Eisenstein z_global = Pv_g - Pu_g;
        if (z_local.norm2() != z_global.norm2()) continue;
        D6Affine T = align(z_local, z_global);
        auto apply = [&](Eisenstein p) {
            Eisenstein q = p - Pu_l;
            if (T.reflect) q = complex_conj(q);
            return q * T.unit + Pu_g;
        };
        Eisenstein P0g = apply(F.P0), P1g = apply(F.P1), P2g = apply(F.P2);

        int parent_id = D.he_face[h ^ 1];
        place(fid, P0g, P1g, P2g, parent_id, Pu_l, Pu_g, T);
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
                            const std::vector<CellPlacement>& cells,
                            int start_cell_id,
                            int dfs_max_seeds)
{
    // Collect candidate seed cells.
    std::vector<int> seeds;
    if (start_cell_id >= 0 && start_cell_id < D.nf && cells[start_cell_id].ok) {
        seeds.push_back(start_cell_id);
    }
    if (dfs_max_seeds > (int)seeds.size() || seeds.empty()) {
        for (int f = 0; f < D.nf; ++f) {
            if (!cells[f].ok) continue;
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

    LatticeUnfolding best = unfold_from_seed(D, cells, seeds[0]);
    for (size_t i = 1; i < seeds.size(); ++i) {
        LatticeUnfolding u = unfold_from_seed(D, cells, seeds[i]);
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
            bool cone = (v < 12);
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

    // Duplicate hex vertices.
    dbg << "\n-- DUPLICATE hex vertices (at >=2 distinct global positions) --\n";
    int n_dup = 0;
    for (const auto& [v, occ] : occurrences) {
        if (v < 12) continue;
        std::set<std::pair<int, int>> distinct;
        for (const auto& [p, f] : occ) distinct.insert({p.first, p.second});
        if (distinct.size() > 1) {
            ++n_dup;
            dbg << "    hex " << v << " at " << distinct.size()
                << " positions:\n";
            for (const auto& [p, f] : occ) {
                dbg << "       face " << f << " pos (" << p.first
                    << ", " << p.second << ")\n";
            }
        }
    }
    if (n_dup == 0) dbg << "    (none)\n";

    // Missing hex vertices.
    dbg << "\n-- MISSING hex vertices (T_sorted vertex never placed) --\n";
    int n_miss = 0;
    for (int v = 12; v < T_sorted.N; ++v) {
        if (occurrences.find(v) == occurrences.end()) {
            ++n_miss;
            dbg << "    hex " << v << "\n";
        }
    }
    if (n_miss == 0) dbg << "    (none)\n";

    dbg << "\nsummary: " << n_dup << " duplicates, " << n_miss
        << " missing; unfolding " << U.n_cells << "/20 faces\n";
}
