#include "fullerenes/delaunay_strip.hh"
#include "fullerenes/eisenstein_raster.hh"
#include "fullerenes/eisenstein_tikz.hh"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ostream>
#include <set>
#include <unordered_map>
#include <utility>

// Lib-side: global namespace.

// Diagnostic flag (toggle via set_unfold_arc_strip_debug; CLI-driven from
// test drivers — no environment variables).
static bool& unfold_arc_strip_debug_flag() { static bool v = false; return v; }
void set_unfold_arc_strip_debug(bool on) { unfold_arc_strip_debug_flag() = on; }
bool unfold_arc_strip_debug() { return unfold_arc_strip_debug_flag(); }


// Index of nbr in T[v]'s neighbour list, or -1.
static int index_of_neighbour(const Triangulation& T, int v, int nbr) {
    const auto& nbrs = T[v];
    int deg = (int)nbrs.size();
    for (int i = 0; i < deg; ++i) if (nbrs[i] == nbr) return i;
    return -1;
}

using StripData = std::unordered_map<int, std::pair<Eisenstein, int>>;

// Register a vertex; first registration wins (don't overwrite k0).
static void register_vertex(StripData& data, int v, Eisenstein pos, int k0) {
    if (data.find(v) != data.end()) return;
    data[v] = {pos, k0};
}

// Compute k0_canon for vertex w given a known neighbour `nbr` of w whose
// position relative to w is along lattice direction `d` (i.e.,
// P_nbr - P_w == unit_direction(d)).
//   k0[w] = (d - index_of(nbr, T[w])) mod 6
// Returns -1 if nbr is not in T[w]'s neighbour list (shouldn't happen
// for valid walks).
static int derive_k0_from_neighbour(const Triangulation& T, int w,
                                     int nbr, int d) {
    int j = index_of_neighbour(T, w, nbr);
    if (j < 0) return -1;
    return ((d - j) % 6 + 6) % 6;
}

// Determine the lattice direction d such that unit_direction(d) == v.
// Returns -1 if v is not a unit lattice direction.
static int direction_of_unit(Eisenstein v) {
    for (int d = 0; d < 6; ++d)
        if (unit_direction(d) == v) return d;
    return -1;
}

// Register every face vertex in `wr.face_path` / `wr.face_pos` into
// `data`, deriving k0_canon for each from a face-neighbour pair.
//
// For each face (a, b, c) with positions (P_a, P_b, P_c):
//   - For corner a: a's neighbour b sits at direction d_ab where
//     unit_direction(d_ab) == P_b - P_a. Then k0[a] = (d_ab - j_b) mod 6,
//     where j_b = index_of(b, T[a]).
// (Symmetric for b, c.)
//
// Aborts on malformed input (face edges should always be unit lattice
// vectors after the walk's frame propagation).
static bool register_walk(const Triangulation& T,
                           const WalkResult& wr,
                           StripData& data)
{
    if (wr.face_path.size() != wr.face_pos.size()) return false;
    for (size_t fi = 0; fi < wr.face_path.size(); ++fi) {
        const auto& F = wr.face_path[fi];
        const auto& P = wr.face_pos[fi];
        // Three corner registrations using the next CCW corner as the
        // reference neighbour.
        for (int k = 0; k < 3; ++k) {
            int a = F[k];
            int b = F[(k + 1) % 3];
            Eisenstein dvec = P[(k + 1) % 3] - P[k];
            int d = direction_of_unit(dvec);
            if (d < 0) {
                std::fprintf(stderr,
                    "fillin::register_walk: face edge not unit-direction "
                    "(face %zu corner %d, vertex %d, dvec=(%d,%d))\n",
                    fi, k, a, dvec.first, dvec.second);
                return false;
            }
            int k0 = derive_k0_from_neighbour(T, a, b, d);
            if (k0 < 0) {
                std::fprintf(stderr,
                    "fillin::register_walk: %d not a neighbour of %d\n",
                    b, a);
                return false;
            }
            register_vertex(data, a, P[k], k0);
        }
    }
    return true;
}

// Internal termination guard: hard ceiling on per-walk steps. Picked
// generously above any plausible iDT edge length on realistic
// fullerenes (largest expected ~ a + b + 1 < 100).
static constexpr int STRIP_WALK_MAX_STEPS = 10000;

// Walk along a single lattice direction (0 = east, 1 = north) for `n`
// steps from u, registering only the on-line vertices visited. Used for
// axis-aligned (b == 0 or a == 0) iDT edges where the geodesic runs
// along T_sorted edges; the "strip" in this case is degenerate (no
// vertex strictly above or below the line) and contains only on-line
// vertices.
//
// Aborts with diagnostic if n exceeds STRIP_WALK_MAX_STEPS.
// Returns the end vertex (last visited) on success, -1 on failure.
static int walk_axis_aligned(const Triangulation& T,
                              int u, int k0_u, int dir, int n,
                              StripData& data)
{
    if (n > STRIP_WALK_MAX_STEPS) {
        char buf[200];
        std::snprintf(buf, sizeof buf,
            "fillin::walk_axis_aligned: n=%d exceeds MAX_STEPS=%d "
            "(u=%d, dir=%d) -- likely malformed input",
            n, STRIP_WALK_MAX_STEPS, u, dir);
        throw std::runtime_error(buf);
    }
    int curr = u, curr_k0 = k0_u;
    Eisenstein curr_pos(0, 0);
    register_vertex(data, curr, curr_pos, curr_k0);
    for (int step = 0; step < n; ++step) {
        int deg = (int)T[curr].size();
        int nbr_idx = ((dir - curr_k0) % 6 + 6) % 6;
        if (deg == 5 && nbr_idx == 5) return -1;       // cone gap
        if (nbr_idx >= deg) return -1;
        int next_v = T[curr][nbr_idx];
        int j = index_of_neighbour(T, next_v, curr);
        if (j < 0) return -1;
        int next_k0 = ((dir + 3 - j) % 6 + 6) % 6;
        curr_pos = curr_pos + unit_direction(dir);
        register_vertex(data, next_v, curr_pos, next_k0);
        curr = next_v; curr_k0 = next_k0;
    }
    return curr;
}

// Try one alignment (i = neighbour index at u to send east in the
// canonical frame, z_try = candidate target lattice position for v).
//
// Builds the initial face (u, v_a=T[u][i], v_b=T[u][(i+1) mod deg]) at
// canonical positions ((0,0), unit(0), unit(1)), runs walk_line to
// z_try, and on success extracts strip data from the full face_path.
//
// Rejects walks that pass through a cone other than u and v: at a
// cone the angle deficit makes "continuing straight" undefined, so a
// chain that touches an intermediate cone is not a valid geodesic
// (it can reach v by accident on cone-clustered fullerenes).
static bool try_walk_via_walk_line(const Triangulation& T,
                                    int u, int v, int i,
                                    Eisenstein z_try,
                                    StripData& out_data)
{
    int deg = (int)T[u].size();
    if (i < 0 || i >= deg) return false;
    int v_a = T[u][i];
    int v_b = T[u][(i + 1) % deg];

    // Initial face must be CCW in T_sorted — by Triangulation invariants
    // T[u][i] and T[u][(i+1) mod deg] form a face (u, T[u][i], T[u][(i+1)
    // mod deg]) in CCW orientation, so this is automatic.
    WalkResult wr = walk_line(T, u, v_a, v_b, /*dir_uv=*/0, z_try);
    if (wr.final_vertex != v) return false;

    StripData local;
    if (!register_walk(T, wr, local)) return false;

    // Sanity: u must be at (0,0) and v at z_try.
    auto it_u = local.find(u);
    auto it_v = local.find(v);
    if (it_u == local.end() || it_u->second.first != Eisenstein(0, 0))
        return false;
    if (it_v == local.end() || it_v->second.first != z_try)
        return false;

    // Reject chains whose LINE segment passes directly through a cone
    // (= an intermediate on-line lattice point that's a cone).  At a
    // cone the angle deficit makes "continuing straight" undefined, so
    // a geodesic cannot pass through one.  unfold_arc_strip does not know
    // F, so the F-interior cone-freedom check in run_frame_walker isn't
    // available here — the per-arc strip retains off-line cone
    // neighbours, and it's up to the face-level pipeline to filter
    // them when assembling F's lmap.
    for (size_t k = 1; k + 1 < wr.walk.size(); ++k) {
        int vid = wr.walk[k].second;
        if (vid < 12) return false;
    }

    out_data = std::move(local);
    return true;
}

Strip unfold_arc_strip(const DelaunayTriangulation& D,
                       const Triangulation& T_sorted,
                       int arc_h)
{
    Strip S;
    int h = arc_h;
    if (h < 0 || h >= D.nh) return S;
    if (!D.alive(h)) return S;
    if (unfold_arc_strip_debug()) {
        std::fprintf(stderr, "unfold_arc_strip arc_h=%d u=%d v=%d L=%g\n",
                     h, D.he_origin[h], D.he_origin[h ^ 1],
                     D.he_length[h]);
    }

    int u = D.he_origin[h];
    int v = D.he_origin[h ^ 1];

    double L = D.he_length[h];
    int N = (int)std::lround(L * L);
    if (N <= 0) return S;

    // Enumerate all sector-0 reps of N: 1 entry for most norms, 2 for
    // split-prime norms (covering both rotation orbits).
    std::vector<Eisenstein> z_candidates = sector0_reps_of_norm(N);
    if (z_candidates.empty()) {
        char buf[200];
        std::snprintf(buf, sizeof buf,
            "fillin::unfold_arc_strip: no Eisenstein integer at norm %d "
            "for edge {%d, %d}", N, u, v);
        throw std::runtime_error(buf);
    }

    int deg_u = (int)T_sorted[u].size();

    // Split-prime norms admit TWO sector-0 reps in distinct rotation
    // orbits.  For symmetric surfaces the walker may land v at a
    // target in EITHER orbit (surface has multiple same-length paths
    // from u to v crossing different T_sorted face chains).  Collect
    // every successful (z_chosen, strip data) pair — the face placement
    // will select the one consistent with F's closure.
    struct Candidate { Eisenstein z_chosen; StripData data; };
    std::vector<Candidate> successes;

    for (size_t ci = 0; ci < z_candidates.size(); ++ci) {
        Eisenstein z_try = z_candidates[ci];
        int a = z_try.first, b = z_try.second;

        // Axis-aligned (b == 0 or a == 0): the line runs along T_sorted
        // edges; the strip is degenerate (no off-line vertex strictly
        // above or below). Walk on-line vertices only.
        if (b == 0 || a == 0) {
            int dir = (b == 0) ? 0 : 1;
            int n   = (b == 0) ? a : b;
            for (int i = 0; i < deg_u; ++i) {
                int k0_u = ((6 - i) % 6 + 6) % 6;
                int nbr_idx = ((dir - k0_u) % 6 + 6) % 6;
                if (nbr_idx != i) continue;     // alignment inconsistent
                StripData local;
                int end_v = walk_axis_aligned(T_sorted, u, k0_u, dir, n, local);
                if (end_v == v) {
                    successes.push_back({z_try, std::move(local)});
                    break;      // one success per z_try is enough
                }
            }
            continue;
        }

        // General/tied: full geodesic walk via walk_line.
        for (int i = 0; i < deg_u; ++i) {
            StripData local;
            if (try_walk_via_walk_line(T_sorted, u, v, i, z_try, local)) {
                successes.push_back({z_try, std::move(local)});
                break;
            }
            if (unfold_arc_strip_debug()) {
                std::fprintf(stderr,
                    "  try_walk u=%d v=%d z_try=(%d,%d) i=%d -> failed\n",
                    u, v, z_try.first, z_try.second, i);
            }
        }
    }
    if (successes.empty()) {
        char buf[200];
        std::snprintf(buf, sizeof buf,
            "fillin::unfold_arc_strip: alignment failed for arc {%d -> %d} "
            "norm=%d (%zu candidate sector-0 reps tried)",
            u, v, N, z_candidates.size());
        throw std::runtime_error(buf);
    }

    // Return the first success as the "primary" strip; ancillary orbits
    // (if any) are queried via compute_strip_with_forced_z_e from
    // place_face when orbit disagreement needs to be resolved.
    const Candidate& chosen = successes[0];
    S.u   = u;
    S.v   = v;
    S.v_position = chosen.z_chosen;
    for (const auto& kv : chosen.data) {
        StripVertex sv;
        sv.vertex_id = kv.first;
        sv.position       = kv.second.first;
        sv.frame_offset  = kv.second.second;
        long w_sign = wedge(chosen.z_chosen, sv.position);
        if (w_sign > 0) {
            S.left.push_back(sv);
        } else if (w_sign < 0) {
            S.right.push_back(sv);
        } else {
            // on-line: in both
            S.left.push_back(sv);
            S.right.push_back(sv);
        }
    }
    S.ok = true;
    return S;
}

std::vector<Strip> unfold_all_arc_strips(const DelaunayTriangulation& D,
                                         const Triangulation& T_sorted)
{
    std::vector<Strip> strips(D.nh);
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int h = 0; h < D.nh; ++h) {
        strips[h] = unfold_arc_strip(D, T_sorted, h);
    }
    return strips;
}

void dump_strip_tikz(const Strip& S,
                     const Triangulation& T_sorted,
                     std::ostream& os,
                     double scale)
{
    using tikz::cart;

    // Side: 1 = above-only, -1 = below-only, 0 = on-line (in both).
    struct Info { Eisenstein pos; int k0; int side; };
    std::unordered_map<int, Info> verts;
    for (const auto& sv : S.left) {
        auto& info = verts[sv.vertex_id];
        info.pos = sv.position; info.k0 = sv.frame_offset; info.side |= 1;
    }
    for (const auto& sv : S.right) {
        auto& info = verts[sv.vertex_id];
        info.pos = sv.position; info.k0 = sv.frame_offset; info.side |= 2;
    }
    auto side_of = [](int s) { return s == 3 ? 0 : (s == 1 ? +1 : -1); };

    os << "\\documentclass[tikz,border=4pt]{standalone}\n"
          "\\usepackage{tikz}\n"
          "\\begin{document}\n"
          "\\begin{tikzpicture}["
       << "scale=" << scale
       << ", every node/.style={font=\\tiny}]\n";

    os << "% iDT edge u=" << S.u << " --> v=" << S.v
       << "  z_e=(" << S.v_position.first << "," << S.v_position.second
       << ")  norm^2=" << S.v_position.norm2()
       << "  strip size = " << verts.size() << "\n";

    // Eisenstein lattice background — extend well beyond the strip so
    // the grid reads as a backdrop, not a tight frame. Margin scales
    // with the strip size so it stays generous on short and long strips.
    tikz::BBox bb;
    for (const auto& kv : verts) bb.bump(kv.second.pos);
    const double margin = std::max(5.0, 0.5 * bb.span());
    bb.xmin -= margin; bb.xmax += margin;
    bb.ymin -= margin; bb.ymax += margin;
    os << "% Eisenstein lattice background\n";
    tikz::emit_grid(os, bb);

    // T_sorted edges between strip vertices.
    os << "% T_sorted edges within the strip\n";
    std::set<std::pair<int, int>> drawn;
    for (const auto& kv : verts) {
        int w = kv.first;
        Eisenstein pw = kv.second.pos;
        for (int n : T_sorted[w]) {
            auto it = verts.find(n);
            if (it == verts.end()) continue;
            std::pair<int, int> e = std::minmax(w, n);
            if (!drawn.insert(e).second) continue;
            auto [x1, y1] = cart(pw);
            auto [x2, y2] = cart(it->second.pos);
            os << "  \\draw[gray!55,line width=0.3pt] ("
               << x1 << "," << y1 << ") -- (" << x2 << "," << y2 << ");\n";
        }
    }

    // Directed iDT geodesic u -> v (drawn UNDER the vertices).
    auto [ux, uy] = cart(Eisenstein(0, 0));
    auto [vx, vy] = cart(S.v_position);
    os << "% iDT geodesic u -> v\n"
       << "  \\draw[->,red!80!black,line width=1pt] ("
       << ux << "," << uy << ") -- (" << vx << "," << vy << ");\n";

    // Strip vertices.
    os << "% strip vertices (square = endpoint cone, circle = interior)\n";
    for (const auto& kv : verts) {
        int v_id = kv.first;
        auto [x, y] = cart(kv.second.pos);
        int s = side_of(kv.second.side);
        const char* fill = (s ==  0) ? "red!75!black"
                          : (s == +1) ? "blue!70!black"
                                      : "green!55!black";
        bool endpoint = (v_id == S.u || v_id == S.v);
        const char* shape = endpoint ? "rectangle" : "circle";
        os << "  \\node[" << shape << ",draw=black,fill=" << fill
           << ",inner sep=1.6pt,label=below right:{\\tiny " << v_id
           << "}] at (" << x << "," << y << ") {};\n";
    }

    os << "\\end{tikzpicture}\n\\end{document}\n";
}

// end of file
