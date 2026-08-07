#include "fullerenes/eisenstein_raster.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>

// =====================================================================
// walk_line
// =====================================================================

namespace {

// Integer gcd (positive result, handling both zeros as 0).
int igcd(int a, int b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

// Primitive walk_line: gcd(endpoint.a, endpoint.b) == 1.  The line
// passes through no intermediate lattice points.
WalkResult walk_line_primitive(const TriangulationView& T,
                               int u, int v, int w,
                               int dir_uv,
                               Eisenstein endpoint)
{
    WalkResult r;

    if (endpoint == Eisenstein(0, 0)) {
        r.walk.push_back({Eisenstein(0, 0), u});
        r.final_vertex = u;
        r.final_face_verts = {u, v, w};
        r.final_face_pos   = { Eisenstein(0, 0),
                               unit_direction(dir_uv),
                               unit_direction(dir_uv + 1) };
        return r;
    }

    int cu = u, cv = v, cw = w;
    Eisenstein PU(0, 0);
    Eisenstein PV = unit_direction(dir_uv);
    Eisenstein PW = unit_direction(dir_uv + 1);

    r.walk.push_back({PU, cu});
    r.face_path.push_back({cu, cv, cw});
    r.face_pos.push_back({PU, PV, PW});

    auto finalize = [&](int fv, Eisenstein FPU, Eisenstein FPV, Eisenstein FPW,
                         int fu, int fv_id, int fw_id) {
        r.final_vertex = fv;
        r.final_face_verts = {fu, fv_id, fw_id};
        r.final_face_pos   = {FPU, FPV, FPW};
    };

    if (PV == endpoint) {
        r.walk.push_back({PV, cv});
        finalize(cv, PU, PV, PW, cu, cv, cw);
        return r;
    }
    if (PW == endpoint) {
        r.walk.push_back({PW, cw});
        finalize(cw, PU, PV, PW, cu, cv, cw);
        return r;
    }

    long cur_t_num = 0, cur_t_den = 1;
    int  entry_edge = -1;           // 0=(U,V), 1=(V,W), 2=(W,U)

    for (int step = 0; step < walk_max_steps; ++step) {
        struct Cand {
            int idx;
            int A, B;
            Eisenstein PA, PB;
            long t_num, den;
        };
        Cand best{};
        bool has_best = false;

        auto try_edge = [&](int idx, int A, int B,
                             Eisenstein PA, Eisenstein PB)
        {
            if (idx == entry_edge) return;
            Eisenstein edge_vec = PB - PA;
            long den_raw = wedge(endpoint, edge_vec);
            if (den_raw == 0) return;
            long t_num_raw = wedge(PA, edge_vec);
            long s_num_raw = -wedge(endpoint, PA);
            long t_num = t_num_raw, s_num = s_num_raw, den = den_raw;
            if (den < 0) { t_num = -t_num; s_num = -s_num; den = -den; }
            // gcd=1 primitive: strict interior crossing.
            if (s_num <= 0 || s_num >= den) return;
            if (t_num * cur_t_den <= cur_t_num * den) return;
            if (!has_best || t_num * best.den < best.t_num * den) {
                best = {idx, A, B, PA, PB, t_num, den};
                has_best = true;
            }
        };

        try_edge(0, cu, cv, PU, PV);
        try_edge(1, cv, cw, PV, PW);
        try_edge(2, cw, cu, PW, PU);

        if (!has_best) {
            r.final_vertex = -1;
            return r;
        }

        int new_third = T.next_on_face(best.B, best.A);
        if (new_third < 0) { r.final_vertex = -1; return r; }
        Eisenstein P_new = best.PB + (best.PA - best.PB).nextCCW();

        cu = best.B; cv = best.A; cw = new_third;
        PU = best.PB; PV = best.PA; PW = P_new;
        entry_edge = 0;
        cur_t_num = best.t_num;
        cur_t_den = best.den;

        r.face_path.push_back({cu, cv, cw});
        r.face_pos.push_back({PU, PV, PW});

        if (PU == endpoint) { r.walk.push_back({PU, cu}); finalize(cu, PU, PV, PW, cu, cv, cw); return r; }
        if (PV == endpoint) { r.walk.push_back({PV, cv}); finalize(cv, PU, PV, PW, cu, cv, cw); return r; }
        if (PW == endpoint) { r.walk.push_back({PW, cw}); finalize(cw, PU, PV, PW, cu, cv, cw); return r; }
    }

    r.final_vertex = -1;
    return r;
}

// At hex vertex V (final vertex of a primitive sub-walk) at lattice
// position P_V, find the next face to enter so the line continues in
// direction `line_dir`.  Returns (next_v, next_w, dir_uv) for a
// primitive sub-walk starting at V.
struct NextFace {
    int  next_v;
    int  next_w;
    int  dir_uv;
    bool ok;
};

NextFace find_next_face_at_hex(const TriangulationView& T,
                               int V,
                               Eisenstein P_V,
                               const std::array<int, 3>& prev_face_verts,
                               const std::array<Eisenstein, 3>& prev_face_pos,
                               Eisenstein line_dir)
{
    NextFace out{-1, -1, -1, false};

    // Identify V's index in prev_face_verts.
    int V_idx = -1;
    for (int i = 0; i < 3; ++i)
        if (prev_face_verts[i] == V) { V_idx = i; break; }
    if (V_idx < 0) return out;

    // Pick any other corner as the "known neighbour" at V.
    int X_idx = (V_idx + 1) % 3;
    int X     = prev_face_verts[X_idx];
    Eisenstein P_X = prev_face_pos[X_idx];
    Eisenstein dXV  = P_X - P_V;    // X relative to V (must be a unit direction)

    // Global lattice direction of X from V.
    int d_X = -1;
    for (int d = 0; d < 6; ++d) {
        if (unit_direction(d) == dXV) { d_X = d; break; }
    }
    if (d_X < 0) return out;

    // X's index in T[V].
    const auto& nbrs = T[V];
    int deg = (int)nbrs.size();
    int j_X = -1;
    for (int j = 0; j < deg; ++j) if (nbrs[j] == X) { j_X = j; break; }
    if (j_X < 0) return out;

    // Frame offset: T[V][0] is at lattice direction (d_X - j_X) mod 6.
    int k0 = ((d_X - j_X) % 6 + 6) % 6;

    // Lattice direction of the line at V: the sextant [unit(j), unit(j+1))
    // containing line_dir, decided by exact integer sector membership
    // (Sector::contains_co is the [R, L) half-open test, i.e. floor-sextant
    // semantics: a line along a unit direction belongs to the sector it
    // OPENS).  line_dir is an Eisenstein pointing in the line's direction
    // (any positive-scale rep).  This replaces floor(atan2/(pi/3)), whose
    // result was 1-ULP-decided exactly on those boundary directions.
    int j_line = -1;
    for (int j = 0; j < 6; ++j)
        if (Sector{unit_direction(j), unit_direction(j + 1)}.contains_co(line_dir)) {
            j_line = j;
            break;
        }
    if (j_line < 0) return out;    // unreachable: the six sectors cover the circle

    // Face index i such that T[V][i] is at direction j_line.
    int i = ((j_line - k0) % 6 + 6) % 6;
    if (deg == 5 && i >= 5) return out;    // cone gap; shouldn't happen for hex
    if (i >= deg) return out;

    out.next_v = nbrs[i];
    out.next_w = nbrs[(i + 1) % deg];
    out.dir_uv = j_line;
    out.ok     = true;
    return out;
}

}  // anonymous namespace

WalkResult walk_line(const TriangulationView& T, int u, int v, int w,
                     int dir_uv, Eisenstein endpoint)
{
    int g = igcd(endpoint.first, endpoint.second);
    if (g <= 1) {
        return walk_line_primitive(T, u, v, w, dir_uv, endpoint);
    }

    Eisenstein sub(endpoint.first / g, endpoint.second / g);
    WalkResult total;

    int cu = u, cv = v, cw = w;
    int cur_dir = dir_uv;
    Eisenstein cur_base(0, 0);

    for (int step = 0; step < g; ++step) {
        WalkResult r = walk_line_primitive(T, cu, cv, cw, cur_dir, sub);
        if (r.final_vertex == -1) {
            total.final_vertex = -1;
            return total;
        }

        // Append walk entries (skip first if not first iteration; it
        // duplicates the previous sub-walk's endpoint).
        size_t start_i = (step == 0) ? 0 : 1;
        for (size_t i = start_i; i < r.walk.size(); ++i) {
            total.walk.push_back({ r.walk[i].first + cur_base,
                                    r.walk[i].second });
        }
        for (auto& f : r.face_path) total.face_path.push_back(f);
        for (auto& fp : r.face_pos) {
            std::array<Eisenstein, 3> shifted = {
                fp[0] + cur_base, fp[1] + cur_base, fp[2] + cur_base
            };
            total.face_pos.push_back(shifted);
        }
        total.final_vertex      = r.final_vertex;
        total.final_face_verts  = r.final_face_verts;
        for (int k = 0; k < 3; ++k)
            total.final_face_pos[k] = r.final_face_pos[k] + cur_base;

        cur_base = cur_base + sub;

        if (step == g - 1) break;

        // Transition to next sub-walk at the intermediate vertex.
        NextFace nf = find_next_face_at_hex(T,
                                            r.final_vertex,
                                            sub,   // P_V in r's local frame
                                            r.final_face_verts,
                                            r.final_face_pos,
                                            sub);
        if (!nf.ok) {
            total.final_vertex = -1;
            return total;
        }

        cu      = r.final_vertex;
        cv      = nf.next_v;
        cw      = nf.next_w;
        cur_dir = nf.dir_uv;
    }

    return total;
}

// =====================================================================
// scan_triangle
// =====================================================================

ScanLines scan_triangle(Eisenstein P0, Eisenstein P1, Eisenstein P2)
{
    // Owner wrapper over the header-inline span body (scan_triangle_into,
    // eisenstein_raster.hh -- the device batch pipelines run the same
    // body): owning storage in, the body's refusal converted to the
    // documented throw at this boundary.
    ScanLines out;
    const int b_min = std::min({P0.second, P1.second, P2.second});
    const int b_max = std::max({P0.second, P1.second, P2.second});
    out.lines.assign((size_t)(b_max - b_min + 1), {});
    const long n = scan_triangle_into(P0, P1, P2, out.lines, out.b_min, out.b_max);
    if (n < 0) {
        const long w = wedge(P1 - P0, P2 - P0);
        char buf[256];
        std::snprintf(buf, sizeof buf,
            "scan_triangle: triangle ((%d,%d), (%d,%d), (%d,%d)) "
            "is %s -- caller invariant violation",
            P0.first, P0.second, P1.first, P1.second, P2.first, P2.second,
            w == 0 ? "degenerate (collinear)" : "CW");
        throw std::runtime_error(buf);
    }
    return out;
}
