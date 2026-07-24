// Eisenstein-paint, Layer II: from a cone-point Alexandrov embedding to a
// full fullerene embedding (see the header banner).  This is the only paint
// translation unit that touches 3D coordinates; the intrinsic chart
// machinery lives in eisenstein_paint.cc.

#include "fullerenes/eisenstein_paint_geometry.hh"

#include "fullerenes/auxiliary.hh"           // IDCounter
#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/barycentric.hh"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace eisenstein_paint {

// Internal Layer-I export (eisenstein_paint.cc): the raw iDT compute
// without the chartability guard -- the Alexandrov flips supersede it.
DelaunayTriangulation detail_compute_dual_idt(const SortedDual& S, const char* who);

namespace {

[[noreturn]] void paint_throw(const char* fmt, ...) {
    char buf[512];
    va_list ap; va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    throw std::runtime_error(buf);
}

}  // namespace

// =====================================================================
// realize_dual: Alexandrov realization of the dual cone metric.
// =====================================================================

DualPolytope realize_dual(const SortedDual& S) {
    DualPolytope P;
    P.D = detail_compute_dual_idt(S, "eisenstein_paint::realize_dual");

    AlexandrovSolver A;
    A.D        = std::move(P.D);
    P.cone_pos = A.solve();
    P.D        = std::move(A.D);
    if (A.valid()) return P;

    using VS = AlexandrovSolver::ValidationStatus;
    const std::string detail = std::string("AlexandrovSolver: ")
                             + AlexandrovSolver::status_str(A.stats_status);
    const Stage st = (A.stats_status == VS::FAIL_NOT_SIMPLE)
                       ? Stage::NON_SIMPLICIAL
                       : Stage::ALEXANDROV;
    throw StageError(st, "eisenstein_paint::realize_dual: " + detail);
}

Deltahedron DualPolytope::deltahedron() const {
    // CCW neighbour rings straight from the DCEL: incident(v) is the CW
    // ring, so reverse it.  Orientation is constructed, never inferred.
    std::vector<std::vector<node_t>> nb(D.nv);
    for (int v = 0; v < D.nv; ++v) {
        if (D.v_out[v] < 0)
            throw std::logic_error("DualPolytope::deltahedron: dead vertex "
                                   + std::to_string(v) + " in T(0)");
        for (int h : D.incident(v)) nb[v].push_back(D.he_origin[h ^ 1]);
        std::reverse(nb[v].begin(), nb[v].end());
        std::vector<node_t> sorted_nb(nb[v].begin(), nb[v].end());
        std::sort(sorted_nb.begin(), sorted_nb.end());
        if (std::adjacent_find(sorted_nb.begin(), sorted_nb.end()) != sorted_nb.end())
            throw std::logic_error(
                "DualPolytope::deltahedron: T(0) is non-simplicial at vertex "
                + std::to_string(v) + " (flat-diagonal multi-edge); the "
                "1-skeleton is not a simple graph -- use "
                "AlexandrovSolver::polytope_tesselation for T-bar(0) instead");
    }
    if ((int)cone_pos.size() != D.nv)
        throw std::logic_error("DualPolytope::deltahedron: cone_pos size "
                               + std::to_string(cone_pos.size()) + " != nv "
                               + std::to_string(D.nv));
    Triangulation T{Graph(Spanify::OwnedDenseGraph<node_t>(nb))};
    return Deltahedron(T, cone_pos);
}

// =====================================================================
// Chart evaluation.
// =====================================================================

void interpolate_cell(const Cell& F,
                      const LatticeMap& lmap,
                      const std::vector<coord3d>& anchors,
                      int n_cones,
                      std::vector<coord3d>& pos3d)
{
    if (!F.ok)
        paint_throw("eisenstein_paint::interpolate_cell: F.ok=false");
    if ((int)anchors.size() != n_cones)
        paint_throw("eisenstein_paint::interpolate_cell: anchors size %zu != n_cones %d",
                    anchors.size(), n_cones);
    const long g = wedge(F.P1 - F.P0, F.P2 - F.P0);
    if (g <= 0)
        paint_throw("eisenstein_paint::interpolate_cell(cell %d): g=%ld <= 0",
                    F.cell_id, g);
    const coord3d& C0 = anchors[F.c0];
    const coord3d& C1 = anchors[F.c1];
    const coord3d& C2 = anchors[F.c2];

    for (const auto& [p, vid] : lmap.entries) {
        if (vid < n_cones) continue;       // cone -- pre-painted
        const IntBary bw = integer_barycentric(p, F.P0, F.P1, F.P2);
        if (bw.n0 < 0 || bw.n1 < 0 || bw.n2 < 0 || bw.denom != g)
            paint_throw("eisenstein_paint::interpolate_cell(cell %d): bad "
                        "barycentric for vertex %d at pos=(%d,%d): "
                        "(n0=%ld n1=%ld n2=%ld, denom=%ld, g=%ld)",
                        F.cell_id, vid, p.first, p.second,
                        bw.n0, bw.n1, bw.n2, bw.denom, g);
        pos3d[vid] = barycentric_combine(reduce_to_lowest_terms(bw), C0, C1, C2);
    }
}

std::vector<coord3d> evaluate(const SurfaceParametrization& A,
                              const std::vector<coord3d>& anchors,
                              const Permutation& perm)
{
    // Cones are written directly; interpolate_cell skips them
    // ("pre-painted" contract).  Non-cone slots are written by
    // interpolation, on-edge vertices idempotently from both adjacent
    // cells.
    const int Nv_sorted = (int)A.T.N;
    std::vector<coord3d> pos3d(Nv_sorted);
    for (int c = 0; c < A.n_cones; ++c) pos3d[c] = anchors[c];
    try {
        for (size_t fi = 0; fi < A.cells.size(); ++fi)
            if (A.cells[fi].ok)
                interpolate_cell(A.cells[fi], A.lmaps[fi], anchors,
                                 A.n_cones, pos3d);
    } catch (const std::exception& e) {
        throw StageError(Stage::INTERPOLATE,
            std::string("interpolate_throw: ") + e.what());
    }

    // perm[u_orig] = u_sorted, so the position written at sorted vertex
    // perm[u_orig] is u_orig's 3D coordinate.
    const int Nv_orig = (int)perm.size();
    std::vector<coord3d> out(Nv_orig);
    for (int u = 0; u < Nv_orig; ++u)
        out[u] = pos3d[perm[u]];
    return out;
}

// =====================================================================
// The cubic polytope + assembly.
// =====================================================================

CubicPolytope realize_cubic(const TriangulationView& T) {
    // AlexandrovIDTCubic solve + cone bookkeeping.  Staging mirrors
    // realize_dual: kis/iDT trouble -> IDT_COMPUTE, non-simplicial cubic
    // iDT -> NON_SIMPLICIAL, any other non-validating outcome -> ALEXANDROV.
    CubicPolytope C;
    AlexandrovIDTCubic AC;
    try {
        C.cone_pos = AC.solve(T);
    } catch (const std::exception& e) {
        throw StageError(Stage::IDT_COMPUTE,
            std::string("eisenstein_paint::realize_cubic: AlexandrovIDTCubic: ")
            + e.what());
    }
    C.cone_triangle = AC.cone_triangle;

    if (!AC.solver.valid()) {
        using VS = AlexandrovSolver::ValidationStatus;
        const std::string detail = std::string("AlexandrovIDTCubic: ")
                                 + AlexandrovSolver::status_str(AC.solver.stats_status);
        const Stage st = (AC.solver.stats_status == VS::FAIL_NOT_SIMPLE)
                           ? Stage::NON_SIMPLICIAL
                           : Stage::ALEXANDROV;
        throw StageError(st, "eisenstein_paint::realize_cubic: " + detail);
    }
    if (C.cone_pos.size() != C.cone_triangle.size())
        throw StageError(Stage::ALEXANDROV,
            "eisenstein_paint::realize_cubic: AlexandrovIDTCubic returned "
            + std::to_string(C.cone_pos.size()) + " positions for "
            + std::to_string(C.cone_triangle.size()) + " cones");
    return C;
}

std::vector<coord3d> pentagon_anchors(const TriangulationView& T,
                                      const Permutation& perm,
                                      const CubicPolytope& C)
{
    // Every dual triangle containing a deg-5 vertex is pentagon-incident
    // and hence a cone of C, so each pentagon accumulates exactly 5
    // contributions; anything else is broken cone bookkeeping.
    std::vector<coord3d> anchors(12, coord3d{0, 0, 0});
    std::vector<int>     count(12, 0);
    for (size_t i = 0; i < C.cone_triangle.size(); ++i)
        for (int k = 0; k < 3; ++k) {
            const int p = C.cone_triangle[i][k];
            if (T.degree(p) == 6) continue;      // hexagon corner
            const int ps = perm[p];              // sorted label, < 12
            if (ps < 0 || ps >= 12)
                paint_throw("eisenstein_paint::pentagon_anchors: deg-%d vertex %d "
                            "has sorted label %d (expected cone label < 12)",
                            T.degree(p), p, ps);
            anchors[ps] += C.cone_pos[i];
            ++count[ps];
        }
    for (int c = 0; c < 12; ++c) {
        if (count[c] != 5)
            paint_throw("eisenstein_paint::pentagon_anchors: cone %d accumulated "
                        "%d cone-vertex contributions (expected 5)", c, count[c]);
        anchors[c] = anchors[c] / 5.0;
    }
    return anchors;
}

std::vector<coord3d> assemble_cubic(const TriangulationView& T,
                                    const std::vector<coord3d>& dual_coords,
                                    const CubicPolytope& C)
{
    const std::vector<tri_t> tris = T.triangles();
    IDCounter<tri_t> tri_numbers;
    for (const tri_t& t : tris) tri_numbers.insert(t.sorted());

    std::vector<coord3d> cubic(tris.size());
    for (size_t U = 0; U < tris.size(); ++U)
        cubic[U] = tris[U].centroid(dual_coords);

    for (size_t i = 0; i < C.cone_triangle.size(); ++i) {
        const auto it = tri_numbers.find(C.cone_triangle[i].sorted());
        if (it == tri_numbers.end())
            paint_throw("eisenstein_paint::assemble_cubic: cone %zu's dual "
                        "triangle (%d,%d,%d) not in T.triangles()",
                        i, C.cone_triangle[i][0], C.cone_triangle[i][1],
                        C.cone_triangle[i][2]);
        cubic[it->second] = C.cone_pos[i];
    }
    return cubic;
}

// =====================================================================
// Facades.
// =====================================================================

namespace {

// Fill R.status from the caught failure; shared by both facades.
template <typename Body>
Status run_facade(Body&& body) {
    Status st;
    try {
        body();
    } catch (const StageError& e) {
        st.stage = e.stage;
        st.why   = e.what();
    } catch (const std::exception& e) {
        st.stage = Stage::UNEXPECTED;
        st.why   = std::string("unexpected: ") + e.what();
    }
    return st;
}

// Cubic-graph bond statistics of `cubic` (indexed in T.triangles()
// order): two cubic vertices are bonded iff their dual triangles share
// an arc.  Derived from the arc -> face-left table (each directed arc
// lies in exactly one CCW triangle).
void fill_bond_stats(const TriangulationView& T,
                     const std::vector<coord3d>& cubic,
                     CubicGeometry& R)
{
    const std::vector<tri_t> tris = T.triangles();
    std::map<arc_t, int> face_of_arc;
    for (int t = 0; t < (int)tris.size(); ++t)
        for (int j = 0; j < 3; ++j)
            face_of_arc[{tris[t][j], tris[t][(j + 1) % 3]}] = t;

    double sum = 0, sum2 = 0;
    int n = 0;
    R.bond_min = 0; R.bond_max = 0;
    for (const auto& [arc, U] : face_of_arc) {
        const auto it = face_of_arc.find({arc.second, arc.first});
        if (it == face_of_arc.end() || it->second <= U) continue;   // count each pair once
        const double L = (cubic[U] - cubic[it->second]).norm();
        if (n == 0) { R.bond_min = R.bond_max = L; }
        else {
            R.bond_min = std::min(R.bond_min, L);
            R.bond_max = std::max(R.bond_max, L);
        }
        sum += L; sum2 += L * L;
        ++n;
    }
    if (n == 0) return;
    R.bond_mean = sum / n;
    const double var = sum2 / n - R.bond_mean * R.bond_mean;
    R.bond_cv = R.bond_mean > 0 ? std::sqrt(std::max(0.0, var)) / R.bond_mean : 0;
}

}  // namespace

DualGeometry dual_geometry(const TriangulationView& T) {
    DualGeometry R;
    R.status = run_facade([&] {
        const SortedDual   S = sorted_dual(T);
        const DualPolytope P = realize_dual(S);
        const SurfaceParametrization A = parametrize(P.D, S);
        R.coords = evaluate(A, P.cone_pos, S.perm);
    });
    if (!R.status.ok()) R.coords.clear();
    return R;
}

CubicGeometry cubic_geometry(const TriangulationView& T) {
    CubicGeometry R;
    R.status = run_facade([&] {
        const SortedDual    S = sorted_dual(T);
        // BOTH metrics are realized: the cubic polytope supplies the
        // anchors, the dual polytope supplies the FLAT cells the charts
        // are evaluated on.  Evaluating on the raw iDT instead is what
        // mangled elongated isomers (see the header banner); the dual
        // solve (12 cones) is cheap next to the cubic one (20..60).
        const DualPolytope  P = realize_dual(S);
        const CubicPolytope C = realize_cubic(T);
        const std::vector<coord3d> anchors = pentagon_anchors(T, S.perm, C);
        const SurfaceParametrization A = parametrize(P.D, S);
        R.dual_coords  = evaluate(A, anchors, S.perm);
        R.cubic_coords = assemble_cubic(T, R.dual_coords, C);
        fill_bond_stats(T, R.cubic_coords, R);
    });
    if (!R.status.ok()) { R.cubic_coords.clear(); R.dual_coords.clear(); }
    return R;
}

}  // namespace eisenstein_paint
