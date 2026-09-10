// eisenstein_atlas.cc -- Tier 3.5: the exact cross-cell atlas (see the header
// banner).  Promoted from claude-projects/curvature-flow/src/cage_transfer.cc
// (stage A of the tab-4 cage transfer), where it was validated by the
// exactness oracle: transferring a candidate onto its own equilateral complex
// reproduces eisenstein_paint::run's dual-vertex positions to ~1e-15.
//
// All arithmetic in this file is integer (long / __int128) or rational over
// integers -- no floating point enters any predicate or position.

#include "fullerenes/eisenstein_atlas.hh"
#include "fullerenes/union_find.hh"   // UnionFind (intrinsic_dual's point gluing)

#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace eisenstein_paint {
namespace {

[[noreturn]] void fail(const std::string& what) {
  throw std::logic_error("eisenstein_atlas: " + what);
}

// Checked narrowing (the header's width contract): denominators and sample
// weights must fit int; a violated bound throws, never silently truncates.
inline int narrow(long x, const char* what) {
  if (x < std::numeric_limits<int>::min() || x > std::numeric_limits<int>::max())
    fail(std::string(what) + " exceeds int width");
  return (int)x;
}

// q = num/den inside the CLOSED cell with frame corners P (all three edge
// half-planes)?  Cross-multiplied exact predicate; __int128 keeps
// den * coordinate exact.  den is narrow-checked ONCE by the caller.
inline bool inside_cell(const std::array<Eisenstein, 3>& P, Eisenstein num, int den) {
  for (int k = 0; k < 3; k++) {
    const Eisenstein Pi = P[k], Pj = P[(k + 1) % 3];
    const Eisenstein ev = Pj - Pi;
    const Eisenstein rel = num - Pi * den;
    if ((__int128)ev.first * rel.second - (__int128)ev.second * rel.first < 0) return false;
  }
  return true;
}

// T applied to the raw rational num/den WITHOUT re-normalizing: trace_segment
// keeps its two endpoints on one common denominator, which a canonicalizing
// constructor would destroy.  den is narrow-checked ONCE by the caller.
inline Eisenstein apply_raw(const LatticeIsometry& T, Eisenstein num, int den) {
  const Eisenstein n = T.reflect ? num.complex_conj() : num;
  return T.t * den + T.u * n;
}

// Corner index k of face t whose CCW arc (t[k], t[k+1]) carries the
// undirected edge {u, v}; -1 when the face does not carry that edge.
inline int corner_index_of_edge(const tri_t& t, int u, int v) {
  for (int k = 0; k < 3; k++) {
    const int x = t[k], y = t[(k + 1) % 3];
    if ((x == u && y == v) || (x == v && y == u)) return k;
  }
  return -1;
}

// Develop face t into a chart from the developed positions pu, pv of the
// undirected edge {eu, ev} it carries at corner k_arc (corner_index_of_edge):
// the edge's corners keep their positions (matched by vertex id), the third
// corner is the unit apex.  Returns the corner positions in face-slot order.
inline std::array<Eisenstein, 3> develop_face_on_edge(const tri_t& t, int k_arc, int eu,
                                                      Eisenstein pu, Eisenstein pv) {
  const Eisenstein px = (t[k_arc] == eu) ? pu : pv;
  const Eisenstein py = (t[k_arc] == eu) ? pv : pu;
  std::array<Eisenstein, 3> p{};
  p[k_arc]           = px;
  p[(k_arc + 1) % 3] = py;
  p[(k_arc + 2) % 3] = unit_apex(px, py);
  return p;
}

// ---------------------------------------------------------------------------
// build_atlas: a composition of value-producing constructions.  Each takes
// exactly the values it depends on, so the dependency order is the argument
// flow -- enforced by the compiler, not by call-order convention.
// ---------------------------------------------------------------------------

// The corner-k = k-th-cycle-origin theorem, gated: every live face has a
// live chart whose CCW corner ids equal the origins of its face cycle.
// cell_metric constructs charts exactly so; everything downstream --
// transition matching and trace_segment's exit selection -- stands on it,
// so the atlas verifies it per cell instead of trusting the cross-module
// invariant silently.
void check_charts_match_cycles(const DelaunayView& D, const ParamTablesView& V) {
  for (int f = 0; f < D.nf; f++) {
    if (D.f_he[f] < 0) continue;
    if (!V.cell_live(f)) fail("live cell " + std::to_string(f) + " not charted");
    const auto cyc = D.face_halfedges(f);
    const auto ids = V.corner_ids(f);
    for (int k = 0; k < 3; k++)
      if (D.he_origin[cyc[k]] != ids[k])
        fail("cell " + std::to_string(f) + " chart corner order departs from its face cycle");
  }
}

// THE transition body: the chart isometry across half-edge h, from the
// frame Pf of face(h) onto the frame Pg of face(twin h), matched on the
// crossed edge's corner SLOTS in each face cycle (Lemma: transitions are
// exact).  Slot matching -- not cone-id matching -- so a repeated-corner
// cell cannot confuse the pairing; the twin cycles the same edge in the
// opposite direction, hence the reversed slot pair on the g side.
//   @throws std::logic_error (from isometry_from_segments) when the two
//           frames do not develop the shared edge into one another --
//           the falsifier the intrinsic selection tests candidates by.
inline LatticeIsometry edge_transition(const DelaunayView& D, int h,
                                       const std::array<Eisenstein, 3>& Pf,
                                       const std::array<Eisenstein, 3>& Pg) {
  const int k  = D.cycle_slot(h);           // = the chart corner slot the
  const int k2 = D.cycle_slot(D.twin(h));   //   edge leaves from (gated)
  return isometry_from_segments(Pf[k], Pf[(k + 1) % 3], Pg[(k2 + 1) % 3], Pg[k2]);
}

// Per-half-edge chart transitions over a parametrization's charts.
// Keyed by the crossing half-edge, so a cell pair sharing several edges
// (legal on a delta-complex) is unambiguous by construction.  Dead slots
// stay identity and are unreachable (@inv: every half-edge on a live
// face cycle is alive).
std::vector<LatticeIsometry> half_edge_transitions(const DelaunayView& D,
                                                   const ParamTablesView& V) {
  std::vector<LatticeIsometry> trans(D.nh);
  for (int h = 0; h < D.nh; h++) {
    if (!D.alive(h)) continue;
    trans[h] = edge_transition(D, h, V.frame_points(D.he_face[h]),
                               V.frame_points(D.he_face[D.twin(h)]));
  }
  return trans;
}

// The T_sorted face/edge combinatorics: faces CCW, and the DENSE
// arcid-indexed tables (arcid = u*dmax + arc index of v in u's list --
// the DCEL build's own vocabulary): arc -> face-left, arc -> undirected
// edge id, plus the edge list and each edge's two incident faces.
struct FaceEdgeTables {
  std::vector<tri_t>             tface;
  std::vector<int32_t>           arc_face;    // [N*dmax]
  std::vector<int32_t>           arc_edge;    // [N*dmax]
  std::vector<std::array<int,2>> tedge;
  std::vector<std::array<int,2>> edge_faces;
};

FaceEdgeTables face_edge_tables(const TriangulationView& T) {
  FaceEdgeTables tab;
  tab.tface = T.triangles();
  tab.arc_face.assign(arc_space(T), -1);   // the lib's named arc-index space
  tab.arc_edge.assign(arc_space(T), -1);
  // @pre (u, v) is a T arc at every use below (consecutive face corners
  // and enumerated adjacencies) -- the guard-free form is correct here.
  auto arcid = [&](int u, int v) { return u * T.dmax + T.arc_ix(u, v); };
  for (int i = 0; i < (int)tab.tface.size(); i++) {
    const tri_t& t = tab.tface[i];
    tab.arc_face[arcid(t[0], t[1])] = i;
    tab.arc_face[arcid(t[1], t[2])] = i;
    tab.arc_face[arcid(t[2], t[0])] = i;
  }
  for (int u = 0; u < T.N; u++)
    for (int v : T[u])
      if (u < v) {
        const int id = (int)tab.tedge.size();
        tab.tedge.push_back({ u, v });
        tab.arc_edge[arcid(u, v)] = id;
        tab.arc_edge[arcid(v, u)] = id;
        const int fuv = tab.arc_face[arcid(u, v)];
        const int fvu = tab.arc_face[arcid(v, u)];
        if (fuv < 0 || fvu < 0)
          fail("T_sorted arc without a face (graph not closed/oriented)");
        tab.edge_faces.push_back({ fuv, fvu });
      }
  return tab;
}

// Anchored edges: two ADJACENT lattice points claimed by one cell.  Their
// unit segment lies inside the closed convex cell, so their relative chart
// positions are ground truth (the Anchoring lemma).  Iteration is cells in
// id order and entries in SCAN order, so anchor selection -- and hence
// every resolved chain -- is deterministic and platform-stable.
struct AnchoredEdges {
  std::vector<int32_t>    anchor_of_edge;   // edge id -> anchors index, or -1
  std::vector<AnchorEdge> anchors;
};

AnchoredEdges anchored_edges(const ParamTablesView& V, const TriangulationView& T,
                             const FaceEdgeTables& tab) {
  AnchoredEdges out;
  out.anchor_of_edge.assign(tab.tedge.size(), -1);
  const Eisenstein half_dirs[3] = { Eisenstein(1, 0), Eisenstein(0, 1), Eisenstein(-1, 1) };
  for (int f = 0; f < V.nf; f++) {
    if (!V.cell_live(f)) continue;
    for (const LatticePoint& e : V.cell_entries(f)) {
      const Eisenstein p = e.pos();
      for (const Eisenstein d : half_dirs) {
        const LatticePoint* q = V.claim(f, p + d);
        if (!q) continue;
        const int u = e.vid, v = q->vid;
        const int ix = T.arc_ix(u, v);
        if (ix < 0)
          fail("cell " + std::to_string(f) + " developed vertices "
               + std::to_string(u) + " and " + std::to_string(v)
               + " -- which are NOT a mesh edge -- onto adjacent lattice positions: "
               "its iDT geodesic triangle does not embed flat (a folded development). "
               "This is a residual non-embedding of certain obtuse iDT faces, seen on "
               "both simplicial and non-simplicial raw iDTs; it is NOT specific to "
               "multi-edges. Parametrize the Alexandrov-realized iDT "
               "(eisenstein_paint::realize_dual) instead of the raw dual_idt "
               "for this isomer");
        const int eid = tab.arc_edge[u * T.dmax + ix];
        if (out.anchor_of_edge[eid] >= 0) continue;
        out.anchor_of_edge[eid] = (int32_t)out.anchors.size();
        out.anchors.push_back({ f, u, v, p, p + d });
      }
    }
  }
  return out;
}

// The routing forest: a multi-source BFS over the edge graph (edges
// adjacent iff they share a face) from the anchored edges, so every
// T_sorted edge has a parent chain of via-face midpoint hops back to
// ground truth.
struct AnchorRouting {
  std::vector<int32_t> parent_edge;   // edge id -> parent edge (-1 root, from an anchor)
  std::vector<int32_t> via_face;      // edge id -> face shared with parent
  std::vector<int32_t> depth;         // edge id -> chain length back to its anchor
};

AnchorRouting route_to_anchors(const TriangulationView& T,
                               const FaceEdgeTables& tab,
                               const std::vector<int32_t>& anchor_of_edge) {
  const int ne = (int)tab.tedge.size();
  AnchorRouting out;
  out.parent_edge.assign(ne, -2);   // -2 unvisited, -1 root
  out.via_face.assign(ne, -1);
  out.depth.assign(ne, 0);
  std::vector<int> queue;
  queue.reserve(ne);
  for (int e = 0; e < ne; e++)
    if (anchor_of_edge[e] >= 0) { out.parent_edge[e] = -1; queue.push_back(e); }
  if (queue.empty()) fail("no anchored edge anywhere (no cell claims two adjacent vertices)");
  auto face_edge = [&](int fi, int k) {
    const tri_t& t = tab.tface[fi];
    const int u = t[k], v = t[(k + 1) % 3];
    return tab.arc_edge[u * T.dmax + T.arc_ix(u, v)];
  };
  for (size_t qi = 0; qi < queue.size(); qi++) {
    const int e = queue[qi];
    for (const int fi : tab.edge_faces[e])
      for (int k = 0; k < 3; k++) {
        const int e2 = face_edge(fi, k);
        if (out.parent_edge[e2] != -2) continue;
        out.parent_edge[e2] = e;
        out.via_face[e2] = fi;
        out.depth[e2] = out.depth[e] + 1;
        queue.push_back(e2);
      }
  }
  for (int e = 0; e < ne; e++)
    if (out.parent_edge[e] == -2) fail("edge graph not covered by anchored-edge BFS");
  return out;
}

}  // namespace

CellAtlas build_atlas(const SurfaceParametrization& P) {
  CellAtlas A;
  A.V = P.view();
  A.D = P.D;
  A.T = P.T;
  check_charts_match_cycles(A.D, A.V);
  A.he_trans = half_edge_transitions(A.D, A.V);
  auto tab   = face_edge_tables(A.T);
  auto anch  = anchored_edges(A.V, A.T, tab);
  auto route = route_to_anchors(A.T, tab, anch.anchor_of_edge);

  A.tface          = std::move(tab.tface);
  A.arc_face       = std::move(tab.arc_face);
  A.arc_edge       = std::move(tab.arc_edge);
  A.tedge          = std::move(tab.tedge);
  A.edge_faces     = std::move(tab.edge_faces);
  A.anchor_of_edge = std::move(anch.anchor_of_edge);
  A.anchors        = std::move(anch.anchors);
  A.bfs_parent_edge = std::move(route.parent_edge);
  A.bfs_via_face    = std::move(route.via_face);
  A.bfs_depth       = std::move(route.depth);
  A.resolved.assign(A.tface.size(), {});
  return A;
}

// The crossing cap: a valid trace crosses each cell a bounded number of
// times (the crossing parameter never decreases and strictly increases
// off the cones), so this is a deep-invariant guard, not a bound any
// valid trace approaches.
inline int max_cell_crossings(int nf) { return 4 * nf + 64; }

CellTrace trace_segment(const CellAtlas& A, int cell,
                        EisensteinRational a, EisensteinRational b) {
  // Common denominator l for BOTH endpoints (raw, deliberately un-reduced:
  // the crossing parameter compares the two wedge values, which is only
  // meaningful on one shared scale).  Narrow-checked ONCE; everything
  // below works on the checked int.
  const int l = narrow(std::lcm(a.den, b.den), "common denominator");
  Eisenstein an = a.num * narrow(l / a.den, "denominator scale");
  Eisenstein bn = b.num * narrow(l / b.den, "denominator scale");
  LatticeIsometry carried;   // identity
  __int128 t_num = 0, t_den = 1;   // progress along the segment (rational)
  const int cap = max_cell_crossings(A.V.nf);
  for (int step = 0; step < cap; step++) {
    const std::array<Eisenstein, 3> P = A.V.frame_points(cell);
    if (inside_cell(P, bn, l)) return { cell, EisensteinRational(bn, l), carried };
    int best = -1;
    __int128 bt_num = 0, bt_den = 1;
    for (int k = 0; k < 3; k++) {
      const Eisenstein Pi = P[k], Pj = P[(k + 1) % 3];
      const Eisenstein ev = Pj - Pi;
      auto w = [&](const Eisenstein& num) -> __int128 {
        const Eisenstein rel = num - Pi * l;
        return (__int128)ev.first * rel.second - (__int128)ev.second * rel.first;
      };
      const __int128 wa = w(an), wb = w(bn);
      if (wb >= 0) continue;                    // b not strictly beyond this edge
      if (wa < 0) continue;                     // a beyond it too (cannot exit through it)
      const __int128 cn = wa, cd = wa - wb;     // t_cross = wa / (wa - wb), cd > 0
      // Progress filter is NON-strict (an equal-parameter crossing is
      // admitted): a segment through a corner may re-cross at the same
      // parameter.  Strict increase -- hence termination well inside the
      // cap -- holds whenever the open segment meets no cone, which every
      // internal caller guarantees (anchored segments, face midlines,
      // closed-face samples).
      if (cn * t_den < t_num * cd) continue;    // behind current progress
      if (best < 0 || cn * bt_den < bt_num * cd) { best = k; bt_num = cn; bt_den = cd; }
    }
    if (best < 0) fail("segment trace found no exit edge (point location stuck)");
    // The exit arc IS the best-th half-edge of the face cycle (corner k =
    // origin of cycle half-edge k, gated by build_atlas), so the crossing
    // transition is he_trans of exactly that half-edge -- no corner-pair
    // search, and parallel edges of a delta-complex cell pair cannot be
    // confused.
    const int h_exit = A.D.face_halfedges(cell)[best];
    const int ncell = A.D.he_face[A.D.twin(h_exit)];
    const LatticeIsometry T = A.he_trans[h_exit];
    an = apply_raw(T, an, l);
    bn = apply_raw(T, bn, l);
    carried = T * carried;
    cell = ncell;
    t_num = bt_num;
    t_den = bt_den;
  }
  fail("segment trace exceeded " + std::to_string(cap) + " cell crossings");
}

// The face edge with the shallowest routing chain back to an anchor
// (depths come from the BFS itself; k-ascending tie-break).
static int shallowest_face_edge(const AtlasView& AV, const tri_t& t) {
  int e_target = -1, best_depth = INT32_MAX;
  for (int k = 0; k < 3; k++) {
    const int e = AV.edge_of(t[k], t[(k + 1) % 3]);
    if (AV.bfs_depth[e] < best_depth) { best_depth = AV.bfs_depth[e]; e_target = e; }
  }
  return e_target;
}

// The root-first parent chain of edge e (root = an anchored edge).
static std::vector<int> anchor_chain(const AtlasView& AV, int e) {
  std::vector<int> chain;
  for (int cur = e; cur >= 0; cur = AV.bfs_parent_edge[cur]) chain.push_back(cur);
  std::reverse(chain.begin(), chain.end());
  return chain;
}

const ResolvedFace& resolve_face(CellAtlas& A, int fi) {
  ResolvedFace& RS = A.resolved[fi];
  if (RS.ok) return RS;

  const AtlasView AV = A.view();
  const tri_t& t = A.tface[fi];
  const std::vector<int> chain = anchor_chain(AV, shallowest_face_edge(AV, t));

  const AnchorEdge& AE = A.anchors[A.anchor_of_edge[chain[0]]];
  int  cell = AE.cell;
  int  eu = AE.u, ev = AE.v;
  Eisenstein pu = AE.pu, pv = AE.pv;
  EisensteinRational mid(pu + pv, 2);
  if (!inside_cell(A.V.frame_points(cell), mid.num,
                   narrow(mid.den, "anchor denominator")))
    fail("anchor midpoint not inside its cell");

  // Midpoint hops: each hop's segment is a midline of the via-face -- inside
  // one flat face -- so its exact trace never wraps a cone.
  for (size_t ci = 1; ci < chain.size(); ci++) {
    const int e_next = chain[ci];
    const int via = A.bfs_via_face[e_next];
    const tri_t& ft = A.tface[via];
    const int k_arc = corner_index_of_edge(ft, eu, ev);
    if (k_arc < 0) fail("hop via-face does not contain the current edge");
    const std::array<Eisenstein, 3> p = develop_face_on_edge(ft, k_arc, eu, pu, pv);
    auto pos_of = [&](int vid) -> Eisenstein {
      for (int k = 0; k < 3; k++)
        if (ft[k] == vid) return p[k];
      fail("hop target edge not on the via-face");
    };
    const auto& E2 = A.tedge[e_next];
    const Eisenstein qa = pos_of(E2[0]), qb = pos_of(E2[1]);
    const EisensteinRational nmid(qa + qb, 2);
    CellTrace out = trace_segment(A, cell, mid, nmid);
    cell = out.cell;
    mid  = out.pos;
    eu = E2[0];
    ev = E2[1];
    pu = out.carried.apply(qa);
    pv = out.carried.apply(qb);
  }

  const int k_arc = corner_index_of_edge(t, eu, ev);
  if (k_arc < 0) fail("resolved edge not on its face");
  const std::array<Eisenstein, 3> p = develop_face_on_edge(t, k_arc, eu, pu, pv);
  // exact verification: unit CCW triangle (catches any orientation-convention drift)
  for (int k = 0; k < 3; k++)
    if ((p[(k + 1) % 3] - p[k]).norm2() != 1) fail("resolved face is not a unit triangle");
  if (wedge(p[1] - p[0], p[2] - p[0]) != 1)
    fail("resolved face is not CCW in the lattice frame");

  RS.ok = true;
  RS.cell = cell;
  RS.P = p;
  RS.anchor = mid;
  return RS;
}

void resolve_all(CellAtlas& A) {
  for (int fi = 0; fi < (int)A.tface.size(); fi++) resolve_face(A, fi);
}

CellPoint locate_sample(CellAtlas& A, int fi, long a, long b, long c, long den) {
  // @pre checked: the sample must lie in the CLOSED face (affine,
  // non-negative weights) -- only there does the anchor-to-sample trace
  // stay in flat cone-free territory.
  if (a < 0 || b < 0 || c < 0 || den <= 0 || a + b + c != den)
    fail("locate_sample weights must be non-negative with a + b + c == den");
  const ResolvedFace& R = resolve_face(A, fi);
  const Eisenstein num = R.P[0] * narrow(a, "sample weight")
                       + R.P[1] * narrow(b, "sample weight")
                       + R.P[2] * narrow(c, "sample weight");
  const EisensteinRational q(num, den);
  CellTrace out = trace_segment(A, R.cell, R.anchor, q);
  return { out.cell, out.pos };
}

// =====================================================================
// The T-free (intrinsic) atlas: development selection by cross-edge
// consistency, and the dual reconstruction that selection carries.
// =====================================================================

namespace {

// Do the two charts currently accepted on the faces of h develop their
// shared edge into one another?  edge_transition throws exactly for the
// mirror (wrong-orbit) development, so a clean return IS the verdict --
// the falsifier that replaces T's frame walkers.
bool transition_consistent(const IntrinsicAtlas& A, int h) {
  const DelaunayView& D = A.D;
  try {
    (void)edge_transition(D, h, A.frame_points(D.he_face[h]),
                          A.frame_points(D.he_face[D.twin(h)]));
    return true;
  } catch (const std::logic_error&) {
    return false;
  }
}

}  // namespace

IntrinsicAtlas build_intrinsic_atlas(const DelaunayView& D) {
  IntrinsicAtlas A;
  A.D = D;
  // A pure flat cone metric: every D vertex IS a cone, and the two counts
  // enter only the scratch-capacity formulas the intrinsic path never reads.
  A.dev = cell_developments(D, D.nv, D.nv);
  A.accepted.assign(D.nf, -1);

  int nlive = 0;
  for (int f = 0; f < D.nf; f++) if (A.dev.n_developments(f) > 0) nlive++;

  // Choose the seed's development so the whole atlas propagates.  A cell with
  // several developments has only ONE that is globally consistent, and an
  // arbitrary pick can wedge the entire BFS.  So prefer an unambiguous
  // (single-development) seed, and otherwise try each of the seed's
  // developments until one places every live cell; the BFS then resolves every
  // remaining ambiguous cell by transition consistency.
  int seed = -1;
  for (int f = 0; f < D.nf; f++) if (A.dev.n_developments(f) == 1) { seed = f; break; }
  if (seed < 0)
    for (int f = 0; f < D.nf; f++) if (A.dev.n_developments(f) > 0) { seed = f; break; }
  if (seed < 0) fail("build_intrinsic_atlas: no placeable cell in iDT");

  for (int k0 = 0; k0 < A.dev.n_developments(seed); k0++) {
    A.accepted.assign(D.nf, -1);
    A.accepted[seed] = k0;
    std::vector<int> order{ seed };
    for (size_t qi = 0; qi < order.size(); qi++) {
      const int f = order[qi];
      for (const int h : D.face_halfedges(f)) {
        const int g = D.he_face[D.twin(h)];
        if (A.placed(g)) continue;
        for (int k = 0; k < A.dev.n_developments(g); k++) {
          A.accepted[g] = k;
          if (transition_consistent(A, h)) { order.push_back(g); break; }
          A.accepted[g] = -1;   // wrong-orbit development; try the next
        }
      }
    }
    if ((int)order.size() == nlive) break;   // consistent development found
  }

  // The transitions the selection makes exact.  Edges touching an unplaced
  // cell (a folded cell on a noisy metric) simply have none; between two
  // placed cells the transition must exist, so a refusal here is loud.
  A.he_trans.assign(D.nh, LatticeIsometry{});
  for (int h = 0; h < D.nh; h++) {
    if (!D.alive(h)) continue;
    const int f = D.he_face[h], g = D.he_face[D.twin(h)];
    if (!A.placed(f) || !A.placed(g)) continue;
    A.he_trans[h] = edge_transition(D, h, A.frame_points(f), A.frame_points(g));
  }

  // Global numbering of the accepted charts' lattice points (the union-find
  // node ids intrinsic_dual glues on).
  A.node_first.assign(D.nf + 1, 0);
  for (int f = 0; f < D.nf; f++)
    A.node_first[f + 1] =
        A.node_first[f] + (A.placed(f) ? A.chart(f).scan.n_entries : 0);
  return A;
}

Triangulation intrinsic_dual(const IntrinsicAtlas& A) {
  const DelaunayView& D = A.D;

  // Unite the copies of each shared lattice point: along every edge, walk its
  // integer points (gcd of the edge vector's components + 1 of them) and map
  // each from one cell's frame into the other's by the chart transition.
  UnionFind uf(A.n_nodes());
  auto node_at = [&](int cell, Eisenstein p) {
    const int n = A.node_at(cell, p);
    if (n < 0) fail("intrinsic_dual: shared lattice point absent from cell scan");
    return n;
  };
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    const int f = D.he_face[h], g = D.he_face[D.twin(h)];
    if (!A.placed(f) || !A.placed(g)) continue;
    const std::array<Eisenstein, 3> Pf = A.frame_points(f);
    const int k = D.cycle_slot(h);       // the crossed edge's corner slot in f
    const Eisenstein Pa = Pf[k], Pb = Pf[(k + 1) % 3];
    const Eisenstein e = Pb - Pa;
    const int segs = std::gcd(std::abs(e.first), std::abs(e.second));
    const Eisenstein step = e / segs;
    const LatticeIsometry Tfg = A.he_trans[h];
    for (int i = 0; i <= segs; i++) {
      const Eisenstein pf = Pa + step * i;
      uf.unite(node_at(f, pf), node_at(g, Tfg.apply(pf)));
    }
  }

  // Dense vertex id per union-find class.
  const std::vector<std::vector<int>> comps = uf.components();
  std::vector<int> vid(A.n_nodes());
  for (int c = 0; c < (int)comps.size(); c++)
    for (const int u : comps[c]) vid[u] = c;
  const int V = (int)comps.size();
  auto gid = [&](int cell, Eisenstein p) { return vid[node_at(cell, p)]; };

  // Locate the apex of the unit triangle (v, w, apex): it is unit_apex(pv, pw)
  // in cell c's frame, but that triangle may straddle one or more cone-iDT
  // edges.  Walk the straight segment v -> apex through the cells: while the
  // apex is outside the current cell, cross the edge the segment truly exits --
  // the apex strictly exterior to it AND the crossing WITHIN that edge's
  // segment (s in [0,1]), not merely on its infinite line -- re-expressing v
  // and the apex by the chart transition, until the apex lands in a cell that
  // claims it.  The oriented multi-sheet fold, boundary by boundary.  The exit
  // arc IS the best-th half-edge of the face cycle (corner k = origin of cycle
  // half-edge k -- true by construction here: cell_developments takes each
  // cell's corners straight off its face cycle).
  struct ApexLoc { int cell; Eisenstein pv, apex; };
  const int cap = max_cell_crossings(D.nf);
  auto locate_apex = [&](int c, Eisenstein pv, Eisenstein apex) -> ApexLoc {
    for (int step = 0; step < cap; step++) {
      if (A.node_at(c, apex) >= 0) return { c, pv, apex };
      const std::array<Eisenstein, 3> P = A.frame_points(c);
      int best = -1;
      for (int k = 0; k < 3; k++) {
        const Eisenstein Pi = P[k], Pj = P[(k + 1) % 3], d = Pj - Pi;
        if (wedge(d, apex - Pi) >= 0) continue;     // apex not strictly exterior to edge k
        const long sden = wedge(apex - pv, d);      // crossing s = snum/sden along [Pi, Pj]
        const long snum = wedge(apex - pv, pv - Pi);
        const bool in_edge = sden > 0 ? (snum >= 0 && snum <= sden)
                                      : sden < 0 && (snum <= 0 && snum >= sden);
        if (in_edge) { best = k; break; }
      }
      if (best < 0) fail("intrinsic_dual: no exit edge toward apex");
      const int h_exit = D.face_halfedges(c)[best];
      const int g = D.he_face[D.twin(h_exit)];
      if (!A.placed(g)) fail("intrinsic_dual: no placed neighbour across exit edge");
      const LatticeIsometry T = A.he_trans[h_exit];
      pv = T.apply(pv);  apex = T.apply(apex);  c = g;
    }
    fail("intrinsic_dual: apex walk exceeded " + std::to_string(cap) + " cell crossings");
  };

  // Seed each vertex with one incident dual edge (v at pv in cell c, neighbour
  // at pw).  The fast path is an in-cell unit neighbour; a cone in a coarse
  // cell can have ALL its neighbours across boundaries, so fall back to
  // stepping a unit in each lattice direction and locating where it lands
  // (skipping the direction into the cone's angle deficit, which reaches no
  // neighbour).
  struct Seed { int cell = -1; Eisenstein pv, pw; };
  std::vector<Seed> seed(V);
  std::vector<std::pair<int, Eisenstein>> occ(V, { -1, Eisenstein{} });
  for (int f = 0; f < D.nf; f++) {
    if (!A.placed(f)) continue;
    const DevelopmentView ch = A.chart(f);
    for (int b = ch.scan.b_min; b <= ch.scan.b_max; b++) {
      const ScanRow row = ch.rows[b - ch.scan.b_min];
      for (int a = row.a_left; a <= row.a_right; a++) {
        const Eisenstein p(a, b);
        const int v = gid(f, p);
        if (occ[v].first < 0) occ[v] = { f, p };
        if (seed[v].cell >= 0) continue;
        for (int d = 0; d < 6; d++)
          if (A.node_at(f, p + unit_direction(d)) >= 0) {
            seed[v] = { f, p, p + unit_direction(d) };
            break;
          }
      }
    }
  }
  for (int v = 0; v < V; v++) {
    const auto [c, pv] = occ[v];
    for (int d = 0; d < 6 && seed[v].cell < 0; d++)
      try {
        const ApexLoc nb = locate_apex(c, pv, pv + unit_direction(d));
        if (gid(nb.cell, nb.apex) != v) seed[v] = { nb.cell, nb.pv, nb.apex };
      } catch (const std::logic_error&) { /* direction into the deficit; try the next */ }
    if (seed[v].cell < 0)
      fail("intrinsic_dual: vertex " + std::to_string(v) + " has no incident edge");
  }

  // Read each vertex's CCW neighbour ring by rotating around it: after
  // neighbour w, the next one CCW is the apex of the unit triangle (v, w, apex)
  // on the CCW side.  The ring closes when the rotation returns to the start
  // neighbour -- after 5 steps at a cone, 6 at a hex: the degree is read off,
  // never assumed.
  std::vector<std::vector<node_t>> nbr(V);
  for (int v = 0; v < V; v++) {
    int c = seed[v].cell;
    Eisenstein pv = seed[v].pv, pw = seed[v].pw;
    const int w0 = gid(c, pw);
    for (;;) {
      nbr[v].push_back(gid(c, pw));
      if ((int)nbr[v].size() > V)
        fail("intrinsic_dual: neighbour ring failed to close at vertex "
             + std::to_string(v));
      const ApexLoc nx = locate_apex(c, pv, unit_apex(pv, pw));
      c = nx.cell;  pv = nx.pv;  pw = nx.apex;
      if (gid(c, pw) == w0) break;
    }
  }

  return Triangulation(Graph(Spanify::OwnedDenseGraph<node_t>(nbr)));
}

}  // namespace eisenstein_paint
