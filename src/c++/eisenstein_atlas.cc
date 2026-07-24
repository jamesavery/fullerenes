// eisenstein_atlas.cc -- Tier 3.5: the exact cross-cell atlas (see the header
// banner).  Promoted from claude-projects/curvature-flow/src/cage_transfer.cc
// (stage A of the tab-4 cage transfer), where it was validated by the
// exactness oracle: transferring a candidate onto its own equilateral complex
// reproduces eisenstein_paint::run's dual-vertex positions to ~1e-15.
//
// All arithmetic in this file is integer (long / __int128) or rational over
// integers -- no floating point enters any predicate or position.

#include "fullerenes/eisenstein_atlas.hh"

#include <algorithm>
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

// q = num/den inside the CLOSED cell (all three edge half-planes)?
// Cross-multiplied exact predicate; __int128 keeps den * coordinate exact.
inline bool inside_cell(const AtlasCell& R, Eisenstein num, long den) {
  for (int k = 0; k < 3; k++) {
    const Eisenstein Pi = R.P[k], Pj = R.P[(k + 1) % 3];
    const Eisenstein ev = Pj - Pi;
    const Eisenstein rel = num - Pi * narrow(den, "denominator");
    if ((__int128)ev.first * rel.second - (__int128)ev.second * rel.first < 0) return false;
  }
  return true;
}

// T applied to the raw rational num/den WITHOUT re-normalizing: trace_segment
// keeps its two endpoints on one common denominator, which a canonicalizing
// constructor would destroy.  (The missing vocabulary word here is a
// "two points on one shared denominator" representation; coin it if a third
// use appears.)
inline Eisenstein apply_raw(const LatticeIsometry& T, Eisenstein num, long den) {
  const Eisenstein n = T.reflect ? num.complex_conj() : num;
  return T.t * narrow(den, "shared denominator") + T.u * n;
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
// build_atlas: a composition of five value-producing constructions.  Each
// takes exactly the values it depends on, so the dependency order is the
// argument flow -- enforced by the compiler, not by call-order convention.
// ---------------------------------------------------------------------------

// The hash index of the parametrization's charts: an AtlasCell is the
// chart's corner identity plus its claim map re-indexed for O(1) position
// lookup.  No chart is recomputed -- the walkers and scans ran once, in
// parametrize().
std::vector<AtlasCell> index_charts(const SurfaceParametrization& P) {
  const DelaunayTriangulation& D = *P.D;
  std::vector<AtlasCell> cells(P.cells.size());
  for (int f = 0; f < (int)P.cells.size(); f++) {
    if (D.f_he[f] < 0) continue;
    const Cell& C = P.cells[f];
    if (!C.ok) fail("live cell " + std::to_string(f) + " not charted");
    AtlasCell& R = cells[f];
    R.ok = true;
    R.corners = C.corners;
    R.P       = C.P;
    const LatticeMap& lm = P.lmaps[f];
    R.claim.reserve(lm.entries.size() * 2);
    for (const auto& [pos, vid] : lm.entries)
      R.claim.emplace(pos, vid);
  }
  return cells;
}

// The chart transition across every D edge, both directions: the unique
// orientation-preserving lattice isometry matching the shared edge's cone
// corners across the two charts (Lemma: transitions are exact).
std::unordered_map<long long, LatticeIsometry>
chart_transitions(const DelaunayTriangulation& D,
                  const std::vector<AtlasCell>& cells) {
  std::unordered_map<long long, LatticeIsometry> trans;
  auto corner_pos = [&](const AtlasCell& R, int cid) -> Eisenstein {
    for (int k = 0; k < 3; k++)
      if (R.corners[k] == cid) return R.P[k];
    fail("cone " + std::to_string(cid) + " is not a corner of its cell");
  };
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    const int f = D.he_face[h], g = D.he_face[D.twin(h)];
    const int a = D.he_origin[h], b = D.dest(h);
    const AtlasCell& Rf = cells[f];
    const AtlasCell& Rg = cells[g];
    trans[CellAtlas::arc_key(g, f)] =
        isometry_from_segments(corner_pos(Rg, a), corner_pos(Rg, b),
                               corner_pos(Rf, a), corner_pos(Rf, b));
    trans[CellAtlas::arc_key(f, g)] =
        isometry_from_segments(corner_pos(Rf, a), corner_pos(Rf, b),
                               corner_pos(Rg, a), corner_pos(Rg, b));
  }
  return trans;
}

// The T_sorted face/edge combinatorics: faces CCW, arc -> face-left,
// undirected edges with ids, and each edge's two incident faces.
struct FaceEdgeTables {
  std::vector<tri_t>                 tface;
  std::unordered_map<long long, int> arc_face;
  std::vector<std::array<int, 2>>    tedge;
  std::unordered_map<long long, int> edge_id;
  std::vector<std::array<int, 2>>    edge_faces;
};

FaceEdgeTables face_edge_tables(const TriangulationView& T) {
  FaceEdgeTables tab;
  tab.tface = T.triangles();
  for (int i = 0; i < (int)tab.tface.size(); i++) {
    const tri_t& t = tab.tface[i];
    tab.arc_face[CellAtlas::arc_key(t[0], t[1])] = i;
    tab.arc_face[CellAtlas::arc_key(t[1], t[2])] = i;
    tab.arc_face[CellAtlas::arc_key(t[2], t[0])] = i;
  }
  for (int u = 0; u < T.N; u++)
    for (int v : T[u])
      if (u < v) {
        const int id = (int)tab.tedge.size();
        tab.tedge.push_back({ u, v });
        tab.edge_id[CellAtlas::edge_key(u, v)] = id;
        const auto itf = tab.arc_face.find(CellAtlas::arc_key(u, v));
        const auto itg = tab.arc_face.find(CellAtlas::arc_key(v, u));
        if (itf == tab.arc_face.end() || itg == tab.arc_face.end())
          fail("T_sorted arc without a face (graph not closed/oriented)");
        tab.edge_faces.push_back({ itf->second, itg->second });
      }
  return tab;
}

// Anchored edges: two ADJACENT lattice points claimed by one cell.  Their
// unit segment lies inside the closed convex cell, so their relative chart
// positions are ground truth (the Anchoring lemma).
struct AnchoredEdges {
  std::vector<int>                        anchor_of_edge;   // edge id -> anchors index, or -1
  std::vector<CellAtlas::AnchorEdge>      anchors;
};

AnchoredEdges anchored_edges(const std::vector<AtlasCell>& cells,
                             const FaceEdgeTables& tab) {
  AnchoredEdges out;
  out.anchor_of_edge.assign(tab.tedge.size(), -1);
  const Eisenstein half_dirs[3] = { Eisenstein(1, 0), Eisenstein(0, 1), Eisenstein(-1, 1) };
  for (int f = 0; f < (int)cells.size(); f++) {
    const AtlasCell& R = cells[f];
    if (!R.ok) continue;
    for (const auto& [p, vid] : R.claim) {
      for (const Eisenstein d : half_dirs) {
        const auto it = R.claim.find(p + d);
        if (it == R.claim.end()) continue;
        const int u = vid, v = it->second;
        const auto eit = tab.edge_id.find(CellAtlas::edge_key(u, v));
        if (eit == tab.edge_id.end())
          fail("cell " + std::to_string(f) + " developed vertices "
               + std::to_string(u) + " and " + std::to_string(v)
               + " -- which are NOT a mesh edge -- onto adjacent lattice positions: "
               "its iDT geodesic triangle does not embed flat (a folded development). "
               "This is a residual non-embedding of certain obtuse iDT faces, seen on "
               "both simplicial and non-simplicial raw iDTs; it is NOT specific to "
               "multi-edges. Parametrize the Alexandrov-realized iDT "
               "(eisenstein_paint::realize_dual) instead of the raw dual_idt "
               "for this isomer");
        if (out.anchor_of_edge[eit->second] >= 0) continue;
        out.anchor_of_edge[eit->second] = (int)out.anchors.size();
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
  std::vector<int> parent_edge;   // edge id -> parent edge (-1 root, from an anchor)
  std::vector<int> via_face;      // edge id -> face shared with parent
};

AnchorRouting route_to_anchors(const FaceEdgeTables& tab,
                               const std::vector<int>& anchor_of_edge) {
  const int ne = (int)tab.tedge.size();
  AnchorRouting out;
  out.parent_edge.assign(ne, -2);   // -2 unvisited, -1 root
  out.via_face.assign(ne, -1);
  std::vector<int> queue;
  queue.reserve(ne);
  for (int e = 0; e < ne; e++)
    if (anchor_of_edge[e] >= 0) { out.parent_edge[e] = -1; queue.push_back(e); }
  if (queue.empty()) fail("no anchored edge anywhere (no cell claims two adjacent vertices)");
  auto face_edge = [&](int fi, int k) {
    const tri_t& t = tab.tface[fi];
    return tab.edge_id.at(CellAtlas::edge_key(t[k], t[(k + 1) % 3]));
  };
  for (size_t qi = 0; qi < queue.size(); qi++) {
    const int e = queue[qi];
    for (const int fi : tab.edge_faces[e])
      for (int k = 0; k < 3; k++) {
        const int e2 = face_edge(fi, k);
        if (out.parent_edge[e2] != -2) continue;
        out.parent_edge[e2] = e;
        out.via_face[e2] = fi;
        queue.push_back(e2);
      }
  }
  for (int e = 0; e < ne; e++)
    if (out.parent_edge[e] == -2) fail("edge graph not covered by anchored-edge BFS");
  return out;
}

}  // namespace

LatticeIsometry CellAtlas::transition(int f_from, int f_to) const {
  const auto it = trans.find(arc_key(f_from, f_to));
  if (it == trans.end())
    throw std::logic_error("eisenstein_atlas: missing cell transition "
                           + std::to_string(f_from) + " -> " + std::to_string(f_to));
  return it->second;
}

CellAtlas build_atlas(const SurfaceParametrization& P) {
  auto cells = index_charts(P);
  auto trans = chart_transitions(*P.D, cells);
  auto tab   = face_edge_tables(P.T);
  auto anch  = anchored_edges(cells, tab);
  auto route = route_to_anchors(tab, anch.anchor_of_edge);

  CellAtlas A;
  A.D              = P.D;
  A.T              = P.T;
  A.cells          = std::move(cells);
  A.trans          = std::move(trans);
  A.tface          = std::move(tab.tface);
  A.arc_face       = std::move(tab.arc_face);
  A.tedge          = std::move(tab.tedge);
  A.edge_id        = std::move(tab.edge_id);
  A.edge_faces     = std::move(tab.edge_faces);
  A.anchor_of_edge = std::move(anch.anchor_of_edge);
  A.anchors        = std::move(anch.anchors);
  A.bfs_parent_edge = std::move(route.parent_edge);
  A.bfs_via_face    = std::move(route.via_face);
  A.resolved.assign(A.tface.size(), {});
  return A;
}

CellTrace trace_segment(const CellAtlas& A, int cell,
                        EisensteinRational a, EisensteinRational b) {
  // Common denominator l for BOTH endpoints (raw, deliberately un-reduced:
  // the crossing parameter compares the two wedge values, which is only
  // meaningful on one shared scale).
  const long l = std::lcm(a.den, b.den);
  Eisenstein an = a.num * narrow(l / a.den, "denominator scale");
  Eisenstein bn = b.num * narrow(l / b.den, "denominator scale");
  LatticeIsometry carried;   // identity
  __int128 t_num = 0, t_den = 1;   // progress along the segment (rational)
  const int cap = 4 * (int)A.cells.size() + 64;
  for (int step = 0; step < cap; step++) {
    const AtlasCell& R = A.cells[cell];
    if (inside_cell(R, bn, l)) return { cell, EisensteinRational(bn, l), carried };
    int best = -1;
    __int128 bt_num = 0, bt_den = 1;
    for (int k = 0; k < 3; k++) {
      const Eisenstein Pi = R.P[k], Pj = R.P[(k + 1) % 3];
      const Eisenstein ev = Pj - Pi;
      auto w = [&](const Eisenstein& num) -> __int128 {
        const Eisenstein rel = num - Pi * (int)l;   // l checked above
        return (__int128)ev.first * rel.second - (__int128)ev.second * rel.first;
      };
      const __int128 wa = w(an), wb = w(bn);
      if (wb >= 0) continue;                    // b not strictly beyond this edge
      if (wa < 0) continue;                     // a beyond it too (cannot exit through it)
      const __int128 cn = wa, cd = wa - wb;     // t_cross = wa / (wa - wb), cd > 0
      if (cn * t_den < t_num * cd) continue;    // behind current progress
      if (best < 0 || cn * bt_den < bt_num * cd) { best = k; bt_num = cn; bt_den = cd; }
    }
    if (best < 0) fail("segment trace found no exit edge (point location stuck)");
    const int ca = R.corners[best], cb = R.corners[(best + 1) % 3];
    const DelaunayTriangulation& D = *A.D;
    int found = -1;
    for (const int h : D.face_halfedges(cell))
      if (D.he_origin[h] == ca && D.dest(h) == cb) { found = h; break; }
    if (found < 0) fail("cell corner edge not found in the DCEL");
    const int ncell = D.he_face[D.twin(found)];
    const LatticeIsometry T = A.transition(cell, ncell);
    an = apply_raw(T, an, l);
    bn = apply_raw(T, bn, l);
    carried = T * carried;
    cell = ncell;
    t_num = bt_num;
    t_den = bt_den;
  }
  fail("segment trace exceeded " + std::to_string(cap) + " cell crossings");
}

const ResolvedFace& resolve_face(CellAtlas& A, int fi) {
  ResolvedFace& RS = A.resolved[fi];
  if (RS.ok) return RS;

  // Target: the face edge with the shallowest BFS chain back to an anchor.
  const tri_t& t = A.tface[fi];
  int e_target = -1;
  {
    int best_depth = INT32_MAX;
    for (int k = 0; k < 3; k++) {
      const int e = A.edge_id.at(CellAtlas::edge_key(t[k], t[(k + 1) % 3]));
      int d = 0, cur = e;
      while (A.bfs_parent_edge[cur] >= 0) { cur = A.bfs_parent_edge[cur]; d++; }
      if (d < best_depth) { best_depth = d; e_target = e; }
    }
  }
  std::vector<int> chain;   // root ... e_target
  for (int cur = e_target; cur >= 0; cur = A.bfs_parent_edge[cur]) chain.push_back(cur);
  std::reverse(chain.begin(), chain.end());

  const CellAtlas::AnchorEdge& AE = A.anchors[A.anchor_of_edge[chain[0]]];
  int  cell = AE.cell;
  int  eu = AE.u, ev = AE.v;
  Eisenstein pu = AE.pu, pv = AE.pv;
  EisensteinRational mid(pu + pv, 2);
  if (!inside_cell(A.cells[cell], mid.num, mid.den)) fail("anchor midpoint not inside its cell");

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

CellPoint locate_sample(CellAtlas& A, int fi, long a, long b, long c, long den) {
  const ResolvedFace& R = resolve_face(A, fi);
  const Eisenstein num = R.P[0] * narrow(a, "sample weight")
                       + R.P[1] * narrow(b, "sample weight")
                       + R.P[2] * narrow(c, "sample weight");
  const EisensteinRational q(num, den);
  CellTrace out = trace_segment(A, R.cell, R.anchor, q);
  return { out.cell, out.pos };
}

}  // namespace eisenstein_paint
