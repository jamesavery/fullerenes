// seam_step.cc -- fold Step 3 primitives: direction algebra, SeamAtlas,
// cone_holonomy.  See include/fullerenes/seam_step.hh and
// doc/ALGORITHM-seam-step.md.

#include "fullerenes/seam_step.hh"
#include "fullerenes/unfold.hh"   // Unfolding::transform_line (the glue isometry)
#include <sstream>


using std::ostringstream;

// The direction algebra (dot2, same_dir, primitive_dir, first_ccw, first_cw)
// and the half-open reflex-safe sector tests (Sector::contains_oc/co) live in
// eisenstein.hh alongside wedge and Sector.

// ---------------------------------------------------------------------------
// SeamAtlas.
// ---------------------------------------------------------------------------

SeamAtlas::SeamAtlas(const std::vector<std::pair<Eisenstein, node_t>>& outline) : O(outline) {
  const int m = (int)O.size();
  if (m < 3) throw std::runtime_error("SeamAtlas: outline has fewer than 3 corners");

  long area2 = 0;                                       // twice the signed shear area
  for (int i = 0; i < m; i++) area2 += wedge(pos(i), pos(i + 1));
  if (area2 == 0) throw std::runtime_error("SeamAtlas: outline has zero area");
  int_turn = area2 < 0 ? -1 : 1;

  for (int i = 0; i < m; i++) {
    const auto key = std::make_pair(lab(i), lab(i + 1));
    if (!mate_ix.emplace(key, i).second) {
      ostringstream msg;
      msg << "SeamAtlas: duplicate outline segment " << key.first << "->" << key.second
          << " (seam is not a tree)";
      throw std::runtime_error(msg.str());
    }
    if (!corner_ix.emplace(pos(i), i).second)
      throw std::runtime_error("SeamAtlas: repeated corner position (outline not simple)");
  }
  for (int i = 0; i < m; i++) mate(i);                  // every segment must have its mate
}

int SeamAtlas::mate(int i) const {
  i = wrap(i);
  const auto it = mate_ix.find({lab(i + 1), lab(i)});
  if (it == mate_ix.end()) {
    ostringstream msg;
    msg << "SeamAtlas: segment " << lab(i) << "->" << lab(i + 1) << " (edge " << i << ") has no mate";
    throw std::runtime_error(msg.str());
  }
  return it->second;
}

EisFrame SeamAtlas::glue(int i) const {
  i = wrap(i);
  const int j = mate(i);
  Eisenstein x0, x0p, w;
  Unfolding::transform_line({pos(i), pos(i + 1)}, {pos(j), pos(j + 1)}, x0, x0p, w);
  const EisFrame g{w, x0p - w * x0};                       // q |-> (q - x0)*w + x0p
  if (!g.is_isometry()) {
    ostringstream msg;
    msg << "SeamAtlas: glue of edge " << i << " is not a lattice isometry (rot norm2 "
        << g.rot.norm2() << ")";
    throw std::runtime_error(msg.str());
  }
  return g;
}

bool SeamAtlas::material(int m, Eisenstein u, bool ccw) const {
  const Eisenstein c = pos(m);
  const Eisenstein a = primitive_dir(pos(m - 1) - c);   // incoming-edge direction
  const Eisenstein b = primitive_dir(pos(m + 1) - c);   // outgoing-edge direction
  // CW outline (interior turn -1): material is the CCW sweep from incoming to
  // outgoing; mirrored for a CCW outline.  Ownership side per the traversal
  // convention (see the header): CCW walks (incoming, outgoing], CW walks
  // [incoming, outgoing).
  if (int_turn == -1) return ccw ? Sector{a, b}.contains_oc(u) : Sector{a, b}.contains_co(u);
  return ccw ? Sector{b, a}.contains_oc(u) : Sector{b, a}.contains_co(u);
}

SeamAtlas::Exit SeamAtlas::first_exit(Eisenstein G, Eisenstein Gn, long min_num, long min_den) const {
  const int ext = -int_turn;
  const int m = (int)O.size();
  const bool at_start = min_num < 0;

  // At the segment's true start, departure from an outline corner is governed
  // by the corner's CLOSED material wedge: a departure along either incident
  // edge is a legal on-boundary walk (its cells are shared boundary points),
  // so the gate ORs the two half-open ownership conventions -- ownership
  // matters for arc enumeration, not for chart existence.  Only a departure
  // strictly outside the wedge has no chart (kind == corner).  The incident
  // edges are non-transversal at the start; beyond it the parameter filter
  // handles everything.
  int skip1 = -1, skip2 = -1;
  if (at_start) {
    const auto at = corner_ix.find(G);
    if (at != corner_ix.end()) {
      if (!material(at->second, Gn - G, true) && !material(at->second, Gn - G, false))
        return {Exit::corner, at->second};
      skip1 = at->second; skip2 = wrap(at->second - 1);
    }
  }

  Exit best{Exit::none, -1};
  for (int i = 0; i < m; i++) {
    if (i == skip1 || i == skip2) continue;
    const Eisenstein s = pos(i), t = pos(i + 1);
    const int tG = Eisenstein::turn(s, t, G), tN = Eisenstein::turn(s, t, Gn);
    if (tN != ext || tG == ext) continue;               // must go from inside-or-on to strictly outside
    const int a = Eisenstein::turn(G, Gn, s), b = Eisenstein::turn(G, Gn, t);
    if (a == b) continue;                               // the step's line misses the edge's span
    const long num = std::labs(wedge(s - G, t - s)), den = std::labs(wedge(Gn - G, t - s));
    if (!at_start && num * min_den <= min_num * den) continue;    // at or before the last crossing
    if (a == 0 || b == 0) return {Exit::corner, i};     // passes exactly through an outline corner
    if (best.kind == Exit::none || num * best.den < best.num * den) best = {Exit::edge, i, num, den};
  }
  return best;
}

// ---------------------------------------------------------------------------
// cone_holonomy.
// ---------------------------------------------------------------------------

// The chart walk of one straight UNIT arc from an arbitrary development
// position P0 (a cone corner, an on-seam point, or an interior cell) in
// direction e: glue at each exit under the strictly-increasing crossing
// parameter, land, look up.  A unit arc contains no interior lattice point,
// hence no interior cone, so the walk stays within flat charts regardless of
// any cones elsewhere in a seam strip (shadow-immune -- unlike a row track).
ArcLanding resolve_arc(const SeamAtlas& atlas, const std::map<Eisenstein, node_t>& grid,
                       Eisenstein P0, Eisenstein e, const std::string& who) {
  const int MAX_GLUE = 16;
  EisFrame F;
  Eisenstein G0 = P0, G1 = P0 + e;
  long tn = -1, td = 1;                                 // last crossing parameter (start sentinel)
  for (int guard = 0;; guard++) {
    const auto ex = atlas.first_exit(G0, G1, tn, td);
    if (ex.kind == SeamAtlas::Exit::none) break;
    if (ex.kind == SeamAtlas::Exit::corner || guard >= MAX_GLUE) {
      ostringstream msg;
      msg << who << ": unit arc from (" << P0.first << "," << P0.second << ") dir ("
          << e.first << "," << e.second << ") "
          << (guard >= MAX_GLUE ? "exceeds the glue bound" : "exits through an outline corner")
          << " [at edge " << ex.index << ", segment (" << G0.first << "," << G0.second
          << ")->(" << G1.first << "," << G1.second << ")]";
      throw std::runtime_error(msg.str());
    }
    const EisFrame g = atlas.glue(ex.index);
    G0 = g(G0); G1 = g(G1); F = g.after(F);
    tn = ex.num; td = ex.den;                           // parameters survive the glue
  }
  const auto it = grid.find(G1);
  if (it == grid.end()) {
    ostringstream msg;
    msg << who << ": unit arc from (" << P0.first << "," << P0.second << ") lands undeveloped at ("
        << G1.first << "," << G1.second << ")";
    throw std::runtime_error(msg.str());
  }
  return {F, G1, it->second};
}

// Resolve one unit ring arc of the cone at corner m (see ConeRing::Entry).
static ConeRing::Entry resolve_unit_arc(const SeamAtlas& atlas, const std::map<Eisenstein, node_t>& grid,
                                        int m, Eisenstein e) {
  const ArcLanding a = resolve_arc(atlas, grid, atlas.pos(m), e,
                                   "ring arc of cone " + std::to_string(atlas.lab(m)));
  return {m, e, a.F, a.cell, a.id};
}

ConeRing cone_holonomy(const SeamAtlas& atlas, const std::map<Eisenstein, node_t>& grid,
                       int m0, Eisenstein anchor) {
  ConeRing out;
  if (!atlas.material(m0, anchor)) {
    ostringstream msg;
    msg << "cone_holonomy: anchor (" << anchor.first << "," << anchor.second
        << ") is not material at corner " << m0 << " (cone " << atlas.lab(m0) << ")";
    throw std::runtime_error(msg.str());
  }

  int m = m0;
  Eisenstein e = anchor;
  const int MAX_RING = 32, MAX_GLUE = 16;               // generous guards, never reached on valid input
  for (int cnt = 0; cnt < MAX_RING; cnt++) {
    out.ring.push_back(resolve_unit_arc(atlas, grid, m, e));

    Eisenstein w_acc(1, 0);
    e = e.nextCCW();
    for (int guard = 0; !atlas.material(m, e); guard++) {
      // Rotating CCW has exited through the OUTGOING edge at m.
      if (guard >= MAX_GLUE) {
        ostringstream msg;
        msg << "cone_holonomy: glue chain exceeds " << MAX_GLUE << " at cone " << atlas.lab(m0);
        throw std::runtime_error(msg.str());
      }
      out.rays.push_back({cnt, w_acc, primitive_dir(atlas.pos(m + 1) - atlas.pos(m))});
      const EisFrame g = atlas.glue(m);
      m = (atlas.mate(m) + 1) % atlas.n();              // the glue maps pos(m) to the mate edge's end
      e = g.rot * e;
      w_acc = g.rot * w_acc;
    }
    if (m == m0 && same_dir(e, anchor)) return out;     // closed onto the anchor
  }
  ostringstream msg;
  msg << "cone_holonomy: ring of cone " << atlas.lab(m0) << " did not close within "
      << MAX_RING << " entries";
  throw std::runtime_error(msg.str());
}

std::vector<ConeRing::Entry> ray_seeds(const SeamAtlas& atlas, const std::map<Eisenstein, node_t>& grid,
                                       int m0, Eisenstein ray, bool ccw, int count) {
  std::vector<ConeRing::Entry> out;
  const int n = atlas.n(), MAX_GLUE = 16;
  int m = m0;
  Eisenstein e = ccw ? first_ccw(ray) : first_cw(ray);
  for (int cnt = 0; cnt < count; cnt++) {
    // Glue until e is material at the current copy (under the traversal
    // convention's ownership side): CCW rotation exits through the OUTGOING
    // edge (landing at the mate's end), CW through the INCOMING edge (landing
    // at the mate's start).  This runs before the first emission too, and
    // steps through angle-only sliver copies.
    for (int guard = 0; !atlas.material(m, e, ccw); guard++) {
      if (guard >= MAX_GLUE) {
        ostringstream msg;
        msg << "ray_seeds: glue bound exceeded at cone " << atlas.lab(m0) << " (corner " << m0 << ")";
        throw std::runtime_error(msg.str());
      }
      if (ccw) {
        const EisFrame g = atlas.glue(m);
        m = (atlas.mate(m) + 1) % n;
        e = g.rot * e;
      } else {
        const int j = ((m - 1) % n + n) % n;
        const EisFrame g = atlas.glue(j);
        m = atlas.mate(j);
        e = g.rot * e;
      }
    }
    out.push_back(resolve_unit_arc(atlas, grid, m, e));
    e = ccw ? e.nextCCW() : e.nextCW();
  }
  return out;
}

std::array<node_t, 6> unit_ring(const SeamAtlas& atlas, const std::map<Eisenstein, node_t>& grid,
                                Eisenstein p) {
  std::array<node_t, 6> ring;
  for (int j = 0; j < 6; j++)
    ring[j] = resolve_arc(atlas, grid, p, Eisenstein::unit[j], "unit_ring").id;
  return ring;
}

