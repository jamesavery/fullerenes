// fold.cc -- BYTE-FAITHFUL namespaced copy of the production library's Folding
// method definitions (and the file-static reduce_dir / in_wedge helpers), from
// src/c++/fold.cc.
//
// Baseline copy for the fold-unfold sub-project: NO algorithm changes. Only
// edits vs production: (a) wrapped in `namespace foldunfold`, (b) include path
// points at the local fullerenes/unfold.hh.
//
// Issues that claude tends to forget, leading to repeated bugs:
//
// 1. Vertices DO NOT have unique positions. Use arc positions instead (pair
//    of coords for src and target). Transform arcs instead of single points, etc.
//    
// 2. Don't use undirected edges, which throw away orientation information. Use arcs
//    (directed edges / half-edges), and directed seam segments.
//
// 3. Upper case U,V,W denote cone vertices (outline corners), of which there are at most 12.
//    U->V is a directed geodesic, not generally edge-aligned and not generally unit length.
//    U--V is an undirected geodesic.
//    Lower case u,v,w denote general graph vertices (can be on outline or interior).
//    u->v is an arc in the full graph, unit length. u--v is an (undirected) edge in the graph.
//
// Wiring notes:
//   - reduce_dir uses std::gcd from <numeric> (NOT the foldunfold::gcd(int,int) <- JA Note: Why? Or rather, why did you write foldunfold::gcd(int,int) if you don't use it?
//     in fu_geometry.cc) -- left as std::gcd, byte-identical to production.
//   - polygon::draw_line is called as polygon::draw_line(...) and resolves to
//     the local foldunfold::polygon static member.
//   - Eisenstein::unit, wedge, Unfolding::transform_line, Graph, GRAPH_DMAX
//     resolve to the foldunfold copies (Unfolding) or the GLOBAL library
//     symbols (Eisenstein::unit, wedge, Graph, GRAPH_DMAX).

// JA Note: There is a general issue of CW vs CCW confusion. We ultimately want neighbours in CW order as seen from the surface
//          (matching Buckygen's orientation). However, the most important thing is to make sure it is consistent. 

#include "fullerenes/unfold.hh"
#include "fullerenes/seam_step.hh"           // SeamAtlas, cone_holonomy (exact corner degrees)
#include <array>
#include <numeric>
#include <sstream>
#include <stdexcept>


// gcd-reduced direction of a lattice vector: v divided by gcd(|x|,|y|).
// NOT necessarily a unit direction -- e.g. reduce_dir({9,3})=={3,1}, norm2=13.
static Eisenstein reduce_dir(Eisenstein v){
  int g = std::gcd(std::abs(v.first), std::abs(v.second));
  return g ? Eisenstein(v.first/g, v.second/g) : v;
}

// Is direction d in the closed CCW angular sector from a to b? Unlike the
// Sector helper, this handles sectors > 180 deg (a cone's interior can be up to
// 300 deg): for a major sector (wedge(a,b) < 0) the interior is everything
// except the minor gap CCW from b to a.
static bool in_wedge(Eisenstein a, Eisenstein b, Eisenstein d){
  long s = wedge(a, b);
  if(s > 0) return  wedge(a, d) >= 0 && wedge(d, b) >= 0;    // minor (<= 180 deg)
  if(s < 0) return !(wedge(b, d) > 0 && wedge(d, a) > 0);    // major: not strictly in the gap
  return true;                                               // a,b colinear (degenerate)
}

/*************************************************************/
/*************************************************************/
/*                 FOLDING IMPLEMENTATION                    */
/*************************************************************/
/*************************************************************/
//
//   - connect_polygon : scan-convert the outline interior in each of the 3
//                       lattice edge-directions, drawing horizontal edges into
//                       oriented neighbour slots.
//   - connect_cross   : This is for connecting deg 6 vertices across seams.
//                       forward/backward raster of each seam segment; where the
//                       two rasters disagree (a split triangle), connect the
//                       edge across the seam via a gluing transform.
//                       //JA Note: There can be up to 5 split edges incident to a deg6 vertex (consider a thin long triangle). Scanning in the 3 lattice edge-directions is necessary here also.
//   - connect_cones   : Connect neigbours to cone vertices U (corners of outline). 
//                       Chain the wedges W(U)_i = V_i -> U -> V_{i+1}, gluing each outgoing
//                       seam segment of W_i to its matching seam segment of W_{i+1},
//                       transforming the arcs with source U to a common frame, from which U's neighbours can be read off in CW order. 
//
//   - identify_nodes  : unchanged (intact pre-LLM algorithm).
//
// Orientation note: neighbours[u][d] is u's neighbour in Eisenstein unit-
// direction d, so single-copy (deg-6) vertices come out CCW for free in connect_polygon. Split
// arcs must be transformed to a single frame for this to work. That is the job of connect_cross and
// connect_cones (currently named assemble_cone).
// JA Note: I edited the above to remove some confusion and align with my intent.

// Assemble every outline-corner label exactly via the glue-chained cone
// holonomy walk (seam_step.hh): corner_cycles[c] is the cone's plane-CCW
// neighbour cycle (fold Step 2), and degrees[c] its exact degree -- a
// corner's wedge can be cut fractionally by the seams, so no local count or
// bare-lookup binning is exact (the walk resolves each ring arc through its
// glue chain). Overwrites the provisional ctor-3 degrees; for graph-backed
// unfoldings the result must agree with the graph. `degrees` is indexed by
// compacted id; outline labels are tracked originals, preserved by the
// compaction. Replaces the retired assemble_cone (unit-cell angle binning,
// which mis-collected pinched and sliver wedges, e.g. on star unfoldings).
void Folding::assemble_corners()
{
  const SeamAtlas atlas(outline);
  std::set<node_t> labels;
  for(const auto& o : outline) labels.insert(o.second);
  for(node_t c : labels){
    int m0 = -1; Eisenstein anchor;
    for(int m = 0; m < atlas.n() && m0 < 0; m++){
      if(atlas.lab(m) != c) continue;
      for(int e = 0; e < 6; e++)
        if(atlas.material(m, Eisenstein::unit[e])){ m0 = m; anchor = Eisenstein::unit[e]; break; }
    }
    if(m0 < 0)
      throw std::runtime_error("Folding: cone " + std::to_string(c) + " has no material unit at any copy");
    const ConeRing ring = cone_holonomy(atlas, final_grid, m0, anchor);
    if(c >= (node_t)degrees.size())
      throw std::logic_error("Folding: outline label " + std::to_string(c) + " outside degrees range");
    degrees[c] = (int)ring.ring.size();
    auto& cyc = corner_cycles[c];
    cyc.clear();
    for(const auto& r : ring.ring) cyc.push_back(r.id);
    corner_rings[c] = ring;
  }
}

// connect_cross connects edges in split triangles
// (for outline segments that do not align with the Eisenstein grid).
void Folding::connect_cross(int i_omega, vector<array<node_t,6>>& neighbours)
{
  map<arc_t,arccoord_t> reverse_arc;
  Eisenstein xu, xv;
  node_t u,v;

  Eisenstein omega     = Eisenstein::unit[i_omega];
  Eisenstein omega_inv = Eisenstein::unit[6-i_omega];

  // Register reverse arcs (rotated coordinates of arc v->u, keyed by {v,u}).
  for(size_t i=0;i<outline.size();i++){
    tie(xu,u) = outline[i];
    tie(xv,v) = outline[(i+1) % outline.size()];

    reverse_arc[{v,u}] = {xu*omega, xv*omega};
  }

  for(size_t i=0;i<outline.size();i++){
    tie(xu,u) = outline[i];
    tie(xv,v) = outline[(i+1) % outline.size()];

    // Coordinates of arc u->v and its glued mate v->u.
    arccoord_t Xuv = {xu*omega,xv*omega}, Xvu = reverse_arc[{u,v}];

    // Affine transform taking the line segment Xuv into Xvu.
    Eisenstein xu0, xu1, T;
    Unfolding::transform_line(Xuv,Xvu, xu0,xu1, T);

    // Rasterize the u->v segment forwards and backwards.
    vector<Eisenstein> segment   (polygon::draw_line(xu*omega,xv*omega)),
                       revsegment(polygon::draw_line(xv*omega,xu*omega));
    reverse(revsegment.begin(),revsegment.end());
    if(segment.size() != revsegment.size()) continue;

    // Forward rasterization rounds one way, backwards the other; where they
    // agree (x==y), we have an on-seam deg-6 vertex (one copy on this U->V lip, the
    // matching copy on the glued V->U lip); where they disagree, a split triangle whose
    // edge crosses the seam.
    for(size_t j=0;j<segment.size();j++){
      const Eisenstein& x(segment[j]), y(revsegment[j]);
      if(x == y){
        // On-seam deg-6 vertex. Its complete oriented ring is assembled by the
        // SEAM STEP (seam_rings, fold Step 3) and installed wholesale after the
        // scans -- the per-scan placement that used to live here mixed the two
        // copies' frames and could not resolve arcs pinched into thin wedges.
        // Identification of the two copies remains identify_nodes' job.
        continue;
      } else {
        Eisenstein yp((y-xu0)*T+xu1); // y mapped through the gluing transform

        // Connect untransformed x-side node to transformed y-side node.
        auto itu = final_grid.find(x *omega_inv);
        auto itv = final_grid.find(yp*omega_inv);
        if(itu == final_grid.end() || itv == final_grid.end()) continue;

        node_t cu = itu->second, cv = itv->second;
        if(cu == cv) continue;

        // Connect ONLY this segment's near side (cu), at the slot of the cross-edge's
        // own direction dcu (cu's frame). cv gets its reciprocal cv->cu when the MATCHING
        // V->U segment is scanned, as the near side in ITS frame -- so we never write cv
        // in this foreign U->V frame (that corrupted single-copy cyclic order). The slot
        // is per-direction, NOT i_omega: a deg-6 vertex can have several cross edges, and
        // forcing them all into i_omega would overwrite each other.
        Eisenstein dcu = (y - x) * omega_inv;          // cross-edge direction at cu
        if(dcu.norm2() != 1) continue;
        neighbours[cu][(6 - dcu.unit_angle()) % 6] = cv;
      }
    }
  }
}

// connect_polygon connects all interior edges of the outline polygon by exact
// scan-conversion in direction i_omega.
void Folding::connect_polygon(int i_omega, vector<array<node_t,6>>& neighbours)
{
  Eisenstein omega     = Eisenstein::unit[i_omega];
  Eisenstein omega_inv = Eisenstein::unit[6-i_omega];

  vector<Eisenstein> outline_coords;
  outline_coords.reserve(outline.size());
  for(const auto& o: outline) outline_coords.push_back(o.first*omega);

  polygon Prot(outline_coords);
  polygon::scanline S(Prot.scanConvert());

  for(size_t y=0;y<S.xs.size();y++){            // For each scanline..
    for(size_t j=0;j<S.xs[y].size()/2;j++){     // ..go through each inner interval
      int x_start = S.xs[y][2*j], x_end = S.xs[y][2*j+1];
      for(int x=x_start;x<x_end;x++)
        connect_arc(neighbours, Eisenstein(x,y+S.minY), Eisenstein(x+1,y+S.minY), omega_inv, i_omega);
    }
  }
}

// Write arc p0--p1 (rotated back by omega_inv) into both endpoints' oriented
// slots. Idempotent: a shared edge written twice with the same value is a no-op.
void Folding::connect_arc(vector<array<node_t,6>>& nb, Eisenstein p0, Eisenstein p1,
                          Eisenstein omega_inv, int i_omega) const
{
  auto a = final_grid.find(p0*omega_inv), b = final_grid.find(p1*omega_inv);
  if(a==final_grid.end() || b==final_grid.end()) return; // <- JA Note: This looks like a silent failure, corrupting the result!
  if(a->second == b->second) return;                     // <- JA Note: This also looks like a silent failure, corrupting the result.
  nb[a->second][i_omega]   = b->second;
  nb[b->second][i_omega+3] = a->second;
}

// One rotation: interior (in-sheet) edges only.  Seam-crossing arcs are the
// seam step's (SeamStep::cross_arcs) and the cone holonomies' -- the raster
// connect_cross is RETIRED from the pipeline (its index-paired draw_line
// detection is invalid on slanted seams; kept compiled for old diagnostics).
void Folding::connect(int i_omega, vector<array<node_t,6>>& neighbours)
{
  if(!(debug_flags & DONT_CONNECT_POLYGON))
    connect_polygon(i_omega, neighbours);
}

// The whole outline is connected into a triangulation by drawing the horizontal
// edges in each of the 3 lattice edge-directions (rotations 0, 1, 2).
void Folding::connect(vector<array<node_t,6>>& neighbours)
{
  connect(0,neighbours);
  connect(1,neighbours);
  connect(2,neighbours);
}

// identify_nodes takes a polygon outline and a map from the Eisenstein grid
// to node ids where some of the (triangulation) graph nodes on the outline are
// split into two or more grid points: to fold up the polygon into a polyhedron,
// we must identify which ones are the same.
//
// HOW:
//     0. For each arc U->V, register its grid position as the reverse to V->U
//
//     A graph G of the node equivalence classes is next built up as follows:
//     1. For each arc U->V: draw U->V forwards and backwards to detect split triangles
//     2. For each node u on the U->V path that is *not* on a split, it must be repeated
//        as a node v on the V->U path, so add edge u-v to the equivalence graph G.
//        (This is an on-seam deg 6 vertex)
//
//     Finally, we produce a canonical representation of each equivalence class:
//     3. Compute connected components of G
//     4. For each connected component, choose smallest id as canonical represention
//
// RESULT:
//     A vector same_as such that for each original node id u, same_as[u] is the smallest
//     node id in its equivalence class.
vector<int> Folding::identify_nodes(const IDCounter<Eisenstein>& grid, const vector< pair<Eisenstein,node_t>>& outline) const
{
  node_t u, v, U, V;
  Eisenstein xu, xv, XU, XV;
  set<edge_t> same_as;

  for(int i_omega=0;i_omega<3;i_omega++){
    map<arc_t,arccoord_t> reverse_arc;

    Eisenstein omega     = Eisenstein::unit[i_omega],
               omega_inv = Eisenstein::unit[6-i_omega];

    if(omega * omega_inv != Eisenstein{1, 0})
      throw std::logic_error("Folding::identify_nodes: omega * omega_inv != 1 -- Eisenstein::unit table corrupt");

    // Register reverse arcs
    for(int i=0;i<outline.size();i++){
      tie(XU,U) = outline[i];
      tie(XV,V) = outline[(i+1) % outline.size()];

      reverse_arc[{V,U}] = {XU*omega,XV*omega};
    }

    // For each segment U->V of the outline, find the reverse one, and identify
    // the nodes on the path that are *not* on split triangles
    for(int i=0;i<outline.size();i++){
      tie(XU,U) = outline[i];
      tie(XV,V) = outline[(i+1) % outline.size()];

      arccoord_t XUV(XU*omega,XV*omega), XVU(reverse_arc[{U,V}]);

      Eisenstein x0,x0p,T;
      Unfolding::transform_line(XUV,XVU, x0,x0p, T);

      vector<Eisenstein>
	segment   (polygon::draw_line(XU*omega, XV*omega)),
	revsegment(polygon::draw_line(XV*omega, XU*omega));

      reverse(revsegment.begin(),revsegment.end());

      if(segment.size() != revsegment.size())
	throw std::logic_error("Folding::identify_nodes: fwd/bwd raster lengths differ on segment "
	  + std::to_string(U) + "->" + std::to_string(V));

      for(int j=0;j<segment.size();j++){
	const Eisenstein& x(segment[j]), y(revsegment[j]);
	if(x == y){
	  Eisenstein xp = (x-x0)*T+x0p;

	  node_t u = grid(x*omega_inv), v = grid(xp*omega_inv);
	  if((u>=0 && v>= 0)  && (u != v))
	    same_as.insert(edge_t{u,v});
	}
      }
    }
  }

  // Find connected components
  // Use grid.nextid (number of unique IDs) not grid.size() (number of map entries,
  // which includes multiple Eisenstein positions mapping to the same node ID).
  vector<int> same(grid.nextid);
  for(size_t i=0;i<grid.nextid;i++) same[i] = i;

  // Build graph from edge set for connected_components
  Graph same_nb(grid.nextid, GRAPH_DMAX);
  for(const edge_t& e: same_as){
    same_nb.push_back(e.first, e.second);
    same_nb.push_back(e.second, e.first);
  }
  Graph S(same_nb);
  vector<vector<node_t> > components(S.connected_components());

  for(auto& c: components){
    node_t canonical = *min_element(c.begin(),c.end());
    for(auto t: c) same[t] = canonical;
  }

  return same;
}


// Exact-arithmetic alternative to identify_nodes, for A/B validation: the
// on-seam lattice points of segment i are exactly pos(i) + s*(Delta/g) for
// s in [1, g); each is one flat vertex split into two copies, the mate being
// its glue image on the label-reversed segment. No rasterization involved.
// Must produce the SAME partition as identify_nodes wherever both apply
// (checked by the test suite); kept separate so the production raster path
// stays untouched.
vector<int> Folding::identify_nodes_exact(const IDCounter<Eisenstein>& grid,
                                          const vector< pair<Eisenstein,node_t>>& outline) const
{
  const SeamAtlas atlas(outline);
  set<edge_t> same_as;

  for(int i = 0; i < atlas.n(); i++){
    const Eisenstein D = atlas.pos(i+1) - atlas.pos(i);
    const int g = std::gcd(std::abs(D.first), std::abs(D.second));
    const Eisenstein d(D.first/g, D.second/g);
    const EisFrame glue = atlas.glue(i);
    for(int s = 1; s < g; s++){
      const Eisenstein p = atlas.pos(i) + d*s;
      const node_t u = grid(p), v = grid(glue(p));
      if(u >= 0 && v >= 0 && u != v) same_as.insert(edge_t{u,v});
    }
  }

  // Canonicalize: connected components of the identification graph, smallest
  // id as representative (same finalization as identify_nodes).
  vector<int> same(grid.nextid);
  for(size_t i=0;i<grid.nextid;i++) same[i] = i;
  Graph same_nb(grid.nextid, GRAPH_DMAX);
  for(const edge_t& e: same_as){
    same_nb.push_back(e.first, e.second);
    same_nb.push_back(e.second, e.first);
  }
  Graph S(same_nb);
  for(auto& c: S.connected_components()){
    node_t canonical = *min_element(c.begin(),c.end());
    for(auto t: c) same[t] = canonical;
  }
  return same;
}

// Assemble cone `cid`'s neighbour cycle directly CCW-ordered.
//
// A cone appears at one plane copy per wedge V_i -> U -> V_{i+1} along the
// outline. We chain the wedges (each outgoing outline seam glues to its
// label-reversed mate, which names the adjacent copy), bring every wedge into
// one common frame at U=origin by composing transform_line's rotation across
// each crossed seam, read each wedge's unit-distance neighbours from final_grid
// restricted to that wedge's angular sector, and bin them by their common-frame
// unit-direction. The angle bins are the CCW neighbour cycle; the cone's
// deficit shows up as the empty direction. final_grid (not arc_coords) is read,
// so this also picks up the Goldberg-Coxeter subdivision vertices.
// JA Note: This looks more complicated than I think it needs to be - and it looks like it transforms the entire wedge?
//          Also: it should only be used to connect the cones, not the on-seam deg 6 vertices.
// RETIRED from fold(): superseded by assemble_corners() (the exact holonomy
// walk). The unit-cell angle binning below mis-collects pinched arcs and
// sliver wedges (bare final_grid lookups read the wrong sheet or nothing).
// Kept only because older diagnostic tools still link it; do not add callers.
vector<node_t> Folding::assemble_cone(node_t cid, const map<pair<node_t,node_t>,size_t>& seg_index) const
{
  const int n   = (int)outline.size();
  const int deg = degrees[cid];

  // First outline appearance (copy) of this cone -- start anywhere; the deficit
  // gap falls out of the angle binning.
  int i0 = -1;
  for(int i = 0; i < n; i++){
    auto it = final_grid.find(outline[i].first);
    if(it != final_grid.end() && it->second == cid){ i0 = i; break; }
  }
  if(i0 < 0)
    throw std::logic_error("assemble_cone: cone " + std::to_string(cid) + " not on the outline");

  map<int,node_t> by_angle;          // common-frame unit-angle (0..5) -> neighbour id
  Eisenstein rot(1,0);               // this wedge's frame -> common frame (unit Eisenstein)
  std::set<int> visited;
  int i = i0;
  const size_t ncopies = node_pos[cid].size();

  for(size_t step = 0; step < ncopies; step++){
    if(visited.count(i)) break;      // chain closed
    visited.insert(i);

    Eisenstein p = outline[i].first;
    Eisenstein a = outline[(i-1+n)%n].first - p;   // ray toward prev (incoming seam)
    Eisenstein b = outline[(i+1)%n].first - p;     // ray toward next (outgoing seam)
    Eisenstein incoming = reduce_dir(a);           // unit dir along the incoming seam

    for(int k = 0; k < 6; k++){
      Eisenstein d = Eisenstein::unit[k];
      // Half-open: the incoming-seam neighbour belongs to the previous wedge, so
      // each shared seam neighbour is binned once -- and the cone's opening cut
      // (first wedge's incoming == last wedge's outgoing, 60 deg apart = the
      // deficit) is not double-counted.
      if(d == incoming) continue;
      if(!in_wedge(a, b, d)) continue;             // wedge interior (handles > 180 deg)
      auto it = final_grid.find(p + d);
      if(it == final_grid.end()) continue;
      if(it->second == cid) continue;              // skip an adjacent copy of U itself (not a neighbour)
      Eisenstein cd = d * rot;                     // direction in the common frame
      by_angle[cd.unit_angle()] = it->second;
    }

    // Advance across the outgoing seam to the adjacent copy, composing its glue.
    node_t U = outline[i].second, V = outline[(i+1)%n].second;
    auto mit = seg_index.find({V,U});              // label-reversed mate segment
    if(mit == seg_index.end()) break;              // open (deficit) seam -- end of chain
    size_t j = mit->second;
    arccoord_t l1 = { outline[i].first,     outline[(i+1)%n].first };
    arccoord_t l2 = { outline[j].first,     outline[(j+1)%n].first };
    Eisenstein x0, x0p, w;
    Unfolding::transform_line(l1, l2, x0, x0p, w);
    rot = rot * w.complex_conj();                  // next wedge's frame -> common frame
    i = (int)((j+1)%n);
  }

  vector<node_t> cycle;
  for(const auto& kv : by_angle) cycle.push_back(kv.second);   // surface-CW (conjugated common frame)
  if((int)cycle.size() != deg)
    throw std::logic_error("assemble_cone: cone " + std::to_string(cid) + " collected "
                           + std::to_string(cycle.size()) + " != degree " + std::to_string(deg));
  return cycle;
}

Triangulation Folding::fold()
{
  node_t N = node_pos.size();

  // Base scan-conversion: oriented slots neighbours[u][d] = u's neighbour in
  // unit-direction d. Correct for single-copy nodes (and, for now, on-seam deg-6).
  vector<array<node_t,6>> neighbours(N);
  for(auto& a: neighbours) a.fill(-1);
  connect(neighbours);

  // Complete the boundary-adjacent flat vertices.  connect_polygon has
  // written every in-sheet arc, so a flat vertex is incomplete exactly when
  // part of its hexagon is developed across a seam; each such arc is one
  // unit-segment chart walk (resolve_arc -- exact for any outline, immune to
  // cones near the seams).  A multi-copy (on-seam) vertex takes its WHOLE row
  // from one copy, so the row is coherent in that copy's frame; a single-copy
  // vertex's frame is the global one, so it fills just its empty slots.
  // Slot s of the scan convention holds direction unit[(6-s)%6].
  {
    const SeamAtlas atlas(outline);
    for(node_t u = 0; u < N; u++){
      if(degrees[u] != 6) continue;
      if(node_pos[u].size() > 1){
        const auto ring = unit_ring(atlas, final_grid, node_pos[u][0]);
        for(int j = 0; j < 6; j++) neighbours[u][(6 - j) % 6] = ring[j];
      } else {
        for(int s = 0; s < 6; s++)
          if(neighbours[u][s] == -1)
            neighbours[u][s] = resolve_arc(atlas, final_grid, node_pos[u][0],
                                           unit_direction(6 - s), "fold completion").id;
      }
    }
  }

  // Build each vertex's oriented neighbour cycle and WRITE it to oriented
  // positions in the graph (neighbours[u*dmax + i]), never push_back.
  //
  // Two sources, one per vertex class (all plane-CCW): cones (degree != 6)
  // from corner_cycles (the ctor's holonomy walks); flat vertices from the
  // slots (in-sheet scan + the completion above).  The slots run CW by index
  // (slot s = direction unit[(6-s)%6]), so they are emitted reversed to give
  // the plane-CCW cycle.
  Graph nbr(N, GRAPH_DMAX);
  for(node_t u = 0; u < N; u++){
    vector<node_t> cw;
    if(degrees[u] != 6){
      cw = corner_cycles.at(u);
    } else {
      for(int d = 5; d >= 0; d--)
        if(neighbours[u][d] != -1) cw.push_back(neighbours[u][d]);
    }
    if((int)cw.size() != degrees[u])
      throw std::logic_error("fold: vertex " + std::to_string(u) + " assembled "
        + std::to_string(cw.size()) + " != degree " + std::to_string(degrees[u]));
    for(size_t i = 0; i < cw.size(); i++)
      nbr.neighbours[(size_t)u*nbr.dmax + i] = cw[i];                  // write to oriented position i
    nbr.deg[u] = (uint8_t)cw.size();
  }

  // 100%-detection: every arc was computed independently from both endpoints
  // (in-sheet scan, chart-walk completion, holonomy) -- the results must be
  // symmetric.
  for(node_t u = 0; u < N; u++)
    for(int i = 0; i < nbr.deg[u]; i++){
      const node_t v = nbr.neighbours[(size_t)u*nbr.dmax + i];
      bool sym = false;
      for(int k = 0; k < nbr.deg[v] && !sym; k++)
        sym = nbr.neighbours[(size_t)v*nbr.dmax + k] == u;
      if(!sym)
        throw std::logic_error("fold: asymmetric arc " + std::to_string(u) + "->" + std::to_string(v));
    }

  Triangulation T(nbr);
  restore_original_labels(T);
  return T;
}

// Build the compacted-id -> original-label permutation and apply it to the folded
// graph and the Folding's id-keyed state. A tracked original (its class's canonical
// grid id < N_orig) is restored to cone_perm.inverse()[canonical]; every other vertex
// takes one of the remaining ids of [0,N). This is a bijection of [0,N) for ANY partial
// subset of originals present in the result (the usual case after a transform). Empty
// cone_perm (outline-only ctors) => identity restore over the tracked range.
void Folding::restore_original_labels(Triangulation& T)
{
  const node_t N = (node_t)node_pos.size();
  const bool have_perm = (node_t)cone_perm.size() >= N_orig && N_orig > 0;
  const Permutation cone_inv = have_perm ? cone_perm.inverse()
                                         : Permutation::identity(N_orig);

  // Pin each present tracked original to its original id (cone_perm.inverse()); the rest
  // (untracked / fold-created) take the remaining ids. Same retain-then-fill as the
  // Folding compaction (fill_free_ids), so there is one implementation of the pattern.
  vector<node_t> pin(N, -1);
  for(const auto& [canonical, c] : compact)
    if(canonical < N_orig) pin[c] = cone_inv[canonical];
  const Permutation relabel(fill_free_ids(pin));

  // One isomorphism across every id-keyed view so the whole Folding shares the
  // original-label id space (fold() is single-shot).
  T.apply_permutation(relabel);
  for(auto& kv : final_grid) kv.second = relabel[kv.second];
  vector<vector<Eisenstein>> np(N);
  vector<int>                dg(N);
  for(node_t c = 0; c < N; c++){ np[relabel[c]] = std::move(node_pos[c]); dg[relabel[c]] = degrees[c]; }
  node_pos = std::move(np);
  degrees  = std::move(dg);
}

vector<node_t> Folding::outline_nodes() const
{
  node_t u, outline_N=0;
  Eisenstein xu;

  vector<node_t> same_nodes(identify_nodes(grid,outline));
  vector<node_t> outline_newnames(outline.size());

  for(int i=0;i<outline.size();i++){
    tie(xu,u) = outline[i];
    outline_newnames[i] = grid(xu);
    outline_N = max(outline_N,u+1);
  }

  vector<node_t> new_nodenames(outline_N,-1);
  for(int i=0;i<outline.size();i++){
    int u = outline[i].second;
    int stored_up = new_nodenames[u];

    if(stored_up != -1 && same_nodes[stored_up] != same_nodes[outline_newnames[i]]){
      std::ostringstream msg;
      msg << "Folding::outline_nodes: same-node disagreement at outline["
          << i << "] = {{" << outline[i].first.first << ","
          << outline[i].first.second << "}," << outline[i].second
          << "} -> u = " << u << ". stored_up = " << stored_up
          << ", outline_newnames[" << i << "] = " << outline_newnames[i];
      throw std::logic_error(msg.str());
    }
    new_nodenames[u] = same_nodes[outline_newnames[i]];
  }

  return new_nodenames;
}


string Folding::to_latex(int K, int L, int label_vertices, bool draw_equilaterally, bool include_headers) const
{
  string s;
  return s;
}

