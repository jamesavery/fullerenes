#include "fullerenes/unfold.hh"
#include <array>
#include <numeric>
#include <sstream>
#include <stdexcept>

// Primitive (unit) direction of a lattice vector: v divided by gcd(|x|,|y|).
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
  if(s > 0) return wedge(a, d) >= 0 && wedge(d, b) >= 0;     // minor (<= 180 deg)
  if(s < 0) return !(wedge(b, d) > 0 && wedge(d, a) > 0);    // major: not strictly in the gap
  return true;                                               // a,b colinear (degenerate)
}

/*************************************************************/
/*************************************************************/
/*                 FOLDING IMPLEMENTATION                    */
/*************************************************************/
/*************************************************************/
//
// Restored pre-LLM polygon scan-conversion engine:
//   - connect_polygon : scan-convert the outline interior in each of the 3
//                       lattice edge-directions, drawing horizontal edges into
//                       oriented neighbour slots.
//   - connect_cross   : forward/backward raster of each seam segment; where the
//                       two rasters disagree (a split triangle), connect the
//                       edge across the seam via the gluing transform.
//   - identify_nodes  : unchanged (intact pre-LLM algorithm).
//
// Orientation note: neighbours[u][d] is u's neighbour in Eisenstein unit-
// direction d, so single-copy (deg-6) vertices come out CCW for free. Split
// vertices (the 12 cones; on-seam deg-6 with gcd>1) appear at several plane
// copies in different frames and are NOT yet correctly cyclically ordered --
// pre-LLM leaned on a post-hoc orient flag that no longer exists. The direct-
// oriented split-vertex assembly via the seam transforms is the next piece
// (see ../claude-projects/fold-unfold/DESIGN.md §6).


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
    // disagree we have a split triangle whose edge crosses the seam.
    for(size_t j=0;j<segment.size();j++){
      const Eisenstein& x(segment[j]), y(revsegment[j]);
      if(x != y){
        Eisenstein yp((y-xu0)*T+xu1); // y mapped through the gluing transform

        // Connect untransformed x-side node to transformed y-side node.
        auto itu = final_grid.find(x *omega_inv);
        auto itv = final_grid.find(yp*omega_inv);
        if(itu == final_grid.end() || itv == final_grid.end()) continue;

        node_t cu = itu->second, cv = itv->second;
        if(cu == cv) continue;

        neighbours[cu][i_omega]   = cv;
        neighbours[cv][i_omega+3] = cu;
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

      for(int x=x_start;x<x_end;x++){
        auto itu = final_grid.find(Eisenstein(x,  y+S.minY)*omega_inv);
        auto itv = final_grid.find(Eisenstein(x+1,y+S.minY)*omega_inv);
        if(itu == final_grid.end() || itv == final_grid.end()) continue;

        node_t u = itu->second, v = itv->second;
        neighbours[u][i_omega]   = v;
        neighbours[v][i_omega+3] = u;
      }
    }
  }
}

// One rotation: interior edges + edges across the seam.
void Folding::connect(int i_omega, vector<array<node_t,6>>& neighbours)
{
  if(!(debug_flags & DONT_CONNECT_POLYGON)) connect_polygon(i_omega, neighbours);
  if(!(debug_flags & DONT_CONNECT_ACROSS))  connect_cross  (i_omega, neighbours);
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

      if(segment.size() != revsegment.size()) continue;

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
  for(const auto& kv : by_angle) cycle.push_back(kv.second);   // CCW by unit-angle
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

  // Outline segment (srcLabel,tgtLabel) -> outline index, for the cone wedge chain.
  map<pair<node_t,node_t>,size_t> seg_index;
  for(size_t s = 0; s < outline.size(); s++)
    seg_index[{outline[s].second, outline[(s+1)%outline.size()].second}] = s;

  // Assemble every node that lies on the outline (all cones, and on-seam deg-6
  // -- i.e. cones OR multi-copy); only strictly-interior single-copy deg-6 use
  // the scan-conversion slots in CCW order.
  Graph nbr(N, GRAPH_DMAX);
  for(node_t u = 0; u < N; u++){
    if(degrees[u] != 6 || node_pos[u].size() > 1){
      for(node_t v : assemble_cone(u, seg_index)) nbr.push_back(u, v);
    } else {
      for(int d = 0; d < 6; d++)
        if(neighbours[u][d] != -1) nbr.push_back(u, neighbours[u][d]);
    }
  }

  Triangulation T(nbr);
  restore_original_labels(T);
  return T;
}

// Build the compacted-id -> original-label permutation and apply it to the folded
// graph and the Folding's id-keyed state. A tracked original (its class's canonical
// grid id < N_orig) is restored to its original id cone_perm.inverse()[canonical];
// every other vertex takes one of the remaining ids of [0,N). This is a bijection of
// [0,N) for ANY partial subset of originals present in the result, which is the usual
// case after a transform (the result keeps only some pre-transform vertices). A
// present original whose original id is out of the result's range (a shrinking
// transform) cannot keep its id and is treated as untracked.
void Folding::restore_original_labels(Triangulation& T)
{
  const node_t N = (node_t)node_pos.size();
  // cone_perm covers the tracked range (size >= N_orig): for a cone-star outline it
  // is the dense [0,Vmax] permutation, larger than N_orig=Ncones. Index its inverse
  // only at canonical < N_orig. Empty => no relabel known => identity restore.
  const bool have_perm = (node_t)cone_perm.size() >= N_orig && N_orig > 0;
  const Permutation cone_inv = have_perm ? cone_perm.inverse()
                                         : Permutation::identity(N_orig);

  Permutation relabel(N);
  vector<bool>   used(N, false);
  vector<node_t> untracked;                  // classes that take a complement id
  for(const auto& [canonical, c] : compact){
    const node_t orig = (canonical < N_orig) ? cone_inv[canonical] : node_t(-1);
    if(orig >= 0 && orig < N){ relabel[c] = orig; used[orig] = true; }  // present tracked original
    else untracked.push_back(c);
  }
  node_t next = 0;                           // fill the complement [0,N) \ used in order
  for(node_t c : untracked){
    while(used[next]) ++next;                // a free id always remains: |untracked| == #unused
    relabel[c] = next; used[next] = true;
  }

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
