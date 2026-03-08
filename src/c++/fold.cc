#include "fullerenes/unfold.hh"
#include <array>
#include <queue>

/*************************************************************/
/*************************************************************/
/*                 FOLDING IMPLEMENTATION                   */
/*************************************************************/
/*************************************************************/



// connect_cross connects edges in split triangles
// (for outline segments that do not align with Eisenstein grid)
void Folding::connect_cross(int i_omega, set<edge_t> &edges, cross_info_t& cross_info)
{
  map<arc_t,arccoord_t> reverse_arc;
  Eisenstein xu, xv;
  node_t u,v;

  Eisenstein omega     = Eisenstein::unit[i_omega];
  Eisenstein omega_inv = Eisenstein::unit[6-i_omega];

  // Register reverse arcs
  for(int i=0;i<outline.size();i++){
    tie(xu,u) = outline[i];
    tie(xv,v) = outline[(i+1) % outline.size()];

    reverse_arc[{v,u}] = {xu*omega, xv*omega}; // Rotated coordinates of arc v->u
  }

  for(int i=0;i<outline.size();i++){
    tie(xu,u) = outline[i];
    tie(xv,v) = outline[(i+1) % outline.size()];

    // First get the coordinates of arc u->v and the reverse v->u
    arccoord_t Xuv = {xu*omega,xv*omega}, Xvu = reverse_arc[{u,v}];

    // What the affine transform that takes the line segment Xuv into Xvu?
    Eisenstein
      xu0,			// Coord of u in u->v
      xu1,			// Coord of u in v->u
      T;			// Rotation transform
    Unfolding::transform_line(Xuv,Xvu, xu0,xu1, T);

    // Alongside u->v, rasterize u->v line segment forwards and backwards
    vector<Eisenstein> segment   (polygon::draw_line(xu*omega,xv*omega)),
                       revsegment(polygon::draw_line(xv*omega,xu*omega));
    reverse(revsegment.begin(),revsegment.end());

    if(segment.size() != revsegment.size()) continue;

    // Go through the nodes of the line segments rasterized back and forth along u->v
    for(int j=0;j<segment.size();j++){
      const Eisenstein& x(segment[j]), y(revsegment[j]);
      // Forward rasterization rounds to the right, backwards to the left.
      // So when x != y, we have a split triangle and an edge that needs to be connected across the boundary
      if(x != y){
	Eisenstein xp((x-xu0)*T+xu1); // Rotate and translate from u->v to v->u coords
	Eisenstein yp((y-xu0)*T+xu1); // Rotate and translate from u->v to v->u coords

	// Connect untransformed u to transformed v
	Eisenstein pu = x*omega_inv, pv = yp*omega_inv;
	auto itu = final_grid.find(pu), itv = final_grid.find(pv);
	if(itu == final_grid.end() || itv == final_grid.end()) continue;

	node_t u = itu->second, v = itv->second;
	if(u == v) continue;
	if(debug_flags & WRITE_FILE) printf("Connect cross arc %d to %d \n",u,v);

	// Direction from u to v at pu: in omega-rotated frame, edge goes x→y.
	// Unrotate to get the absolute Eisenstein direction at pu.
	Eisenstein d_u_to_v = (y - x) * omega_inv;
	// Direction from v to u at pv: in omega-rotated frame on v's side,
	// u is at xp, v at yp. Direction yp→xp = (x-y)*T. Unrotate.
	Eisenstein d_v_to_u = (x - y) * T * omega_inv;

	int dir_u = d_u_to_v.unit_angle();
	int dir_v = d_v_to_u.unit_angle();

	edge_t e(u, v);
	// edge_t always has first < second. Map directions accordingly.
	if(cross_info.find(e) == cross_info.end()){
	  if(u <= v)
	    cross_info[e] = {dir_u, dir_v, pu, pv};
	  else
	    cross_info[e] = {dir_v, dir_u, pv, pu};
	}

	edges.insert(e);
      }
    }
  }
}

// connect_interior finds all edges where both endpoints are in the grid
// by directly checking the 3 positive Eisenstein neighbor directions.
void Folding::connect_interior(set<edge_t> &edges)
{
  static const Eisenstein dirs[3] = {{1,0}, {0,1}, {-1,1}};

  for(const auto& kv: final_grid){
    const Eisenstein& p = kv.first;
    node_t u = kv.second;

    for(int d = 0; d < 3; d++){
      auto it = final_grid.find(p + dirs[d]);
      if(it != final_grid.end()){
        node_t v = it->second;
        if(u != v) edges.insert(edge_t(u,v));
      }
    }
  }
}

// The whole outline is connected into a triangulation / cubic graph dual
// by finding interior edges via direct neighbor lookup and cross-boundary
// edges via line rasterization in each of the 3 Eisenstein directions.
void Folding::connect(set<edge_t> &edges, cross_info_t& cross_info)
{
  if(!(debug_flags & DONT_CONNECT_POLYGON))  connect_interior(edges);
  if(!(debug_flags & DONT_CONNECT_ACROSS)){
    connect_cross(0, edges, cross_info);
    connect_cross(1, edges, cross_info);
    connect_cross(2, edges, cross_info);
  }
}

// identify_nodes takes a polygon outline and a map from the eisenstein grid
// to node id's where some of the (triangulation) graph nodes on the outline are
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

    assert((omega * omega_inv == Eisenstein{1,0}));

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
  neighbours_t same_nb(grid.nextid, GRAPH_DMAX);
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


Triangulation Folding::fold()
{
  node_t N = node_pos.size();
  set<edge_t> edge_set;
  cross_info_t cross_info;
  connect(edge_set, cross_info);

  // Step 1: Build per-position direction map.
  // pos_nb[p][d] = the compacted node at grid position p + unit[d], or -1.
  map<Eisenstein, array<node_t, 6>> pos_nb;
  for(const auto& [p, u] : final_grid){
    array<node_t, 6> dirs;
    dirs.fill(-1);
    for(int d = 0; d < 6; d++){
      auto it = final_grid.find(p + Eisenstein::unit[d]);
      if(it != final_grid.end() && it->second != u)
        dirs[d] = it->second;
    }
    pos_nb[p] = dirs;
  }

  // Step 2: Compute rotation from each grid position to its node's canonical frame.
  // The rotation comes from the boundary identification transform T.
  // Re-run the identification loop (same as identify_nodes) to extract rotations.
  map<Eisenstein, int> pos_rotation;
  for(const auto& [p, u] : final_grid)
    pos_rotation[p] = 0;

  // Build adjacency list for position identifications: pos -> [(other_pos, T_angle)]
  // where T_angle is the rotation from pos's frame to other_pos's frame.
  map<Eisenstein, vector<pair<Eisenstein, int>>> id_adj;

  for(int i_omega = 0; i_omega < 3; i_omega++){
    Eisenstein omega     = Eisenstein::unit[i_omega];
    Eisenstein omega_inv = Eisenstein::unit[6-i_omega];

    map<arc_t, arccoord_t> reverse_arc;
    for(size_t i = 0; i < outline.size(); i++){
      auto [xu, u] = outline[i];
      auto [xv, v] = outline[(i+1) % outline.size()];
      reverse_arc[{v,u}] = {xu*omega, xv*omega};
    }

    for(size_t i = 0; i < outline.size(); i++){
      auto [xu, u] = outline[i];
      auto [xv, v] = outline[(i+1) % outline.size()];

      arccoord_t XUV(xu*omega, xv*omega), XVU(reverse_arc[{u,v}]);
      Eisenstein x0, x0p, T;
      Unfolding::transform_line(XUV, XVU, x0, x0p, T);

      if(!T.isUnit()) continue;
      int T_angle = T.unit_angle();

      auto segment    = polygon::draw_line(xu*omega, xv*omega);
      auto revsegment = polygon::draw_line(xv*omega, xu*omega);
      reverse(revsegment.begin(), revsegment.end());

      if(segment.size() != revsegment.size()) continue;

      for(size_t j = 0; j < segment.size(); j++){
        if(segment[j] == revsegment[j]){
          Eisenstein xp = (segment[j] - x0) * T + x0p;
          Eisenstein p1 = segment[j] * omega_inv;
          Eisenstein p2 = xp * omega_inv;

          if(p1 != p2 && final_grid.count(p1) && final_grid.count(p2) &&
             final_grid.at(p1) == final_grid.at(p2)){
            // p1 and p2 are identified copies of the same node.
            // Direction d at p1 corresponds to direction (d + T_angle) % 6 at p2.
            id_adj[p1].push_back({p2, T_angle});
            id_adj[p2].push_back({p1, (6 - T_angle) % 6});
          }
        }
      }
    }
  }

  // BFS from each canonical position to propagate rotations.
  for(node_t u = 0; u < N; u++){
    if(node_pos[u].size() <= 1) continue;

    const Eisenstein& p0 = node_pos[u][0];
    pos_rotation[p0] = 0;

    queue<Eisenstein> bfs;
    bfs.push(p0);
    set<Eisenstein> visited;
    visited.insert(p0);

    while(!bfs.empty()){
      Eisenstein p = bfs.front(); bfs.pop();
      int rot_p = pos_rotation[p];

      for(auto& [q, T_angle] : id_adj[p]){
        if(visited.count(q)) continue;
        if(final_grid.at(q) != u) continue;

        // Direction d at p corresponds to direction (d + T_angle) at q.
        // pos_rotation maps to canonical frame: d_canonical = (d + pos_rotation[p]) % 6.
        // From q: d_canonical = (d_q + pos_rotation[q]) % 6.
        // d_q = d_p + T_angle, so: (d_p + T_angle + pos_rotation[q]) = (d_p + rot_p)
        // => pos_rotation[q] = rot_p - T_angle
        pos_rotation[q] = (rot_p - T_angle + 6) % 6;
        visited.insert(q);
        bfs.push(q);
      }
    }
  }

  // Step 3: Build unified per-node direction map using rotations.
  vector<array<node_t, 6>> dir_nb(N);
  for(auto& a : dir_nb) a.fill(-1);

  for(node_t u = 0; u < N; u++){
    for(const auto& p : node_pos[u]){
      int rot = pos_rotation[p];
      const auto& nb = pos_nb[p];
      for(int d = 0; d < 6; d++){
        if(nb[d] < 0) continue;
        int d_canonical = (d + rot) % 6;
        if(dir_nb[u][d_canonical] == -1)
          dir_nb[u][d_canonical] = nb[d];
      }
    }
  }

  // Step 4: Add cross-boundary directions to dir_nb from cross_info.
  // cross_info maps each cross-boundary edge to (dir_first, dir_second, pos_first, pos_second)
  // where the directions are in each position's local Eisenstein frame.
  // Rotate them to the canonical frame using pos_rotation.
  for(const auto& [e, info] : cross_info){
    auto [dir_first_raw, dir_second_raw, pos_f, pos_s] = info;
    node_t u = e.first, v = e.second;

    int dir_uv = (dir_first_raw  + pos_rotation[pos_f]) % 6;
    int dir_vu = (dir_second_raw + pos_rotation[pos_s]) % 6;

    dir_nb[u][dir_uv] = v;
    dir_nb[v][dir_vu] = u;
  }

  // Step 5: Build neighbour lists directly from dir_nb in CW order.
  // Eisenstein directions 0-5 go CCW, so scanning 5→4→3→2→1→0 gives CW order,
  // which is what Triangulation's next_on_face (= prev) convention requires.
  neighbours_t neighbours(N, GRAPH_DMAX);
  for(node_t u = 0; u < N; u++){
    for(int d = 5; d >= 0; d--){
      if(dir_nb[u][d] >= 0)
        neighbours.push_back(u, dir_nb[u][d]);
    }
  }

  Triangulation T(neighbours);
  return T;
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
      fprintf(stderr,"outline[%d] = {{%d,%d},%d} -> u = %d. stored_up = %d, outline_newnames[%d] = %d\n",i,outline[i].first.first,outline[i].first.second,outline[i].second,u,stored_up,i,outline_newnames[i]);
      abort();
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
