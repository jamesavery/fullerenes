#include "fullerenes/unfold.hh"

/************************************************************/ 
/************************************************************/ 
/*                 FOLDING IMPLEMENTATION                   */
/************************************************************/ 
/************************************************************/ 


// TODO: Use Eisenstein layout to insert arcs correctly oriented directly
//       Replace omega by i_omega={0,1,2}
//
//       u--v--w places n[v][i_omega] = w, n[v][i_omega+3] = u, so
//       u--v    places n[u][i_omega] = v, n[v][i_omega+3] = u
// HOW TO IMPLEMENT:
//       void connect_xxx(Graph &n, int i_omega)
//       Fill in n


// connect_cross connects edges in split triangles
// (for outline segments that do not align with Eisenstein grid)
void Folding::connect_cross(int i_omega, set<edge_t> &edges)
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
    assert(segment.size() == revsegment.size());

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
	if(debug_flags & WRITE_FILE) printf("Connect cross arc %d to %d \n",u,v);
	edges.insert(edge_t(u,v));
      }
    }
  }
}

// connect_polygon connects all the inner edges in the outline polygon
// by exact scan-conversion and rasterization
void Folding::connect_polygon(int i_omega, set<edge_t> &edges)
{
  Eisenstein
    omega     = Eisenstein::unit[i_omega],
    omega_inv = Eisenstein::unit[6-i_omega];

  vector<Eisenstein> outline_coords(omega*get_keys(outline));
  polygon P(outline_coords);
  polygon::scanline S(P.scanConvert());

  for(int y=0;y<S.xs.size();y++){        // For each y..
    for(int j=0;j<S.xs[y].size()/2;j++){ // ..go through each inner interval
      int x_start = S.xs[y][2*j], x_end = S.xs[y][2*j+1];

      for(int x=x_start;x<x_end;x++){
	Eisenstein pu = Eisenstein(x,y+S.minY)*omega_inv, pv = Eisenstein(x+1,y+S.minY)*omega_inv;
	auto itu = final_grid.find(pu), itv = final_grid.find(pv);
	if(itu == final_grid.end() || itv == final_grid.end()) continue;

	node_t u = itu->second, v = itv->second;
	if(debug_flags & WRITE_FILE) printf("Connect polygon arc %d to %d \n",u,v);
	edges.insert(edge_t(u,v));
      }
    }
  }
}

// The inner edges and the cross-outline edges are all the edges
void Folding::connect(int i_omega, set<edge_t> &edges)
{
  if(!(debug_flags & DONT_CONNECT_POLYGON))  connect_polygon(i_omega, edges);
  if(!(debug_flags & DONT_CONNECT_ACROSS))   connect_cross(i_omega, edges);
}


// The whole outline is connected into a triangulation / cubic graph dual
// by rotating 0, 60, and 120 degrees and "drawing" the horizontal edges
void Folding::connect(set<edge_t> &edges)
{
  connect(0, edges);
  connect(1, edges);
  connect(2, edges);
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
    if(debug_flags & WRITE_FILE) cout << "omega: " << i_omega << ": " << omega << "; " << omega_inv << endl;
    // Register reverse arcs
    for(int i=0;i<outline.size();i++){
      tie(XU,U) = outline[i];
      tie(XV,V) = outline[(i+1) % outline.size()];

      reverse_arc[{V,U}] = {XU*omega,XV*omega}; 
    }

    // cout << "reverse_keys   = " << get_keys(reverse_arc) << endl;
    // cout << "reverse_values = " << get_values(reverse_arc) << endl;    
    
    // For each segment U->V of the outline, find the reverse one, and identify
    // the nodes on the path that are *not* on split triangles
    for(int i=0;i<outline.size();i++){
      tie(XU,U) = outline[i];
      tie(XV,V) = outline[(i+1) % outline.size()];

      arccoord_t XUV(XU*omega,XV*omega), XVU(reverse_arc[{U,V}]);

      if(debug_flags & WRITE_FILE){
        cout << "\noutline["<<i<<"]\n";
        cout << "{U,V} = " << make_pair(U,V) << endl;
        cout << "{XU,XV}   = " << make_pair(XU,XV) << endl;
        cout << "XUV = " << XUV << endl;
        cout << "XVU = " << XVU << endl;
      }
      
      Eisenstein x0,x0p,T;
      Unfolding::transform_line(XUV,XVU, x0,x0p, T);

      //      cout << "{x0,x0p} = " << make_pair(x0,x0p) << endl;
      //TODO: Handle horizontal lines. 
      vector<Eisenstein>
	segment   (polygon::draw_line(XU*omega, XV*omega)), 
	revsegment(polygon::draw_line(XV*omega, XU*omega));

      if(debug_flags & WRITE_FILE)
        cout << "segment = " << segment    << endl
	     << "regveg  = " << revsegment << endl;
      reverse(revsegment.begin(),revsegment.end());
      assert(segment.size() == revsegment.size());

      // 
      for(int j=0;j<segment.size();j++){
	const Eisenstein& x(segment[j]), y(revsegment[j]);
	if(x == y){
	  if(debug_flags & WRITE_FILE) cout << "x = " << x << ", u = " << grid(x*omega_inv) << endl;
	  //	  cout << "{x,y} = " << make_pair(x,y) << endl;
	  //	  cout << "{xw,yw} = " << make_pair(x*omega_inv,y*omega_inv) << endl;	  
	  Eisenstein xp = (x-x0)*T+x0p;
	  //	  cout << "xp = " << xp << endl;

	  node_t u = grid(x*omega_inv), v = grid(xp*omega_inv);
	  if((u>=0 && v>= 0)  && (u != v)){
	    if(debug_flags & WRITE_FILE) cout << "same_as {u,v} = " << edge_t{u,v} << endl;
	    same_as.insert(edge_t{u,v});
	  } else {
	    if(debug_flags & WRITE_FILE) cout << "u==v at "<<(x*omega_inv)<<"/"<<(xp*omega_inv)<<": " << edge_t{u,v} << endl;
	  }
	    //	  assert(u>=0 && v>=0);	
	}
      }
    }
  }
  if(debug_flags & WRITE_FILE) cout << "same_as = " << same_as << "\n\n";
  
  // Find connected components
  vector<int> same(grid.size());
  for(int i=0;i<grid.size();i++) same[i] = i;

  Graph S(same_as);
  vector<vector<node_t> > components(S.connected_components());
  if(debug_flags & WRITE_FILE) cout << "S = " << S << endl;

  for(auto& c: components){
    node_t canonical = *min_element(c.begin(),c.end());
    sort(c.begin(), c.end());
    if(debug_flags & WRITE_FILE) cout << "Component " << c << " has canonical element " << canonical << endl;
      
    for(auto t: c) same[t] = canonical;
  }


  if(debug_flags&WRITE_FILE) 
    debug_file << "samecomponents = " << components << ";\n";

  if(debug_flags & WRITE_FILE)
    cout << "same_components = " << components.size()<<", " << components << "\n\n"
         << "same    = " << same << "\n\n";
  return same;
}

Triangulation Folding::fold()
{
  node_t N = node_pos.size();
  set<edge_t> edge_set;

  connect(edge_set);

  // Build neighbour lists from the collected edge set.
  // Using a set avoids both duplicates and clobbering that the old
  // 6-slot positional scheme suffered from with multi-position nodes.
  neighbours_t neighbours(N);
  for(const auto& e: edge_set){
    neighbours[e.first].push_back(e.second);
    neighbours[e.second].push_back(e.first);
  }

  Triangulation T(neighbours,true);
  // TODO:
  // 1. CubicPair -> also dual
  // 2. uv-map for both
  // 3. include for polyhedron
  // 4. multiplication -> keep old nodes in same spot, interpolation of uv-map, 3D coordinates, 2D coordinates, anything really
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
  // printf("\\newcommand{\\outline}{");
  // for(int i=0;i<outline.size();i++){
  //   auto o = outline[i];
  //   printf("%d/%d/%d%s",o.second,o.first.first,o.first.second,i+1==outline.size()?"":",");
  // }
  // printf("}\n\n");
  // printf("\\newcommand{\\xzero}{0}\n"
  // 	 "\\newcommand{\\yzero}{0}\n"
  // 	 "\\newcommand{\\nx}{16}\n"
  // 	 "\\newcommand{\\ny}{16}\n\n");
  
  // latex_scanconversion("A",{1,0}, p);
  // latex_scanconversion("B",{0,1}, p);
  // latex_scanconversion("C",{1,-1},p);    

  // int N = -1;
  // matrix<int> node_grid          = node_grid_from_outline(outline,N);
  // vector<Eisenstein> inner_faces = cubic_inner_faces(p, node_grid);

  // printf("\\newcommand{\\innernodes}{");
  // for(int i=0;i<inner_faces.size();i++)
  //   printf("%d/%d%s",inner_faces[i].first,inner_faces[i].second,i+1<inner_faces.size()?",":"");
  // printf("}\n");

  // printf("\\newcommand{\\innerhexagons}{");
  // for(int i=0;i<inner_faces.size();i++){
  //   Eisenstein xy = inner_faces[i], omega = {0,1}, omega_n = {1,0};
  //   double cx, cy;

  //   for(int j=0;j<6;j++, omega_n *= omega){
  //     Eisenstein C = (xy) + (xy+omega_n) + (xy+omega_n*omega);

  //     printf("%.2f/%.2f%s",C.first/3.0, C.second/3.0, j+1<6?"/":"");
  //   }
  //   if(i+1<inner_faces.size()) printf(",\n");
  // }
  // printf("}\n");    
  return s;
}
