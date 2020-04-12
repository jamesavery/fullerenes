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
vector<edge_t> Folding::connect_cross(const Eisenstein& w)
{
  vector<edge_t> edges;

  // Register reverse arcs
  map<dedge_t,dedgecoord_t> reverse_arc;
  vector< pair<Eisenstein,Eisenstein> > EC, ECp; // Split edge Eisenstein coordinates 

  Eisenstein xu, xv;
  node_t u,v;
  for(int i=0;i<outline.size();i++){
    tie(xu,u) = outline[i];
    tie(xv,v) = outline[(i+1) % outline.size()];

    reverse_arc[{v,u}] = {xu*w, xv*w}; // Rotated coordinates of arc v->u
  }

  Eisenstein iw = w.invertn();
  for(int i=0;i<outline.size();i++){
    tie(xu,u) = outline[i];
    tie(xv,v) = outline[(i+1) % outline.size()];

    // First get the coordinates of arc u->v and the reverse v->u
    dedgecoord_t Xuv = {xu*w,xv*w}, Xvu = reverse_arc[{u,v}];

    // What the affine transform that takes the line segment Xuv into Xvu?
    Eisenstein
      xu0,			// Coord of u in u->v
      xu1,			// Coord of u in v->u
      T;			// Rotation transform
    Unfolding::transform_line(Xuv,Xvu, xu0,xu1, T);

    // Alongside u->v, rasterize u->v line segment forwards and backwards
    vector<Eisenstein> segment   (polygon::draw_line(xu*w,xv*w)), 
                       revsegment(polygon::draw_line(xv*w,xu*v));
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
	node_t u = grid(x*iw), v = grid(yp*iw);
	edges.push_back(edge_t(u,v));

	// Debugging
	assert(u>=0 && v>=0);	
	EC.push_back({x,y});
	ECp.push_back({xp,yp});
      }
      
    }
  }

  // Debugging
  polygon P(w*get_keys(outline));
  if(debug_flags & WRITE_FILE)
    debug_file
      << "points = " << P.allpoints() << ";\n"
      << "outlinecoords = " << P.outline << ";\n" 
      << "EC = " << EC << ";\n"
      << "ECp = " << ECp << ";\n"
      ;

  return edges;
}

// connect_polygon connects all the inner edges in the outline polygon
// by exact scan-conversion and rasterization
vector<edge_t> Folding::connect_polygon(const Eisenstein& w)
{
  vector<edge_t> edges;  
  Eisenstein iw(w.invertn());
  vector<Eisenstein> outline_coords(w*get_keys(outline));
  polygon P(outline_coords);
  
  polygon::scanline S(P.scanConvert());

  for(int y=0;y<S.xs.size();y++){        // For each y..
    for(int j=0;j<S.xs[y].size()/2;j++){ // ..go through each inner interval
      int x_start = S.xs[y][2*j], x_end = S.xs[y][2*j+1];

      for(int x=x_start;x<x_end;x++){
	node_t u = grid(Eisenstein(x,y+S.minY)), v = grid(Eisenstein(x+1,y+S.minY)*iw);
	edges.push_back({u,v});

	assert(u>=0 && v>=0);
      }      
    }
  }
  return edges;
}

// The inner edges and the cross-outline edges are all the edges
vector<edge_t> Folding::connect(const Eisenstein& w)
{
  vector<edge_t> 
    polygon_edges(connect_polygon(w)), 
    cross_edges(connect_cross(w)),
    edges;
  
  if(!(debug_flags & DONT_CONNECT_POLYGON)) copy(polygon_edges.begin(),polygon_edges.end(),back_inserter(edges));
  if(!(debug_flags & DONT_CONNECT_ACROSS))  copy(cross_edges.begin(),cross_edges.end(),back_inserter(edges));

  return edges;
}


// The whole outline is connected into a triangulation / cubic graph dual
// by rotating 0, 60, and 120 degrees and "drawing" the horizontal edges
vector<edge_t> Folding::connect()
{
  Eisenstein I = {1,0}, CW = {1,-1}, CCW = {0,1};
  
  vector<edge_t> 
    e(connect(I)),
    eCW(connect(CW)),
    eCCW(connect(CCW)), edges;
  
  copy(e.begin(),   e.end(),back_inserter(edges));

  if(!(debug_flags & DONT_ROTATE)){
    copy(eCW.begin(), eCW.end(),back_inserter(edges));
    copy(eCCW.begin(),eCCW.end(),back_inserter(edges));
  }

  return edges;
}


// TODO: Shouldn't this do the 3 rotations?
vector<int> Folding::identify_nodes() const
{
  node_t u, v, U, V;
  Eisenstein xu, xv, XU, XV;  
  set<edge_t> same_as;

  map<dedge_t,dedgecoord_t> reverse_arc;


  Eisenstein I = {1,0}, CW = {1,-1}, CCW = {0,1};
  Eisenstein omegas[3] = {I,CW,CW};

  for(int i_omega=0;i_omega<3;i_omega++){
    Eisenstein omega = omegas[i_omega], omega_inv = omega.invertn();
    
    // Register reverse arcs
    for(int i=0;i<outline.size();i++){
      tie(xu,u) = outline[i];
      tie(xv,v) = outline[(i+1) % outline.size()];

      reverse_arc[{v,u}] = {xu*omega,xv*omega}; 
    }

    // For each segment U->V of the outline, find the reverse one, and identify
    // the nodes on the path that are *not* on split triangles
    for(int i=0;i<outline.size();i++){
      tie(XU,U) = outline[i];
      tie(XV,V) = outline[(i+1) % outline.size()];

      dedgecoord_t XUV(XU*omega,XV*omega), XVU(reverse_arc[{U,V}]);

      Eisenstein x0,x0p,T;
      Unfolding::transform_line(XUV,XVU, x0,x0p, T);

      //TODO: Handle horizontal lines. (?)
      vector<Eisenstein> segment(polygon::draw_line(XU,XV)), 
	revsegment(polygon::draw_line(XV,XU));
      reverse(revsegment.begin(),revsegment.end());
      assert(segment.size() == revsegment.size());

      // 
      for(int j=0;j<segment.size();j++){
	const Eisenstein& x(segment[j]), y(revsegment[j]);
	if(x == y){
	  Eisenstein xp = (x-x0)*T+x0p;
	  node_t u = grid(x*omega_inv), v = grid(xp*omega_inv);
	  same_as.insert({u,v});

	  assert(u>=0 && v>=0);	
	}
      }
    }
  }
  // Find connected components
  vector<int> same(grid.size());
  for(int i=0;i<grid.size();i++) same[i] = i;

  Graph S(same_as);
  vector<vector<node_t> > components(S.connected_components());

  for(auto& c: components){
    node_t canonical = *min_element(c.begin(),c.end());

    for(auto t: c) same[t] = canonical;
  }


  if(debug_flags&WRITE_FILE) 
    debug_file << "samecomponents = " << components << ";\n";


  return same;
}

PlanarGraph Folding::fold()
{
  vector<int> same_as;
  vector<edge_t> edge_list(connect());

  if(!(debug_flags & DONT_IDENTIFY_NODES)){
    same_as = identify_nodes();
    for(int i=0;i<edge_list.size();i++){
      node_t u = edge_list[i].first, v = edge_list[i].second;
      edge_list[i].first = same_as[u];
      edge_list[i].second = same_as[v];
    }
  }
  
  // Compactify node names
  set<edge_t> edges;
  IDCounter<int> new_names;

  if(!(debug_flags & DONT_IDENTIFY_NODES)) {
    for(int i=0;i<edge_list.size();i++){
      edge_t e = edge_list[i];
      node_t u = e.first, v = e.second;

      int newu = new_names.insert(u);
      int newv = new_names.insert(v);

      edges.insert(edge_t(newu,newv));
    }
  } else {
    edges = set<edge_t>(edge_list.begin(),edge_list.end());
    for(int i=0;i<grid.size();i++) new_names.insert(i);
  }

  int N = new_names.size();
  vector<coord2d> layout2d;

  if(debug_flags & DO_NONPLANAR_LAYOUT){
    layout2d.resize(N);
    for(IDCounter<Eisenstein>::const_iterator xi(grid.begin()); xi!=grid.end(); xi++){
      Eisenstein xy(xi->first);
      node_t      u(new_names(xi->second));

      if(u>=0) layout2d[u] = coord2d(xy.first,xy.second);
    }
  }

  PlanarGraph G(edges,layout2d);

  if(debug_flags & WRITE_FILE) debug_file << "G = " << G << ";\n";

  return G;
}

vector<node_t> Folding::outline_nodes() const
{
  node_t u, outline_N=0;
  Eisenstein xu;
  
  vector<node_t> same_nodes(identify_nodes());
  vector<node_t> outline_newnames(outline.size());

  for(int i=0;i<outline.size();i++){
    tie(xu,u) = outline[i];
    outline_newnames[i] = grid(xu);
    outline_N = max(outline_N,u+1);
  }
  
  vector<node_t> new_nodenames(outline_N,-1);
  for(int i=0;i<outline.size();i++){
    cerr << new_nodenames << endl;

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
