#include "graph.hh"

node_t observe_point = 8;

struct Polyhedron : public Graph {
  vector<coord3d> points;
  coord3d centre;
  vector<face_t> faces;


  Polyhedron(const Graph& G, const vector<coord3d>& points, const int face_max = INT_MAX) : 
    Graph(G), points(points), centre(centre3d(points)), faces(G.compute_faces_flat(face_max))
  {
    // TODO: Currently assumes layout2d is given -- but it is not!
  }

  double surface_area() const {
    assert(layout2d.size() == N);
    vector<face_t> tris(triangulation(faces).compute_faces_flat(4,layout2d));
    double A = 0;  
    for(size_t i=0;i<tris.size();i++){
      const face_t& tri(tris[i]);
      Tri3D T(points[tri[0]],points[tri[1]],points[tri[2]]);
      A += T.area();
    }    
    return A;
  }

  double volume() const {
    assert(layout2d.size() == N);
    vector<face_t> tris(triangulation(faces).compute_faces_flat(4,layout2d));
    double V = 0;

    // Sort triangle vertices CCW according to 2D layout (N.B.: may have to reverse triangles generated from face 0! How?
#if 0
    for(size_t i=0;i<tris.size();i++){
      coord2d c((layout2d[tris[i][0]]+layout2d[tris[i][0]]+layout2d[tris[i][0]])/3.0);
      sort_ccw_point CCW(layout2d,c);
      sort(tris[i].begin(),tris[i].end(),CCW);
      //      cout << tris[i] << endl;
    }
#endif
    
    // Now generate tetrahedra and either add or subtract volume according to which direction the face is pointing
    coord3d zero;
    for(size_t i=0;i<tris.size();i++){
      const face_t& tri(tris[i]);
      Tri3D T(points[tri[0]],points[tri[1]],points[tri[2]]);

      Tetra3D tet(T.a,T.b,T.c,centre);
      if(T.back_face(centre)) V += tet.volume();
      else V += tet.volume();
    }
    return V;
  }

  // TODO: Handle coplanar triangles.
  bool include_point(const node_t& u, const coord3d& interior_point) const {
   const vector<node_t>& ns(neighbours[u]);
      const coord3d& ux(points[u]);
      bool include_point = true;

      if(observe_point == u){
	cout << "centre = " << interior_point << ";\n";
	cout << "uindex = " << u << ";\n";
	cout << "upoint = " << points[u] << ";\n";
	cout << "uneighbours = {"; for(int i=0;i<ns.size();i++) cout << points[ns[i]] << (i+1<ns.size()? ", ":"};\n");
      }
      for(unsigned int i=0;i<ns.size();i++)
	for(unsigned int j=i+1;j<ns.size();j++)
	  for(unsigned int k=j+1;k<ns.size();k++){
	    const coord3d vs[3] = {points[ns[i]],points[ns[j]],points[ns[k]]};
	    int intersects = -1;
	    // Test whether Tetrahedron(u,v1,v2,v3) is contained in the hull 
	    Tri3D T(vs);
	    const Tri3D::segment_t line(interior_point,(vs[0]+vs[1]+vs[2])/3.0);
	    
	    if(observe_point == u){
	      cout << "cpoint"<<i<<j<<k<<" = " << ((vs[0]+vs[1]+vs[2])/3.0) << ";\n";
	      cout << "Tetra"<<i<<j<<k<<" = {"<<Tri3D(vs) << ", ";
	    }
	    for(int p=0;p<3;p++){
	      coord3d uvv[3] = vs;
	      uvv[p] = ux;
	      if(observe_point == u)
		cout << Tri3D(uvv) << (p+1<3? ", ":"};\n");
	      if(Tri3D(uvv).intersects(line))
		intersects = p;
	    }
	    if(intersects >= 0)
	      fprintf(stderr,"%d - (%d,%d,%d) - intersects %d\n",u,i,j,k,intersects);
	    if(observe_point == u && intersects>=0)
	      cout << "intersects = " << intersects << ";\n";
	    if(intersects >= 0){
	      include_point = false;
	      cout << "itets = Tetras" << i << j << k << ";\n";
	      cout << "cpoint = cpoint" << i << j << k << ";\n";
	    }

	  }
      return include_point;
  }



  void remove_vertex(const node_t& u) {
    vector<node_t> ns(neighbours[u]);
    sort_ccw_point CCW(layout2d,layout2d[u]);
    sort(ns.begin(),ns.end(),CCW);    
    
    neighbours[u].clear();
    for(unsigned int i=0;i<ns.size();i++){
      const node_t& vleft(ns[(i+ns.size()-1)%ns.size()]), v(ns[i]), vright(ns[(i+1)%ns.size()]);
      set<node_t> nvs(neighbours[v].begin(),neighbours[v].end());
      nvs.erase(u);
      nvs.insert(vleft);
      nvs.insert(vright);
      
      neighbours[v] = vector<node_t>(nvs.begin(), nvs.end());
    }
  }

  Polyhedron& convex_hull_inplace() {
    set<node_t> workset;
    for(node_t u=0;u<N;u++) workset.insert(u); // n log n, can be made linear

    int step=0;
    while(!workset.empty()){
      node_t u = *workset.begin(); workset.erase(workset.begin()); // pop one element from the working set
     
      if(!include_point(u,centre)){
	fprintf(stderr,"Discarding node %d.\n",u);
	for(unsigned int i=0;i<neighbours[u].size();i++){
	  fprintf(stderr,"Inserting neighbour %d back into work set.\n",neighbours[u][i]);
	  workset.insert(neighbours[u][i]);
	}
	remove_vertex(u);
	update_from_neighbours();
	cout << "step" << (++step) << " = " << *this << ";\n";
      } else {
	fprintf(stderr,"Keeping node %d for now.\n",u);
      }
    }
    // TODO: Move to method
    edge_set.clear();
    for(node_t u=0;u<N;u++){
      for(unsigned int i=0;i<neighbours[u].size();i++)
	edge_set.insert(edge_t(u,neighbours[u][i]));
    }
    update_auxiliaries();

    return *this;
  }

  Polyhedron convex_hull() const {
    Polyhedron hull(*this);
    return hull.convex_hull_inplace();
  }

  friend ostream& operator<<(ostream& s, const Polyhedron& P){
    vector<node_t> reachable_points;
    for(node_t u=0;u<P.N;u++) if(P.neighbours[u].size()!=0) reachable_points.push_back(u);
    s << "{{";
    for(unsigned int i=0;i<reachable_points.size();i++) s << reachable_points[i] << (i+1<reachable_points.size()?", ":"},{");
    for(unsigned int i=0;i<P.points.size();i++) s << P.points[i] << (i+1<P.points.size()?", ":"},");
    s << static_cast<Graph>(P) << "}";
    return s;
  }
};

