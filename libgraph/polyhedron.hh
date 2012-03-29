#include "graph.hh"
#include <fstream>
#include <sstream>

node_t observe_point = 8;

struct Polyhedron : public Graph {
  int face_max;
  vector<coord3d> points;
  coord3d centre;
  vector<face_t> faces;

  Polyhedron(const int face_max=INT_MAX) : face_max(face_max) {}
  Polyhedron(const Graph& G, const vector<coord3d>& points, const int face_max = INT_MAX) : 
    Graph(G), face_max(face_max), points(points), centre(centre3d(points)), faces(G.compute_faces_flat(face_max))
  {
    //    if(layout2d.size() != N)
    //layout2d = polar_angles();
    layout2d = tutte_layout();
  }
  Polyhedron(const string& path) {
    ifstream file(path.c_str());
    file >> *this;
    file.close();
  }

  vector<coord2d> polar_angles() const {
    vector<coord2d> angles(N);
    for(node_t u=0;u<N;u++){
      const coord3d& x(points[u]);
      const double r = x.norm();
      angles[u].first  = acos(x[2]/r);
      angles[u].second = atan2(x[1],x[0]);
    }
    return angles;
  }

  double surface_area() const {
    assert(layout2d.size() == N);
    double A = 0;  

    vector<face_t> tris(triangulation(faces));

    for(size_t i=0;i<tris.size();i++){
      const face_t& tri(tris[i]);
      Tri3D T(points[tri[0]],points[tri[1]],points[tri[2]]);
      A += T.area();
    } 

    return A;
  }

  double volume() const {
    vector<face_t> tris(triangulation(faces));
    double V = 0;

    // Now generate tetrahedra and either add or subtract volume according to which direction the face is pointing
    coord3d zero(-10,0,0);
    for(size_t i=0;i<tris.size();i++){
      const face_t& t(tris[i]);
      Tri3D T(points[t[0]],points[t[1]],points[t[2]]);
      double dV = Tetra3D(T.a,T.b,T.c,zero).volume();
      V += (T.back_face(zero)? 1 : -1)*dV;
    }
    return fabs(V);
  }

  double volume2() const {
    vector<face_t> tris(triangulation(faces));
    double V = 0;

    // Now generate tetrahedra and either add or subtract volume according to which direction the face is pointing
    for(size_t i=0;i<tris.size();i++){
      const face_t& t(tris[i]);
      Tri3D T(points[t[0]],points[t[1]],points[t[2]]);

      V += ((T.a).dot(T.n))*T.area()/T.n.norm();
    }
    return fabs(V/3.0);
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

  friend istream& operator>>(istream& f, Polyhedron& P){
    string s;
    node_t u=0,v;
    coord3d x;

    while(getline(f,s)){
      stringstream l(s);
      l >> x;
      if(l.fail()) continue; // Invalid line
      while(!l.fail()){
	l >> v;
	P.edge_set.insert(edge_t(u,v-1)); // File format numbers nodes from 1 
      }
      P.points.push_back(x);
      u++;
    }
    P.update_auxiliaries();
    P.layout2d = P.tutte_layout();
    P.faces = P.compute_faces_flat(P.face_max,P.layout2d);
    P.centre = P.centre3d(P.points);
    return f;
  }
};

