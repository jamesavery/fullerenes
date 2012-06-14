#include "planargraph.hh"

bool PlanarGraph::this_is_a_fullerene() const {
  for(node_t u=0;u<N;u++)
    if(neighbours[u].size() != 3){ 
      fprintf(stderr,"Graph is not cubic: vertex %d has %d neighbours.\n",u,int(neighbours[u].size())); 
      return false;
    }
    
  facemap_t faces(PlanarGraph(*this).compute_faces(7));
  int n_faces = 0;
  for(facemap_t::const_iterator f(faces.begin()); f!=faces.end();f++)
    n_faces += f->second.size();

  const int E = 3*N/2;
  const int F = 2+E-N;
    
  if(E != edge_set.size()){
    fprintf(stderr,"Graph is not planar: wrong number of edges: %d != %d\n",int(edge_set.size()),E);
    return false;
  }

  if(F != n_faces){
    fprintf(stderr,"Graph is not planar: wrong number of faces: %d != %d\n",n_faces,F);
    return false;
  }

  if(faces[5].size() != 12){
    fprintf(stderr,"Graph is not fullerene: wrong number of pentagons: %d != 12\n",int(faces[5].size()));
    return false;
  }

  if(faces[6].size() != (F-12)){
    fprintf(stderr,"Graph is not fullerene: wrong number of hexagons: %d != %d\n",int(faces[6].size()),F-12);
    return false;
  }

  return true;
}


PlanarGraph PlanarGraph::dual_graph(unsigned int Fmax) const {
  PlanarGraph dual;
  unsigned int Nfaces = edge_set.size()-N+2;
  dual.N = Nfaces;
  dual.neighbours.resize(Nfaces);
  dual.edges.resize(Nfaces*(Nfaces-1)/2);
  
  //  cerr << "dual_graph(" << Fmax << ")\n";
  const vector<face_t> allfaces(compute_faces_flat(Fmax));

  if(Nfaces != allfaces.size()){
    fprintf(stderr,"%d != %d faces: Graph is not polyhedral.\n",Nfaces,int(allfaces.size()));
    cout << "errgraph = " << *this << endl;
  }

  // Construct mapping e -> faces containing e (these are mutually adjacent)
  //  cerr << "dual_graph::construct facenodes\n";
  map< edge_t, set<int> > facenodes;
  for(unsigned int i=0;i<allfaces.size(); i++){
    const face_t& face(allfaces[i]);
    //  cerr << "Face "<<i<<": " << face << endl;
    for(unsigned int j=0;j<face.size();j++)
      facenodes[edge_t(face[j],face[(j+1)%face.size()])].insert(i);
  }
  //  cerr << "dual_graph::test planarity\n";
  for(map<edge_t,set<int> >::const_iterator fs(facenodes.begin());fs!=facenodes.end();fs++){
    const edge_t&   e(fs->first);
    const set<int>& connects(fs->second);
    if(connects.size() != 2)
      fprintf(stderr,"Edge (%d,%d) connects %d faces: Graph is not planar.\n",e.first,e.second,int(connects.size()));
  }
  
  // Insert edge between each pair of faces that share an edge
  //  cerr << "dual_graph::construct graph\n";
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const set<int>& adjacent_faces(facenodes[*e]);
    for(set<int>::const_iterator f(adjacent_faces.begin()); f!= adjacent_faces.end(); f++){
      set<int>::const_iterator g(f);
      for(++g; g!= adjacent_faces.end(); g++)
	dual.edge_set.insert(edge_t(*f,*g));
    }
  }
  //fprintf(stderr,"%d nodes, and %d edges in dual graph.\n",int(dual.N), int(dual.edge_set.size()));

  dual.update_auxiliaries();

  // If original graph was planar with 2D layout, there's a corresponding layout for the dual graph
  if(layout2d.size() == N){
    //    cerr << "dual_graph::compute layout.\n";
    dual.layout2d = vector<coord2d>(Nfaces);
#pragma omp parallel for
    for(int i=0;i<Nfaces;i++)
      dual.layout2d[i] = allfaces[i].centroid(layout2d);
  }
  return dual;
}



// NB: TODO: What happens, for example, if a triangle is comprised of three smaller triangles?
// This produces "phantom" faces! Fix and use the oriented version instead.
facemap_t PlanarGraph::compute_faces(unsigned int Nmax, bool planar_layout) const 
{
  facemap_t facemap;
  // TODO: This is a much better and faster method, but needs to be debugged.
  if(planar_layout && layout2d.size() == N) return compute_faces_oriented();

  cerr << "Non-oriented face computation (loop search)\n";
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const node_t s = e->first, t = e->second;

    const vector<node_t>& ns(neighbours[t]);

    for(unsigned int i=0;i<ns.size();i++)
      if(ns[i] != s) {
	const node_t u = ns[i];

	face_t face(shortest_cycle(s,t,u,Nmax));  
	//	cerr << face << endl;
	if(face.size() > 0 && face.size() <= Nmax){
	  facemap[face.size()].insert(face);
	} //else {
	  //	  fprintf(stderr,"Erroneous face starting at (%d -> %d -> %d) found: ",s,t,u); 
	  //	  cerr << face << endl;
	  
	//	}
      }
  }
  return facemap;
}

facemap_t PlanarGraph::compute_faces_oriented() const 
{
  assert(layout2d.size() == N);
  facemap_t facemap;
  cerr << "Computing faces using 2D orientation.\n";
  set<dedge_t> workset;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const node_t s = e->first, t = e->second;
    workset.insert(dedge_t(s,t));
    workset.insert(dedge_t(t,s));
  }

  // Outer face must exist and be ordered CW, rest of faces CCW
  // assert(outer_face.size() > 0);
  if(outer_face.size() < 3)
    outer_face = find_outer_face();

  if(outer_face.size() < 3){
    cerr << "Invaid outer face: " << outer_face << endl;
    assert(outer_face.size() < 3);
  }

  //  cerr << "Outer face: " << outer_face << endl;

  coord2d centre(centre2d(layout2d));
  
  for(node_t u=0;u<N;u++)
    if(!outer_face.contains(u) && !outer_face.point_inside(layout2d,u)){
      cerr << "Point " << u << "/" << layout2d[u] << " is outside outer face " << outer_face << endl;
      cerr << "Winding number: " << outer_face.winding_number(layout2d,u) << endl;
      abort();
    }
  //  cerr << "Outer face is OK: All vertices are inside face.\n";

  // Add outer face to output, remove directed edges from work set
  facemap[outer_face.size()].insert(outer_face);
  for(unsigned int i=0;i<outer_face.size();i++){
    const node_t u = outer_face[i], v = outer_face[(i+1)%outer_face.size()];
    //    printf("Removing directed edge (%d,%d)\n",u,v);
    workset.erase(dedge_t(u,v));
  }

  // Now visit every other edge once in each direction.
  while(!workset.empty()){
    dedge_t e = *workset.begin(); workset.erase(workset.begin());

    face_t face;
    face.push_back(e.first);
    face.push_back(e.second);

    //    printf("%d->%d\n",e.first,e.second);
    while(e.second != face[0]){
      const node_t u = e.first, v = e.second;
      const vector<node_t>& ns(neighbours[v]);

      coord2d vu(layout2d[u]-layout2d[v]);
      double angle_min = -M_PI;

      node_t w=-1;
      for(unsigned int i=0;i<ns.size();i++) {
	//	printf("%d : %d (%d->%d) angle %g\n",i,ns[i],u,v,vu.line_angle(layout[ns[i]]-layout[v]));
	if(ns[i] != u) { // Find and use first unvisited edge in order of angle to u->v
	  set<dedge_t>::iterator ei(workset.find(dedge_t(v,ns[i])));

	  if(ei != workset.end()){ // directed edge is not yet visited
	    coord2d vw(layout2d[ns[i]]-layout2d[v]);
	    double angle = vu.line_angle(vw);

	    if(angle>= angle_min){
	      angle_min = angle;
	      w = ns[i];
	    } 
	  } 
	}
      }
      if(w == -1) abort(); // There is no face!

      e = dedge_t(v,w);
      workset.erase(e);
      
      if(e.second != face[0]) face.push_back(e.second);

    }
    //    cout << "face = " << face << endl;
    facemap[face.size()].insert(face);
  }
  return facemap;
}



vector<face_t> PlanarGraph::compute_faces_flat(unsigned int Nmax, bool planar_layout) const 
{
  vector<face_t> faces;
  facemap_t facemap(compute_faces(Nmax,planar_layout));

  for(facemap_t::const_iterator fs(facemap.begin()); fs != facemap.end(); fs++)
    copy(fs->second.begin(),fs->second.end(),inserter(faces,faces.end()));

  // Make sure that outer face is at position 0
  if(planar_layout){
    const node_t s(outer_face[0]), t(outer_face[1]), r(outer_face[2]);
    for(int i=0;i<faces.size();i++){
      const face_t f(faces[i]);
      if(f[0] == s && f[1] == t && f[2] == r){ // swap faces[i] with faces[0]
	faces[i] = faces[0];
	faces[0] = f;
      }
    }
  } else outer_face = face_t(faces[0]);

  return faces;
}


vector<tri_t> PlanarGraph::triangulation(int face_max) const
{
  vector<face_t> faces(compute_faces_flat(face_max));  
  return triangulation(faces);
}

vector<tri_t> PlanarGraph::centroid_triangulation(const vector<face_t>& faces) const 
{
  // Test whether faces already form a triangulation
  bool is_tri = true; for(int i=0;i<faces.size();i++) if(faces[i].size() != 3) is_tri = false;
  if(is_tri){
    //    cerr << "Faces already form a triangulation.\n";
    vector<tri_t> tris(faces.begin(),faces.end());
    return orient_triangulation(tris);
  }

  // Triangulate by inserting extra vertex at face centroid and connecting
  // each face vertex to this midpoint.
  vector<tri_t> tris;
  for(int i=0;i<faces.size();i++){
    const node_t v_new = N+i;
    const face_t& f(faces[i]);

    for(int j=0;j<f.size();j++)
      tris.push_back(tri_t(f[j],v_new,f[(j+1)%f.size()]));
  }
  
  return orient_triangulation(tris);
}
  

vector<tri_t> PlanarGraph::triangulation(const vector<face_t>& faces) const
{
  // Test whether faces already form a triangulation
  bool is_tri = true; for(int i=0;i<faces.size();i++) if(faces[i].size() != 3) is_tri = false;
  if(is_tri){
    //    cerr << "Faces already form a triangulation.\n";
    vector<tri_t> tris(faces.begin(),faces.end());
    return orient_triangulation(tris);
  } else {
    for(int i=0;i<faces.size();i++) 
      if(faces[i].size() != 3){
	fprintf(stderr,"Face %d has %d sides: ",i,int(faces[i].size())); cerr << faces[i] << endl;
      }
  }

  vector<tri_t> tris;
  // First, break up the faces into a non-consistent triangulation
  for(size_t i=0;i<faces.size();i++){
    face_t f(faces[i]);
    assert(f.size() >= 3);
    for(size_t j=1;j<f.size()-1;j++)
      tris.push_back(tri_t(f[0],f[j],f[j+1]));
  }
  
  return orient_triangulation(tris);
}


vector<tri_t>& PlanarGraph::orient_triangulation(vector<tri_t>& tris) const
{
  // Now, pick an orientation for triangle 0. We choose the one it
  // already has. This determines the orientation of the remaining triangles!
  map<dedge_t,bool> done;
  for(int i=0;i<3;i++) done[dedge_t(tris[0][i],tris[0][(i+1)%3])] = true;

  for(int i=1;i<tris.size();i++){
    tri_t& t(tris[i]);
    if(done[dedge_t(t[0],t[1])] || done[dedge_t(t[1],t[2])] || done[dedge_t(t[2],t[0])]){
      node_t u = t[2]; t[2] = t[1]; t[1] = u;
    }
    
    // if(done[dedge_t(t[0],t[1])]){ node_t u = t[1]; t[1] = t[0]; t[0] = u; }
    // if(done[dedge_t(t[1],t[2])]){ node_t u = t[2]; t[2] = t[1]; t[1] = u; }
    // if(done[dedge_t(t[2],t[0])]){ node_t u = t[0]; t[0] = t[2]; t[2] = u; }
    done[dedge_t(t[0],t[1])] = true;
    done[dedge_t(t[1],t[2])] = true;
    done[dedge_t(t[2],t[0])] = true;
  }
  // Check consistency of orientation. It is consistent if and only if
  // each edge has been used exactly once in each direction.
  bool consistent = true;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    if(!done[dedge_t(e->first,e->second)]){
      fprintf(stderr,"A: Directed edge %d->%d is missing: triangulation is not consistently oriented.\n",e->first,e->second);
      consistent = false;
    }
    if(!done[dedge_t(e->second,e->first)]){
      fprintf(stderr,"B: Directed edge %d->%d is missing: triangulation is not consistently oriented.\n",e->second,e->first);
      consistent = false;
    }
  }

  if(!consistent){
    cerr << "(*** Inconsistent triangulation: ***)\n";
    cerr << "tris = {"; for(int i=0;i<tris.size();i++) cerr << tris[i] << (i+1<tris.size()? ", ":"};\n");
    cerr << "outerface = " << outer_face << ";\n";
  }
  assert(consistent == true);
  return tris;
}

// Finds the vertices belonging to the outer face in a symmetric planar
// layout centered at (0,0). Returns the face in CW order.
face_t PlanarGraph::find_outer_face() const 	
{
  assert(layout2d.size() == N);

  vector<double> radii(N);

  node_t u_farthest = 0;
  double rmax = 0;
  for(node_t u=0;u<N;u++){
    radii[u] = layout2d[u].norm();
    if(radii[u] > rmax){ rmax = radii[u]; u_farthest = u; }
  }
  
  face_t outer_face;
  int i = 0;
  for(node_t t = u_farthest, u = u_farthest, v = -1; v != u_farthest && i <= N; i++){
    const vector<node_t>& ns(neighbours[u]);
    double r = 0;
    for(int i=0;i<ns.size();i++)
      if(ns[i] != t && ns[i] != u && radii[ns[i]] > r){ r = radii[ns[i]]; v = ns[i]; }
    outer_face.push_back(u);
    t = u;
    u = v;
  }
  // fprintf(stderr,"(u_farthest,rmax) = (%d,%f); i = %d\n",u_farthest,rmax,i);
  // cerr << "Outer face: " << outer_face << endl;
  // cerr << "Radii: "; for(int i=0;i<outer_face.size();i++) cerr << " " << radii[outer_face[i]]; cerr << "\n";

  assert(i<N);

  sort_ccw_point CCW(layout2d,outer_face.centroid(layout2d));
  sort(outer_face.begin(),outer_face.end(),CCW);
  reverse(outer_face.begin(),outer_face.end());  

  cerr << "Found outer face: " << outer_face << endl;
  return outer_face;
}

vector<double> PlanarGraph::edge_lengths() const 
{
  assert(layout2d.size() == N);

  vector<double> lengths(edge_set.size());
  unsigned int i = 0;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end();e++, i++)
    lengths[i] = (layout2d[e->first]-layout2d[e->second]).norm();

  return lengths;
}

coord2d PlanarGraph::width_height() const {
  double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY;
  for(node_t u=0;u<N;u++){
    double x = layout2d[u].first, y = layout2d[u].second;
    if(x<xmin) xmin = x;
    if(x>xmax) xmax = x;
    if(y<ymin) ymin = y;
    if(y>ymax) ymax = y;
  }
  return coord2d(xmax-xmin,ymax-ymin);
}

void PlanarGraph::scale(const coord2d& x) {
  for(node_t u=0;u<N;u++) layout2d[u] *= x;
}

void PlanarGraph::move(const coord2d& x) {
  for(node_t u=0;u<N;u++) layout2d[u] += x;
}


ostream& operator<<(ostream& s, const PlanarGraph& g) 
{
  s << g.name<< "Graph[Range[0,"<<(g.N-1)<<"],\n\tUndirectedEdge@@#&/@{";
  for(set<edge_t>::const_iterator e(g.edge_set.begin()); e!=g.edge_set.end(); ){    
    s << "{" << e->first << "," << e->second << "}";
    if(++e != g.edge_set.end())
      s << ", ";
    else
      s << "}";
  }

  if(g.layout2d.size() == g.N){
    s << g.name << ",\n\tVertexCoordinates->{";
    for(unsigned int i=0;i<g.N;i++){
      coord2d xy(g.layout2d[i]);
      s << xy << (i+1<g.N?", ":"}");
    }
  } // else { fprintf(stderr,"No layout, man!\n"); }
  s << "\n]";

  return s;
}
