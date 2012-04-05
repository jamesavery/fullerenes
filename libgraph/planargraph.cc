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

  const vector<face_t> allfaces(compute_faces_flat(Fmax));

  if(Nfaces != allfaces.size()){
    fprintf(stderr,"%d != %d faces: Graph is not polyhedral.\n",Nfaces,int(allfaces.size()));
    cout << "errgraph = " << *this << endl;
  }

  // Construct mapping e -> faces containing e (these are mutually adjacent)
  map< edge_t, set<int> > facenodes;
  for(unsigned int i=0;i<allfaces.size(); i++){
    const face_t& face(allfaces[i]);
    //  cerr << "Face "<<i<<": " << face << endl;
    for(unsigned int j=0;j<face.size();j++)
      facenodes[edge_t(face[j],face[(j+1)%face.size()])].insert(i);
  }

  for(map<edge_t,set<int> >::const_iterator fs(facenodes.begin());fs!=facenodes.end();fs++){
    const edge_t&   e(fs->first);
    const set<int>& connects(fs->second);
    if(connects.size() != 2)
      fprintf(stderr,"Edge (%d,%d) connects %d faces: Graph is not planar.\n",e.first,e.second,int(connects.size()));
  }
  
  // Insert edge between each pair of faces that share an edge
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const set<int>& adjacent_faces(facenodes[*e]);
    for(set<int>::const_iterator f(adjacent_faces.begin()); f!= adjacent_faces.end(); f++){
      set<int>::const_iterator g(f);
      for(++g; g!= adjacent_faces.end(); g++)
	dual.edge_set.insert(edge_t(*f,*g));
    }
  }
  //  fprintf(stderr,"%d nodes, and %d edges in dual graph.\n",int(dual.N), int(dual.edge_set.size()));

  dual.update_auxiliaries();

  // If original graph was planar with 2D layout, there's a corresponding layout for the dual grap
  if(layout2d.size() == N){
    dual.layout2d = vector<coord2d>(Nfaces);
    for(unsigned int i=0;i<Nfaces;i++){
      face_t face(allfaces[i]);
      coord2d centre = 0;
      for(unsigned int j=0;j<face.size();j++) centre += layout2d[face[j]];
      dual.layout2d[i] = centre / face.size();
    }
  }

  // for(unsigned int u=0;u<Nfaces;u++){
  //   printf("%d: ",u);
  //   for(int i=0;i<dual.neighbours[u].size();i++)
  //     printf("%d ",dual.neighbours[u][i]);
  //   printf("\n");
  // }  
  return dual;
}

// NB: TODO: What happens, for example, if a triangle is comprised of three smaller triangles?
// This produces "phantom" faces! Fix and use the oriented version instead.
facemap_t PlanarGraph::compute_faces(unsigned int Nmax) const 
{
  facemap_t facemap;

  // TODO: This is a much better and faster method, but needs to be debugged.
  if(layout2d.size() == N) return compute_faces_oriented();

  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const node_t s = e->first, t = e->second;

    const vector<node_t>& ns(neighbours[t]);
    for(unsigned int i=0;i<ns.size();i++)
      if(ns[i] != s) {
	const node_t u = ns[i];

	face_t face(shortest_cycle(s,t,u,Nmax));  
	if(face.size() <= Nmax){
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
  facemap_t facemap;
  cout << "Computing faces using 2D orientation.\n";
  set<dedge_t> workset;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const node_t s = e->first, t = e->second;
    workset.insert(dedge_t(s,t));
    workset.insert(dedge_t(t,s));
  }

  // Outer face must exist and be ordered CW, rest of faces CCW
  assert(outer_face.size() > 0);

  cerr << "Outer face: " << outer_face << endl;

  coord2d centre(centre2d(layout2d));

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

      node_t w=0;
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
	    } // else {
	    //   fprintf(stderr,"\t[%d->%d already used.]\n",v,ns[i]);
	    // }
	  }
	}
      }
      e = dedge_t(v,w);
      workset.erase(e);

      if(e.second != face[0]) face.push_back(e.second);
    }
    //    cout << "face = " << face << endl;
    facemap[face.size()].insert(face);
  }
  return facemap;
}



vector<face_t> PlanarGraph::compute_faces_flat(unsigned int Nmax) const 
{
  vector<face_t> result;
  facemap_t facemap(compute_faces(Nmax));
  for(facemap_t::const_iterator fs(facemap.begin()); fs != facemap.end(); fs++)
    copy(fs->second.begin(),fs->second.end(),inserter(result,result.end()));

  return result;
}



vector<face_t> PlanarGraph::triangulation(int face_max) const
{
  vector<face_t> faces(compute_faces_flat(face_max));  
  return triangulation(faces);
}
  
vector<face_t> PlanarGraph::triangulation(const vector<face_t>& faces) const
{
  assert(layout2d.size() == N);
  vector<face_t> tris;
  
  
  for(size_t i=0;i<faces.size();i++){
    face_t f(faces[i]);

    for(size_t j=1;j<f.size()-1;j++){
      face_t t(3); 
      t[0] = f[0]; t[1] = f[j]; t[2] = f[j+1];

      coord2d c(t.centroid(layout2d));
      sort_ccw_point CCW(layout2d,c);
      
      sort(t.begin(),t.end(),CCW);
      if(i == 0) reverse(t.begin(), t.end()); // TODO: Show normals!

      tris.push_back(t);
    }
  }

  return tris;
}

// Finds the vertices belonging to the outer face in a symmetric planar
// layout centered at (0,0). Returns the face in CW order.
face_t PlanarGraph::find_outer_face() const 	
{
  assert(layout2d.size() == N);

  vector<double> radii(N);

  for(node_t u=0;u<N;u++)
    radii[u] = layout2d[u].norm();
  
  double rmax = *max_element(radii.begin(),radii.end());

  face_t outer_face;
  for(node_t u=0;u<N;u++)
    if(fabs(radii[u]-rmax) < 1e-10)
      outer_face.push_back(u);

  sort_ccw_point CCW(layout2d,centre2d(layout2d));
  sort(outer_face.begin(),outer_face.end(),CCW);
  reverse(outer_face.begin(),outer_face.end());  

  return outer_face;
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
  }
  s << "\n]";

  return s;
}
