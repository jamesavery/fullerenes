#include "planargraph.hh"
#include <queue>
#include <list>
using namespace std;

bool PlanarGraph::is_cubic() const {
  for(node_t u=0;u<N;u++)
    if(neighbours[u].size() != 3)
      return false;
  return true;
}

bool PlanarGraph::is_triangulation() const { // NB: A bit expensive
  vector<face_t> faces(compute_faces_flat(INT_MAX, true));

  for(int i=0;i<faces.size();i++) if(faces[i].size() != 3) return false;
  return true;
}

bool PlanarGraph::is_a_fullerene() const {
  if(!is_cubic()){
    fprintf(stdout,"Graph is not cubic.\n"); 
    return false;
  }
    
  facemap_t faces(compute_faces(6,true));
  int n_faces = 0;
  for(facemap_t::const_iterator f(faces.begin()); f!=faces.end();f++)
    n_faces += f->second.size();

  const int E = 3*N/2;
  const int F = 2+E-N;

  set<edge_t> edge_set = undirected_edges(); // TODO: Do with neighbours - this is a bit slow.
  if(E != edge_set.size()){
    fprintf(stdout,"Graph is not planar cubic: wrong number of edges: %d != %d\n",int(edge_set.size()),E);
    return false;
  }

  if(F != n_faces){
    fprintf(stdout,"Graph is not planar cubic: wrong number of faces: %d != %d\n",n_faces,F);
    cout << "faces = " << get_values(faces) << ";\n";
    return false;
  }

  if(faces[5].size() != 12){
    fprintf(stdout,"Graph is not fullerene: wrong number of pentagons: %d != 12\n",int(faces[5].size()));
    return false;
  }

  if(faces[6].size() != (F-12)){
    fprintf(stdout,"Graph is not fullerene: wrong number of hexagons: %d != %d\n",int(faces[6].size()),F-12);
    return false;
  }

  return true;
}

// the following is a naive approach that iterates over all pairs of edges
// for some purposes it would be sufficient to ensure that each face intersects itself an even number of times (while figures of eight are problematic)
bool PlanarGraph::layout_is_crossingfree() const 
{
  assert(layout2d.size() == N);
  set<edge_t> es = undirected_edges(); // TODO: In new planargraph, this is unnecessary
  for (set<edge_t>::iterator e1(es.begin()); e1!=es.end(); e1++){
    for (set<edge_t>::iterator e2(e1); e2!=es.end(); e2++){
      if (e1->first == e2->first || e1->second == e2->first || e1->first == e2->second || e1->second == e2->second) continue; // equal edges and edges that share a vertex
      const double e1ax = layout2d[e1->first].first,
                   e1ay = layout2d[e1->first].second,
                   e1bx = layout2d[e1->second].first,
                   e1by = layout2d[e1->second].second,
                   e2ax = layout2d[e2->first].first,
                   e2ay = layout2d[e2->first].second,
                   e2bx = layout2d[e2->second].first,
                   e2by = layout2d[e2->second].second;
      const double a1 = (e1ay - e1by)/(e1ax - e1bx);
      const double b1 = e1ay - a1 * e1ax;
      if ((e2ay > a1*e2ax+b1 && e2by > a1*e2bx+b1) || (e2ay < a1*e2ax+b1 && e2by < a1*e2bx+b1)) continue; // both points of the second edge lie on the same side of the first edge
      const double a2 = (e2ay - e2by)/(e2ax - e2bx);
      const double b2 = e2ay - a2 * e2ax;
      if ((e1ay > a2*e1ax+b2 && e1by > a2*e1bx+b2) || (e1ay < a2*e1ax+b2 && e1by < a2*e1bx+b2)) continue; // both points of the first edge lie on the same side of the second edge
      cerr << "edges " << *e1 << " and " << *e2 << " intersect." << endl;
      return false;
    }
  }
  return true;
}




PlanarGraph PlanarGraph::dual_graph(unsigned int Fmax, bool planar_layout) const {
  // TODO: Simplify
  PlanarGraph dual;
  set<edge_t> edge_set = undirected_edges(); // TODO: In new planargraph, this is unnecessary
  unsigned int Nfaces = edge_set.size()-N+2;
  dual.N = Nfaces;
  dual.neighbours.resize(Nfaces);

  //  cerr << "dual_graph(" << Fmax << ")\n";
  const vector<face_t> allfaces(compute_faces_flat(Fmax,planar_layout));

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
  set<edge_t> dual_edges;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const set<int>& adjacent_faces(facenodes[*e]);
    for(set<int>::const_iterator f(adjacent_faces.begin()); f!= adjacent_faces.end(); f++){
      set<int>::const_iterator g(f);
      for(++g; g!= adjacent_faces.end(); g++)
        dual_edges.insert(edge_t(*f,*g));
    }
  }
  //fprintf(stderr,"%d nodes, and %d edges in dual graph.\n",int(dual.N), int(dual.edge_set.size()));

  dual = Graph(dual_edges);

  // If original graph was planar with 2D layout, there's a corresponding layout for the dual graph
  // (but it is not planar -- might not want to use this!)
  if(planar_layout && layout2d.size() == N){
    //    cerr << "dual_graph::compute layout.\n";
    dual.layout2d = vector<coord2d>(Nfaces);

    for(int i=0;i<Nfaces;i++)
      dual.layout2d[i] = allfaces[i].centroid(layout2d);
  }
  return dual;
}



// NB: TODO: What happens, for example, if a triangle is comprised of three smaller triangles?
// This produces "phantom" faces! Fix and use the oriented version instead.
facemap_t PlanarGraph::compute_faces(unsigned int Nmax, bool planar_layout) const 
{
  set<edge_t> edge_set = undirected_edges();

  facemap_t facemap;
  // TODO: This is a much better and faster method, but requires a planar layout
  if(planar_layout && layout2d.size() == N) return compute_faces_oriented();

  cerr << " Non-oriented face computation (loop search). This is not reliable!\n";
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

face_t PlanarGraph::get_face_oriented(node_t s, node_t t) const 
{
  face_t face;
  face.push_back(s);
  face.push_back(t);
  
  node_t u = s, v = t;
  //    printf("%d->%d\n",e.first,e.second);
  while(v != s){
    const vector<node_t>& ns(neighbours[v]);

    
    //coord2d vu = coord2d::displacement(layout2d[u],layout2d[v],layout_is_spherical);
    coord2d vu = layout2d[u]-layout2d[v];
    double angle_max = -M_PI;

    node_t w=-1;
    for(unsigned int i=0;i<ns.size();i++) {
      //	printf("%d : %d (%d->%d) angle %g\n",i,ns[i],u,v,vu.line_angle(layout[ns[i]]-layout[v]));
      if(ns[i] != u) { // Find and use first unvisited edge in order of angle to u->v
        //	coord2d vw = coord2d::displacement(layout2d[ns[i]],layout2d[v],layout_is_spherical);
        coord2d vw = layout2d[ns[i]]-layout2d[v];
        double angle = vu.line_angle(vw);
        
        if(angle>= angle_max){
          angle_max = angle;
          w = ns[i];
        } 
      } 
    }
    if(w == -1) abort(); // There is no face!

    u = v; v = w;
      
    if(w != s) face.push_back(w);
  }
  return face;
}

facemap_t PlanarGraph::compute_faces_oriented() const 
{
  assert(layout2d.size() == N);
  facemap_t facemap;
  int faces_found =0;
  //  cout << "Computing faces using 2D orientation." << endl;
  set<dedge_t> workset;
  set<edge_t> edge_set = undirected_edges();

  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const node_t s = e->first, t = e->second;
    workset.insert(dedge_t(s,t));
    workset.insert(dedge_t(t,s));
  }

  // If layout is planar, outer face must exist and be ordered CW,
  // rest of faces CCW. If layout is spherical / periodic, all faces
  // should be ordered CCW.
  if(!layout_is_spherical){
    if(outer_face.size() < 3)
    outer_face = find_outer_face();

    if(outer_face.size() < 3){
      cerr << "Invalid outer face: " << outer_face << endl;
      assert(outer_face.size() < 3);
    }

    for(node_t u=0;u<N;u++)
      if(!outer_face.contains(u) && !outer_face.point_inside(layout2d,u)){
        cerr << "Point " << u << "/" << layout2d[u] << " is outside outer face " << outer_face << endl;
        for(int i=0;i<outer_face.size();i++) cerr << "\t" << layout2d[outer_face[i]] << endl;
        cerr << "Winding number: " << outer_face.winding_number(layout2d,u) << endl;
        abort();
      }
    //    cout << "compute_faces_oriented: Outer face "<<outer_face<<" is OK: All vertices are inside face." << endl;
    facemap[outer_face.size()].insert(outer_face);
    faces_found++;
    // Add outer face to output, remove directed edges from work set
    for(unsigned int i=0;i<outer_face.size();i++){
      const node_t u = outer_face[i], v = outer_face[(i+1)%outer_face.size()];
      //    printf("Removing directed edge (%d,%d)\n",u,v);
      workset.erase(dedge_t(u,v));
    }
  }

  // Now visit every other edge once in each direction.
  while(!workset.empty()){
    dedge_t e = *workset.begin(); 
    face_t face(get_face_oriented(e.first,e.second));
    facemap[face.size()].insert(face);
    faces_found++;

    for(int i=0;i<face.size();i++)
      workset.erase(dedge_t(face[i],face[(i+1)%face.size()]));
  }
  assert(faces_found == N/2+2 || faces_found == (N-2)*2); // assuming the graph is either cubic or a triangulation // FIXME remove when things are more stable
  return facemap;
}

// sort neighbour list CW
void PlanarGraph::orient_neighbours() 
{
  assert(layout2d.size() == N);
  for(node_t u=0;u<N;u++){
    sort_ccw_point CCW(layout2d,layout2d[u]);
    sort(neighbours[u].begin(),neighbours[u].end(),CCW);
    reverse(neighbours[u].begin(),neighbours[u].end());
  }
}

vector<face_t> PlanarGraph::compute_faces_flat(unsigned int Nmax, bool planar_layout) const 
{
  vector<face_t> faces;
  facemap_t facemap(compute_faces(Nmax,planar_layout));

  for(facemap_t::const_iterator fs(facemap.begin()); fs != facemap.end(); fs++)
    copy(fs->second.begin(),fs->second.end(),inserter(faces,faces.end()));

  // Check that faces are orientable: Every edge must appear in two faces
  map<edge_t,int> edgecount;
  for(int i=0;i<faces.size();i++)
    for(int j=0;j<faces[i].size();j++)
      edgecount[edge_t(faces[i][j],faces[i][(j+1)%faces[i].size()])]++;

  for(map<edge_t,int>::const_iterator e(edgecount.begin()); e!=edgecount.end();e++)
    if(e->second != 2){
      cerr << "compute_faces_flat: Graph not orientable - edge "<< e->first << " appears in " << e->second <<" faces, not two.\n";
      cerr << "faces = {"; for(int i=0;i<faces.size();i++) cerr << faces[i] << (i+1<faces.size()?", ":"};\n");
      cerr << "G = " << *this << ";\n";

      abort();
    }


  // Make sure that outer face is at position 0
  if(planar_layout){
    if(outer_face.size() < 3)
      outer_face = find_outer_face();

    const set<node_t> of(outer_face.begin(),outer_face.end());
    for(int i=0;i<faces.size();i++){
      const face_t &f(faces[i]);
      const set<node_t> sf(f.begin(),f.end());

      if(of==sf){ // swap faces[i] with faces[0]
       faces[i] = faces[0];
       faces[0] = outer_face;
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
    //    cerr << "centroid_triangulation: Faces already form a triangulation.\n";
    vector<tri_t> tris(faces.begin(),faces.end());
    return orient_triangulation(tris);
  } else {
    //    cerr << "centroid_triangulation: Not a triangulation. Building centroid triangulation!\n";
    // cerr << "Original faces:\n";
    // cerr << "faces = {"; for(int i=0;i<faces.size();i++) cerr << faces[i] << (i+1<faces.size()?", ":"};\n");
    // cerr << "layout = {"; for(int i=0;i<layout2d.size();i++) cerr << layout2d[i] << (i+1<layout2d.size()?", ":"};\n");
    // cerr << "G = " << *this << ";\n";
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
    //cerr << "PlanarGraph::triangulation: Faces already form a triangulation.\n";
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
  // Check that triangles are orientable: Every edge must appear in two faces
  map<edge_t,int> edgecount;
  for(int i=0;i<tris.size();i++)
    for(int j=0;j<3;j++){
      edgecount[edge_t(tris[i][j],tris[i][(j+1)%3])]++;
      if(edgecount[edge_t(tris[i][j],tris[i][(j+1)%3])]>2)
        cerr << tris[i] << " bad!\n";
    }
  for(map<edge_t,int>::const_iterator e(edgecount.begin()); e!=edgecount.end();e++)
    if(e->second != 2){
      cerr << "Triangulation not orientable: Edge "<< e->first << " appears in " << e->second <<" tris, not two.\n";
      abort();
    }

  // Now, pick an orientation for triangle 0. We choose the one it
  // already has. This determines the orientation of the remaining triangles!
  map<dedge_t,bool> done;
  for(int i=0;i<3;i++){
    done[dedge_t(tris[0][i],tris[0][(i+1)%3])] = true;
  }

  queue<int> workset; 
  for(int i=1;i<tris.size();i++) workset.push(i);

  while(!workset.empty()){
    int i = workset.front(); workset.pop();
    tri_t& t(tris[i]);


    // Is this triangle connected to any already processed triangle?
    bool seen = false, rev_seen = false;
    for(int j=0;j<3;j++){  seen |= done[dedge_t(t[j],t[(j+1)%3])]; rev_seen |= done[dedge_t(t[(j+1)%3],t[j])]; }
    if(!seen && !rev_seen) {
      workset.push(i);
      continue;
    }

    if(seen){
      node_t u = t[2]; t[2] = t[1]; t[1] = u;
    }
    
    done[dedge_t(t[0],t[1])] = true;
    done[dedge_t(t[1],t[2])] = true;
    done[dedge_t(t[2],t[0])] = true;
  }
  // Check consistency of orientation. It is consistent if and only if
  // each edge has been used exactly once in each direction.
  bool consistent = true;
  set<edge_t> edge_set = undirected_edges();

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
  sort(outer_face.begin(),outer_face.end(),CCW); // sort CCW
  reverse(outer_face.begin(),outer_face.end()); // reverse to get CW

  //  cout << "Found outer face: " << outer_face << endl;
  return outer_face;
}

vector<double> PlanarGraph::edge_lengths() const 
{
  assert(layout2d.size() == N);
  set<edge_t> edge_set = undirected_edges();

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
  set<edge_t> edge_set = g.undirected_edges();

  s << "Graph[Range["<<g.N<<"],\n\tUndirectedEdge@@#&/@{";
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end(); ){    
    s << "{" << (e->first+1) << "," << (e->second+1) << "}";
    if(++e != edge_set.end())
      s << ", ";
    else
      s << "}";
  }

  if(g.layout2d.size() == g.N){
    s << ",\n\tVertexCoordinates->{";
    for(unsigned int i=0;i<g.N;i++){
      coord2d xy(g.layout2d[i]);
      s << xy << (i+1<g.N?", ":"}");
    }
  } // else { fprintf(stderr,"No layout, man!\n"); }
  s << "\n]";

  return s;
}


// **********************************************************************
//		       COMBINATORIAL PROPERTIES
// **********************************************************************

void perfmatch_dfs(map<dedge_t,int>& faceEdge, const vector<face_t>& faces, 
		   map<dedge_t,int>& matrix, vector<bool>& faceSum, vector<bool>& visited, const dedge_t& e) 
{
  int frev = faceEdge[reverse(e)];
  if(visited[frev]) return;
  visited[frev] = true;

  const face_t &f(faces[frev]);
  for(int i=0;i<f.size();i++)
    perfmatch_dfs(faceEdge,faces,matrix,faceSum,visited,dedge_t(f[i],f[(i+1)%f.size()]));

  // NB: How to handle outer face?
  if(!faceSum[frev]) { //not odd sum of CW edges
    int fe = faceEdge[e];
    faceSum[frev] = !faceSum[frev];
    faceSum[fe] = !faceSum[fe];
    matrix[e] *= -1;
    matrix[reverse(e)] *= -1;
  }

}

#ifdef HAS_LAPACK
#ifdef HAS_MKL
#include <mkl.h>
#else
extern "C" void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
#endif

double lu_det(const vector<double> &A, int N)
{
  int info = 0;
  double *result = new double[N*N];
  int    *ipiv   = new int[N];
  double prod = 1.0;
  memcpy(result,&A[0],N*N*sizeof(double));
  dgetrf_(&N,&N, result, &N, ipiv, &info);
  {
    int i;
    for(i=0;i<N;i++) prod *= result[(N+1)*i];
  }
  free(result);
  free(ipiv);
  return fabs(prod);
}


size_t PlanarGraph::count_perfect_matchings() const 
{
  map<dedge_t,int> faceEdge;
  vector<face_t> faces(compute_faces_flat(max_degree(), true));
  vector<bool> faceSum(faces.size()), visited(faces.size());  

  map<dedge_t,int> A;
  set<edge_t> edge_set = undirected_edges();
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end(); e++){
    A[*e] = 1;
    A[reverse(*e)] = -1;
  }
  
  for(int i=0;i<faces.size();i++){
    const face_t &f(faces[i]);
    for(int j=0;j<f.size();j++){
      const dedge_t e(f[j],f[(j+1)%f.size()]);
      faceEdge[e] = i;
      if(A[e] == 1) faceSum[i] = !faceSum[i];
    }
  }

  perfmatch_dfs(faceEdge,faces,A,faceSum,visited,*edge_set.begin());

  vector<double> Af(N*N);
  for(map<dedge_t,int>::const_iterator a(A.begin()); a!=A.end(); a++)
    Af[a->first.first*N+a->first.second] = a->second;

  return round(sqrtl(fabs(lu_det(Af,N))));
}
#else
size_t PlanarGraph::count_perfect_matchings() const 
{
  cerr << "count_perfect_matchings() requires LAPACK.\n";
  return 0;
}
#endif


vector<coord3d> PlanarGraph::zero_order_geometry(double scalerad) const
{
  assert(layout2d.size() == N);
  vector<coord2d> angles(spherical_projection());

  // Spherical projection
  vector<coord3d> coordinates(N);
  for(int i=0;i<N;i++){
    double theta = angles[i].first, phi = angles[i].second;
    double x = cos(theta)*sin(phi), y = sin(theta)*sin(phi), z = cos(phi);
    coordinates[i] = coord3d(x,y,z);
  }

  // Move to centroid
  coord3d cm;
  for(node_t u=0;u<N;u++) cm += coordinates[u];
  cm /= double(N);
  coordinates -= cm;

  // Scale spherical projection
  double Ravg = 0;
  for(node_t u=0;u<N;u++)
    for(int i=0;i<3;i++) Ravg += (coordinates[u]-coordinates[neighbours[u][i]]).norm();
  Ravg /= (3.0*N);
  
  coordinates *= scalerad*1.5/Ravg;

  return coordinates;
}
