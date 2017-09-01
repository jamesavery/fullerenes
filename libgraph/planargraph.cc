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
  vector<face_t> faces(compute_faces());

  for(int i=0;i<faces.size();i++) if(faces[i].size() != 3) return false;
  return true;
}

bool PlanarGraph::is_a_fullerene(bool verbose) const {
  if(!is_cubic()){
    if(verbose) fprintf(stdout,"Graph is not cubic.\n");
    return false;
  }

  vector<face_t> faces(compute_faces(6,true));
  int n_faces = faces.size();
  int n_edges = count_edges();
  
  const int E = 3*N/2;
  const int F = 2+E-N;

  if(E != n_edges){
    if(verbose) fprintf(stdout,"Graph is not planar cubic: wrong number of edges: %d != %d\n",n_edges,E);
    return false;
  }

  if(F != n_faces){
    if(verbose){
      fprintf(stdout,"Graph is not planar cubic: wrong number of faces: %d != %d\n",n_faces,F);
      cout << "faces = " << faces << ";\n";
    }
    return false;
  }

  int Np=0, Nh=0;
  for(const face_t &f: faces){
    if(f.size()==5) Np++;
    if(f.size()==6) Nh++;
  }
  
  if(Np != 12){
    if(verbose) fprintf(stdout,"Graph is not fullerene: wrong number of pentagons: %d != 12\n",Np);
    return false;
  }

  if(Nh != (F-12)){
    if(verbose) fprintf(stdout,"Graph is not fullerene: wrong number of hexagons: %d != %d\n",Nh,F-12);
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


// checks if the planar graph stays connected after removing v.  this function
// implies and relies on the condition that the graph has at most one face
// larger than a triangle.  If there is more than one larger face than a
// triangle, the function may return 'false', even though the correct answer is
// 'true'.
bool PlanarGraph::is_cut_vertex(const node_t v) const {
  assert(is_oriented); // we need oriented (sorted) neighbours of v (direction doesn't matter)
  const vector<node_t> &nv = neighbours[v];
  const int n_neighbours = nv.size();
  if(n_neighbours < 2) return false;

  int n_edges = 0;
  for(int i=0; i<n_neighbours; i++){
    const int v1=nv[i], v2=nv[(i+1)%n_neighbours];
        // and by counting this way we don't count edges between non-neighbours,
        // thus avoid the separating triangle problem
    if( edge_exists(edge_t(v1,v2)) ){
      n_edges++;
    }
  }
  // in a ring of n vertices where each except one adjacent face are triangles,
  // the induced graph is connected exactly when there are at least n-1
  // triangles
  return n_edges < n_neighbours-1;
}


PlanarGraph PlanarGraph::dual_graph(unsigned int Fmax, bool planar_layout, unordered_map<dedge_t,int> * FM) const
{
  if(is_oriented){

    // Each directed edge uniquely identifies a face
    vector<dedge_t>            face_reps = compute_face_representations(Fmax);
    unordered_map<dedge_t,int> face_numbers(face_reps.size());
    for(int i=0;i<face_reps.size();i++) face_numbers[face_reps[i]] = i;

    // if FM is not a null pointer, copy face numbers back
    if(FM) *FM=face_numbers;

    // cerr << "face_reps  = " << face_reps << ";\n";
    // cerr << "face_reps' = " << get_keys(face_numbers) << ";\n";
    // cerr << "face_nums  = " << get_values(face_numbers) << ";\n";
    //    cerr << "faces      = " << compute_faces_oriented(Fmax) << ";\n";
    
    PlanarGraph dual(face_numbers.size());
    
    for(const auto &ei: face_numbers){
      // e_f is minimal directed edge representation of face f, i_f is its face number.
      dedge_t e_f = ei.first; node_t i_f = ei.second;
      // cerr << "Processing face " << i_f << ": " << e_f << " -> " << get_face_oriented(e_f,Fmax) << ";\n";
      
      // Now iterate along face f's directed edges in CCW order:
      // This visits each face neighbour in CCW order.
      node_t u=e_f.first, v=e_f.second, w=-1, i=0;
      do {
        // e_g is MDE-representation of opposite face along edge e_f
        dedge_t e_g = get_face_representation({v,u},Fmax);
        // cerr << "Processing face: " << dedge_t{u,v} << " -> " << get_face_oriented({u,v},Fmax)<<" "
        //      << " opposite " << e_g << " -> " << get_face_oriented(e_g,Fmax) << ";\n";
        dual.neighbours[i_f].push_back(face_numbers[e_g]);

        w = prev(v,u); u = v; v = w; // CCW node neighbour order + CCW face order
        assert(++i <= Fmax);        // Face larger than Fmax or corrupted graph
      } while (u != e_f.first);
    }
    assert(dual.is_consistently_oriented());
    dual.is_oriented = true;
    
    if(planar_layout && layout2d.size() == N){
      dual.layout2d = vector<coord2d>(face_numbers.size());

      for(const auto &ei: face_numbers){
        dedge_t e_f = ei.first; node_t i_f = ei.second;
        face_t f = get_face_oriented(e_f,Fmax);
        
        dual.layout2d[i_f] = f.centroid(layout2d);
      }
    }

    return dual;
    // Proper implementation of general oriented planar graph dual ends
  } else {
    // TODO: ********** Get rid of all this junk below! ******************
    IDCounter<face_t> face_numbers;
    PlanarGraph dual;
    set<edge_t> edge_set = undirected_edges(); // TODO: In new planargraph, this is unnecessary
    int Nfaces = edge_set.size()-N+2;
    dual.N = Nfaces;
    dual.neighbours.resize(Nfaces);

    //  cerr << "dual_graph(" << Fmax << ")\n";
    vector<face_t> allfaces = compute_faces(Fmax,planar_layout);

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
    // (but it is not planar, because the outer face is placed in (0,0))
    if(planar_layout && layout2d.size() == N){
      //    cerr << "dual_graph::compute layout.\n";
      dual.layout2d = vector<coord2d>(allfaces.size());

      for(int i=0;i<allfaces.size();i++)
        dual.layout2d[i] = allfaces[i].centroid(layout2d);
    }
    return dual;
  }    
}


Graph PlanarGraph::leapfrog_dual() const
{
  assert(is_oriented);
  vector<face_t> faces = compute_faces_oriented();

  Graph lf(N+faces.size(),true);

  // Start with all the existing nodes
  for(node_t u=0;u<N;u++) lf.neighbours[u] = neighbours[u];

  // Now connect new face-center nodes in oriented order
  for(int i=0;i<faces.size();i++){
    const face_t &f  = faces[i];
    node_t c = N+i;                // Face-center node

    cerr << "new node " << c << " at face " << f << "\n";
    
    for(int j=0;j<f.size();j++){
      node_t u = f[j], v = f[(j+1)%f.size()];

      lf.insert_edge(dedge_t{v,c},u,-1);
    }
  }

  return lf;
}


vector<face_t> PlanarGraph::compute_faces(unsigned int Nmax, bool planar_layout) const
{
  // TODO: This should supercede using the planar embedding for orientation
  if(is_oriented) return compute_faces_oriented();
  
  // TODO: Clean up.
  if(planar_layout && layout2d.size() == N) return compute_faces_layout_oriented();

  
  // TODO: This should never be used
  cerr << " Non-oriented face computation (loop search). This is not reliable!\n";
  // abort();
  cerr << "This shouldn't happen but we'll accept it for now." << endl;
  set<edge_t> edge_set = undirected_edges();
  set<face_t> faces;
  for(const edge_t &e: edge_set){
    const node_t s = e.first, t = e.second;
    
    const vector<node_t>& nt(neighbours[t]);
    
    for(unsigned int i=0;i<nt.size();i++)
      if(nt[i] != s) {
        const node_t u = nt[i];

        face_t face(shortest_cycle(s,t,u,Nmax));
        //	cerr << face << endl;
        if(face.size() > 0 && face.size() <= Nmax)
          faces.insert(face.normalized());
        } //else {
          //	  fprintf(stderr,"Erroneous face starting at (%d -> %d -> %d) found: ",s,t,u);
          //	  cerr << face << endl;
        //	}
  }


  // // Make sure that outer face is at position 0
  // if(planar_layout){
  //   if(outer_face.size() < 3)
  //     outer_face = find_outer_face();

  //   const set<node_t> of(outer_face.begin(),outer_face.end());
  //   for(int i=0;i<faces.size();i++){
  //     const face_t &f(faces[i]);
  //     const set<node_t> sf(f.begin(),f.end());

  //     if(of==sf){ // swap faces[i] with faces[0]
  //      faces[i] = faces[0];
  //      faces[0] = outer_face;
  //     }
  //   }
  // } else outer_face = face_t(faces[0]);
  
  vector<face_t> face_vector(faces.begin(),faces.end());

  return face_vector;
}

 
face_t PlanarGraph::get_face_layout_oriented(node_t s, node_t t) const
{
  face_t face;
  face.push_back(s);
  face.push_back(t);

  node_t u = s, v = t;
  //    printf("%d->%d\n",e.first,e.second);
  while(v != s){
    const vector<node_t>& ns(neighbours[v]);


    coord2d vu = layout2d[u]-layout2d[v];
    double angle_max = -M_PI;

    node_t w=-1;
    for(unsigned int i=0;i<ns.size();i++) {
      //        printf("%d : %d (%d->%d) angle %g\n",i,ns[i],u,v,vu.line_angle(layout[ns[i]]-layout[v]));
      if(ns[i] != u) { // Find and use first unvisited edge in order of angle to u->v
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



vector<face_t> PlanarGraph::compute_faces_layout_oriented() const
{
  assert(layout2d.size() == N);
  //  cout << "Computing faces using 2D orientation." << endl;
  set<dedge_t> workset;
  set<edge_t> edge_set = undirected_edges();

  vector<face_t> faces;
  set<face_t> face_set;
  
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const node_t s = e->first, t = e->second;
    workset.insert(dedge_t(s,t));
    workset.insert(dedge_t(t,s));
  }

  // If layout is planar, outer face must exist and be ordered CW,
  // rest of faces CCW. If layout is spherical / periodic, all faces
  // should be ordered CCW.

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
    faces.push_back(outer_face);
    // Add outer face to output, remove directed edges from work set
    for(unsigned int i=0;i<outer_face.size();i++){
      const node_t u = outer_face[i], v = outer_face[(i+1)%outer_face.size()];
      //    printf("Removing directed edge (%d,%d)\n",u,v);
      workset.erase(dedge_t(u,v));
    }

  // Now visit every other edge once in each direction.
  while(!workset.empty()){
    dedge_t e = *workset.begin();
    face_t face(get_face_layout_oriented(e.first,e.second));
    face_set.insert(face);

    for(int i=0;i<face.size();i++)
      workset.erase(dedge_t(face[i],face[(i+1)%face.size()]));
  }

  copy(face_set.begin(), face_set.end(), std::back_inserter(faces));
  
  return faces;
}

// sort neighbour list CW
void PlanarGraph::orient_neighbours()
{
  if(is_oriented) return;
  
  assert(layout2d.size() == N);
  for(node_t u=0;u<N;u++){
    sort_ccw_point CCW(layout2d,layout2d[u]);
    sort(neighbours[u].begin(),neighbours[u].end(),CCW);
    reverse(neighbours[u].begin(),neighbours[u].end());
  }
  is_oriented = true;
}


vector<tri_t> PlanarGraph::triangulation(int face_max) const
{
  vector<face_t> faces(compute_faces(face_max));
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

    if(f.size() > 3)
      for(int j=0;j<f.size();j++)
        tris.push_back({f[j],v_new,f[(j+1)%f.size()]});
    else
      tris.push_back({f[0],f[1],f[2]});
  }

  return tris;                        // TODO: Make sure triangulation is oriented.
  //return orient_triangulation(tris);
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
      cerr << "tris = " << tris << "+1;\n";
      cerr << "g    = " << *this << ";\n";
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
//                       COMBINATORIAL PROPERTIES
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
  vector<face_t> faces(compute_faces());
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

// TODO: Where does this belong?
// Assumes file is at position of a graph start
Graph PlanarGraph::read_hog_planarcode(FILE *planarcode_file)
{
  // Read the number N of vertices per graph.
  int number_length=1, N=0;
  fread(reinterpret_cast<unsigned char*>(&N), 1, 1, planarcode_file);
  if(N == 0){
    fread(reinterpret_cast<unsigned char*>(&N), 2, 1, planarcode_file);
    number_length=2;
  }
  
  Graph g(N,true);
  for(node_t u=0; u<N && !feof(planarcode_file); ++u){
    int v=0;
    do{
      int n_read = fread(reinterpret_cast<char*>(&v), number_length, 1, planarcode_file);
      if(n_read != 1 && !feof(planarcode_file)){
        perror("Error reading HoG PlanarCode file: ");
        abort();
      }
      if(v!=0) g.neighbours[u].push_back(v-1); // In oriented order
    } while(v!=0 && !feof(planarcode_file));
  }
  // Check graph
  for(node_t u=0;u<N;u++){
    for(auto v: g.neighbours[u]){
      bool found_vu = false;
      
      for(node_t w: g.neighbours[v])
        if(w == u) found_vu = true;
      if(!found_vu){
        fprintf(stderr,"Graph is not symmetric: (u,v) = (%d,%d) has\n",u,v);
        cerr << "neighbours["<<u<<"] = " << g.neighbours[u] <<";\n";
        cerr << "neighbours["<<v<<"] = " << g.neighbours[v] <<";\n";
        abort();
      }
    }
  }

  return g;
}

// TODO: Read only a range
vector<Graph> PlanarGraph::read_hog_planarcodes(FILE *planarcode_file) {
  const int header_size = 15;
  vector<Graph> graph_list;

  //the actual parsing of the selected graph:
  //go to the beginning of the selected graph
  fseek(planarcode_file,  header_size, SEEK_SET);

  //  int i = 1;
  while(!feof(planarcode_file)){
    //    cerr << "Reading graph " << (i++) << ".\n";
    Graph g = read_hog_planarcode(planarcode_file);
    //    cerr << "Got graph on " << g.N << " vertices.\n";
    if(g.N != 0){
      graph_list.push_back(g);
    }
  }
    
  return graph_list;
}

 
 
// In an oriented planar graph, the directed edge starting in the smallest node
// is a unique representation of the face.
dedge_t PlanarGraph::get_face_representation(dedge_t e, int Fmax) const
{
  assert(is_oriented);


  int i=0;  
  dedge_t e_min = e;
  node_t u = e.first, v = e.second;

  while(v!=e.first){
    node_t w = prev(v,u);        // Previous neighbour in v following u in CCW (corner u-v-w in face)
    u=v; v=w; 

    if(u<e_min.first) e_min = {u,v};
    
    assert(w != -1);
    assert(++i<=Fmax); // Fmax is a back-stop to avoid infinite loops in a corrupted graph
  }
  return e_min;
}

// In an oriented planar graph, the directed edge starting in the smallest node
// is a unique representation of the face.
vector<dedge_t> PlanarGraph::compute_face_representations(int Fmax) const
{
  assert(is_oriented);

  unordered_set<dedge_t> faces(2*count_edges());
  
  for(node_t u=0;u<N;u++)
    for(node_t v: neighbours[u]){
      // For each directed edge, find the representative edge of the specified face
      // and assign an identifier
      faces.insert(get_face_representation({u,v},Fmax));
    }

  return vector<dedge_t>(faces.begin(),faces.end());
}


face_t PlanarGraph::get_face_oriented(const dedge_t &e, int Fmax) const
{
  assert(is_oriented);

  int i=0;
  node_t u = e.first, v=e.second;
  face_t f = vector<int>{{u}};
  
  while(v!=e.first){
    node_t w = prev(v,u);        // Previous neighbour to u in v defines corner u-v-w in face

    f.push_back(v);
    u=v; v=w; i++;
    assert(w != -1);
    assert(i<=Fmax);                // Fmax is a back-stop to avoid infinite loops in a corrupted graph
  }
  return f;
}

vector<face_t> PlanarGraph::compute_faces_oriented(int Fmax) const
{
  vector<dedge_t> face_representations = compute_face_representations(Fmax);

  vector<face_t> faces(face_representations.size());
  for(int i=0;i<face_representations.size();i++) faces[i] = get_face_oriented(face_representations[i],Fmax);

  // cerr << "facereps = " << face_representations << ";\n"
  //      << "faces    = " << faces << ";\n";
  
  return faces;
}
