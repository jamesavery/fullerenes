#include <queue>

#include "fullerenes/spiral.hh"
#include "fullerenes/planargraph.hh"
#include "fullerenes/layout2d.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/cubicgraph.hh"

using namespace std;

PlanarGraph::PlanarGraph(const spiral_nomenclature &fsn){
  switch(fsn.construction_scheme){
    case spiral_nomenclature::CS_NONE:
      // cerr << "none" << endl;
      assert(false);
      break;
    case spiral_nomenclature::CUBIC:
      // cerr << "CUBIC" << endl;
      *this = CubicGraph(fsn);
      break;
    case spiral_nomenclature::TRIANGULATION:
      // cerr << "TRIANGULATION" << endl;
      *this = Triangulation(fsn);
      break;
    case spiral_nomenclature::LEAPFROG:
      // cerr << "LEAPFROG" << endl;
      Triangulation T(fsn);
      *this = T.inverse_leapfrog_dual();
      break;
  }  
}

// Every polyhedral graph G can be represented by a triangulation.
//  1. If G is a triangulation, it is G
//  2. If G is cubic, it is its dual
//  3. If G is non-cubic and non-triangulation, it is G's leapfrog dual
PlanarGraph PlanarGraph::enveloping_triangulation(construction_scheme_t &scheme) const
{
  if(is_triangulation()){
    scheme = spiral_nomenclature::TRIANGULATION;
    return *this;
  } else if(is_cubic()){
    scheme = spiral_nomenclature::CUBIC;
    return dual_graph();
  } else {
    scheme = spiral_nomenclature::LEAPFROG;
    return leapfrog_dual();
  }
}

PlanarGraph PlanarGraph::enveloping_triangulation(const construction_scheme_t &scheme) const
{
  switch(scheme){
  case spiral_nomenclature::TRIANGULATION:
    return *this;
  case spiral_nomenclature::CUBIC:
    return dual_graph();
  case spiral_nomenclature::LEAPFROG:
  default:
    return leapfrog_dual();
  }
}

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

  vector<face_t> faces(compute_faces_oriented());
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
    if(verbose){
      fprintf(stdout,"Graph is not fullerene: wrong number of pentagons: %d != 12\n",Np);
      cout << "faces = " << faces << "\n";
    }
    return false;
  }

  if(Nh != (F-12)){
    if(verbose) fprintf(stdout,"Graph is not fullerene: wrong number of hexagons: %d != %d\n",Nh,F-12);
    return false;
  }

  return true;
}

// layout_is_crossingfree moved to layout2d.cc


// checks if the planar graph stays connected after removing v.  this function
// implies and relies on the condition that the graph has at most one face
// larger than a triangle.  If there is more than one larger face than a
// triangle, the function may return 'false', even though the correct answer is
// 'true'.
bool PlanarGraph::is_cut_vertex(const node_t v) const {
  // Requires oriented (sorted) neighbours of v (direction doesn't matter)
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



PlanarGraph PlanarGraph::dual_graph(unsigned int Fmax) const
{
  // Each directed edge uniquely identifies a face
  vector<arc_t>            face_reps = compute_face_representations(Fmax);
  unordered_map<arc_t,int> face_numbers(face_reps.size());
  for(int i=0;i<face_reps.size();i++) face_numbers[face_reps[i]] = i;

  PlanarGraph dual(face_numbers.size());

  for(const auto &ei: face_numbers){
    auto [e_f,i_f] = ei;

    node_t u=e_f.first, v=e_f.second, w=-1, i=0;
    do {
      arc_t e_g = get_face_representation({v,u},Fmax);
      dual.neighbours[i_f].push_back(face_numbers[e_g]);

      w = prev(v,u); u = v; v = w;
      assert(++i <= Fmax);
    } while (u != e_f.first);
  }
  assert(dual.is_consistently_oriented());

  return dual;
}

// the dual of the LF, ie a Triangulation is returned
PlanarGraph PlanarGraph::leapfrog_dual() const
{
  vector<face_t> faces = compute_faces_oriented();
  size_t Nf = faces.size();

  PlanarGraph lf(Graph(N+Nf));

  // Start with all the existing nodes
  for(node_t u=0;u<N;u++) lf.neighbours[u] = neighbours[u];

  // Now connect new face-center nodes in oriented order
  for(int i=0;i<Nf;i++){
    const face_t &f  = faces[i];
    node_t c = N+i;                // Face-center node

    for(int j=0;j<f.size();j++){
      node_t u = f[j], v = f[(j+1)%f.size()];

      lf.insert_edge(arc_t{v,c},u,-1);
    }
  }

  return lf;
}


vector<face_t> PlanarGraph::compute_faces(unsigned int Fmax) const {
  return compute_faces_oriented(Fmax);
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
  map<arc_t,bool> done;
  for(int i=0;i<3;i++){
    done[arc_t(tris[0][i],tris[0][(i+1)%3])] = true;
  }

  queue<int> workset;
  for(int i=1;i<tris.size();i++) workset.push(i);

  while(!workset.empty()){
    int i = workset.front(); workset.pop();
    tri_t& t(tris[i]);


    // Is this triangle connected to any already processed triangle?
    bool seen = false, rev_seen = false;
    for(int j=0;j<3;j++){  seen |= done[arc_t(t[j],t[(j+1)%3])]; rev_seen |= done[arc_t(t[(j+1)%3],t[j])]; }
    if(!seen && !rev_seen) {
      workset.push(i);
      continue;
    }

    if(seen){
      node_t u = t[2]; t[2] = t[1]; t[1] = u;
    }

    done[arc_t(t[0],t[1])] = true;
    done[arc_t(t[1],t[2])] = true;
    done[arc_t(t[2],t[0])] = true;
  }
  // Check consistency of orientation. It is consistent if and only if
  // each edge has been used exactly once in each direction.
  bool consistent = true;
  vector<edge_t> edges = undirected_edges();

  for(edge_t e: edges){
    if(!done[arc_t(e.first,e.second)]){
      fprintf(stderr,"A: Directed edge %d->%d is missing: triangulation is not consistently oriented.\n",e.first,e.second);
      consistent = false;
    }
    if(!done[arc_t(e.second,e.first)]){
      fprintf(stderr,"B: Directed edge %d->%d is missing: triangulation is not consistently oriented.\n",e.second,e.first);
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

// find_outer_face, edge_lengths, width_height, scale, move moved to layout2d.cc


ostream& operator<<(ostream& s, const PlanarGraph& g)
{
  vector<edge_t> edges = g.undirected_edges();

  s << "Graph[Range["<<g.N<<"],\n\tUndirectedEdge@@#&/@{";
  for(int i=0;i<edges.size();i++){
    s << "{" << (edges[i].first+1) << "," << (edges[i].second+1) << "}";
    if(i+1<edges.size())
      s << ", ";
    else
      s << "}";
  }
  s << "\n]";

  return s;
}


// **********************************************************************
//                       COMBINATORIAL PROPERTIES
// **********************************************************************

void perfmatch_dfs(map<arc_t,int>& faceEdge, const vector<face_t>& faces,
                   map<arc_t,int>& matrix, vector<bool>& faceSum, vector<bool>& visited, const arc_t& e)
{
  int frev = faceEdge[reverse(e)];
  if(visited[frev]) return;
  visited[frev] = true;

  const face_t &f(faces[frev]);
  for(int i=0;i<f.size();i++)
    perfmatch_dfs(faceEdge,faces,matrix,faceSum,visited,arc_t(f[i],f[(i+1)%f.size()]));

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
  map<arc_t,int> faceEdge;
  assert(is_consistently_oriented());
  vector<face_t> faces(compute_faces());
  vector<bool> faceSum(faces.size()), visited(faces.size());

  map<arc_t,int> A;
  vector<edge_t> edges = undirected_edges();
  for(edge_t e: edges){
    A[e] = 1;
    A[reverse(e)] = -1;
  }

  for(int i=0;i<faces.size();i++){
    const face_t &f(faces[i]);
    for(int j=0;j<f.size();j++){
      const arc_t e(f[j],f[(j+1)%f.size()]);
      faceEdge[e] = i;
      if(A[e] == 1) faceSum[i] = !faceSum[i];
    }
  }

  perfmatch_dfs(faceEdge,faces,A,faceSum,visited,edges[0]);

  vector<double> Af(N*N);
  for(map<arc_t,int>::const_iterator a(A.begin()); a!=A.end(); a++)
    Af[a->first.first*N+a->first.second] = a->second;

  return round(sqrtl(fabs(lu_det(Af,N))));
}
#else
size_t PlanarGraph::count_perfect_matchings() const
{
   cerr << "count_perfect_matchings() requires LAPACK.\n";
   //    cerr << "count_perfect_matchings() is temporarily out of service.\n";
  return 0;
}
#endif


vector<coord3d> PlanarGraph::zero_order_geometry(double scalerad) const
{
  vector<coord2d> flat_layout = tutte_layout();
  vector<coord2d> angles(layout2d::spherical_projection(*this, flat_layout));

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

 
// In an oriented planar graph, the directed edge starting in the smallest node
// is a unique representation of the face.
arc_t PlanarGraph::get_face_representation(arc_t e, int Fmax) const
{
  assert(is_consistently_oriented());


  int i=0;  
  arc_t e_min = e;
  node_t u = e.first, v = e.second;

  while(v!=e.first){
    node_t w = next_on_face(u,v);
    u=v; v=w; 

    if(u<e_min.first) e_min = {u,v};
    
    assert(w != -1);
    assert(++i<=Fmax); // Fmax is a back-stop to avoid infinite loops in a corrupted graph
  }
  return e_min;
}

// In an oriented planar graph, the directed edge starting in the smallest node
// is a unique representation of the face.
vector<arc_t> PlanarGraph::compute_face_representations(int Fmax) const
{
  assert(is_consistently_oriented());

  unordered_set<arc_t> faces(2*count_edges());
  
  for(node_t u=0;u<N;u++)
    for(node_t v: neighbours[u]){
      // For each directed edge, find the representative edge of the specified face
      // and assign an identifier
      faces.insert(get_face_representation({u,v},Fmax));
    }

  return vector<arc_t>(faces.begin(),faces.end());
}


face_t PlanarGraph::get_face_oriented(const arc_t &e, int Fmax) const
{
  assert(is_consistently_oriented());

  int i=0;
  node_t u = e.first, v=e.second;
  face_t f = vector<int>(1,u);
  
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
  vector<arc_t> face_representations = compute_face_representations(Fmax);

  vector<face_t> faces(face_representations.size());
  for(int i=0;i<face_representations.size();i++) faces[i] = get_face_oriented(face_representations[i],Fmax);

  // cerr << "facereps = " << face_representations << ";\n"
  //      << "faces    = " << faces << ";\n";
  
  return faces;
}


// permutation of vertex numbers (ie, replace v by vertex_numbers[v], to get numbered vertices)
// where permutations are as returned by PG.leapfrog_dual().get_spiral()
// locants are vertices that should have small vertex numbers (as far as permitted by symmetry equivalent canonical spirals)
vector<node_t> PlanarGraph::vertex_numbers(vector<vector<node_t>> &permutations, const vector<node_t> &locants) const{
  assert(!is_cubic());
  vector<node_t> vertex_numbers_inv(N,INT_MAX);
  for(int p=0; p<permutations.size(); p++){
    const vector<node_t> &perm=permutations[p];
    vector<node_t> vertex_numbers_tmp;
    // strip face-vertices, keep only vertex-vertices
    for(int i=0; i<perm.size(); i++){
      if(perm[i] < N) vertex_numbers_tmp.push_back(perm[i]);
    }
    assert(vertex_numbers_tmp.size() == N);
    //invert
    vector<node_t> vertex_numbers_inv_tmp(N);
    for(int i=0; i<vertex_numbers_tmp.size(); i++) vertex_numbers_inv_tmp[vertex_numbers_tmp[i]] = i;
    // copy to vertex_numbers_inv?
    if(locants.size()==0){
      vertex_numbers_inv = vertex_numbers_inv_tmp;
      break;
    }
    // compare two vectors, but only at chosen positions
    for(int l=0; l<locants.size(); l++){
      if(vertex_numbers_inv_tmp[locants[l]] > vertex_numbers_inv[locants[l]]) break;
      if(vertex_numbers_inv_tmp[locants[l]] < vertex_numbers_inv[locants[l]]){
        vertex_numbers_inv = vertex_numbers_inv_tmp;
        break;
      }
    }
  }
  //invert
  vector<node_t> vertex_numbers(N);
  for(int i=0; i<vertex_numbers.size(); i++) vertex_numbers[vertex_numbers_inv[i]] = i;
  return vertex_numbers; 
}


