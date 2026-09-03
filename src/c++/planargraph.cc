#include <queue>
#include <stdexcept>

#include "fullerenes/spiral.hh"
#include "fullerenes/planargraph.hh"
#include "fullerenes/layout2d.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/cubicgraph.hh"
#include "fullerenes/dense_linalg.hh"

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
PlanarGraph PlanarGraphView::enveloping_triangulation(construction_scheme_t &scheme) const
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

PlanarGraph PlanarGraphView::enveloping_triangulation(const construction_scheme_t &scheme) const
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

bool PlanarGraphView::is_cubic() const {
  for(node_t u=0;u<N;u++)
    if(degree(u) != 3)
      return false;
  return true;
}

bool PlanarGraphView::is_triangulation() const { // NB: A bit expensive
  vector<face_t> faces(compute_faces());

  for(int i=0;i<faces.size();i++) if(faces[i].size() != 3) return false;
  return true;
}

bool PlanarGraphView::is_a_fullerene(bool verbose) const {
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
bool PlanarGraphView::is_cut_vertex(const node_t v) const {
  // Requires oriented (sorted) neighbours of v (direction doesn't matter)
  auto nv = nbrs(v);
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



PlanarGraph PlanarGraphView::dual_graph(unsigned int Fmax) const
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
      dual.push_back(i_f, face_numbers[e_g]);

      w = prev(v,u); u = v; v = w;
      if(++i > Fmax) throw std::logic_error("dual_graph: face exceeds Fmax (corrupted/non-oriented graph)");
    } while (u != e_f.first);
  }
  if(!dual.is_consistently_oriented())
    throw std::logic_error("dual_graph: produced dual is not consistently oriented");

  return dual;
}

// the dual of the LF, ie a Triangulation is returned
PlanarGraph PlanarGraphView::leapfrog_dual() const
{
  vector<face_t> faces = compute_faces_oriented();
  size_t Nf = faces.size();

  PlanarGraph lf(N+Nf);

  // Start with all the existing nodes
  for(node_t u=0;u<N;u++) lf.assign_row(u, (*this)[u]);

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


vector<face_t> PlanarGraphView::compute_faces(unsigned int Fmax) const {
  return compute_faces_oriented(Fmax);
}



vector<tri_t> PlanarGraphView::triangulation(int face_max) const
{
  vector<face_t> faces(compute_faces(face_max));
  return triangulation(faces);
}

// A face too small to cut into triangles.  Both triangulations below index f[2],
// so this is the guard that keeps a corrupt face list from reading past its end
// -- and the fan's loop bound j+1 < f.size() from underflowing on an empty one.
static void require_polygon(const face_t& f, const char* operation)
{
  if(f.size() < 3)
    throw std::logic_error(string(operation) + ": face with " + to_string(f.size())
                           + " vertices cannot be triangulated (needs at least 3)");
}

// Split face i of more than 3 vertices at its centroid, the new vertex N+i.  The
// documented contract, including why no orientation-repair pass follows, is at
// the declaration in graphview.hh.
//
// Each triangle keeps its face's arc f[j] -> f[j+1] and adds the two spokes
// f[j+1] -> c -> f[j]; the next triangle traverses that same spoke as
// f[j+1] -> c, so every added arc appears once in each direction and the face
// arcs are untouched.  Hence: oriented faces in, oriented surface out -- and a
// triangle face, which is its own centroid triangulation, comes out unchanged,
// so the two cases agree on the boundary arcs rather than winding oppositely
// (which is what made a mixed triangle/polygon mesh come out inconsistent before
// 2026-08-09).
vector<tri_t> PlanarGraphView::centroid_triangulation(const vector<face_t>& faces) const
{
  vector<tri_t> tris;
  for(size_t i=0;i<faces.size();i++){
    const face_t& f(faces[i]);
    require_polygon(f, "PlanarGraphView::centroid_triangulation");
    const node_t c = N+int(i);            // face i's centroid vertex

    if(f.size() == 3) tris.push_back({f[0],f[1],f[2]});
    else for(size_t j=0;j<f.size();j++)
      tris.push_back({f[j],f[(j+1)%f.size()],c});
  }
  return tris;
}


// Fan from f[0].  Same contract, no new vertices: the interior diagonal f[0]-f[j]
// is traversed f[j] -> f[0] by triangle j-1 and f[0] -> f[j] by triangle j, so
// again every added arc appears once in each direction.  A triangle face is the
// j == 1 term alone, so it too comes out unchanged and needs no special case.
vector<tri_t> PlanarGraphView::triangulation(const vector<face_t>& faces) const
{
  vector<tri_t> tris;
  for(const face_t& f: faces){
    require_polygon(f, "PlanarGraphView::triangulation");
    for(size_t j=1;j+1<f.size();j++) tris.push_back(tri_t(f[0],f[j],f[j+1]));
  }
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

// |det| of a flat row-major N×N matrix via the in-house BLAS-free LU
// (LinAlg::det).  Formerly routed through LAPACK dgetrf_, which the deployed
// OpenBLAS silently miscomputes for N ≳ 60 — exactly the FKT regime below.
double lu_det(const vector<double> &A, int N)
{
  return fabs(LinAlg::det(matrix<double>(N, N, A)));
}


size_t PlanarGraphView::count_perfect_matchings() const
{
  map<arc_t,int> faceEdge;
  require_oriented_surface(*this, "PlanarGraphView::count_perfect_matchings");
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


vector<coord3d> PlanarGraphView::tutte_sphere_geometry(double scalerad) const
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
    for(int i=0;i<3;i++) Ravg += (coordinates[u]-coordinates[nbrs(u)[i]]).norm();
  Ravg /= (3.0*N);

  coordinates *= scalerad*1.5/Ravg;

  return coordinates;
}

 
// In an oriented planar graph, the directed edge starting in the smallest node
// is a unique representation of the face.
arc_t PlanarGraphView::get_face_representation(arc_t e, int Fmax) const
{
  // Precondition (consistent orientation) is a maintained global invariant and is
  // checked once (throwing) by the batch callers (compute_face_representations,
  // dual_graph). It is NOT re-checked here: this accessor is called O(E) times per
  // batch, and is_consistently_oriented() is itself O(E) -- re-checking would make
  // face computation O(E^2). The Fmax back-stop below is the always-compiled guard
  // that turns a corrupted/non-oriented graph into a throw rather than an infinite loop.
  int i=0;
  arc_t e_min = e;
  node_t u = e.first, v = e.second;

  while(v!=e.first){
    node_t w = next_on_face(u,v);
    u=v; v=w;

    if(u<e_min.first) e_min = {u,v};

    if(w == -1)    throw std::logic_error("get_face_representation: dangling arc (corrupted/non-oriented graph)");
    if(++i > Fmax) throw std::logic_error("get_face_representation: face exceeds Fmax (corrupted/non-oriented graph)");
  }
  return e_min;
}

// In an oriented planar graph, the directed edge starting in the smallest node
// is a unique representation of the face.
vector<arc_t> PlanarGraphView::compute_face_representations(int Fmax) const
{
  // @pre: graph is consistently oriented. Checked once here (O(E)) so the O(E)
  // per-arc tracer below can trust it; a violation throws rather than aborting.
  if(!is_consistently_oriented())
    throw std::logic_error("compute_face_representations: graph is not consistently oriented");

  unordered_set<arc_t> faces(2*count_edges());
  
  for(node_t u=0;u<N;u++)
    for(node_t v: nbrs(u)){
      // For each directed edge, find the representative edge of the specified face
      // and assign an identifier
      faces.insert(get_face_representation({u,v},Fmax));
    }

  return vector<arc_t>(faces.begin(),faces.end());
}


face_t PlanarGraphView::get_face_oriented(const arc_t &e, int Fmax) const
{
  // See get_face_representation: orientation precondition is checked once by the
  // batch caller (compute_faces_oriented), not per-arc, to keep face computation
  // linear rather than O(E^2). The Fmax back-stop below throws (always compiled)
  // on a corrupted/non-oriented graph rather than looping forever.
  int i=0;
  node_t u = e.first, v=e.second;
  face_t f = vector<int>(1,u);

  while(v!=e.first){
    node_t w = prev(v,u);        // Previous neighbour to u in v defines corner u-v-w in face

    f.push_back(v);
    u=v; v=w; i++;
    if(w == -1)  throw std::logic_error("get_face_oriented: dangling arc (corrupted/non-oriented graph)");
    if(i > Fmax) throw std::logic_error("get_face_oriented: face exceeds Fmax (corrupted/non-oriented graph)");
  }
  return f;
}

vector<face_t> PlanarGraphView::compute_faces_oriented(int Fmax) const
{
  // @pre (consistent orientation) is enforced once by compute_face_representations,
  // called first below, before any face is traced.
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
vector<node_t> PlanarGraphView::vertex_numbers(vector<vector<node_t>> &permutations, const vector<node_t> &locants) const{
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


