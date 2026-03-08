#include <iomanip>
#include <limits>

#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/layout2d.hh"


double Polyhedron::diameter() const {
  double dmax = -INFINITY;
  for(int i=0;i<N;i++)
    for(int j=i+1;j<N;j++){
      double d = (points[i]-points[j]).norm();
      if(d > dmax) dmax = d;
    }
  return dmax;
};

// Helper: build the extended point list for centroid triangulation.
// Appends each face's centroid to the vertex list.
static vector<coord3d> centroid_points(std::span<const coord3d> points, const vector<face_t>& faces) {
  vector<coord3d> pts(points.begin(), points.end());
  for (const auto& f : faces) pts.push_back(f.centroid(points));
  return pts;
}

double Polyhedron::surface_area() const {
  vector<tri_t>  tris(centroid_triangulation(faces));
  vector<coord3d> pts(centroid_points(points, faces));

  double A = 0;
  for (const auto& tri : tris)
    A += Tri3D(pts[tri[0]], pts[tri[1]], pts[tri[2]]).area();
  return A;
}

// Signed volume via the divergence theorem: V = (1/3) ∮ r·n̂ dA.
// For a flat triangle, a·n̂ = b·n̂ = c·n̂ (constant on plane), so using
// any vertex gives the exact integral over the triangle.
double Polyhedron::volume_divergence() const {
  vector<tri_t>   tris(centroid_triangulation(faces));
  vector<coord3d> pts(centroid_points(points, faces));

  double V = 0;
  for (const auto& t : tris) {
    Tri3D T(pts[t[0]], pts[t[1]], pts[t[2]]);
    V += T.a.dot(T.n) * T.area() / T.n.norm();
  }
  return fabs(V / 3.0);
}

// Signed volume via signed tetrahedra from the origin:
// V = (1/6) Σ a·(b×c) for each triangle (a,b,c).
double Polyhedron::volume_tetra() const {
  vector<tri_t>   tris(centroid_triangulation(faces));
  vector<coord3d> pts(centroid_points(points, faces));

  double V = 0;
  for (const auto& t : tris) {
    const coord3d &a = pts[t[0]], &b = pts[t[1]], &c = pts[t[2]];
    V += a.dot(b.cross(c));
  }
  return fabs(V / 6.0);
}

Polyhedron Polyhedron::incremental_convex_hull() const {
  list<tri_t> output;
  typedef list<tri_t>::iterator triit;
  srandom(42); // Seed random numbers with constant for reproducible behaviour

  // 1. Create initial tetrahedron from 4 non-coplanar points.
  Tri3D T(points, tri_t(0, 1, 2));

  // Find the point furthest from the (0,1,2)-plane
  double distmax = 0;
  node_t v_best = 3;
  for (node_t u = 3; u < N; u++) {
    double dist = T.distance(points[u]);
    if (dist > distmax) { distmax = dist; v_best = u; }
  }

  output.push_back(tri_t(0, 1, 2));
  output.push_back(tri_t(0, 1, v_best));
  output.push_back(tri_t(0, 2, v_best));
  output.push_back(tri_t(1, 2, v_best));

  // Orient all initial faces outward (away from tetrahedron centroid)
  coord3d c((points[0] + points[1] + points[2] + points[v_best]) / 4.0);
  for (tri_t& t : output)
    if (!Tri3D(points, t).back_face(c)) t.flip();

  // 2. Incrementally add each remaining vertex
  for (node_t u = 3; u < N; u++) {
    if (u == v_best) continue;

    // Small perturbation for numerical robustness (breaks degeneracies from coplanar points)
    long r = random();
    coord3d perturbation(r & 0xff, (r >> 8) & 0xff, (r >> 16) & 0xff);
    coord3d p = points[u] * (coord3d(1, 1, 1) + perturbation * 1e-13);

    // 2.1 Find all faces visible from p
    list<triit> visible;
    map<arc_t, bool> is_visible;
    coord3d centre;
    for (triit t = output.begin(); t != output.end(); t++) {
      if (!Tri3D(points, *t).back_face(p)) {
        visible.push_back(t);
        for (int i = 0; i < 3; i++)
          is_visible[{t->u(i), t->u((i + 1) % 3)}] = true;
        centre += t->centroid(points);
      }
    }
    if (visible.empty()) continue;
    centre /= visible.size();

    // 2.2 Build horizon edges (boundary between visible and invisible faces)
    //     and delete visible faces
    list<edge_t> horizon;
    for (triit tv : visible) {
      for (int j = 0; j < 3; j++) {
        arc_t e{(*tv)[j], (*tv)[(j + 1) % 3]};
        arc_t e_rev{e.second, e.first};
        if (is_visible[e] != is_visible[e_rev])
          horizon.push_back(edge_t(e));
      }
      output.erase(tv);
    }

    // 2.3 Add new faces connecting u to each horizon edge, oriented outward
    for (const edge_t& e : horizon) {
      tri_t t(u, e.first, e.second);
      if (!Tri3D(points, t).back_face(centre)) t.flip();
      output.push_back(t);
    }

    if (output.size() > (size_t)N * N * 10) {
      cerr << "Convex hull error: output size " << output.size() << " exceeds safety limit\n";
    }
  }

  // 3. Construct oriented graph directly from hull triangles
  // 3.1 Remap vertices to contiguous indices
  set<node_t> used_nodes;
  for (const tri_t& t : output)
    for (int i = 0; i < 3; i++)
      used_nodes.insert(t.u(i));

  map<node_t, node_t> nodemap;
  vector<coord3d> remaining_points(used_nodes.size());
  node_t idx = 0;
  for (node_t u : used_nodes) {
    nodemap[u] = idx;
    remaining_points[idx] = points[u];
    idx++;
  }
  node_t M = remaining_points.size();

  // 3.2 Build arc_next map: for each directed edge (a,b) in a triangle (a,b,c),
  //     arc_next[(a,b)] = c. Also record one neighbour per vertex for fan start.
  map<arc_t, node_t> arc_next;
  vector<node_t> first_neighbour(M, -1);
  for (const tri_t& t : output) {
    node_t a = nodemap[t[0]], b = nodemap[t[1]], c = nodemap[t[2]];
    arc_next[{a, b}] = c;
    arc_next[{b, c}] = a;
    arc_next[{c, a}] = b;
    if (first_neighbour[a] < 0) first_neighbour[a] = b;
    if (first_neighbour[b] < 0) first_neighbour[b] = c;
    if (first_neighbour[c] < 0) first_neighbour[c] = a;
  }

  // 3.3 Build oriented neighbour lists by fan traversal.
  //     For vertex u, follow arc_next[(u, v)] around the fan.
  //     This produces CCW order (outward normals).
  neighbours_t nb(M, GRAPH_DMAX);
  for (node_t u = 0; u < M; u++) {
    node_t v = first_neighbour[u];
    node_t v0 = v;
    do {
      nb.push_back(u, v);
      v = arc_next.at({u, v});
    } while (v != v0);
  }

  // 3.4 Assemble oriented PlanarGraph and face list
  PlanarGraph g{Graph(nb)};
  vector<face_t> faces;
  faces.reserve(output.size());
  for (const tri_t& t : output)
    faces.push_back(tri_t(nodemap[t[0]], nodemap[t[1]], nodemap[t[2]]));
  return Polyhedron(g, remaining_points, 3, faces);
}

struct sort_ccw_coord3d {
  std::span<const coord3d> points;
  coord3d X, Y, n;

  // Points are only neighbour displacements from origin-node
  sort_ccw_coord3d(std::span<const coord3d> points)
    : points(points) {

    // TODO: More numerically robust method    
    coord3d xc(0,0,0);
    for(const coord3d &p: points) xc += p;
    xc /= points.size();

    n = /*x0*/-xc; n /= n.norm();	//

    coord3d x1 = points[0] /* - x0*/;
    X = x1 - n*x1.dot(n); X /= X.norm();
    Y = n.cross(x1);
  }

  bool operator()(const node_t& s, const node_t& t) const {
    coord3d xs(points[s]/*-x0*/), xt(points[t]/*-x0*/);
    coord2d Xs(X.dot(xs), Y.dot(xs)), Xt(X.dot(xt), Y.dot(xt));
    
    double angs = atan2(Xs.first,Xs.second), angt = atan2(Xt.first,Xt.second);
    return angs <= angt;
  } 
};



// Orient neighbours using planar_orient + 3D volume sign check.
// This is a free-standing helper for Polyhedron construction from unoriented data.
static void orient_polyhedron_neighbours(Polyhedron& P)
{
  // First get a consistent planar orientation
  layout2d::planar_orient(P);

  // Check volume sign to ensure outward-pointing normals (CCW-on-outside convention)
  double V=0;
  for(node_t u=0;u<P.N;u++){
    auto nu = P.nbrs(u);
    const coord3d ux(P.points[u]);

    for(int i=0;i<nu.size();i++){
      Tri3D T(ux, P.points[nu[i]],P.points[nu[(i+1)%nu.size()]]);
      V += ((T.a).dot(T.n))*T.area()/T.n.norm();
    }
  }

  if(V<0){ // Calculated normals are pointing inwards - reverse order.
    P.flip_all_orientations();
  }
}

Polyhedron::Polyhedron(const PlanarGraph& G, const vector<coord3d>& points_, const int face_max_, const vector<face_t> faces_) :
  PlanarGraph(G), face_max(face_max_), owned_points(points_), faces(faces_)
{
  repoint_coords();
  if(!is_consistently_oriented()) orient_polyhedron_neighbours(*this);

  if(faces.size() == 0){
    faces = compute_faces(face_max);
    face_max = 0;
    for(int i=0;i<faces.size();i++) if(faces[i].size() > face_max) face_max = faces[i].size();
  }
}

Polyhedron::Polyhedron(const PlanarGraph& G, std::span<coord3d> points_, const int face_max_, const vector<face_t> faces_) :
  PlanarGraph(G), face_max(face_max_), points(points_), faces(faces_)
{
  if(!is_consistently_oriented()) orient_polyhedron_neighbours(*this);

  if(faces.size() == 0){
    faces = compute_faces(face_max);
    face_max = 0;
    for(int i=0;i<faces.size();i++) if(faces[i].size() > face_max) face_max = faces[i].size();
  }
}

Polyhedron::Polyhedron(const vector<coord3d>& xs, double tolerance) 
{
  double bondlength = INFINITY;

  for(int i=0;i<xs.size();i++){
    for(int j=i+1;j<xs.size();j++){
      double d = (xs[i]-xs[j]).norm();
      if(d < bondlength) bondlength = d;
    }
  }
     
  neighbours_t nb(xs.size(), GRAPH_DMAX);
  for(int i=0;i<xs.size();i++){
    for(int j=i+1;j<xs.size();j++){
      double d = (xs[i]-xs[j]).norm();
      if(d <= bondlength*tolerance) {
        nb.push_back(i, j);
        nb.push_back(j, i);
      }
    }
  }

  *this = Polyhedron(PlanarGraph(Graph(nb)), xs);
}


matrix3d Polyhedron::inertia_matrix() const
{
  matrix3d I;

  for(int k=0;k<points.size();k++){
    const coord3d& x(points[k]);
    const long double xx(x.dot(x));
    for(int i=0;i<3;i++){
      I(i,i) += xx;

      for(int j=0;j<3;j++)
        I(i,j) -= x[i]*x[j];
    }
  }
  return I;
}

matrix3d Polyhedron::principal_axes() const
{
  const matrix3d I(inertia_matrix());
  pair<coord3d,matrix3d> ES(I.eigensystem());

  matrix3d Id;
  Id(0,0) = 1; 
  Id(1,1) = 1; 
  Id(2,2) = 1; 
/*
  cerr << "Inertial frame:\n " 
       << " inertia_matrix = " << I << ";\n"
       << " lambda  = " << ES.first << ";\n"
       << " vectors = " << ES.second << ";\n";
*/

  for(int i=0;i<3;i++) 
    if(std::isnan(ES.first[i])){
      cerr << "Warning: Inertial frame returned NaN. Setting inertial frame transformation to identity.\n";
      return Id;
    }
  
  if((ES.second*ES.second.transpose() - Id).norm() > 1e-2){
    cerr << "Warning: Inertial frame transform is not unitary. Setting inertial frame transformation to identity.\n";
    return Id;
  }

  return ES.second;
}

coord3d Polyhedron::width_height_depth() const {
  double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY,zmin=INFINITY,zmax=-INFINITY;
  for(node_t u=0;u<N;u++){
    const coord3d& x(points[u]);
    if(x[0]<xmin) xmin = x[0];
    if(x[0]>xmax) xmax = x[0];
    if(x[1]<ymin) ymin = x[1];
    if(x[1]>ymax) ymax = x[1];
    if(x[2]<zmin) zmin = x[2];
    if(x[2]>zmax) zmax = x[2];
  }
  return coord3d(xmax-xmin,ymax-ymin,zmax-zmin);
}



Polyhedron Polyhedron::dual() const 
{
  PlanarGraph d(dual_graph());

  vector<coord3d> coordinates(d.N);
  for(node_t u=0;u<d.N;u++){
    const face_t& f = faces[u];
    coord3d avg;
    for(auto v: f) avg += points[v];
    coordinates[u] = avg/double(f.size());
  }
 
  return Polyhedron(d,coordinates);
}


Polyhedron Polyhedron::leapfrog_dual() const 
{
  assert(is_consistently_oriented());
  size_t Nf = faces.size();

  Polyhedron Plf(Graph(N+Nf));
  Plf.owned_points.resize(N+Nf);
  Plf.repoint_coords();
   
  // Start with all the existing nodes
  for(node_t u=0;u<N;u++){
    Plf.assign_row(u, (*this)[u]);
    Plf.points[u]     = points[u];
  }

  // The result is a deltahedron: a polyhedron consisting of
  // only triangles. The first Nv points are the original vertices,
  // the last Nf are the midpoints of the original faces.
  int n_tris = 0;
  for(auto f: faces) n_tris += f.size();
  Plf.faces = vector<face_t>(n_tris);

  // Now connect new face-center nodes in oriented order
  for(int i=0;i<faces.size();i++){
    const face_t &f  = faces[i];
    node_t c = N+i;                // Face-center node
    
    // cerr << "new node " << c << " at face " << f << "\n";
    coord3d xc = {0,0,0};
    size_t   d = f.size();
    for(int j=0;j<d;j++){
      node_t u = f[j], v = f[(j+1)%f.size()];

      // Center node position is middle of face
      xc += points[u]/d;

      // Add edge mumble mumble
      Plf.insert_edge(arc_t{v,c},u,-1);

      // Add triangle
      Plf.faces[c] = tri_t{u,v,c};
    }
    Plf.points[c] = xc;
  }

  return Plf;
}

Polyhedron Polyhedron::fullerene_polyhedron(FullereneGraph G)
{
  Polyhedron P(G,G.zero_order_geometry(),6);
  P.set_points(G.optimized_geometry(P.points));

  P.move_to_origin();		// Center of mass at (0,0,0)
  P.align_with_axes();		// Align with principal axes
  
  return P;
}

bool Polyhedron::optimize(int opt_method, double ftol)
{
  if(is_a_fullerene()){
    FullereneGraph g(*this);
    set_points(g.optimized_geometry(points,opt_method,ftol));
    return true;
  } if(is_cubic()) {
    bool optimize_angles = true;
    return optimize_other(optimize_angles);
  } else if(is_triangulation()) {
    bool optimize_angles = false;
    return optimize_other(optimize_angles);
  }else{
     Triangulation LFD = leapfrog_dual();

    // inverse_tranges for faster lookup:
    // indices in 'triangles' at which there are triangles containing this vertex
    vector<vector<int>> inverse_triangle_list(LFD.N);
    auto lfd_tris = LFD.triangles();
    for (int i=0; i<(int)lfd_tris.size(); i++){
      const tri_t& tri = lfd_tris[i];
      inverse_triangle_list[tri[0]].push_back(i);
      inverse_triangle_list[tri[1]].push_back(i);
      inverse_triangle_list[tri[2]].push_back(i);
    }

    // generate and optimize LF of initial polyhedron
    PlanarGraph LF = LFD.dual_graph();
    Polyhedron P(LF,LF.zero_order_geometry());
    bool optimize_angles = true;
    bool opt_success = P.optimize_other(optimize_angles);

    // for each face in LF which corresponds to a vertex in the initial graph,
    // find the average coordinates of all vertices (ie the face centre)
    for (int i=0; i<N; i++){
      const int face_size(inverse_triangle_list[i].size());
      coord3d face_centre;
      for (int j=0; j<face_size; j++){
        face_centre += P.points[inverse_triangle_list[i][j]];
      }
      face_centre /= face_size;
      points[i] = face_centre/sqrt(3);
    }
    
    return opt_success;
  }
}

bool Polyhedron::is_triangulation() const {
  for(int i=0;i<faces.size();i++) if(faces[i].size()!=3) return false;
  return true;
}

// TODO: Add function for checking if forcefield convergence is achieved

bool Polyhedron::is_invalid() const {
  bool has_nans = false;
  for(auto p: points){
    if(std::isnan(p[0])||std::isnan(p[1])||std::isnan(p[2])) has_nans = true;
  }
  return has_nans;
}

pair<coord3d,coord3d> Polyhedron::bounding_box() const {
  coord3d xmin{INFINITY,INFINITY,INFINITY},xmax{-INFINITY,-INFINITY,-INFINITY};
  for(const auto& p: points)
    for(int i=0;i<3;i++){
      if(p[i] < xmin[i]) xmin[i] = p[i];
      if(p[i] > xmax[i]) xmax[i] = p[i];
    }
  return make_pair(xmin,xmax);
}

// Regular dodecahedron with unit bond length.
// Vertices: (±1,±1,±1), (0,±φ,±1/φ), (±1/φ,0,±φ), (±φ,±1/φ,0) scaled by φ/2,
// permuted to match the buckygen CW-oriented adjacency of FullereneGraph::C20().
double Polyhedron::dodecahedron_points[20][3] = {
  { 0.80901699437494745, 0.80901699437494745, 0.80901699437494745},
  { 0.00000000000000000, 1.30901699437494750, 0.50000000000000000},
  {-0.80901699437494745, 0.80901699437494745, 0.80901699437494745},
  {-0.50000000000000000, 0.00000000000000000, 1.30901699437494750},
  { 0.50000000000000000, 0.00000000000000000, 1.30901699437494750},
  { 0.80901699437494745,-0.80901699437494745, 0.80901699437494745},
  { 1.30901699437494750,-0.50000000000000000, 0.00000000000000000},
  { 1.30901699437494750, 0.50000000000000000, 0.00000000000000000},
  { 0.80901699437494745, 0.80901699437494745,-0.80901699437494745},
  { 0.00000000000000000, 1.30901699437494750,-0.50000000000000000},
  {-1.30901699437494750, 0.50000000000000000, 0.00000000000000000},
  {-1.30901699437494750,-0.50000000000000000, 0.00000000000000000},
  {-0.80901699437494745, 0.80901699437494745,-0.80901699437494745},
  {-0.80901699437494745,-0.80901699437494745, 0.80901699437494745},
  { 0.00000000000000000,-1.30901699437494750, 0.50000000000000000},
  { 0.00000000000000000,-1.30901699437494750,-0.50000000000000000},
  { 0.80901699437494745,-0.80901699437494745,-0.80901699437494745},
  { 0.50000000000000000, 0.00000000000000000,-1.30901699437494750},
  {-0.50000000000000000, 0.00000000000000000,-1.30901699437494750},
  {-0.80901699437494745,-0.80901699437494745,-0.80901699437494745}
};
