#include <iomanip>
#include <limits>

#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"


double Polyhedron::diameter() const {
  double dmax = -INFINITY;
  for(int i=0;i<N;i++)
    for(int j=i+1;j<N;j++){
      double d = (points[i]-points[j]).norm();
      if(d > dmax) dmax = d;
    }
  return dmax;
};

double Polyhedron::surface_area() const {
  double A = 0;

  // vector<tri_t> tris(triangulation(faces));
  
  // for(size_t i=0;i<tris.size();i++){
  //   const tri_t& tri(tris[i]);
  //   Tri3D T(points[tri[0]],points[tri[1]],points[tri[2]]);
  //   A += T.area();
  // }
  // return A;

  vector<tri_t> tris(centroid_triangulation(faces));
  vector<coord3d> centroid_points(points.begin(),points.end());
  for(int i=0;i<faces.size();i++) centroid_points.push_back(faces[i].centroid(points));
  
  for(size_t i=0;i<tris.size();i++){
    const tri_t& tri(tris[i]);
    Tri3D T(centroid_points[tri[0]],centroid_points[tri[1]],centroid_points[tri[2]]);
    A += T.area();
  }
  return A;
}

double Polyhedron::volume_tetra() const {
  vector<tri_t>   tris(centroid_triangulation(faces));
  vector<coord3d> centroid_points(points.begin(),points.end());
  for(int i=0;i<faces.size();i++) centroid_points.push_back(faces[i].centroid(points));

  double V = 0,Vm=0,Vp=0;

  // Now generate tetrahedra and either add or subtract volume according to which direction the face is pointing
  coord3d zero(0,0,0);
  for(size_t i=0;i<tris.size();i++){
    const face_t& t(tris[i]);
    Tri3D T(centroid_points[t[0]],centroid_points[t[1]],centroid_points[t[2]]);
    double dV = Tetra3D(T.a,T.b,T.c,zero).volume();
    V += (T.back_face(zero)? 1 : -1)*dV;
    if(T.back_face(zero)) Vp += dV;
    else Vm += dV;
    //    if(!T.back_face(zero))
      //      cerr << "Tri " << t << " / " << T << " is a front face.\n";
  }
  fprintf(stderr,"V = %f - %f = %f\n",Vp,Vm,V);
  return fabs(V);
}

double Polyhedron::volume_divergence() const {
  vector<tri_t>   tris(centroid_triangulation(faces));
  vector<coord3d> centroid_points(points.begin(),points.end());
  for(int i=0;i<faces.size();i++) centroid_points.push_back(faces[i].centroid(points));

  //  cerr << "points = {"; for(int i=0;i<centroid_points.size();i++) cerr << centroid_points[i] << (i+1<centroid_points.size()? ", ":"};\n");

  double V = 0;

  // Now generate tetrahedra and either add or subtract volume according to which direction the face is pointing
  for(size_t i=0;i<tris.size();i++){
    const face_t& t(tris[i]);
    Tri3D T(centroid_points[t[0]],centroid_points[t[1]],centroid_points[t[2]]);

    V += ((T.a).dot(T.n))*T.area()/T.n.norm();
  }
  return fabs(V/3.0);
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
  neighbours_t nb(M);
  for (node_t u = 0; u < M; u++) {
    node_t v = first_neighbour[u];
    node_t v0 = v;
    do {
      nb[u].push_back(v);
      v = arc_next.at({u, v});
    } while (v != v0);
  }

  // 3.4 Assemble oriented PlanarGraph and face list
  PlanarGraph g(Graph(nb, true));
  vector<face_t> faces;
  faces.reserve(output.size());
  for (const tri_t& t : output)
    faces.push_back(tri_t(nodemap[t[0]], nodemap[t[1]], nodemap[t[2]]));
  g.outer_face = faces[0];

  return Polyhedron(g, remaining_points, 3, faces);
}

struct sort_ccw_coord3d {
  const vector<coord3d> &points;
  coord3d X, Y, n;

  // Points are only neighbour displacements from origin-node
  sort_ccw_coord3d(const vector<coord3d>& points)
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



void Polyhedron::orient_neighbours()
{
  
  if(layout2d.size() != N)
    layout2d = tutte_layout();
  PlanarGraph::orient_neighbours();

  // TODO: Sort this out
  //  else if(points.size() == N){
  // Orient neighbours locally (CW or CCW depending on luck)
  // for(node_t u=0;u<N;u++){
  //   vector<node_t> &ns(neighbours[u]);
  //   size_t degree = ns.size();
      
  //   int ns_index[degree];
  //   for(int i=0;i<degree;i++) ns_index[i] = i;
				     
  //   vector<coord3d> neighbour_points(degree);
  //   const coord3d &x0 = points[u];

  //   for(int i=0;i<degree;i++) neighbour_points[i] = points[ns[i]] - x0;
  //   sort_ccw_coord3d CCW(neighbour_points);

  //   sort(ns_index,ns_index+degree,CCW);
  //   vector<node_t> ns_sorted(degree);
  //   for(int i=0;i<degree;i++) ns_sorted[i] = ns[ns_index[i]];
  //   ns = ns_sorted;
  // }


  // // Choose that first node is correct, then flip to consistency
  // // TODO: How? For now, cheat slowly.
  // if(is_consistently_oriented())
  //   is_oriented = true;
  // else {
  //   layout2d = tutte_layout();
  //   PlanarGraph::orient_neighbours();
  // }
  //}
  
  // Calculate volume
  double V=0;
  for(node_t u=0;u<N;u++){
    const face_t nu(neighbours[u]);
    const coord3d ux(points[u]);

    for(int i=0;i<nu.size();i++){
      Tri3D T(ux, points[nu[i]],points[nu[(i+1)%nu.size()]]);
      V += ((T.a).dot(T.n))*T.area()/T.n.norm();
    }
  }

  if(V<0){ // Calculated normals are pointing inwards - reverse order.
    //    printf("Inverted normals - reversing neighbours lists.\n");
    for(node_t u=0;u<N;u++) reverse(neighbours[u].begin(), neighbours[u].end());
  }
}

Polyhedron::Polyhedron(const PlanarGraph& G, const vector<coord3d>& points_, const int face_max_, const vector<face_t> faces_) : 
  PlanarGraph(G), face_max(face_max_), points(points_), faces(faces_)
{
  if(!is_oriented) orient_neighbours();

  //  for(node_t u=0;u<N;u++) points[u] = points_[u];
  
  if(faces.size() == 0){
    faces = compute_faces(face_max);
    assert(outer_face.size() <= face_max);
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
     
  set<edge_t> edges;
  for(int i=0;i<xs.size();i++){
    for(int j=i+1;j<xs.size();j++){
      double d = (xs[i]-xs[j]).norm();
      if(d <= bondlength*tolerance) {
        edges.insert(edge_t(i,j));
      }
    }
  }
  
  (*this) = Polyhedron(PlanarGraph(edges), xs);
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
  assert(is_oriented);
  size_t Nf = faces.size();
  
  Polyhedron Plf(Graph(N+Nf,true));
  Plf.points.reserve(N+Nf);
   
  // Start with all the existing nodes
  for(node_t u=0;u<N;u++){
    Plf.neighbours[u] = neighbours[u];
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
  if(G.layout2d.empty())
    G.layout2d       = G.tutte_layout();
  
  Polyhedron P(G,G.zero_order_geometry(),6);
  P.points = G.optimized_geometry(P.points);

  P.move_to_origin();		// Center of mass at (0,0,0)
  P.align_with_axes();		// Align with principal axes
  
  return P;
}

bool Polyhedron::optimize(int opt_method, double ftol)
{
  if(is_a_fullerene()){
    FullereneGraph g(*this);
    points = g.optimized_geometry(points,opt_method,ftol);
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
    for (int i=0; i<LFD.triangles.size(); i++){
      const tri_t& tri = LFD.triangles[i];
      inverse_triangle_list[tri[0]].push_back(i);
      inverse_triangle_list[tri[1]].push_back(i);
      inverse_triangle_list[tri[2]].push_back(i);
    }

    // generate and optimize LF of initial polyhedron
    PlanarGraph LF = LFD.dual_graph();
    LF.layout2d = LF.tutte_layout();
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
    if(std::isnan(p[0])||std::isnan(p[1])||std::isnan(p[1])) has_nans = true;
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
