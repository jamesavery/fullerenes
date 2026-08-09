#include <iomanip>
#include <limits>

#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/layout2d.hh"


template<>
double PolyhedronView<double>::diameter() const {
  double dmax = -INFINITY;
  for(int i=0;i<N;i++)
    for(int j=i+1;j<N;j++){
      double d = (points[i]-points[j]).norm();
      if(d > dmax) dmax = d;
    }
  return dmax;
};

// The surface every integral below runs over -- documented at the declaration.
template<>
PolyhedronView<double>::CentroidSurface PolyhedronView<double>::centroid_surface() const {
  const vector<face_t> fs(faces());

  CentroidSurface S{vector<coord3d>(points.begin(), points.end()), centroid_triangulation(fs)};
  for (const face_t& f : fs) S.points.push_back(f.centroid(points));
  return S;
}

template<>
double PolyhedronView<double>::surface_area() const {
  const CentroidSurface S(centroid_surface());

  double A = 0;
  for (const tri_t& t : S.tris)
    A += Tri3D(S.points[t[0]], S.points[t[1]], S.points[t[2]]).area();
  return A;
}

// Signed volume via the divergence theorem: V = (1/3) ∮ r·n̂ dA.
// For a flat triangle, a·n̂ = b·n̂ = c·n̂ (constant on plane), so using
// any vertex gives the exact integral over the triangle.
template<>
double PolyhedronView<double>::volume_divergence() const {
  const CentroidSurface S(centroid_surface());

  double V = 0;
  for (const tri_t& t : S.tris) {
    Tri3D T(S.points[t[0]], S.points[t[1]], S.points[t[2]]);
    V += T.a.dot(T.n) * T.area() / T.n.norm();
  }
  return fabs(V / 3.0);
}

// Signed volume via signed tetrahedra from the origin:
// V = (1/6) Σ a·(b×c) for each triangle (a,b,c).
template<>
double PolyhedronView<double>::volume_tetra() const {
  const CentroidSurface S(centroid_surface());

  double V = 0;
  for (const tri_t& t : S.tris) {
    const coord3d &a = S.points[t[0]], &b = S.points[t[1]], &c = S.points[t[2]];
    V += a.dot(b.cross(c));
  }
  return fabs(V / 6.0);
}

template<>
Polyhedron PolyhedronView<double>::incremental_convex_hull() const {
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
  Graph nb(M, GRAPH_DMAX);
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
  return Polyhedron(g, remaining_points, 3);
}

struct sort_ccw_coord3d {
  std::span<const coord3d> points;
  coord3d X, Y, n;

  // Points are only neighbour displacements from origin-node
  sort_ccw_coord3d(std::span<const coord3d> points)
    : points(points) {

    // TODO: More numerically robust method
    const coord3d xc(mean(points));

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

Polyhedron::Polyhedron(const PlanarGraphView& G, const vector<coord3d>& points_, const int face_max_) :
  base_t(G), face_max(face_max_)
{
  owned_points = points_;
  repoint();
  if(!is_consistently_oriented()) orient_polyhedron_neighbours(*this);

  if(face_max == INT_MAX){
    auto fs = faces();
    face_max = 0;
    for(auto& f: fs) if(int(f.size()) > face_max) face_max = f.size();
  }
}

Polyhedron::Polyhedron(const PlanarGraphView& G, std::span<coord3d> points_, const int face_max_) :
  base_t(G), face_max(face_max_)
{
  // Copy coordinates into owned storage
  owned_points.assign(points_.begin(), points_.end());
  repoint();
  if(!is_consistently_oriented()) orient_polyhedron_neighbours(*this);

  if(face_max == INT_MAX){
    auto fs = faces();
    face_max = 0;
    for(auto& f: fs) if(int(f.size()) > face_max) face_max = f.size();
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
     
  Graph nb(xs.size(), GRAPH_DMAX);
  for(int i=0;i<xs.size();i++){
    for(int j=i+1;j<xs.size();j++){
      double d = (xs[i]-xs[j]).norm();
      if(d <= bondlength*tolerance) {
        nb.push_back(i, j);
        nb.push_back(j, i);
      }
    }
  }

  *this = Polyhedron(PlanarGraph(nb), xs);
}


// A second moment is a Gram matrix -- an integral of x x^T against a nonnegative
// measure -- and is therefore positive semidefinite by construction.  A negative
// eigenvalue is not a small error but a proof that the input was not the mass
// distribution it claimed to be, and the sign test costs no tolerance.
//
// It is the ONLY detector for a MIS-ORIENTED face list.  Half the surface wound
// inwards makes det_sum cancel while the first moment does not, so the enclosed
// volume survives the `vol > 0` guard as a near-zero number, the centroid ratio
// blows up, and the covariance still comes out FINITE -- every other guard here
// and in principal_frame() passes, and the caller gets a confidently wrong frame.
// Measured on the mis-oriented cube in tests/volume-moments-test.cc: unguarded
// volume 2.5e-4 (> 0), covariance norm 188 (finite), lambda_min = -188.
//
// matrix3d::eigenvalues() returns descending, so [2] is lambda_min.
static bool is_positive_semidefinite(const matrix3d& M) { return M.eigenvalues()[2] >= 0; }

// The moments of the solid enclosed by an oriented triangle surface.  The
// arithmetic is documented at the declaration in graphview.hh; the apex `o` is
// taken at the vertex mean, a pure translation (the moments are exactly
// translation-covariant) bought for nothing: summing 10^5 tetrahedra whose apex
// sits 10^3 units outside the body loses ~5 digits to cancellation before the
// answer appears.
MassMoments volume_moments(std::span<const coord3d> P, const vector<tri_t>& tris)
{
  MassMoments m;
  if(P.size() < 4 || tris.empty()) return m;

  const coord3d o(mean(P));

  double   det_sum = 0;
  coord3d  m1{0,0,0};
  matrix3d m2{};
  for(const tri_t& t: tris){
    const coord3d a = P[t[0]]-o, b = P[t[1]]-o, c = P[t[2]]-o;
    const double  det = a.dot(b.cross(c));
    const coord3d s = a+b+c;
    det_sum += det;
    m1 += s*det;
    m2 += (s.outer(s) + a.outer(a) + b.outer(b) + c.outer(c)) * (det/120.0);
  }
  // An inward-oriented surface gives det_sum < 0, which would negate the
  // covariance and invert its eigenvalue ordering (the longest axis would sort
  // last).  Normalise the sign here so no caller has to know the winding; the
  // centroid is a ratio and is unaffected.
  const double sgn = (det_sum < 0)? -1.0 : 1.0;
  const double vol = sgn*det_sum/6.0;
  if(!(vol > 0)) return m;
  const coord3d  cr = m1 * (sgn/(24.0*vol));   // centroid relative to the apex o
  const matrix3d cov = m2*sgn - cr.outer(cr)*vol;
  if(!is_positive_semidefinite(cov)) return MassMoments{MassMoments::Code::NotPositiveSemidefinite};

  m.code       = MassMoments::Code::Ok;
  m.mass       = vol;
  m.centroid   = o + cr;
  m.covariance = cov;
  return m;
}

// The point-mass sibling: unit mass at each point.  Same triple, same
// postcondition, no topology -- documented at the declaration in graphview.hh.
MassMoments vertex_moments(std::span<const coord3d> P)
{
  MassMoments m;
  if(P.empty()) return m;

  const coord3d c(mean(P));

  matrix3d cov{};
  for(const coord3d& p: P){ const coord3d x(p-c); cov += x.outer(x); }
  if(!is_positive_semidefinite(cov)) return MassMoments{MassMoments::Code::NotPositiveSemidefinite};

  m.code       = MassMoments::Code::Ok;
  m.mass       = (double)P.size();
  m.centroid   = c;
  m.covariance = cov;
  return m;
}

// Inertia tensor from a CENTRAL second-moment matrix M: I = tr(M) Id - M.
// Shared by both mass models so the identity is written once.
static matrix3d inertia_from_second_moment(const matrix3d& M)
{
  return matrix3d::unit_matrix()*M.trace() - M;
}

// The faces are triangulated the way the rest of the class triangulates them --
// centroid_surface(), the same surface surface_area() / volume_tetra() /
// volume_divergence() integrate over -- so this method's `mass` (the enclosed
// volume, under this weighting) and volume_tetra() are the same integral rather
// than two slightly different answers to the same question.  A face of a real cage is not exactly planar,
// and the two triangulation conventions genuinely differ: fan-triangulating
// (v0,vi,vi+1) instead moves the C60 volume by 3.3e-5 relative AND makes it
// depend on which vertex the face list happens to start at, which centroid
// triangulation does not.
template<>
MassMoments PolyhedronView<double>::volume_moments() const
{
  const CentroidSurface S(centroid_surface());
  return ::volume_moments(std::span<const coord3d>(S.points.data(), S.points.size()), S.tris);
}

// The moments of the chosen mass distribution -- the one place the two models
// differ, so that everything downstream is written once.
static MassMoments moments_of(const PolyhedronView<double>& P, MassModel m)
{
  return m == MassModel::Atoms
       ? ::vertex_moments(std::span<const coord3d>(P.points.data(), P.points.size()))
       : P.volume_moments();
}

// Both models are CENTRAL: the second moment is taken about the mass
// distribution's own centre (the volume centroid / the vertex mean), so the
// tensor -- and every frame built on it -- is invariant under translating the
// polyhedron.  Before 2026-08-07 the vertex sum was taken about the ORIGIN, so
// an off-centre cage got a frame that depended on where it happened to sit.
//
// A degenerate or non-PSD distribution has a zeroed covariance, so this yields
// the ZERO tensor -- a plausible-looking wrong answer that no guard downstream
// catches.  principal_axes() on the zero tensor returns the identity, but NOT
// because either of its checks fires: the eigenvalues are 0, not NaN, and the
// returned matrix is perfectly unitary.  The identity comes from
// matrix3d::eigensystem()'s ISOTROPIC branch (geometry.hh), which returns the
// standard basis whenever the deviatoric part vanishes -- and the zero matrix is
// isotropic.  The fact is therefore only reachable through the named outcome:
// moments_of(...).code, or principal_frame()'s translation of it.
template<>
matrix3d PolyhedronView<double>::inertia_matrix(MassModel m) const
{
  return inertia_from_second_moment(moments_of(*this, m).covariance);
}

// The principal frame, with the reason attached when it is the identity fallback
// rather than a frame.  Until 2026-08-08 the two guards below wrote a sentence to
// cerr and returned the identity, which no caller could branch on and no sweep
// could count; the codes are that information on the return channel.
template<>
InertialFrame PolyhedronView<double>::principal_frame(MassModel m) const
{
  const MassMoments M(moments_of(*this, m));
  switch(M.code){
  case MassMoments::Code::Degenerate:
    return {InertialFrame::Code::DegenerateMass};
  case MassMoments::Code::NotPositiveSemidefinite:
    return {InertialFrame::Code::NotPositiveSemidefinite};
  case MassMoments::Code::Ok: break;
  }

  const matrix3d I(inertia_from_second_moment(M.covariance));
  const pair<coord3d,matrix3d> ES(I.eigensystem());

  for(int i=0;i<3;i++)
    if(std::isnan(ES.first[i])) return {InertialFrame::Code::NonFiniteTensor};

  if((ES.second*ES.second.transpose() - matrix3d::unit_matrix()).norm() > 1e-2)
    return {InertialFrame::Code::NotUnitary};

  return {InertialFrame::Code::Ok, ES.second};
}

template<>
matrix3d PolyhedronView<double>::principal_axes(MassModel m) const
{
  return principal_frame(m).axes;
}

template<>
coord3d PolyhedronView<double>::width_height_depth() const {
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



template<>
Polyhedron PolyhedronView<double>::dual() const
{
  PlanarGraph d(dual_graph());
  auto fs = faces();

  vector<coord3d> coordinates(d.N);
  for(node_t u=0;u<d.N;u++){
    const face_t& f = fs[u];
    coord3d avg;
    for(auto v: f) avg += points[v];
    coordinates[u] = avg/double(f.size());
  }
 
  return Polyhedron(d,coordinates);
}


template<>
Polyhedron PolyhedronView<double>::leapfrog_dual() const
{
  assert(is_consistently_oriented());
  auto fs = faces();
  size_t Nf = fs.size();

  Polyhedron Plf;
  static_cast<Owned<PolyhedronView<double>>&>(Plf) = Owned<PolyhedronView<double>>(int(N+Nf));
  // points already allocated by Owned<PolyhedronView<double>>(N+Nf)

  // Start with all the existing nodes
  for(node_t u=0;u<N;u++){
    Plf.assign_row(u, (*this)[u]);
    Plf.points[u]     = points[u];
  }

  // Now connect new face-center nodes in oriented order.
  // The result is a deltahedron: all triangular faces.
  for(int i=0;i<int(Nf);i++){
    const face_t &f  = fs[i];
    node_t c = N+i;                // Face-center node

    coord3d xc = {0,0,0};
    size_t   d = f.size();
    for(int j=0;j<int(d);j++){
      node_t u = f[j], v = f[(j+1)%f.size()];

      // Center node position is middle of face
      xc += points[u]/d;

      // Add edge
      Plf.insert_edge(arc_t{v,c},u,-1);
    }
    Plf.points[c] = xc;
  }

  Plf.face_max = 3;
  return Plf;
}

Polyhedron Polyhedron::fullerene_polyhedron(FullereneGraph G)
{
  // Start geometry: the Eisenstein-paint warm start (deterministic,
  // shape-respecting) with the legacy Tutte-on-sphere start as fallback
  // for the rare isomers where the paint pipeline reports failure.
  vector<coord3d> x0;
  try {
    x0 = G.eisenstein_paint_geometry();
  } catch (const std::exception& e) {
    cerr << "fullerene_polyhedron: " << e.what()
         << " -- falling back to zero_order_geometry\n";
    x0 = G.zero_order_geometry();
  }

  Polyhedron P(G,x0,6);
  P.owned_points = G.optimized_geometry(P.points);
  P.repoint();

  P.move_to_origin();		// Center of mass at (0,0,0)
  P.align_with_axes(MassModel::Atoms);

  return P;
}

template<>
bool PolyhedronView<double>::optimize(int opt_method, double ftol)
{
  if(is_a_fullerene()){
    FullereneGraph g(*this);
    auto new_pts = g.optimized_geometry(points,opt_method,ftol);
    std::copy(new_pts.begin(), new_pts.end(), points.begin());
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
  for(auto& f: faces(face_max)) if(f.size()!=3) return false;
  return true;
}

// TODO: Add function for checking if forcefield convergence is achieved

template<>
bool PolyhedronView<double>::is_invalid() const {
  bool has_nans = false;
  for(auto p: points){
    if(std::isnan(p[0])||std::isnan(p[1])||std::isnan(p[2])) has_nans = true;
  }
  return has_nans;
}

template<>
pair<coord3d,coord3d> PolyhedronView<double>::bounding_box() const {
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
