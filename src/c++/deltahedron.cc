#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"

Deltahedron::Deltahedron(const Triangulation& T, const vector<coord3d>& points)
  : Triangulation(T), points(points)
{
  assert((int)points.size() == N);
}

Deltahedron::Deltahedron(const Polyhedron& P)
  : Triangulation(static_cast<const Graph&>(P)), points(P.points)
{
  assert(P.is_triangulation());
}

vector<face_t> Deltahedron::compute_dual_faces() const {
  vector<face_t> faces(triangles.size());
  for(size_t i = 0; i < triangles.size(); i++)
    faces[i] = face_t(triangles[i]);
  return faces;
}

void Deltahedron::smooth(double q) {
  vector<coord3d> new_points(N);
  for(node_t u = 0; u < N; u++){
    coord3d avg;
    for(node_t v : neighbours[u]) avg += points[v];
    avg /= neighbours[u].size();
    new_points[u] = points[u]*(1.0-q) + avg*q;
  }
  points = new_points;
}

Deltahedron Deltahedron::GCtransform(unsigned k, unsigned l) const {
  assert(l == 0 && "General GC(k,l) with l!=0 not yet implemented");
  int m = k - 1;
  int n = k;  // = m+1

  // 1. Get topology + face grids from existing halma_transform
  vector<map<edge_t,node_t>> face_grids;
  Triangulation T_new = Triangulation::halma_transform(m, &face_grids);

  // 2. Compute 3D coordinates for all vertices via barycentric interpolation
  //
  // Each face grid maps edge_t{a,b} -> node_id, where {a,b} are stored as
  // {min,max} by the edge_t constructor. The three triangle corners sit at
  // {0,0}, {0,n}, {n,n} where n = k = m+1.
  //
  // The barycentric coordinates of grid point {a,b} are:
  //   lambda_0 = (n-b)/n    weight of T[0]
  //   lambda_1 = (b-a)/n    weight of T[1]
  //   lambda_2 = a/n        weight of T[2]
  //
  // With the scaling built in:
  //   point({a,b}) = (n-b)*P[T[0]] + (b-a)*P[T[1]] + a*P[T[2]]
  //
  // Edge/corner vertices appear in multiple face grids but all give the
  // same coordinate (barycentric on shared edges/corners), so overwrites
  // are harmless.

  vector<coord3d> new_points(T_new.N);

  for(int i = 0; i < (int)triangles.size(); i++){
    const face_t& T = triangles[i];
    const auto& grid = face_grids[i];

    for(const auto& [ab, node_id] : grid){
      int a = ab.first, b = ab.second;  // stored as {min, max}
      new_points[node_id] = points[T[0]]*(n-b) + points[T[1]]*(b-a) + points[T[2]]*a;
    }
  }

  return Deltahedron(T_new, new_points);
}
