#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/unfold.hh"
#include <cmath>

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
  if(l == 0){
    // Halma path: direct subdivision via face grids, preserves node IDs
    int m = k - 1;
    int n = k;  // = m+1

    vector<map<edge_t,node_t>> face_grids;
    Triangulation T_new = Triangulation::halma_transform(m, &face_grids);

    // Barycentric interpolation in each face grid.
    // Grid corners at {0,0}, {0,n}, {n,n} map to T[0], T[1], T[2].
    // point({a,b}) = (n-b)*P[T[0]] + (b-a)*P[T[1]] + a*P[T[2]]
    // Weights sum to n=k, so corner vertices get k*P_original.
    vector<coord3d> new_points(T_new.N);

    for(int i = 0; i < (int)triangles.size(); i++){
      const face_t& T = triangles[i];
      const auto& grid = face_grids[i];

      for(const auto& [ab, node_id] : grid){
        int a = ab.first, b = ab.second;
        new_points[node_id] = points[T[0]]*(n-b) + points[T[1]]*(b-a) + points[T[2]]*a;
      }
    }

    return Deltahedron(T_new, new_points);
  }

  // General (k,l) path: unfold to Eisenstein plane, scale, fold back.
  // Then assign 3D coordinates via barycentric interpolation within
  // each original face's scaled triangle in the Eisenstein plane.
  Eisenstein w(k, l);
  int T = w.norm2();  // k^2 + kl + l^2
  double sqrtT = sqrt((double)T);

  // 1. Unfold, scale by w, fold to get new topology
  Unfolding u_orig(static_cast<const Triangulation&>(*this));
  Unfolding gcu(u_orig * w);
  Folding gcf(gcu);
  Triangulation T_new(gcf.fold());

  // 2. Assign 3D coordinates by iterating over original faces.
  //    Each face (u,v,w) has Eisenstein corners in the scaled plane
  //    from gcu.arc_coords. We rasterize the scaled triangle's bounding
  //    box, test containment, and interpolate via barycentric coordinates.
  //
  //    For lattice point e inside triangle (EU, EV, EW):
  //      d1 = EV-EU, d2 = EW-EU, det = T (CCW faces)
  //      sT = d2.b*rel.a - d2.a*rel.b
  //      tT = -d1.b*rel.a + d1.a*rel.b
  //      rT = T - sT - tT
  //      point = (rT*Pu + sT*Pv + tT*Pw) / sqrt(T)
  //
  //    Corner vertices get sqrt(T)*P_original, matching the halma
  //    convention (k = sqrt(k^2) for l=0).
  //    Boundary points get identical coords from adjacent faces.

  vector<coord3d> new_points(T_new.N);

  for(int i = 0; i < (int)triangles.size(); i++){
    const tri_t& tri = triangles[i];
    node_t nu = tri[0], nv = tri[1], nw = tri[2];

    // Scaled Eisenstein corners of this face
    Eisenstein EU = gcu.arc_coords.at({nu, nv}).first;
    Eisenstein EV = gcu.arc_coords.at({nu, nv}).second;
    Eisenstein EW = EU + (EV - EU).nextCCW();

    // Side vectors and determinant (should equal T for CCW faces)
    Eisenstein d1 = EV - EU;
    Eisenstein d2 = EW - EU;

    // Bounding box of the scaled triangle
    int amin = min({EU.first, EV.first, EW.first});
    int amax = max({EU.first, EV.first, EW.first});
    int bmin = min({EU.second, EV.second, EW.second});
    int bmax = max({EU.second, EV.second, EW.second});

    for(int ea = amin; ea <= amax; ea++){
      for(int eb = bmin; eb <= bmax; eb++){
        Eisenstein e(ea, eb);

        // Containment test: all turns >= 0 for CCW triangle
        if(Eisenstein::turn(EU, EV, e) < 0) continue;
        if(Eisenstein::turn(EV, EW, e) < 0) continue;
        if(Eisenstein::turn(EW, EU, e) < 0) continue;

        // Look up the output node ID
        auto it = gcf.final_grid.find(e);
        if(it == gcf.final_grid.end()) continue;
        node_t node_id = it->second;

        // Barycentric weights (integer, sum to T)
        Eisenstein rel = e - EU;
        int sT = d2.second * rel.first - d2.first * rel.second;
        int tT = -d1.second * rel.first + d1.first * rel.second;
        int rT = T - sT - tT;

        new_points[node_id] = (points[nu]*rT + points[nv]*sT + points[nw]*tT) / sqrtT;
      }
    }
  }

  return Deltahedron(T_new, new_points);
}
