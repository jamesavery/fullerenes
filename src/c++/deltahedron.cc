#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/unfold.hh"
#include "fullerenes/buckinverse.hh"
#include <cmath>
#include <climits>
#include <numeric>
#include <queue>
#include <set>
#include <map>
#include <array>
#include <chrono>
#include <stdexcept>

// ============================================================
// VertexHData: per-vertex convexity geometry and derivatives.
//
// Level 1 (basic): h, n_hat, centroid, N_len — for h queries.
// Level 2 (derivs): adds r_perp, per-neighbor De/g/w/Dex/Q —
//   for E_conv gradient, Hv product, and Hessian assembly.
//
// h = (x_v - centroid) · n_hat,  n_hat = N/|N|,
// N = Sum_j (e_j × e_{j+1}),  e_j = x[n_j] - x_v.
// ============================================================
struct VertexHData {
  int vertex;           // vertex index
  int d;                // degree
  double h;             // signed convexity height (h>0 = convex)
  coord3d n_hat;        // unit outward fan normal
  coord3d centroid;     // neighbor centroid
  double N_len;         // |N_fan|
  bool valid;           // false if degenerate (N_len < eps)

  // Derivative-level data (populated by compute_derivs)
  coord3d r_perp;       // (x_v - centroid) - h*n_hat
  struct Nbr {
    int id;             // neighbor vertex index
    coord3d De;         // e_{j-1} - e_{j+1}
    coord3d g;          // dh/dx_{n_j} = -n_hat/d + r_perp × De / N_len
    coord3d w;          // r_perp × De / N_len (normal-derivative part of g)
    matrix3d Dex;       // [De]×  (cross-product matrix)
    matrix3d Q;         // P · [De]× / N_len  (d(n_hat)/dx_{n_j})
  };
  vector<Nbr> nb;
  bool has_derivs = false;

  // Compute basic h data for vertex v.
  static VertexHData compute_h(const DeltahedronView<double>& D, std::span<const coord3d> x, int v) {
    VertexHData vd;
    vd.vertex = v;
    vd.d = D.degree(v);
    vd.has_derivs = false;

    // Neighbor centroid
    vd.centroid = coord3d(0,0,0);
    for(int j = 0; j < vd.d; j++) vd.centroid += x[D[v][j]];
    vd.centroid /= (double)vd.d;

    // Fan normal (unnormalized)
    coord3d N_fan(0,0,0);
    for(int j = 0; j < vd.d; j++){
      coord3d e1 = x[D[v][j]] - x[v];
      coord3d e2 = x[D[v][(j+1)%vd.d]] - x[v];
      N_fan += e1.cross(e2);
    }
    vd.N_len = N_fan.norm();
    vd.valid = (vd.N_len > 1e-15);
    if(!vd.valid){ vd.h = 0; vd.n_hat = coord3d(0,0,0); return vd; }

    vd.n_hat = N_fan / vd.N_len;
    vd.h = (x[v] - vd.centroid).dot(vd.n_hat);
    return vd;
  }

  // Compute full derivative data (extends basic h data).
  static VertexHData compute_derivs(const DeltahedronView<double>& D, std::span<const coord3d> x, int v) {
    VertexHData vd = compute_h(D, x, v);
    if(!vd.valid) return vd;
    vd.has_derivs = true;

    vd.r_perp = (x[v] - vd.centroid) - vd.n_hat * vd.h;
    matrix3d P = matrix3d::unit_matrix() - vd.n_hat.outer(vd.n_hat);

    vd.nb.resize(vd.d);
    for(int j = 0; j < vd.d; j++){
      vd.nb[j].id = D[v][j];
      coord3d ej_prev = x[D[v][(j+vd.d-1)%vd.d]] - x[v];
      coord3d ej_next = x[D[v][(j+1)%vd.d]]       - x[v];
      vd.nb[j].De  = ej_prev - ej_next;
      vd.nb[j].g   = vd.n_hat * (-1.0/vd.d) + vd.r_perp.cross(vd.nb[j].De) / vd.N_len;
      vd.nb[j].w   = vd.r_perp.cross(vd.nb[j].De) / vd.N_len;
      vd.nb[j].Dex = matrix3d::cross_matrix(vd.nb[j].De);
      vd.nb[j].Q   = P * vd.nb[j].Dex * (1.0/vd.N_len);
    }
    return vd;
  }

  // --- Operations on derivative data ---

  // Gather: (dh/dx) · v  where dh/dx_v = n_hat, dh/dx_{n_j} = g_j.
  double dhdx_dot(const vector<coord3d>& v) const {
    double s = n_hat.dot(v[vertex]);
    for(int j = 0; j < d; j++) s += nb[j].g.dot(v[nb[j].id]);
    return s;
  }

  // Scatter: out[v] += alpha * n_hat;  out[n_j] += alpha * g_j.
  // (Gradient of a scalar function f(h): df/dx = dfdh * dh/dx.)
  void scatter_dhdx(double alpha, vector<coord3d>& out) const {
    out[vertex] = out[vertex] + n_hat * alpha;
    for(int j = 0; j < d; j++)
      out[nb[j].id] = out[nb[j].id] + nb[j].g * alpha;
  }

  // Correction-term Hv: Hv += dEdh * d²h/dx² * v.
  // d²h/dx² has blocks H_{vv}=0, H_{v,j}=Q_j, H_{j,k}=7-term formula.
  void scatter_d2h_hv(double dEdh, const vector<coord3d>& v, vector<coord3d>& Hv) const {
    if(fabs(dEdh) < 1e-15) return;
    matrix3d R = matrix3d::cross_matrix(r_perp);

    // H_{v,j} blocks: Hv[v] += Q_j * dEdh * v[n_j], Hv[n_j] += Q_j^T * dEdh * v[v]
    for(int j = 0; j < d; j++){
      matrix3d Hvj = nb[j].Q * dEdh;
      Hv[vertex]   = Hv[vertex]   + Hvj * v[nb[j].id];
      Hv[nb[j].id] = Hv[nb[j].id] + Hvj.transpose() * v[vertex];
    }

    // H_{j,k} blocks (7-term formula)
    for(int j = 0; j < d; j++){
      for(int k = 0; k < d; k++){
        matrix3d Hjk = nb[k].Q * (-1.0/d)                                        // A: -(1/d) Q_k
                     + nb[j].Dex * (1.0 / (d * N_len))                            // B: [De_j]× / (d·|N|)
                     + nb[j].De.cross(n_hat).outer(nb[k].g) * (1.0/N_len)         // C: (De_j × n_hat)⊗g_k / |N|
                     + nb[j].Dex * nb[k].Q * (h / N_len);                         // D: h·[De_j]×·Q_k / |N|
        if(k == (j+d-1)%d) Hjk += R * ( 1.0/N_len);                              // E: +[r_perp]× for k=j-1
        if(k == (j+1)%d)   Hjk += R * (-1.0/N_len);                              // F: -[r_perp]× for k=j+1
        Hjk += nb[j].w.outer(n_hat.cross(nb[k].De)) * (-1.0/N_len);              // G: -w_j ⊗ (n_hat × De_k)/|N|

        Hv[nb[j].id] = Hv[nb[j].id] + (Hjk * dEdh) * v[nb[k].id];
      }
    }
  }

  // Correction-term Hessian blocks: calls add_block(row_vertex, col_vertex, 3x3 matrix).
  // Same formula as scatter_d2h_hv but assembles explicit blocks.
  template<typename AddBlockFn>
  void scatter_d2h_blocks(double dEdh, AddBlockFn&& add_block) const {
    if(fabs(dEdh) < 1e-15) return;
    matrix3d R = matrix3d::cross_matrix(r_perp);

    // H_{v,j} blocks
    for(int j = 0; j < d; j++){
      matrix3d Hvj = nb[j].Q * dEdh;
      add_block(vertex, nb[j].id, Hvj);
      add_block(nb[j].id, vertex, Hvj.transpose());
    }

    // H_{j,k} blocks
    for(int j = 0; j < d; j++){
      for(int k = 0; k < d; k++){
        matrix3d Hjk = nb[k].Q * (-1.0/d)
                     + nb[j].Dex * (1.0 / (d * N_len))
                     + nb[j].De.cross(n_hat).outer(nb[k].g) * (1.0/N_len)
                     + nb[j].Dex * nb[k].Q * (h / N_len);
        if(k == (j+d-1)%d) Hjk += R * ( 1.0/N_len);
        if(k == (j+1)%d)   Hjk += R * (-1.0/N_len);
        Hjk += nb[j].w.outer(n_hat.cross(nb[k].De)) * (-1.0/N_len);

        add_block(nb[j].id, nb[k].id, Hjk * dEdh);
      }
    }
  }

  // Rank-1 Hessian blocks: d²E/dh² * (dh/dx)⊗(dh/dx).
  // Calls add_block(row_vertex, col_vertex, 3x3 matrix).
  template<typename AddBlockFn>
  void scatter_rank1_blocks(double d2Edh2, AddBlockFn&& add_block) const {
    if(fabs(d2Edh2) < 1e-15) return;
    // v-v block
    add_block(vertex, vertex, n_hat.outer(n_hat) * d2Edh2);
    // v-j and j-v blocks
    for(int j = 0; j < d; j++){
      matrix3d Bvj = n_hat.outer(nb[j].g) * d2Edh2;
      add_block(vertex, nb[j].id, Bvj);
      add_block(nb[j].id, vertex, Bvj.transpose());
    }
    // j-k blocks
    for(int j = 0; j < d; j++)
      for(int k = 0; k < d; k++)
        add_block(nb[j].id, nb[k].id, nb[j].g.outer(nb[k].g) * d2Edh2);
  }
};

// ============================================================
// EdgeBondData: per-edge bond geometry for E_bond.
//
// E_bond = (k/2) * (r - L)^2  per edge.
// Hessian: M = k[(1-L/r)I + (L/r^3)d⊗d].
// ============================================================
struct EdgeBondData {
  int u, v;
  coord3d d;          // x[u] - x[v]
  double r, dev;      // |d|, r - L
  bool valid;

  static EdgeBondData compute(std::span<const coord3d> x, int u, int v, double L) {
    EdgeBondData ed;
    ed.u = u; ed.v = v;
    ed.d = x[u] - x[v];
    ed.r = ed.d.norm();
    ed.valid = (ed.r > 1e-15);
    ed.dev = ed.r - L;
    return ed;
  }

  double energy(double k) const { return 0.5 * k * dev * dev; }

  // out[u] += k*(dev/r)*d,  out[v] -= same
  void scatter_gradient(double k, vector<coord3d>& out) const {
    coord3d g = d * (k * dev / r);
    out[u] = out[u] + g;
    out[v] = out[v] - g;
  }

  // Hv[u] += M*(v[u]-v[v]),  Hv[v] -= same
  // M*dv = k*(1-L/r)*dv + k*(L/r^3)*(d.dv)*d
  void scatter_hv(double k, double L, const vector<coord3d>& v, vector<coord3d>& Hv) const {
    coord3d dv = v[u] - v[this->v];
    coord3d Mdv = dv * (k * (1 - L/r)) + d * (k * L / (r*r*r) * d.dot(dv));
    Hv[u] = Hv[u] + Mdv;
    Hv[this->v] = Hv[this->v] - Mdv;
  }

  // Hessian blocks: H(u,u) = H(v,v) = +M,  H(u,v) = H(v,u) = -M
  template<typename F>
  void scatter_blocks(double k, double L, F&& add_block) const {
    matrix3d M = matrix3d::unit_matrix() * (k * (1 - L/r))
               + d.outer(d) * (k * L / (r*r*r));
    add_block(u, v,       -M);
    add_block(this->v, u, -M);
    add_block(u, u,        M);
    add_block(this->v, this->v, M);
  }
};

// ============================================================
// CornerAngleData: per-triangle-corner angle geometry.
//
// Angle theta at vertex b between arms va = x[a]-x[b], vc = x[d]-x[b].
// E_angle = (k/2) * (theta - pi/3)^2  per corner.
// Shared by E_angle (per triangle) and E_curv (per vertex fan).
// ============================================================
struct CornerAngleData {
  int a, b, d;              // arm1 end, apex, arm2 end
  double theta, dev;        // angle, deviation from pi/3
  coord3d p, q;             // dtheta/d(va), dtheta/d(vc)
  double ra, rc, S, alpha;  // arm lengths, sin(theta), 1 - dev*cos/sin
  coord3d ua, uc;           // unit arm vectors
  bool valid;

  static CornerAngleData compute(std::span<const coord3d> x, int a, int b, int d) {
    CornerAngleData ca;
    ca.a = a; ca.b = b; ca.d = d;

    coord3d va = x[a] - x[b], vc = x[d] - x[b];
    ca.ra = va.norm(); ca.rc = vc.norm();
    ca.valid = (ca.ra > 1e-15 && ca.rc > 1e-15);
    if(!ca.valid) return ca;

    ca.ua = va / ca.ra;
    ca.uc = vc / ca.rc;
    double C = max(-1.0, min(1.0, ca.ua.dot(ca.uc)));
    ca.theta = acos(C);
    ca.S = sin(ca.theta);
    if(ca.S < 1e-10){ ca.valid = false; return ca; }

    ca.dev = ca.theta - M_PI / 3.0;
    ca.alpha = 1.0 - ca.dev * C / ca.S;

    coord3d::dangle(va, vc, ca.p, ca.q);
    return ca;
  }

  // out[a] += w*p,  out[d] += w*q,  out[b] -= w*(p+q)
  void scatter_gradient(double w, vector<coord3d>& out) const {
    out[a] = out[a] + p * w;
    out[d] = out[d] + q * w;
    out[b] = out[b] - (p + q) * w;
  }

  // Hv scatter with weight k.
  // Builds arm-space blocks Haa/Hcc/Hac, scatters k*(block*v) to {a,d,b}.
  // Includes both rank-1 (k*alpha*p⊗p etc) and correction terms.
  void scatter_hv(double k, const vector<coord3d>& v, vector<coord3d>& Hv) const {
    matrix3d sym_ac = ua.outer(uc) + uc.outer(ua);

    matrix3d Haa = p.outer(p) * (k * alpha)
      + (sym_ac + matrix3d::unit_matrix() * (ua.dot(uc))
         - ua.outer(ua) * (3*ua.dot(uc))) * (k * dev / (ra*ra * S));

    matrix3d Hcc = q.outer(q) * (k * alpha)
      + (sym_ac + matrix3d::unit_matrix() * (ua.dot(uc))
         - uc.outer(uc) * (3*ua.dot(uc))) * (k * dev / (rc*rc * S));

    matrix3d Hac = p.outer(q) * (k * alpha)
      + (ua.outer(ua) + uc.outer(uc) - ua.outer(uc) * (ua.dot(uc))
         - matrix3d::unit_matrix()) * (k * dev / (ra*rc * S));

    matrix3d Hab = Haa + Hac;
    matrix3d Hcb = Hcc + Hac.transpose();

    Hv[a] = Hv[a] + Haa * v[a] + Hac * v[d] + (-Hab) * v[b];
    Hv[d] = Hv[d] + Hac.transpose() * v[a] + Hcc * v[d] + (-Hcb) * v[b];
    Hv[b] = Hv[b] + (-Hab.transpose()) * v[a] + (-Hcb.transpose()) * v[d] + (Hab + Hcb) * v[b];
  }

  // Hessian blocks via add_block callback, weight k.
  template<typename F>
  void scatter_blocks(double k, F&& add_block) const {
    matrix3d sym_ac = ua.outer(uc) + uc.outer(ua);
    double C = ua.dot(uc);
    matrix3d I = matrix3d::unit_matrix();

    matrix3d Haa = p.outer(p) * (k * alpha)
      + (sym_ac + I * C - ua.outer(ua) * (3*C)) * (k * dev / (ra*ra * S));

    matrix3d Hcc = q.outer(q) * (k * alpha)
      + (sym_ac + I * C - uc.outer(uc) * (3*C)) * (k * dev / (rc*rc * S));

    matrix3d Hac = p.outer(q) * (k * alpha)
      + (ua.outer(ua) + uc.outer(uc) - ua.outer(uc) * C - I) * (k * dev / (ra*rc * S));

    add_block(a, a, Haa);
    add_block(a, d, Hac);
    add_block(d, a, Hac.transpose());
    add_block(d, d, Hcc);

    matrix3d Hab = Haa + Hac;
    matrix3d Hcb = Hcc + Hac.transpose();
    add_block(a, b, -Hab);
    add_block(b, a, -Hab.transpose());
    add_block(d, b, -Hcb);
    add_block(b, d, -Hcb.transpose());
    add_block(b, b, Hab + Hcb);
  }
};

Deltahedron::Deltahedron(const TriangulationView& T, const vector<coord3d>& pts)
  : base_t(static_cast<const GraphView&>(T), std::vector<coord3d>(pts))
{
  assert((int)points.size() == N);
}

Deltahedron::Deltahedron(const TriangulationView& T, std::span<const coord3d> pts)
  : base_t(static_cast<const GraphView&>(T), std::vector<coord3d>(pts.begin(), pts.end()))
{
  assert((int)points.size() == N);
}

Deltahedron::Deltahedron(const Polyhedron& P)
  : base_t(static_cast<const GraphView&>(P), std::vector<coord3d>(P.points.begin(), P.points.end()))
{
  assert(P.is_triangulation());
}

template<>
vector<face_t> DeltahedronView<double>::compute_dual_faces() const {
  auto tris = triangles();
  vector<face_t> faces(tris.size());
  for(size_t i = 0; i < tris.size(); i++)
    faces[i] = face_t(tris[i]);
  return faces;
}

template<>
double DeltahedronView<double>::max_angle_relerr() const {
  double max_re = 0;
  const double target = M_PI / 3.0;
  for (const auto& t : triangles()) {
    for (int c = 0; c < 3; c++) {
      coord3d ea = points[t[(c+1)%3]] - points[t[c]];
      coord3d eb = points[t[(c+2)%3]] - points[t[c]];
      double na = ea.norm(), nb = eb.norm();
      if (na < 1e-15 || nb < 1e-15) return 1.0;  // degenerate → max error
      double cos_th = ea.dot(eb) / (na * nb);
      cos_th = max(-1.0, min(1.0, cos_th));
      double th = acos(cos_th);
      double re = fabs(th - target) / target;
      if (re > max_re) max_re = re;
    }
  }
  return max_re;
}

template<>
int DeltahedronView<double>::count_concave() const {
  int n_concave = 0;
  for(int v = 0; v < N; v++){
    auto vd = VertexHData::compute_h(*this, points, v);
    if(vd.valid && vd.h < 0) n_concave++;
  }
  return n_concave;
}

template<>
void DeltahedronView<double>::smooth(double q) {
  vector<coord3d> new_points(N);
  for(node_t u = 0; u < N; u++){
    coord3d avg;
    for(node_t v : (*this)[u]) avg += points[v];
    avg /= degree(u);
    new_points[u] = points[u]*(1.0-q) + avg*q;
  }
  std::copy(new_points.begin(), new_points.end(), points.begin());
}

// Taubin lambda|mu low-pass: a forward umbrella step (lambda>0) then a shrink-
// correcting step (mu<-lambda), per pass -- the standard volume-preserving filter
// that suppresses high-frequency mesh noise without the shrinkage of plain smoothing.
template<>
void DeltahedronView<double>::taubin_smooth(int iters, double lambda, double mu) {
  for(int i = 0; i < iters; i++){ smooth(lambda); smooth(mu); }
}

// Discrete Gaussian curvature: the angle defect K(v) = 2*pi - sum of the triangle
// corner angles at v. A degenerate triangle (a zero-length incident edge) has an
// undefined corner angle and is skipped, contributing nothing to the defect. The
// tris-taking overload lets a caller that already has the triangulation avoid
// recomputing it (triangles() is O(E)).
template<>
std::vector<double> DeltahedronView<double>::angle_defects(const std::vector<tri_t>& tris) const {
  std::vector<double> K(N, 2.0*M_PI);
  for(const tri_t& t : tris)
    for(int c = 0; c < 3; c++){
      coord3d e1 = points[t[(c+1)%3]] - points[t[c]], e2 = points[t[(c+2)%3]] - points[t[c]];
      if(e1.norm2() > 1e-300 && e2.norm2() > 1e-300) K[t[c]] -= coord3d::angle(e1, e2);
    }
  return K;
}

template<>
std::vector<double> DeltahedronView<double>::angle_defects() const { return angle_defects(triangles()); }

template<>
Deltahedron DeltahedronView<double>::halma_transform(int m) const {
    // Halma path: direct subdivision via face grids, preserves node IDs
    int n = m + 1;  // = k in GC(k,0) terminology

    vector<map<edge_t,node_t>> face_grids;
    Triangulation T_new = TriangulationView::halma_transform(m, &face_grids);

    // Barycentric interpolation in each face grid.
    // Grid corners at {0,0}, {0,n}, {n,n} map to T[0], T[1], T[2].
    // point({a,b}) = (n-b)*P[T[0]] + (b-a)*P[T[1]] + a*P[T[2]]
    // Weights sum to n=k, so corner vertices get k*P_original.
    vector<coord3d> new_points(T_new.N);
    auto tris = triangles();

    for(int i = 0; i < (int)tris.size(); i++){
      const face_t& T = tris[i];
      const auto& grid = face_grids[i];

      for(const auto& [ab, node_id] : grid){
        int a = ab.first, b = ab.second;
        new_points[node_id] = points[T[0]]*(n-b) + points[T[1]]*(b-a) + points[T[2]]*a;
      }
    }

    return Deltahedron(T_new, new_points);  
}

template<>
Deltahedron DeltahedronView<double>::GCtransform(unsigned k, unsigned l) const {
  if(l==0 || k==0) return halma_transform(max(k,l) - 1);

  // General (k,l) path: unfold to Eisenstein plane, scale, fold back.
  // Then assign 3D coordinates via barycentric interpolation within
  // each original face's scaled triangle in the Eisenstein plane.
  Eisenstein w(k, l);
  int T = w.norm2();  // k^2 + kl + l^2
  double sqrtT = sqrt((double)T);

  // 1. Unfold, scale by w, fold to get new topology
  Unfolding u_orig(static_cast<const Triangulation&>(*this));
  Unfolding gcu(u_orig * w);
  Folding   gcf(gcu);
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
  auto tris = triangles();

  for(int i = 0; i < (int)tris.size(); i++){
    const tri_t& tri = tris[i];
    node_t nu = tri[0], nv = tri[1], nw = tri[2];

    // Scaled Eisenstein corners of this face. nu,nv are ORIGINAL labels (from
    // triangles()); gcu.arc_coords is keyed by the cones-first labels, so map through
    // gcu.cone_perm (which survives the *w copy). final_grid (read below) is already
    // back in original labels, matching T_new.
    auto [EU,EV] = gcu.arc_coords.at({gcu.cone_perm[nu], gcu.cone_perm[nv]});
    Eisenstein EW = EU + (EV - EU).nextCCW();

    // Side vectors and determinant (should equal T for CCW faces)
    Eisenstein d1 = EV - EU, d2 = EW - EU;

    // Bounding box of the scaled triangle
    int amin = min({EU.first, EV.first, EW.first});
    int amax = max({EU.first, EV.first, EW.first});
    int bmin = min({EU.second, EV.second, EW.second});
    int bmax = max({EU.second, EV.second, EW.second});

    for(int ea = amin; ea <= amax; ea++){
      for(int eb = bmin; eb <= bmax; eb++){
        Eisenstein e = {ea, eb};

        // Containment test: all turns >= 0 for CCW triangle
        if(Eisenstein::turn(EU, EV, e) < 0) continue;
        if(Eisenstein::turn(EV, EW, e) < 0) continue;
        if(Eisenstein::turn(EW, EU, e) < 0) continue;

        // Look up the output node ID for each Eisenstein point e inside the triangle T
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

// =====================================================================
// Incremental geometry from extension path
// =====================================================================

using buckinverse::ExtensionPath;
using buckinverse::ExtensionStep;
using buckinverse::ExpKind;
using buckinverse::ReducibleDual;

// Thomas algorithm for tridiagonal system with coord3d values.
// diag[i] = diagonal entries, rhs[i] = right-hand sides.
// Sub- and super-diagonal entries are all -1.
// Returns solution in-place in rhs.
static void solveTridiagonal(const vector<double>& diag, vector<coord3d>& rhs) {
    int n = (int)rhs.size();
    if (n == 0) return;
    if (n == 1) { rhs[0] = rhs[0] / diag[0]; return; }

    // Forward sweep
    vector<double> d(diag);
    vector<double> c(n - 1, -1.0);  // super-diagonal

    for (int i = 1; i < n; i++) {
        double w = -1.0 / d[i - 1];  // sub-diagonal / d'[i-1]
        d[i] -= w * c[i - 1];
        rhs[i] = rhs[i] - rhs[i - 1] * w;
    }

    // Back substitution
    rhs[n - 1] = rhs[n - 1] / d[n - 1];
    for (int i = n - 2; i >= 0; i--)
        rhs[i] = (rhs[i] - rhs[i + 1] * c[i]) / d[i];
}

// For F-ring expansion: translate + rotate tp-side to make room for the new ring.
// In an equilateral antiprism, consecutive rings are separated by height h along
// the axis AND rotated by pi/5 (36 degrees). Before this call, path and tp are
// adjacent rings. After, they're separated by 2h with 2*(pi/5) rotation, so the
// strip placed at the midpoint has correct edge lengths to both sides.
static void shiftForFRing(const ReducibleDual& rd, const ExtensionStep& step,
                          vector<coord3d>& points) {
    const auto& path = step.path;  // ring that stays in place
    const auto& tp = step.tp;      // outer vertices that get shifted

    // Compute axis: path centroid -> tp centroid
    coord3d c_path(0,0,0), c_tp(0,0,0);
    for (int i = 0; i < 5; i++) {
        c_path += points[path[i]];
        c_tp   += points[tp[i]];
    }
    c_path = c_path / 5.0;
    c_tp   = c_tp   / 5.0;

    coord3d shift = c_tp - c_path;  // one ring spacing along axis
    double h = shift.norm();
    coord3d axis = shift / h;       // unit axis direction

    // Determine rotation direction from existing geometry.
    // In the antiprism, tp[0] connects to path[0] and path[1] (or path[4]).
    // The bisector of path[0]/path[1] should be near tp[0]'s current angle.
    // After inserting a ring, tp needs to rotate by another pi/5 in the same
    // direction that path->tp already rotates.
    // Detect direction: project path[0] and tp[0] onto the plane perpendicular
    // to axis, compute the signed angle from path[0] to tp[0].
    coord3d r_p0 = points[path[0]] - c_path;
    r_p0 = r_p0 - axis * r_p0.dot(axis);  // radial component
    coord3d r_t0 = points[tp[0]] - c_tp;
    r_t0 = r_t0 - axis * r_t0.dot(axis);  // radial component
    // Signed angle from r_p0 to r_t0 around axis
    double cross_z = (r_p0.cross(r_t0)).dot(axis);
    double theta = (cross_z > 0) ? -M_PI / 5.0 : M_PI / 5.0;  // rotate FURTHER in same direction
    double ct = cos(theta), st = sin(theta);

    // BFS from tp vertices to find all vertices on the tp-side.
    // Stop at path vertices (they form the boundary and don't move).
    set<int> path_set(path.begin(), path.end());
    set<int> visited;
    queue<int> q;
    for (int v : tp) {
        if (!visited.count(v)) {
            visited.insert(v);
            q.push(v);
        }
    }
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int slot = 0; slot < ReducibleDual::D_MAX; slot++) {
            if (!(rd.V[u].active & (1 << slot))) continue;
            int nb = rd.V[u].nbr[slot];
            if (path_set.count(nb) || visited.count(nb)) continue;
            visited.insert(nb);
            q.push(nb);
        }
    }

    // Apply translation + rotation (Rodrigues' formula) to all tp-side vertices.
    // Rotation center is c_tp (rotate around axis through tp centroid), then shift.
    for (int v : visited) {
        // Rotate around axis through c_tp by theta
        coord3d p = points[v] - c_tp;
        coord3d p_rot = p * ct + axis.cross(p) * st + axis * axis.dot(p) * (1 - ct);
        // Translate by one ring spacing
        points[v] = c_tp + p_rot + shift;
    }
}

// Compute strip vertex coordinates for one expansion step.
// points[] is indexed by ReducibleDual vertex IDs (full graph numbering).
// Places each strip vertex at the centroid of its 4 non-strip neighbors.
static void computeStripCoords(const ExtensionStep& step, vector<coord3d>& points) {
    const auto& strip = step.strip;
    const auto& path = step.path;
    const auto& tp = step.tp;
    int n = (int)strip.size();

    if (step.kind.type == ExpKind::L_type) {
        for (int j = 0; j < n; j++)
            points[strip[j]] = (points[path[j]] + points[path[j + 1]]
                              + points[tp[j]] + points[tp[j + 1]]) * 0.25;

    } else if (step.kind.type == ExpKind::B_type) {
        assert(step.kind.i == 0 && step.kind.j == 0 && n == 3);
        points[strip[0]] = (points[path[0]] + points[path[1]]
                          + points[tp[0]] + points[tp[1]]) * 0.25;
        points[strip[1]] = (points[path[1]] + points[path[2]]
                          + points[path[3]] + points[tp[1]]) * 0.25;
        points[strip[2]] = (points[path[3]] + points[path[4]]
                          + points[tp[1]] + points[tp[2]]) * 0.25;

    } else {
        // F-ring: antiprism placement on the local cylinder.
        // For each strip[i], place at the bisector angle of path[i] and
        // path[(i+1)%5], at the cylinder radius and midpoint height.
        // This correctly handles any index offset between path and tp arrays.
        assert(n == 5);

        // 1. Ring centroids
        coord3d c_path(0,0,0), c_tp(0,0,0);
        for (int i = 0; i < 5; i++) {
            c_path += points[path[i]];
            c_tp   += points[tp[i]];
        }
        c_path = c_path / 5.0;
        c_tp   = c_tp   / 5.0;

        // 2. Cylinder axis (path center toward tp center).  Coincident ring
        // centroids leave the axis (and the whole antiprism frame) undefined;
        // that means the upstream geometry is invalid, so fail loudly rather
        // than smooth over it (a Laplacian fallback would only hide the bug).
        coord3d axis = c_tp - c_path;
        double axis_len = axis.norm();
        if (axis_len < 1e-9)
            throw std::runtime_error(
                "shiftForFRing: path and tp ring centroids coincide "
                "(degenerate cylinder axis; upstream geometry invalid)");
        axis = axis / axis_len;

        // 3. Cylinder radius from path ring
        double r = 0;
        for (int i = 0; i < 5; i++) {
            coord3d v = points[path[i]] - c_path;
            coord3d radial = v - axis * v.dot(axis);
            r += radial.norm();
        }
        r /= 5.0;

        // 4. Local coordinate frame perpendicular to axis
        coord3d v0 = points[path[0]] - c_path;
        coord3d e1 = v0 - axis * v0.dot(axis);
        e1 = e1 / e1.norm();
        coord3d e2 = axis.cross(e1);

        // 5. Strip center at midpoint between ring centroids
        coord3d c_strip = (c_path + c_tp) * 0.5;

        // 6. Place each strip[i] at the bisector angle of path[i] and path[(i+1)%5]
        for (int i = 0; i < 5; i++) {
            int ip1 = (i + 1) % 5;
            // Compute radial directions for path[i] and path[ip1]
            coord3d vi   = points[path[i]] - c_path;
            coord3d vip1 = points[path[ip1]] - c_path;
            coord3d ri   = vi - axis * vi.dot(axis);
            coord3d rip1 = vip1 - axis * vip1.dot(axis);
            // Bisector direction (sum of normalized radials)
            coord3d bisector = ri / ri.norm() + rip1 / rip1.norm();
            double blen = bisector.norm();
            if (blen < 1e-12) {
                // path[i] and path[ip1] are diametrically opposite; use perpendicular
                bisector = axis.cross(ri);
                blen = bisector.norm();
            }
            bisector = bisector / blen;
            points[strip[i]] = c_strip + bisector * r;
        }
        return;
    }
}

// Lift strip vertices outward to match the surface height of their neighbours.
// After computeStripCoords + rd.expand, each strip vertex sits at the quad
// centroid which is inside the surface (chord vs arc). We compute the fan
// normal (points outward from CCW convention) and push each strip vertex
// outward along it until its h value (height above neighbour centroid along
// the fan normal) matches the average h of its non-strip neighbours.
//
// TODO: Cleanup (cleaner code, factor into higher-level abstractions  + add clearer descriptions of the method).
static void liftStripToSurface(const ExtensionStep& step,
                                const ReducibleDual& rd,
                                vector<coord3d>& points) {
    // Collect strip vertex IDs into a set for fast lookup.
    set<int> strip_set(step.strip.begin(), step.strip.end());

    for (int s : step.strip) {
        // Gather CCW-ordered neighbours from active bitmask
        vector<int> nbrs;
        uint8_t m = rd.V[s].active;
        for (; m; m &= m - 1)
            nbrs.push_back(rd.V[s].nbr[__builtin_ctz(m)]);
        int d = (int)nbrs.size();
        if (d < 3) continue;

        // Compute h for each non-strip neighbour: its height above ITS own
        // neighbour centroid along ITS own fan normal. This tells us how far
        // "above" the local surface each boundary vertex sits.
        double h_target = 0;
        int n_boundary = 0;
        for (int nb : nbrs) {
            if (strip_set.count(nb)) continue;  // skip other strip vertices

            // Gather nb's neighbours
            vector<int> nb_nbrs;
            uint8_t mb = rd.V[nb].active;
            for (; mb; mb &= mb - 1)
                nb_nbrs.push_back(rd.V[nb].nbr[__builtin_ctz(mb)]);
            int db = (int)nb_nbrs.size();
            if (db < 3) continue;

            coord3d c_nb(0,0,0);
            for (int nn : nb_nbrs) c_nb += points[nn];
            c_nb /= (double)db;

            coord3d nf(0,0,0);
            for (int j = 0; j < db; j++) {
                coord3d e1 = points[nb_nbrs[j]] - points[nb];
                coord3d e2 = points[nb_nbrs[(j+1)%db]] - points[nb];
                nf += e1.cross(e2);
            }
            double nl = nf.norm();
            if (nl < 1e-15) continue;

            double h_nb = (points[nb] - c_nb).dot(nf / nl);
            h_target += h_nb;
            n_boundary++;
        }
        if (n_boundary == 0) continue;
        h_target /= n_boundary;

        // Now compute h for the strip vertex itself
        coord3d centroid(0,0,0);
        for (int nb : nbrs) centroid += points[nb];
        centroid /= (double)d;

        coord3d n_fan(0,0,0);
        for (int j = 0; j < d; j++) {
            coord3d e1 = points[nbrs[j]] - points[s];
            coord3d e2 = points[nbrs[(j+1)%d]] - points[s];
            n_fan += e1.cross(e2);
        }
        double n_len = n_fan.norm();
        if (n_len < 1e-15) continue;
        coord3d n_hat = n_fan / n_len;

        double h_current = (points[s] - centroid).dot(n_hat);

        // Push outward along the fan normal to match target h
        points[s] = points[s] + n_hat * (h_target - h_current);
    }
}


// TODO: Move seed data to its own file that all can use.
// =====================================================================
// Precomputed seed dual geometry (force-field optimized face centroids)
// Generated by gen_seeds from IsomerDB + FullereneGraph::optimized_geometry
// =====================================================================

// C20 dual (icosahedron): 12 vertices, all degree-5
static const Graph C20_seed_neighbours = {
    {5, 3, 2, 1, 4},
    {11, 4, 0, 2, 10},
    {10, 1, 0, 3, 8},
    {8, 2, 0, 5, 6},
    {5, 0, 1, 11, 7},
    {6, 3, 0, 4, 7},
    {8, 3, 5, 7, 9},
    {6, 5, 4, 11, 9},
    {10, 2, 3, 6, 9},
    {10, 8, 6, 7, 11},
    {11, 1, 2, 8, 9},
    {9, 7, 4, 1, 10}
};
static const vector<coord3d> C20_seed_points = {
    {-7.08588823223712971e-01, -5.36369371716198806e-01, -4.58486006315860395e-01},
    {-2.63305120431595041e-01, -8.59841516687480545e-01, 4.37427275138782212e-01},
    {2.97738250992196696e-01, -8.41051885497480334e-01, -4.51645105251351542e-01},
    {9.38550249517450946e-03, 8.54836290281664989e-03, -9.99919236306042403e-01},
    {-8.98401741631351780e-01, -2.18538989986770127e-02, 4.38630093623971318e-01},
    {-7.29869668249228720e-01, 5.14840562074550201e-01, -4.49698904060023441e-01},
    {2.63305120431593931e-01, 8.59841516687480989e-01, -4.37427275138781824e-01},
    {-2.97738250992197417e-01, 8.41051885497480889e-01, 4.51645105251350709e-01},
    {8.98401741631351780e-01, 2.18538989986784421e-02, -4.38630093623970596e-01},
    {7.08588823223713082e-01, 5.36369371716198140e-01, 4.58486006315860617e-01},
    {7.29869668249228165e-01, -5.14840562074551311e-01, 4.49698904060022442e-01},
    {-9.38550249517353975e-03, -8.54836290281626304e-03, 9.99919236306042736e-01}
};

// C28 dual (Td): 16 vertices, 12 deg-5 + 4 deg-6
static const Graph C28_seed_neighbours = {
    {6, 4, 3, 2, 1, 5},
    {9, 5, 0, 2, 8},
    {8, 1, 0, 3, 7},
    {7, 2, 0, 4, 15},
    {15, 3, 0, 6, 14},
    {6, 0, 1, 9, 12},
    {14, 4, 0, 5, 12},
    {8, 2, 3, 15, 10},
    {11, 9, 1, 2, 7, 10},
    {12, 5, 1, 8, 11},
    {11, 8, 7, 15, 13},
    {12, 9, 8, 10, 13},
    {14, 6, 5, 9, 11, 13},
    {14, 12, 11, 10, 15},
    {15, 4, 6, 12, 13},
    {13, 10, 7, 3, 4, 14}
};
static const vector<coord3d> C28_seed_points = {
    {3.34261161509050497e-01, -1.63608600073426391e+00, 1.46018963749153186e-01},
    {1.12114359230029703e+00, -1.10236932111558161e+00, -4.45486668492173443e-01},
    {7.72128376690640161e-01, -6.23262654666507099e-01, -1.39551906183534102e+00},
    {-3.39369586004024670e-01, -5.48090254453619141e-01, -1.28221523958660777e+00},
    {-4.20658530314941270e-01, -8.94700466192063892e-01, -2.20532690226496852e-01},
    {-2.74654620369771785e-01, -1.04803812254568718e+00, 8.79052135236750409e-01},
    {6.60115089027072499e-01, -6.07977804813272571e-01, 4.47231716409980551e-01},
    {1.95501520122909694e-01, -1.47531034525772098e+00, -9.53441502794699414e-01},
    {2.83252171975811096e-01, 3.81672051566855852e-01, -1.32464833702157248e+00},
    {8.91710597131574856e-01, -6.58364910568587109e-03, -4.68522048021769888e-01},
    {3.94355331258453290e-01, 9.93890241770606364e-01, -3.93648153055207339e-01},
    {8.00289010343711915e-01, 5.00215199814869060e-01, 5.25818631607869458e-01},
    {4.56949224383211006e-02, -1.59385214877641561e-02, 1.17241701193989933e+00},
    {-2.76110774878330190e-01, 8.05949691469507123e-01, 4.83318161599257334e-01},
    {-8.44848972266442844e-01, -1.58317454188156764e-01, 5.08657885380037755e-01},
    {-3.93953111470164208e-01, 2.00700090168516948e-01, -4.51402369252568603e-01}
};

// C30 dual (D5h): 17 vertices, 12 deg-5 + 5 deg-6
static const Graph C30_seed_neighbours = {
    {11, 3, 2, 1, 4},
    {5, 4, 0, 2, 6},
    {6, 1, 0, 3, 7},
    {7, 2, 0, 11, 9},
    {11, 0, 1, 5, 8},
    {16, 8, 4, 1, 6, 15},
    {15, 5, 1, 2, 7, 13},
    {13, 6, 2, 3, 9, 12},
    {9, 11, 4, 5, 16, 10},
    {12, 7, 3, 11, 8, 10},
    {12, 9, 8, 16, 14},
    {9, 3, 0, 4, 8},
    {13, 7, 9, 10, 14},
    {15, 6, 7, 12, 14},
    {15, 13, 12, 10, 16},
    {16, 5, 6, 13, 14},
    {14, 10, 8, 5, 15}
};
static const vector<coord3d> C30_seed_points = {
    {-5.84787895500947119e-01, -4.00918129916895627e-01, -7.05184461461710876e-01},
    {-1.25581455441115430e-01, -8.05923823815222096e-01, -3.14079386022112106e-01},
    {1.32211672555218540e-01, -3.29799539673337860e-01, -7.98549938197729126e-01},
    {-2.92175720934314598e-01, 2.59671810981700280e-01, -7.81750720323426074e-01},
    {-7.09293764098534729e-01, -5.10713463629080766e-01, 2.13909964515144029e-03},
    {-6.92430546345066400e-02, -5.07452121425790992e-01, 3.45923624716787614e-01},
    {4.50835906935334596e-01, -3.95641167383027959e-01, -1.48929415089563449e-01},
    {3.47874988529713769e-01, 2.62933153184990498e-01, -4.37966195251790036e-01},
    {-4.93630448124039667e-01, 8.20192292292472175e-02, 3.62722842591090389e-01},
    {-2.35837320127705918e-01, 5.58143513371131550e-01, -1.21747709584526562e-01},
    {-8.94170733034767712e-02, 6.43423563408358934e-01, 5.84759028601800401e-01},
    {-8.12254682504155667e-01, 1.47860856938937107e-01, -2.86897680517075149e-01},
    {4.30661888266364479e-01, 7.55234517451122134e-01, 8.99059887954494347e-02},
    {8.55049281755898005e-01, 1.65763166796083800e-01, 7.31067709211465494e-02},
    {5.84787924532465420e-01, 4.00919172707515981e-01, 7.05185720414510220e-01},
    {5.97256153759563868e-01, -3.10361117345800575e-01, 5.57577323096763555e-01},
    {1.35438451021443352e-02, -1.51507571596593664e-02, 8.73795808764026960e-01}
};

// Get precomputed seed data by type
static const Graph& seedNeighbours(buckinverse::SeedType s) {
    switch (s) {
        case buckinverse::SeedType::C20: return C20_seed_neighbours;
        case buckinverse::SeedType::C28: return C28_seed_neighbours;
        case buckinverse::SeedType::C30: return C30_seed_neighbours;
        default: assert(false && "Unknown seed type"); return C20_seed_neighbours;
    }
}

static const vector<coord3d>& seedPoints(buckinverse::SeedType s) {
    switch (s) {
        case buckinverse::SeedType::C20: return C20_seed_points;
        case buckinverse::SeedType::C28: return C28_seed_points;
        case buckinverse::SeedType::C30: return C30_seed_points;
        default: assert(false && "Unknown seed type"); return C20_seed_points;
    }
}

// Match seed vertices via canonical spiral.
// Computes canonical spirals for both the precomputed seed graph and the
// ep seed_state graph. Since both are isomorphic, they produce the same
// canonical spiral, giving permutations that compose into the mapping.
// Returns mapping: precomputed_vertex -> ep_vertex_id.
//TODO: Is this still needed?
static vector<int> matchSeedViaSpiralImpl(
    const Graph& precomp,
    const ExtensionPath& ep)
{
    int seed_N = (int)precomp.size();

    // 1. Canonical spiral of precomputed seed graph
    Triangulation T_pre(precomp);
    vector<int> spiral_pre;
    jumplist_t jumps_pre;
    vector<vector<node_t>> perms_pre;
    T_pre.get_spiral(spiral_pre, jumps_pre, perms_pre);
    // perms_pre[0][canonical_i] = precomputed_vertex_id

    // 2. Compact the scattered ep vertex IDs to 0..seed_N-1
    vector<int> ep_ids;
    ep_ids.reserve(seed_N);
    for (const auto& sv : ep.seed_state)
        ep_ids.push_back(sv.id);
    sort(ep_ids.begin(), ep_ids.end());

    map<int, int> to_compact;
    for (int i = 0; i < seed_N; i++)
        to_compact[ep_ids[i]] = i;

    // Build compact adjacency
    Graph compact_adj(seed_N, 6);
    for (const auto& sv : ep.seed_state) {
        int ci = to_compact[sv.id];
        uint8_t m = sv.active;
        for (int p = 0; p < 6; p++) {
            if (m & (1 << p))
                compact_adj.push_back(ci, to_compact[sv.nbr[p]]);
        }
    }

    // 3. Canonical spiral of ep seed graph
    Triangulation T_ep(compact_adj);
    vector<int> spiral_ep;
    jumplist_t jumps_ep;
    vector<vector<node_t>> perms_ep;
    T_ep.get_spiral(spiral_ep, jumps_ep, perms_ep);
    // perms_ep[0][canonical_i] = compact_ep_id

    // 4. Compose: precomputed vertex perms_pre[0][i] maps to ep vertex ep_ids[perms_ep[0][i]]
    vector<int> mapping(seed_N);
    for (int i = 0; i < seed_N; i++)
        mapping[perms_pre[0][i]] = ep_ids[perms_ep[0][i]];

    return mapping;
}

// Fallback: BFS-based graph isomorphism between precomputed seed and ep seed state.
// Returns mapping: precomputed_vertex -> ep_vertex_id.
static vector<int> matchSeedViaBFS(
    const Graph& precomp,
    const ExtensionPath& ep)
{
    int seed_N = (int)precomp.size();

    // Build adjacency from seed_state for fast lookup
    map<int, vector<int>> ep_adj;
    for (const auto& sv : ep.seed_state) {
        vector<int> nbrs;
        uint8_t m = sv.active;
        for (int p = 0; p < 6; p++) {
            if (m & (1 << p))
                nbrs.push_back(sv.nbr[p]);
        }
        ep_adj[sv.id] = nbrs;
    }

    // Try matching each precomputed degree-5 vertex with each ep degree-5 vertex
    for (int start_p = 0; start_p < seed_N; start_p++) {
        if ((int)precomp[start_p].size() != 5) continue;

        for (const auto& sv : ep.seed_state) {
            if ((int)ep_adj[sv.id].size() != 5) continue;

            // Try all 5 rotations of the CW neighbor list
            for (int rot = 0; rot < 5; rot++) {
                vector<int> mapping(seed_N, -1);
                bool valid = true;

                queue<int> q;
                mapping[start_p] = sv.id;
                q.push(start_p);

                while (!q.empty() && valid) {
                    int p_v = q.front(); q.pop();
                    int ep_v = mapping[p_v];
                    const auto& p_nbrs = precomp[p_v];
                    const auto& e_nbrs = ep_adj[ep_v];

                    if (p_nbrs.size() != e_nbrs.size()) { valid = false; break; }

                    int offset = -1;
                    if (p_v == start_p) {
                        offset = rot;
                    } else {
                        for (int j = 0; j < (int)p_nbrs.size(); j++) {
                            if (mapping[p_nbrs[j]] >= 0) {
                                int target = mapping[p_nbrs[j]];
                                for (int k = 0; k < (int)e_nbrs.size(); k++) {
                                    if (e_nbrs[k] == target) {
                                        offset = (k - j + (int)e_nbrs.size()) % (int)e_nbrs.size();
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                        if (offset < 0) { valid = false; break; }
                    }

                    int deg = (int)p_nbrs.size();
                    for (int j = 0; j < deg; j++) {
                        int p_nbr = p_nbrs[j];
                        int e_nbr = e_nbrs[(j + offset) % deg];

                        if (mapping[p_nbr] < 0) {
                            mapping[p_nbr] = e_nbr;
                            q.push(p_nbr);
                        } else if (mapping[p_nbr] != e_nbr) {
                            valid = false;
                            break;
                        }
                    }
                }

                if (valid) {
                    bool complete = true;
                    for (int i = 0; i < seed_N; i++)
                        if (mapping[i] < 0) { complete = false; break; }
                    if (complete) return mapping;
                }
            }
        }
    }

    assert(false && "Failed to match seed vertices");
    return {};
}

// Load precomputed seed geometry into the points array, mapped to the
// extension path's vertex IDs.
static void computeSeedGeometry(const ExtensionPath& ep, vector<coord3d>& points) {
    const auto& precomp_nbrs = seedNeighbours(ep.seed);
    const auto& precomp_pts = seedPoints(ep.seed);

    // Match via BFS isomorphism (faster than spiral for these small graphs)
    vector<int> mapping = matchSeedViaBFS(precomp_nbrs, ep);

    // Normalize to unit sphere (precomputed coords are at physical scale ~1.5 A)
    double scale = 0;
    for (const auto& p : precomp_pts)
        scale = max(scale, p.norm());

    for (int i = 0; i < (int)precomp_pts.size(); i++) {
        points[mapping[i]] = precomp_pts[i] / scale;
    }
}

Deltahedron Deltahedron::fromExtensionPath(const ExtensionPath& ep) {
    int full_N = ep.full_N;
    vector<coord3d> points(full_N);

    // 1. Compute seed geometry
    computeSeedGeometry(ep, points);

    // 2. Create ReducibleDual and load seed state
    ReducibleDual rd(full_N);
    for (const auto& sv : ep.seed_state) {
        for (int i = 0; i < ReducibleDual::D_MAX; i++)
            rd.V[sv.id].nbr[i] = sv.nbr[i];
        rd.V[sv.id].active = sv.active;
        rd.n_alive++;
        if (rd.degree(sv.id) == 5) rd.deg5.insert(sv.id);
    }

    // 3. For each expansion step: compute coords, expand topology, enforce convexity
    for (const auto& step : ep.steps) {
        computeStripCoords(step, points);
        rd.expand(step);
        liftStripToSurface(step, rd, points);
    }

    // 4. Extract compact Graph with renumbered coordinates
    vector<int> remap(full_N, -1);
    int id = 0;
    for (int u = 0; u < full_N; u++)
        if (rd.alive(u)) remap[u] = id++;

    Graph adj(id, 6);
    vector<coord3d> compact_points(id);
    for (int u = 0; u < full_N; u++) {
        if (!rd.alive(u)) continue;
        compact_points[remap[u]] = points[u];
        uint8_t m = rd.V[u].active;
        for (; m; m &= m - 1)
            adj.push_back(remap[u], remap[rd.V[u].nbr[__builtin_ctz(m)]]);
    }

    return Deltahedron(Triangulation(adj), compact_points);
}

// Helper: extract compact Deltahedron from ReducibleDual state + full points array.
// Also fills remap[full_N] with mapping from full IDs to compact 0..n-1 IDs (-1 if dead).
// Returns the compact Deltahedron and the number of alive vertices.
static Deltahedron extractCompact(
    const ReducibleDual& rd, int full_N,
    const vector<coord3d>& points,
    vector<int>& remap)
{
    remap.assign(full_N, -1);
    int id = 0;
    for (int u = 0; u < full_N; u++)
        if (rd.alive(u)) remap[u] = id++;

    Graph adj(id, 6);
    vector<coord3d> compact_pts(id);
    for (int u = 0; u < full_N; u++) {
        if (!rd.alive(u)) continue;
        compact_pts[remap[u]] = points[u];
        uint8_t m = rd.V[u].active;
        for (; m; m &= m - 1)
            adj.push_back(remap[u], remap[rd.V[u].nbr[__builtin_ctz(m)]]);
    }

    return Deltahedron(Triangulation(adj), compact_pts);
}

// Helper: extract a small patch sub-Deltahedron containing only the vertices
// near a just-expanded strip. Free vertices = strip + path + tp; boundary
// (fixed) = one-ring of (path ∪ tp) that is not in {strip, path, tp}.
// Returns the sub-Deltahedron plus:
//   free_mask[m]: true for free vertices in the sub-graph
//   interior_mask[m]: true for vertices whose full neighbor set is in the patch
//                     (i.e. patch degree == graph degree; false for truncated boundary)
//   remap[full_N]: full→sub ID mapping (-1 if not in patch)
static Deltahedron extractPatch(
    const ReducibleDual& rd, int full_N,
    const ExtensionStep& step,
    const vector<coord3d>& points,
    vector<bool>& free_mask,
    vector<bool>& interior_mask,
    vector<int>& remap,
    bool wider = false)
{
    // 1. Collect core free vertices (strip + path + tp)
    set<int> core_free;
    for (int v : step.strip) core_free.insert(v);
    for (int v : step.path)  core_free.insert(v);
    for (int v : step.tp)    core_free.insert(v);

    // 2. Expand to one-ring of (path ∪ tp) for boundary ring 1
    set<int> ring1;
    for (int v : step.path) {
        uint8_t m = rd.V[v].active;
        for (; m; m &= m - 1) {
            int nb = rd.V[v].nbr[__builtin_ctz(m)];
            if (!core_free.count(nb)) ring1.insert(nb);
        }
    }
    for (int v : step.tp) {
        uint8_t m = rd.V[v].active;
        for (; m; m &= m - 1) {
            int nb = rd.V[v].nbr[__builtin_ctz(m)];
            if (!core_free.count(nb)) ring1.insert(nb);
        }
    }

    // Free set: core + (if wider, ring1 is also free)
    set<int> free_set = core_free;
    if (wider) free_set.insert(ring1.begin(), ring1.end());

    // Patch set: free_set + boundary ring
    set<int> patch_set = free_set;
    if (wider) {
        // Ring 2: one-ring of ring1 not already in patch
        for (int v : ring1) {
            uint8_t m = rd.V[v].active;
            for (; m; m &= m - 1)
                patch_set.insert(rd.V[v].nbr[__builtin_ctz(m)]);
        }
    } else {
        patch_set.insert(ring1.begin(), ring1.end());
    }

    // 3. Remap patch vertices to 0..m-1
    remap.assign(full_N, -1);
    int id = 0;
    for (int v : patch_set) remap[v] = id++;
    int m = id;

    // 4. Build adjacency for the sub-graph (only edges within patch)
    Graph adj(m, 6);
    for (int u : patch_set) {
        uint8_t mask = rd.V[u].active;
        for (; mask; mask &= mask - 1) {
            int nb = rd.V[u].nbr[__builtin_ctz(mask)];
            if (remap[nb] >= 0)
                adj.push_back(remap[u], remap[nb]);
        }
    }

    // 5. Build free_mask
    free_mask.assign(m, false);
    for (int v : free_set)
        if (remap[v] >= 0) free_mask[remap[v]] = true;

    // 6. Build interior_mask: true iff all graph neighbors are in the patch
    //    (patch degree == graph degree). False for boundary vertices whose
    //    neighbor sets are truncated — E_conv would compute bogus h values.
    interior_mask.assign(m, false);
    for (int u : patch_set) {
        int graph_deg = rd.degree(u);
        int patch_deg = adj.degree(remap[u]);
        interior_mask[remap[u]] = (patch_deg == graph_deg);
    }

    // 7. Extract coordinates
    vector<coord3d> patch_pts(m);
    for (int v : patch_set) patch_pts[remap[v]] = points[v];

    return Deltahedron(Triangulation(adj), patch_pts);
}

Deltahedron Deltahedron::fromExtensionPathOptimized(const ExtensionPath& ep, FILE* log, StepCallback diag,
                                                     OptMethod method, double step_tol, double final_tol, long long max_work_per_step,
                                                     double step_angle_tol, double final_angle_tol,
                                                     OptMethod final_method, double patch_grad_tol,
                                                     bool global_post_patch_reflect) {
    int full_N = ep.full_N;
    vector<coord3d> points(full_N);
    PipelineDiag pd;  // Pipeline diagnostics accumulator

    // Timing accumulators (only used when opt_log is set)
    using clk = chrono::steady_clock;
    double ms_seed = 0, ms_place = 0, ms_reflect = 0, ms_patch = 0, ms_relax = 0, ms_final = 0;
    int total_relax_iters = 0, total_patch_iters = 0;
    int acc_energy = 0, acc_grad = 0, acc_hv = 0;
    auto t0 = clk::now();

    // 1. Compute seed geometry
    computeSeedGeometry(ep, points);
    pd.set_seed_type((int)ep.seed);  // 0=C20, 1=C28, 2=C30

    // 2. Create ReducibleDual and load seed state
    ReducibleDual rd(full_N);
    for (const auto& sv : ep.seed_state) {
        for (int i = 0; i < ReducibleDual::D_MAX; i++)
            rd.V[sv.id].nbr[i] = sv.nbr[i];
        rd.V[sv.id].active = sv.active;
        rd.n_alive++;
        if (rd.degree(sv.id) == 5) rd.deg5.insert(sv.id);
    }
    ms_seed = chrono::duration<double,milli>(clk::now() - t0).count();

    // Diagnostic: seed geometry snapshot
    if (diag) {
        vector<int> seed_remap;
        Deltahedron D_seed = extractCompact(rd, full_N, points, seed_remap);
        diag(0, "seed", D_seed);
    }

    // 3. For each expansion step: place strip, expand topology, optimize
    int step_idx = 0;
    for (const auto& step : ep.steps) {
        auto ts = clk::now();

        // a. Place strip vertices
        if (step.kind.type == ExpKind::F_type)
            shiftForFRing(rd, step, points);  // shift tp-side to make room
        computeStripCoords(step, points);

        // b. Expand topology + enforce outward convexity
        rd.expand(step);
        if (step.kind.type != ExpKind::F_type)
            liftStripToSurface(step, rd, points);

        auto t_place = clk::now();
        ms_place += chrono::duration<double,milli>(t_place - ts).count();

        // F-ring placement is exact — skip reflect, patch, and relaxation.
        if (step.kind.type == ExpKind::F_type) {
            pd.flags |= PipelineDiag::HAS_F_RING;
            if (diag) {
                vector<int> diag_remap;
                Deltahedron D_diag = extractCompact(rd, full_N, points, diag_remap);
                diag(step_idx + 1, "placed", D_diag);
                diag(step_idx + 1, "reflected", D_diag);
                diag(step_idx + 1, "patched", D_diag);
                diag(step_idx + 1, "relaxed", D_diag);
            }
            if (log)
                fprintf(log, "  step %2d: N=%3d F-ring (exact placement, no optimization)\n",
                        step_idx, (int)(rd.N()));
            step_idx++;
            continue;
        }

        // Diagnostic: after strip placement + lift
        if (diag) {
            vector<int> diag_remap;
            Deltahedron D_diag = extractCompact(rd, full_N, points, diag_remap);
            diag(step_idx + 1, "placed", D_diag);
        }

        // c. Reflect all concave vertices on the full graph (before patch extraction).
        //    This is O(N) and saves the optimizer from spending expensive iterations
        //    fighting concavity from inverted strip vertices or pre-existing concavities.
        {
            vector<int> refl_remap;
            Deltahedron D = extractCompact(rd, full_N, points, refl_remap);

            D.reflect_all_concave(D.points);

            // Copy reflected coords back to full array
            for (int u = 0; u < full_N; u++)
                if (refl_remap[u] >= 0)
                    points[u] = D.points[refl_remap[u]];
        }
        auto t_refl = clk::now();
        ms_reflect += chrono::duration<double,milli>(t_refl - t_place).count();

        // Diagnostic: after reflect
        if (diag) {
            vector<int> diag_remap;
            Deltahedron D_diag = extractCompact(rd, full_N, points, diag_remap);
            diag(step_idx + 1, "reflected", D_diag);
        }

        // d-f. Patch reflect-optimize loop: optimize patch without hard convexity
        //       constraint (which causes it to get stuck at h=0 boundary), then
        //       reflect any concavities on the full graph and re-optimize.
        //       Cycle guard: if reflection count doesn't decrease, break early
        //       (the full-graph optimizer will handle remaining concavities).
        {
            int prev_refl = INT_MAX;
            int patch_round;
            for (patch_round = 0; patch_round < 10; patch_round++) {
                // d. Extract small patch sub-graph (O(1) vertices)
                vector<bool> free_mask, interior_mask;
                vector<int> patch_remap;
                Deltahedron patch = extractPatch(rd, full_N, step, points,
                                                 free_mask, interior_mask, patch_remap,
                                                 global_post_patch_reflect);

                // e. Trust-region optimize on patch only (no hard convexity constraint;
                //    E_conv softplus still provides a soft bias toward convexity)
                patch.opt_log = log;
                patch.optimize_patch(patch.points, free_mask, interior_mask, 0, 150, patch_grad_tol, false);
                patch.opt_log = nullptr;
                total_patch_iters += patch.iterations_used;

                // f. Copy patch free-vertex coords back to full array
                for (int u = 0; u < full_N; u++)
                    if (patch_remap[u] >= 0 && free_mask[patch_remap[u]])
                        points[u] = patch.points[patch_remap[u]];

                // g. Reflect concave on full graph; break if convex or cycling
                vector<int> refl_remap;
                Deltahedron D = extractCompact(rd, full_N, points, refl_remap);
                int n_refl = D.reflect_all_concave(D.points);
                for (int u = 0; u < full_N; u++)
                    if (refl_remap[u] >= 0)
                        points[u] = D.points[refl_remap[u]];

                if (n_refl == 0) break;  // patch produced convex result

                if (log)
                    fprintf(log, "    step %2d patch_round %d: reflected=%d, conc=%d, ang=%.2e\n",
                            step_idx, patch_round, n_refl, D.count_concave(), D.max_angle_relerr());

                // Cycle guard: if not making progress, stop looping
                if (n_refl >= prev_refl) {
                    pd.flags |= PipelineDiag::PATCH_CYCLING;
                    break;
                }
                prev_refl = n_refl;
            }
            pd.set_max_patch_rounds(patch_round + 1);
        }

        auto t_patch = clk::now();
        ms_patch += chrono::duration<double,milli>(t_patch - t_refl).count();

        // Diagnostic: after patch optimize, before full-graph relaxation
        if (diag) {
            vector<int> diag_remap;
            Deltahedron D_diag = extractCompact(rd, full_N, points, diag_remap);
            diag(step_idx + 1, "patched", D_diag);
        }

        // g. Full-graph reflect-optimize loop.
        //    Reflect concave vertices, then optimize.  If cycling is detected
        //    (concavity not decreasing), escalate to convex hull projection.
        vector<int> full_remap;
        Deltahedron D = extractCompact(rd, full_N, points, full_remap);

        D.opt_k_flat = 0;  // skip E_flat in intermediate steps
        D.opt_k_conv = 0;   // pure quality; reflect/hull handles convexity
        D.opt_skip_post_reflect = true;  // we handle convexity in the loop
        D.opt_method = method;
        int n_steps = (int)ep.steps.size();
        bool log_this = log && (step_idx <= 2 || step_idx == n_steps/2 || step_idx >= n_steps - 2);

        {
          int prev_conc = INT_MAX;
          bool used_hull = false;
          int round;
          for(round = 0; round < 10; round++){
            int n_fix;
            if(!used_hull){
              n_fix = D.reflect_all_concave(D.points);
            } else {
              n_fix = D.project_onto_convex_hull(D.points);
            }
            if(round > 0 && n_fix == 0) break;  // stable in convex basin

            // Optimize
            if (log_this) D.opt_log = log;
            OptResult step_result = D.optimize(D.points, 0, step_tol, {}, max_work_per_step, step_angle_tol);
            if (log_this) D.opt_log = nullptr;
            pd.flags |= D.opt_diag_flags;  // accumulate optimizer-level flags
            if(step_result == OptResult::STAGNATED) {
              pd.flags |= PipelineDiag::STAG_STEP;
              pd.inc_stag_steps();
            }
            total_relax_iters += D.iterations_used;
            acc_energy += D.n_energy_evals;
            acc_grad += D.n_grad_evals;
            acc_hv += D.n_hv_evals;

            int conc = D.count_concave();
            if(log && (log_this || conc > 0))
              fprintf(log, "    step %2d round %d: %s=%d, %d iters, ang=%.2e, conc=%d\n",
                      step_idx, round, used_hull ? "projected" : "reflected",
                      n_fix, D.iterations_used, D.max_angle_relerr(), conc);

            if(conc == 0) break;  // success

            // Cycle guard: if concavity not decreasing, escalate to hull projection
            if(conc >= prev_conc){
              if(!used_hull){
                pd.flags |= PipelineDiag::REFLECT_CYCLING_STEP;
                used_hull = true;  // escalate — retry with hull projection
                pd.flags |= PipelineDiag::HULL_USED_STEP;
                if(log) fprintf(log, "    step %2d: reflect cycling (%d→%d), escalating to hull projection\n",
                                step_idx, prev_conc, conc);
              } else {
                pd.flags |= PipelineDiag::HULL_CYCLING_STEP;
                break;  // hull projection also cycling — give up
              }
            }
            prev_conc = conc;
          }
          pd.set_max_reflect_rounds_step(round + 1);
        }

        // Convexity invariant: per-step must produce convex geometry
        if(D.count_concave() > 0) {
          pd.flags |= PipelineDiag::CONVEXITY_FAIL_STEP;
          pd.inc_cvx_fail_steps();
          fprintf(stderr, "CONVEXITY FAILURE: step %d (N=%d) ended with %d concave vertices\n",
                  step_idx, D.N, D.count_concave());
        }

        // h. Write optimized coordinates back to full array
        for (int u = 0; u < full_N; u++)
            if (full_remap[u] >= 0)
                points[u] = D.points[full_remap[u]];

        auto t_relax = clk::now();
        ms_relax += chrono::duration<double,milli>(t_relax - t_patch).count();

        // Diagnostic: after full-graph relaxation
        if (diag) {
            vector<int> diag_remap;
            Deltahedron D_diag = extractCompact(rd, full_N, points, diag_remap);
            D_diag.iterations_used = D.iterations_used;
            diag(step_idx + 1, "relaxed", D_diag);
        }

        if (log) {
            fprintf(log, "  step %2d: N=%3d place=%.1f refl=%.1f patch=%.1f(%d) relax=%.1f(%d) ms\n",
                    step_idx, D.N,
                    chrono::duration<double,milli>(t_place - ts).count(),
                    chrono::duration<double,milli>(t_refl - t_place).count(),
                    chrono::duration<double,milli>(t_patch - t_refl).count(), total_patch_iters,
                    chrono::duration<double,milli>(t_relax - t_patch).count(), D.iterations_used);
        }
        step_idx++;
    }

    // 4. Final extraction + full optimization
    auto t_final = clk::now();
    vector<int> remap;
    Deltahedron D = extractCompact(rd, full_N, points, remap);
    // Diagnostic: before final optimization
    if (diag) diag((int)ep.steps.size() + 1, "patched", D);

    D.opt_k_flat = 0;  // geometry is already 3D from seed expansion; E_flat fights equilateral

    // Final reflect-optimize loop with hull projection escalation.
    // Reflect concave vertices, then optimize.  If cycling is detected
    // (concavity not decreasing), escalate to convex hull projection.
    D.opt_k_conv = 0;
    D.opt_convex_constraint = false;
    D.opt_skip_post_reflect = true;  // we handle convexity in the loop
    D.opt_method = final_method;
    int total_final_iters = 0;

    {
      int prev_conc = INT_MAX;
      bool used_hull = false;
      int round;
      for(round = 0; round < 10; round++){
        int n_fix;
        if(!used_hull){
          n_fix = D.reflect_all_concave(D.points);
        } else {
          n_fix = D.project_onto_convex_hull(D.points);
        }
        if(round > 0 && n_fix == 0) break;  // stable in convex basin

        if(diag && round == 0) diag((int)ep.steps.size() + 1, "reflected", D);

        // Optimize
        if (log) D.opt_log = log;
        OptResult fr = D.optimize(D.points, 0, final_tol, {}, max_work_per_step, final_angle_tol);
        D.opt_log = nullptr;
        pd.flags |= D.opt_diag_flags;  // accumulate optimizer-level flags
        if(fr == OptResult::STAGNATED) pd.flags |= PipelineDiag::STAG_FINAL;
        total_final_iters += D.iterations_used;
        acc_energy += D.n_energy_evals;
        acc_grad += D.n_grad_evals;
        acc_hv += D.n_hv_evals;

        int conc = D.count_concave();
        if(log) fprintf(log, "  final round %d: %s=%d, %d iters, ang=%.2e, conc=%d\n",
                        round, used_hull ? "projected" : "reflected",
                        n_fix, D.iterations_used, D.max_angle_relerr(), conc);

        if(conc == 0) break;  // success

        // Cycle guard: if concavity not decreasing, escalate to hull projection
        if(conc >= prev_conc){
          if(!used_hull){
            pd.flags |= PipelineDiag::REFLECT_CYCLING_FINAL;
            used_hull = true;  // escalate — retry with hull projection
            pd.flags |= PipelineDiag::HULL_USED_FINAL;
            if(log) fprintf(log, "  final: reflect cycling (%d→%d), escalating to hull projection\n",
                            prev_conc, conc);
          } else {
            pd.flags |= PipelineDiag::HULL_CYCLING_FINAL;
            break;  // hull projection also cycling — give up
          }
        }
        prev_conc = conc;
      }
      pd.set_final_reflect_rounds(round + 1);
    }

    if(D.count_concave() > 0 && log) {
      fprintf(log, "  final loop ended with %d concave, attempting hull+reflect cleanup\n",
              D.count_concave());
    }

    // Final constrained Steihaug polish: h>=0 trust region prevents regression
    // to concave.  k_conv=0 so the Hessian is pure quality — the constraint
    // is the only convexity mechanism.
    // Use hull projection + reflect to ensure all vertices start convex.
    // Hull projection handles deep concavities (large graphs); reflect handles
    // remaining h < 0 cases (small graphs where all vertices are on the hull).
    D.project_onto_convex_hull(D.points);
    D.reflect_all_concave(D.points);
    D.opt_convex_constraint = true;
    D.opt_method = OptMethod::STEIHAUG;
    if (log) D.opt_log = log;
    OptResult final_result = D.optimize(D.points, 0, final_tol, {}, max_work_per_step, final_angle_tol);
    D.opt_log = nullptr;
    pd.flags |= D.opt_diag_flags;  // accumulate optimizer-level flags from constrained phase
    if(final_result == OptResult::STAGNATED) pd.flags |= PipelineDiag::STAG_CONSTRAINED;
    if(final_result == OptResult::BUDGET_EXHAUSTED) pd.flags |= PipelineDiag::BUDGET_CONSTRAINED;
    total_final_iters += D.iterations_used;
    acc_energy += D.n_energy_evals;
    acc_grad += D.n_grad_evals;
    acc_hv += D.n_hv_evals;
    D.opt_convex_constraint = false;

    if(log) fprintf(log, "  final constrained: %d iters, ang=%.2e, conc=%d, %s\n",
                    D.iterations_used, D.max_angle_relerr(), D.count_concave(),
                    opt_result_name(final_result));

    // Override result if convexity loop couldn't fix concavity
    if(D.count_concave() > 0)
      final_result = OptResult::CONVEXITY_STUCK;

    D.iterations_used = total_final_iters;
    D.final_opt_result = final_result;
    pd.set_final_result(final_result);
    D.diag = pd;
    ms_final = chrono::duration<double,milli>(clk::now() - t_final).count();

    // Diagnostic: after final optimization
    if (diag) diag((int)ep.steps.size() + 1, "final", D);

    // Set accumulated eval counters
    D.n_energy_evals = acc_energy;
    D.n_grad_evals = acc_grad;
    D.n_hv_evals = acc_hv;

    if (log) {
        double ms_total = chrono::duration<double,milli>(clk::now() - t0).count();
        fprintf(log, "  totals: seed=%.0f place=%.0f refl=%.0f patch=%.0f(%d) relax=%.0f(%d) final=%.0f(%d) total=%.0f ms result=%s\n",
                ms_seed, ms_place, ms_reflect, ms_patch, total_patch_iters,
                ms_relax, total_relax_iters, ms_final, D.iterations_used, ms_total,
                opt_result_name(D.final_opt_result));
    }

    D.iterations_used += total_relax_iters;
    return D;
}

// =====================================================================
// Deltahedron geometry optimization toward equilateral triangles
// =====================================================================

// Smallest eigenpair of a 3x3 symmetric PSD matrix.
// Returns {lambda_min, eigenvector}, or {0, {}} if degenerate.
static pair<double,coord3d> smallest_eigenpair_3x3(const matrix3d& A)
{
  coord3d lambda = A.eigenvalues();

  // Find smallest eigenvalue (clamp negative FP noise for PSD)
  int imin = 0;
  if(lambda[1] < lambda[0]) imin = 1;
  if(lambda[2] < lambda[imin]) imin = 2;
  double lmin = max(0.0, lambda[imin]);
  if(lmin == 0) return {0.0, coord3d()};

  // Compute eigenvector as cross product of two rows of (A - lmin*I).
  // For a rank-2 matrix, two of the three rows are linearly independent,
  // and their cross product gives the null-space direction.
  matrix3d B; // B = A - lmin*I
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      B(i,j) = A(i,j) - (i==j ? lmin : 0);

  coord3d row0(B(0,0), B(0,1), B(0,2));
  coord3d row1(B(1,0), B(1,1), B(1,2));
  coord3d row2(B(2,0), B(2,1), B(2,2));

  // Try all three row-pair cross products, pick the largest
  coord3d c01 = row0.cross(row1);
  coord3d c02 = row0.cross(row2);
  coord3d c12 = row1.cross(row2);
  double n01 = c01.norm(), n02 = c02.norm(), n12 = c12.norm();

  coord3d n;
  if(n01 >= n02 && n01 >= n12) n = c01 / n01;
  else if(n02 >= n12)          n = c02 / n02;
  else                         n = c12 / n12;

  if(n.norm() < 0.5) return {0.0, coord3d()};  // degenerate
  return {lmin, n};
}

// Compute total energy and gradient for the equilateral-triangle force field.
// Four terms: E_bond (edge lengths), E_angle (triangle angles), E_curv (Gaussian curvature),
// E_flat (face centroid ring flatness for deg<=6 vertices with all deg<=6 neighbors).
// Returns total energy. If grad is non-null, gradient is accumulated (zeroed at start).
static double deltahedron_energy_and_gradient(
    const DeltahedronView<double>& D,
    const vector<edge_t>& edges,
    std::span<const coord3d> x,
    vector<coord3d>* grad,
    double L,           // target edge length
    double k_bond,
    double k_angle,
    double k_curv,
    double k_flat,
    double k_conv, double sigma_conv = 0,  // sigma_conv>0: softplus E_conv, 0: quadratic
    const vector<bool>& conv_mask = {})  // if non-empty, restrict E_conv to these vertices
{
  const int N = D.N;
  double energy = 0.0;

  // Zero gradient
  if(grad)
    for(int i = 0; i < N; i++) (*grad)[i] = coord3d(0,0,0);

  // --- E_bond = (k/2) Sum_edges (r - L)^2 ---
  for(const edge_t& e : edges){
    auto ed = EdgeBondData::compute(x, e.first, e.second, L);
    if(!ed.valid) continue;
    energy += ed.energy(k_bond);
    if(grad) ed.scatter_gradient(k_bond, *grad);
  }

  // --- E_angle = (k/2) Sum_corners (theta - pi/3)^2 ---
  for(const auto& tri : D.triangles()){
    for(int c = 0; c < 3; c++){
      auto ca = CornerAngleData::compute(x, tri[(c+2)%3], tri[c], tri[(c+1)%3]);
      if(!ca.valid) continue;
      energy += 0.5 * k_angle * ca.dev * ca.dev;
      if(grad) ca.scatter_gradient(k_angle * ca.dev, *grad);
    }
  }

  // --- E_curv = (k/2) Sum_v (K_target - angle_sum)^2 ---
  // K(v) = 2pi - angle_sum(v), K_target = 2pi - deg*pi/3.
  // dev = K - K_target = deg*pi/3 - angle_sum.
  for(int v = 0; v < N; v++){
    int d = D.degree(v);

    // Compute per-fan-angle data
    double angle_sum = 0;
    vector<CornerAngleData> fan(d);
    for(int i = 0; i < d; i++){
      fan[i] = CornerAngleData::compute(x, D[v][i], v, D[v][(i+1)%d]);
      if(fan[i].valid) angle_sum += fan[i].theta;
    }

    double dev = d * M_PI / 3.0 - angle_sum;
    energy += 0.5 * k_curv * dev * dev;

    if(grad){
      double w = -k_curv * dev;  // d(dev)/d(angle) = -1
      for(int i = 0; i < d; i++)
        if(fan[i].valid) fan[i].scatter_gradient(w, *grad);
    }
  }

  // E_flat: face centroid ring flatness (from flatness.tex).
  // For each qualifying vertex v, the centroids of surrounding triangles should
  // be approximately coplanar. The smallest eigenvalue lambda_0 of X^T X
  // (centroid-corrected) measures deviation from planarity.
  // E_flat = (k_flat/2) * Sum_v lambda_0(v)
  if(k_flat > 0){
    for(int v = 0; v < N; v++){
      int d = D.degree(v);
      if(d > 6) continue;

      // 1. Compute face centroids and their mean
      vector<coord3d> fc(d);
      coord3d c_bar(0,0,0);
      for(int i = 0; i < d; i++){
        int ni  = D[v][i];
        int ni1 = D[v][(i+1) % d];
        fc[i] = (x[v] + x[ni] + x[ni1]) / 3.0;
        c_bar += fc[i];
      }
      c_bar /= (double)d;

      // 2. Centroid-corrected coordinates
      for(int i = 0; i < d; i++) fc[i] -= c_bar;

      // 3. Build A = X^T X (3x3 symmetric)
      matrix3d A;
      for(int a = 0; a < 3; a++)
        for(int b = a; b < 3; b++){
          double s = 0;
          for(int i = 0; i < d; i++) s += fc[i][a] * fc[i][b];
          A(a,b) = s;
          A(b,a) = s;
        }

      // 4. Smallest eigenvalue and eigenvector (robust solver)
      auto [lambda0, n] = smallest_eigenpair_3x3(A);
      // Skip when already flat: relative threshold
      double trA = A(0,0) + A(1,1) + A(2,2);
      if(lambda0 < 1e-12 * trA) continue;

      // Scale-invariant measure: lambda0 / trA (dimensionless, in [0, 1/3])
      // E_flat = (k_flat/2) * sum_v (lambda0 / trA)
      // Full gradient via quotient rule:
      //   f_i = (k/trA) * [(X_i.n)n - (lambda0/trA) X_i]
      //   dE/dx[n_j] = (1/3)(f_j + f_{j-1})
      //   dE/dx[v] = 0 (v appears in all centroids equally, cancels via centroid correction)
      double scale = k_flat / trA;
      double ratio = lambda0 / trA;
      energy += 0.5 * scale * lambda0;

      if(grad){
        for(int j = 0; j < d; j++){
          int jprev = (j + d - 1) % d;
          coord3d fj    = n * fc[j].dot(n)    - fc[j]    * ratio;
          coord3d fjpre = n * fc[jprev].dot(n) - fc[jprev] * ratio;
          (*grad)[D[v][j]] += (fj + fjpre) * (scale / 3.0);
        }
      }
    }
  }

  // E_conv: convexity penalty.
  // sigma_conv > 0: softplus  E = k * sigma * log(1 + exp(-h/sigma)).
  //                 Smooth (C-inf), linear for h << -sigma, zero for h >> sigma.
  // sigma_conv == 0: quadratic  E = k * h^2 for h < 0.
  //                 Simpler, C1 at h=0 (discontinuous Hessian).
  if(k_conv > 0){
    for(int v = 0; v < N; v++){
      int d = D.degree(v);
      if(d > 6) continue;
      if(!conv_mask.empty() && !conv_mask[v]) continue;

      auto vd = grad ? VertexHData::compute_derivs(D, x, v)
                     : VertexHData::compute_h(D, x, v);
      if(!vd.valid) continue;

      double dEdh;
      if(sigma_conv > 0){
        // Softplus: E = k * sigma * log(1 + exp(-h/sigma))
        double z = vd.h / sigma_conv;
        double sig = (z > 20) ? 0.0 : (z < -20) ? 1.0 : 1.0 / (1.0 + exp(z));
        if(sig < 1e-15) continue;  // h >> sigma: numerically zero

        // log(1+exp(-z)) overflows for z << -20; use -z asymptote there
        if(z < -20)
          energy += k_conv * (-vd.h);  // k * sigma * (-z) = -k * h
        else
          energy += k_conv * sigma_conv * log1p(exp(-z));

        dEdh = -k_conv * sig;
      } else {
        // Quadratic: E = k * h^2 for h < 0
        if(vd.h >= 0) continue;
        energy += k_conv * vd.h * vd.h;
        dEdh = 2.0 * k_conv * vd.h;
      }

      if(grad)
        vd.scatter_dhdx(dEdh, *grad);
    }
  }

  return energy;
}

// Dot product of two vector<coord3d> arrays (sum of component-wise products).
static double vec_dot(const vector<coord3d>& a, const vector<coord3d>& b){
  double s = 0;
  for(int i = 0; i < (int)a.size(); i++)
    s += a[i].dot(b[i]);
  return s;
}

static double vec_norm(const vector<coord3d>& a){
  return sqrt(vec_dot(a, a));
}

// a[i] += alpha * b[i]
static void vec_axpy(vector<coord3d>& a, double alpha, const vector<coord3d>& b){
  for(int i = 0; i < (int)a.size(); i++)
    a[i] = a[i] + b[i] * alpha;
}

static void vec_scale(vector<coord3d>& a, double alpha){
  for(int i = 0; i < (int)a.size(); i++)
    a[i] = a[i] * alpha;
}

static void vec_zero(vector<coord3d>& a){
  for(int i = 0; i < (int)a.size(); i++)
    a[i] = coord3d(0,0,0);
}

// ---------- Hessian-vector product (matrix-free) ----------
//
// Computes Hv = H * v for the same energy terms used in optimize():
// E_bond, E_angle, E_curv, and E_flat.  Does NOT assemble H.
//
// Caller must zero-initialize Hv before calling.
//
static void deltahedron_hv_product(
    const DeltahedronView<double>& D,
    const vector<edge_t>& edges,
    std::span<const coord3d> x,
    const vector<coord3d>& v,
    vector<coord3d>& Hv,
    double L,
    double k_bond, double k_angle, double k_curv, double k_flat,
    const vector<bool>& fixed = {},
    double k_conv = 0, double sigma_conv = 0)
{
  const int N = D.N;
  const matrix3d I = matrix3d::unit_matrix();
  const double theta0 = M_PI / 3.0;

  // Helper: apply 3x3 block M to vertex pair (i,j), accumulating M*v[j] into Hv[i]
  auto apply_block = [&](int vi, int vj, const matrix3d& M){
    Hv[vi] = Hv[vi] + M * v[vj];
  };

  // --- E_bond Hv ---
  if(k_bond > 0){
    for(const edge_t& e : edges){
      auto ed = EdgeBondData::compute(x, e.first, e.second, L);
      if(!ed.valid) continue;
      ed.scatter_hv(k_bond, L, v, Hv);
    }
  }

  // --- E_angle Hv ---
  if(k_angle > 0){
    for(const auto& tri : D.triangles()){
      for(int c = 0; c < 3; c++){
        auto ca = CornerAngleData::compute(x, tri[(c+2)%3], tri[c], tri[(c+1)%3]);
        if(!ca.valid) continue;
        ca.scatter_hv(k_angle, v, Hv);
      }
    }
  }

  // --- E_curv Hv ---
  // H_curv = k_curv * sum_v [dK⊗dK + dev * d²K/dx²]
  // Rank-1: gather dK.v over fan, scatter dK.
  // Correction: sum of per-fan-angle Hessian blocks (same as E_angle).
  if(k_curv > 0){
    for(int vertex = 0; vertex < N; vertex++){
      int deg = D.degree(vertex);

      // Compute per-fan-angle data
      double angle_sum = 0;
      vector<CornerAngleData> fan(deg);
      for(int i = 0; i < deg; i++){
        fan[i] = CornerAngleData::compute(x, D[vertex][i], vertex,
                                          D[vertex][(i+1)%deg]);
        if(fan[i].valid) angle_sum += fan[i].theta;
      }
      double dev = deg * M_PI / 3.0 - angle_sum;

      // Rank-1: Hv += k_curv * (dK.v) * dK
      // dK = d(angle_sum)/dx, scatter via fan[i].scatter_gradient
      double dK_dot_v = 0;
      for(int i = 0; i < deg; i++)
        if(fan[i].valid){
          // dK.v contribution from fan angle i: p.v[a] + q.v[d] - (p+q).v[b]
          dK_dot_v += fan[i].p.dot(v[fan[i].a] - v[vertex])
                    + fan[i].q.dot(v[fan[i].d] - v[vertex]);
        }
      for(int i = 0; i < deg; i++)
        if(fan[i].valid) fan[i].scatter_gradient(k_curv * dK_dot_v, Hv);

      // Correction: Hv += (-k_curv * dev) * d²(angle_sum)/dx² * v
      if(fabs(dev) > 1e-15){
        double w2 = -k_curv * dev;
        for(int i = 0; i < deg; i++)
          if(fan[i].valid) fan[i].scatter_hv(w2, v, Hv);
      }
    }
  }

  // --- E_conv Hv ---
  // sigma_conv > 0: softplus E = k*sigma*log(1+exp(-h/sigma)).
  // sigma_conv == 0: quadratic E = k*h² for h<0, 0 for h>=0.
  // Exact Hv = d²E/dh² * (dh/dx ⊗ dh/dx) * v  +  dE/dh * d²h/dx² * v.
  if(k_conv > 0){
    for(int vertex = 0; vertex < N; vertex++){
      int d = D.degree(vertex);
      if(d > 6) continue;

      auto vd = VertexHData::compute_derivs(D, x, vertex);
      if(!vd.valid) continue;

      double dEdh, d2Edh2;
      if(sigma_conv > 0){
        double z = vd.h / sigma_conv;
        double sig = (z > 20) ? 0.0 : (z < -20) ? 1.0 : 1.0 / (1.0 + exp(z));
        if(sig < 1e-15) continue;
        dEdh   = -k_conv * sig;
        d2Edh2 = (k_conv / sigma_conv) * sig * (1.0 - sig);
      } else {
        if(vd.h >= 0) continue;
        dEdh   = 2.0 * k_conv * vd.h;
        d2Edh2 = 2.0 * k_conv;
      }

      // Rank-1: Hv += d²E/dh² * (dh/dx · v) * dh/dx
      vd.scatter_dhdx(d2Edh2 * vd.dhdx_dot(v), Hv);

      // Correction: Hv += dE/dh * d²h/dx² * v
      vd.scatter_d2h_hv(dEdh, v, Hv);
    }
  }

  // --- E_flat Hv --- (phase 1 only, via finite-difference)
  // Hv_flat = (g_flat(x + eps*v) - g_flat(x - eps*v)) / (2*eps)
  if(k_flat > 0){
    double eps = 1e-7 * L;
    vector<coord3d> x_plus(N), x_minus(N), g_plus(N), g_minus(N);
    for(int i = 0; i < N; i++){
      x_plus[i]  = x[i] + v[i] * eps;
      x_minus[i] = x[i] - v[i] * eps;
    }
    // Compute E_flat gradient at x+eps*v and x-eps*v
    // Use the full energy_and_gradient but only with k_flat, everything else zero
    deltahedron_energy_and_gradient(D, edges, x_plus,  &g_plus,  L, 0, 0, 0, k_flat, 0);
    deltahedron_energy_and_gradient(D, edges, x_minus, &g_minus, L, 0, 0, 0, k_flat, 0);

    double inv2eps = 1.0 / (2.0 * eps);
    for(int i = 0; i < N; i++)
      Hv[i] = Hv[i] + (g_plus[i] - g_minus[i]) * inv2eps;
  }

  // Zero out fixed vertices
  if(!fixed.empty())
    for(int i = 0; i < N; i++)
      if(fixed[i]) Hv[i] = coord3d(0,0,0);
}

// Compute energy only (no gradient), for line search.
static double deltahedron_energy_only(
    const DeltahedronView<double>& D,
    const vector<edge_t>& edges,
    std::span<const coord3d> x,
    double L, double k_bond, double k_angle, double k_curv, double k_flat, double k_conv,
    double sigma_conv = 0,
    const vector<bool>& conv_mask = {})
{
  return deltahedron_energy_and_gradient(D, edges, x, nullptr, L, k_bond, k_angle, k_curv, k_flat, k_conv, sigma_conv, conv_mask);
}

// Compute signed convexity height h for all vertices.
// h > 0 = convex, h < 0 = concave.  Fixed or high-degree vertices get h = 1.0.
static void compute_h_values(const DeltahedronView<double>& D, std::span<const coord3d> x,
                              vector<double>& h, const vector<bool>& fixed = {})
{
  int N = D.N;
  h.resize(N);
  for(int v = 0; v < N; v++){
    int d = D.degree(v);
    if(d > 6 || (!fixed.empty() && fixed[v])){
      h[v] = 1.0;
      continue;
    }
    auto vd = VertexHData::compute_h(D, x, v);
    h[v] = vd.h;  // 0 if degenerate
  }
}

// Check convexity constraint: h(v) >= -tau*L for all free/interior vertices.
// Returns true if all constraints satisfied.
static bool check_convexity(const DeltahedronView<double>& D, std::span<const coord3d> x,
                            const vector<bool>& free_mask, double L, double tau = 0.05,
                            const vector<bool>& interior_mask = {})
{
  for(int v = 0; v < D.N; v++){
    bool is_interior = !interior_mask.empty() && interior_mask[v];
    if(!free_mask[v] && !is_interior) continue;
    auto vd = VertexHData::compute_h(D, x, v);
    if(vd.valid && vd.h < -tau * L) return false;
  }
  return true;
}

// Assemble exact analytical Hessian for the patch optimizer.
// E_bond: exact.  E_angle: exact.  E_conv: exact (including d²h/dx² correction).
// H must be ndof x ndof (ndof = 3*nfree), zero-initialized by caller.
static void assemble_patch_hessian(
    const DeltahedronView<double>& D,
    const vector<edge_t>& edges,
    std::span<const coord3d> x,
    vector<vector<double>>& H,
    const vector<int>& free_idx,
    double L,
    double k_bond, double k_angle, double k_conv, double sigma,
    const vector<bool>& conv_mask)
{
  int nfree = (int)free_idx.size();
  const matrix3d I = matrix3d::unit_matrix();

  // Reverse map: vertex -> index in free_idx (-1 if not free)
  vector<int> fmap(D.N, -1);
  for(int k = 0; k < nfree; k++) fmap[free_idx[k]] = k;

  // Scatter 3x3 block M into H at vertex-pair (vi, vj).
  auto add_block = [&](int vi, int vj, const matrix3d& M){
    int ki = fmap[vi], kj = fmap[vj];
    if(ki < 0 || kj < 0) return;
    for(int p = 0; p < 3; p++)
      for(int q = 0; q < 3; q++)
        H[3*ki+p][3*kj+q] += M(p,q);
  };

  // --- E_bond = (k/2) Sum_edges (r - L)^2 ---
  if(k_bond > 0){
    for(const edge_t& e : edges){
      auto ed = EdgeBondData::compute(x, e.first, e.second, L);
      if(!ed.valid) continue;
      ed.scatter_blocks(k_bond, L, add_block);
    }
  }

  // --- E_angle = (k/2) Sum_corners (theta - pi/3)^2 ---
  if(k_angle > 0){
    for(const auto& tri : D.triangles()){
      for(int c = 0; c < 3; c++){
        auto ca = CornerAngleData::compute(x, tri[(c+2)%3], tri[c], tri[(c+1)%3]);
        if(!ca.valid) continue;
        ca.scatter_blocks(k_angle, add_block);
      }
    }
  }

  // --- E_conv = k * Sum_v softplus(-h/sigma) ---
  // Exact Hessian: d²E/dh² * (dh/dx)⊗(dh/dx) + dE/dh * d²h/dx²
  // where d²E/dh² = (k/sigma²) * sig * (1 - sig),  dE/dh = -(k/sigma) * sig.
  if(k_conv > 0){
    for(int v = 0; v < D.N; v++){
      int d = D.degree(v);
      if(d > 6) continue;
      if(!conv_mask.empty() && !conv_mask[v]) continue;

      auto vd = VertexHData::compute_derivs(D, x, v);
      if(!vd.valid) continue;

      double z = vd.h / sigma;
      double sig = (z > 20) ? 0.0 : (z < -20) ? 1.0 : 1.0 / (1 + exp(z));

      double d2Edh2 = (k_conv / sigma) * sig * (1 - sig);
      double dEdh   = -k_conv * sig;
      if(fabs(d2Edh2) < 1e-15 && fabs(dEdh) < 1e-15) continue;

      // Rank-1: d²E/dh² * (dh/dx)⊗(dh/dx)
      vd.scatter_rank1_blocks(d2Edh2, add_block);

      // Correction: dE/dh * d²h/dx²
      vd.scatter_d2h_blocks(dEdh, add_block);
    }
  }
}

bool Deltahedron::optimize_patch(std::span<const coord3d> initial_geometry,
                                 const vector<bool>& free_mask,
                                 const vector<bool>& interior_mask,
                                 double target_L, int max_iter, double grad_tol,
                                 bool convex_constraint)
{
  assert((int)initial_geometry.size() == N);
  assert((int)free_mask.size() == N);
  std::copy(initial_geometry.begin(), initial_geometry.end(), points.begin());

  vector<edge_t> edges = undirected_edges();

  // Target edge length from boundary-only edges (both endpoints fixed).
  double L = target_L;
  if(L <= 0){
    double sum = 0; int count = 0;
    for(const edge_t& e : edges){
      if(free_mask[e.first] || free_mask[e.second]) continue;
      sum += coord3d::dist(points[e.first], points[e.second]);
      count++;
    }
    if(count == 0)  // fallback: all edges
      for(const edge_t& e : edges){
        sum += coord3d::dist(points[e.first], points[e.second]);
        count++;
      }
    L = sum / count;
  }

  // Force constants: bond + angle + convexity bias, no curvature/flatness.
  const double k_bond = 1.0, k_angle = 1.0;
  const double k_curv = 0, k_flat = 0, k_conv = 5.0;
  // sigma_conv matches assemble_patch_hessian's softplus E_conv formulation.
  const double sigma_conv = 0.2 * L;

  // "No new concavities" constraint: compute h at entry for checked vertices.
  // Reject any step that makes a vertex with h_start > 0 have h_trial < 0.
  // Build a mask of checked vertices (free + interior).
  vector<bool> checked_mask(N, false);
  for(int v = 0; v < N; v++)
    checked_mask[v] = free_mask[v] || (!interior_mask.empty() && interior_mask[v]);
  vector<double> h_start;
  compute_h_values(*this, points, h_start);

  // Collect free vertex indices
  vector<int> free_idx;
  for(int i = 0; i < N; i++)
    if(free_mask[i]) free_idx.push_back(i);
  int nfree = (int)free_idx.size();
  int ndof = 3 * nfree;
  if(ndof == 0) return true;

  // Lambdas wrapping energy/gradient computation.
  // interior_mask restricts E_conv to vertices with full neighbor sets in the
  // patch — boundary vertices have truncated degree and produce bogus h values.
  auto energy_fn = [&](std::span<const coord3d> x) -> double {
    return deltahedron_energy_only(*this, edges, x, L, k_bond, k_angle, k_curv, k_flat, k_conv, sigma_conv, interior_mask);
  };
  auto grad_fn = [&](std::span<const coord3d> x, vector<coord3d>& g) -> double {
    return deltahedron_energy_and_gradient(*this, edges, x, &g, L, k_bond, k_angle, k_curv, k_flat, k_conv, sigma_conv, interior_mask);
  };

  // Extract free-vertex gradient components into flat vector
  auto extract_grad = [&](const vector<coord3d>& g, vector<double>& gf){
    gf.resize(ndof);
    for(int k = 0; k < nfree; k++)
      for(int c = 0; c < 3; c++)
        gf[3*k+c] = g[free_idx[k]][c];
  };

  vector<coord3d> grad(N);
  vector<double> gf(ndof);
  int consec_rejects = 0;

  // Trust region parameters
  double Delta_max = L;          // max trust region radius
  double Delta = 0.5 * L;       // initial trust region radius

  for(int iter = 0; iter < max_iter; iter++){
    iterations_used = iter + 1;

    double E = grad_fn(points, grad);
    extract_grad(grad, gf);

    double gnorm = 0;
    for(double v : gf) gnorm += v*v;
    gnorm = sqrt(gnorm);
    if(gnorm < grad_tol) return true;

    // Trust region stall detection: 5 consecutive rejected steps means
    // we've hit a gradient floor and can't make further progress.
    if(consec_rejects >= 5) return false;

    // Assemble exact analytical Hessian
    vector<vector<double>> H(ndof, vector<double>(ndof, 0.0));
    const double sigma = 0.2 * L;
    assemble_patch_hessian(*this, edges, points, H, free_idx,
                           L, k_bond, k_angle, k_conv, sigma, interior_mask);

    // --- Solve (H + lambda*I) delta = -g via Gaussian elimination ---
    auto solve_linear = [&](const vector<vector<double>>& M,
                            const vector<double>& rhs_in,
                            vector<double>& x_out) -> bool {
      vector<vector<double>> A = M;
      vector<double> rhs = rhs_in;
      x_out.resize(ndof);
      for(int col = 0; col < ndof; col++){
        int best = col;
        for(int row = col+1; row < ndof; row++)
          if(fabs(A[row][col]) > fabs(A[best][col])) best = row;
        swap(A[col], A[best]);
        swap(rhs[col], rhs[best]);
        if(fabs(A[col][col]) < 1e-30) return false;
        for(int row = col+1; row < ndof; row++){
          double f = A[row][col] / A[col][col];
          for(int j = col; j < ndof; j++) A[row][j] -= f * A[col][j];
          rhs[row] -= f * rhs[col];
        }
      }
      for(int i = ndof-1; i >= 0; i--){
        x_out[i] = rhs[i];
        for(int j = i+1; j < ndof; j++) x_out[i] -= A[i][j] * x_out[j];
        x_out[i] /= A[i][i];
      }
      return true;
    };

    auto solve_shifted = [&](double lambda, vector<double>& delta_out) -> bool {
      vector<vector<double>> A = H;
      for(int i = 0; i < ndof; i++) A[i][i] += lambda;
      vector<double> rhs(ndof);
      for(int i = 0; i < ndof; i++) rhs[i] = -gf[i];
      return solve_linear(A, rhs, delta_out);
    };

    auto vec_norm = [](const vector<double>& v) -> double {
      double s = 0; for(double x : v) s += x*x; return sqrt(s);
    };

    auto vec_dot = [](const vector<double>& a, const vector<double>& b) -> double {
      double s = 0; for(size_t i = 0; i < a.size(); i++) s += a[i]*b[i]; return s;
    };

    // LDL^T to count negative pivots and find most-negative pivot
    int n_neg_eig = 0;
    double min_pivot = 1e30;
    {
      vector<vector<double>> Lf(ndof, vector<double>(ndof, 0.0));
      vector<double> Df(ndof, 0.0);
      for(int i = 0; i < ndof; i++){
        double sum = H[i][i];
        for(int k = 0; k < i; k++) sum -= Lf[i][k]*Lf[i][k]*Df[k];
        Df[i] = sum;
        if(sum < min_pivot) min_pivot = sum;
        if(sum < 0) n_neg_eig++;
        double dabs = max(fabs(sum), 1e-30);
        for(int j = i+1; j < ndof; j++){
          double s = H[j][i];
          for(int k = 0; k < i; k++) s -= Lf[j][k]*Lf[i][k]*Df[k];
          Lf[j][i] = s / dabs;
        }
      }
    }

    // --- Trust-region subproblem: bisect on lambda to find ||delta(lambda)|| ~ Delta ---
    vector<double> delta(ndof);
    double dnorm = 0;
    double lambda = 0;
    const char* step_type = "newton";

    // Try unconstrained Newton first
    bool solved = solve_shifted(0, delta);
    dnorm = solved ? vec_norm(delta) : 1e30;
    double slope = solved ? vec_dot(gf, delta) : 1;

    if(!solved || dnorm > Delta || slope > 0){
      // Bisect on lambda to find (H+lambda*I)^{-1}g with ||delta|| ~ Delta
      double lambda_lo = 0, lambda_hi = gnorm / Delta + 1.0;
      for(int probe = 0; probe < 10; probe++){
        if(solve_shifted(lambda_hi, delta)){
          dnorm = vec_norm(delta);
          slope = vec_dot(gf, delta);
          if(dnorm <= Delta && slope <= 0) break;
        }
        lambda_hi *= 4;
      }
      for(int bis = 0; bis < 20; bis++){
        lambda = 0.5 * (lambda_lo + lambda_hi);
        if(!solve_shifted(lambda, delta)){
          lambda_lo = lambda; continue;
        }
        dnorm = vec_norm(delta);
        slope = vec_dot(gf, delta);
        if(dnorm > Delta || slope > 0) lambda_lo = lambda;
        else lambda_hi = lambda;
      }
      lambda = lambda_hi;
      solve_shifted(lambda, delta);
      dnorm = vec_norm(delta);
      step_type = "reg-newton";
    }

    // Predicted reduction from quadratic model: -(g'*delta + 0.5*delta'*H*delta)
    double pred = 0;
    for(int i = 0; i < ndof; i++){
      pred -= gf[i] * delta[i];
      for(int j = 0; j < ndof; j++)
        pred -= 0.5 * delta[i] * H[i][j] * delta[j];
    }

    // Trial point
    vector<coord3d> x_trial(points.begin(), points.end());
    for(int k = 0; k < nfree; k++)
      for(int c = 0; c < 3; c++)
        x_trial[free_idx[k]][c] = points[free_idx[k]][c] + delta[3*k+c];

    double Et = energy_fn(x_trial);
    double actual = E - Et;

    // "No new concavities" constraint: reject if any checked vertex that was
    // convex at entry (h_start > 0) becomes concave in the trial geometry.
    vector<double> h_trial;
    compute_h_values(*this, x_trial, h_trial);
    bool convex = true;
    if(convex_constraint)
      for(int v = 0; v < N && convex; v++)
        if(checked_mask[v] && h_start[v] > 0 && h_trial[v] < 0)
          convex = false;

    // Trust region update
    double rho = (pred > 0) ? actual / pred : -1;
    bool accepted = (rho > 0.1) && convex;
    const char* tr_action;
    if(accepted){
      consec_rejects = 0;
      for(int k = 0; k < nfree; k++)
        for(int c = 0; c < 3; c++)
          points[free_idx[k]][c] = x_trial[free_idx[k]][c];
      if(rho > 0.75 && dnorm > 0.5 * Delta){
        Delta = min(2.0 * Delta, Delta_max);
        tr_action = "expand";
      } else {
        tr_action = "keep";
      }
    } else {
      consec_rejects++;
      Delta *= 0.25;
      if(Delta < 1e-14 * L) Delta = 1e-14 * L;
      tr_action = convex ? "shrink" : "conv-shrink";
    }

    if(opt_log){
      // Compute h_min over interior vertices of the current (possibly trial) geometry
      std::span<const coord3d> xlog = accepted ? std::span<const coord3d>(points) : std::span<const coord3d>(x_trial);
      double h_min_log = 1e30;
      double ang_log = 0;
      for(int v = 0; v < N; v++){
        bool is_int = !interior_mask.empty() && interior_mask[v];
        if(!free_mask[v] && !is_int) continue;
        int d = degree(v);
        coord3d cen(0,0,0);
        for(int j = 0; j < d; j++) cen += xlog[(*this)[v][j]];
        cen /= (double)d;
        coord3d nf(0,0,0);
        for(int j = 0; j < d; j++){
          coord3d e1 = xlog[(*this)[v][j]] - xlog[v];
          coord3d e2 = xlog[(*this)[v][(j+1)%d]] - xlog[v];
          nf += e1.cross(e2);
        }
        double nl = nf.norm();
        if(nl < 1e-15) continue;
        double h = (xlog[v] - cen).dot(nf / nl);
        if(h < h_min_log) h_min_log = h;
      }
      ang_log = max_angle_relerr();
      fprintf(opt_log, "    patch %3d: E=%.6e |g|=%.3e |d|=%.2e D=%.2e rho=%.2f %-14s %-11s neg=%d h=%+.4f ang=%.2e\n",
              iter, E, gnorm, dnorm, Delta, rho, step_type,
              tr_action, n_neg_eig, h_min_log, ang_log);
    }
  }
  return false;  // didn't converge
}

template<>
int DeltahedronView<double>::reflect_concave(std::span<coord3d> pts, double threshold,
                                  const vector<bool>& fixed,
                                  vector<bool>* reflected_mask) const
{
  bool has_fixed = !fixed.empty();
  int count = 0;
  for(int v = 0; v < N; v++){
    if(has_fixed && fixed[v]) continue;
    int d = degree(v);
    if(d > 6) continue;

    auto vd = VertexHData::compute_h(*this, pts, v);
    if(!vd.valid) continue;
    if(vd.h < -threshold){
      pts[v] = vd.centroid + vd.n_hat * (-vd.h);
      if(reflected_mask) (*reflected_mask)[v] = true;
      count++;
    }
  }
  return count;
}

template<>
int DeltahedronView<double>::reflect_all_concave(std::span<coord3d> pts, double threshold,
                                      const vector<bool>& fixed,
                                      vector<bool>* reflected_mask) const
{
  int total = 0;
  for(int pass = 0; pass < 20; pass++){
    int n = reflect_concave(pts, threshold, fixed, reflected_mask);
    if(n == 0) break;
    total += n;
    if(pass == 19)
      fprintf(stderr, "WARNING: reflect_all_concave hit 20-pass limit (N=%d)\n", N);
  }
  return total;
}

// ============================================================
// Convex hull projection: compute the convex hull of the point set,
// then project every concave vertex onto the nearest hull face.
// ============================================================

// Incremental convex hull returning just the triangle list.
// Each triangle is an array of 3 vertex indices into the pts array.
// Triangles are oriented with outward normals.
static vector<array<int,3>> build_convex_hull(std::span<coord3d> pts)
{
  int n = (int)pts.size();
  if(n < 4) return {};

  // Find 4 non-coplanar points for initial tetrahedron
  int p0 = 0, p1 = -1, p2 = -1, p3 = -1;

  // Find p1: furthest from p0
  double best = 0;
  for(int i = 1; i < n; i++){
    double d = (pts[i] - pts[p0]).dot(pts[i] - pts[p0]);
    if(d > best){ best = d; p1 = i; }
  }
  if(p1 < 0) return {};

  // Find p2: furthest from line p0-p1
  coord3d u01 = pts[p1] - pts[p0];
  double u01_len2 = u01.dot(u01);
  best = 0;
  for(int i = 0; i < n; i++){
    if(i == p0 || i == p1) continue;
    coord3d v = pts[i] - pts[p0];
    double proj = v.dot(u01) / u01_len2;
    coord3d perp = v - u01 * proj;
    double d = perp.dot(perp);
    if(d > best){ best = d; p2 = i; }
  }
  if(p2 < 0) return {};

  // Find p3: furthest from plane p0-p1-p2
  coord3d normal = (pts[p1] - pts[p0]).cross(pts[p2] - pts[p0]);
  double nlen = normal.norm();
  if(nlen < 1e-15) return {};
  normal /= nlen;
  best = 0;
  for(int i = 0; i < n; i++){
    if(i == p0 || i == p1 || i == p2) continue;
    double d = fabs((pts[i] - pts[p0]).dot(normal));
    if(d > best){ best = d; p3 = i; }
  }
  if(p3 < 0) return {};

  // Build initial tetrahedron with outward-facing triangles
  coord3d centroid = (pts[p0] + pts[p1] + pts[p2] + pts[p3]) / 4.0;

  auto make_outward = [&](int a, int b, int c) -> array<int,3> {
    coord3d fn = (pts[b] - pts[a]).cross(pts[c] - pts[a]);
    coord3d fc = (pts[a] + pts[b] + pts[c]) / 3.0 - centroid;
    if(fn.dot(fc) < 0) return {a, c, b};  // flip
    return {a, b, c};
  };

  vector<array<int,3>> tris;
  tris.push_back(make_outward(p0, p1, p2));
  tris.push_back(make_outward(p0, p1, p3));
  tris.push_back(make_outward(p0, p2, p3));
  tris.push_back(make_outward(p1, p2, p3));

  // Track which vertices are on the hull
  vector<bool> on_hull(n, false);
  on_hull[p0] = on_hull[p1] = on_hull[p2] = on_hull[p3] = true;

  // Incrementally add each remaining vertex
  for(int i = 0; i < n; i++){
    if(on_hull[i]) continue;
    const coord3d& p = pts[i];

    // Find visible faces (point is in front of the face)
    vector<int> visible;
    for(int f = 0; f < (int)tris.size(); f++){
      coord3d fn = (pts[tris[f][1]] - pts[tris[f][0]]).cross(pts[tris[f][2]] - pts[tris[f][0]]);
      if(fn.dot(p - pts[tris[f][0]]) > 1e-14 * fn.norm())
        visible.push_back(f);
    }
    if(visible.empty()) continue;  // inside hull

    on_hull[i] = true;

    // Build directed edge -> face map for visible faces
    map<pair<int,int>, int> edge_face;
    for(int fi : visible){
      auto& t = tris[fi];
      for(int j = 0; j < 3; j++)
        edge_face[{t[j], t[(j+1)%3]}] = fi;
    }

    // Horizon edges: directed edges of visible faces whose reverse is NOT in a visible face
    vector<pair<int,int>> horizon;
    for(auto& [edge, fi] : edge_face){
      if(edge_face.find({edge.second, edge.first}) == edge_face.end())
        horizon.push_back(edge);
    }

    // Remove visible faces (in reverse order to preserve indices)
    sort(visible.rbegin(), visible.rend());
    for(int fi : visible)
      tris.erase(tris.begin() + fi);

    // Add new faces connecting point i to each horizon edge
    for(auto& [a, b] : horizon)
      tris.push_back({i, a, b});  // a,b is CCW from outside of invisible neighbor,
                                   // so i,a,b has outward normal (away from hull interior)
  }

  return tris;
}

template<>
int DeltahedronView<double>::project_onto_convex_hull(std::span<coord3d> pts) const
{
  // 1. Identify concave vertices
  vector<int> concave;
  for(int v = 0; v < N; v++){
    auto vd = VertexHData::compute_h(*this, pts, v);
    if(vd.valid && vd.h < 0) concave.push_back(v);
  }
  if(concave.empty()) return 0;

  // 2. Build convex hull
  auto tris = build_convex_hull(pts);
  if(tris.empty()) return 0;

  // 3. For each concave vertex, project onto nearest hull face. closest_point_on_triangle
  // reports dist2 = +inf for a degenerate (zero-area) face; only move the vertex if some
  // non-degenerate face was found, else leave it in place (never snap it to the origin).
  for(int v : concave){
    double best_dist2 = std::numeric_limits<double>::infinity();
    coord3d best_pt = pts[v];
    for(auto& tri : tris){
      ClosestPoint cp = closest_point_on_triangle(pts[v], tri_t(tri[0],tri[1],tri[2]), pts);
      if(cp.dist2 < best_dist2){ best_dist2 = cp.dist2; best_pt = cp.pos; }
    }
    pts[v] = best_pt;
  }

  return (int)concave.size();
}

const char* opt_result_name(OptResult r) {
  switch(r) {
    case OptResult::CONVERGED:        return "CONVERGED";
    case OptResult::STAGNATED:        return "STAGNATED";
    case OptResult::BUDGET_EXHAUSTED: return "BUDGET_EXHAUSTED";
    case OptResult::CONVEXITY_STUCK:  return "CONVEXITY_STUCK";
  }
  return "UNKNOWN";
}

const char* PipelineDiag::flag_name(uint32_t f) {
  switch(f) {
    case REFLECT_CYCLING_STEP:  return "reflect_cycling_step";
    case HULL_USED_STEP:        return "hull_used_step";
    case HULL_CYCLING_STEP:     return "hull_cycling_step";
    case CONVEXITY_FAIL_STEP:   return "cvx_fail_step";
    case PATCH_CYCLING:         return "patch_cycling";
    case STAG_STEP:             return "stag_step";
    case REFLECT_CYCLING_FINAL: return "reflect_cycling_final";
    case HULL_USED_FINAL:       return "hull_used_final";
    case HULL_CYCLING_FINAL:    return "hull_cycling_final";
    case STAG_FINAL:            return "stag_final";
    case STAG_CONSTRAINED:      return "stag_constrained";
    case BUDGET_CONSTRAINED:    return "budget_constrained";
    case NEG_CURVATURE:         return "neg_curvature";
    case TR_BOUNDARY:           return "tr_boundary";
    case STEP_REJECTED:         return "step_rejected";
    case CVX_REJECTED:          return "cvx_rejected";
    case LBFGS_RESET:           return "lbfgs_reset";
    case HAS_F_RING:            return "has_f_ring";
  }
  return "unknown";
}

OptResult Deltahedron::optimize(std::span<const coord3d> initial_geometry, double target_L,
                                double grad_tol, const vector<bool>& fixed,
                                long long max_work, double angle_tol)
{
  assert((int)initial_geometry.size() == N);
  assert(fixed.empty() || (int)fixed.size() == N);
  std::copy(initial_geometry.begin(), initial_geometry.end(), points.begin());
  const bool has_fixed = !fixed.empty();
  opt_diag_flags = 0;  // Reset per-call optimizer diagnostics

  // Cache edge list (avoid recomputing on every energy evaluation)
  vector<edge_t> edges = undirected_edges();

  // Target edge length: use provided value, or mean of initial edge lengths.
  double L = target_L;
  if(L <= 0){
    double sum = 0; int count = 0;
    for(const edge_t& e : edges){
      sum += coord3d::dist(points[e.first], points[e.second]);
      count++;
    }
    L = sum / count;
  }

  // Force constants (fixed throughout)
  const double k_bond   = 1.0;
  const double k_angle  = 1.0;
  const double k_curv   = 2.0;
  const double k_conv   = opt_k_conv;  // quadratic one-sided convexity penalty

  // Two-phase optimization:
  //   Phase 1: k_flat active — settle into flat/equilateral
  //   Phase 2: k_flat off — pure equilateral convergence
  double k_flat = opt_k_flat;
  int phase = (k_flat > 0) ? 1 : 2;
  double phase1_grad_norm0 = 0;

  const double c1 = 1e-4;        // Armijo parameter

  // Shared helpers
  auto zero_fixed_grad = [&](vector<coord3d>& grad){
    if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) grad[i] = coord3d(0,0,0);
  };

  // Evaluation counters
  n_energy_evals = 0;
  n_grad_evals = 0;
  n_hv_evals = 0;

  // Work budget: total_work = n_energy + 2*n_grad + 7*n_hv
  // All three primitives are O(Nv) per call with constant cost ratios.
  // Empirical cost ratios measured at Nv=32,52,102 (consistent within 10%):
  //   energy_eval : gradient_eval : hv_product ≈ 1 : 2 : 7
  // Budget in units of "energy evaluations".  Default: 400*Nv^2.
  // Each CG/LBFGS iteration costs ~17 energy evals (line search) + 1 gradient (=2),
  // so ~19 work per iteration.  400*Nv^2/19 ≈ 21*Nv^2 iterations.
  // Real wall time scales as budget * Nv (since each eval is O(Nv)).
  if(max_work <= 0) max_work = 400LL * N * N;
  const long long phase1_work_budget = max_work / 4;
  auto total_work_fn = [&]() -> long long {
    return (long long)n_energy_evals + 2LL * n_grad_evals + 7LL * n_hv_evals;
  };

  auto compute_eg = [&](vector<coord3d>& grad) -> double {
    n_grad_evals++;
    double E = deltahedron_energy_and_gradient(*this, edges, points, &grad, L,
                                                k_bond, k_angle, k_curv, k_flat, k_conv);
    zero_fixed_grad(grad);
    return E;
  };

  auto compute_e_only = [&](const vector<coord3d>& x_trial) -> double {
    n_energy_evals++;
    return deltahedron_energy_only(*this, edges, x_trial, L, k_bond, k_angle, k_curv, k_flat, k_conv);
  };

  auto compute_gmax = [&](const vector<coord3d>& grad) -> double {
    double gmax = 0;
    for(int i = 0; i < N; i++){
      if(has_fixed && fixed[i]) continue;
      gmax = max(gmax, grad[i].norm());
    }
    return gmax * L;
  };

  auto edge_cv = [&]() -> double {
    double s = 0, s2 = 0;
    for(const edge_t& e : edges){
      double d = coord3d::dist(points[e.first], points[e.second]);
      s += d; s2 += d*d;
    }
    int ne = (int)edges.size();
    double mu = s / ne;
    return sqrt(max(0.0, s2/ne - mu*mu)) / mu;
  };

  // Phase transition logic (shared by all methods)
  // Phase 1 ends when: work budget quarter exhausted, or gradient drops 100x.
  // Returns true if phase changed.
  long long phase1_work_start = 0;
  auto check_phase_transition = [&](int iter, const vector<coord3d>& grad) -> bool {
    if(phase != 1 || iter % 50 != 49) return false;
    bool advance = (total_work_fn() - phase1_work_start >= phase1_work_budget);
    if(!advance && phase1_grad_norm0 > 0){
      double gn = vec_norm(grad);
      if(gn < phase1_grad_norm0 * 0.01) advance = true;
    }
    if(advance){
      k_flat = 0;
      phase = 2;
      return true;
    }
    return false;
  };

  // Backtracking Armijo line search (shared by CG and L-BFGS)
  auto line_search = [&](double E, const vector<coord3d>& grad, const vector<coord3d>& dir,
                         vector<coord3d>& x_trial) -> double {
    double slope = vec_dot(grad, dir);
    double alpha = 1.0;
    for(int ls = 0; ls < 60; ls++){
      for(int i = 0; i < N; i++) x_trial[i] = points[i] + dir[i] * alpha;
      double E_trial = compute_e_only(x_trial);
      if(E_trial <= E + c1 * alpha * slope) break;
      double denom = 2.0 * (E_trial - E - slope * alpha);
      if(denom > 1e-30){
        double alpha_q = -slope * alpha * alpha / denom;
        alpha = max(0.1 * alpha, min(0.5 * alpha, alpha_q));
      } else {
        alpha *= 0.5;
      }
    }
    return alpha;
  };

  // Log ~20 times over the work budget. Estimate iters from work/N (each iter ~ N work).
  const int log_interval = opt_log ? max(1, (int)(max_work / (20LL * N))) : 0;
  const char* method_name = (opt_method == OptMethod::CG) ? "CG" :
                            (opt_method == OptMethod::LBFGS) ? "LBFGS" : "ST";

  vector<coord3d> grad(N);
  double E = compute_eg(grad);
  phase1_grad_norm0 = vec_norm(grad);

  if(opt_log)
    fprintf(opt_log, "  %s start: E=%.6f |g|=%.4e L=%.4f cv=%.4f ph=%d tol=%.2e\n",
            method_name, E, phase1_grad_norm0, L, edge_cv(), phase, grad_tol);

  bool converged = false;

  // Stagnation detection for angle-based convergence: if energy hasn't
  // decreased meaningfully for stag_window consecutive iterations, the
  // optimizer is stuck at a local minimum and can't reduce angle error
  // further.  Break to avoid burning the entire work budget.
  const int stag_window = 50;
  double stag_E_ref = E;   // energy at start of current window
  int stag_count = 0;      // iterations since last meaningful decrease

  // ==================== CG ====================
  if(opt_method == OptMethod::CG){
    vector<coord3d> grad_old(N), dir(N), x_trial(N);
    for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0);
    if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir[i] = coord3d(0,0,0);

    for(int iter = 0; ; iter++){
      iterations_used = iter + 1;
      if(total_work_fn() >= max_work) break;

      // Phase transition
      if(check_phase_transition(iter, grad)){
        E = compute_eg(grad);
        for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0);
        if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir[i] = coord3d(0,0,0);
      }

      // Convergence check
      double gmax = compute_gmax(grad);
      if(gmax < grad_tol){ converged = true; break; }
      if(angle_tol > 0 && max_angle_relerr() < angle_tol && count_concave() == 0){ converged = true; break; }
      if(stag_count >= stag_window) break;  // no progress for 50 iterations

      // Ensure descent direction
      double slope = vec_dot(grad, dir);
      if(slope > 0){
        for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0);
        if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir[i] = coord3d(0,0,0);
      }

      // Line search and update
      double alpha = line_search(E, grad, dir, x_trial);
      for(int i = 0; i < N; i++) points[i] = points[i] + dir[i] * alpha;

      grad_old = grad;
      double E_old = E;
      E = compute_eg(grad);

      // Stagnation tracking
      if(E_old - E > 1e-15 * max(1.0, fabs(E_old))){ stag_count = 0; stag_E_ref = E; }
      else stag_count++;

      // Logging
      if(log_interval > 0 && iter % log_interval == 0)
        fprintf(opt_log, "  CG %4d: E=%.6f |g|=%.4e gmax*L=%.4e a=%.3e cv=%.4f ang=%.2e ph=%d\n",
                iter, E, vec_norm(grad), compute_gmax(grad), alpha, edge_cv(),
                max_angle_relerr(), phase);

      // Polak-Ribiere beta
      double gg_old = vec_dot(grad_old, grad_old);
      double beta = 0.0;
      if(gg_old > 1e-30){
        vector<coord3d> gdiff(N);
        for(int i = 0; i < N; i++) gdiff[i] = grad[i] - grad_old[i];
        beta = max(0.0, vec_dot(grad, gdiff) / gg_old);
      }
      for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0) + dir[i] * beta;
      if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir[i] = coord3d(0,0,0);
    }
  }

  // ==================== L-BFGS ====================
  else if(opt_method == OptMethod::LBFGS){
    const int m = 10;  // history depth
    deque<vector<coord3d>> S, Y;
    deque<double> rho_hist;
    vector<coord3d> dir(N), x_trial(N), grad_old(N);

    for(int iter = 0; ; iter++){
      iterations_used = iter + 1;
      if(total_work_fn() >= max_work) break;

      // Phase transition
      if(check_phase_transition(iter, grad)){
        E = compute_eg(grad);
        S.clear(); Y.clear(); rho_hist.clear();
      }

      // Convergence check
      double gmax = compute_gmax(grad);
      if(gmax < grad_tol){ converged = true; break; }
      if(angle_tol > 0 && max_angle_relerr() < angle_tol && count_concave() == 0){ converged = true; break; }
      if(stag_count >= stag_window) break;  // no progress for 50 iterations

      // Two-loop recursion to compute search direction
      dir = grad;  // q = grad
      int hist_size = (int)S.size();
      vector<double> alpha_hist(hist_size);

      // Forward loop
      for(int i = hist_size - 1; i >= 0; i--){
        alpha_hist[i] = rho_hist[i] * vec_dot(S[i], dir);
        vec_axpy(dir, -alpha_hist[i], Y[i]);
      }

      // Initial Hessian approximation: gamma * I
      if(hist_size > 0){
        double ys = vec_dot(Y.back(), S.back());
        double yy = vec_dot(Y.back(), Y.back());
        if(yy > 1e-30) vec_scale(dir, ys / yy);
      }

      // Backward loop
      for(int i = 0; i < hist_size; i++){
        double beta = rho_hist[i] * vec_dot(Y[i], dir);
        vec_axpy(dir, alpha_hist[i] - beta, S[i]);
      }

      // dir = -H*g (negate for descent direction)
      vec_scale(dir, -1.0);
      if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir[i] = coord3d(0,0,0);

      // Safeguard: if not a descent direction, reset to -grad
      double slope = vec_dot(grad, dir);
      if(slope > 0){
        opt_diag_flags |= PipelineDiag::LBFGS_RESET;
        for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0);
        if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir[i] = coord3d(0,0,0);
        S.clear(); Y.clear(); rho_hist.clear();
      }

      // Save old state for history update
      grad_old = grad;

      // Line search and update
      double alpha = line_search(E, grad, dir, x_trial);
      for(int i = 0; i < N; i++) points[i] = points[i] + dir[i] * alpha;
      double E_old = E;
      E = compute_eg(grad);

      // Stagnation tracking
      if(E_old - E > 1e-15 * max(1.0, fabs(E_old))){ stag_count = 0; stag_E_ref = E; }
      else stag_count++;

      // Update L-BFGS history: s = alpha*dir, y = grad - grad_old
      {
        vector<coord3d> s(N), y(N);
        for(int i = 0; i < N; i++){
          s[i] = dir[i] * alpha;
          y[i] = grad[i] - grad_old[i];
        }
        double ys = vec_dot(y, s);
        if(ys > 1e-10 * vec_dot(s, s)){  // curvature condition
          S.push_back(std::move(s));
          Y.push_back(std::move(y));
          rho_hist.push_back(1.0 / ys);
          if((int)S.size() > m){ S.pop_front(); Y.pop_front(); rho_hist.pop_front(); }
        }
      }

      // Logging
      if(log_interval > 0 && iter % log_interval == 0){
        // Compute h_min for log
        double h_min_log = 1e30;
        for(int v = 0; v < N; v++){
          if(has_fixed && fixed[v]) continue;
          int d = degree(v);
          if(d > 6) continue;
          coord3d cen(0,0,0);
          for(int j = 0; j < d; j++) cen += points[(*this)[v][j]];
          cen /= (double)d;
          coord3d nf(0,0,0);
          for(int j = 0; j < d; j++){
            coord3d e1 = points[(*this)[v][j]] - points[v];
            coord3d e2 = points[(*this)[v][(j+1)%d]] - points[v];
            nf += e1.cross(e2);
          }
          double nl = nf.norm();
          if(nl < 1e-15) continue;
          double h = (points[v] - cen).dot(nf / nl);
          if(h < h_min_log) h_min_log = h;
        }
        fprintf(opt_log, "  LB %4d: E=%.6f |g|=%.4e gmax*L=%.4e a=%.3e cv=%.4f ang=%.2e hmin=%+.4f ph=%d h=%d\n",
                iter, E, vec_norm(grad), compute_gmax(grad), alpha, edge_cv(),
                max_angle_relerr(), h_min_log, phase, (int)S.size());
      }
    }
  }

  // ==================== Steihaug-Toint ====================
  else if(opt_method == OptMethod::STEIHAUG){
    double Delta_max = L;
    double Delta = 0.5 * L;
    const int max_inner = min(3 * N, 200);

    // Temp vectors for inner CG
    vector<coord3d> z(N), r_cg(N), d_cg(N), Hd(N), x_trial(N);

    // Convexity constraint: h_current tracks which vertices are convex.
    // Updated at each accepted step; only convex→concave transitions are rejected.
    // Works with or without E_conv: constraint prevents new concavities,
    // E_conv (if active) pushes existing concave vertices toward h=0.
    vector<double> h_current, h_trial;
    if(opt_convex_constraint){
      compute_h_values(*this, points, h_current, fixed);
    }

    for(int iter = 0; ; iter++){
      iterations_used = iter + 1;
      if(total_work_fn() >= max_work) break;

      // Phase transition
      if(check_phase_transition(iter, grad)){
        E = compute_eg(grad);
        Delta = 0.5 * L;  // reset trust region on phase change
        if(opt_convex_constraint) compute_h_values(*this, points, h_current, fixed);
      }

      // Convergence check
      double gmax = compute_gmax(grad);
      if(gmax < grad_tol){ converged = true; break; }
      if(angle_tol > 0 && max_angle_relerr() < angle_tol && count_concave() == 0){ converged = true; break; }
      if(stag_count >= stag_window) break;  // no progress for 50 iterations

      // --- Steihaug CG to solve trust-region subproblem ---
      // Approximately solve: min_z  g^T z + 0.5 z^T H z,  ||z|| <= Delta
      vec_zero(z);
      for(int i = 0; i < N; i++) r_cg[i] = grad[i] * (-1.0);  // r = -g
      if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) r_cg[i] = coord3d(0,0,0);
      d_cg = r_cg;  // d = r
      double rr = vec_dot(r_cg, r_cg);
      double gnorm = sqrt(rr);
      int inner_iters = 0;

      for(int j = 0; j < max_inner; j++){
        inner_iters = j + 1;

        // Hv product
        vec_zero(Hd);
        n_hv_evals++;
        deltahedron_hv_product(*this, edges, points, d_cg, Hd, L,
                               k_bond, k_angle, k_curv, k_flat, fixed, k_conv);

        double kappa = vec_dot(d_cg, Hd);

        if(kappa <= 1e-15 * rr){
          // Negative or zero curvature: step to trust-region boundary along d
          opt_diag_flags |= PipelineDiag::NEG_CURVATURE;
          double zz = vec_dot(z, z);
          double zd = vec_dot(z, d_cg);
          double dd = vec_dot(d_cg, d_cg);
          // Solve ||z + tau*d||^2 = Delta^2
          double a_coeff = dd;
          double b_coeff = 2.0 * zd;
          double c_coeff = zz - Delta * Delta;
          double disc = b_coeff * b_coeff - 4.0 * a_coeff * c_coeff;
          double tau = (-b_coeff + sqrt(max(0.0, disc))) / (2.0 * a_coeff);
          vec_axpy(z, tau, d_cg);
          break;
        }

        double alpha_cg = rr / kappa;

        // Check trust-region boundary
        {
          vector<coord3d> z_new(N);
          for(int i = 0; i < N; i++) z_new[i] = z[i] + d_cg[i] * alpha_cg;
          double z_new_norm = vec_norm(z_new);
          if(z_new_norm >= Delta){
            // Truncate to boundary
            opt_diag_flags |= PipelineDiag::TR_BOUNDARY;
            double zz = vec_dot(z, z);
            double zd = vec_dot(z, d_cg);
            double dd = vec_dot(d_cg, d_cg);
            double a_coeff = dd;
            double b_coeff = 2.0 * zd;
            double c_coeff = zz - Delta * Delta;
            double disc = b_coeff * b_coeff - 4.0 * a_coeff * c_coeff;
            double tau = (-b_coeff + sqrt(max(0.0, disc))) / (2.0 * a_coeff);
            vec_axpy(z, tau, d_cg);
            break;
          }
          z = z_new;
        }

        // Update residual
        vector<coord3d> r_new(N);
        for(int i = 0; i < N; i++) r_new[i] = r_cg[i] - Hd[i] * alpha_cg;
        if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) r_new[i] = coord3d(0,0,0);

        double rr_new = vec_dot(r_new, r_new);
        if(sqrt(rr_new) < 0.01 * gnorm) break;  // inner convergence

        double beta_cg = rr_new / rr;
        for(int i = 0; i < N; i++) d_cg[i] = r_new[i] + d_cg[i] * beta_cg;
        if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) d_cg[i] = coord3d(0,0,0);
        r_cg = r_new;
        rr = rr_new;
      }

      // --- Evaluate trust-region step ---
      double znorm = vec_norm(z);

      // Predicted reduction: -(g.z) - 0.5*(z.Hz)
      vec_zero(Hd);
      n_hv_evals++;
      deltahedron_hv_product(*this, edges, points, z, Hd, L,
                             k_bond, k_angle, k_curv, k_flat, fixed, k_conv);
      double pred = -(vec_dot(grad, z) + 0.5 * vec_dot(z, Hd));

      // Trial point
      for(int i = 0; i < N; i++) x_trial[i] = points[i] + z[i];
      double E_trial = compute_e_only(x_trial);
      double actual = E - E_trial;

      double rho = (pred > 1e-30) ? actual / pred : -1;

      // Accept based on energy reduction.  If opt_convex_constraint is on,
      // also reject steps that make a currently-convex vertex concave.
      bool accepted = (rho > 0.1);
      if(accepted && opt_convex_constraint){
        compute_h_values(*this, x_trial, h_trial, fixed);
        for(int v = 0; v < N && accepted; v++)
          if(h_current[v] > 0 && h_trial[v] < 0){
            accepted = false;
            opt_diag_flags |= PipelineDiag::CVX_REJECTED;
          }
      }

      if(accepted){
        double E_old = E;
        std::copy(x_trial.begin(), x_trial.end(), points.begin());
        E = compute_eg(grad);
        if(opt_convex_constraint) compute_h_values(*this, points, h_current, fixed);
        if(rho > 0.75 && znorm > 0.5 * Delta) Delta = min(2.0 * Delta, Delta_max);
        // Stagnation tracking
        if(E_old - E > 1e-15 * max(1.0, fabs(E_old))){ stag_count = 0; stag_E_ref = E; }
        else stag_count++;
      } else {
        Delta *= 0.25;
        if(Delta < 1e-14 * L) Delta = 1e-14 * L;
        stag_count++;  // rejected step = no progress
        opt_diag_flags |= PipelineDiag::STEP_REJECTED;
      }

      // Logging
      if(log_interval > 0 && iter % log_interval == 0)
        fprintf(opt_log, "  ST %4d: E=%.6f |g|=%.4e gmax*L=%.4e |z|=%.2e D=%.2e rho=%.2f ang=%.2e in=%d ph=%d %s\n",
                iter, E, vec_norm(grad), compute_gmax(grad), znorm, Delta, rho,
                max_angle_relerr(), inner_iters, phase, accepted ? "acc" : "REJ");
    }
  }

  // Final stats
  final_gmax_L = compute_gmax(grad);
  bool stagnated = !converged && stag_count >= stag_window;
  if(opt_log)
    fprintf(opt_log, "  %s done: %d iters, E=%.6f gmax*L=%.4e cv=%.4f %s\n",
            method_name, iterations_used, E, final_gmax_L, edge_cv(),
            converged ? "CONVERGED" : stagnated ? "STAGNATED" : "budget");

  // Post-optimization strict convexity cleanup.
  // Hull projection can disturb angle quality, so CG polish after projecting.
  // Skipped when opt_skip_post_reflect is set (caller will handle convexity).
  if(!opt_skip_post_reflect)
  {
    int total_projected = project_onto_convex_hull(points);
    if(total_projected > 0){
      // Projection moved vertices — brief CG (Polak-Ribiere) polish to recover angle quality.
      E = compute_eg(grad);
      zero_fixed_grad(grad);
      vector<coord3d> dir_r(N), grad_old_r(N), x_trial_r(N);
      for(int i = 0; i < N; i++) dir_r[i] = grad[i] * (-1.0);
      if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir_r[i] = coord3d(0,0,0);
      for(int iter = 0; iter < 50; iter++){
        if(compute_gmax(grad) < grad_tol) break;
        grad_old_r = grad;
        double alpha = line_search(E, grad, dir_r, x_trial_r);
        for(int i = 0; i < N; i++) points[i] = points[i] + dir_r[i] * alpha;
        E = compute_eg(grad);
        // Polak-Ribiere beta
        double gg_old = vec_dot(grad_old_r, grad_old_r);
        double beta = 0.0;
        if(gg_old > 1e-30){
          vector<coord3d> gdiff(N);
          for(int i = 0; i < N; i++) gdiff[i] = grad[i] - grad_old_r[i];
          beta = max(0.0, vec_dot(grad, gdiff) / gg_old);
        }
        for(int i = 0; i < N; i++) dir_r[i] = grad[i] * (-1.0) + dir_r[i] * beta;
        if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir_r[i] = coord3d(0,0,0);
      }
      // Final projection in case CG polish re-introduced barely-concave vertices
      project_onto_convex_hull(points);
      if(opt_log)
        fprintf(opt_log, "  Post-project polish: projected=%d ang=%.4e\n",
                total_projected, max_angle_relerr());
    }
  }

  final_angle_relerr = max_angle_relerr();
  final_n_concave = count_concave();
  final_opt_result = converged ? OptResult::CONVERGED
                   : stagnated ? OptResult::STAGNATED
                               : OptResult::BUDGET_EXHAUSTED;

  return final_opt_result;
}

template<>
double DeltahedronView<double>::gradient_check(std::span<const coord3d> geometry, double target_L, double eps) const
{
  vector<coord3d> x(geometry.begin(), geometry.end());
  vector<edge_t> edges = undirected_edges();

  double L = target_L;
  if(L <= 0){
    for(const edge_t& e : edges)
      L += coord3d::dist(x[e.first], x[e.second]);
    L /= edges.size();
  }

  // Same force constants as optimize()
  const double k_bond  = 1.0;
  const double k_angle = 1.0;
  const double k_curv  = 2.0;
  const double k_flat  = 2.0;
  const double k_conv  = 10.0;

  // Analytic gradient
  vector<coord3d> grad(N, coord3d(0,0,0));
  double E0 = deltahedron_energy_and_gradient(*this, edges, x, &grad, L,
                                               k_bond, k_angle, k_curv, k_flat, k_conv);

  // Finite-difference gradient
  double max_rel_err = 0;
  for(int i = 0; i < N; i++){
    for(int c = 0; c < 3; c++){
      vector<coord3d> xp = x, xm = x;
      xp[i][c] += eps;
      xm[i][c] -= eps;
      double Ep = deltahedron_energy_and_gradient(*this, edges, xp, nullptr, L,
                                                   k_bond, k_angle, k_curv, k_flat, k_conv);
      double Em = deltahedron_energy_and_gradient(*this, edges, xm, nullptr, L,
                                                   k_bond, k_angle, k_curv, k_flat, k_conv);
      double fd = (Ep - Em) / (2 * eps);
      double an = grad[i][c];
      double denom = max(1.0, max(abs(fd), abs(an)));
      double rel_err = abs(fd - an) / denom;
      if(rel_err > max_rel_err) max_rel_err = rel_err;
    }
  }
  return max_rel_err;
}

template<>
double DeltahedronView<double>::hessian_check(std::span<const coord3d> geometry,
                                  const vector<bool>& free_mask,
                                  const vector<bool>& interior_mask,
                                  double target_L, double eps, bool verbose) const
{
  vector<coord3d> x(geometry.begin(), geometry.end());
  vector<edge_t> edges = undirected_edges();

  // Same force constants as optimize_patch
  const double k_bond = 1.0, k_angle = 1.0, k_conv = 5.0;
  const double k_curv = 0, k_flat = 0;

  // Target L from boundary edges (same logic as optimize_patch)
  double L = target_L;
  if(L <= 0){
    double sum = 0; int count = 0;
    for(const edge_t& e : edges){
      if(free_mask[e.first] || free_mask[e.second]) continue;
      sum += coord3d::dist(x[e.first], x[e.second]);
      count++;
    }
    if(count == 0)
      for(const edge_t& e : edges){
        sum += coord3d::dist(x[e.first], x[e.second]);
        count++;
      }
    L = sum / count;
  }
  const double sigma = 0.2 * L;

  // Collect free vertex indices
  vector<int> free_idx;
  for(int i = 0; i < N; i++)
    if(free_mask[i]) free_idx.push_back(i);
  int nfree = (int)free_idx.size();
  int ndof = 3 * nfree;

  // Analytical Hessian
  vector<vector<double>> H_an(ndof, vector<double>(ndof, 0.0));
  assemble_patch_hessian(*this, edges, x, H_an, free_idx,
                         L, k_bond, k_angle, k_conv, sigma, interior_mask);

  // FD Hessian: differentiate gradient w.r.t. each free DOF
  auto grad_fn = [&](const vector<coord3d>& xx, vector<double>& gf){
    vector<coord3d> g(N, coord3d(0,0,0));
    deltahedron_energy_and_gradient(*this, edges, xx, &g, L,
                                    k_bond, k_angle, k_curv, k_flat, k_conv, sigma, interior_mask);
    gf.resize(ndof);
    for(int k = 0; k < nfree; k++)
      for(int c = 0; c < 3; c++)
        gf[3*k+c] = g[free_idx[k]][c];
  };

  vector<vector<double>> H_fd(ndof, vector<double>(ndof, 0.0));
  for(int j = 0; j < ndof; j++){
    int v = free_idx[j/3];
    int c = j % 3;
    vector<coord3d> xp = x, xm = x;
    xp[v][c] += eps;
    xm[v][c] -= eps;
    vector<double> gp, gm;
    grad_fn(xp, gp);
    grad_fn(xm, gm);
    for(int i = 0; i < ndof; i++)
      H_fd[i][j] = (gp[i] - gm[i]) / (2 * eps);
  }

  // Compare
  double max_rel_err = 0;
  int worst_i = 0, worst_j = 0;
  for(int i = 0; i < ndof; i++){
    for(int j = 0; j < ndof; j++){
      double denom = max(1.0, max(abs(H_an[i][j]), abs(H_fd[i][j])));
      double rel_err = abs(H_an[i][j] - H_fd[i][j]) / denom;
      if(rel_err > max_rel_err){
        max_rel_err = rel_err;
        worst_i = i; worst_j = j;
      }
    }
  }

  if(verbose){
    fprintf(stderr, "Hessian check: ndof=%d, L=%.4f, eps=%.1e, max_rel_err=%.3e\n",
            ndof, L, eps, max_rel_err);
    fprintf(stderr, "  worst at H[%d][%d]: analytical=%.8e  FD=%.8e\n",
            worst_i, worst_j, H_an[worst_i][worst_j], H_fd[worst_i][worst_j]);
    // Print top-10 worst entries
    vector<tuple<double,int,int>> errs;
    for(int i = 0; i < ndof; i++)
      for(int j = 0; j < ndof; j++){
        double denom = max(1.0, max(abs(H_an[i][j]), abs(H_fd[i][j])));
        errs.push_back({abs(H_an[i][j] - H_fd[i][j]) / denom, i, j});
      }
    sort(errs.rbegin(), errs.rend());
    fprintf(stderr, "  top-10 errors:\n");
    for(int k = 0; k < min(10, (int)errs.size()); k++){
      auto [err, i, j] = errs[k];
      fprintf(stderr, "    H[%d][%d]: an=% .6e  fd=% .6e  err=%.3e\n",
              i, j, H_an[i][j], H_fd[i][j], err);
    }
  }

  return max_rel_err;
}
