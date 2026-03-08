#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/unfold.hh"
#include "fullerenes/buckinverse.hh"
#include <cmath>
#include <climits>
#include <numeric>
#include <queue>
#include <set>
#include <chrono>

Deltahedron::Deltahedron(const Triangulation& T, std::span<const coord3d> pts)
  : Triangulation(T), points(vector<coord3d>(pts.begin(), pts.end()))
{
  assert((int)points.size() == N);
}

Deltahedron::Deltahedron(const Polyhedron& P)
  : Triangulation(static_cast<const Graph&>(P)),
    points(vector<coord3d>(P.points.begin(), P.points.end()))
{
  assert(P.is_triangulation());
}

vector<face_t> Deltahedron::compute_dual_faces() const {
  auto tris = triangles();
  vector<face_t> faces(tris.size());
  for(size_t i = 0; i < tris.size(); i++)
    faces[i] = face_t(tris[i]);
  return faces;
}

double Deltahedron::max_angle_relerr() const {
  double max_re = 0;
  const double target = M_PI / 3.0;
  for (const auto& t : triangles()) {
    for (int c = 0; c < 3; c++) {
      coord3d ea = points[t[(c+1)%3]] - points[t[c]];
      coord3d eb = points[t[(c+2)%3]] - points[t[c]];
      double cos_th = ea.dot(eb) / (ea.norm() * eb.norm());
      cos_th = max(-1.0, min(1.0, cos_th));
      double th = acos(cos_th);
      double re = fabs(th - target) / target;
      if (re > max_re) max_re = re;
    }
  }
  return max_re;
}

int Deltahedron::count_concave() const {
  int n_concave = 0;
  for (int v = 0; v < N; v++) {
    int deg = degree(v);
    coord3d centroid(0,0,0);
    for (int nb : nbrs(v)) centroid = centroid + points[nb];
    centroid = centroid * (1.0 / deg);
    coord3d fan_normal(0,0,0);
    for (int j = 0; j < deg; j++) {
      coord3d e1 = points[nbrs(v)[j]] - points[v];
      coord3d e2 = points[nbrs(v)[(j+1)%deg]] - points[v];
      fan_normal = fan_normal + e1.cross(e2);
    }
    double nn = fan_normal.norm();
    if (nn < 1e-15) continue;
    coord3d nhat = fan_normal / nn;
    double h = (points[v] - centroid).dot(nhat);
    if (h < 0) n_concave++;
  }
  return n_concave;
}

void Deltahedron::smooth(double q) {
  vector<coord3d> new_points(N);
  for(node_t u = 0; u < N; u++){
    coord3d avg;
    for(node_t v : nbrs(u)) avg += points[v];
    avg /= degree(u);
    new_points[u] = points[u]*(1.0-q) + avg*q;
  }
  points = std::move(new_points);
}

Deltahedron Deltahedron::halma_transform(int m) const {
    // Halma path: direct subdivision via face grids, preserves node IDs
    int n = m + 1;  // = k in GC(k,0) terminology

    vector<map<edge_t,node_t>> face_grids;
    Triangulation T_new = Triangulation::halma_transform(m, &face_grids);

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

Deltahedron Deltahedron::GCtransform(unsigned k, unsigned l) const {
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

    // Scaled Eisenstein corners of this face
    auto [EU,EV] = gcu.arc_coords.at({nu,nv});
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

// Solve 5x5 cyclic tridiagonal system for F-ring.
// All diagonal entries = 6, sub/super-diagonal = -1, corner entries = -1.
// Uses direct Gaussian elimination (5x5 is tiny).
static void solveCyclicTridiag5(vector<coord3d>& rhs) {
    // Matrix: 6 on diagonal, -1 on sub/super, -1 in corners (0,4) and (4,0)
    // Solve by Sherman-Morrison: write A = T + u*v^T where T is standard
    // tridiagonal (with A[0][0]=7, A[4][4]=7) and u*v^T corrects the corners.
    // Or just do direct 5x5 solve — it's fast enough.

    // Direct Gaussian elimination with partial pivoting on 5x5
    double A[5][5] = {
        { 6, -1,  0,  0, -1},
        {-1,  6, -1,  0,  0},
        { 0, -1,  6, -1,  0},
        { 0,  0, -1,  6, -1},
        {-1,  0,  0, -1,  6}
    };
    coord3d b[5];
    for (int i = 0; i < 5; i++) b[i] = rhs[i];

    // Forward elimination
    for (int col = 0; col < 5; col++) {
        // Find pivot
        int pivot = col;
        for (int row = col + 1; row < 5; row++)
            if (fabs(A[row][col]) > fabs(A[pivot][col])) pivot = row;
        if (pivot != col) {
            swap(A[col], A[pivot]);
            swap(b[col], b[pivot]);
        }
        // Eliminate below
        for (int row = col + 1; row < 5; row++) {
            double factor = A[row][col] / A[col][col];
            for (int j = col; j < 5; j++)
                A[row][j] -= factor * A[col][j];
            b[row] = b[row] - b[col] * factor;
        }
    }

    // Back substitution
    for (int i = 4; i >= 0; i--) {
        for (int j = i + 1; j < 5; j++)
            b[i] = b[i] - b[j] * A[i][j];
        b[i] = b[i] / A[i][i];
    }

    for (int i = 0; i < 5; i++) rhs[i] = b[i];
}

// Compute strip vertex coordinates for one expansion step.
// points[] is indexed by ReducibleDual vertex IDs (full graph numbering).
// TODO: Cleanup (cleaner code, factor into higher-level abstractions  + add clearer descriptions of the method).
static void computeStripCoords(const ExtensionStep& step, vector<coord3d>& points) {
    const auto& strip = step.strip;
    const auto& path = step.path;
    const auto& tp = step.tp;
    int n = (int)strip.size();

    if (step.kind.type == ExpKind::L_type) {
        // L-type: quad centroid initial placement.
        for (int j = 0; j < n; j++)
            points[strip[j]] = (points[path[j]] + points[path[j + 1]]
                              + points[tp[j]] + points[tp[j + 1]]) * 0.25;

    } else if (step.kind.type == ExpKind::B_type) {
        // B(0,0): quad centroid for each of the 3 strip vertices.
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

        // 2. Cylinder axis (path center toward tp center)
        coord3d axis = c_tp - c_path;
        double axis_len = axis.norm();
        if (axis_len < 1e-9) {
	  // I'm pretty sure this should never happen. TODO: Verify.
	  fprintf(stderr,"ERROR: Rings coincide. This should never happen.\n");
	  abort();
	  /*
            // Degenerate: rings coincide. Fall back to Laplacian.
            vector<coord3d> rhs(5);
            for (int i = 0; i < 5; i++) {
                int ip1 = (i + 1) % 5;
                rhs[i] = points[path[i]] + points[path[ip1]]
                       + points[tp[i]] + points[tp[ip1]];
            }
            solveCyclicTridiag5(rhs);
            for (int i = 0; i < 5; i++)
                points[strip[i]] = rhs[i];
            return;  // skip projection below
	  */
        }
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
static const neighbours_t C20_seed_neighbours = {
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
    {-1.166968544300538e+00, -8.833420725694727e-01, -7.550764378709098e-01},
    {-4.336347350479605e-01, -1.416065749126212e+00, 7.203946333542167e-01},
    {4.903421527455298e-01, -1.385121171636910e+00, -7.438099299798976e-01},
    {1.545645795507564e-02, 1.407807989188461e-02, -1.646757705195113e+00},
    {-1.479569901277310e+00, -3.599142029230822e-02, 7.223755901917346e-01},
    {-1.202015638761476e+00, 8.478869660586799e-01, -7.406045616538722e-01},
    {4.336352872061797e-01, 1.416065471679303e+00, -7.203950152649671e-01},
    {-4.903422681984767e-01, 1.385120809704024e+00, 7.438101271351236e-01},
    {1.479569795174217e+00, 3.599103708487350e-02, -7.223758222997361e-01},
    {1.166968569761731e+00, 8.833421297215243e-01, 7.550759932463559e-01},
    {1.202015986893285e+00, -8.478857583613777e-01, 7.406053704820970e-01},
    {-1.545716215025896e-02, -1.407832215400697e-02, 1.646757757854968e+00}
};

// C28 dual (Td): 16 vertices, 12 deg-5 + 4 deg-6
static const neighbours_t C28_seed_neighbours = {
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
    {-9.405871929312047e-01, -1.591437541535750e+00, -6.424042055924868e-01},
    {-1.052186933899916e+00, 8.221058471004183e-03, -1.847241739500054e+00},
    {-2.100705949588757e+00, 1.371615510463313e-02, -3.253969728510948e-01},
    {-1.573984731240869e+00, -8.826169230448901e-01, 1.123777687025010e+00},
    {2.049663394302654e-02, -1.816781750144524e+00, 1.103578219245478e+00},
    {4.865624763842553e-01, -8.933243418540859e-01, -1.866649081554954e+00},
    {1.032443532675570e+00, -1.822129193467265e+00, -3.650576037518890e-01},
    {-1.562563067140939e+00, 9.008049711206201e-01, 1.125262888656033e+00},
    {-9.201960031842845e-01, 1.604356014272783e+00, -6.399226210993696e-01},
    {4.978435493587666e-01, 8.901454619450764e-01, -1.865261867604557e+00},
    {4.360413802100851e-02, 1.814635206188322e+00, 1.106524177648666e+00},
    {1.055595166598490e+00, 1.809286918474982e+00, -3.620881547467743e-01},
    {1.837068606500722e+00, -1.122787591788580e-02, -6.747295752041378e-01},
    {1.582280767772851e+00, 9.130387113224273e-01, 1.087163432326433e+00},
    {1.570643194380010e+00, -9.349496985618820e-01, 1.085463797982192e+00},
    {2.369060856185360e-02, -1.729409781154545e-03, 1.956994082833905e+00}
};

// C30 dual (D5h): 17 vertices, 12 deg-5 + 5 deg-6
static const neighbours_t C30_seed_neighbours = {
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
    {-1.559070459063749e+00, -1.068747933944833e+00, -1.879973982189866e+00},
    {-4.682161425554466e-01, -2.046776877036546e+00, -9.492046811648749e-01},
    {1.532096024990294e-01, -9.070467001411456e-01, -2.112729254739929e+00},
    {-8.625543883353577e-01, 5.089098250300979e-01, -2.075384368963578e+00},
    {-1.868119958019254e+00, -1.335183568767094e+00, -1.929375039951397e-01},
    {-2.037636619568821e-01, -1.492201823618911e+00, 1.017290557863881e+00},
    {1.325606589502483e+00, -1.163536525491193e+00, -4.378322309024945e-01},
    {1.023052451276474e+00, 7.730257842827236e-01, -1.287868315997610e+00},
    {-1.451572499786121e+00, 2.413203433573819e-01, 1.066446125811059e+00},
    {-6.933261514500439e-01, 1.641274678779233e+00, -3.581780451280046e-01},
    {-5.135808727388225e-02, 1.654043304436030e+00, 1.603517492377413e+00},
    {-2.111865666607679e+00, 2.443439039860929e-01, -8.890610449743985e-01},
    {1.194817677758050e+00, 1.926052064891783e+00, 4.155148946807128e-01},
    {2.215295968842247e+00, 5.135870518177549e-01, 3.722701229064814e-01},
    {1.559076114414663e+00, 1.068846066695498e+00, 1.880080905318174e+00},
    {1.599833943239003e+00, -6.313650866734757e-01, 1.533561078157907e+00},
    {1.989553219992953e-01, 7.347900093475675e-02, 2.294516632610892e+00}
};

// Get precomputed seed data by type
static const neighbours_t& seedNeighbours(buckinverse::SeedType s) {
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
    const neighbours_t& precomp,
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
    neighbours_t compact_adj(seed_N, GRAPH_DMAX);
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
    const neighbours_t& precomp,
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

    // Try both orientations: the input graph may be either reflection
    // of the planar embedding (Whitney's theorem: unique up to reflection).
    for (int reflect = 0; reflect < 2; reflect++) {
        if (reflect == 1) {
            for (auto& kv : ep_adj)
                reverse(kv.second.begin(), kv.second.end());
        }

        for (int start_p = 0; start_p < seed_N; start_p++) {
            if ((int)precomp[start_p].size() != 5) continue;

            for (const auto& sv : ep.seed_state) {
                if ((int)ep_adj[sv.id].size() != 5) continue;

                // Try all rotations of the neighbor list
                int deg_start = (int)precomp[start_p].size();
                for (int rot = 0; rot < deg_start; rot++) {
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
    }

    return {};  // No match found
}

// Load precomputed seed geometry into the points array, mapped to the
// extension path's vertex IDs.
static void computeSeedGeometry(const ExtensionPath& ep, vector<coord3d>& points) {
    const auto& precomp_nbrs = seedNeighbours(ep.seed);
    const auto& precomp_pts = seedPoints(ep.seed);

    // Try BFS match against precomputed seed (fast path)
    vector<int> mapping = matchSeedViaBFS(precomp_nbrs, ep);

    if (!mapping.empty()) {
        // Normalize to unit sphere (precomputed coords are at physical scale ~1.5 A)
        double scale = 0;
        for (const auto& p : precomp_pts)
            scale = max(scale, p.norm());

        for (int i = 0; i < (int)precomp_pts.size(); i++)
            points[mapping[i]] = precomp_pts[i] / scale;
        return;
    }

    // Fallback: the ep's seed is a different isomer than the precomputed one.
    // Build the seed graph from ep.seed_state and compute geometry from scratch.
    int seed_N = (int)precomp_nbrs.size();
    vector<int> ep_ids;
    ep_ids.reserve(seed_N);
    for (const auto& sv : ep.seed_state)
        ep_ids.push_back(sv.id);
    sort(ep_ids.begin(), ep_ids.end());

    map<int, int> to_compact;
    for (int i = 0; i < seed_N; i++)
        to_compact[ep_ids[i]] = i;

    neighbours_t compact_adj(seed_N, GRAPH_DMAX);
    for (const auto& sv : ep.seed_state) {
        int ci = to_compact[sv.id];
        uint8_t m = sv.active;
        for (int p = 0; p < 6; p++)
            if (m & (1 << p))
                compact_adj.push_back(ci, to_compact[sv.nbr[p]]);
    }

    Graph seed_graph(compact_adj);
    PlanarGraph pg(seed_graph);
    vector<coord3d> seed_pts = pg.zero_order_geometry();

    // Normalize to unit sphere
    double scale = 0;
    for (const auto& p : seed_pts)
        scale = max(scale, p.norm());
    for (auto& p : seed_pts)
        p /= scale;

    // Map compact indices back to ep vertex IDs
    for (int i = 0; i < seed_N; i++)
        points[ep_ids[i]] = seed_pts[i];
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

    neighbours_t adj(id, GRAPH_DMAX);
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

    neighbours_t adj(id, GRAPH_DMAX);
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
    vector<int>& remap)
{
    // 1. Collect free vertices (strip + path + tp)
    set<int> free_set;
    for (int v : step.strip) free_set.insert(v);
    for (int v : step.path)  free_set.insert(v);
    for (int v : step.tp)    free_set.insert(v);

    // 2. Expand to one-ring of (path ∪ tp) for the fixed boundary
    set<int> patch_set = free_set;
    for (int v : step.path) {
        uint8_t m = rd.V[v].active;
        for (; m; m &= m - 1)
            patch_set.insert(rd.V[v].nbr[__builtin_ctz(m)]);
    }
    for (int v : step.tp) {
        uint8_t m = rd.V[v].active;
        for (; m; m &= m - 1)
            patch_set.insert(rd.V[v].nbr[__builtin_ctz(m)]);
    }

    // 3. Remap patch vertices to 0..m-1
    remap.assign(full_N, -1);
    int id = 0;
    for (int v : patch_set) remap[v] = id++;
    int m = id;

    // 4. Build adjacency for the sub-graph (only edges within patch)
    neighbours_t adj(m, GRAPH_DMAX);
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
        int patch_deg = (int)adj[remap[u]].size();
        interior_mask[remap[u]] = (patch_deg == graph_deg);
    }

    // 7. Extract coordinates
    vector<coord3d> patch_pts(m);
    for (int v : patch_set) patch_pts[remap[v]] = points[v];

    return Deltahedron(Triangulation(adj), patch_pts);
}

Deltahedron Deltahedron::fromExtensionPathOptimized(const ExtensionPath& ep, int max_iter_per_step, FILE* log, StepCallback diag,
                                                     OptMethod method, double step_tol, double final_tol, long long max_work_per_step,
                                                     double step_angle_tol, double final_angle_tol) {
    int full_N = ep.full_N;
    vector<coord3d> points(full_N);

    // Timing accumulators (only used when opt_log is set)
    using clk = chrono::steady_clock;
    double ms_seed = 0, ms_place = 0, ms_reflect = 0, ms_patch = 0, ms_cg = 0, ms_final = 0;
    int total_cg_iters = 0, total_patch_iters = 0;
    int acc_energy = 0, acc_grad = 0, acc_hv = 0;
    auto t0 = clk::now();

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

        // a. Place strip vertices (quad centroid initial guess)
        computeStripCoords(step, points);

        // b. Expand topology + enforce outward convexity
        rd.expand(step);
        liftStripToSurface(step, rd, points);

        auto t_place = clk::now();
        ms_place += chrono::duration<double,milli>(t_place - ts).count();

        // Diagnostic: after strip placement + lift
        if (diag) {
            vector<int> diag_remap;
            Deltahedron D_diag = extractCompact(rd, full_N, points, diag_remap);
            diag(step_idx + 1, "placed", D_diag);
        }

        // c. Reflect concave free vertices on the full graph (before patch extraction).
        //    This is O(N) and saves the optimizer from spending expensive CG iterations
        //    fighting concavity from inverted strip vertices.
        {
            set<int> free_set;
            for (int v : step.strip) free_set.insert(v);
            for (int v : step.path)  free_set.insert(v);
            for (int v : step.tp)    free_set.insert(v);

            vector<int> refl_remap;
            Deltahedron D = extractCompact(rd, full_N, points, refl_remap);

            // Build fixed mask: only free_set vertices may be reflected
            vector<bool> refl_fixed(D.N, true);
            for (int u : free_set)
                if (refl_remap[u] >= 0) refl_fixed[refl_remap[u]] = false;

            for (int pass = 0; pass < 3; pass++)
                if (D.reflect_concave(D.points, 0, refl_fixed) == 0) break;

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

        // d. Extract small patch sub-graph (O(1) vertices)
        vector<bool> free_mask, interior_mask;
        vector<int> patch_remap;
        Deltahedron patch = extractPatch(rd, full_N, step, points,
                                         free_mask, interior_mask, patch_remap);

        // e. Trust-region optimize on patch only (strip + path + tp free)
        patch.opt_log = log;
        patch.optimize_patch(patch.points, free_mask, interior_mask);
        patch.opt_log = nullptr;
        total_patch_iters += patch.iterations_used;

        // f. Copy patch free-vertex coords back to full array
        for (int u = 0; u < full_N; u++)
            if (patch_remap[u] >= 0 && free_mask[patch_remap[u]])
                points[u] = patch.points[patch_remap[u]];

        auto t_patch = clk::now();
        ms_patch += chrono::duration<double,milli>(t_patch - t_refl).count();

        // Diagnostic: after patch optimize, BEFORE full-graph CG (the key snapshot)
        if (diag) {
            vector<int> diag_remap;
            Deltahedron D_diag = extractCompact(rd, full_N, points, diag_remap);
            diag(step_idx + 1, "patched", D_diag);
        }

        // g. Full-graph CG relaxation to release accumulated strain.
        vector<int> full_remap;
        Deltahedron D = extractCompact(rd, full_N, points, full_remap);

        int iter_budget = max_iter_per_step > 0 ? max_iter_per_step : 3*D.N;
        D.opt_k_flat = 0;  // skip E_flat in intermediate steps (final CG uses default)
        D.opt_method = method;
        // Log CG details for a few representative steps
        int n_steps = (int)ep.steps.size();
        bool log_this = log && (step_idx <= 2 || step_idx == n_steps/2 || step_idx >= n_steps - 2);
        if (log_this) D.opt_log = log;
        D.optimize(D.points, 0, iter_budget, step_tol, {}, max_work_per_step, step_angle_tol);
        if (log_this) D.opt_log = nullptr;
        total_cg_iters += D.iterations_used;
        acc_energy += D.n_energy_evals;
        acc_grad += D.n_grad_evals;
        acc_hv += D.n_hv_evals;

        // h. Write optimized coordinates back to full array
        for (int u = 0; u < full_N; u++)
            if (full_remap[u] >= 0)
                points[u] = D.points[full_remap[u]];

        auto t_cg = clk::now();
        ms_cg += chrono::duration<double,milli>(t_cg - t_patch).count();

        // Diagnostic: after CG
        if (diag) {
            vector<int> diag_remap;
            Deltahedron D_diag = extractCompact(rd, full_N, points, diag_remap);
            diag(step_idx + 1, "cg", D_diag);
        }

        if (log) {
            fprintf(log, "  step %2d: N=%3d place=%.1f refl=%.1f patch=%.1f(%d) cg=%.1f(%d) ms\n",
                    step_idx, D.N,
                    chrono::duration<double,milli>(t_place - ts).count(),
                    chrono::duration<double,milli>(t_refl - t_place).count(),
                    chrono::duration<double,milli>(t_patch - t_refl).count(), patch.iterations_used,
                    chrono::duration<double,milli>(t_cg - t_patch).count(), D.iterations_used);
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
    D.opt_method = method;
    D.optimize(D.points, 0, 12*D.N, final_tol, {}, max_work_per_step, final_angle_tol);
    ms_final = chrono::duration<double,milli>(clk::now() - t_final).count();

    // Diagnostic: after final optimization
    if (diag) diag((int)ep.steps.size() + 1, "final", D);

    // Accumulate eval counters from all per-step + final optimize calls
    acc_energy += D.n_energy_evals;
    acc_grad += D.n_grad_evals;
    acc_hv += D.n_hv_evals;
    D.n_energy_evals = acc_energy;
    D.n_grad_evals = acc_grad;
    D.n_hv_evals = acc_hv;

    if (log) {
        double ms_total = chrono::duration<double,milli>(clk::now() - t0).count();
        fprintf(log, "  totals: seed=%.0f place=%.0f refl=%.0f patch=%.0f(%d) cg=%.0f(%d) final=%.0f(%d) total=%.0f ms\n",
                ms_seed, ms_place, ms_reflect, ms_patch, total_patch_iters,
                ms_cg, total_cg_iters, ms_final, D.iterations_used, ms_total);
    }

    D.iterations_used += total_cg_iters;
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
    const Deltahedron& D,
    const vector<edge_t>& edges,
    std::span<const coord3d> x,
    vector<coord3d>* grad,
    double L,           // target edge length
    double k_bond,
    double k_angle,
    double k_curv,
    double k_flat,
    double k_conv,
    const vector<bool>& conv_mask = {})  // if non-empty, restrict E_conv to these vertices
{
  const int N = D.N;
  double energy = 0.0;

  // Zero gradient
  if(grad)
    for(int i = 0; i < N; i++) (*grad)[i] = coord3d(0,0,0);

  // --- E_bond: harmonic springs on edge lengths ---
  // E_bond = (k_bond/2) * Sum_edges (|x_u - x_v| - L)^2
  for(const edge_t& e : edges){
    int u = e.first, v = e.second;
    coord3d diff = x[u] - x[v];
    double r = diff.norm();
    if(r < 1e-15) continue;  // degenerate edge, skip
    double dev = r - L;
    energy += 0.5 * k_bond * dev * dev;

    if(grad){
      coord3d g = diff * (k_bond * dev / r);
      (*grad)[u] += g;
      (*grad)[v] -= g;
    }
  }

  // --- E_angle: harmonic springs on triangle corner angles ---
  // E_angle = (k_angle/2) * Sum_triangles Sum_corners (angle - pi/3)^2
  const double target_angle = M_PI / 3.0;
  for(const auto& tri : D.triangles()){
    int v[3] = {tri[0], tri[1], tri[2]};

    for(int c = 0; c < 3; c++){
      // Angle at vertex v[c] between v[c]->v[c-1] and v[c]->v[c+1]
      int b = v[c];
      int a = v[(c+2) % 3];  // previous
      int d = v[(c+1) % 3];  // next

      coord3d va = x[a] - x[b];
      coord3d vc = x[d] - x[b];

      double theta = coord3d::angle(va, vc);
      double dev = theta - target_angle;
      energy += 0.5 * k_angle * dev * dev;

      if(grad){
        coord3d da, dc;
        coord3d::dangle(va, vc, da, dc);
        double w = k_angle * dev;
        (*grad)[a] += da * w;
        (*grad)[d] += dc * w;
        (*grad)[b] -= (da + dc) * w;
      }
    }
  }

  // --- E_curv: harmonic springs on discrete Gaussian curvature ---
  // K(v) = 2pi - sum of face angles at v
  // K_target(v) = 2pi - deg(v)*pi/3
  // E_curv = (k_curv/2) * Sum_v (K(v) - K_target(v))^2
  //        = (k_curv/2) * Sum_v (deg(v)*pi/3 - angle_sum(v))^2
  for(int v = 0; v < N; v++){
    int d = D.degree(v);
    double angle_sum = 0.0;

    // Compute angle sum (and store derivatives if needed)
    vector<coord3d> da_list, dc_list;
    if(grad){ da_list.resize(d); dc_list.resize(d); }

    for(int i = 0; i < d; i++){
      int ni   = D.nbrs(v)[i];
      int ni1  = D.nbrs(v)[(i+1) % d];
      coord3d va = x[ni]  - x[v];
      coord3d vc = x[ni1] - x[v];

      angle_sum += coord3d::angle(va, vc);
      if(grad) coord3d::dangle(va, vc, da_list[i], dc_list[i]);
    }

    double target_sum = d * M_PI / 3.0;
    double dev = target_sum - angle_sum;  // = K(v) - K_target(v)
    energy += 0.5 * k_curv * dev * dev;

    if(grad){
      // Gradient: dE/d(...) = k_curv * dev * d(dev)/d(...)
      // dev = target_sum - angle_sum, so d(dev)/d(...) = -d(angle_sum)/d(...)
      double w = -k_curv * dev;  // negative because dev = target - sum
      for(int i = 0; i < d; i++){
        int ni  = D.nbrs(v)[i];
        int ni1 = D.nbrs(v)[(i+1) % d];

        (*grad)[ni]  += da_list[i] * w;
        (*grad)[ni1] += dc_list[i] * w;
        (*grad)[v]   -= (da_list[i] + dc_list[i]) * w;
      }
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
        int ni  = D.nbrs(v)[i];
        int ni1 = D.nbrs(v)[(i+1) % d];
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
          (*grad)[D.nbrs(v)[j]] += (fj + fjpre) * (scale / 3.0);
        }
      }
    }
  }

  // E_conv: smooth convexity bias via softplus (currently disabled in optimize(),
  // replaced by periodic vertex reflection; retained for gradient_check() and
  // potential future use).
  // For each qualifying vertex, compute signed height h above neighbor centroid
  // plane (positive = convex). Penalty: k_conv * sigma * log(1 + exp(-h/sigma)).
  // Nearly zero for h > 0, linear in -h for h < 0, smooth everywhere.
  // Exact gradient including normal-rotation term for each neighbour.
  if(k_conv > 0){
    const double sigma = 0.2 * L;  // transition width ~ 20% of edge length
    for(int v = 0; v < N; v++){
      int d = D.degree(v);
      if(d > 6) continue;
      if(!conv_mask.empty() && !conv_mask[v]) continue;  // skip truncated boundary verts

      // Neighbor centroid
      coord3d centroid(0,0,0);
      for(int i = 0; i < d; i++) centroid += x[D.nbrs(v)[i]];
      centroid /= (double)d;

      // Fan normal (unnormalized, outward for convex)
      coord3d N_fan(0,0,0);
      for(int i = 0; i < d; i++){
        coord3d e1 = x[D.nbrs(v)[i]] - x[v];
        coord3d e2 = x[D.nbrs(v)[(i+1)%d]] - x[v];
        N_fan += e1.cross(e2);
      }
      double N_len = N_fan.norm();
      if(N_len < 1e-15) continue;
      coord3d n_hat = N_fan / N_len;

      double h = (x[v] - centroid).dot(n_hat);

      // Softplus: sigma * log(1 + exp(-h/sigma))
      double z = h / sigma;
      double sp, sig;
      if(z > 20){         // convex: softplus ≈ 0
        sp = sigma * exp(-z);
        sig = 0;
      } else if(z < -20){ // concave: softplus ≈ -h
        sp = -h;
        sig = 1;
      } else {
        sp = sigma * log(1 + exp(-z));
        sig = 1.0 / (1 + exp(z));  // sigmoid(-z)
      }

      energy += k_conv * sp;

      if(grad){
        // dE/dh = -k_conv * sigmoid(-h/sigma)
        double dEdh = -k_conv * sig;

        // Exact gradient: dh/dx[v] = n_hat (same as frozen-normal,
        // since dN/dx[v] = 0 by telescoping sum of edge differences)
        (*grad)[v] += n_hat * dEdh;

        // Exact gradient for neighbours:
        // dh/dx[n_j] = -n_hat/d + r_perp x (e_{j-1} - e_{j+1}) / |N|
        // where r_perp = (x[v] - centroid) - h*n_hat, e_k = x[n_k] - x[v]
        coord3d r_perp = (x[v] - centroid) - n_hat * h;
        for(int j = 0; j < d; j++){
          coord3d ej_prev = x[D.nbrs(v)[(j+d-1)%d]] - x[v];
          coord3d ej_next = x[D.nbrs(v)[(j+1)%d]]   - x[v];
          coord3d dhdx_nj = n_hat * (-1.0/d)
                          + r_perp.cross(ej_prev - ej_next) / N_len;
          (*grad)[D.nbrs(v)[j]] += dhdx_nj * dEdh;
        }
      }
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
    const Deltahedron& D,
    const vector<edge_t>& edges,
    std::span<const coord3d> x,
    const vector<coord3d>& v,
    vector<coord3d>& Hv,
    double L,
    double k_bond, double k_angle, double k_curv, double k_flat,
    const vector<bool>& fixed = {})
{
  const int N = D.N;
  const matrix3d I = matrix3d::unit_matrix();
  const double theta0 = M_PI / 3.0;

  // Helper: apply 3x3 block M to vertex pair (i,j), accumulating M*v[j] into Hv[i]
  auto apply_block = [&](int vi, int vj, const matrix3d& M){
    Hv[vi] = Hv[vi] + M * v[vj];
  };

  // --- E_bond Hv ---
  // Per edge: M = k*[(1-L/r)*I + (L/r^3)*d⊗d]
  // Hv[u] += M*(v[u]-v[v]),  Hv[v] -= M*(v[u]-v[v])
  if(k_bond > 0){
    for(const edge_t& e : edges){
      int u = e.first, vv = e.second;
      coord3d d = x[u] - x[vv];
      double r = d.norm();
      if(r < 1e-15) continue;

      coord3d dv = v[u] - v[vv];
      // M*dv = k*(1 - L/r)*dv + k*(L/r^3)*(d.dv)*d
      coord3d Mdv = dv * (k_bond * (1 - L/r)) + d * (k_bond * L / (r*r*r) * d.dot(dv));
      Hv[u] = Hv[u] + Mdv;
      Hv[vv] = Hv[vv] - Mdv;
    }
  }

  // --- E_angle Hv ---
  // Per triangle corner: compute arm-space blocks Haa, Hcc, Hac and scatter.
  if(k_angle > 0){
    for(const auto& tri : D.triangles()){
      for(int c = 0; c < 3; c++){
        int b = tri[c], a = tri[(c+2)%3], dd = tri[(c+1)%3];
        coord3d va = x[a] - x[b], vc = x[dd] - x[b];
        double ra = va.norm(), rc = vc.norm();
        if(ra < 1e-15 || rc < 1e-15) continue;

        coord3d ua = va / ra, uc = vc / rc;
        double C = max(-1.0, min(1.0, ua.dot(uc)));
        double theta = acos(C);
        double S = sin(theta);
        if(S < 1e-10) continue;

        double dev = theta - theta0;
        double alpha = 1.0 - dev * C / S;

        coord3d p, q;
        coord3d::dangle(va, vc, p, q);

        matrix3d sym_ac = ua.outer(uc) + uc.outer(ua);

        matrix3d Haa = p.outer(p) * (k_angle * alpha)
          + (sym_ac + I * C - ua.outer(ua) * (3*C)) * (k_angle * dev / (ra*ra * S));

        matrix3d Hcc = q.outer(q) * (k_angle * alpha)
          + (sym_ac + I * C - uc.outer(uc) * (3*C)) * (k_angle * dev / (rc*rc * S));

        matrix3d Hac = p.outer(q) * (k_angle * alpha)
          + (ua.outer(ua) + uc.outer(uc) - ua.outer(uc) * C - I) * (k_angle * dev / (ra*rc * S));

        // Scatter: same pattern as assemble_patch_hessian, but M*v instead of storing M
        matrix3d Hab = Haa + Hac;
        matrix3d Hcb = Hcc + Hac.transpose();

        // H(a,a)*v[a] + H(a,dd)*v[dd] + H(a,b)*v[b]
        Hv[a] = Hv[a] + Haa * v[a] + Hac * v[dd] + (-Hab) * v[b];
        // H(dd,a)*v[a] + H(dd,dd)*v[dd] + H(dd,b)*v[b]
        Hv[dd] = Hv[dd] + Hac.transpose() * v[a] + Hcc * v[dd] + (-Hcb) * v[b];
        // H(b,a)*v[a] + H(b,dd)*v[dd] + H(b,b)*v[b]
        Hv[b] = Hv[b] + (-Hab.transpose()) * v[a] + (-Hcb.transpose()) * v[dd] + (Hab + Hcb) * v[b];
      }
    }
  }

  // --- E_curv Hv ---
  // E_curv = (k_curv/2) * sum_v (K_target - angle_sum)^2
  // H_curv = k_curv * sum_v [dK⊗dK + dev * d^2K/dx^2]
  //
  // For the rank-1 term: dK_v/dx is the Gaussian curvature gradient at vertex v.
  // dK_v/dx[ni] = -da_i - dc_{i-1}  (sum of angle derivatives at v)
  // dK_v/dx[v]  = sum_i (da_i + dc_i)
  //
  // For the d^2K/dx^2 term: same arm-space angle Hessian blocks as E_angle,
  // but summed per-vertex (over the fan of angles at v) instead of per-triangle.
  if(k_curv > 0){
    for(int vertex = 0; vertex < N; vertex++){
      int deg = D.degree(vertex);

      // First pass: compute angle sum, derivatives, and curvature deviation
      double angle_sum = 0;
      vector<coord3d> da_list(deg), dc_list(deg);
      for(int i = 0; i < deg; i++){
        int ni  = D.nbrs(vertex)[i];
        int ni1 = D.nbrs(vertex)[(i+1) % deg];
        coord3d va = x[ni] - x[vertex], vc = x[ni1] - x[vertex];
        angle_sum += coord3d::angle(va, vc);
        coord3d::dangle(va, vc, da_list[i], dc_list[i]);
      }
      double dev = deg * M_PI / 3.0 - angle_sum;

      // Rank-1 term: Hv += k_curv * (dK . v) * dK
      // dK/dx[ni] = -sum of da_i terms pointing to ni - sum of dc_j terms pointing to ni
      // dK/dx[v] = sum_i (da_i + dc_i)  (the chain rule d(va)/d(x_b) = -I part)
      //
      // But dK.v = -dev_weight * sum_i (da_i . (v[ni] - v[vertex]) + dc_i . (v[ni1] - v[vertex]))
      // where dev_weight = ... actually let's compute it directly.
      //
      // The gradient of E_curv at vertex `vertex` contributes:
      //   g[ni]  += w * da_i   (for each i where ni = neighbours[vertex][i])
      //   g[ni1] += w * dc_i
      //   g[vertex] -= w * (da_i + dc_i)
      // where w = -k_curv * dev.
      //
      // So the Jacobian dK/dx (without the -k_curv*dev factor) has:
      //   (dK/dx)[ni]     += -da_i       (and also -dc_{i-1})
      //   (dK/dx)[vertex] += (da_i + dc_i)
      //
      // For the rank-1 Hv: factor = k_curv * (dK . v), then scatter factor * dK.
      // dK = d(angle_sum)/dx  (gradient of angle sum, not energy gradient)

      // Compute dK . v = sum_i [ da_i . v[ni] + dc_i . v[ni1] - (da_i + dc_i) . v[vertex] ]
      double dK_dot_v = 0;
      for(int i = 0; i < deg; i++){
        int ni  = D.nbrs(vertex)[i];
        int ni1 = D.nbrs(vertex)[(i+1) % deg];
        dK_dot_v += da_list[i].dot(v[ni] - v[vertex]) + dc_list[i].dot(v[ni1] - v[vertex]);
      }

      // Scatter: Hv += k_curv * dK_dot_v * dK
      double w = k_curv * dK_dot_v;
      for(int i = 0; i < deg; i++){
        int ni  = D.nbrs(vertex)[i];
        int ni1 = D.nbrs(vertex)[(i+1) % deg];
        Hv[ni]     = Hv[ni]     + da_list[i] * w;
        Hv[ni1]    = Hv[ni1]    + dc_list[i] * w;
        Hv[vertex] = Hv[vertex] - (da_list[i] + dc_list[i]) * w;
      }

      // Curvature correction term: k_curv * dev * d^2(angle_sum)/dx^2 * v
      // This is the same angle Hessian blocks as E_angle, but summed over the
      // fan of angles at vertex `vertex`, weighted by k_curv * dev.
      if(fabs(dev) > 1e-15){
        double w2 = -k_curv * dev;  // negative because dev = target - sum, d(dev)/d(angle) = -1

        for(int i = 0; i < deg; i++){
          int ni  = D.nbrs(vertex)[i];
          int ni1 = D.nbrs(vertex)[(i+1) % deg];
          int b = vertex;

          coord3d va = x[ni] - x[b], vc = x[ni1] - x[b];
          double ra = va.norm(), rc = vc.norm();
          if(ra < 1e-15 || rc < 1e-15) continue;

          coord3d ua = va / ra, uc = vc / rc;
          double Cos = max(-1.0, min(1.0, ua.dot(uc)));
          double theta = acos(Cos);
          double Sin = sin(theta);
          if(Sin < 1e-10) continue;

          double angle_dev = theta - theta0;
          double alph = 1.0 - angle_dev * Cos / Sin;

          coord3d p_a, q_c;
          coord3d::dangle(va, vc, p_a, q_c);

          matrix3d sym = ua.outer(uc) + uc.outer(ua);

          matrix3d Haa = p_a.outer(p_a) * alph
            + (sym + I * Cos - ua.outer(ua) * (3*Cos)) * (angle_dev / (ra*ra * Sin));

          matrix3d Hcc = q_c.outer(q_c) * alph
            + (sym + I * Cos - uc.outer(uc) * (3*Cos)) * (angle_dev / (rc*rc * Sin));

          matrix3d Hac_m = p_a.outer(q_c) * alph
            + (ua.outer(ua) + uc.outer(uc) - ua.outer(uc) * Cos - I) * (angle_dev / (ra*rc * Sin));

          // Scale by w2 and scatter (same pattern)
          matrix3d Hab_m = Haa + Hac_m;
          matrix3d Hcb_m = Hcc + Hac_m.transpose();

          Hv[ni]  = Hv[ni]  + (Haa * v[ni] + Hac_m * v[ni1] + (-Hab_m) * v[b]) * w2;
          Hv[ni1] = Hv[ni1] + (Hac_m.transpose() * v[ni] + Hcc * v[ni1] + (-Hcb_m) * v[b]) * w2;
          Hv[b]   = Hv[b]   + ((-Hab_m.transpose()) * v[ni] + (-Hcb_m.transpose()) * v[ni1] + (Hab_m + Hcb_m) * v[b]) * w2;
        }
      }
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
    const Deltahedron& D,
    const vector<edge_t>& edges,
    std::span<const coord3d> x,
    double L, double k_bond, double k_angle, double k_curv, double k_flat, double k_conv,
    const vector<bool>& conv_mask = {})
{
  return deltahedron_energy_and_gradient(D, edges, x, nullptr, L, k_bond, k_angle, k_curv, k_flat, k_conv, conv_mask);
}

// Check convexity constraint: h(v) >= -tau*L for all free vertices.
// Returns true if all constraints satisfied.
static bool check_convexity(const Deltahedron& D, const vector<coord3d>& x,
                            const vector<bool>& free_mask, double L, double tau = 0.05)
{
  for(int v = 0; v < D.N; v++){
    if(!free_mask[v]) continue;
    int d = D.degree(v);

    coord3d centroid(0,0,0);
    for(int j = 0; j < d; j++) centroid += x[D.nbrs(v)[j]];
    centroid /= (double)d;

    coord3d n_fan(0,0,0);
    for(int j = 0; j < d; j++){
      coord3d e1 = x[D.nbrs(v)[j]] - x[v];
      coord3d e2 = x[D.nbrs(v)[(j+1)%d]] - x[v];
      n_fan += e1.cross(e2);
    }
    double n_len = n_fan.norm();
    if(n_len < 1e-15) continue;

    double h = (x[v] - centroid).dot(n_fan / n_len);
    if(h < -tau * L) return false;
  }
  return true;
}

// Assemble exact analytical Hessian for the patch optimizer.
// E_bond: exact.  E_angle: exact.  E_conv: exact (including d²h/dx² correction).
// H must be ndof x ndof (ndof = 3*nfree), zero-initialized by caller.
static void assemble_patch_hessian(
    const Deltahedron& D,
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

  // --- E_bond = (k/2) Sum_edges (|x_u - x_v| - L)^2 ---
  // Exact Hessian per edge: M = k*[(1-L/r)*I + (L/r^3)*d⊗d]
  // H(u,u) = H(v,v) = +M,  H(u,v) = H(v,u) = -M
  if(k_bond > 0){
    for(const edge_t& e : edges){
      int u = e.first, v = e.second;
      coord3d d = x[u] - x[v];
      double r = d.norm();
      if(r < 1e-15) continue;

      matrix3d M = I * (k_bond * (1 - L/r)) + d.outer(d) * (k_bond * L / (r*r*r));

      add_block(u, u, M);
      add_block(v, v, M);
      add_block(u, v, -M);
      add_block(v, u, -M);
    }
  }

  // --- E_angle = (k/2) Sum_corners (theta - pi/3)^2 ---
  // Exact Hessian.  Angle theta at vertex b, arms va = x_a - x_b, vc = x_c - x_b.
  // Unit vectors ua = va/ra, uc = vc/rc.  C = cos(theta), S = sin(theta).
  // First derivatives: p = dtheta/d(va), q = dtheta/d(vc) (from coord3d::dangle).
  //
  // Arm-space Hessian (second derivatives w.r.t. va, vc):
  //   alpha = 1 - dev*cot(theta)
  //   Haa = k*alpha * p⊗p + (k*dev)/(ra^2*S) * [ua⊗uc + uc⊗ua + C*I - 3C*ua⊗ua]
  //   Hcc = k*alpha * q⊗q + (k*dev)/(rc^2*S) * [uc⊗ua + ua⊗uc + C*I - 3C*uc⊗uc]
  //   Hac = k*alpha * p⊗q + (k*dev)/(ra*rc*S) * [ua⊗ua + uc⊗uc - C*ua⊗uc - I]
  //
  // Chain rule to vertex coordinates (d(va)/d(x_b) = -I, d(vc)/d(x_b) = -I):
  //   H(a,b) = -(Haa + Hac),  H(c,b) = -(Hcc + Hca),  H(b,b) = sum of all four
  if(k_angle > 0){
    const double theta0 = M_PI / 3.0;
    for(const auto& tri : D.triangles()){
      for(int c = 0; c < 3; c++){
        int b = tri[c], a = tri[(c+2)%3], dd = tri[(c+1)%3];
        coord3d va = x[a] - x[b], vc = x[dd] - x[b];
        double ra = va.norm(), rc = vc.norm();
        if(ra < 1e-15 || rc < 1e-15) continue;

        coord3d ua = va / ra, uc = vc / rc;
        double C = max(-1.0, min(1.0, ua.dot(uc)));
        double theta = acos(C);
        double S = sin(theta);
        if(S < 1e-10) continue;

        double dev = theta - theta0;
        double alpha = 1.0 - dev * C / S;

        coord3d p, q;
        coord3d::dangle(va, vc, p, q);

        // Arm-space Hessian blocks (rank-1 + correction)
        matrix3d sym_ac = ua.outer(uc) + uc.outer(ua);  // symmetric part, shared by Haa and Hcc

        matrix3d Haa = p.outer(p) * (k_angle * alpha)
          + (sym_ac + I * C - ua.outer(ua) * (3*C)) * (k_angle * dev / (ra*ra * S));

        matrix3d Hcc = q.outer(q) * (k_angle * alpha)
          + (sym_ac + I * C - uc.outer(uc) * (3*C)) * (k_angle * dev / (rc*rc * S));

        matrix3d Hac = p.outer(q) * (k_angle * alpha)
          + (ua.outer(ua) + uc.outer(uc) - ua.outer(uc) * C - I) * (k_angle * dev / (ra*rc * S));

        // Scatter to vertex coordinates via chain rule
        add_block(a, a, Haa);
        add_block(a, dd, Hac);
        add_block(dd, a, Hac.transpose());
        add_block(dd, dd, Hcc);

        matrix3d Hab = Haa + Hac;               // H(a,b) = -Hab
        matrix3d Hcb = Hcc + Hac.transpose();   // H(c,b) = -Hcb
        add_block(a, b, -Hab);
        add_block(b, a, -Hab.transpose());
        add_block(dd, b, -Hcb);
        add_block(b, dd, -Hcb.transpose());
        add_block(b, b, Hab + Hcb);
      }
    }
  }

  // --- E_conv = k * Sum_v softplus(-h/sigma) ---
  // Exact Hessian: d²E/dh² * (dh/dx)⊗(dh/dx) + dE/dh * d²h/dx²
  // where d²E/dh² = (k/sigma²) * sig * (1 - sig),  dE/dh = -(k/sigma) * sig.
  //
  // h = (x_v - centroid) · n_hat,  n_hat = N/|N|,  N = Σ (e_i × e_{i+1}).
  // Key fact: dN/dx_v = 0 (telescoping), so dn_hat/dx_v = 0 and d²h/dx_v² = 0.
  //
  // Hessian blocks of h:
  //   H_{vv} = 0
  //   H_{v,j} = P · [Δe_j]× / |N|                     (dn_hat/dx_{n_j})
  //   H_{j,k} = -(1/d)Qk + {[Δe_j]×/d + (Δe_j × n_hat)⊗g_k + h·[Δe_j]×·Qk
  //             + δ_{k,j-1}[r_perp]× - δ_{k,j+1}[r_perp]×} / |N|
  //             - w_j ⊗ (n_hat × Δe_k) / |N|²
  // where Qk = P·[Δe_k]×/|N|, Δe_j = e_{j-1} - e_{j+1}, [·]× = cross_matrix.
  if(k_conv > 0){
    using Mx = matrix3d;  // shorthand

    for(int v = 0; v < D.N; v++){
      int d = D.degree(v);
      if(d > 6) continue;
      if(!conv_mask.empty() && !conv_mask[v]) continue;

      coord3d centroid(0,0,0);
      for(int i = 0; i < d; i++) centroid += x[D.nbrs(v)[i]];
      centroid /= (double)d;

      coord3d N_fan(0,0,0);
      for(int i = 0; i < d; i++){
        coord3d e1 = x[D.nbrs(v)[i]] - x[v];
        coord3d e2 = x[D.nbrs(v)[(i+1)%d]] - x[v];
        N_fan += e1.cross(e2);
      }
      double N_len = N_fan.norm();
      if(N_len < 1e-15) continue;
      coord3d n_hat = N_fan / N_len;

      double h = (x[v] - centroid).dot(n_hat);
      double z = h / sigma;
      double sig;
      if(z > 20) sig = 0;
      else if(z < -20) sig = 1;
      else sig = 1.0 / (1 + exp(z));

      double d2Edh2 = (k_conv / sigma) * sig * (1 - sig);
      double dEdh   = -k_conv * sig;
      if(abs(d2Edh2) < 1e-15 && abs(dEdh) < 1e-15) continue;

      coord3d r_perp = (x[v] - centroid) - n_hat * h;
      const Mx P = I - n_hat.outer(n_hat);
      const Mx R = Mx::cross_matrix(r_perp);  // [r_perp]×

      // Precompute per-neighbor data
      struct NbrData {
        int id;          // global vertex index
        coord3d De;      // Δe_j = e_{j-1} - e_{j+1}
        coord3d g;       // dh/dx_{n_j}
        coord3d w;       // r_perp × Δe_j / |N|
        Mx Dex;          // [Δe_j]×
        Mx Q;            // P · [Δe_j]× / |N|
      };
      vector<NbrData> nb(d);
      for(int j = 0; j < d; j++){
        nb[j].id = D.nbrs(v)[j];
        coord3d ej_prev = x[D.nbrs(v)[(j+d-1)%d]] - x[v];
        coord3d ej_next = x[D.nbrs(v)[(j+1)%d]]   - x[v];
        nb[j].De  = ej_prev - ej_next;
        nb[j].g   = n_hat * (-1.0/d) + r_perp.cross(nb[j].De) / N_len;
        nb[j].w   = r_perp.cross(nb[j].De) / N_len;
        nb[j].Dex = Mx::cross_matrix(nb[j].De);
        nb[j].Q   = P * nb[j].Dex * (1.0/N_len);
      }

      // --- Rank-1 term: d²E/dh² * (dh/dx)⊗(dh/dx) ---
      if(abs(d2Edh2) > 1e-15){
        // v-v block (g_v = n_hat)
        add_block(v, v, n_hat.outer(n_hat) * d2Edh2);
        // v-j and j-v blocks
        for(int j = 0; j < d; j++){
          Mx Bvj = n_hat.outer(nb[j].g) * d2Edh2;
          add_block(v, nb[j].id, Bvj);
          add_block(nb[j].id, v, Bvj.transpose());
        }
        // j-k blocks
        for(int j = 0; j < d; j++)
          for(int k = 0; k < d; k++)
            add_block(nb[j].id, nb[k].id, nb[j].g.outer(nb[k].g) * d2Edh2);
      }

      // --- Correction term: dE/dh * d²h/dx² ---
      if(abs(dEdh) > 1e-15){
        // H_{vv} = 0  (nothing to add)

        // H_{v,j} = P · [Δe_j]× / |N|
        for(int j = 0; j < d; j++){
          Mx Hvj = nb[j].Q * dEdh;
          add_block(v, nb[j].id, Hvj);
          add_block(nb[j].id, v, Hvj.transpose());
        }

        // H_{j,k} blocks
        for(int j = 0; j < d; j++){
          for(int k = 0; k < d; k++){
            // Term A: -(1/d) Qk
            Mx Hjk = nb[k].Q * (-1.0/d);
            // Term B: [Δe_j]× / (d·|N|)
            Hjk += nb[j].Dex * (1.0 / (d * N_len));
            // Term C: (Δe_j × n_hat) ⊗ g_k / |N|
            Hjk += nb[j].De.cross(n_hat).outer(nb[k].g) * (1.0/N_len);
            // Term D: h · [Δe_j]× · Qk / |N|
            Hjk += nb[j].Dex * nb[k].Q * (h / N_len);
            // Term E: δ_{k,j-1} [r_perp]× / |N|
            if(k == (j+d-1)%d) Hjk += R * (1.0/N_len);
            // Term F: -δ_{k,j+1} [r_perp]× / |N|
            if(k == (j+1)%d)   Hjk += R * (-1.0/N_len);
            // Term G: -w_j ⊗ (n_hat × Δe_k) / |N|
            // (from product rule: (r_perp × Δe_j) · d(1/|N|)/dx_{n_k}, and w_j = (r_perp × Δe_j)/|N|)
            Hjk += nb[j].w.outer(n_hat.cross(nb[k].De)) * (-1.0/N_len);

            add_block(nb[j].id, nb[k].id, Hjk * dEdh);
          }
        }
      }
    }
  }
}

bool Deltahedron::optimize_patch(std::span<const coord3d> initial_geometry,
                                 const vector<bool>& free_mask,
                                 const vector<bool>& interior_mask,
                                 double target_L, int max_iter, double grad_tol)
{
  assert((int)initial_geometry.size() == N);
  assert((int)free_mask.size() == N);
  points = vector<coord3d>(initial_geometry.begin(), initial_geometry.end());

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
  const double conv_tau = 0.05;  // hard convexity threshold: h >= -tau*L

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
    return deltahedron_energy_only(*this, edges, x, L, k_bond, k_angle, k_curv, k_flat, k_conv, interior_mask);
  };
  auto grad_fn = [&](std::span<const coord3d> x, vector<coord3d>& g) -> double {
    return deltahedron_energy_and_gradient(*this, edges, x, &g, L, k_bond, k_angle, k_curv, k_flat, k_conv, interior_mask);
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
  double E_prev = 1e30;
  int stall_count = 0;

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

    // Trust region stall detection: exit if Delta has shrunk to near-zero
    if(Delta < 1e-12 * L){
      stall_count++;
      if(stall_count >= 3) return false;
    } else stall_count = 0;

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
    bool convex = check_convexity(*this, x_trial, free_mask, L, conv_tau);

    // Trust region update
    double rho = (pred > 0) ? actual / pred : -1;
    bool accepted = (rho > 0.1) && convex;
    const char* tr_action;
    if(accepted){
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
      Delta *= 0.25;
      if(Delta < 1e-14 * L) Delta = 1e-14 * L;
      tr_action = convex ? "shrink" : "conv-shrink";
    }

    if(opt_log){
      fprintf(opt_log, "    patch %3d: E=%.6e |g|=%.3e |d|=%.2e D=%.2e rho=%.2f %-14s %s neg=%d\n",
              iter, E, gnorm, dnorm, Delta, rho, step_type,
              tr_action, n_neg_eig);
    }
  }
  return false;  // didn't converge
}

int Deltahedron::reflect_concave(std::span<coord3d> pts, double threshold,
                                  const vector<bool>& fixed) const
{
  bool has_fixed = !fixed.empty();
  int count = 0;
  for(int v = 0; v < N; v++){
    if(has_fixed && fixed[v]) continue;
    int d = degree(v);
    if(d > 6) continue;

    coord3d centroid(0,0,0);
    for(int j = 0; j < d; j++) centroid += pts[nbrs(v)[j]];
    centroid /= (double)d;

    coord3d n_fan(0,0,0);
    for(int j = 0; j < d; j++){
      coord3d e1 = pts[nbrs(v)[j]] - pts[v];
      coord3d e2 = pts[nbrs(v)[(j+1)%d]] - pts[v];
      n_fan += e1.cross(e2);
    }
    double n_len = n_fan.norm();
    if(n_len < 1e-15) continue;
    coord3d n_hat = n_fan / n_len;

    double h = (pts[v] - centroid).dot(n_hat);
    if(h < -threshold){
      pts[v] = centroid + n_hat * (-h);
      count++;
    }
  }
  return count;
}

bool Deltahedron::optimize(std::span<const coord3d> initial_geometry, double target_L, int max_iter, double grad_tol,
                           const vector<bool>& fixed, long long max_work, double angle_tol)
{
  assert((int)initial_geometry.size() == N);
  assert(fixed.empty() || (int)fixed.size() == N);
  points = vector<coord3d>(initial_geometry.begin(), initial_geometry.end());
  const bool has_fixed = !fixed.empty();

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
  const double k_conv   = 0;

  // Two-phase optimization:
  //   Phase 1: k_flat active — settle into flat/equilateral
  //   Phase 2: k_flat off — pure equilateral convergence
  const int phase_budget = max_iter / 2;
  double k_flat = opt_k_flat;
  int phase = (k_flat > 0) ? 1 : 2;
  double phase1_grad_norm0 = 0;

  const double refl_threshold = 0.05 * L;
  const double c1 = 1e-4;        // Armijo parameter

  // Shared helpers
  auto zero_fixed_grad = [&](vector<coord3d>& grad){
    if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) grad[i] = coord3d(0,0,0);
  };

  // Evaluation counters
  n_energy_evals = 0;
  n_grad_evals = 0;
  n_hv_evals = 0;

  // Work budget: total_work = n_energy + N*n_grad + N*n_hv
  auto total_work_fn = [&]() -> long long {
    return (long long)n_energy_evals + (long long)N * n_grad_evals + (long long)N * n_hv_evals;
  };
  if(max_work > 0 && max_iter <= 0) max_iter = INT_MAX;

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

  // Periodic reflection (shared by all methods)
  auto do_reflect = [&]() -> bool {
    return reflect_concave(points, refl_threshold, fixed) > 0;
  };

  // Phase transition logic (shared by all methods)
  // Returns true if phase changed.
  auto check_phase_transition = [&](int iter, int phase_start_iter, const vector<coord3d>& grad) -> bool {
    if(phase != 1 || iter <= phase_start_iter || iter % 50 != 49) return false;
    bool advance = (iter - phase_start_iter >= phase_budget);
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

  const int log_interval = opt_log ? max(1, max_iter / 20) : 0;
  const char* method_name = (opt_method == OptMethod::CG) ? "CG" :
                            (opt_method == OptMethod::LBFGS) ? "LBFGS" : "ST";

  vector<coord3d> grad(N);
  double E = compute_eg(grad);
  phase1_grad_norm0 = vec_norm(grad);

  if(opt_log)
    fprintf(opt_log, "  %s start: E=%.6f |g|=%.4e L=%.4f cv=%.4f ph=%d tol=%.2e\n",
            method_name, E, phase1_grad_norm0, L, edge_cv(), phase, grad_tol);

  bool converged = false;
  int phase_start = 0;

  // ==================== CG ====================
  if(opt_method == OptMethod::CG){
    vector<coord3d> grad_old(N), dir(N), x_trial(N);
    for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0);
    if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir[i] = coord3d(0,0,0);

    for(int iter = 0; iter < max_iter; iter++){
      iterations_used = iter + 1;
      if(max_work > 0 && total_work_fn() >= max_work) break;

      // Periodic reflection
      bool reflected = false;
      if(iter > 0){
        if(do_reflect()){
          E = compute_eg(grad);
          reflected = true;
        }
      }

      // Phase transition
      if(check_phase_transition(iter, phase_start, grad)){
        phase_start = iter;
        E = compute_eg(grad);
        for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0);
        if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir[i] = coord3d(0,0,0);
      }

      // Convergence check
      double gmax = compute_gmax(grad);
      if(gmax < grad_tol){ converged = true; break; }
      if(angle_tol > 0 && max_angle_relerr() < angle_tol && count_concave() == 0){ converged = true; break; }

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
      E = compute_eg(grad);

      // Logging
      if(log_interval > 0 && (iter % log_interval == 0 || iter == max_iter - 1))
        fprintf(opt_log, "  CG %4d: E=%.6f |g|=%.4e gmax*L=%.4e a=%.3e cv=%.4f ph=%d%s\n",
                iter, E, vec_norm(grad), compute_gmax(grad), alpha, edge_cv(), phase,
                reflected ? " R" : "");

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

    for(int iter = 0; iter < max_iter; iter++){
      iterations_used = iter + 1;
      if(max_work > 0 && total_work_fn() >= max_work) break;

      // Periodic reflection
      bool reflected = false;
      if(iter > 0){
        if(do_reflect()){
          E = compute_eg(grad);
          reflected = true;
          // Clear L-BFGS history after reflection (geometry changed externally)
          S.clear(); Y.clear(); rho_hist.clear();
        }
      }

      // Phase transition
      if(check_phase_transition(iter, phase_start, grad)){
        phase_start = iter;
        E = compute_eg(grad);
        S.clear(); Y.clear(); rho_hist.clear();
      }

      // Convergence check
      double gmax = compute_gmax(grad);
      if(gmax < grad_tol){ converged = true; break; }
      if(angle_tol > 0 && max_angle_relerr() < angle_tol && count_concave() == 0){ converged = true; break; }

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
        for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0);
        if(has_fixed) for(int i = 0; i < N; i++) if(fixed[i]) dir[i] = coord3d(0,0,0);
        S.clear(); Y.clear(); rho_hist.clear();
      }

      // Save old state for history update
      grad_old = grad;

      // Line search and update
      double alpha = line_search(E, grad, dir, x_trial);
      for(int i = 0; i < N; i++) points[i] = points[i] + dir[i] * alpha;
      E = compute_eg(grad);

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
      if(log_interval > 0 && (iter % log_interval == 0 || iter == max_iter - 1))
        fprintf(opt_log, "  LB %4d: E=%.6f |g|=%.4e gmax*L=%.4e a=%.3e cv=%.4f ph=%d h=%d%s\n",
                iter, E, vec_norm(grad), compute_gmax(grad), alpha, edge_cv(), phase,
                (int)S.size(), reflected ? " R" : "");
    }
  }

  // ==================== Steihaug-Toint ====================
  else if(opt_method == OptMethod::STEIHAUG){
    double Delta_max = L;
    double Delta = 0.5 * L;
    const int max_inner = min(3 * N, 200);

    // Temp vectors for inner CG
    vector<coord3d> z(N), r_cg(N), d_cg(N), Hd(N), x_trial(N);

    for(int iter = 0; iter < max_iter; iter++){
      iterations_used = iter + 1;
      if(max_work > 0 && total_work_fn() >= max_work) break;

      // Periodic reflection
      bool reflected = false;
      if(iter > 0){
        if(do_reflect()){
          E = compute_eg(grad);
          reflected = true;
        }
      }

      // Phase transition
      if(check_phase_transition(iter, phase_start, grad)){
        phase_start = iter;
        E = compute_eg(grad);
        Delta = 0.5 * L;  // reset trust region on phase change
      }

      // Convergence check
      double gmax = compute_gmax(grad);
      if(gmax < grad_tol){ converged = true; break; }
      if(angle_tol > 0 && max_angle_relerr() < angle_tol && count_concave() == 0){ converged = true; break; }

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
                               k_bond, k_angle, k_curv, k_flat, fixed);

        double kappa = vec_dot(d_cg, Hd);

        if(kappa <= 1e-15 * rr){
          // Negative or zero curvature: step to trust-region boundary along d
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
                             k_bond, k_angle, k_curv, k_flat, fixed);
      double pred = -(vec_dot(grad, z) + 0.5 * vec_dot(z, Hd));

      // Trial point
      for(int i = 0; i < N; i++) x_trial[i] = points[i] + z[i];
      double E_trial = compute_e_only(x_trial);
      double actual = E - E_trial;

      double rho = (pred > 1e-30) ? actual / pred : -1;

      // Accept based on energy reduction only; convexity is handled
      // by periodic reflect_concave (same as CG and L-BFGS).
      bool accepted = (rho > 0.1);

      if(accepted){
        points = x_trial;
        E = compute_eg(grad);
        if(rho > 0.75 && znorm > 0.5 * Delta) Delta = min(2.0 * Delta, Delta_max);
      } else {
        Delta *= 0.25;
        if(Delta < 1e-14 * L) Delta = 1e-14 * L;
      }

      // Logging
      if(log_interval > 0 && (iter % log_interval == 0 || iter == max_iter - 1))
        fprintf(opt_log, "  ST %4d: E=%.6f |g|=%.4e gmax*L=%.4e |z|=%.2e D=%.2e rho=%.2f in=%d ph=%d %s%s\n",
                iter, E, vec_norm(grad), compute_gmax(grad), znorm, Delta, rho, inner_iters,
                phase, accepted ? "acc" : "REJ", reflected ? " R" : "");
    }
  }

  // Final stats
  final_gmax_L = compute_gmax(grad);
  if(opt_log)
    fprintf(opt_log, "  %s done: %d iters, E=%.6f gmax*L=%.4e cv=%.4f %s\n",
            method_name, iterations_used, E, final_gmax_L, edge_cv(),
            converged ? "CONVERGED" : "budget");

  // Post-optimization strict convexity cleanup
  for(int pass = 0; pass < 3; pass++)
    if(reflect_concave(points, 0, fixed) == 0) break;

  final_angle_relerr = max_angle_relerr();
  final_n_concave = count_concave();

  return converged;
}

double Deltahedron::gradient_check(std::span<const coord3d> geometry, double target_L, double eps) const
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

double Deltahedron::hessian_check(std::span<const coord3d> geometry,
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
                                    k_bond, k_angle, k_curv, k_flat, k_conv, interior_mask);
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
