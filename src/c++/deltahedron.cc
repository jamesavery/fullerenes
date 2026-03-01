#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/unfold.hh"
#include "fullerenes/buckinverse.hh"
#include <cmath>
#include <numeric>
#include <queue>

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

    // Forward sweep (modify diag and rhs in-place via copies)
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
static void computeStripCoords(const ExtensionStep& step, vector<coord3d>& points) {
    const auto& strip = step.strip;
    const auto& path = step.path;
    const auto& tp = step.tp;
    int n = (int)strip.size();

    if (step.kind.type == ExpKind::L_type) {
        // L-type: n strip verts, n+1 path/tp verts
        // Tridiagonal system: diag = 5 at endpoints, 6 at interior
        vector<double> diag(n);
        vector<coord3d> rhs(n);

        for (int j = 0; j < n; j++) {
            diag[j] = (j == 0 || j == n - 1) ? 5.0 : 6.0;
            rhs[j] = points[path[j]] + points[path[j + 1]]
                    + points[tp[j]] + points[tp[j + 1]];
        }
        solveTridiagonal(diag, rhs);
        for (int j = 0; j < n; j++)
            points[strip[j]] = rhs[j];

    } else if (step.kind.type == ExpKind::B_type) {
        // B(0,0): 3 strip verts, 5 path verts, 3 tp verts
        assert(step.kind.i == 0 && step.kind.j == 0 && n == 3);
        vector<double> diag = {5.0, 6.0, 5.0};
        vector<coord3d> rhs(3);
        rhs[0] = points[path[0]] + points[path[1]] + points[tp[0]] + points[tp[1]];
        rhs[1] = points[path[1]] + points[path[2]] + points[path[3]] + points[tp[1]];
        rhs[2] = points[path[3]] + points[path[4]] + points[tp[1]] + points[tp[2]];

        solveTridiagonal(diag, rhs);
        for (int j = 0; j < 3; j++)
            points[strip[j]] = rhs[j];

    } else {
        // F-ring: 5 strip verts in a cycle, all deg-6
        assert(n == 5);
        vector<coord3d> rhs(5);
        for (int i = 0; i < 5; i++) {
            int ip1 = (i + 1) % 5;
            rhs[i] = points[path[i]] + points[path[ip1]]
                   + points[tp[i]] + points[tp[ip1]];
        }
        solveCyclicTridiag5(rhs);
        for (int i = 0; i < 5; i++)
            points[strip[i]] = rhs[i];
    }

    // Reproject strip vertices to the unit sphere.
    // The tridiagonal solve gives correct angular positions but places
    // strip vertices interior to the sphere (centroid of points on a sphere
    // is inside it). Since all vertices should lie on a unit sphere centered
    // at origin, normalize each strip vertex to unit length.
    for (int j = 0; j < n; j++) {
        double r = points[strip[j]].norm();
        if (r > 1e-12)
            points[strip[j]] = points[strip[j]] / r;
    }
}

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
static vector<int> matchSeedViaSpiralImpl(
    const neighbours_t& precomp,
    const ExtensionPath& ep)
{
    int seed_N = (int)precomp.size();

    // 1. Canonical spiral of precomputed seed graph
    Triangulation T_pre(precomp, true);
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
    neighbours_t compact_adj(seed_N);
    for (const auto& sv : ep.seed_state) {
        int ci = to_compact[sv.id];
        uint8_t m = sv.active;
        for (int p = 0; p < 6; p++) {
            if (m & (1 << p))
                compact_adj[ci].push_back(to_compact[sv.nbr[p]]);
        }
    }

    // 3. Canonical spiral of ep seed graph
    Triangulation T_ep(compact_adj, true);
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

    // 3. For each expansion step: compute coords then expand topology
    for (const auto& step : ep.steps) {
        computeStripCoords(step, points);
        rd.expand(step);
    }

    // 4. Extract compact Graph with renumbered coordinates
    vector<int> remap(full_N, -1);
    int id = 0;
    for (int u = 0; u < full_N; u++)
        if (rd.alive(u)) remap[u] = id++;

    neighbours_t adj(id);
    vector<coord3d> compact_points(id);
    for (int u = 0; u < full_N; u++) {
        if (!rd.alive(u)) continue;
        compact_points[remap[u]] = points[u];
        uint8_t m = rd.V[u].active;
        for (; m; m &= m - 1)
            adj[remap[u]].push_back(remap[rd.V[u].nbr[__builtin_ctz(m)]]);
    }

    return Deltahedron(Triangulation(adj, true), compact_points);
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

    neighbours_t adj(id);
    vector<coord3d> compact_pts(id);
    for (int u = 0; u < full_N; u++) {
        if (!rd.alive(u)) continue;
        compact_pts[remap[u]] = points[u];
        uint8_t m = rd.V[u].active;
        for (; m; m &= m - 1)
            adj[remap[u]].push_back(remap[rd.V[u].nbr[__builtin_ctz(m)]]);
    }

    return Deltahedron(Triangulation(adj, true), compact_pts);
}

Deltahedron Deltahedron::fromExtensionPathOptimized(const ExtensionPath& ep, int max_iter_per_step) {
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

    // 3. For each expansion step: place strip, expand topology, optimize
    for (const auto& step : ep.steps) {

        // a. Place strip vertices via Laplacian solve
        computeStripCoords(step, points);

        // b. Expand topology
        rd.expand(step);

        // c. Extract compact Deltahedron with current coordinates
        vector<int> remap;
        Deltahedron D = extractCompact(rd, full_N, points, remap);

        // d. Optimize geometry (CG relaxation toward equilateral triangles)
        //    Use limited iterations with loose tolerance for intermediate steps.
        D.optimize(D.points, 0, max_iter_per_step, 1e-8);

        // e. Write optimized coordinates back to full array
        for (int u = 0; u < full_N; u++)
            if (remap[u] >= 0)
                points[u] = D.points[remap[u]];
    }

    // 4. Final extraction
    vector<int> remap;
    return extractCompact(rd, full_N, points, remap);
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
    const vector<coord3d>& x,
    vector<coord3d>* grad,
    double L,           // target edge length
    double k_bond,
    double k_angle,
    double k_curv,
    double k_flat,
    double k_conv)
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
  for(const auto& tri : D.triangles){
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
    int d = (int)D.neighbours[v].size();
    double angle_sum = 0.0;

    // Compute angle sum (and store derivatives if needed)
    vector<coord3d> da_list, dc_list;
    if(grad){ da_list.resize(d); dc_list.resize(d); }

    for(int i = 0; i < d; i++){
      int ni   = D.neighbours[v][i];
      int ni1  = D.neighbours[v][(i+1) % d];
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
        int ni  = D.neighbours[v][i];
        int ni1 = D.neighbours[v][(i+1) % d];

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
      int d = (int)D.neighbours[v].size();
      if(d > 6) continue;

      // 1. Compute face centroids and their mean
      vector<coord3d> fc(d);
      coord3d c_bar(0,0,0);
      for(int i = 0; i < d; i++){
        int ni  = D.neighbours[v][i];
        int ni1 = D.neighbours[v][(i+1) % d];
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
          (*grad)[D.neighbours[v][j]] += (fj + fjpre) * (scale / 3.0);
        }
      }
    }
  }

  // E_conv: smooth convexity bias via softplus.
  // For each qualifying vertex, compute signed height h above neighbor centroid
  // plane (positive = convex). Penalty: k_conv * sigma * log(1 + exp(-h/sigma)).
  // Nearly zero for h > 0, linear in -h for h < 0, smooth everywhere.
  // Exact gradient including normal-rotation term for each neighbour.
  if(k_conv > 0){
    const double sigma = 0.2 * L;  // transition width ~ 20% of edge length
    for(int v = 0; v < N; v++){
      int d = (int)D.neighbours[v].size();
      if(d > 6) continue;

      // Neighbor centroid
      coord3d centroid(0,0,0);
      for(int i = 0; i < d; i++) centroid += x[D.neighbours[v][i]];
      centroid /= (double)d;

      // Fan normal (unnormalized, outward for convex)
      coord3d N_fan(0,0,0);
      for(int i = 0; i < d; i++){
        coord3d e1 = x[D.neighbours[v][i]] - x[v];
        coord3d e2 = x[D.neighbours[v][(i+1)%d]] - x[v];
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
          coord3d ej_prev = x[D.neighbours[v][(j+d-1)%d]] - x[v];
          coord3d ej_next = x[D.neighbours[v][(j+1)%d]]   - x[v];
          coord3d dhdx_nj = n_hat * (-1.0/d)
                          + r_perp.cross(ej_prev - ej_next) / N_len;
          (*grad)[D.neighbours[v][j]] += dhdx_nj * dEdh;
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

// Compute energy only (no gradient), for line search.
static double deltahedron_energy_only(
    const Deltahedron& D,
    const vector<edge_t>& edges,
    const vector<coord3d>& x,
    double L, double k_bond, double k_angle, double k_curv, double k_flat, double k_conv)
{
  return deltahedron_energy_and_gradient(D, edges, x, nullptr, L, k_bond, k_angle, k_curv, k_flat, k_conv);
}

bool Deltahedron::optimize(const vector<coord3d>& initial_geometry, double target_L, int max_iter, double grad_tol)
{
  assert((int)initial_geometry.size() == N);
  points = initial_geometry;

  // Cache edge list (avoid recomputing on every energy evaluation)
  vector<edge_t> edges = undirected_edges();

  // Target edge length: use provided value, or mean of initial edge lengths
  double L = target_L;
  if(L <= 0){
    for(const edge_t& e : edges)
      L += coord3d::dist(points[e.first], points[e.second]);
    L /= edges.size();
  }

  // Force constants (fixed throughout)
  const double k_bond   = 1.0;
  const double k_angle  = 1.0;
  const double k_curv   = 2.0;

  // CG parameters
  const double c1 = 1e-4;        // Armijo parameter

  // Three-phase optimization:
  //   Phase 1: k_conv + k_flat active — flip concave regions (approximate conv gradient)
  //   Phase 2: k_flat active, k_conv off — settle into flat/equilateral (exact gradient)
  //            Early exit when gradient norm drops 100x from phase start.
  //   Phase 3: k_flat off, k_conv off — pure equilateral convergence (exact gradient)
  //            Converges when gradient norm < grad_tol.
  // Each phase gets up to max_iter/3 iterations.
  const int phase_budget = max_iter / 3;
  double k_flat = 2.0;
  double k_conv = 2.0;
  int phase = 1;
  double phase2_grad_norm0 = 0;   // gradient norm at start of phase 2
  int phase1_prev_concave = N;    // concave vertex count at last phase 1 check
  int phase1_stall_count = 0;     // consecutive checks with no concavity improvement

  // Helper: count concave vertices (deg<=6 with negative convexity height)
  auto count_concave = [&]() -> int {
    int n_concave = 0;
    for(int v = 0; v < N; v++){
      int dv = (int)neighbours[v].size();
      if(dv > 6) continue;
      coord3d cent(0,0,0);
      for(int j = 0; j < dv; j++) cent += points[neighbours[v][j]];
      cent /= (double)dv;
      coord3d nfan(0,0,0);
      for(int j = 0; j < dv; j++){
        coord3d e1 = points[neighbours[v][j]] - points[v];
        coord3d e2 = points[neighbours[v][(j+1)%dv]] - points[v];
        nfan += e1.cross(e2);
      }
      double nl = nfan.norm();
      if(nl < 1e-15) continue;
      if((points[v] - cent).dot(nfan / nl) < 0) n_concave++;
    }
    return n_concave;
  };

  // Helper: reset CG with current force constants
  auto reset_cg = [&](double& E, vector<coord3d>& grad, vector<coord3d>& dir){
    E = deltahedron_energy_and_gradient(*this, edges, points, &grad, L,
                                         k_bond, k_angle, k_curv, k_flat, k_conv);
    for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0);
  };

  vector<coord3d> grad(N), grad_old(N), dir(N), x_trial(N);

  double E;
  reset_cg(E, grad, dir);

  bool converged = false;
  int phase_start = 0;

  for(int iter = 0; iter < max_iter; iter++){
    iterations_used = iter + 1;

    // Phase transition checks (every 50 iterations)
    if(phase < 3 && iter > phase_start && iter % 50 == 49){
      bool advance = (iter - phase_start >= phase_budget);

      if(!advance && phase == 1){
        int nc = count_concave();
        if(nc == 0){
          advance = true;  // all convex — move on
        } else if(nc >= phase1_prev_concave){
          phase1_stall_count++;
          if(phase1_stall_count >= 2)
            advance = true;  // no improvement for 100 iters — convexity term settled
        } else {
          phase1_stall_count = 0;  // reset on improvement
        }
        phase1_prev_concave = nc;
      }

      // Phase 2 early exit: gradient dropped 100x from phase start
      if(!advance && phase == 2 && phase2_grad_norm0 > 0){
        double gn = sqrt(vec_dot(grad, grad));
        if(gn < phase2_grad_norm0 * 0.01) advance = true;
      }

      if(advance){
        if(phase == 1){
          k_conv = 0;
          phase = 2;
          phase_start = iter;
          reset_cg(E, grad, dir);
          phase2_grad_norm0 = sqrt(vec_dot(grad, grad));
        } else {
          k_flat = 0;
          phase = 3;
          phase_start = iter;
          reset_cg(E, grad, dir);
        }
      }
    }

    double grad_norm2 = vec_dot(grad, grad);
    if(grad_norm2 < grad_tol * grad_tol){
      converged = true;
      break;
    }

    // Backtracking line search (Armijo condition)
    double slope = vec_dot(grad, dir);  // should be negative
    if(slope > 0){
      // Search direction is not a descent direction; reset to steepest descent
      for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0);
      slope = -grad_norm2;
    }

    double alpha = 1.0;
    for(int ls = 0; ls < 60; ls++){
      for(int i = 0; i < N; i++) x_trial[i] = points[i] + dir[i] * alpha;
      double E_trial = deltahedron_energy_only(*this, edges, x_trial, L,
                                               k_bond, k_angle, k_curv, k_flat, k_conv);
      if(E_trial <= E + c1 * alpha * slope) break;
      alpha *= 0.5;
    }

    // Update position
    for(int i = 0; i < N; i++) points[i] = points[i] + dir[i] * alpha;

    // Save old gradient
    grad_old = grad;

    // New energy and gradient
    E = deltahedron_energy_and_gradient(*this, edges, points, &grad, L,
                                        k_bond, k_angle, k_curv, k_flat, k_conv);

    // Polak-Ribiere beta (with restart: beta = max(0, ...))
    double gg_old = vec_dot(grad_old, grad_old);
    double beta = 0.0;
    if(gg_old > 1e-30){
      // PR formula: beta = dot(g, g - g_old) / dot(g_old, g_old)
      vector<coord3d> gdiff(N);
      for(int i = 0; i < N; i++) gdiff[i] = grad[i] - grad_old[i];
      beta = max(0.0, vec_dot(grad, gdiff) / gg_old);
    }

    // Update search direction
    for(int i = 0; i < N; i++)
      dir[i] = grad[i] * (-1.0) + dir[i] * beta;
  }

  // Post-optimization fix: reflect any deeply concave vertices through their
  // neighbor centroid plane. The equilateral energy is roughly symmetric around
  // h=0, so if we landed in a concave local minimum at h=-d, there should be a
  // convex one near h=+d. After reflecting, re-run a short optimization with
  // k_conv active to settle into the convex basin.
  {
    bool any_reflected = false;
    for(int v = 0; v < N; v++){
      int d = (int)neighbours[v].size();
      if(d > 6) continue;

      coord3d centroid(0,0,0);
      for(int j = 0; j < d; j++) centroid += points[neighbours[v][j]];
      centroid /= (double)d;

      coord3d n_fan(0,0,0);
      for(int j = 0; j < d; j++){
        coord3d e1 = points[neighbours[v][j]] - points[v];
        coord3d e2 = points[neighbours[v][(j+1)%d]] - points[v];
        n_fan += e1.cross(e2);
      }
      double n_len = n_fan.norm();
      if(n_len < 1e-15) continue;
      coord3d n_hat = n_fan / n_len;

      double h = (points[v] - centroid).dot(n_hat);
      if(h < -0.1 * L){  // deeply concave: more than 10% of edge length
        // Reflect through centroid plane: place vertex at centroid + |h| * n_hat
        points[v] = centroid + n_hat * (-h);
        any_reflected = true;
      }
    }

    if(any_reflected){
      // Re-optimize with k_conv active, using remaining iteration budget
      // or at least max_iter/3 more iterations.
      int reopt_budget = max(max_iter / 3, max_iter - iterations_used);
      k_conv = 2.0;
      k_flat = 2.0;
      phase = 1;
      phase_start = iterations_used;
      phase1_prev_concave = N;
      phase1_stall_count = 0;

      reset_cg(E, grad, dir);

      converged = false;
      int reopt_start = iterations_used;
      for(int iter = 0; iter < reopt_budget; iter++){
        iterations_used = reopt_start + iter + 1;

        // Phase transitions (same logic as main loop)
        if(phase < 3 && iter > 0 && iter % 50 == 49){
          bool advance = (iter >= phase_budget);

          if(!advance && phase == 1){
            int nc = count_concave();
            if(nc == 0) advance = true;
            else if(nc >= phase1_prev_concave){
              phase1_stall_count++;
              if(phase1_stall_count >= 2) advance = true;
            } else {
              phase1_stall_count = 0;
            }
            phase1_prev_concave = nc;
          }

          if(!advance && phase == 2 && phase2_grad_norm0 > 0){
            double gn = sqrt(vec_dot(grad, grad));
            if(gn < phase2_grad_norm0 * 0.01) advance = true;
          }

          if(advance){
            if(phase == 1){
              k_conv = 0; phase = 2;
              reset_cg(E, grad, dir);
              phase2_grad_norm0 = sqrt(vec_dot(grad, grad));
            } else {
              k_flat = 0; phase = 3;
              reset_cg(E, grad, dir);
            }
          }
        }

        double grad_norm2 = vec_dot(grad, grad);
        if(grad_norm2 < grad_tol * grad_tol){
          converged = true;
          break;
        }

        double slope = vec_dot(grad, dir);
        if(slope > 0){
          for(int i = 0; i < N; i++) dir[i] = grad[i] * (-1.0);
          slope = -grad_norm2;
        }

        double alpha = 1.0;
        for(int ls = 0; ls < 60; ls++){
          for(int i = 0; i < N; i++) x_trial[i] = points[i] + dir[i] * alpha;
          double E_trial = deltahedron_energy_only(*this, edges, x_trial, L,
                                                   k_bond, k_angle, k_curv, k_flat, k_conv);
          if(E_trial <= E + c1 * alpha * slope) break;
          alpha *= 0.5;
        }

        for(int i = 0; i < N; i++) points[i] = points[i] + dir[i] * alpha;
        grad_old = grad;
        E = deltahedron_energy_and_gradient(*this, edges, points, &grad, L,
                                            k_bond, k_angle, k_curv, k_flat, k_conv);

        double gg_old = vec_dot(grad_old, grad_old);
        double beta = 0.0;
        if(gg_old > 1e-30){
          vector<coord3d> gdiff(N);
          for(int i = 0; i < N; i++) gdiff[i] = grad[i] - grad_old[i];
          beta = max(0.0, vec_dot(grad, gdiff) / gg_old);
        }
        for(int i = 0; i < N; i++)
          dir[i] = grad[i] * (-1.0) + dir[i] * beta;
      }
    }
  }

  return converged;
}

double Deltahedron::gradient_check(const vector<coord3d>& geometry, double target_L, double eps) const
{
  vector<coord3d> x = geometry;
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
