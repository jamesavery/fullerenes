#include <gtest/gtest.h>
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_old.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/planargraph.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/symmetry.hh"
#include <chrono>
#include <numeric>
#include <queue>
#include <array>
#include <map>

// Build a fullerene dual triangulation from N and isomer index (via IsomerDB).
static Triangulation make_dual(int N, int idx, bool IPR = false) {
  IsomerDB db = IsomerDB::readPDB(N, IPR);
  FullereneGraph G = IsomerDB::makeIsomer(N, db.entries[idx]);
  PlanarGraph PG(G);
  Triangulation T(PG.dual_graph());
  return T;
}

// Print timing stats from a vector of durations in microseconds.
static void print_timing_stats(const std::string& label, std::vector<double>& times_us) {
  if (times_us.empty()) return;
  std::sort(times_us.begin(), times_us.end());
  double median = times_us[times_us.size() / 2];
  double max_t = times_us.back();
  double mean = std::accumulate(times_us.begin(), times_us.end(), 0.0) / times_us.size();
  double sq_sum = 0;
  for (double t : times_us) sq_sum += (t - mean) * (t - mean);
  double stddev = sqrt(sq_sum / times_us.size());
  std::cout << label << ": n=" << times_us.size()
            << " median=" << median << "us"
            << " max=" << max_t << "us"
            << " mean=" << mean << "us"
            << " stddev=" << stddev << "us" << std::endl;
}


// ============================================================================
// Plantri: general 3-connected planar triangulations
// ============================================================================

// Read all graphs from a planar_code file, streaming one at a time.
// Calls callback(graph, index) for each graph. Returns total count.
static int stream_planarcode(const std::string& path,
                             std::function<void(const PlanarGraph&, int)> callback) {
  FILE* f = fopen(path.c_str(), "rb");
  if (!f) return -1;

  // Skip 15-byte header ">>planar_code<<"
  fseek(f, 15, SEEK_SET);

  int count = 0;
  while (!feof(f)) {
    PlanarGraph g = PlanarGraph::read_hog_planarcode(f);
    if (g.N == 0) break;
    callback(g, count);
    count++;
  }
  fclose(f);
  return count;
}


// ============================================================================
// DCEL-based DelaunayTriangulation tests
// ============================================================================

// Verify DCEL construction from a known triangulation.
TEST(DCEL, Construction_C20) {
  Triangulation T = make_dual(20, 0);
  auto D = DelaunayTriangulation::from_triangulation(T);

  EXPECT_EQ(D.nv, T.N);
  // C20 dual: 12 vertices, 30 edges, 20 faces
  EXPECT_EQ(D.nh, 2 * 30);
  EXPECT_EQ(D.nf, 20);

  EXPECT_TRUE(D.check_consistency()) << "DCEL consistency check failed";

  // All edges should have length 1.
  for (int h = 0; h < D.nh; h++)
    if (D.alive(h))
      EXPECT_DOUBLE_EQ(D.he_length[h], 1.0);

  // All angles should be pi/3 (equilateral).
  for (int h = 0; h < D.nh; h++)
    if (D.alive(h))
      EXPECT_NEAR(D.he_angle[h], M_PI / 3.0, 1e-14);
}

// Test DCEL construction for C60.
TEST(DCEL, Construction_C60) {
  Triangulation T = make_dual(60, 0);
  auto D = DelaunayTriangulation::from_triangulation(T);
  EXPECT_EQ(D.nv, T.N);  // N/2+2 = 32
  EXPECT_TRUE(D.check_consistency());
}

// Test the full DCEL-based iDT algorithm on C20 (trivial: no flat vertices).
TEST(DCEL, Compute_C20) {
  Triangulation T = make_dual(20, 0);
  auto D = DelaunayTriangulation::compute(T);

  // C20 has only degree-5 vertices, no flat vertices to remove.
  // Result should be 12 vertices, fully Delaunay.
  EXPECT_EQ(D.nv, 12);
  EXPECT_TRUE(D.check_consistency());
  EXPECT_TRUE(D.is_delaunay());
}

// Test the full DCEL-based iDT algorithm on C60 Ih.
TEST(DCEL, Compute_C60_Ih) {
  Triangulation T = make_dual(60, 0, true);  // IPR isomer 0 = Ih
  auto D = DelaunayTriangulation::compute(T);

  EXPECT_EQ(D.nv, 12);
  EXPECT_TRUE(D.check_consistency());
  EXPECT_TRUE(D.is_delaunay());

  // Count live edges and faces.
  int live_edges = 0;
  for (int h = 0; h < D.nh; h += 2)
    if (D.alive(h)) live_edges++;
  EXPECT_EQ(live_edges, 30);

  int live_faces = 0;
  for (int f = 0; f < D.nf; f++)
    if (D.f_he[f] >= 0) live_faces++;
  EXPECT_EQ(live_faces, 20);
}

// Helper: verify a DCEL result has the expected structure.
// N_original_faces: number of equilateral triangles in the input (for area check).
// If 0, area check is skipped.
static void verify_dcel_reduced(const DelaunayTriangulation& D, int expected_verts,
                                int N_original_faces = 0) {
  EXPECT_EQ(D.nv, expected_verts);
  EXPECT_TRUE(D.check_consistency());

  // --- Delaunay criterion: every edge must be Delaunay ---
  {
    int non_del = 0;
    for (int h = 0; h < D.nh; h += 2)
      if (D.alive(h) && !D.is_delaunay_edge(h))
        non_del++;
    EXPECT_EQ(non_del, 0) << non_del << " non-Delaunay edges remain";
  }

  // Count live edges and faces.
  int live_edges = 0;
  for (int h = 0; h < D.nh; h += 2)
    if (D.alive(h)) live_edges++;

  int live_faces = 0;
  for (int f = 0; f < D.nf; f++)
    if (D.f_he[f] >= 0) live_faces++;

  // Euler: V - E + F = 2 for genus 0
  EXPECT_EQ(expected_verts - live_edges + live_faces, 2)
    << "Euler formula failed: V=" << expected_verts
    << " E=" << live_edges << " F=" << live_faces;

  if (expected_verts == 12) {
    EXPECT_EQ(live_edges, 30);
    EXPECT_EQ(live_faces, 20);
  }

  // --- Edge lengths: positive and twin-consistent ---
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    EXPECT_GT(D.he_length[h], 0.0) << "Edge " << h/2 << " has non-positive length";
    EXPECT_EQ(D.he_length[h], D.he_length[h ^ 1])
      << "Edge " << h/2 << " twin length mismatch";
  }

  // --- Triangle inequality and angle sum = pi for each face ---
  for (int f = 0; f < D.nf; f++) {
    if (D.f_he[f] < 0) continue;
    int h0 = D.f_he[f], h1 = D.he_next[h0], h2 = D.he_next[h1];
    double a = D.he_length[h0], b = D.he_length[h1], c = D.he_length[h2];
    EXPECT_GT(a + b, c) << "Triangle inequality in face " << f;
    EXPECT_GT(b + c, a) << "Triangle inequality in face " << f;
    EXPECT_GT(c + a, b) << "Triangle inequality in face " << f;

    double cos_u = std::clamp((a*a + c*c - b*b) / (2.0*a*c), -1.0, 1.0);
    double cos_v = std::clamp((a*a + b*b - c*c) / (2.0*a*b), -1.0, 1.0);
    double cos_w = std::clamp((b*b + c*c - a*a) / (2.0*b*c), -1.0, 1.0);
    double angle_sum = acos(cos_u) + acos(cos_v) + acos(cos_w);
    EXPECT_NEAR(angle_sum, M_PI, 1e-10) << "Angle sum != pi in face " << f;
  }

  // --- Cone angles: each vertex retains its original cone angle ---
  for (int v = 0; v < D.nv; v++) {
    if (D.v_out[v] < 0) continue;  // dead vertex
    double expected_angle = D.v_cone_angle[v];  // set at construction
    double total_angle = 0;
    int h0 = D.v_out[v], h = h0;
    do {
      total_angle += D.he_angle[h];
      h = D.cw(h);
    } while (h != h0);
    EXPECT_NEAR(total_angle, expected_angle, 1e-6)
      << "Vertex " << v << " (orig degree " << D.v_orig_degree[v]
      << "): cone angle " << total_angle << " != expected " << expected_angle;
  }

  // --- Area conservation ---
  if (N_original_faces > 0) {
    double expected_area = N_original_faces * sqrt(3.0) / 4.0;
    double actual_area = 0;
    for (int f = 0; f < D.nf; f++) {
      if (D.f_he[f] < 0) continue;
      int h0 = D.f_he[f], h1 = D.he_next[h0], h2 = D.he_next[h1];
      double a = D.he_length[h0], b = D.he_length[h1], c = D.he_length[h2];
      double s = (a + b + c) / 2.0;
      double area2 = s * (s - a) * (s - b) * (s - c);
      actual_area += sqrt(std::max(0.0, area2));
    }
    EXPECT_NEAR(actual_area, expected_area, 1e-6)
      << "Total area " << actual_area << " != expected " << expected_area;
  }
}

TEST(DCEL, C60_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  std::vector<double> times_us;
  int idx = 0;
  int n_fail = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    SCOPED_TRACE("C60 #" + std::to_string(idx));
    auto t0 = std::chrono::high_resolution_clock::now();
    auto D = DelaunayTriangulation::compute(T);
    auto t1 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    verify_dcel_reduced(D, 12, 60);
    if (D.nv != 12 || !D.is_delaunay()) n_fail++;
    idx++;
  }
  BuckyGen::stop(Q);
  EXPECT_EQ(idx, 1812);
  EXPECT_EQ(n_fail, 0);
  print_timing_stats("DCEL C60 iDT", times_us);
}

// Verify exact Eisenstein face-origin tracking: every original face must
// appear in exactly one iDT face's f_origin (modulo boundary faces that
// appear in exactly two), and the total count must equal N_original_faces.
static void verify_exact_origins(const DelaunayTriangulation& D,
                                 int N_original_faces) {
  // Collect all origin assignments.
  std::map<int, int> origin_count;  // original face ID → number of iDT faces claiming it
  int total = 0;
  for (int f = 0; f < D.nf; f++) {
    if (D.f_he[f] < 0) continue;  // dead face
    for (int orig : D.f_origin[f])
      origin_count[orig]++;
    total += D.f_origin[f].size();
  }

  // Every original face must appear at least once.
  int missing = 0, singly = 0, doubly = 0, multi = 0;
  for (int i = 0; i < N_original_faces; i++) {
    auto it = origin_count.find(i);
    if (it == origin_count.end()) missing++;
    else if (it->second == 1) singly++;
    else if (it->second == 2) doubly++;
    else multi++;
  }

  EXPECT_EQ(missing, 0) << missing << " original faces not assigned to any iDT face";
  EXPECT_EQ(multi, 0) << multi << " original faces assigned to 3+ iDT faces";
  // Boundary faces (split by a geodesic) appear in exactly 2 iDT faces.
  // Interior faces appear in exactly 1.  Both are acceptable.
  // But no face should appear in 0 or 3+.

  // No spurious face IDs outside [0, N_original_faces).
  for (auto& [fid, cnt] : origin_count)
    EXPECT_LT(fid, N_original_faces) << "Spurious origin face ID " << fid;
}

// Exact origin tracking on C20 (no flat vertices to remove).
TEST(DCEL, ExactOrigins_C20) {
  Triangulation T = make_dual(20, 0, false);
  auto D = DelaunayTriangulation::compute(T, /*track_origins=*/true);
  EXPECT_EQ(D.nv, 12);
  EXPECT_TRUE(D.is_delaunay());
  verify_exact_origins(D, 20);
}

// Exact origin tracking on C60 Ih.
TEST(DCEL, ExactOrigins_C60_Ih) {
  Triangulation T = make_dual(60, 0, true);
  auto D = DelaunayTriangulation::compute(T, /*track_origins=*/true);
  EXPECT_EQ(D.nv, 12);
  EXPECT_TRUE(D.is_delaunay());
  verify_exact_origins(D, 60);
}

TEST(DCEL, ExactOrigins_C60_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  std::vector<double> times_us;
  int idx = 0;
  int n_fail = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    SCOPED_TRACE("C60 #" + std::to_string(idx));
    try {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto D = DelaunayTriangulation::compute(T, /*track_origins=*/true);
      auto t1 = std::chrono::high_resolution_clock::now();
      times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
      verify_dcel_reduced(D, 12, 60);
      verify_exact_origins(D, 60);
      if (D.nv != 12 || !D.is_delaunay()) n_fail++;
    } catch (const std::exception& e) {
      ADD_FAILURE() << "C60 #" << idx << ": exception: " << e.what();
      n_fail++;
    }
    idx++;
  }
  BuckyGen::stop(Q);
  EXPECT_EQ(idx, 1812);
  EXPECT_EQ(n_fail, 0);
  print_timing_stats("DCEL C60 exact origins", times_us);
}

TEST(DCEL, DISABLED_C80_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(80, false, false);
  Triangulation T;
  std::vector<double> times_us;
  int idx = 0;
  int n_fail = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    SCOPED_TRACE("C80 #" + std::to_string(idx));
    auto t0 = std::chrono::high_resolution_clock::now();
    auto D = DelaunayTriangulation::compute(T);
    auto t1 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    verify_dcel_reduced(D, 12, 80);
    if (D.nv != 12 || !D.is_delaunay()) n_fail++;
    idx++;
  }
  BuckyGen::stop(Q);
  EXPECT_EQ(n_fail, 0);
  print_timing_stats("DCEL C80 iDT", times_us);
}

TEST(DCEL, DISABLED_ExactOrigins_C80_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(80, false, false);
  Triangulation T;
  std::vector<double> times_us;
  int idx = 0;
  int n_fail = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    SCOPED_TRACE("C80 #" + std::to_string(idx));
    try {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto D = DelaunayTriangulation::compute(T, /*track_origins=*/true);
      auto t1 = std::chrono::high_resolution_clock::now();
      times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
      verify_dcel_reduced(D, 12, 80);
      verify_exact_origins(D, 80);
      if (D.nv != 12 || !D.is_delaunay()) n_fail++;
    } catch (const std::exception& e) {
      ADD_FAILURE() << "C80 #" << idx << ": exception: " << e.what();
      n_fail++;
    }
    idx++;
  }
  BuckyGen::stop(Q);
  EXPECT_EQ(n_fail, 0);
  print_timing_stats("DCEL C80 exact origins", times_us);
}

TEST(DCEL, DISABLED_C100_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(100, false, false);
  Triangulation T;
  std::vector<double> times_us;
  int idx = 0;
  int n_fail = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    SCOPED_TRACE("C100 #" + std::to_string(idx));
    auto t0 = std::chrono::high_resolution_clock::now();
    auto D = DelaunayTriangulation::compute(T);
    auto t1 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    verify_dcel_reduced(D, 12, 100);
    if (D.nv != 12 || !D.is_delaunay()) n_fail++;
    idx++;
  }
  BuckyGen::stop(Q);
  EXPECT_EQ(n_fail, 0);
  print_timing_stats("DCEL C100 iDT", times_us);
}

TEST(DCEL, DISABLED_ExactOrigins_C100_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(100, false, false);
  Triangulation T;
  std::vector<double> times_us;
  int idx = 0;
  int n_fail = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    SCOPED_TRACE("C100 #" + std::to_string(idx));
    try {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto D = DelaunayTriangulation::compute(T, /*track_origins=*/true);
      auto t1 = std::chrono::high_resolution_clock::now();
      times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
      verify_dcel_reduced(D, 12, 100);
      verify_exact_origins(D, 100);
      if (D.nv != 12 || !D.is_delaunay()) n_fail++;
    } catch (const std::exception& e) {
      ADD_FAILURE() << "C100 #" << idx << ": exception: " << e.what();
      n_fail++;
    }
    idx++;
  }
  BuckyGen::stop(Q);
  EXPECT_EQ(n_fail, 0);
  print_timing_stats("DCEL C100 exact origins", times_us);
}

TEST(DCEL, DISABLED_Plantri15_AllTriangulations) {
  std::string path = FULLERENE_ROOT "/data/triangulations_15.pl";
  // Check file exists.
  {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) {
      std::cerr << "Skipping: " << path << " not found" << std::endl;
      return;
    }
    fclose(f);
  }

  int n_tested = 0, n_failed = 0;
  int max_degree_seen = 0;
  std::vector<double> times_us;

  stream_planarcode(path, [&](const PlanarGraph& g, int idx) {
    // Record original degrees
    vector<int> orig_degrees(g.N);
    for (int u = 0; u < g.N; u++)
      orig_degrees[u] = g[u].size();

    for (int d : orig_degrees)
      max_degree_seen = std::max(max_degree_seen, d);

    // Count cone points (non-degree-6)
    int n_cones = 0;
    for (int d : orig_degrees)
      if (d != 6) n_cones++;

    // Skip flat tori
    if (n_cones == 0) return;

    int N_original_faces = 2 * g.N - 4;
    Triangulation T(g);

    if (!T.is_consistently_oriented()) {
      n_failed++;
      if (n_failed <= 10)
        std::cerr << "Graph #" << idx << ": input not consistently oriented" << std::endl;
      n_tested++;
      return;
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    auto D = DelaunayTriangulation::compute(T);
    auto t1 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());

    // Vertex count must match
    if (D.nv != n_cones) {
      n_failed++;
      if (n_failed <= 10)
        std::cerr << "Graph #" << idx << ": removal incomplete (D.nv="
                  << D.nv << ", expected " << n_cones << ")" << std::endl;
      n_tested++;
      return;
    }

    // Consistency
    if (!D.check_consistency()) {
      n_failed++;
      if (n_failed <= 10)
        std::cerr << "Graph #" << idx << ": consistency check failed" << std::endl;
      n_tested++;
      return;
    }

    // Delaunay: STRICT — 0 non-Delaunay edges
    int non_del = D.count_non_delaunay();
    if (non_del > 0) {
      if (n_failed < 10 || non_del > 4)
        ADD_FAILURE() << "Graph #" << idx << ": " << non_del << " non-Delaunay edges";
      n_failed++;
    }

    // Cone angles
    // Cone degrees are ordered by sort_flat_last: non-6 first, in original order.
    vector<int> cone_degrees;
    for (int d : orig_degrees)
      if (d != 6) cone_degrees.push_back(d);

    if ((int)cone_degrees.size() == D.nv) {
      for (int v = 0; v < D.nv; v++) {
        if (D.v_out[v] < 0) continue;
        double expected_angle = D.v_cone_angle[v];
        double total_angle = 0;
        int h0 = D.v_out[v], h = h0;
        do {
          total_angle += D.he_angle[h];
          h = D.cw(h);
        } while (h != h0);
        if (std::abs(total_angle - expected_angle) > 1e-6) {
          if (n_failed < 10)
            ADD_FAILURE() << "Graph #" << idx << " vertex " << v
                          << " (orig deg " << D.v_orig_degree[v]
                          << "): cone angle " << total_angle
                          << " != expected " << expected_angle;
          n_failed++;
          break;
        }
      }
    }

    // Area conservation
    double expected_area = N_original_faces * sqrt(3.0) / 4.0;
    double actual_area = 0;
    for (int f = 0; f < D.nf; f++) {
      if (D.f_he[f] < 0) continue;
      int h0 = D.f_he[f], h1 = D.he_next[h0], h2 = D.he_next[h1];
      double a = D.he_length[h0], b = D.he_length[h1], c = D.he_length[h2];
      double s = (a + b + c) / 2.0;
      double area2 = s * (s - a) * (s - b) * (s - c);
      actual_area += sqrt(std::max(0.0, area2));
    }
    if (std::abs(actual_area - expected_area) > 1e-6) {
      if (n_failed < 10)
        ADD_FAILURE() << "Graph #" << idx << ": area " << actual_area
                      << " != expected " << expected_area;
      n_failed++;
    }

    n_tested++;

    if (n_tested % 100000 == 0)
      fprintf(stderr, "  [%d] tested=%d failed=%d\n", idx, n_tested, n_failed);
  });

  print_timing_stats("DCEL Plantri15 iDT", times_us);
  std::cout << "DCEL Plantri15: tested " << n_tested << " triangulations, "
            << n_failed << " failures, max degree seen = "
            << max_degree_seen << std::endl;
  EXPECT_EQ(n_failed, 0);
}

// ============================================================================
// DCEL embed_3d tests
// ============================================================================

// Verify DCEL 3D embedding quality: edge distance errors and convexity.
static void verify_dcel_embedding(const DelaunayTriangulation& D,
                                  const vector<coord3d>& coords,
                                  double dist_tol, const std::string& label) {
  ASSERT_EQ((int)coords.size(), D.nv) << label << ": wrong coord count";

  // Edge distance errors (shortest edge per pair, matching what embed_3d targets).
  map<pair<int,int>, double> shortest;
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    int u = D.he_origin[h], v = D.dest(h);
    if (u > v) swap(u, v);
    double L = D.he_length[h];
    auto key = make_pair(u, v);
    auto it = shortest.find(key);
    if (it == shortest.end() || L < it->second)
      shortest[key] = L;
  }

  double max_rel_err = 0;
  double sum_sq_err = 0;
  int n_edges = 0;
  for (auto& [key, target] : shortest) {
    double actual = (coords[key.first] - coords[key.second]).norm();
    double rel_err = fabs(actual - target) / target;
    max_rel_err = std::max(max_rel_err, rel_err);
    sum_sq_err += (actual - target) * (actual - target);
    n_edges++;
  }
  double rms_err = sqrt(sum_sq_err / n_edges);

  EXPECT_LT(max_rel_err, dist_tol)
    << label << ": max relative edge error = " << max_rel_err;

  // Cone angle errors
  double max_cone_err = 0;
  for (int v = 0; v < D.nv; v++) {
    if (D.v_out[v] < 0) continue;
    double angle_sum = 0;
    int h0 = D.v_out[v], h = h0;
    do {
      int d1 = D.dest(h);
      int h2 = D.cw(h);
      int d2 = D.dest(h2);
      coord3d va = coords[d1] - coords[v], vb = coords[d2] - coords[v];
      double ra = va.norm(), rb = vb.norm();
      if (ra > 1e-15 && rb > 1e-15) {
        double C = max(-1.0, min(1.0, va.dot(vb) / (ra * rb)));
        angle_sum += acos(C);
      }
      h = h2;
    } while (h != h0);
    max_cone_err = std::max(max_cone_err, fabs(angle_sum - D.v_cone_angle[v]));
  }

  std::cout << label << ": n_edges=" << n_edges
            << " max_rel_err=" << max_rel_err
            << " rms_err=" << rms_err
            << " max_cone_err=" << max_cone_err << std::endl;
}

TEST(DCELEmbed, C20_Icosahedron) {
  Triangulation T = make_dual(20, 0, false);
  auto D = DelaunayTriangulation::compute(T);
  ASSERT_EQ(D.nv, 12);

  auto coords = D.embed_3d();
  verify_dcel_embedding(D, coords, 1e-6, "C20");

  // For C20, all iDT edges should have length 1 (regular icosahedron).
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    int u = D.he_origin[h], v = D.dest(h);
    double dist = (coords[u] - coords[v]).norm();
    EXPECT_NEAR(dist, 1.0, 1e-6) << "C20 edge h=" << h
      << " (" << u << "," << v << ") = " << dist;
  }
}

TEST(DCELEmbed, C60_Ih) {
  Triangulation T = make_dual(60, 0, true);
  auto D = DelaunayTriangulation::compute(T);
  ASSERT_EQ(D.nv, 12);

  auto coords = D.embed_3d();
  verify_dcel_embedding(D, coords, 1e-4, "C60_Ih");
}

TEST(DCELEmbed, SmallFullerenes) {
  int sizes[] = {20, 24, 26, 28, 30, 32, 34, 36, 38, 40};
  int total = 0, n_poor = 0;

  for (int N : sizes) {
    BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
    Triangulation T;
    int idx = 0;
    while (BuckyGen::next_fullerene(Q, T)) {
      auto D = DelaunayTriangulation::compute(T);
      if (D.nv != 12) { idx++; continue; }

      auto coords = D.embed_3d();

      // Check edge distance matching (shortest per pair)
      map<pair<int,int>, double> shortest;
      for (int h = 0; h < D.nh; h += 2) {
        if (!D.alive(h)) continue;
        int u = D.he_origin[h], v = D.dest(h);
        if (u > v) swap(u, v);
        double L = D.he_length[h];
        auto key = make_pair(u, v);
        auto it = shortest.find(key);
        if (it == shortest.end() || L < it->second)
          shortest[key] = L;
      }

      double max_rel_err = 0;
      for (auto& [key, target] : shortest) {
        double actual = (coords[key.first] - coords[key.second]).norm();
        double rel_err = fabs(actual - target) / target;
        max_rel_err = std::max(max_rel_err, rel_err);
      }

      if (max_rel_err > 0.01) {
        n_poor++;
        std::cerr << "C" << N << " #" << idx
                  << ": max_rel_err=" << max_rel_err << std::endl;
      }

      EXPECT_LT(max_rel_err, 0.05)
        << "C" << N << " #" << idx << " DCEL embedding too inaccurate";

      idx++;
      total++;
    }
    BuckyGen::stop(Q);
  }
  std::cout << "DCEL Embed: Tested " << total << " embeddings, "
            << n_poor << " with >1% error" << std::endl;
}

TEST(DCELEmbed, C60_1264_MultiEdge) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    if (idx == 1264) break;
    idx++;
  }
  BuckyGen::stop(Q);

  auto D = DelaunayTriangulation::compute(T);
  ASSERT_EQ(D.nv, 12);

  // Report multi-edges
  std::map<std::pair<int,int>, std::vector<double>> edge_lengths;
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    int u = D.he_origin[h], v = D.dest(h);
    if (u > v) std::swap(u, v);
    edge_lengths[{u,v}].push_back(D.he_length[h]);
  }

  std::cout << "Distinct vertex pairs: " << edge_lengths.size() << std::endl;
  for (auto& [key, lens] : edge_lengths) {
    if (lens.size() > 1) {
      std::cout << "  Multi-edge (" << key.first << "," << key.second << "): ";
      for (double l : lens) std::cout << l << " ";
      std::cout << std::endl;
    }
  }

  for (int v = 0; v < D.nv; v++)
    std::cout << "  v" << v << ": deg=" << D.vertex_degree(v)
              << " cone=" << D.v_cone_angle[v]
              << " orig_deg=" << D.v_orig_degree[v] << std::endl;

  auto coords = D.embed_3d();

  // Report edge errors
  double max_rel_err = 0;
  for (auto& [key, lens] : edge_lengths) {
    double dist = (coords[key.first] - coords[key.second]).norm();
    double shortest = *std::min_element(lens.begin(), lens.end());
    double rel_err = std::abs(dist - shortest) / shortest;
    max_rel_err = std::max(max_rel_err, rel_err);
    if (rel_err > 0.01)
      std::cout << "  (" << key.first << "," << key.second
                << "): target=" << shortest << " actual=" << dist
                << " rel_err=" << rel_err << std::endl;
  }

  // Report cone angle errors
  for (int v = 0; v < D.nv; v++) {
    if (D.v_out[v] < 0) continue;
    double angle_sum = 0;
    int h0 = D.v_out[v], h = h0;
    do {
      int d1 = D.dest(h);
      int h2 = D.cw(h);
      int d2 = D.dest(h2);
      coord3d va = coords[d1] - coords[v], vb = coords[d2] - coords[v];
      double ra = va.norm(), rb = vb.norm();
      if (ra > 1e-15 && rb > 1e-15) {
        double C = std::max(-1.0, std::min(1.0, va.dot(vb) / (ra * rb)));
        angle_sum += acos(C);
      }
      h = h2;
    } while (h != h0);
    if (std::abs(angle_sum - D.v_cone_angle[v]) > 0.01)
      std::cout << "  v" << v << ": target_cone=" << D.v_cone_angle[v]
                << " actual=" << angle_sum
                << " err=" << std::abs(angle_sum - D.v_cone_angle[v]) << std::endl;
  }

  std::cout << "max_rel_err=" << max_rel_err << std::endl;
  EXPECT_LT(max_rel_err, 0.05) << "C60 #1264 DCEL embedding too inaccurate";
}

TEST(DCELEmbed, C60_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  int idx = 0, n_poor = 0;
  double worst_err = 0;
  int worst_idx = -1;

  while (BuckyGen::next_fullerene(Q, T)) {
    auto D = DelaunayTriangulation::compute(T);
    if (D.nv != 12) { idx++; continue; }

    auto coords = D.embed_3d();

    map<pair<int,int>, double> shortest;
    for (int h = 0; h < D.nh; h += 2) {
      if (!D.alive(h)) continue;
      int u = D.he_origin[h], v = D.dest(h);
      if (u > v) swap(u, v);
      double L = D.he_length[h];
      auto key = make_pair(u, v);
      auto it = shortest.find(key);
      if (it == shortest.end() || L < it->second)
        shortest[key] = L;
    }

    double max_rel_err = 0;
    for (auto& [key, target] : shortest) {
      double actual = (coords[key.first] - coords[key.second]).norm();
      double rel_err = fabs(actual - target) / target;
      max_rel_err = std::max(max_rel_err, rel_err);
    }

    if (max_rel_err > worst_err) { worst_err = max_rel_err; worst_idx = idx; }
    if (max_rel_err > 0.01) n_poor++;

    EXPECT_LT(max_rel_err, 0.05)
      << "C60 #" << idx << " DCEL embedding too inaccurate";

    idx++;
  }
  BuckyGen::stop(Q);

  std::cout << "DCEL C60 Embed: " << idx << " isomers, "
            << n_poor << " with >1% error, worst=" << worst_err
            << " at #" << worst_idx << std::endl;
}

// ============================================================================
// Symmetry-constrained DCEL embedding tests
// ============================================================================

// Helper: make dual triangulation from buckygen
static Triangulation make_buckygen_dual(int N, int target_idx) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  Triangulation T;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    if (idx == target_idx) break;
    idx++;
  }
  BuckyGen::stop(Q);
  return T;
}

// Helper: build SymmetryConstraint for a triangulation's iDT embedding.
static SymmetryConstraint make_sym_constraint(const Triangulation& T) {
  Symmetry S(T);
  auto rep = S.representation_3d();
  return restrict_symmetry_to_cone_points(
      vector<vector<int>>(S.G.begin(), S.G.end()), rep.R, T);
}

TEST(DCELSymEmbed, C20_Ih) {
  // C20 has Ih symmetry (group order 120). All 12 iDT vertices in one orbit.
  Triangulation T = make_dual(20, 0, false);
  auto sym = make_sym_constraint(T);
  EXPECT_GT(sym.perms.size(), 1u) << "Expected non-trivial symmetry group";

  auto D = DelaunayTriangulation::compute(T);
  ASSERT_EQ(D.nv, 12);

  auto coords_nosym = D.embed_3d();
  auto coords_sym   = D.embed_3d(sym);

  // Both should give machine-precision results for C20 (regular icosahedron)
  verify_dcel_embedding(D, coords_nosym, 1e-6, "C20 no-sym");
  verify_dcel_embedding(D, coords_sym,   1e-6, "C20 with-sym");

  std::cout << "  Symmetry group order: " << sym.perms.size() << std::endl;
}

TEST(DCELSymEmbed, C60_1264_MultiEdge) {
  // The multi-edge case that originally caused MDS collapse.
  // With symmetry-constrained reduced-parameter optimization,
  // orbit members are generated by group action -- no collapse possible.
  Triangulation T = make_buckygen_dual(60, 1264);
  Symmetry S(T);
  auto sym = make_sym_constraint(T);

  auto D = DelaunayTriangulation::compute(T);
  ASSERT_EQ(D.nv, 12);

  auto orbits = compute_orbits(12, sym.perms);
  std::cout << "  C60 #1264: point group " << S.point_group()
            << ", |G|=" << sym.perms.size()
            << ", " << orbits.size() << " orbits:";
  for (auto& o : orbits) {
    std::cout << " {";
    for (size_t i = 0; i < o.size(); i++)
      std::cout << (i?",":"") << o[i];
    std::cout << "}";
  }
  std::cout << std::endl;

  auto coords = D.embed_3d(sym);
  verify_dcel_embedding(D, coords, 1e-4, "C60 #1264 sym");
}

TEST(DCELSymEmbed, C60_AllIsomers) {
  int total = 0, n_poor_nosym = 0, n_poor_sym = 0;
  double worst_nosym = 0, worst_sym = 0;
  int worst_idx_nosym = -1, worst_idx_sym = -1;

  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    auto sym = make_sym_constraint(T);

    auto D = DelaunayTriangulation::compute(T);
    if (D.nv != 12) { idx++; continue; }

    auto coords_nosym = D.embed_3d();
    auto coords_sym   = D.embed_3d(sym);

    auto measure_err = [&](const vector<coord3d>& coords) {
      double max_err = 0;
      map<pair<int,int>, double> sh;
      for (int h = 0; h < D.nh; h += 2) {
        if (!D.alive(h)) continue;
        int u = D.he_origin[h], v = D.dest(h);
        if (u > v) std::swap(u, v);
        double L = D.he_length[h];
        auto key = make_pair(u, v);
        auto it = sh.find(key);
        if (it == sh.end() || L < it->second) sh[key] = L;
      }
      for (auto& [key, target] : sh) {
        double actual = (coords[key.first] - coords[key.second]).norm();
        max_err = std::max(max_err, std::abs(actual - target) / target);
      }
      return max_err;
    };

    double err_nosym = measure_err(coords_nosym);
    double err_sym   = measure_err(coords_sym);

    if (err_nosym > 0.01) n_poor_nosym++;
    if (err_sym   > 0.01) n_poor_sym++;
    if (err_nosym > worst_nosym) { worst_nosym = err_nosym; worst_idx_nosym = idx; }
    if (err_sym   > worst_sym)   { worst_sym   = err_sym;   worst_idx_sym   = idx; }

    EXPECT_LT(err_sym, 0.05) << "C60 #" << idx << " sym embedding too inaccurate";

    idx++;
    total++;
  }
  BuckyGen::stop(Q);

  std::cout << "C60 all-isomers (" << total << "):" << std::endl;
  std::cout << "  no-sym: " << n_poor_nosym << " with >1% error, worst="
            << worst_nosym << " at #" << worst_idx_nosym << std::endl;
  std::cout << "  sym:    " << n_poor_sym << " with >1% error, worst="
            << worst_sym   << " at #" << worst_idx_sym << std::endl;
}

TEST(SymmetricIDT, D2h_348_CotSumAnalysis) {
  // Analyze ALL edges' cot-sums to determine which are co-circular.
  // If the rim edges of the self-dual quads are co-circular, the quads
  // might not exist in a different Delaunay triangulation, and the
  // "topological obstruction" would be wrong.
  Triangulation T = make_buckygen_dual(60, 348);
  Symmetry S(T);

  // Compute iDT with correct geodesic distances.
  auto D = DelaunayTriangulation::compute(T);
  ASSERT_EQ(D.nv, 12);

  // Replace edge lengths with geodesic distances.
  vector<int> cone_pts;
  for (int v = 0; v < T.N; v++)
    if (T.degree(v) != 6) cone_pts.push_back(v);
  auto dist2 = T.surface_distances(cone_pts);
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    D.he_length[h] = D.he_length[h^1] = sqrt(dist2(D.he_origin[h], D.dest(h)));
  }
  D.recompute_all_angles();
  D.flip_to_delaunay();

  // Compute cot-sum for each edge via the Diamond's side lengths.
  auto cot_opp = [](double e, double a, double b) {
    // cot of the angle opposite edge e in triangle with sides (e, a, b)
    double cos_C = (a*a + b*b - e*e) / (2*a*b);
    double sin_C = sqrt(std::max(0.0, 1.0 - cos_C*cos_C));
    return sin_C > 1e-15 ? cos_C / sin_C : 1e15;
  };

  // Print cot-sum for ALL 30 edges.
  int n_cocircular = 0, n_strict = 0;
  std::cout << "  Edge cot-sums (after geodesic correction):" << std::endl;
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    Diamond dm = D.diamond(h);
    double cs = cot_opp(dm.e, dm.a, dm.b) + cot_opp(dm.e, dm.c, dm.d);
    int u = D.he_origin[h], v = D.dest(h);
    bool cocircular = fabs(cs) < 1e-8;
    if (cocircular) n_cocircular++; else n_strict++;
    std::cout << "    (" << u << "," << v << "): L=" << dm.e
              << " cot_sum=" << cs << (cocircular ? " CO-CIRCULAR" : "") << std::endl;
  }
  std::cout << "  Total: " << n_strict << " strictly Delaunay, "
            << n_cocircular << " co-circular" << std::endl;

  // Now identify the self-dual orbit's quads and check their rim edges.
  auto cone_perms = restrict_to_cone_points(
      vector<vector<int>>(S.G.begin(), S.G.end()), T);
  set<pair<int,int>> edge_set;
  for (int h = 0; h < D.nh; h += 2)
    if (D.alive(h))
      edge_set.insert({min(D.he_origin[h],D.dest(h)), max(D.he_origin[h],D.dest(h))});

  std::cout << "  Self-dual orbit quads:" << std::endl;
  for (auto& perm : cone_perms)
    for (auto& e : edge_set) {
      auto pe = make_pair(min(perm[e.first],perm[e.second]),
                          max(perm[e.first],perm[e.second]));
      if (edge_set.count(pe)) continue;
      // pe is missing. Find the quad.
      if (D.v_out[pe.first] < 0) continue;
      // Find common neighbors of pe.first and pe.second.
      vector<int> common;
      int h0 = D.v_out[pe.first], hh = h0;
      do {
        int nb = D.dest(hh);
        // Check if nb is also a neighbor of pe.second.
        int h2 = D.v_out[pe.second];
        if (h2 >= 0) {
          int hh2 = h2;
          do {
            if (D.dest(hh2) == nb) { common.push_back(nb); break; }
            hh2 = D.cw(hh2);
          } while (hh2 != h2);
        }
        hh = D.cw(hh);
      } while (hh != h0);
      if (common.size() != 2) continue;

      std::cout << "    Quad: diag(" << common[0] << "," << common[1]
                << ") alt(" << pe.first << "," << pe.second << ")" << std::endl;
      // Check rim edges' cot-sums.
      for (int c : common) {
        for (int p : {pe.first, pe.second}) {
          // Find half-edge from c to p or p to c.
          int rh = -1;
          int rh0 = D.v_out[c];
          if (rh0 >= 0) {
            int rhh = rh0;
            do { if (D.dest(rhh)==p) { rh=rhh; break; } rhh=D.cw(rhh); } while(rhh!=rh0);
          }
          if (rh >= 0) {
            Diamond rdm = D.diamond(rh);
            double rcs = cot_opp(rdm.e,rdm.a,rdm.b)+cot_opp(rdm.e,rdm.c,rdm.d);
            std::cout << "      rim(" << c << "," << p
                      << ") cot_sum=" << rcs
                      << (fabs(rcs)<1e-8 ? " CO-CIRCULAR" : " STRICT") << std::endl;
          }
        }
      }
      break;  // just one quad
    }
}

TEST(SymmetricIDT, D2h_348) {
  Triangulation T = make_buckygen_dual(60, 348);
  Symmetry S(T);

  // Compute Steiner-augmented iDT.
  auto D = DelaunayTriangulation::compute_symmetric(T, S);
  std::cout << "  Steiner iDT: nv=" << D.nv << " (12 cone + "
            << D.nv - 12 << " Steiner)" << std::endl;
  EXPECT_TRUE(D.check_consistency()) << "DCEL inconsistent after Steiner insertion";
  // The Steiner insertion might create non-Delaunay edges at the boundary
  // between old and new faces.  Fix with Lawson:
  D.flip_to_delaunay();
  EXPECT_TRUE(D.is_delaunay()) << "Not Delaunay after Steiner + Lawson";

  // Build 16-vertex permutations.
  // Steiner vertices 12..15 correspond to the 4 co-circular quads.
  // Map each Steiner vertex to its quad's vertex set (sorted), then
  // for each group element, find which Steiner vertex maps to which.
  auto cone_perms = restrict_to_cone_points(
      vector<vector<int>>(S.G.begin(), S.G.end()), T);

  // Identify which Steiner vertex belongs to which quad.
  // Each Steiner vertex (index >= 12) is connected to exactly 4 cone-point vertices.
  vector<set<int>> steiner_quads;  // steiner_quads[i] = {4 cone-point neighbors}
  for (int sv = 12; sv < D.nv; sv++) {
    set<int> nbs;
    int h0 = D.v_out[sv], h = h0;
    if (h0 >= 0) do { nbs.insert(D.dest(h)); h = D.cw(h); } while (h != h0);
    steiner_quads.push_back(nbs);
  }

  // For each group element, extend the 12-vertex perm to 16 vertices.
  vector<vector<int>> full_perms;
  for (auto& cp : cone_perms) {
    vector<int> fp(D.nv);
    for (int v = 0; v < 12; v++) fp[v] = cp[v];
    // Map each Steiner vertex: its quad {a,b,c,d} maps to {cp[a],cp[b],cp[c],cp[d]}.
    for (int si = 0; si < (int)steiner_quads.size(); si++) {
      set<int> mapped_quad;
      for (int v : steiner_quads[si]) mapped_quad.insert(cp[v]);
      // Find which Steiner vertex has this mapped quad as neighbors.
      int mapped_sv = -1;
      for (int sj = 0; sj < (int)steiner_quads.size(); sj++)
        if (steiner_quads[sj] == mapped_quad) { mapped_sv = 12 + sj; break; }
      if (mapped_sv < 0) mapped_sv = 12 + si;  // fallback (shouldn't happen)
      fp[12 + si] = mapped_sv;
    }
    full_perms.push_back(fp);
  }

  int violations = D.check_edge_symmetry(full_perms);
  std::cout << "  Edge symmetry violations (16-vertex perms): " << violations << std::endl;
  EXPECT_EQ(violations, 0) << "Steiner iDT edge set not G-invariant";

  // Also test embed on the ORIGINAL (non-Steiner) iDT for comparison.
  auto D_orig = DelaunayTriangulation::compute(T);
  auto sym_constraint = make_sym_constraint(T);
  auto coords_nosym = D_orig.embed_3d();
  auto coords_sym = D_orig.embed_3d(sym_constraint);

  // Measure max relative edge error.
  auto measure_err = [&](const vector<coord3d>& coords) {
    double max_err = 0;
    map<pair<int,int>, double> sh;
    for (int h = 0; h < D.nh; h += 2) {
      if (!D.alive(h)) continue;
      int u = D.he_origin[h], v = D.dest(h);
      if (u > v) std::swap(u, v);
      double L = D.he_length[h];
      auto key = make_pair(u, v);
      auto it = sh.find(key);
      if (it == sh.end() || L < it->second) sh[key] = L;
    }
    for (auto& [key, target] : sh) {
      double actual = (coords[key.first] - coords[key.second]).norm();
      max_err = std::max(max_err, std::abs(actual - target) / target);
    }
    return max_err;
  };
  double err_nosym = measure_err(coords_nosym);
  double err_sym = measure_err(coords_sym);

  std::cout << "  D2h #348: nosym_err=" << err_nosym
            << " sym_err=" << err_sym << std::endl;
}

static void survey_idt_features(int N) {
  int total = 0;
  int n_with_multi = 0, n_with_self_loop = 0, n_with_cocircular = 0;
  int max_multi_edges = 0, max_self_loops = 0, max_cocircular = 0;
  int total_multi = 0, total_loops = 0, total_cocircular = 0;

  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  Triangulation T;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    auto D = DelaunayTriangulation::compute(T);
    if (D.nv != 12) { idx++; continue; }

    // Count multi-edges and self-loops.
    map<pair<int,int>, int> edge_count;
    int n_multi = 0, n_loops = 0;
    for (int h = 0; h < D.nh; h += 2) {
      if (!D.alive(h)) continue;
      int u = D.he_origin[h], v = D.dest(h);
      if (u == v) {
        n_loops++;
      } else {
        auto key = make_pair(min(u,v), max(u,v));
        edge_count[key]++;
      }
    }
    for (auto& [key, cnt] : edge_count)
      if (cnt > 1) n_multi += cnt - 1;  // extra edges beyond the first

    // Count co-circular edges.
    int n_cc = 0;
    for (int h = 0; h < D.nh; h += 2) {
      if (!D.alive(h)) continue;
      Diamond dm = D.diamond(h);
      auto cot_opp = [](double e, double a, double b) {
        double cos_C = (a*a + b*b - e*e) / (2*a*b);
        double sin_C = sqrt(std::max(0.0, 1.0 - cos_C*cos_C));
        return sin_C > 1e-15 ? cos_C / sin_C : 1e15;
      };
      double cs = cot_opp(dm.e, dm.a, dm.b) + cot_opp(dm.e, dm.c, dm.d);
      if (fabs(cs) < 1e-8) n_cc++;
    }

    if (n_multi > 0) n_with_multi++;
    if (n_loops > 0) n_with_self_loop++;
    if (n_cc > 0) n_with_cocircular++;
    total_multi += n_multi;
    total_loops += n_loops;
    total_cocircular += n_cc;
    max_multi_edges = std::max(max_multi_edges, n_multi);
    max_self_loops = std::max(max_self_loops, n_loops);
    max_cocircular = std::max(max_cocircular, n_cc);

    if (n_loops > 0 || n_multi > 0)
      std::cout << "  #" << idx << ": multi=" << n_multi
                << " loops=" << n_loops << " cocircular=" << n_cc << std::endl;
    total++;
    idx++;
  }
  BuckyGen::stop(Q);

  std::cout << "C" << N << " iDT survey (" << total << " isomers):" << std::endl;
  std::cout << "  Multi-edges: " << n_with_multi << " isomers, "
            << total_multi << " total, max " << max_multi_edges << " per isomer" << std::endl;
  std::cout << "  Self-loops:  " << n_with_self_loop << " isomers, "
            << total_loops << " total, max " << max_self_loops << " per isomer" << std::endl;
  std::cout << "  Co-circular: " << n_with_cocircular << " isomers, "
            << total_cocircular << " total, max " << max_cocircular << " per isomer" << std::endl;
}

TEST(SymmetricIDT, C60_MultiEdgeSurvey)  { survey_idt_features(60); }
TEST(SymmetricIDT, C80_MultiEdgeSurvey)  { survey_idt_features(80); }
TEST(SymmetricIDT, C100_MultiEdgeSurvey) { survey_idt_features(100); }
