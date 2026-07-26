#include <gtest/gtest.h>
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_embed3d.hh"
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

// Round-trip the .idt text serializer (to_ascii / from_ascii) on a genuine multi-edge
// delta-complex, through a temp file. The reconstruction must preserve the metric + topology and
// pass check_consistency (from_ascii throws otherwise), and its v_cone_angle cache must be coherent.
TEST(DCEL, IdtRoundTrip) {
  // C60 isomer #1264: DelaunayTriangulation::compute yields a 12-cone delta-complex WITH a
  // multi-edge (two geodesics between one cone pair) -- exercises the non-simplicial path.
  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) { if (idx == 1264) break; idx++; }
  BuckyGen::stop(Q);
  DelaunayTriangulation D = DelaunayTriangulation::compute(T);
  ASSERT_EQ(D.nv, 12);
  ASSERT_TRUE(D.check_consistency());
  ASSERT_FALSE(D.is_simplicial()) << "expected a multi-edge delta-complex for C60 #1264";

  FILE* f = std::tmpfile();
  ASSERT_NE(f, nullptr);
  ASSERT_TRUE(DelaunayTriangulation::to_ascii(D, f));
  std::rewind(f);
  DelaunayTriangulation D2 = DelaunayTriangulation::from_ascii(f);   // validates via check_consistency
  std::fclose(f);

  EXPECT_EQ(D2.nv, D.nv);
  EXPECT_FALSE(D2.is_simplicial()) << "the multi-edge must survive the round-trip";
  EXPECT_TRUE(D2.check_consistency());

  auto sorted_lengths = [](const DelaunayTriangulation& G) {
    std::vector<double> L;
    for (int h = 0; h < G.nh; h += 2) if (G.alive(h)) L.push_back(G.he_length[h]);
    std::sort(L.begin(), L.end());
    return L;
  };
  auto sorted_defects = [](const DelaunayTriangulation& G) {
    std::vector<double> K;
    for (int v = 0; v < G.nv; v++) K.push_back(2 * M_PI - G.vertex_angle_sum(v));
    std::sort(K.begin(), K.end());
    return K;
  };
  EXPECT_EQ(sorted_lengths(D2), sorted_lengths(D));   // exact double round-trip (%.17g)
  std::vector<double> Ka = sorted_defects(D), Kb = sorted_defects(D2);
  ASSERT_EQ(Ka.size(), Kb.size());
  for (size_t i = 0; i < Ka.size(); i++) EXPECT_NEAR(Ka[i], Kb[i], 1e-12);

  // v_cone_angle cache coherence: curvature() (reads the cache) == 2pi - vertex_angle_sum (direct).
  std::vector<double> Kcache = D2.curvature();
  for (int v = 0; v < D2.nv; v++)
    EXPECT_NEAR(Kcache[v], 2 * M_PI - D2.vertex_angle_sum(v), 1e-12);
}

// A malformed .idt record is rejected, not loaded as a silently-wrong triangulation.
TEST(DCEL, IdtRejectsMalformed) {
  // An out-of-range origin vertex (o0=9 with nv=3) trips from_ascii's per-field range guard.
  FILE* f = std::tmpfile();
  ASSERT_NE(f, nullptr);
  std::fprintf(f, "iDT-DCEL 1\n3 1 3\n3\n3\n3\n9 1 2 5 1.0\n1 2 4 1 1.0\n2 0 0 3 1.0\n");
  std::rewind(f);
  EXPECT_THROW(DelaunayTriangulation::from_ascii(f), std::runtime_error);
  std::fclose(f);
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

// Origin tracking tests (ExactOrigins_C20, C60_Ih, C60_AllIsomers,
// DISABLED_ExactOrigins_C80_AllIsomers, DISABLED_ExactOrigins_C100_AllIsomers)
// were removed when the f_origin / OriginTracker machinery was archived.
// See src/c++/attic/delaunay_origin_tracking.cc.attic for the original code.

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

  auto coords = embed_delaunay_3d(D);
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

  auto coords = embed_delaunay_3d(D);
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

      auto coords = embed_delaunay_3d(D);

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

  // Bisect multi-edges before embedding.
  int n_bisected = D.bisect_multi_edges();
  if (n_bisected > 0) {
    std::cout << "Bisected " << n_bisected << " multi-edges, nv=" << D.nv << std::endl;
    EXPECT_TRUE(D.check_consistency()) << "Inconsistent after bisection";
  }

  auto coords = embed_delaunay_3d(D);

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
    D.bisect_multi_edges();

    auto coords = embed_delaunay_3d(D);

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

static void test_embed_all(int N, double threshold = 0.05) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  Triangulation T;
  int idx = 0, n_poor = 0, total = 0;
  double worst_err = 0;
  int worst_idx = -1;
  vector<double> all_errors;

  while (BuckyGen::next_fullerene(Q, T)) {
    auto D = DelaunayTriangulation::compute(T);
    if (D.nv != 12) { idx++; continue; }
    D.bisect_multi_edges();

    auto coords = embed_delaunay_3d(D);

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
      max_rel_err = std::max(max_rel_err, fabs(actual - target) / target);
    }

    if (max_rel_err > worst_err) { worst_err = max_rel_err; worst_idx = idx; }
    if (max_rel_err > 0.01) n_poor++;
    if (max_rel_err > 0.25) {
      // Detailed diagnosis of complete failures.
      int n_multi = 0;
      map<pair<int,int>, int> edge_count;
      for (int h = 0; h < D.nh; h += 2) {
        if (!D.alive(h)) continue;
        int eu = D.he_origin[h], ev = D.dest(h);
        edge_count[{min(eu,ev), max(eu,ev)}]++;
      }
      for (auto& [k, c] : edge_count) if (c > 1) n_multi += c - 1;

      // Get RSPI (pentagon positions in the spiral).
      vector<int> spiral;
      jumplist_t jumps;
      T.get_spiral(spiral, jumps);
      std::cout << "  C" << N << " #" << idx << ": err=" << max_rel_err
                << " nv=" << D.nv << " multi=" << n_multi
                << " RSPI=[";
      bool first = true;
      for (size_t si = 0; si < spiral.size(); si++)
        if (spiral[si] == 5) { if (!first) std::cout << ","; std::cout << si+1; first = false; }
      std::cout << "]";
      if (!jumps.empty()) {
        std::cout << " jumps=[";
        first = true;
        for (auto& j : jumps) { if (!first) std::cout << ","; std::cout << "(" << j.first << "," << j.second << ")"; first = false; }
        std::cout << "]";
      }
      std::cout << std::endl;
    } else if (max_rel_err > 0.05) {
      std::cout << "  C" << N << " #" << idx << ": err=" << max_rel_err
                << " nv=" << D.nv << std::endl;
    }
    all_errors.push_back(max_rel_err);
    total++;
    idx++;
  }
  BuckyGen::stop(Q);

  sort(all_errors.begin(), all_errors.end());
  int n = all_errors.size();
  std::cout << "C" << N << " Embed: " << total << " isomers, "
            << n_poor << " with >1% error, worst=" << worst_err
            << " at #" << worst_idx << std::endl;
  if (n > 0) {
    std::cout << "  Percentiles:";
    for (double p : {90.0, 95.0, 99.0, 99.9, 99.99, 100.0}) {
      int i = std::min((int)(p/100.0 * n), n-1);
      std::cout << " p" << p << "=" << all_errors[i];
    }
    std::cout << std::endl;
  }
  EXPECT_EQ(n_poor, 0) << "Some isomers have >1% embedding error";
}

TEST(DCELEmbed, C70_AllIsomers)  { test_embed_all(70); }
TEST(DCELEmbed, C80_AllIsomers)  { test_embed_all(80); }
TEST(DCELEmbed, C100_AllIsomers) { test_embed_all(100); }

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

  auto coords_nosym = embed_delaunay_3d(D);
  auto coords_sym   = embed_delaunay_3d(D, sym);

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

  auto coords = embed_delaunay_3d(D, sym);
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

    // Bisect multi-edges for nosym path only. The sym path uses
    // full-space optimization from symmetric init and doesn't need bisection.
    auto D_bisected = D;
    D_bisected.bisect_multi_edges();

    auto coords_nosym = embed_delaunay_3d(D_bisected);
    auto coords_sym   = embed_delaunay_3d(D, sym);

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


// ============================================================================
// DCELOwnership: the owner's span/storage aliasing contract.
//
// Each test states one claim about the Owned-view pattern (delaunay.hh
// banner): copies are independent under in-place mutation, every owned field
// survives a copy, moves steal storage and leave the empty complex, growth
// re-binds every span.  These are the surfaces the compute() byte gates
// cannot see (reduction never grows, and no gate copies-then-mutates-both).
// ============================================================================

// The equality relation on owned DCEL state: counts, the nine arrays (folded
// over the canonical field tuple, so a field added in stage 3 is compared
// automatically), free lists, and the tracker's points.
static bool same_dcel_state(const DelaunayTriangulation& A,
                            const DelaunayTriangulation& B) {
  if (A.nv != B.nv || A.nh != B.nh || A.nf != B.nf) return false;
  bool eq = true;
  auto ta = A.to_tuple(), tb = B.to_tuple();
  auto spans_equal = [](auto sa, auto sb) {
    return sa.size() == sb.size() && std::equal(sa.begin(), sa.end(), sb.begin());
  };
  [&]<std::size_t... I>(std::index_sequence<I...>) {
    ((eq = eq && spans_equal(std::get<I>(ta), std::get<I>(tb))), ...);
  }(std::make_index_sequence<DelaunayView::n_fields>{});
  if (!eq) return false;
  if (A.free_edges != B.free_edges || A.free_faces != B.free_faces) return false;
  if (A.tracker.points.size() != B.tracker.points.size()) return false;
  for (size_t i = 0; i < A.tracker.points.size(); i++) {
    const auto& p = A.tracker.points[i];
    const auto& q = B.tracker.points[i];
    if (p.label != q.label || p.face != q.face ||
        p.b[0] != q.b[0] || p.b[1] != q.b[1] || p.b[2] != q.b[2]) return false;
  }
  return true;
}

// The repoint() postcondition (@inv sized) as a named predicate.
static bool array_sizes_match_counts(const DelaunayTriangulation& D) {
  return (int)D.he_next.size()   == D.nh && (int)D.he_origin.size()     == D.nh
      && (int)D.he_face.size()   == D.nh && (int)D.he_length.size()     == D.nh
      && (int)D.he_angle.size()  == D.nh
      && (int)D.v_out.size()     >= D.nv && (int)D.v_cone_angle.size()  >= D.nv
      && (int)D.v_orig_degree.size() >= D.nv
      && (int)D.f_he.size()      == D.nf;
}

// A copy must be independent under in-place mutation.  Flips are the ONE
// mutation shape that can expose span aliasing: they rewire without ever
// allocating, so a copy whose spans alias the source keeps aliasing it for
// the whole operation (growth would self-heal via repoint).  This is exactly
// the Alexandrov trial-step shape (T_trial = T; flip; accept/reject).
TEST(DCELOwnership, CopyIsIndependentUnderInPlaceMutation) {
  Triangulation T = make_dual(60, 0);
  auto D  = DelaunayTriangulation::compute(T);
  auto D2 = D;
  ASSERT_TRUE(same_dcel_state(D, D2));
  ASSERT_NE(D.he_next.data(), D2.he_next.data());   // deep copy, repointed

  int flipped = -1;
  for (int h : D2.edges())
    if (D2.flip_edge(h)) { flipped = h; break; }
  ASSERT_GE(flipped, 0) << "no flippable edge on the C60 iDT";

  auto D_fresh = DelaunayTriangulation::compute(T);  // compute is deterministic
  EXPECT_TRUE(same_dcel_state(D, D_fresh)) << "mutating the copy changed the source";
  EXPECT_FALSE(same_dcel_state(D2, D)) << "the flip did not change the copy";
}

// Every owned field -- arrays, free lists, tracker points -- survives a copy,
// on a complex where all of them are populated (metric compute with
// track_removed: dead slots, recycled free-list entries, tracked points).
// Pins the documented "copies snapshot the tracker" postcondition.
TEST(DCELOwnership, CopyPreservesEveryField) {
  Triangulation T = make_dual(60, 0);
  auto unit = [](node_t, node_t) { return 1.0; };
  std::vector<int> n2o;
  auto D = DelaunayTriangulation::compute(T, unit, 1e-6, &n2o, /*track_removed=*/true);
  ASSERT_FALSE(D.free_edges.empty());               // reduction recycled slots
  ASSERT_FALSE(D.tracker.points.empty());           // removals were tracked

  auto D2 = D;
  EXPECT_TRUE(same_dcel_state(D, D2));
  EXPECT_NE(D.he_next.data(), D2.he_next.data());
  EXPECT_TRUE(array_sizes_match_counts(D2));
}

// A move steals the storage (same data pointers) and leaves the source as
// the documented empty complex: counts zero, spans re-bound (never dangling).
TEST(DCELOwnership, MoveStealsStorageAndEmptiesSource) {
  auto D = DelaunayTriangulation::compute(make_dual(60, 0));
  const int* p = D.he_next.data();
  auto D2 = std::move(D);
  EXPECT_EQ(D2.he_next.data(), p);                  // storage stolen, repointed
  EXPECT_TRUE(D2.check_consistency());
  EXPECT_EQ(D.nv, 0);
  EXPECT_EQ(D.nh, 0);
  EXPECT_EQ(D.nf, 0);
  EXPECT_TRUE(D.he_next.empty());                   // re-bound to emptied vector
}

TEST(DCELOwnership, SelfAssignmentIsIdentity) {
  auto D  = DelaunayTriangulation::compute(make_dual(60, 0));
  auto D2 = D;
  DelaunayTriangulation& alias = D;
  D = alias;                                        // self copy-assign
  EXPECT_TRUE(same_dcel_state(D, D2));
  D = std::move(alias);                             // self move-assign
  EXPECT_TRUE(same_dcel_state(D, D2));
  EXPECT_TRUE(D.check_consistency());
}

// Growth must re-bind every span.  A FRESH build has exact-size arrays and
// EMPTY free lists, so face-splitting is forced through all three ensure_*
// helpers (vertex, edge, and face growth).  Note a post-compute complex
// cannot exercise this: reduction stocks the free lists, and subsequent
// surgery allocates from them without growing.
TEST(DCELOwnership, GrowthRepointsEveryArray) {
  Triangulation T = make_dual(20, 0);
  auto D = DelaunayTriangulation::from_triangulation(T);
  ASSERT_TRUE(D.free_edges.empty());                // fresh build: nothing recycled
  const int nv0 = D.nv, nh0 = D.nh, nf0 = D.nf;

  int P = D.split_face(D.f_he[0], {1.0, 1.0, 1.0});
  EXPECT_EQ(P, nv0);                                // the appended vertex
  EXPECT_EQ(D.nv, nv0 + 1);
  EXPECT_EQ(D.nh, nh0 + 6);                         // three new edges
  EXPECT_EQ(D.nf, nf0 + 2);                         // 1 face -> 3 (one slot recycled)
  EXPECT_TRUE(array_sizes_match_counts(D));         // every span re-bound
  EXPECT_TRUE(D.check_consistency());

  // And the recycled-slot surgery path stays consistent on a reduced complex
  // (bisect allocates from the free lists; no growth, byte-gated elsewhere).
  // Multi-edge isomers are found by predicate, not by index: enumeration
  // orders differ between IsomerDB and BuckyGen, so a hardcoded index names
  // different isomers in different tests.
  bool found_multi = false;
  for (int idx = 0; idx < 1812; idx++) {
    auto R = DelaunayTriangulation::compute(make_dual(60, idx));
    if (R.is_simplicial()) continue;
    found_multi = true;
    ASSERT_FALSE(R.free_edges.empty());
    EXPECT_GT(R.bisect_multi_edges(), 0);
    EXPECT_TRUE(array_sizes_match_counts(R));
    EXPECT_TRUE(R.check_consistency());
    break;
  }
  ASSERT_TRUE(found_multi) << "no multi-edge C60 iDT found";
}

// The capacity formulas are BUILD-TIME quantities: dcel_capacities takes the
// pre-reduction vertex count nv0, never the live nv of a reduced complex.
TEST(DCELOwnership, CapacitiesAreBuildTime) {
  Triangulation T = make_dual(60, 0);
  const DcelCapacities cap = dcel_capacities(T.N);
  auto B = DelaunayTriangulation::from_triangulation(T);
  EXPECT_EQ(cap.nh_cap, (long)B.he_next.size());
  EXPECT_EQ(cap.nf_cap, (long)B.f_he.size());

  auto D = DelaunayTriangulation::compute(T);
  EXPECT_EQ(D.nv, 12);                              // reduced to the cones
  EXPECT_EQ((long)D.he_next.size(), cap.nh_cap);    // arrays keep build size
  EXPECT_NE(cap.nh_cap, dcel_capacities(D.nv).nh_cap) << "live nv is the wrong argument";
}
