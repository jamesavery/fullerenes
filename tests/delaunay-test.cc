#include <gtest/gtest.h>
#include "fullerenes/delaunay.hh"
#include "delaunay-counting-metric.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/planargraph.hh"
#include "fullerenes/buckygen-wrapper.hh"
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
// DCELExact: the exact-integer predicate layer (paper sec:exactness).
//
// On lattice triangles H = 3*tau^2 (the Heron-Eisenstein identity), so the
// Delaunay form F is sqrt(3) times the INTEGER s_upper*tau_lower +
// s_lower*tau_upper.  The falsifiers pin lattice_tau, the form's argument
// pairing, and the named diamonds; the corpus sweeps assert the exact
// predicates agree with the banded float predicates on every diamond both
// BEFORE reduction (where non-Delaunay edges exist) and after (a
// disagreement is not a test bug, it is the documented sliver and must be
// examined); the trip falsifiers pin every new exact failure mode loudly.
// ============================================================================

TEST(DCELExact, LatticeTauCases) {
  // The unit triangle: H = 3, tau = 1.
  EXPECT_EQ(lattice_tau(1, 1, 1), 1);
  // The (1,1,3) triangle (unit sides, sqrt(3) base): H = 3, tau = 1.
  EXPECT_EQ(lattice_tau(3, 1, 1), 1);
  // A larger lattice triangle: sides^2 (4, 4, 4) -> H = 48 = 3*16, tau = 4.
  EXPECT_EQ(lattice_tau(4, 4, 4), 4);
  // Degenerate: sides 2, 1, 1 (squared 4, 1, 1) are collinear -> tau = 0.
  EXPECT_EQ(lattice_tau(4, 1, 1), 0);
  // Triangle inequality violated: H < 0 -> not a triangle at all.
  EXPECT_EQ(lattice_tau(9, 1, 1), -1);
  // A right isoceles triangle (squared sides 2, 1, 1): H = 4, not divisible
  // by 3 -- a Euclidean triangle that fits NO lattice (so the unit-square
  // diamond below is outside the exact regime's domain).
  EXPECT_EQ(lattice_tau(2, 1, 1), -1);
  // Larger H = 3*(perfect square) cases: (7,7,7) -> H = 147 = 3*49; and
  // (5,5,5) -> H = 75 = 3*25 -- tau exists although 5 is NOT a Loeschian
  // norm: the condition is necessary, not sufficient (the doc's caveat).
  EXPECT_EQ(lattice_tau(7, 7, 7), 7);
  EXPECT_EQ(lattice_tau(5, 5, 5), 5);
  // H divisible by 3 but H/3 not a perfect square: (2,2,3) -> H = 15.
  EXPECT_EQ(lattice_tau(2, 2, 3), -1);
}

TEST(DCELExact, SplitPrimeOrbits) {
  // Norm 7 is a split prime: two sector-0 representatives, (2,1) and (1,2),
  // in DISTINCT unit orbits -- they are not lattice-equivalent, which is
  // why the fan development must scan representatives (the orbit-anchor
  // discovery).  Pin the phenomenon with the (7, 1, 3) lattice triangle:
  // its CCW apex (apex_factor(7,1,3,+1) = 2+w) lands on the lattice over
  // exactly ONE of the two base choices.  (An equilateral triangle would
  // NOT discriminate: its apex is a 60-degree rotation, lattice over any
  // base.)
  auto reps = sector0_reps_of_norm(7);
  ASSERT_EQ(reps.size(), 2u);
  int n_placed = 0;
  for (Eisenstein z : reps)
    if (place_third_eis_total({0, 0}, z, 7, 1, 3, +1)) n_placed++;
  EXPECT_EQ(n_placed, 1) << "expected exactly one orbit to admit the CCW apex";
}

TEST(DCELExact, EquilateralDiamonds) {
  // The unit rhombus (two unit equilateral triangles glued along the unit
  // diagonal): apex angles 60 degrees, cot-sum 2/sqrt(3) > 0 -- strictly
  // Delaunay and strictly convex (60 + 60 = 120 < 180 at each endpoint).
  DiamondSq unit{1, 1, 1, 1, 1};
  ASSERT_TRUE(unit.delaunay_form().has_value());
  EXPECT_GT(*unit.delaunay_form(), 0);
  EXPECT_TRUE(unit.is_delaunay());
  EXPECT_FALSE(unit.is_cocircular());
  EXPECT_TRUE(unit.is_convex());
  // The same rhombus around its LONG diagonal sqrt(3): apex angles 120,
  // cot-sum negative -- must flip (and may: 30 + 30 = 60 < 180, convex).
  DiamondSq long_diag{3, 1, 1, 1, 1};
  EXPECT_LT(*long_diag.delaunay_form(), 0);
  EXPECT_FALSE(long_diag.is_delaunay());
  EXPECT_TRUE(long_diag.is_convex());
  // ... and its exact flipped diagonal is the unit one (and vice versa: the
  // unit rhombus flips to the long diagonal).
  ASSERT_TRUE(long_diag.flipped_length_sq().has_value());
  EXPECT_EQ(*long_diag.flipped_length_sq(), 1);
  EXPECT_EQ(*unit.flipped_length_sq(), 3);
  // A genuinely cocircular LATTICE diamond: two (2, sqrt(3), 1) right
  // triangles glued along the diameter e = 2 (squared sides 4, 3, 1;
  // tau = 2).  Right apex angles, cot-sum exactly 0: F' = 0 in integers.
  DiamondSq cocirc{4, 3, 1, 3, 1};
  ASSERT_TRUE(cocirc.delaunay_form().has_value());
  EXPECT_EQ(*cocirc.delaunay_form(), 0);
  EXPECT_TRUE(cocirc.is_delaunay());
  EXPECT_TRUE(cocirc.is_cocircular());
  // The unit square around a diagonal sqrt(2) is concyclic in the
  // EUCLIDEAN sense (the float predicate says so), but its faces are
  // right-isoceles (2,1,1) triangles that fit no lattice -- it is OUTSIDE
  // the exact regime's domain, and the tau form refuses it rather than
  // answering for a diamond the regime cannot contain.
  DiamondSq square{2, 1, 1, 1, 1};
  EXPECT_FALSE(square.delaunay_form().has_value());
  EXPECT_FALSE(square.is_cocircular());
  EXPECT_FALSE(square.is_delaunay());
  EXPECT_FALSE(square.is_convex());
  EXPECT_TRUE((Diamond{std::sqrt(2.0), 1, 1, 1, 1}
                   .is_cocircular(delaunay_detail::tie_cocircular_tol)));
  // The ARGUMENT PAIRING regression guard (the one transposition the
  // symmetric diamonds above cannot catch): an asymmetric diamond, upper
  // triangle (e,a,b) = (1,1,1) vs lower (e,c,d) = (1,4,3).  cot(B) =
  // 1/sqrt(3); the lower triangle has s_lower = 4+3-1 = 6, H_lower =
  // 2(4+12+3)-(1+16+9) = 12, tau_lower = 2, so F/sqrt(3) = s_upper*
  // tau_lower + s_lower*tau_upper = 1*2 + 6*1 = 8 > 0.  A transposed
  // pairing would give 1*1 + 6*2 = 13 -- also positive, so pin the VALUE.
  DiamondSq asym{1, 1, 1, 4, 3};
  EXPECT_EQ(*asym.delaunay_form(), 8);
  // Outside the predicates' domain: both triangles degenerate (tau = 0).
  // A degenerate diamond is NONE of Delaunay / cocircular / convex --
  // matching the float path's cot_degenerate convention and the historical
  // integer cocircularity predicate.  A live complex never contains such a
  // diamond (positive angles are a class invariant), so this is the carry-
  // corruption verdict, not a geometry verdict.
  DiamondSq degenerate{4, 1, 1, 1, 1};
  EXPECT_FALSE(degenerate.delaunay_form().has_value());
  EXPECT_FALSE(degenerate.is_delaunay());
  EXPECT_FALSE(degenerate.is_cocircular());
  EXPECT_FALSE(degenerate.is_convex());
}

// The corpus verdict of the exact switch (paper sec:exactness): the exact
// pipeline must make the SAME decisions as the banded pipeline everywhere
// the tolerance slivers are empty -- and over THIS corpus (all C20-C60 +
// C80-IPR + C128#0) they are.  Integer state (topology, counts, free
// lists) must be byte-equal; the float state (he_length = sqrt of the
// exact integer vs the float flip/development formulas) may differ only at
// rounding level.  max_lsq records how much envelope headroom the corpus
// actually uses (exact_lsq_max scope-honesty).
static void expect_exact_matches_banded(const Triangulation& T, const std::string& label,
                                        long long& n_isomers, double& worst_rel,
                                        long long& max_lsq) {
  auto E = DelaunayTriangulation::compute(T);          // exact metric
  auto B = DelaunayTriangulation::compute_banded(T);   // banded float metric
  n_isomers++;

  ASSERT_EQ(E.nv, B.nv) << label;
  ASSERT_EQ(E.nh, B.nh) << label;
  ASSERT_EQ(E.nf, B.nf) << label;
  // ALL NINE fields via the canonical tuple fold (the field-list contract):
  // integer fields byte-equal -- v_orig_degree included, the exact
  // flatness predicate's own input -- and the three double fields within
  // the rounding band (same reals, two arithmetics).
  {
    auto te = E.to_tuple(), tb = B.to_tuple();
    bool ok = true;
    [&]<std::size_t... I>(std::index_sequence<I...>) {
      auto cmp = [&](auto se, auto sb) {
        ASSERT_EQ(se.size(), sb.size()) << label;
        for (std::size_t i = 0; i < se.size() && ok; i++) {
          using T = std::remove_cvref_t<decltype(se[0])>;
          if constexpr (std::is_same_v<T, double>) {
            double d = std::abs(se[i] - sb[i]);
            if (!(d <= 1e-9 * std::max(1.0, std::abs(sb[i])))) {
              ok = false;
              ADD_FAILURE() << label << " double field diverges at " << i;
            }
          } else if (se[i] != sb[i]) {
            ok = false;
            ADD_FAILURE() << label << " integer field diverges at " << i;
          }
        }
      };
      (cmp(std::get<I>(te), std::get<I>(tb)), ...);
    }(std::make_index_sequence<DelaunayView::n_fields>{});
    if (!ok) return;
  }
  auto el = E.free_edges.live(), bl = B.free_edges.live();
  ASSERT_TRUE(el.size() == bl.size() && std::equal(el.begin(), el.end(), bl.begin()))
      << label << " free_edges";
  auto ef = E.free_faces.live(), bf = B.free_faces.live();
  ASSERT_TRUE(ef.size() == bf.size() && std::equal(ef.begin(), ef.end(), bf.begin()))
      << label << " free_faces";

  // Float state: same reals computed two ways (exact sqrt-of-integer vs the
  // float flip/development formulas).  he_length must agree to rounding
  // level, and the exact value must be EXACTLY the root of an integer --
  // the discriminator that pins compute(T) actually running the exact
  // regime (the banded flipped_length drifts off sqrt-of-integer).
  for (int h = 0; h < E.nh; h++) {
    if (!E.alive(h)) continue;
    double le = E.he_length[h], lb = B.he_length[h];
    double rel = std::abs(le - lb) / lb;
    worst_rel = std::max(worst_rel, rel);
    ASSERT_LT(rel, 1e-9) << label << " he_length h=" << h;
    long long n = (long long)std::llround(le * le);
    max_lsq = std::max(max_lsq, n);
    ASSERT_EQ(le, std::sqrt((double)n)) << label
        << " exact he_length is not the root of an integer, h=" << h;
  }

  // The canonical-tesselation seed must agree; the integer Gauss-Bonnet
  // must hold on both results; and the equilateral curvature label must
  // agree with the integer cone excess (kappa = eps*pi/3) at every cone --
  // the identity the exact flatness predicate rests on.
  EXPECT_EQ(E.cocircular_edges(), B.cocircular_edges()) << label;
  EXPECT_EQ(E.total_cone_excess(), 12) << label;
  EXPECT_EQ(B.total_cone_excess(), 12) << label;
  constexpr double pi_3 = std::numbers::pi_v<double> / 3.0;
  for (int v = 0; v < E.nv; v++)
    if (E.v_out[v] >= 0)
      EXPECT_NEAR(E.DelaunayView::curvature(v), E.cone_excess(v) * pi_3, 1e-9)
          << label << " v=" << v;
}

TEST(DCELExact, PipelineMatchesBandedOnCorpus) {
  long long n_isomers = 0, max_lsq = 0;
  double worst_rel = 0;
  for (int N = 20; N <= 60; N += 2) {
    if (N == 22) continue;
    BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
    Triangulation T;
    int idx = 0;
    while (BuckyGen::next_fullerene(Q, T)) {
      expect_exact_matches_banded(T, "C" + std::to_string(N) + " #" + std::to_string(idx),
                                  n_isomers, worst_rel, max_lsq);
      if (::testing::Test::HasFatalFailure()) { BuckyGen::stop(Q); return; }
      idx++;
    }
    BuckyGen::stop(Q);
  }
  // The IPR C80s and the non-simplicial C128 (multi-edge paths).  NOTE the
  // tie-break (and so lattice_first_tie_side) fires on NO isomer of the
  // C20-C60 + C80-IPR + C128 corpus (ephemeral instrumented probe of the
  // review round: 0 calls in those 5,778 isomers; for the large-N leg
  // below the byte-equality is consistent with no-firing OR
  // firing-with-agreeing-orders).  Its dynamic coverage lives in
  // ConstructedTieExercisesTieBreak, whose article makes it fire and whose
  // A/B asserts the regimes still agree byte-for-byte through the firing.
  {
    BuckyGen::buckygen_queue Q = BuckyGen::start(80, true, false);
    Triangulation T;
    int idx = 0;
    while (BuckyGen::next_fullerene(Q, T)) {
      expect_exact_matches_banded(T, "C80-IPR #" + std::to_string(idx),
                                  n_isomers, worst_rel, max_lsq);
      if (::testing::Test::HasFatalFailure()) { BuckyGen::stop(Q); return; }
      idx++;
    }
    BuckyGen::stop(Q);
  }
  {
    BuckyGen::buckygen_queue Q = BuckyGen::start(128, false, false);
    Triangulation T;
    if (BuckyGen::next_fullerene(Q, T))
      expect_exact_matches_banded(T, "C128 #0", n_isomers, worst_rel, max_lsq);
    BuckyGen::stop(Q);
    if (::testing::Test::HasFatalFailure()) return;
  }
  // A thin large-N leg (first 50 isomers of C100 / C140 / C160): the
  // consumers re-based by the regime switch (paint, the sub-project
  // validators) run far past C60, so the A/B evidence should not stop
  // there.
  for (int N : {100, 140, 160}) {
    BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
    Triangulation T;
    int idx = 0;
    while (idx < 50 && BuckyGen::next_fullerene(Q, T)) {
      expect_exact_matches_banded(T, "C" + std::to_string(N) + " #" + std::to_string(idx),
                                  n_isomers, worst_rel, max_lsq);
      if (::testing::Test::HasFatalFailure()) { BuckyGen::stop(Q); return; }
      idx++;
    }
    BuckyGen::stop(Q);
  }
  // worst_rel > 0 certifies the two pipelines really computed the lengths
  // by different arithmetic (a 0 would mean the A/B compared one
  // computation with itself).
  EXPECT_GT(worst_rel, 0.0);
  std::cout << "DCELExact pipeline verdict: " << n_isomers
            << " isomers, integer state byte-equal, worst he_length rel dev = "
            << worst_rel << ", max Lsq = " << max_lsq
            << " (envelope " << exact_lsq_max << ")" << std::endl;
}

// Exact vs banded predicate agreement over every diamond of every C60 iDT,
// in BOTH decision directions: the iDT's own diamonds (all Delaunay by
// @post), and their FLIPPED counterparts {f^2, b, d, a, c} -- the diamond
// around the other diagonal, which is non-Delaunay wherever the original
// is strictly Delaunay, so the must-flip direction is genuinely exercised
// on real corpus geometry.  (The pre-reduction complex is all-unit --
// every diamond {1,1,1,1,1} -- and proves nothing.)  A disagreement is the
// documented tolerance sliver and must be examined, never silenced.
static void expect_predicates_agree(const DelaunayTriangulation& D, const std::string& label,
                                    long long& n_edges, long long& n_disagree) {
  auto agree = [&](const DiamondSq& ds, const Diamond& dm, int h, const char* which) {
    n_edges++;
    if (ds.is_delaunay() != dm.is_delaunay() || ds.is_convex() != dm.is_convex()) {
      if (++n_disagree <= 5)
        ADD_FAILURE() << label << " h=" << h << " (" << which
                      << "): exact/banded disagree (sliver hit)"
                      << " del " << ds.is_delaunay() << "/" << dm.is_delaunay()
                      << " cvx " << ds.is_convex() << "/" << dm.is_convex();
    }
  };
  for (int h : D.edges()) {
    Diamond dm = D.diamond(h);
    auto ds = dm.squared();
    ASSERT_TRUE(ds.has_value()) << label << " h=" << h
        << ": equilateral-derived length fails integrality";
    agree(*ds, dm, h, "idt");
    // The flipped diamond is a GENUINE lattice object only when the
    // original is convex (a performable flip); on non-convex originals the
    // formal flipped tuple is unrealizable and the exact regime refuses it
    // (tau <= 0) while the float regime evaluates it as Euclidean data --
    // a domain difference, not a sliver.
    if (!ds->is_convex()) continue;
    auto fsq = ds->flipped_length_sq();
    if (fsq && *fsq > 0) {
      DiamondSq fds{*fsq, ds->b, ds->d, ds->a, ds->c};
      Diamond   fdm{dm.flipped_length(), dm.b, dm.d, dm.a, dm.c};
      agree(fds, fdm, h, "flipped");
    }
  }
}

TEST(DCELExact, CorpusPredicateAgreement) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  long long n_edges = 0, n_disagree = 0;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    auto D = DelaunayTriangulation::compute(T);
    expect_predicates_agree(D, "C60 #" + std::to_string(idx), n_edges, n_disagree);
    if (::testing::Test::HasFatalFailure()) { BuckyGen::stop(Q); return; }
    idx++;
  }
  BuckyGen::stop(Q);
  EXPECT_EQ(n_disagree, 0);
  std::cout << "DCELExact corpus: " << idx << " isomers, " << n_edges
            << " diamonds (idt + flipped), exact == banded on every predicate"
            << std::endl;
}

// The constructed-tie falsifier (the recorded MUST follow-up): the crushed
// (3,1) tube -- 37 vertices, 6 cones (degrees 3,3,4,4,5,5), 31 flats,
// waist vertex 8 whose pop-time self-loop diamond is EXACTLY tied
// (integer-incircle-verified by its generator, the delaunay sub-project's
// gen_crushed_tube.py; graph pre-sorted flats-last).  On this article the
// cocircular tie-break -- and with it lattice_first_tie_side, the exact
// side order -- actually FIRES, which no fullerene of the standard corpus
// makes it do.  The firing is observed through the counting metric
// (delaunay-counting-metric.hh); the reduction must complete Ok to the
// 6-cone Delaunay complex.
//
// The tie geometry is label-independent, but WHICH diagonal the driver
// holds on the tied cocircular side-quads at pop(v) is flip-history-
// hence label-dependent: under this driver 58 of the generator's
// labeling seeds 0..128 fire (delaunay/tools/sweep_tie_exact.sh).  The
// embedded article is the smallest firing seed, regenerated by
//   python3 claude-projects/delaunay/tools/gen_crushed_tube.py <out> 2
TEST(DCELExact, ConstructedTieExercisesTieBreak) {
  static const char* tube =
      "[[13,2,1,7,31],[7,0,2,13],[0,13,1],[36,4,5],[36,6,5,3],[3,4,6,25,36],"
      "[26,24,25,5,4,36],[0,1,13,35,17,31],[15,29,23,20,22,11],"
      "[34,16,28,15,18,12],[32,19,14,33,30,21],[28,18,15,8,22,29],"
      "[27,17,34,9,18,16],[2,0,31,35,7,1],[20,32,21,33,10,19],"
      "[9,28,29,8,11,18],[34,27,12,18,28,9],[31,7,35,34,12,27],"
      "[16,12,9,15,11,28],[23,22,20,14,10,32],[8,23,32,14,19,22],"
      "[14,32,10,30,24,33],[29,11,8,20,19,23],[8,29,22,19,32,20],"
      "[33,21,30,25,6,26],[24,30,26,36,5,6],[33,24,6,36,25,30],"
      "[35,31,17,12,16,34],[9,16,18,11,29,15],[15,28,11,22,23,8],"
      "[21,10,33,26,25,24],[13,0,7,17,27,35],[20,23,19,10,21,14],"
      "[14,21,24,26,30,10],[35,27,16,9,12,17],[13,31,27,34,17,7],"
      "[26,6,4,3,5,25]]";
  FILE* f = std::tmpfile();
  ASSERT_NE(f, nullptr);
  std::fputs(tube, f);
  std::rewind(f);
  PlanarGraph g = PlanarGraph::from_ascii(f);
  std::fclose(f);
  Triangulation T(g);

  auto D = DelaunayTriangulation::from_triangulation(T);
  HostDelaunayWorkspace ws({.nv0 = D.nv, .k_max = D.nh, .nh_explicit = D.nh});
  std::vector<long long> Lsq((std::size_t)D.nh_cap, 1);
  int n_tie = 0, last_side = -1, last_float_side = -1;
  auto st = D.DelaunayView::remove_flat_vertices(
      ws, CountingExactMetric{{std::span<long long>(Lsq)}, &n_tie,
                              &last_side, &last_float_side});
  EXPECT_EQ(st, DelaunayView::Status::Ok);
  EXPECT_GE(n_tie, 1) << "the tie-break never fired: the article no longer ties";
  EXPECT_EQ(D.nv, 6);
  EXPECT_TRUE(D.is_delaunay());
  EXPECT_TRUE(D.check_consistency());
  std::cout << "DCELExact constructed tie: first_tie_side fired " << n_tie
            << " time(s); exact side " << last_side << " (float rule would say "
            << last_float_side << ", rounding noise at an exact tie); reduced to "
            << D.nv << " cones" << std::endl;

  // The debt item's byte-equality clause: on this tie-firing article the
  // exact and banded pipelines must still produce the SAME reduction
  // (integer state byte-equal, floats at rounding level) -- the tie-break
  // flip is energy-neutral and the final Lawson canonicalizes, so a
  // regime-dependent tie resolution must not leak into the output.
  long long n_isomers = 0, max_lsq = 0;
  double worst_rel = 0;
  expect_exact_matches_banded(T, "crushed (3,1) tube", n_isomers, worst_rel, max_lsq);
}

// ============================================================================
// DCELExactStatus: the REACHABLE exact-regime failure modes trip the latch
// loudly (the falsifier pattern of the DCELStatus suite).  Simple carry
// corruptions surface at the SWEEP (the tau-domain guards make the corrupt
// diamond non-Delaunay-and-unflippable), so the development's and
// tie-side's own trips sit behind it as defense in depth -- reachable only
// through corruptions that keep every diamond lattice-valid, which no
// simple injection produces.  A corrupt carry must never crash and never
// return a silently wrong reduction.
// ============================================================================

// A C60 iDT complex plus a workspace and a hand-corruptible carry.
struct ExactRig {
  DelaunayTriangulation D;
  HostDelaunayWorkspace ws;
  std::vector<long long> Lsq;
  ExactRig()
    : D(DelaunayTriangulation::compute(make_dual(60, 0))),
      ws({.nv0 = D.nv, .k_max = D.nh, .nh_explicit = D.nh}),
      Lsq((std::size_t)D.nh_cap, 0) {
    for (int h = 0; h < D.nh; h++)
      if (D.alive(h)) Lsq[h] = (long long)std::llround(D.he_length[h] * D.he_length[h]);
  }
  ExactIntegerMetric metric() { return {std::span<long long>(Lsq)}; }
};

TEST(DCELExactStatus, CorruptCarryRejectsFlipQuietly) {
  // A zero squared length was the historical SIGFPE value (division by the
  // base norm).  The tau-domain guards close that hole STRUCTURALLY: the
  // corrupted diamond is degenerate, so convexity fails and the flip is
  // rejected with the complex untouched -- no crash, no trip, no write.
  ExactRig rig;
  int h = *rig.D.edges().begin();
  rig.Lsq[h] = rig.Lsq[h ^ 1] = 0;
  EXPECT_FALSE(rig.D.DelaunayView::flip_edge(h, rig.metric()));
  EXPECT_EQ(rig.D.status, DelaunayView::Status::Ok);
}

TEST(DCELExactStatus, CorruptCarryTripsSweep) {
  // The same corruption under the SWEEP is loud: the zero-carry edge is
  // non-Delaunay under the exact predicate but unflippable, which is the
  // drain's lem:ndimpliesconvex contradiction -> InvariantViolated.
  ExactRig rig;
  int h = *rig.D.edges().begin();
  rig.Lsq[h] = rig.Lsq[h ^ 1] = 0;
  rig.ws.sweep_in_stack.clear_all();
  rig.D.DelaunayView::flip_to_delaunay(rig.ws.lawson_stack, rig.ws.sweep_in_stack,
                                       rig.metric());
  EXPECT_EQ(rig.D.status, DelaunayView::Status::InvariantViolated);
  ASSERT_NE(rig.D.status_site, nullptr);
  EXPECT_NE(std::string(rig.D.status_site).find("lawson_sweep"), std::string::npos);
}

TEST(DCELExactStatus, OverflowingCarryTripsLoudly) {
  // A whole diamond scaled past the envelope still passes convexity (it is
  // a genuine scaled lattice shape), so the flip reaches the length
  // computation -- where flipped_length_sq's INPUT envelope check refuses
  // it and the no-lattice-diagonal trip fires.  Never a silent overflow.
  ExactRig rig;
  int h = *rig.D.edges().begin();
  auto A = rig.D.diamond_arcs(h);
  for (int arc : {A.e, A.a, A.b, A.c, A.d})
    rig.Lsq[arc] = rig.Lsq[arc ^ 1] = exact_lsq_max + 1;   // scaled unit rhombus
  rig.D.DelaunayView::flip_edge(h, rig.metric());
  EXPECT_EQ(rig.D.status, DelaunayView::Status::InvariantViolated);
  ASSERT_NE(rig.D.status_site, nullptr);
  EXPECT_NE(std::string(rig.D.status_site).find("flip_edge"), std::string::npos);
}

TEST(DCELExactStatus, FlippedDiagonalPastEnvelopeTrips) {
  // The OUTPUT envelope branch: a tall in-envelope diamond whose flipped
  // diagonal exceeds the envelope.  With M = 498169 (4M - 1 = 3*815^2),
  // the diamond {1, M, M, M, M} has tau = 815 both sides, is convex
  // (P = Q = 1), all five inputs pass the input guard, and the flipped
  // diagonal is norm2((-407,815) - (408,-815)) / 1 = 3*815^2 = 1992675
  // > exact_lsq_max -- the policy's envelope trip must fire.
  ExactRig rig;
  int h = *rig.D.edges().begin();
  auto A = rig.D.diamond_arcs(h);
  rig.Lsq[A.e] = rig.Lsq[A.e ^ 1] = 1;
  for (int arc : {A.a, A.b, A.c, A.d})
    rig.Lsq[arc] = rig.Lsq[arc ^ 1] = 498169;
  DiamondSq tall{1, 498169, 498169, 498169, 498169};
  ASSERT_EQ(*tall.flipped_length_sq(), 1992675);
  rig.D.DelaunayView::flip_edge(h, rig.metric());
  EXPECT_EQ(rig.D.status, DelaunayView::Status::InvariantViolated);
  ASSERT_NE(rig.D.status_site, nullptr);
  EXPECT_NE(std::string(rig.D.status_site).find("exceeds the exact envelope"),
            std::string::npos);
}

TEST(DCELExactStatus, CorruptCarryTripsDriver) {
  // A lying carry entry (he_length says 1, the carry says sqrt(5)) makes
  // the corrupted edge non-Delaunay-and-unflippable: the removal driver
  // must trip loudly (at the initial sweep or the development), never
  // return a partial reduction.
  auto D = DelaunayTriangulation::from_triangulation(make_dual(60, 0));
  HostDelaunayWorkspace ws({.nv0 = D.nv, .k_max = D.nh, .nh_explicit = D.nh});
  std::vector<long long> Lsq((std::size_t)D.nh_cap, 1);
  Lsq[0] = Lsq[1] = 5;
  auto st = D.DelaunayView::remove_flat_vertices(ws, ExactIntegerMetric{std::span<long long>(Lsq)});
  EXPECT_EQ(st, DelaunayView::Status::InvariantViolated);
}

TEST(DCELExactStatus, ExactDriverRejectsNonExactMetric) {
  // The owner entry derives its carry and must throw, loudly, on a metric
  // that is not an exact one (the fabricated-carry hole).
  auto D = DelaunayTriangulation::from_intrinsic_metric(
      make_dual(60, 0), [](node_t, node_t) { return 1.1; });
  EXPECT_THROW(D.remove_flat_vertices_exact(), std::runtime_error);
}

TEST(DCELExactStatus, ExactDriverRejectsInsertedVertices) {
  // split_face inserts a flat vertex with a non-equilateral orig_degree
  // label: float curvature says flat, integer cone excess says cone.  The
  // exact driver must refuse the mislabeled complex rather than silently
  // leave the vertex in place.
  auto D = DelaunayTriangulation::compute(make_dual(60, 0));
  int f = -1;
  for (int i = 0; i < D.nf; i++) if (D.f_he[i] >= 0) { f = i; break; }
  ASSERT_GE(f, 0);
  D.split_face(D.f_he[f], {1.0, 1.0, 1.0});
  EXPECT_THROW(D.remove_flat_vertices_exact(), std::runtime_error);
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

// The equality relation on DCEL state: counts, capacities, status, the nine
// arrays (folded over the canonical field tuple), free-list sequences, and
// the tracker's points.  The tuple fold covers only the nine field arrays;
// everything outside the tuple is compared here BY NAME -- extend this list
// when the state grows.
static bool same_dcel_state(const DelaunayTriangulation& A,
                            const DelaunayTriangulation& B) {
  if (A.nv != B.nv || A.nh != B.nh || A.nf != B.nf) return false;
  if (A.nh_cap != B.nh_cap || A.nf_cap != B.nf_cap) return false;
  if (A.status != B.status) return false;
  bool eq = true;
  auto ta = A.to_tuple(), tb = B.to_tuple();
  auto spans_equal = [](auto sa, auto sb) {
    return sa.size() == sb.size() && std::equal(sa.begin(), sa.end(), sb.begin());
  };
  [&]<std::size_t... I>(std::index_sequence<I...>) {
    ((eq = eq && spans_equal(std::get<I>(ta), std::get<I>(tb))), ...);
  }(std::make_index_sequence<DelaunayView::n_fields>{});
  if (!eq) return false;
  if (!spans_equal(A.free_edges.live(), B.free_edges.live())) return false;
  if (!spans_equal(A.free_faces.live(), B.free_faces.live())) return false;
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
      && (int)D.f_he.size()      == D.nf
      && D.free_edges.capacity() >= D.nh_cap / 2    // DcelCapacities::free_edges_cap
      && D.free_faces.capacity() >= D.nf_cap;       // DcelCapacities::free_faces_cap
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


// ============================================================================
// DCELStatus: the view's failure channel as a falsifier group.  Each test
// names the guard it must trip; the owner's growth shadows make these
// structurally unreachable from owner paths, so a BARE view is driven.
// ============================================================================

// A bare view over an owner's arrays with the growth capacity clamped to the
// current size and an empty free list: allocation must refuse, loudly.
TEST(DCELStatus, AllocEdgeAtCapacityIsLoud) {
  auto D = DelaunayTriangulation::from_triangulation(make_dual(20, 0));
  auto D_ref = D;
  DelaunayView V = D;          // slice: same spans, same counts/caps
  ASSERT_TRUE(V.free_edges.empty());
  ASSERT_EQ(V.nh, V.nh_cap);   // fresh build is exactly full

  EXPECT_EQ(V.alloc_edge(), -1);
  EXPECT_EQ(V.status, DelaunayView::Status::CapacityExceeded);
  EXPECT_NE(V.status_site, nullptr);
  EXPECT_TRUE(same_dcel_state(D, D_ref)) << "a refused alloc must not mutate";
}

// An undersized fan workspace trips CapacityExceeded instead of overrunning.
TEST(DCELStatus, UndersizedFanWorkspaceIsLoud) {
  auto D = DelaunayTriangulation::from_triangulation(make_dual(20, 0));
  HostDelaunayWorkspace ws({.nv0 = D.nv, .k_max = 2, .nh_explicit = D.nh});
  DelaunayView& V = D;
  V.extract_fan(0, ws);        // every C20 vertex has degree 5 > 2
  EXPECT_EQ(V.status, DelaunayView::Status::CapacityExceeded);
}

// A metric violating the triangle inequality makes the sweep fail loudly --
// the throw contract downstream Armijo searches use as their reject signal.
TEST(DCELStatus, NonRealisableMetricThrowsFromSweep) {
  Triangulation T = make_dual(20, 0);
  auto bad = [](node_t u, node_t v) { return (u == 0 && v == 1) || (u == 1 && v == 0) ? 10.0 : 1.0; };
  auto D = DelaunayTriangulation::from_intrinsic_metric(T, bad);
  EXPECT_THROW(D.lawson_sweep(), std::runtime_error);
}

// After a trip the latch is sticky: every subsequent mutation is a no-op and
// the state is byte-identical (first-failure-wins is structural).
TEST(DCELStatus, StickyStatusEarlyOutsEverySubsequentMutation) {
  auto D = DelaunayTriangulation::from_triangulation(make_dual(20, 0));
  DelaunayView V = D;
  ASSERT_EQ(V.alloc_edge(), -1);                       // trip it
  ASSERT_NE(V.status, DelaunayView::Status::Ok);
  auto site = V.status_site;

  auto D_ref = D;
  EXPECT_FALSE(V.flip_edge(0));
  EXPECT_EQ(V.alloc_face(), -1);
  EXPECT_TRUE(same_dcel_state(D, D_ref));
  EXPECT_EQ(V.status_site, site) << "a later trip clobbered the first diagnostic";
}

// wire_triangle refuses a dead arc handle (the historical .at() guard).
TEST(DCELStatus, WireTriangleRejectsDeadArc) {
  auto D = DelaunayTriangulation::from_triangulation(make_dual(20, 0));
  auto D_ref = D;
  DelaunayView V = D;
  EXPECT_EQ(V.wire_triangle(-1, 0, 1), -1);
  EXPECT_EQ(V.status, DelaunayView::Status::InvariantViolated);
  EXPECT_TRUE(same_dcel_state(D, D_ref));
}

// ============================================================================
// DCELWorkspace: the layout's exactness and the wrappers' coverage.
// ============================================================================

// Every workspace span lies inside [arena, arena + bytes()) and no two spans
// overlap: bytes() and make() cannot drift (one layout list, checked).
TEST(DCELWorkspace, LayoutFitsAndSpansDisjoint) {
  DelaunayWorkspace::Layout l{.nv0 = 32, .k_max = 32, .nh_explicit = 180};
  std::vector<std::byte> arena(l.bytes());
  DelaunayWorkspace ws = l.make(std::span<std::byte>(arena));
  ASSERT_EQ(ws.k_max, 32) << "carve failed on an exactly-sized arena";

  struct Block { const void* lo; const void* hi; };
  std::vector<Block> blocks;
  auto add = [&](auto sp) {
    if (!sp.empty()) blocks.push_back({sp.data(), sp.data() + sp.size()});
  };
  add(ws.fan.nb); add(ws.fan.spoke_he); add(ws.fan.inner_rim);
  add(ws.fan.spokes); add(ws.fan.rims); add(ws.fan.cum); add(ws.fan.P);
  add(ws.poly); add(ws.rpoly);
  add(ws.tri.diagonals); add(ws.tri.triangles);
  // (The SpanStack/BitSpan members expose only their live prefixes, so the
  // fit/disjointness canary covers the eleven raw spans -- all fields come
  // from the ONE layout list, so drift in any field moves these too.)
  const void* lo = arena.data();
  const void* hi = arena.data() + l.bytes();
  for (auto& b : blocks) {
    EXPECT_GE(b.lo, lo);
    EXPECT_LE(b.hi, hi);
  }
  for (size_t i = 0; i < blocks.size(); i++)
    for (size_t j = i + 1; j < blocks.size(); j++) {
      bool disjoint = blocks[i].hi <= blocks[j].lo || blocks[j].hi <= blocks[i].lo;
      EXPECT_TRUE(disjoint) << "workspace spans " << i << " and " << j << " overlap";
    }

  // An undersized arena fails loudly: the empty workspace trips on first use.
  std::vector<std::byte> small(l.bytes() / 2);
  DelaunayWorkspace bad = l.make(std::span<std::byte>(small));
  EXPECT_EQ(bad.k_max, 0);
}

// The single-vertex removal wrapper: deg-3 merge and deg>=4 ear path, both
// leaving a consistent complex.
TEST(DCELWorkspace, RemoveFlatVertexWrapper) {
  // deg >= 4 path: remove one flat vertex of a fresh C60 dual build
  // (pre-reduction, every deg-6 vertex is flat).
  Triangulation T = make_dual(60, 0);
  auto D = DelaunayTriangulation::from_triangulation(T);
  int flat = -1;
  for (int v = 0; v < D.nv; v++)
    if (D.is_flat(v)) { flat = v; break; }
  ASSERT_GE(flat, 0);
  const int deg = D.vertex_degree(flat);
  ASSERT_GE(deg, 4);
  const int nf0 = D.nf;
  D.remove_flat_vertex(flat);
  EXPECT_LT(D.v_out[flat], 0) << "vertex not removed";
  EXPECT_TRUE(D.check_consistency());
  (void)nf0;
}
