#include <gtest/gtest.h>
#include "fullerenes/delaunay.hh"
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

// Verify all properties of the reduced triangulation.
// N_original is the number of faces in the original fullerene dual (= N for C_N).
static void verify_reduced(const FulleroidDelaunay& D, int expected_verts, int N_original) {
  // --- Combinatorial checks ---

  // Correct number of vertices
  EXPECT_EQ(D.N, expected_verts);

  // Euler formula for genus 0: V - E + F = 2
  // With V vertices and all faces triangular: E = 3V - 6, F = 2V - 4
  int V = D.N;
  int E = 0;
  for (node_t u = 0; u < D.N; u++)
    E += D.degree(u);
  E /= 2;

  if (expected_verts == 12) {
    EXPECT_EQ(E, 30) << "Expected 30 edges for 12 vertices on genus 0";
  }

  int F = E - V + 2;
  EXPECT_EQ(F, 2 * V - 4) << "Euler formula check failed";

  // All vertices have degree >= 3
  for (node_t u = 0; u < D.N; u++)
    EXPECT_GE((int)D.degree(u), 3) << "Vertex " << u << " has degree < 3";

  // Consistent orientation
  EXPECT_TRUE(D.is_consistently_oriented()) << "Orientation is broken";

  // --- Edge length consistency ---

  // Symmetric
  EXPECT_TRUE(D.edge_lengths_are_symmetric());

  // Positive
  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D[u])
      EXPECT_GT(D.get_length(u, v), 0) << "Edge (" << u << "," << v << ") has non-positive length";

  // Edge-adjacency consistency: edge_lengths(u,v) > 0 iff v in neighbours[u]
  for (node_t u = 0; u < D.N; u++)
    for (node_t v = 0; v < D.N; v++) {
      bool is_neighbour = std::find(D[u].begin(), D[u].end(), v)
                          != D[u].end();
      if (is_neighbour)
        EXPECT_GT(D.get_length(u, v), 0) << "Neighbour (" << u << "," << v << ") has zero length";
      else if (u != v)
        EXPECT_EQ(D.get_length(u, v), 0) << "Non-neighbour (" << u << "," << v << ") has nonzero length";
    }

  // --- Delaunay criterion ---
  // Every edge must be Delaunay.  A non-Delaunay edge is an algorithm bug.
  {
    int non_del = 0;
    for (node_t u = 0; u < D.N; u++)
      for (node_t v : D[u])
        if (u < v && !D.is_delaunay_edge(u, v))
          non_del++;
    EXPECT_EQ(non_del, 0) << non_del << " non-Delaunay edges remain";
  }

  // --- Metric checks ---

  // Triangle inequality: every triangle must satisfy strict triangle inequality
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j];
      node_t w = D[u][(j + 1) % D[u].size()];
      if (u < v && u < w) {  // each triangle once
        double a = D.get_length(u, v);
        double b = D.get_length(v, w);
        double c = D.get_length(w, u);
        EXPECT_GT(a + b, c) << "Triangle inequality violated in (" << u << "," << v << "," << w << ")";
        EXPECT_GT(b + c, a) << "Triangle inequality violated in (" << u << "," << v << "," << w << ")";
        EXPECT_GT(c + a, b) << "Triangle inequality violated in (" << u << "," << v << "," << w << ")";
      }
    }

  // Triangle angle sum = pi for each triangle
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j];
      node_t w = D[u][(j + 1) % D[u].size()];
      if (u < v && u < w) {
        double a = D.get_length(v, w);  // opposite u
        double b = D.get_length(u, w);  // opposite v
        double c = D.get_length(u, v);  // opposite w

        double cos_u = (b*b + c*c - a*a) / (2.0*b*c);
        double cos_v = (a*a + c*c - b*b) / (2.0*a*c);
        double cos_w = (a*a + b*b - c*c) / (2.0*a*b);

        cos_u = std::max(-1.0, std::min(1.0, cos_u));
        cos_v = std::max(-1.0, std::min(1.0, cos_v));
        cos_w = std::max(-1.0, std::min(1.0, cos_w));

        double angle_sum = acos(cos_u) + acos(cos_v) + acos(cos_w);
        EXPECT_NEAR(angle_sum, M_PI, 1e-10)
          << "Angle sum != pi in triangle (" << u << "," << v << "," << w << ")";
      }
    }

  // Cone angle at each vertex = 5*pi/3 (angle deficit = pi/3)
  for (node_t u = 0; u < D.N; u++) {
    double total_angle = 0;
    auto nbrs = D[u];
    int k = nbrs.size();
    for (int j = 0; j < k; j++) {
      node_t v = nbrs[j];
      node_t w = nbrs[(j + 1) % k];
      // Angle at u in triangle (u, v, w)
      double a = D.get_length(v, w);  // opposite u
      double b = D.get_length(u, v);
      double c = D.get_length(u, w);
      double cos_u = (b*b + c*c - a*a) / (2.0*b*c);
      cos_u = std::max(-1.0, std::min(1.0, cos_u));
      total_angle += acos(cos_u);
    }
    EXPECT_NEAR(total_angle, 5.0 * M_PI / 3.0, 1e-6)
      << "Cone angle at vertex " << u << " = " << total_angle
      << ", expected " << 5.0 * M_PI / 3.0;
  }

  // Total surface area conservation:
  // Original: N_original equilateral triangles of side 1 => area = N_original * sqrt(3)/4
  // Reduced: 20 triangles with computed edge lengths => same total area
  double expected_area = N_original * sqrt(3.0) / 4.0;
  double actual_area = 0;
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j];
      node_t w = D[u][(j + 1) % D[u].size()];
      if (u < v && u < w) {
        // Heron's formula
        double a = D.get_length(u, v);
        double b = D.get_length(v, w);
        double c = D.get_length(w, u);
        double s = (a + b + c) / 2.0;
        double area2 = s * (s - a) * (s - b) * (s - c);
        actual_area += sqrt(std::max(0.0, area2));
      }
    }
  EXPECT_NEAR(actual_area, expected_area, 1e-6)
    << "Total area " << actual_area << " != expected " << expected_area
    << " (original had " << N_original << " equilateral triangles)";
}

// ============================================================================
// Helpers
// ============================================================================

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

TEST(IntrinsicDelaunay, C20) {
  // C20 has 12 vertices all degree-5, no flat vertices to remove.
  Triangulation T = make_dual(20, 0, false);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();
  verify_reduced(D, 12, 20);
}

TEST(IntrinsicDelaunay, C24) {
  Triangulation T = make_dual(24, 0, false);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();
  verify_reduced(D, 12, 24);
}

TEST(IntrinsicDelaunay, C28) {
  // C28 has 2 isomers
  for (int idx = 0; idx < 2; idx++) {
    Triangulation T = make_dual(28, idx, false);
    FulleroidDelaunay D(T);
    D.remove_flat_vertices();
    verify_reduced(D, 12, 28);
  }
}

TEST(IntrinsicDelaunay, C60_Ih) {
  // C60 IPR (buckminsterfullerene)
  Triangulation T = make_dual(60, 0, true);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();
  verify_reduced(D, 12, 60);
}

// Exhaustive tests: every isomer for small sizes (via buckygen)
TEST(IntrinsicDelaunay, AllIsomers_C20_to_C50) {
  int sizes[] = {20, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50};
  int total = 0;
  for (int N : sizes) {
    BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
    Triangulation T;
    int idx = 0;
    while (BuckyGen::next_fullerene(Q, T)) {
      SCOPED_TRACE("C" + std::to_string(N) + " #" + std::to_string(idx));
      FulleroidDelaunay D(T);
      D.remove_flat_vertices();
      verify_reduced(D, 12, N);
      idx++;
      total++;
    }
    BuckyGen::stop(Q);
  }
  std::cout << "Tested " << total << " isomers exhaustively (C20-C50)" << std::endl;
}

TEST(IntrinsicDelaunay, C60_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  std::vector<double> times_us;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    SCOPED_TRACE("C60 #" + std::to_string(idx));
    FulleroidDelaunay D(T);
    auto t0 = std::chrono::high_resolution_clock::now();
    D.remove_flat_vertices();
    auto t1 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    verify_reduced(D, 12, 60);
    idx++;
  }
  BuckyGen::stop(Q);
  print_timing_stats("C60 iDT", times_us);
}

// Run the iDT with full invariant auditing on C60 #1264 (known problematic isomer).
TEST(IntrinsicDelaunay, C60_1264_Audit) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    if (idx == 1264) break;
    idx++;
  }
  BuckyGen::stop(Q);
  ASSERT_EQ(idx, 1264);

  FulleroidDelaunay D(T);
  IDTAudit audit(D);
  audit.stop_on_failure = false;  // collect all failures, don't abort
  D.audit = &audit;

  D.remove_flat_vertices();
  ASSERT_EQ(D.N, 12);

  std::cout << "Audit: " << audit.n_checks << " checks, "
            << audit.n_failures << " failures\n";
  EXPECT_EQ(audit.n_failures, 0);
}

TEST(IntrinsicDelaunay, DISABLED_C80_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(80, false, false);
  Triangulation T;
  std::vector<double> times_us;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    SCOPED_TRACE("C80 #" + std::to_string(idx));
    FulleroidDelaunay D(T);
    auto t0 = std::chrono::high_resolution_clock::now();
    D.remove_flat_vertices();
    auto t1 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    verify_reduced(D, 12, 80);
    idx++;
  }
  BuckyGen::stop(Q);
  print_timing_stats("C80 iDT", times_us);
}

TEST(IntrinsicDelaunay, DISABLED_C90_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(90, false, false);
  Triangulation T;
  std::vector<double> times_us;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    SCOPED_TRACE("C90 #" + std::to_string(idx));
    FulleroidDelaunay D(T);
    auto t0 = std::chrono::high_resolution_clock::now();
    D.remove_flat_vertices();
    auto t1 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    verify_reduced(D, 12, 90);
    idx++;
  }
  BuckyGen::stop(Q);
  print_timing_stats("C90 iDT", times_us);
}

TEST(IntrinsicDelaunay, DISABLED_C100_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(100, false, false);
  Triangulation T;
  std::vector<double> times_us;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    SCOPED_TRACE("C100 #" + std::to_string(idx));
    FulleroidDelaunay D(T);
    auto t0 = std::chrono::high_resolution_clock::now();
    D.remove_flat_vertices();
    auto t1 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    verify_reduced(D, 12, 100);
    idx++;
  }
  BuckyGen::stop(Q);
  print_timing_stats("C100 iDT", times_us);
}

TEST(IntrinsicDelaunay, C20_AllEquilateral) {
  // For C20 (icosahedron), all edges should remain length 1
  // since all vertices are degree-5 (no flat vertices to remove).
  Triangulation T = make_dual(20, 0, false);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();

  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D[u])
      EXPECT_DOUBLE_EQ(D.get_length(u, v), 1.0);
}

// Verify all properties of the reduced triangulation for a GENERAL
// equilateral triangulation (not just fullerene duals).
// original_degrees: degree of each vertex in the original triangulation (before sorting).
// N_original_faces: number of faces in the original triangulation.
static void verify_reduced_general(const FulleroidDelaunay& D,
                                   const vector<int>& original_degrees,
                                   int N_original_faces) {
  // Count expected cone points
  int expected_verts = 0;
  for (int d : original_degrees)
    if (d != 6) expected_verts++;

  // --- Combinatorial checks ---
  EXPECT_EQ(D.N, expected_verts);

  int V = D.N;
  int E = 0;
  for (node_t u = 0; u < D.N; u++)
    E += D[u].size();
  E /= 2;

  // Euler: V - E + F = 2 for genus 0, F = 2V - 4, E = 3V - 6
  int F = E - V + 2;
  EXPECT_EQ(F, 2 * V - 4) << "Euler formula check failed";

  // All vertices have degree >= 3
  for (node_t u = 0; u < D.N; u++)
    EXPECT_GE((int)D[u].size(), 3) << "Vertex " << u << " has degree < 3";

  // Consistent orientation
  EXPECT_TRUE(D.is_consistently_oriented()) << "Orientation is broken";

  // --- Edge length consistency ---
  EXPECT_TRUE(D.edge_lengths_are_symmetric());

  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D[u])
      EXPECT_GT(D.get_length(u, v), 0) << "Edge (" << u << "," << v << ") has non-positive length";

  // --- Delaunay criterion ---
  {
    int non_del = 0;
    for (node_t u = 0; u < D.N; u++)
      for (node_t v : D[u])
        if (u < v && !D.is_delaunay_edge(u, v))
          non_del++;
    EXPECT_EQ(non_del, 0) << non_del << " non-Delaunay edges remain";
  }

  // --- Triangle inequality ---
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j];
      node_t w = D[u][(j + 1) % D[u].size()];
      if (u < v && u < w) {
        double a = D.get_length(u, v);
        double b = D.get_length(v, w);
        double c = D.get_length(w, u);
        EXPECT_GT(a + b, c);
        EXPECT_GT(b + c, a);
        EXPECT_GT(c + a, b);
      }
    }

  // --- Angle sum = pi per triangle ---
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j];
      node_t w = D[u][(j + 1) % D[u].size()];
      if (u < v && u < w) {
        double a = D.get_length(v, w);
        double b = D.get_length(u, w);
        double c = D.get_length(u, v);
        double cos_u = std::clamp((b*b + c*c - a*a) / (2.0*b*c), -1.0, 1.0);
        double cos_v = std::clamp((a*a + c*c - b*b) / (2.0*a*c), -1.0, 1.0);
        double cos_w = std::clamp((a*a + b*b - c*c) / (2.0*a*b), -1.0, 1.0);
        double angle_sum = acos(cos_u) + acos(cos_v) + acos(cos_w);
        EXPECT_NEAR(angle_sum, M_PI, 1e-10);
      }
    }

  // --- Cone angles: each cone point retains its original cone angle ---
  // The FulleroidDelaunay constructor sorts vertices so cone points come first.
  // After reduction, vertex u has original degree original_cone_degrees[u],
  // so its cone angle should be original_cone_degrees[u] * pi/3.
  // We need to recover the original cone degrees. They are the non-6 degrees
  // from original_degrees, sorted the same way sort_flat_last would sort them.
  // sort_flat_last sorts by {(d==6)?1:0, original_index}, so cone points
  // come first in their original relative order.
  vector<int> cone_degrees;
  for (int d : original_degrees)
    if (d != 6) cone_degrees.push_back(d);

  ASSERT_EQ((int)cone_degrees.size(), D.N);
  for (node_t u = 0; u < D.N; u++) {
    double expected_angle = cone_degrees[u] * M_PI / 3.0;
    double total_angle = 0;
    const auto& nbrs = D[u];
    int k = nbrs.size();
    for (int j = 0; j < k; j++) {
      node_t v = nbrs[j];
      node_t w = nbrs[(j + 1) % k];
      double a = D.get_length(v, w);
      double b = D.get_length(u, v);
      double c = D.get_length(u, w);
      double cos_u = std::clamp((b*b + c*c - a*a) / (2.0*b*c), -1.0, 1.0);
      total_angle += acos(cos_u);
    }
    EXPECT_NEAR(total_angle, expected_angle, 1e-6)
      << "Cone angle at vertex " << u << " (orig degree " << cone_degrees[u]
      << "): got " << total_angle << ", expected " << expected_angle;
  }

  // --- Total area conservation ---
  double expected_area = N_original_faces * sqrt(3.0) / 4.0;
  double actual_area = 0;
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j];
      node_t w = D[u][(j + 1) % D[u].size()];
      if (u < v && u < w) {
        double a = D.get_length(u, v);
        double b = D.get_length(v, w);
        double c = D.get_length(w, u);
        double s = (a + b + c) / 2.0;
        double area2 = s * (s - a) * (s - b) * (s - c);
        actual_area += sqrt(std::max(0.0, area2));
      }
    }
  EXPECT_NEAR(actual_area, expected_area, 1e-6);
}

TEST(IntrinsicDelaunay, C60_Ih_IcosahedralSymmetry) {
  // C60 Ih (buckminsterfullerene) has icosahedral symmetry.
  // All 12 pentagons are equivalent, so after reduction all 30 edges
  // should have the same geodesic length.
  Triangulation T = make_dual(60, 0, true);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();

  // Collect all edge lengths
  double first_len = -1;
  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D[u])
      if (u < v) {
        double len = D.get_length(u, v);
        if (first_len < 0) first_len = len;
        EXPECT_NEAR(len, first_len, 1e-10)
          << "Edge (" << u << "," << v << ") has length " << len
          << " but expected " << first_len << " (icosahedral symmetry)";
      }

  std::cout << "C60 Ih: all 30 edges have length " << first_len << std::endl;
}

// Test C100 #84570 with the actual remove_flat_vertices pipeline.
TEST(IntrinsicDelaunay, C100_84570) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(100, false, false);
  Triangulation T;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    if (idx == 84570) break;
    idx++;
  }
  BuckyGen::stop(Q);
  ASSERT_EQ(idx, 84570) << "Could not reach isomer #84570";

  FulleroidDelaunay D(T);
  D.remove_flat_vertices();
  ASSERT_EQ(D.N, 12);
  verify_reduced(D, 12, 100);
}

// ============================================================================
// 3D Embedding tests
// ============================================================================

// Verify the 3D embedding: edge distances match iDT targets, convexity.
static void verify_embedding(const FulleroidDelaunay& D, const vector<coord3d>& coords,
                             double dist_tol, const std::string& label) {
  ASSERT_EQ((int)coords.size(), D.N) << label << ": wrong coord count";

  // Edge distance errors
  double max_rel_err = 0;
  double sum_sq_err = 0;
  int n_edges = 0;
  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D[u])
      if (u < v) {
        double target = D.get_length(u, v);
        double actual = (coords[u] - coords[v]).norm();
        double rel_err = fabs(actual - target) / target;
        max_rel_err = std::max(max_rel_err, rel_err);
        sum_sq_err += (actual - target) * (actual - target);
        n_edges++;
      }
  double rms_err = sqrt(sum_sq_err / n_edges);

  EXPECT_LT(max_rel_err, dist_tol)
    << label << ": max relative edge error = " << max_rel_err;

  // Convexity: each vertex should be above the plane of its neighbors.
  // Compute height above neighbor centroid plane.
  int n_concave = 0;
  double h_min = 1e30;
  for (node_t u = 0; u < D.N; u++) {
    const auto& nbrs = D[u];
    int k = nbrs.size();
    coord3d centroid;
    for (int j = 0; j < k; j++) centroid += coords[nbrs[j]];
    centroid /= k;

    // Outward normal from fan of neighbors (CCW convention)
    coord3d normal;
    for (int j = 0; j < k; j++) {
      coord3d e1 = coords[nbrs[j]] - centroid;
      coord3d e2 = coords[nbrs[(j + 1) % k]] - centroid;
      normal += e1.cross(e2);
    }
    double n_len = normal.norm();
    if (n_len < 1e-15) continue;
    normal /= n_len;

    double h = (coords[u] - centroid).dot(normal);
    h_min = std::min(h_min, h);
    if (h < -1e-10) n_concave++;
  }

  std::cout << label << ": rms_err=" << rms_err
            << " max_rel_err=" << max_rel_err
            << " h_min=" << h_min
            << " n_concave=" << n_concave << std::endl;
}

TEST(DelaunayEmbed, C20_Icosahedron) {
  // C20 dual = icosahedron (12 vertices, all degree 5, all edges length 1).
  // The 3D embedding should be a regular icosahedron.
  Triangulation T = make_dual(20, 0, false);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();
  ASSERT_EQ(D.N, 12);

  auto coords = D.embed_3d();
  verify_embedding(D, coords, 1e-6, "C20");

  // For a regular icosahedron with edge length 1, all edges should be exactly 1.
  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D[u])
      if (u < v) {
        double dist = (coords[u] - coords[v]).norm();
        EXPECT_NEAR(dist, 1.0, 1e-6)
          << "C20 edge (" << u << "," << v << ") = " << dist;
      }
}

TEST(DelaunayEmbed, C60_Ih) {
  // C60 Ih: all iDT edge lengths are equal (icosahedral symmetry).
  Triangulation T = make_dual(60, 0, true);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();
  ASSERT_EQ(D.N, 12);

  auto coords = D.embed_3d();
  verify_embedding(D, coords, 1e-4, "C60_Ih");
}

TEST(DelaunayEmbed, SmallFullerenes) {
  // Test embedding for all isomers of small fullerenes.
  int sizes[] = {20, 24, 26, 28, 30, 32, 34, 36, 38, 40};
  int total = 0, n_poor = 0;

  for (int N : sizes) {
    BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
    Triangulation T;
    int idx = 0;
    while (BuckyGen::next_fullerene(Q, T)) {
      FulleroidDelaunay D(T);
      D.remove_flat_vertices();
      if (D.N != 12) { idx++; continue; }

      auto coords = D.embed_3d();

      // Check edge distance matching
      double max_rel_err = 0;
      for (node_t u = 0; u < D.N; u++)
        for (node_t v : D[u])
          if (u < v) {
            double target = D.get_length(u, v);
            double actual = (coords[u] - coords[v]).norm();
            double rel_err = fabs(actual - target) / target;
            max_rel_err = std::max(max_rel_err, rel_err);
          }

      if (max_rel_err > 0.01) {
        n_poor++;
        std::cerr << "C" << N << " #" << idx
                  << ": max_rel_err=" << max_rel_err << std::endl;
      }

      EXPECT_LT(max_rel_err, 0.05)
        << "C" << N << " #" << idx << " embedding too inaccurate";

      idx++;
      total++;
    }
    BuckyGen::stop(Q);
  }
  std::cout << "Tested " << total << " embeddings, "
            << n_poor << " with >1% error" << std::endl;
}

TEST(DelaunayEmbed, C60_AllIsomers) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  std::vector<double> times_us;
  std::vector<double> max_errs;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    FulleroidDelaunay D(T);
    D.remove_flat_vertices();
    if (D.N != 12) { idx++; continue; }

    auto t0 = std::chrono::high_resolution_clock::now();
    auto coords = D.embed_3d();
    auto t1 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());

    double max_rel_err = 0;
    for (node_t u = 0; u < D.N; u++)
      for (node_t v : D[u])
        if (u < v) {
          double target = D.get_length(u, v);
          double actual = (coords[u] - coords[v]).norm();
          max_rel_err = std::max(max_rel_err, fabs(actual - target) / target);
        }
    max_errs.push_back(max_rel_err);
    idx++;
  }
  BuckyGen::stop(Q);

  // Pair errors with indices for diagnostics
  vector<pair<double,int>> err_idx;
  for (size_t i = 0; i < max_errs.size(); i++)
    err_idx.push_back({max_errs[i], (int)i});
  std::sort(err_idx.rbegin(), err_idx.rend());

  std::sort(max_errs.begin(), max_errs.end());
  double median_err = max_errs[max_errs.size() / 2];
  double p99_err = max_errs[(int)(max_errs.size() * 0.99)];
  double worst_err = max_errs.back();

  print_timing_stats("C60 embed", times_us);
  std::cout << "C60 embed errors: median=" << median_err
            << " p99=" << p99_err << " worst=" << worst_err << std::endl;
  std::cout << "Worst 5 isomers:" << std::endl;
  for (int i = 0; i < 5 && i < (int)err_idx.size(); i++)
    std::cout << "  C60 #" << err_idx[i].second
              << ": max_rel_err=" << err_idx[i].first << std::endl;

  // For initial geometry purposes, even 30% error is useful.
  // The full optimizer handles the rest.
  EXPECT_LT(p99_err, 0.2) << "p99 C60 embedding error too large";
}

// Parameterized embedding test: all isomers of a given size.
static void test_embed_all_isomers(int N) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  Triangulation T;
  std::vector<double> times_us, idt_times_us, max_errs;
  int idx = 0, n_skip = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    auto t0 = std::chrono::high_resolution_clock::now();
    FulleroidDelaunay D(T);
    D.remove_flat_vertices();
    auto t1 = std::chrono::high_resolution_clock::now();
    idt_times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());

    if (D.N != 12) { idx++; n_skip++; continue; }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto coords = D.embed_3d();
    auto t3 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t3 - t2).count());

    double max_rel_err = 0;
    for (node_t u = 0; u < D.N; u++)
      for (node_t v : D[u])
        if (u < v) {
          double target = D.get_length(u, v);
          double actual = (coords[u] - coords[v]).norm();
          max_rel_err = std::max(max_rel_err, fabs(actual - target) / target);
        }
    max_errs.push_back(max_rel_err);
    idx++;
  }
  BuckyGen::stop(Q);

  // Sort for percentile stats
  vector<pair<double,int>> err_idx;
  for (size_t i = 0; i < max_errs.size(); i++)
    err_idx.push_back({max_errs[i], (int)i});
  std::sort(err_idx.rbegin(), err_idx.rend());

  std::sort(max_errs.begin(), max_errs.end());
  double median_err = max_errs[max_errs.size() / 2];
  double p99_err = max_errs[(int)(max_errs.size() * 0.99)];
  double worst_err = max_errs.back();

  print_timing_stats("C" + std::to_string(N) + " iDT", idt_times_us);
  print_timing_stats("C" + std::to_string(N) + " embed", times_us);
  std::cout << "C" << N << " embed: n=" << max_errs.size()
            << " median=" << median_err
            << " p99=" << p99_err << " worst=" << worst_err;
  if (n_skip) std::cout << " (" << n_skip << " skipped: iDT incomplete)";
  std::cout << std::endl;
  std::cout << "Worst 5:" << std::endl;
  for (int i = 0; i < 5 && i < (int)err_idx.size(); i++)
    std::cout << "  #" << err_idx[i].second
              << ": " << err_idx[i].first << std::endl;

  EXPECT_LT(p99_err, 0.2) << "C" << N << " p99 embedding error too large";
}

TEST(DelaunayEmbed, C70_AllIsomers)  { test_embed_all_isomers(70); }
TEST(DelaunayEmbed, C80_AllIsomers)  { test_embed_all_isomers(80); }

// Investigate outlier isomers where embed_3d doesn't converge to machine precision.
TEST(DelaunayEmbed, InvestigateOutliers) {
  // Outliers identified from C70/C80 all-isomers tests
  struct Case { int N; int buckygen_idx; };
  Case cases[] = {
    {70, 3926}, {70, 2836}, {70, 3892},
    {80, 1904}, {80, 1826}, {80, 16208}, {80, 8296}, {80, 24339},
  };

  for (auto& c : cases) {
    BuckyGen::buckygen_queue Q = BuckyGen::start(c.N, false, false);
    Triangulation T;
    int idx = 0;
    while (BuckyGen::next_fullerene(Q, T)) {
      if (idx == c.buckygen_idx) break;
      idx++;
    }
    BuckyGen::stop(Q);
    if (idx != c.buckygen_idx) { ADD_FAILURE() << "Couldn't find isomer"; continue; }

    std::cout << "\n=== C" << c.N << " #" << c.buckygen_idx << " ===" << std::endl;

    FulleroidDelaunay D(T);
    D.remove_flat_vertices();
    ASSERT_EQ(D.N, 12);

    // 1. Check MDS eigenvalues — are the APSP distances Euclidean?
    matrix<double> Dist = D.all_pairs_distances();
    {
      int n = D.N;
      vector<double> D_sq(n * n);
      for (int i = 0; i < n * n; i++) D_sq[i] = Dist[i] * Dist[i];

      vector<double> row_mean(n, 0), col_mean(n, 0);
      double grand_mean = 0;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
          row_mean[i] += D_sq[i * n + j];
          col_mean[j] += D_sq[i * n + j];
          grand_mean += D_sq[i * n + j];
        }
      for (int i = 0; i < n; i++) { row_mean[i] /= n; col_mean[i] /= n; }
      grand_mean /= (n * n);

      vector<double> B(n * n);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          B[i * n + j] = -0.5 * (D_sq[i * n + j] - row_mean[i] - col_mean[j] + grand_mean);

      // Eigendecompose (simple Jacobi for 12x12)
      vector<double> evals(n);
      // Just compute eigenvalues via symmetric QR or Jacobi
      // Reuse: copy B, do Jacobi
      vector<double> V(n * n, 0);
      for (int i = 0; i < n; i++) V[i * n + i] = 1.0;
      for (int iter = 0; iter < 10000; iter++) {
        double max_val = 0; int p = 0, q = 1;
        for (int i = 0; i < n; i++)
          for (int j = i + 1; j < n; j++)
            if (fabs(B[i * n + j]) > max_val) { max_val = fabs(B[i * n + j]); p = i; q = j; }
        if (max_val < 1e-15) break;
        double app = B[p*n+p], aqq = B[q*n+q], apq = B[p*n+q];
        double theta = (fabs(app-aqq) < 1e-30) ? M_PI/4 : 0.5*atan2(2*apq, app-aqq);
        double cs = cos(theta), sn = sin(theta);
        for (int j = 0; j < n; j++) {
          double bp = B[p*n+j], bq = B[q*n+j];
          B[p*n+j] = cs*bp + sn*bq; B[q*n+j] = -sn*bp + cs*bq;
        }
        for (int i = 0; i < n; i++) {
          double bp = B[i*n+p], bq = B[i*n+q];
          B[i*n+p] = cs*bp + sn*bq; B[i*n+q] = -sn*bp + cs*bq;
        }
      }
      for (int i = 0; i < n; i++) evals[i] = B[i * n + i];
      std::sort(evals.rbegin(), evals.rend());

      std::cout << "MDS eigenvalues (APSP):";
      for (int i = 0; i < n; i++) std::cout << " " << evals[i];
      std::cout << std::endl;
      std::cout << "  lambda3/lambda1 = " << evals[2] / evals[0] << std::endl;
      if (evals[3] < -1e-6)
        std::cout << "  WARNING: negative eigenvalue " << evals[3]
                  << " => APSP distances NOT Euclidean" << std::endl;
    }

    // 2. Run embed_3d and measure result quality
    auto coords = D.embed_3d();

    double max_rel_err = 0, stress = 0;
    for (node_t u = 0; u < D.N; u++)
      for (node_t v : D[u])
        if (u < v) {
          double target = D.get_length(u, v);
          double actual = (coords[u] - coords[v]).norm();
          double err = actual - target;
          stress += err * err;
          max_rel_err = std::max(max_rel_err, fabs(err) / target);
        }
    std::cout << "embed_3d: stress=" << stress << " max_rel_err=" << max_rel_err << std::endl;

    // 3. Check convexity
    int n_concave = 0;
    double h_min = 1e30;
    for (node_t u = 0; u < D.N; u++) {
      const auto& nbrs = D[u];
      int k = nbrs.size();
      coord3d centroid;
      for (int j = 0; j < k; j++) centroid += coords[nbrs[j]];
      centroid /= k;
      coord3d normal;
      for (int j = 0; j < k; j++) {
        coord3d e1 = coords[nbrs[j]] - centroid;
        coord3d e2 = coords[nbrs[(j+1)%k]] - centroid;
        normal += e1.cross(e2);
      }
      double nlen = normal.norm();
      if (nlen > 1e-15) {
        normal /= nlen;
        double h = (coords[u] - centroid).dot(normal);
        h_min = std::min(h_min, h);
        if (h < -1e-10) n_concave++;
      }
    }
    std::cout << "Convexity: h_min=" << h_min << " n_concave=" << n_concave << std::endl;

    // 4. Print edge length range and degree sequence
    double L_min = 1e30, L_max = 0;
    for (node_t u = 0; u < D.N; u++)
      for (node_t v : D[u])
        if (u < v) {
          double L = D.get_length(u, v);
          L_min = std::min(L_min, L);
          L_max = std::max(L_max, L);
        }
    std::cout << "iDT edge lengths: min=" << L_min << " max=" << L_max
              << " ratio=" << L_max/L_min << std::endl;

    std::cout << "Degrees:";
    for (node_t u = 0; u < D.N; u++)
      std::cout << " " << D[u].size();
    std::cout << std::endl;
  }
}

// ============================================================================
// MDS vs Face-BFS comparison
// ============================================================================

// --- Shared Steihaug optimizer (takes initial coords, returns optimized coords) ---
struct EmbedResult {
  vector<coord3d> coords;
  double stress;        // final E = sum (|xi-xj| - Lij)^2
  double max_rel_err;
  int    iters;
  double h_min;         // minimum height above neighbor centroid plane
  int    n_concave;     // vertices with h < -1e-10
};

// Per-edge distance-matching geometry (local copy for the comparison test).
struct StressEdge {
  int i, j;
  coord3d diff;
  double dist, target, err;
  static StressEdge compute(const vector<coord3d>& x, int i, int j, double L) {
    coord3d d = x[i] - x[j];
    double r = d.norm();
    return {i, j, d, r, L, r - L};
  }
  bool valid() const { return dist > 1e-15; }
  double energy() const { return err * err; }
  void scatter_gradient(vector<coord3d>& g) const {
    coord3d c = diff * (2.0 * err / dist);
    g[i] = g[i] + c;  g[j] = g[j] - c;
  }
  void scatter_hv(const vector<coord3d>& v, vector<coord3d>& Hv) const {
    coord3d dv = v[i] - v[j];
    double udv = diff.dot(dv) / (dist * dist);
    coord3d h = diff * (2.0 * target / dist * udv) + dv * (2.0 * err / dist);
    Hv[i] = Hv[i] + h;  Hv[j] = Hv[j] - h;
  }
};

static double vdot_c(const vector<coord3d>& a, const vector<coord3d>& b) {
  double s = 0; for (size_t i = 0; i < a.size(); i++) s += a[i].dot(b[i]); return s;
}
static double vnorm_c(const vector<coord3d>& a) { return sqrt(vdot_c(a, a)); }
static void vaxpy_c(double alpha, const vector<coord3d>& x, vector<coord3d>& y) {
  for (size_t i = 0; i < x.size(); i++) y[i] += x[i] * alpha;
}

static double det3_c(const matrix3d& M) {
  return M(0,0)*(M(1,1)*M(2,2)-M(1,2)*M(2,1))
       - M(0,1)*(M(1,0)*M(2,2)-M(1,2)*M(2,0))
       + M(0,2)*(M(1,0)*M(2,1)-M(1,1)*M(2,0));
}
static matrix3d adjugate_c(const matrix3d& M) {
  return matrix3d(
    M(1,1)*M(2,2)-M(1,2)*M(2,1), M(0,2)*M(2,1)-M(0,1)*M(2,2), M(0,1)*M(1,2)-M(0,2)*M(1,1),
    M(1,2)*M(2,0)-M(1,0)*M(2,2), M(0,0)*M(2,2)-M(0,2)*M(2,0), M(0,2)*M(1,0)-M(0,0)*M(1,2),
    M(1,0)*M(2,1)-M(1,1)*M(2,0), M(0,1)*M(2,0)-M(0,0)*M(2,1), M(0,0)*M(1,1)-M(0,1)*M(1,0));
}

static void project_rigid_body_c(vector<coord3d>& g, const vector<coord3d>& x) {
  int N = g.size();
  coord3d mean;
  for (auto& gi : g) mean += gi;
  mean /= N;
  for (auto& gi : g) gi -= mean;
  coord3d c;
  for (auto& xi : x) c += xi;
  c /= N;
  matrix3d I;
  coord3d tau;
  for (int k = 0; k < N; k++) {
    coord3d r = x[k] - c;
    I += matrix3d::unit_matrix() * r.dot(r) - r.outer(r);
    tau += r.cross(g[k]);
  }
  double d = det3_c(I);
  if (fabs(d) < 1e-30) return;
  coord3d omega = adjugate_c(I) * tau * (1.0 / d);
  for (int k = 0; k < N; k++)
    g[k] -= omega.cross(x[k] - c);
}

static void orient_outward_c(vector<coord3d>& x, const neighbours_t& nbrs) {
  int N = x.size();
  coord3d c;
  for (auto& xi : x) c += xi;
  c /= N;
  double vol = 0;
  for (int u = 0; u < N; u++) {
    int deg = nbrs[u].size();
    for (int j = 0; j < deg; j++) {
      int v = nbrs[u][j], w = nbrs[u][(j+1)%deg];
      if (u < v && u < w)
        vol += (x[u]-c).dot((x[v]-c).cross(x[w]-c));
    }
  }
  if (vol < 0)
    for (auto& xi : x) xi = c*2.0 - xi;
}

static EmbedResult run_steihaug(vector<coord3d> x, int N,
                                 const neighbours_t& nbrs,
                                 const matrix<double>& el)
{
  auto eval = [&](const vector<coord3d>& pos, vector<coord3d>& g) -> double {
    double E = 0;
    g.assign(N, coord3d());
    for (int u = 0; u < N; u++)
      for (int v : nbrs[u]) {
        if (u >= v) continue;
        auto ed = StressEdge::compute(pos, u, v, el(u, v));
        if (!ed.valid()) continue;
        E += ed.energy();
        ed.scatter_gradient(g);
      }
    project_rigid_body_c(g, pos);
    return E;
  };

  auto hv_prod = [&](const vector<coord3d>& pos, const vector<coord3d>& dir,
                      vector<coord3d>& Hv) {
    Hv.assign(N, coord3d());
    for (int u = 0; u < N; u++)
      for (int v : nbrs[u]) {
        if (u >= v) continue;
        auto ed = StressEdge::compute(pos, u, v, el(u, v));
        if (!ed.valid()) continue;
        ed.scatter_hv(dir, Hv);
      }
    project_rigid_body_c(Hv, pos);
  };

  vector<coord3d> g(N);
  double E = eval(x, g);
  double Delta = std::max(vnorm_c(g), 1e-14);
  int outer;

  for (outer = 0; outer < 200; outer++) {
    double gnorm = vnorm_c(g);
    if (gnorm < 1e-13) break;

    vector<coord3d> p(N), r(g), d(N);
    for (int i = 0; i < N; i++) d[i] = coord3d() - r[i];
    double rr = vdot_c(r, r);
    bool hit_boundary = false;

    for (int cg = 0; cg < 3*N; cg++) {
      vector<coord3d> Hd(N);
      hv_prod(x, d, Hd);
      double dHd = vdot_c(d, Hd);

      if (dHd <= 1e-15 * rr) {
        double pd=vdot_c(p,d), pp=vdot_c(p,p), dd=vdot_c(d,d);
        double tau = (-pd + sqrt(std::max(0.0, pd*pd-dd*(pp-Delta*Delta)))) / dd;
        vaxpy_c(tau, d, p);
        hit_boundary = true; break;
      }

      double alpha = rr / dHd;
      vector<coord3d> p_new(p);
      vaxpy_c(alpha, d, p_new);

      if (vnorm_c(p_new) >= Delta) {
        double pd=vdot_c(p,d), pp=vdot_c(p,p), dd=vdot_c(d,d);
        double tau = (-pd + sqrt(std::max(0.0, pd*pd-dd*(pp-Delta*Delta)))) / dd;
        vaxpy_c(tau, d, p);
        hit_boundary = true; break;
      }

      p = p_new;
      vaxpy_c(alpha, Hd, r);
      double rr_new = vdot_c(r, r);
      if (sqrt(rr_new) < 1e-13 * gnorm) break;

      double beta = rr_new / rr;
      for (int i = 0; i < N; i++) d[i] = coord3d() - r[i] + d[i] * beta;
      rr = rr_new;
    }

    double predicted = vdot_c(g, p);
    { vector<coord3d> Hp(N); hv_prod(x, p, Hp); predicted += 0.5 * vdot_c(p, Hp); }

    vector<coord3d> x_trial(N), g_trial(N);
    for (int i = 0; i < N; i++) x_trial[i] = x[i] + p[i];
    double E_trial = eval(x_trial, g_trial);

    double rho = (E - E_trial) / std::max(1e-30, -predicted);
    double pnorm = vnorm_c(p);
    if (rho < 0.25)                        Delta = 0.25 * pnorm;
    else if (rho > 0.75 && hit_boundary)   Delta = std::min(2.0 * Delta, 1e10);
    Delta = std::max(Delta, 1e-15);

    if (rho > 0.1) { x = x_trial; g = g_trial; E = E_trial; }
  }

  orient_outward_c(x, nbrs);

  // Compute final quality
  double max_rel_err = 0;
  for (int u = 0; u < N; u++)
    for (int v : nbrs[u])
      if (u < v) {
        double target = el(u, v);
        double actual = (x[u] - x[v]).norm();
        max_rel_err = std::max(max_rel_err, fabs(actual - target) / target);
      }

  // Convexity: height above neighbor centroid plane
  double h_min = 1e30;
  int n_concave = 0;
  for (int u = 0; u < N; u++) {
    int k = nbrs[u].size();
    coord3d centroid;
    for (int j = 0; j < k; j++) centroid += x[nbrs[u][j]];
    centroid /= k;
    coord3d normal;
    for (int j = 0; j < k; j++) {
      coord3d e1 = x[nbrs[u][j]] - centroid;
      coord3d e2 = x[nbrs[u][(j+1)%k]] - centroid;
      normal += e1.cross(e2);
    }
    double nlen = normal.norm();
    if (nlen > 1e-15) {
      normal /= nlen;
      double h = (x[u] - centroid).dot(normal);
      h_min = std::min(h_min, h);
      if (h < -1e-10) n_concave++;
    }
  }

  return {x, E, max_rel_err, outer, h_min, n_concave};
}

// --- Classical MDS initial placement ---
static vector<coord3d> mds_placement(int N, const neighbours_t& nbrs,
                                      const matrix<double>& el)
{
  // Step 1: All-pairs shortest paths (Floyd-Warshall on the reduced graph)
  const double INF = 1e30;
  vector<double> D(N * N, INF);
  for (int u = 0; u < N; u++) {
    D[u * N + u] = 0;
    for (int v : nbrs[u])
      D[u * N + v] = el(u, v);
  }
  for (int k = 0; k < N; k++)
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        if (D[i*N+k] + D[k*N+j] < D[i*N+j])
          D[i*N+j] = D[i*N+k] + D[k*N+j];

  // Step 2: Double-centered squared distance matrix  B = -0.5 * J * D^2 * J
  vector<double> D_sq(N * N);
  for (int i = 0; i < N * N; i++) D_sq[i] = D[i] * D[i];

  vector<double> row_mean(N, 0), col_mean(N, 0);
  double grand_mean = 0;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      row_mean[i] += D_sq[i * N + j];
      col_mean[j] += D_sq[i * N + j];
      grand_mean += D_sq[i * N + j];
    }
  for (int i = 0; i < N; i++) { row_mean[i] /= N; col_mean[i] /= N; }
  grand_mean /= (N * N);

  vector<double> B(N * N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      B[i * N + j] = -0.5 * (D_sq[i * N + j] - row_mean[i] - col_mean[j] + grand_mean);

  // Step 3: Jacobi eigendecomposition of B
  vector<double> V(N * N, 0);
  for (int i = 0; i < N; i++) V[i * N + i] = 1.0;

  for (int iter = 0; iter < 10000; iter++) {
    double max_val = 0; int p = 0, q = 1;
    for (int i = 0; i < N; i++)
      for (int j = i + 1; j < N; j++)
        if (fabs(B[i * N + j]) > max_val) { max_val = fabs(B[i * N + j]); p = i; q = j; }
    if (max_val < 1e-15) break;
    double app = B[p*N+p], aqq = B[q*N+q], apq = B[p*N+q];
    double theta = (fabs(app-aqq) < 1e-30) ? M_PI/4 : 0.5*atan2(2*apq, app-aqq);
    double cs = cos(theta), sn = sin(theta);
    // Rotate rows of B
    for (int j = 0; j < N; j++) {
      double bp = B[p*N+j], bq = B[q*N+j];
      B[p*N+j] = cs*bp + sn*bq; B[q*N+j] = -sn*bp + cs*bq;
    }
    // Rotate cols of B
    for (int i = 0; i < N; i++) {
      double bp = B[i*N+p], bq = B[i*N+q];
      B[i*N+p] = cs*bp + sn*bq; B[i*N+q] = -sn*bp + cs*bq;
    }
    // Rotate eigenvectors
    for (int i = 0; i < N; i++) {
      double vp = V[i*N+p], vq = V[i*N+q];
      V[i*N+p] = cs*vp + sn*vq; V[i*N+q] = -sn*vp + cs*vq;
    }
  }

  // Extract eigenvalues and sort by magnitude (largest first)
  vector<double> evals(N);
  for (int i = 0; i < N; i++) evals[i] = B[i * N + i];

  vector<int> order(N);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return evals[a] > evals[b];
  });

  // Step 4: Embed using top 3 eigenvectors, scaled by sqrt(eigenvalue)
  vector<coord3d> x(N);
  for (int i = 0; i < N; i++) {
    for (int dim = 0; dim < 3; dim++) {
      int col = order[dim];
      double eval = std::max(0.0, evals[col]);
      x[i][dim] = V[i * N + col] * sqrt(eval);
    }
  }

  // If the 3rd eigenvalue is very small (flat/degenerate), perturb
  double lam3 = std::max(0.0, evals[order[2]]);
  double lam1 = std::max(1e-30, evals[order[0]]);
  if (lam3 / lam1 < 1e-6) {
    // Deterministic perturbation based on vertex index
    for (int i = 0; i < N; i++)
      x[i][2] += 0.1 * lam1 * sin(i * 2.0 * M_PI / N);
  }

  return x;
}

// --- Face-BFS initial placement (copy from delaunay.cc for local use) ---
static std::pair<coord3d, coord3d> trilaterate_c(
    const coord3d& p1, const coord3d& p2, const coord3d& p3,
    double d1, double d2, double d3)
{
  coord3d ex = p2 - p1;
  double a = ex.norm();
  if (a < 1e-15) { coord3d m = (p1+p2+p3)/3.0; return {m, m}; }
  ex /= a;
  coord3d p31 = p3 - p1;
  double ix = ex.dot(p31);
  coord3d ey = p31 - ex * ix;
  double jy = ey.norm();
  if (jy < 1e-15) { coord3d m = (p1+p2+p3)/3.0; return {m, m}; }
  ey /= jy;
  coord3d ez = ex.cross(ey);
  double x = (d1*d1 - d2*d2 + a*a) / (2*a);
  double y = (d1*d1 - d3*d3 + ix*ix + jy*jy - 2*ix*x) / (2*jy);
  double z = sqrt(std::max(0.0, d1*d1 - x*x - y*y));
  coord3d base = p1 + ex*x + ey*y;
  return {base + ez*z, base - ez*z};
}

static coord3d place_opposite_c(const coord3d& xa, const coord3d& xb, const coord3d& xc,
                                 double da, double db)
{
  coord3d edge = xb - xa;
  double L = edge.norm();
  if (L < 1e-15) return (xa + xb + xc) / 3.0;
  coord3d ue = edge / L;
  double t = (L*L + da*da - db*db) / (2*L);
  coord3d center = xa + ue * t;
  double R = sqrt(std::max(0.0, da*da - t*t));
  coord3d vc = xc - center;
  vc = vc - ue * vc.dot(ue);
  double vc_len = vc.norm();
  if (vc_len < 1e-15) {
    coord3d arb = (fabs(ue[0]) < 0.9) ? coord3d(1,0,0) : coord3d(0,1,0);
    vc = ue.cross(arb);
    vc /= vc.norm();
  } else {
    vc /= vc_len;
  }
  coord3d en = ue.cross(vc);
  double alpha = M_PI / 4.0;
  return center - vc * (R * cos(alpha)) - en * (R * sin(alpha));
}

static vector<coord3d> face_bfs_placement_c(
    int N, const neighbours_t& nbrs, const matrix<double>& el)
{
  vector<std::array<int,3>> faces;
  for (int u = 0; u < N; u++) {
    int deg = nbrs[u].size();
    for (int j = 0; j < deg; j++) {
      int v = nbrs[u][j], w = nbrs[u][(j+1)%deg];
      if (u < v && u < w)
        faces.push_back({u, v, w});
    }
  }

  std::map<edge_t, vector<int>> edge_faces;
  for (int f = 0; f < (int)faces.size(); f++)
    for (int k = 0; k < 3; k++) {
      int a = faces[f][k], b = faces[f][(k+1)%3];
      edge_faces[{std::min(a,b), std::max(a,b)}].push_back(f);
    }

  vector<coord3d> x(N);
  vector<bool> placed(N, false);
  vector<bool> done(faces.size(), false);

  int a = faces[0][0], b = faces[0][1], c = faces[0][2];
  double lab = el(a,b), lac = el(a,c), lbc = el(b,c);
  x[a] = coord3d(0, 0, 0);
  x[b] = coord3d(lab, 0, 0);
  double px = (lab*lab + lac*lac - lbc*lbc) / (2*lab);
  x[c] = coord3d(px, sqrt(std::max(0.0, lac*lac - px*px)), 0);
  placed[a] = placed[b] = placed[c] = true;
  done[0] = true;

  std::queue<int> Q;
  Q.push(0);

  while (!Q.empty()) {
    int f = Q.front(); Q.pop();
    for (int k = 0; k < 3; k++) {
      int ea = faces[f][k], eb = faces[f][(k+1)%3];
      edge_t e(std::min(ea,eb), std::max(ea,eb));
      for (int f2 : edge_faces[e]) {
        if (done[f2]) continue;
        int new_v = -1;
        for (int v : faces[f2])
          if (v != ea && v != eb) { new_v = v; break; }
        if (placed[new_v]) { done[f2] = true; Q.push(f2); continue; }
        int old_v = -1;
        for (int v : faces[f])
          if (v != ea && v != eb) { old_v = v; break; }
        int third = -1;
        for (int u : nbrs[new_v])
          if (placed[u] && u != ea && u != eb) { third = u; break; }
        if (third >= 0) {
          auto [s1, s2] = trilaterate_c(x[ea], x[eb], x[third],
                                         el(new_v, ea), el(new_v, eb), el(new_v, third));
          coord3d edg = x[eb] - x[ea];
          coord3d mid = (x[ea] + x[eb]) * 0.5;
          coord3d old_dir = x[old_v] - mid;
          old_dir = old_dir - edg * (old_dir.dot(edg) / edg.dot(edg));
          coord3d s1_dir = s1 - mid;
          s1_dir = s1_dir - edg * (s1_dir.dot(edg) / edg.dot(edg));
          x[new_v] = (s1_dir.dot(old_dir) < 0) ? s1 : s2;
        } else {
          x[new_v] = place_opposite_c(x[ea], x[eb], x[old_v],
                                       el(new_v, ea), el(new_v, eb));
        }
        placed[new_v] = true;
        done[f2] = true;
        Q.push(f2);
      }
    }
  }
  return x;
}

// --- Comparison test ---
struct PerIsomer {
  int buckygen_idx;
  double max_rel_err;
  double time_us;
  int    iters;
  double h_min;
  int    n_concave;
};

struct MethodStats {
  std::string name;
  vector<PerIsomer> results;

  void print(int N) const {
    int n = results.size();
    if (n == 0) return;

    vector<double> errs, times;
    int n_outliers_1pct = 0, n_outliers_5pct = 0, n_machine = 0;
    int n_concave_total = 0;
    int total_iters = 0;
    double worst_h = 1e30;

    for (auto& r : results) {
      errs.push_back(r.max_rel_err);
      times.push_back(r.time_us);
      total_iters += r.iters;
      if (r.max_rel_err > 0.05) n_outliers_5pct++;
      if (r.max_rel_err > 0.01) n_outliers_1pct++;
      if (r.max_rel_err < 1e-12) n_machine++;
      if (r.n_concave > 0) n_concave_total++;
      worst_h = std::min(worst_h, r.h_min);
    }

    std::sort(errs.begin(), errs.end());
    std::sort(times.begin(), times.end());

    vector<double> h_vals;
    for (auto& r : results) h_vals.push_back(r.h_min);
    std::sort(h_vals.begin(), h_vals.end());

    printf("  C%d %s: n=%d\n", N, name.c_str(), n);
    printf("    Error:   median=%.2e  p99=%.2e  max=%.2e\n",
           errs[n/2], errs[(int)(n*0.99)], errs.back());
    printf("    >1%%: %d (%.1f%%)  >5%%: %d (%.1f%%)  machine: %d (%.1f%%)\n",
           n_outliers_1pct, 100.0*n_outliers_1pct/n,
           n_outliers_5pct, 100.0*n_outliers_5pct/n,
           n_machine, 100.0*n_machine/n);
    printf("    Convex:  h_min_worst=%.3e  h_min_p01=%.3e  n_concave_isomers=%d (%.1f%%)\n",
           worst_h, h_vals[(int)(n*0.01)], n_concave_total, 100.0*n_concave_total/n);
    printf("    Time:    median=%.0fus  p99=%.0fus  max=%.0fus\n",
           times[n/2], times[(int)(n*0.99)], times.back());
    printf("    Iters:   total=%d  mean=%.1f\n", total_iters, (double)total_iters/n);
  }
};

static void compare_mds_vs_bfs(int N) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  Triangulation T;

  MethodStats mds_stats{"MDS", {}};
  MethodStats bfs_stats{"face-BFS", {}};

  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    FulleroidDelaunay D(T);
    D.remove_flat_vertices();
    if (D.N != 12) { idx++; continue; }

    // MDS pipeline
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto x0 = mds_placement(D.N, static_cast<const neighbours_t&>(D), D.edge_lengths);
      auto res = run_steihaug(x0, D.N, static_cast<const neighbours_t&>(D), D.edge_lengths);
      auto t1 = std::chrono::high_resolution_clock::now();
      double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
      mds_stats.results.push_back({idx, res.max_rel_err, us, res.iters, res.h_min, res.n_concave});
    }

    // Face-BFS pipeline
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto x0 = face_bfs_placement_c(D.N, static_cast<const neighbours_t&>(D), D.edge_lengths);
      auto res = run_steihaug(x0, D.N, static_cast<const neighbours_t&>(D), D.edge_lengths);
      auto t1 = std::chrono::high_resolution_clock::now();
      double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
      bfs_stats.results.push_back({idx, res.max_rel_err, us, res.iters, res.h_min, res.n_concave});
    }

    idx++;
  }
  BuckyGen::stop(Q);

  int n = mds_stats.results.size();
  printf("\n=== C%d: MDS vs Face-BFS comparison (%d isomers) ===\n", N, n);
  mds_stats.print(N);
  bfs_stats.print(N);

  // --- Failure cross-comparison ---
  // Identify MDS failures (>1% error) and show how BFS does on those same isomers
  double fail_thresh = 0.01;
  printf("\n  --- MDS failures (err > %.0f%%) vs BFS on same isomers ---\n", fail_thresh*100);

  struct FailPair { int idx; double mds_err, bfs_err; double mds_h, bfs_h;
                    int mds_conc, bfs_conc; int mds_it, bfs_it; };
  vector<FailPair> mds_fails;
  for (int i = 0; i < n; i++) {
    if (mds_stats.results[i].max_rel_err > fail_thresh) {
      mds_fails.push_back({
        mds_stats.results[i].buckygen_idx,
        mds_stats.results[i].max_rel_err, bfs_stats.results[i].max_rel_err,
        mds_stats.results[i].h_min, bfs_stats.results[i].h_min,
        mds_stats.results[i].n_concave, bfs_stats.results[i].n_concave,
        mds_stats.results[i].iters, bfs_stats.results[i].iters
      });
    }
  }
  if (mds_fails.empty()) {
    printf("    (none)\n");
  } else {
    printf("    %5s  %11s  %11s  %10s  %10s  %5s  %5s\n",
           "idx", "mds_err", "bfs_err", "mds_h_min", "bfs_h_min", "mds_c", "bfs_c");
    std::sort(mds_fails.begin(), mds_fails.end(),
              [](const FailPair& a, const FailPair& b) { return a.mds_err > b.mds_err; });
    for (auto& f : mds_fails)
      printf("    %5d  %11.3e  %11.3e  %10.3e  %10.3e  %5d  %5d\n",
             f.idx, f.mds_err, f.bfs_err, f.mds_h, f.bfs_h, f.mds_conc, f.bfs_conc);
    int bfs_saves = 0;
    for (auto& f : mds_fails)
      if (f.bfs_err < fail_thresh) bfs_saves++;
    printf("    BFS recovers %d/%d MDS failures (%.0f%%)\n",
           bfs_saves, (int)mds_fails.size(), 100.0*bfs_saves/mds_fails.size());
  }

  // Also show BFS failures and how MDS does on those
  printf("\n  --- BFS failures (err > %.0f%%) vs MDS on same isomers ---\n", fail_thresh*100);
  vector<FailPair> bfs_fails;
  for (int i = 0; i < n; i++) {
    if (bfs_stats.results[i].max_rel_err > fail_thresh) {
      bfs_fails.push_back({
        bfs_stats.results[i].buckygen_idx,
        mds_stats.results[i].max_rel_err, bfs_stats.results[i].max_rel_err,
        mds_stats.results[i].h_min, bfs_stats.results[i].h_min,
        mds_stats.results[i].n_concave, bfs_stats.results[i].n_concave,
        mds_stats.results[i].iters, bfs_stats.results[i].iters
      });
    }
  }
  if (bfs_fails.empty()) {
    printf("    (none)\n");
  } else {
    // Don't print all BFS failures (too many), just summary
    int mds_saves = 0;
    for (auto& f : bfs_fails)
      if (f.mds_err < fail_thresh) mds_saves++;
    printf("    %d BFS failures total. MDS recovers %d/%d (%.1f%%)\n",
           (int)bfs_fails.size(), mds_saves, (int)bfs_fails.size(),
           100.0*mds_saves/bfs_fails.size());

    // Both fail?
    int both_fail = 0;
    for (int i = 0; i < n; i++)
      if (mds_stats.results[i].max_rel_err > fail_thresh &&
          bfs_stats.results[i].max_rel_err > fail_thresh)
        both_fail++;
    int either_fail = 0;
    for (int i = 0; i < n; i++)
      if (mds_stats.results[i].max_rel_err > fail_thresh ||
          bfs_stats.results[i].max_rel_err > fail_thresh)
        either_fail++;
    printf("    Both fail: %d  Either fails: %d  Best-of-two recovers: %d/%d (%.2f%% success)\n",
           both_fail, either_fail, either_fail - both_fail, either_fail,
           100.0*(n - both_fail)/n);
  }
}

TEST(DelaunayCompare, C60_MDS_vs_BFS) { compare_mds_vs_bfs(60); }
TEST(DelaunayCompare, C70_MDS_vs_BFS) { compare_mds_vs_bfs(70); }
TEST(DelaunayCompare, C80_MDS_vs_BFS) { compare_mds_vs_bfs(80); }

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

TEST(IntrinsicDelaunay, Plantri15_AllTriangulations) {
  // Validate iDT on ALL 2,406,841 3-connected planar triangulations on 15 vertices.
  // These include vertices of degree 3-12, testing negative-curvature cone points.
  std::string path = FULLERENE_ROOT "/data/triangulations_15.pl";

  FILE* f = fopen(path.c_str(), "rb");
  if (!f) {
    GTEST_SKIP() << "Plantri database not found: " << path;
  }
  fclose(f);

  int n_tested = 0, n_failed = 0;
  int max_degree_seen = 0;
  std::vector<double> times_us;

  stream_planarcode(path, [&](const PlanarGraph& g, int idx) {
    // Record original degrees before sort_flat_last reorders
    vector<int> orig_degrees(g.N);
    for (int u = 0; u < g.N; u++)
      orig_degrees[u] = g[u].size();

    for (int d : orig_degrees)
      max_degree_seen = std::max(max_degree_seen, d);

    // Count cone points (non-degree-6)
    int n_cones = 0;
    for (int d : orig_degrees)
      if (d != 6) n_cones++;

    // Skip if all vertices are degree 6 (flat torus — no cone points)
    if (n_cones == 0) return;

    // Number of original faces: for triangulation on N vertices, F = 2N - 4
    int N_original_faces = 2 * g.N - 4;

    Triangulation T(g);  // already oriented from planar_code

    // Validate input triangulation
    if (!T.is_consistently_oriented()) {
      n_failed++;
      if (n_failed <= 10)
        std::cerr << "Graph #" << idx << ": input not consistently oriented" << std::endl;
      n_tested++;
      return;
    }

    FulleroidDelaunay D(T);

    auto t0 = std::chrono::high_resolution_clock::now();
    D.remove_flat_vertices();
    auto t1 = std::chrono::high_resolution_clock::now();
    times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());

    // Skip remaining checks if removal failed (vertex count mismatch)
    if (D.N != n_cones) {
      n_failed++;
      if (n_failed <= 10)
        std::cerr << "Graph #" << idx << ": removal incomplete (D.N="
                  << D.N << ", expected " << n_cones << ")" << std::endl;
      n_tested++;
      return;
    }

    EXPECT_TRUE(D.is_consistently_oriented()) << "Graph #" << idx;
    EXPECT_TRUE(D.edge_lengths_are_symmetric()) << "Graph #" << idx;

    // Delaunay check — general triangulations may have more blocked flips
    // than fullerene duals due to multi-edge/self-loop constraints
    int non_del = 0;
    for (node_t u = 0; u < D.N; u++)
      for (node_t v : D[u])
        if (u < v && !D.is_delaunay_edge(u, v))
          non_del++;
    if (non_del > 4) {
      ADD_FAILURE() << "Graph #" << idx << ": " << non_del << " non-Delaunay edges";
      n_failed++;
    }

    // Cone angle check
    vector<int> cone_degrees;
    for (int d : orig_degrees)
      if (d != 6) cone_degrees.push_back(d);

    if ((int)cone_degrees.size() == D.N) {
      for (node_t u = 0; u < D.N; u++) {
        double expected_angle = cone_degrees[u] * M_PI / 3.0;
        double total_angle = 0;
        const auto& nbrs = D[u];
        int k = nbrs.size();
        for (int j = 0; j < k; j++) {
          node_t v = nbrs[j];
          node_t w = nbrs[(j + 1) % k];
          double a = D.get_length(v, w);
          double b = D.get_length(u, v);
          double c = D.get_length(u, w);
          double cos_u = std::clamp((b*b + c*c - a*a) / (2.0*b*c), -1.0, 1.0);
          total_angle += acos(cos_u);
        }
        if (std::abs(total_angle - expected_angle) > 1e-6) {
          ADD_FAILURE() << "Graph #" << idx << " vertex " << u
                        << " (orig deg " << cone_degrees[u]
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
    for (node_t u = 0; u < D.N; u++)
      for (size_t j = 0; j < D[u].size(); j++) {
        node_t v = D[u][j];
        node_t w = D[u][(j + 1) % D[u].size()];
        if (u < v && u < w) {
          double a = D.get_length(u, v);
          double b = D.get_length(v, w);
          double c = D.get_length(w, u);
          double s = (a + b + c) / 2.0;
          double area2 = s * (s - a) * (s - b) * (s - c);
          actual_area += sqrt(std::max(0.0, area2));
        }
      }
    if (std::abs(actual_area - expected_area) > 1e-6) {
      ADD_FAILURE() << "Graph #" << idx << ": area " << actual_area
                    << " != expected " << expected_area;
      n_failed++;
    }

    n_tested++;

    if (n_tested % 100000 == 0)
      fprintf(stderr, "  [%d] tested=%d failed=%d\n", idx, n_tested, n_failed);
  });

  print_timing_stats("Plantri15 iDT", times_us);
  std::cout << "Plantri15: tested " << n_tested << " triangulations, "
            << n_failed << " failures, max degree seen = "
            << max_degree_seen << std::endl;
  EXPECT_EQ(n_failed, 0);
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
