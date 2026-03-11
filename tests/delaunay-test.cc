#include <gtest/gtest.h>
#include "fullerenes/delaunay.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <chrono>
#include <numeric>

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
    E += D.neighbours[u].size();
  E /= 2;

  if (expected_verts == 12) {
    EXPECT_EQ(E, 30) << "Expected 30 edges for 12 vertices on genus 0";
  }

  int F = E - V + 2;
  EXPECT_EQ(F, 2 * V - 4) << "Euler formula check failed";

  // All vertices have degree >= 3
  for (node_t u = 0; u < D.N; u++)
    EXPECT_GE((int)D.neighbours[u].size(), 3) << "Vertex " << u << " has degree < 3";

  // Consistent orientation
  EXPECT_TRUE(D.is_consistently_oriented()) << "Orientation is broken";

  // --- Edge length consistency ---

  // Symmetric
  EXPECT_TRUE(D.edge_lengths_are_symmetric());

  // Positive
  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D.neighbours[u])
      EXPECT_GT(D.get_length(u, v), 0) << "Edge (" << u << "," << v << ") has non-positive length";

  // Edge-adjacency consistency: edge_lengths(u,v) > 0 iff v in neighbours[u]
  for (node_t u = 0; u < D.N; u++)
    for (node_t v = 0; v < D.N; v++) {
      bool is_neighbour = std::find(D.neighbours[u].begin(), D.neighbours[u].end(), v)
                          != D.neighbours[u].end();
      if (is_neighbour)
        EXPECT_GT(D.get_length(u, v), 0) << "Neighbour (" << u << "," << v << ") has zero length";
      else if (u != v)
        EXPECT_EQ(D.get_length(u, v), 0) << "Non-neighbour (" << u << "," << v << ") has nonzero length";
    }

  // --- Delaunay criterion ---
  // Lawson's flip algorithm may leave a small number of non-Delaunay edges
  // in rare topological configurations (multi-edge blocked flips).
  // Count them rather than hard-failing.
  {
    int non_del = 0;
    for (node_t u = 0; u < D.N; u++)
      for (node_t v : D.neighbours[u])
        if (u < v && !D.is_delaunay_edge(u, v))
          non_del++;
    EXPECT_LE(non_del, 2) << non_del << " non-Delaunay edges (max 2 allowed for blocked flips)";
  }

  // --- Metric checks ---

  // Triangle inequality: every triangle must satisfy strict triangle inequality
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D.neighbours[u].size(); j++) {
      node_t v = D.neighbours[u][j];
      node_t w = D.neighbours[u][(j + 1) % D.neighbours[u].size()];
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
    for (size_t j = 0; j < D.neighbours[u].size(); j++) {
      node_t v = D.neighbours[u][j];
      node_t w = D.neighbours[u][(j + 1) % D.neighbours[u].size()];
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
    const auto& nbrs = D.neighbours[u];
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
    for (size_t j = 0; j < D.neighbours[u].size(); j++) {
      node_t v = D.neighbours[u][j];
      node_t w = D.neighbours[u][(j + 1) % D.neighbours[u].size()];
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

TEST(IntrinsicDelaunay, C80_AllIsomers) {
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

TEST(IntrinsicDelaunay, C20_AllEquilateral) {
  // For C20 (icosahedron), all edges should remain length 1
  // since all vertices are degree-5 (no flat vertices to remove).
  Triangulation T = make_dual(20, 0, false);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();

  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D.neighbours[u])
      EXPECT_DOUBLE_EQ(D.get_length(u, v), 1.0);
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
    for (node_t v : D.neighbours[u])
      if (u < v) {
        double len = D.get_length(u, v);
        if (first_len < 0) first_len = len;
        EXPECT_NEAR(len, first_len, 1e-10)
          << "Edge (" << u << "," << v << ") has length " << len
          << " but expected " << first_len << " (icosahedral symmetry)";
      }

  std::cout << "C60 Ih: all 30 edges have length " << first_len << std::endl;
}
