// Round-trip and edge-case tests for PlanarGraph::to_ascii / from_ascii.
//
// to_ascii writes the neighbour list as nested JSON arrays:
//     [[1,2,3],[0,2,3],...]
// from_ascii is the inverse, reading one record from the current stream
// position. The pair is the canonical ASCII oriented-graph dump used by
// the claude-projects validator harness.

#include "fullerenes/planargraph.hh"
#include "fullerenes/graph.hh"

#include <gtest/gtest.h>
#include <cstdio>
#include <stdexcept>

namespace {

bool graphs_equal(const PlanarGraph& A, const PlanarGraph& B) {
  if (A.N != B.N) return false;
  for (node_t u = 0; u < A.N; u++) {
    auto a = A[u], b = B[u];
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
      if (a[i] != b[i]) return false;
    }
  }
  return true;
}

PlanarGraph roundtrip(const PlanarGraph& G) {
  FILE* f = std::tmpfile();
  EXPECT_NE(f, nullptr);
  EXPECT_TRUE(PlanarGraph::to_ascii(G, f));
  std::rewind(f);
  PlanarGraph H = PlanarGraph::from_ascii(f);
  std::fclose(f);
  return H;
}

TEST(PlanarGraphAscii, EmptyGraph) {
  PlanarGraph G(0);
  PlanarGraph H = roundtrip(G);
  EXPECT_EQ(H.N, 0);
}

TEST(PlanarGraphAscii, K3) {
  PlanarGraph G(Graph{{1,2}, {0,2}, {0,1}});
  PlanarGraph H = roundtrip(G);
  EXPECT_TRUE(graphs_equal(G, H));
}

TEST(PlanarGraphAscii, K4) {
  PlanarGraph G(Graph{{1,2,3}, {0,2,3}, {0,1,3}, {0,1,2}});
  PlanarGraph H = roundtrip(G);
  EXPECT_TRUE(graphs_equal(G, H));
}

// CCW orientation is a permutation of the neighbour list. Make sure the
// exact order survives the round-trip rather than being canonicalised.
TEST(PlanarGraphAscii, OrientationPreserved) {
  PlanarGraph G(Graph{{2,3,1}, {3,2,0}, {0,1,3}, {1,0,2}});
  PlanarGraph H = roundtrip(G);
  EXPECT_TRUE(graphs_equal(G, H));
  EXPECT_EQ(int(H[0][0]), 2);
  EXPECT_EQ(int(H[0][1]), 3);
  EXPECT_EQ(int(H[0][2]), 1);
}

// A dodecahedron-shaped neighbour list: 20 vertices, every row of degree 3.
// Real-fullerene-sized record going through the parser.
TEST(PlanarGraphAscii, DodecahedronRoundTrip) {
  PlanarGraph G(Graph{
    { 1, 4, 5}, { 0, 2, 6}, { 1, 3, 7}, { 2, 4, 8}, { 0, 3, 9},
    { 0,10,14}, { 1,10,11}, { 2,11,12}, { 3,12,13}, { 4,13,14},
    { 5, 6,15}, { 6, 7,16}, { 7, 8,17}, { 8, 9,18}, { 9, 5,19},
    {10,16,19}, {11,15,17}, {12,16,18}, {13,17,19}, {14,15,18},
  });
  PlanarGraph H = roundtrip(G);
  EXPECT_TRUE(graphs_equal(G, H));
  EXPECT_EQ(H.N, 20);
}

// Whitespace-tolerant parsing.
TEST(PlanarGraphAscii, AcceptsWhitespace) {
  FILE* f = std::tmpfile();
  std::fputs("  [ [1, 2 , 3] , [0,2,3] , [0,1,3] , [0,1,2] ]  ", f);
  std::rewind(f);
  PlanarGraph H = PlanarGraph::from_ascii(f);
  std::fclose(f);
  PlanarGraph K(Graph{{1,2,3}, {0,2,3}, {0,1,3}, {0,1,2}});
  EXPECT_TRUE(graphs_equal(K, H));
}

// Read two graphs in sequence from one stream (newline-separated).
TEST(PlanarGraphAscii, MultiRecordStreaming) {
  PlanarGraph G1(Graph{{1,2}, {0,2}, {0,1}});
  PlanarGraph G2(Graph{{1,2,3}, {0,2,3}, {0,1,3}, {0,1,2}});

  FILE* f = std::tmpfile();
  EXPECT_TRUE(PlanarGraph::to_ascii(G1, f));
  std::fputc('\n', f);
  EXPECT_TRUE(PlanarGraph::to_ascii(G2, f));
  std::fputc('\n', f);
  std::rewind(f);

  PlanarGraph H1 = PlanarGraph::from_ascii(f);
  PlanarGraph H2 = PlanarGraph::from_ascii(f);
  std::fclose(f);

  EXPECT_TRUE(graphs_equal(G1, H1));
  EXPECT_TRUE(graphs_equal(G2, H2));
}

TEST(PlanarGraphAscii, EofBeforeRecordThrows) {
  FILE* f = std::tmpfile();
  std::rewind(f);
  EXPECT_THROW(PlanarGraph::from_ascii(f), std::runtime_error);
  std::fclose(f);
}

TEST(PlanarGraphAscii, MalformedHeaderThrows) {
  FILE* f = std::tmpfile();
  std::fputs("not a graph", f);
  std::rewind(f);
  EXPECT_THROW(PlanarGraph::from_ascii(f), std::runtime_error);
  std::fclose(f);
}

TEST(PlanarGraphAscii, UnterminatedRecordThrows) {
  FILE* f = std::tmpfile();
  std::fputs("[[1,2,3", f);
  std::rewind(f);
  EXPECT_THROW(PlanarGraph::from_ascii(f), std::runtime_error);
  std::fclose(f);
}

TEST(PlanarGraphAscii, BadSeparatorThrows) {
  FILE* f = std::tmpfile();
  std::fputs("[[1,2,3];[0,2,3]]", f);
  std::rewind(f);
  EXPECT_THROW(PlanarGraph::from_ascii(f), std::runtime_error);
  std::fclose(f);
}

}  // namespace
