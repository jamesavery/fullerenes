// The compact binary geometry format (.geo, GEO-FORMAT.md): the worked example byte for
// byte, round trips for every coordinate type, capacities, appending with crash
// leftovers, the checksum, the reader's rejections (including every twin block of the
// tetrahedron), file-name options, and non-simplicial delta-complexes through
// DelaunayTriangulation.

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/geo-format.hh"
#include "fullerenes/polyhedron.hh"

#include <gtest/gtest.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <vector>

#include <stdlib.h>
#include <unistd.h>

using Code  = mesh_io_error::Code;
using Rows  = std::vector<std::vector<int>>;
using Bytes = std::vector<uint8_t>;

namespace {

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

// A polyhedron from points around the origin and its edges; each row is sorted
// counter-clockwise about the vertex's outward direction, starting at the first
// neighbour in edge order.
Polyhedron star_polyhedron(const std::vector<coord3d>& X, const std::vector<std::pair<int,int>>& edges) {
  const int N = X.size();
  Rows rows(N);
  for (auto [u, v] : edges) { rows[u].push_back(v); rows[v].push_back(u); }
  for (int u = 0; u < N; u++) {
    const int first = rows[u][0];
    const coord3d n = X[u] / X[u].norm();
    coord3d e1 = X[first] - X[u];
    e1 -= n * e1.dot(n);
    e1 /= e1.norm();
    const coord3d e2 = n.cross(e1);
    const auto angle = [&](int v) {
      if (v == first) return 0.0;
      const coord3d d = X[v] - X[u];
      const double a = atan2(d.dot(e2), d.dot(e1));
      return a < 0 ? a + 2 * M_PI : a;
    };
    std::sort(rows[u].begin(), rows[u].end(), [&](int a, int b) { return angle(a) < angle(b); });
  }
  Graph G(size_t(N), 10);
  for (int u = 0; u < N; u++)
    for (int v : rows[u]) G.push_back(u, v);
  return Polyhedron(PlanarGraphView(G.N, G.dmax, G.neighbours, G.deg), X);
}

Polyhedron tetrahedron() {
  return star_polyhedron({{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}},
                         {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}});
}

// Degrees 4, 4, 4, 3, 3.
Polyhedron bipyramid() {
  std::vector<coord3d> X;
  for (int i = 0; i < 3; i++) X.push_back({cos(2 * M_PI * i / 3), sin(2 * M_PI * i / 3), 0.0});
  X.push_back({0, 0, 1.5});
  X.push_back({0, 0, -1.5});
  return star_polyhedron(X, {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}, {0, 4}, {1, 4}, {2, 4}});
}

Polyhedron scaled(Polyhedron P, double s, coord3d shift = {0, 0, 0}) {
  for (auto& p : P.points) p = p * s + shift;
  return P;
}

Rows rows_of(const GraphView& G) {
  Rows rows(G.N);
  for (int u = 0; u < G.N; u++) {
    const auto r = G.nbrs(u);
    rows[u].assign(r.begin(), r.end());
  }
  return rows;
}

// The rows of a decoded record (GEO-FORMAT.md sec. 7.6).
Rows rows_of(const geo_record& r) {
  const auto off = geo::row_offsets(r.degree);
  std::vector<int> origin(off[r.n]);
  for (int v = 0; v < r.n; v++) std::fill(origin.begin() + off[v], origin.begin() + off[v + 1], v);
  Rows rows(r.n);
  for (int v = 0; v < r.n; v++)
    for (int h = off[v]; h < off[v + 1]; h++) rows[v].push_back(origin[r.twin[h]]);
  return rows;
}

geo_options fixed_options(int width) {
  geo_options opt;
  opt.type  = geo_type::FIXED;
  opt.width = uint8_t(width);
  return opt;
}

geo_options declared(geo_options opt, int deg_min, int deg_bits) {
  opt.deg_min  = deg_min;
  opt.deg_bits = deg_bits;
  return opt;
}

// ---------------------------------------------------------------------------
// Raw file access
// ---------------------------------------------------------------------------

FILE* temp() {
  FILE* f = std::tmpfile();
  if (!f) throw std::runtime_error("tmpfile failed");
  return f;
}

Bytes bytes_of(FILE* f) {
  fflush(f);
  fseeko(f, 0, SEEK_END);
  Bytes b(size_t(ftello(f)));
  fseeko(f, 0, SEEK_SET);
  if (fread(b.data(), 1, b.size(), f) != b.size()) throw std::runtime_error("short read");
  return b;
}

FILE* file_with(const Bytes& b) {
  FILE* f = temp();
  if (fwrite(b.data(), 1, b.size(), f) != b.size()) throw std::runtime_error("short write");
  fflush(f);
  return f;
}

void expect_points(std::span<const coord3d> got, std::span<const coord3d> want) {
  ASSERT_EQ(got.size(), want.size());
  for (size_t v = 0; v < want.size(); v++)
    for (int k = 0; k < 3; k++) EXPECT_EQ(got[v][k], want[v][k]) << "vertex " << v << ", axis " << k;
}

// Exactly `code`, from `call`.
template <class F>
void expect_code(Code code, F&& call) {
  try {
    call();
    ADD_FAILURE() << "expected mesh_io_error code " << int(code) << ", nothing was thrown";
  } catch (const mesh_io_error& e) {
    EXPECT_EQ(int(e.code), int(code)) << e.what();
  }
}

double record_scale_at(const Bytes& b, size_t offset) {
  uint32_t v = 0;
  for (int k = 0; k < 4; k++) v |= uint32_t(b[offset + k]) << (8 * k);
  return double(std::bit_cast<float>(v));
}

// GEO-FORMAT.md sec. 12.
const Bytes worked_example = {
  0x10, 0x04, 0x00, 0x00, 0x08, 0x05, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x6e, 0xdf, 0xb3, 0x90, 0x3f, 0x6b, 0x5f, 0xcb,
  0x05, 0x02, 0x01, 0x3c, 0xff, 0xff, 0xff, 0xff, 0x01, 0x01, 0x01, 0xff, 0x01, 0x01, 0x01, 0xff,
  0x63, 0xb9, 0xa7, 0x00,
};

geo_options worked_example_options() {
  geo_options opt = fixed_options(8);
  opt.triangulation = true;
  return opt;
}

}  // namespace

// ===========================================================================
// Writing and reading
// ===========================================================================

TEST(GeoFormat, WorkedExampleByteForByte) {
  const Polyhedron P = tetrahedron();
  ASSERT_EQ(rows_of(P), (Rows{{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}}));

  FILE* f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(P, f, false, worked_example_options()));
  EXPECT_EQ(bytes_of(f), worked_example);
  EXPECT_TRUE(Polyhedron::verify_geo(f));

  const Polyhedron Q = Polyhedron::from_geo(f);
  EXPECT_EQ(rows_of(Q), rows_of(P));
  const double sigma = double(std::bit_cast<float>(uint32_t(0x3c010205)));
  for (int v = 0; v < 4; v++)
    for (int k = 0; k < 3; k++) {
      EXPECT_EQ(Q.points[v][k], (P.points[v][k] > 0 ? 127 : -127) * sigma);
      EXPECT_LE(std::fabs(Q.points[v][k] - P.points[v][k]), sigma / 2);
    }
  fclose(f);
}

TEST(GeoFormat, FloatRoundTrips) {
  const Polyhedron P = Polyhedron::C20();
  for (geo_type type : {geo_type::F64, geo_type::F32}) {
    geo_options opt;
    opt.type = type;
    FILE* f = temp();
    ASSERT_TRUE(Polyhedron::to_geo(P, f, false, opt));
    const geo_header H = Polyhedron::read_geo_header(f);
    EXPECT_EQ(H.count, 1u);
    EXPECT_EQ(H.opt.N, 20u);
    EXPECT_EQ(H.opt.deg_min, 3);
    EXPECT_EQ(H.opt.deg_bits, 0);
    EXPECT_EQ(bytes_of(f).size(), 32 + H.record_size());
    EXPECT_EQ(H.record_size() % (type == geo_type::F64 ? 8 : 4), 0u);

    const Polyhedron Q = Polyhedron::from_geo(f);
    EXPECT_EQ(rows_of(Q), rows_of(P));
    for (int v = 0; v < P.N; v++)
      for (int k = 0; k < 3; k++)
        EXPECT_EQ(Q.points[v][k], type == geo_type::F64 ? P.points[v][k] : double(float(P.points[v][k])));
    fclose(f);
  }
}

TEST(GeoFormat, FixedPointErrorBound) {
  const Polyhedron P = scaled(Polyhedron::C20(), 3.7, {100.25, -50.5, 3.0});
  for (int width : {2, 3, 5, 8, 12, 16, 21, 26, 30})
    for (bool offset : {false, true}) {
      geo_options opt = fixed_options(width);
      opt.offset = offset;
      FILE* f = temp();
      ASSERT_TRUE(Polyhedron::to_geo(P, f, false, opt));
      const double sigma = record_scale_at(bytes_of(f), 32);
      const Polyhedron Q = Polyhedron::from_geo(f);
      EXPECT_EQ(rows_of(Q), rows_of(P));

      // The per-record scale is the smallest f32 covering the extent.
      const double q = double((uint64_t(1) << (width - 1)) - 1);
      double s = 0;
      std::vector<double> o(3, 0.0);
      if (offset)
        for (int k = 0; k < 3; k++) {
          double lo = INFINITY, hi = -INFINITY;
          for (const auto& p : P.points) { lo = std::min(lo, p[k]); hi = std::max(hi, p[k]); }
          o[k] = double(float(lo / 2 + hi / 2));
        }
      for (const auto& p : P.points)
        for (int k = 0; k < 3; k++) s = std::max(s, std::fabs(p[k] - o[k]));
      EXPECT_GE(sigma, s / q);
      EXPECT_LT(double(std::nextafter(float(sigma), 0.0f)), s / q);

      for (int v = 0; v < P.N; v++)
        for (int k = 0; k < 3; k++)
          EXPECT_LE(std::fabs(Q.points[v][k] - P.points[v][k]), sigma / 2 * (1 + 1e-9) + 1e-12)
              << "width " << width << " offset " << offset;
      fclose(f);
    }
}

TEST(GeoFormat, FileScaleRange) {
  const Polyhedron P = Polyhedron::C20();
  double m = 0;
  for (const auto& p : P.points)
    for (int k = 0; k < 3; k++) m = std::max(m, std::fabs(p[k]));

  geo_options fits = fixed_options(8);
  fits.scale = float(m / 127 * 1.01);
  FILE* f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(P, f, false, fits));
  const Polyhedron Q = Polyhedron::from_geo(f);
  for (int v = 0; v < P.N; v++)
    for (int k = 0; k < 3; k++)
      EXPECT_LE(std::fabs(Q.points[v][k] - P.points[v][k]), double(fits.scale) / 2 * (1 + 1e-9));
  EXPECT_EQ(Polyhedron::read_geo_header(f).opt.scale, fits.scale);
  fclose(f);

  geo_options narrow = fixed_options(8);
  narrow.scale = float(m / 127 * 0.9);
  f = temp();
  expect_code(Code::ValueOutOfRange, [&] { Polyhedron::to_geo(P, f, false, narrow); });
  fclose(f);
}

TEST(GeoFormat, AppendAndRandomAccess) {
  const geo_options opt = declared({}, 3, 0);
  FILE* f = temp();
  std::vector<Polyhedron> Ps;
  for (int i = 0; i < 5; i++) {
    Ps.push_back(scaled(Polyhedron::C20(), 1 + i / 10.0));
    ASSERT_TRUE(Polyhedron::to_geo(Ps.back(), f, true, opt));
  }
  EXPECT_EQ(Polyhedron::read_geo_header(f).count, 5u);
  EXPECT_TRUE(Polyhedron::verify_geo(f));
  for (int i : {3, 0, 4, 1, 2}) {
    const Polyhedron Q = Polyhedron::from_geo(f, i);
    expect_points(Q.points, Ps[i].points);
  }
  expect_code(Code::IndexOutOfRange, [&] { Polyhedron::from_geo(f, 5); });

  // A single-polyhedron file is appendable with its derived capacities.
  FILE* g = temp();
  ASSERT_TRUE(Polyhedron::to_geo(Ps[0], g));
  ASSERT_TRUE(Polyhedron::to_geo(Ps[1], g, true));
  EXPECT_EQ(Polyhedron::read_geo_header(g).count, 2u);
  EXPECT_TRUE(Polyhedron::verify_geo(g));
  fclose(g);
  fclose(f);
}

TEST(GeoFormat, AppendRules) {
  const Polyhedron C20 = Polyhedron::C20();

  FILE* f = temp();
  EXPECT_THROW(Polyhedron::to_geo(C20, f, true, geo_options{}), std::invalid_argument)
      << "a file created by appending must declare its degree range";
  geo_options capacity = declared({}, 3, 0);
  capacity.record_n = true;
  EXPECT_THROW(Polyhedron::to_geo(C20, f, true, capacity), std::invalid_argument)
      << "record_n needs a declared N";
  EXPECT_EQ(bytes_of(f).size(), 0u) << "a refused append leaves the file untouched";

  ASSERT_TRUE(Polyhedron::to_geo(C20, f, true, declared({}, 3, 0)));
  expect_code(Code::HeaderMismatch, [&] { Polyhedron::to_geo(C20, f, true, declared(fixed_options(12), 3, 0)); });
  expect_code(Code::HeaderMismatch, [&] { Polyhedron::to_geo(C20, f, true, declared({}, 3, 1)); });
  geo_options other_N = declared({}, 3, 0);
  other_N.N = 21;
  expect_code(Code::HeaderMismatch, [&] { Polyhedron::to_geo(C20, f, true, other_N); });
  expect_code(Code::CapacityExceeded, [&] { Polyhedron::to_geo(tetrahedron(), f, true, {}); });
  EXPECT_EQ(Polyhedron::read_geo_header(f).count, 1u);
  fclose(f);

  f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(tetrahedron(), f, true, declared({}, 3, 0)));
  geo_options four = declared({}, 3, 0);
  four.N = 4;
  expect_code(Code::CapacityExceeded, [&] { Polyhedron::to_geo(bipyramid(), f, true, four); });
  fclose(f);

  // Streams: append mode and read-only are refused.
  char path[] = "/tmp/geo-io-test-XXXXXX";
  const int fd = mkstemp(path);
  ASSERT_GE(fd, 0);
  close(fd);
  FILE* a = fopen(path, "a+b");
  EXPECT_THROW(Polyhedron::to_geo(C20, a, true, declared({}, 3, 0)), std::invalid_argument);
  fclose(a);
  FILE* r = fopen(path, "rb");
  EXPECT_THROW(Polyhedron::to_geo(C20, r, true, declared({}, 3, 0)), std::invalid_argument);
  EXPECT_THROW(Polyhedron::to_geo(C20, r), std::invalid_argument);
  fclose(r);
  unlink(path);

  expect_code(Code::NullFile, [&] { Polyhedron::from_geo(nullptr); });
}

TEST(GeoFormat, VertexCountPerRecordAndGraphOnly) {
  const Polyhedron T = tetrahedron(), B = bipyramid();

  geo_options graph_only = declared({}, 3, 1);
  graph_only.type = geo_type::NONE;
  graph_only.record_n = true;
  graph_only.triangulation = true;
  graph_only.N = 5;
  FILE* f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(T, f, true, graph_only));
  ASSERT_TRUE(Polyhedron::to_geo(B, f, true, graph_only));
  expect_code(Code::CapacityExceeded, [&] { Polyhedron::to_geo(scaled(Polyhedron::C20(), 1), f, true, graph_only); });

  const geo_record r0 = geo::read_record(f, 0), r1 = geo::read_record(f, 1);
  EXPECT_EQ(r0.n, 4);
  EXPECT_EQ(r1.n, 5);
  EXPECT_TRUE(r0.x.empty());
  EXPECT_EQ(rows_of(r0), rows_of(T));
  EXPECT_EQ(rows_of(r1), rows_of(B));
  expect_code(Code::UnsupportedFormat, [&] { Polyhedron::from_geo(f, 0); });
  EXPECT_TRUE(Polyhedron::verify_geo(f));
  fclose(f);

  geo_options floats = declared({}, 3, 1);
  floats.type = geo_type::F32;
  floats.record_n = true;
  floats.N = 5;
  f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(T, f, true, floats));
  ASSERT_TRUE(Polyhedron::to_geo(B, f, true, floats));
  for (int i = 0; i < 2; i++) {
    const Polyhedron& P = i == 0 ? T : B;
    const Polyhedron Q = Polyhedron::from_geo(f, i);
    EXPECT_EQ(rows_of(Q), rows_of(P));
    for (int v = 0; v < P.N; v++)
      for (int k = 0; k < 3; k++) EXPECT_EQ(Q.points[v][k], double(float(P.points[v][k])));
  }

  // Unused slots must be zero: tetrahedron record, vertex slot 4 (record bytes 48..59),
  // then the bit stream at record byte 60: n (3 bits), degree slots 0..4 (1 bit each).
  const Bytes good = bytes_of(f);
  EXPECT_EQ(good[32 + 60], 0x04);
  Bytes bad = good;
  bad[32 + 48] = 1;
  FILE* g = file_with(bad);
  expect_code(Code::MalformedFile, [&] { Polyhedron::from_geo(g, 0); });
  fclose(g);
  bad = good;
  bad[32 + 60] = 0x84;
  g = file_with(bad);
  expect_code(Code::MalformedFile, [&] { Polyhedron::from_geo(g, 0); });
  fclose(g);
  bad = good;
  bad[32 + 60] = 0x07;
  g = file_with(bad);
  expect_code(Code::MalformedFile, [&] { Polyhedron::from_geo(g, 0); });
  fclose(g);
  fclose(f);
}

TEST(GeoFormat, CoordinatesOnly) {
  const Polyhedron P = Polyhedron::C20();
  geo_options opt = fixed_options(16);
  opt.graph = false;
  FILE* f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(P, f, false, opt));
  EXPECT_EQ(Polyhedron::read_geo_header(f).record_size(), 4u + 3 * 20 * 2);
  expect_code(Code::UnsupportedFormat, [&] { Polyhedron::from_geo(f); });

  const Polyhedron Q = Polyhedron::from_geo(f, static_cast<const PlanarGraphView&>(P));
  EXPECT_EQ(rows_of(Q), rows_of(P));
  for (int v = 0; v < P.N; v++)
    for (int k = 0; k < 3; k++) EXPECT_NEAR(Q.points[v][k], P.points[v][k], 1e-4);
  const Polyhedron T = tetrahedron();
  EXPECT_THROW(Polyhedron::from_geo(f, static_cast<const PlanarGraphView&>(T)), std::invalid_argument);
  fclose(f);

  f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(P, f));
  EXPECT_THROW(Polyhedron::from_geo(f, static_cast<const PlanarGraphView&>(P)), std::invalid_argument);
  fclose(f);
}

TEST(GeoFormat, BatchWriterDerivesCapacities) {
  const std::vector<Polyhedron> Ps = {tetrahedron(), bipyramid(), scaled(tetrahedron(), 2)};
  geo_options opt;
  opt.record_n = true;
  FILE* f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(Ps, f, opt));
  const geo_header H = Polyhedron::read_geo_header(f);
  EXPECT_EQ(H.count, 3u);
  EXPECT_EQ(H.opt.N, 5u);
  EXPECT_EQ(H.opt.deg_min, 3);
  EXPECT_EQ(H.opt.deg_bits, 1);
  EXPECT_EQ(H.edge_capacity(), 9u);
  EXPECT_TRUE(Polyhedron::verify_geo(f));
  for (int i = 0; i < 3; i++) {
    const Polyhedron Q = Polyhedron::from_geo(f, i);
    EXPECT_EQ(rows_of(Q), rows_of(Ps[i]));
    expect_points(Q.points, Ps[i].points);
  }
  fclose(f);

  f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(Polyhedron::C20(), f));
  const Bytes before = bytes_of(f);
  expect_code(Code::CapacityExceeded, [&] { Polyhedron::to_geo(Ps, f, geo_options{}); });
  EXPECT_EQ(bytes_of(f), before) << "a refused write leaves the file as it was";
  EXPECT_THROW(Polyhedron::to_geo(std::span<const Polyhedron>{}, f, geo_options{}), std::invalid_argument);
  geo_options empty = declared({}, 3, 0);
  empty.N = 20;
  ASSERT_TRUE(Polyhedron::to_geo(std::span<const Polyhedron>{}, f, empty));
  EXPECT_EQ(bytes_of(f).size(), 32u);
  EXPECT_EQ(Polyhedron::read_geo_header(f).checksum, 0u);
  EXPECT_TRUE(Polyhedron::verify_geo(f));
  ASSERT_TRUE(Polyhedron::to_geo(Polyhedron::C20(), f, true, empty));
  EXPECT_EQ(Polyhedron::read_geo_header(f).count, 1u);
  fclose(f);
}

TEST(GeoFormat, WriterRejectsBadInput) {
  const Polyhedron C20 = Polyhedron::C20();
  FILE* f = temp();
  EXPECT_THROW(Polyhedron::to_geo(C20, f, false, fixed_options(1)), std::invalid_argument);
  EXPECT_THROW(Polyhedron::to_geo(C20, f, false, fixed_options(31)), std::invalid_argument);
  geo_options offset_float;
  offset_float.offset = true;
  EXPECT_THROW(Polyhedron::to_geo(C20, f, false, offset_float), std::invalid_argument);
  geo_options scale_float;
  scale_float.scale = 1;
  EXPECT_THROW(Polyhedron::to_geo(C20, f, false, scale_float), std::invalid_argument);
  geo_options nothing;
  nothing.type = geo_type::NONE;
  nothing.graph = false;
  EXPECT_THROW(Polyhedron::to_geo(C20, f, false, nothing), std::invalid_argument);
  geo_options half_declared;
  half_declared.deg_min = 3;
  EXPECT_THROW(Polyhedron::to_geo(C20, f, false, half_declared), std::invalid_argument);

  geo_options triangles;
  triangles.triangulation = true;
  expect_code(Code::InvalidTopology, [&] { Polyhedron::to_geo(C20, f, false, triangles); });

  Polyhedron broken = C20;
  broken.points[7][1] = NAN;
  expect_code(Code::ValueOutOfRange, [&] { Polyhedron::to_geo(broken, f, false, {}); });
  expect_code(Code::ValueOutOfRange, [&] { Polyhedron::to_geo(broken, f, false, fixed_options(12)); });
  fclose(f);
}

// ===========================================================================
// Crash consistency and integrity
// ===========================================================================

TEST(GeoFormat, InterruptedAppendIsIgnoredThenOverwritten) {
  const geo_options opt = declared({}, 3, 0);
  FILE* f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(scaled(Polyhedron::C20(), 1), f, true, opt));
  ASSERT_TRUE(Polyhedron::to_geo(scaled(Polyhedron::C20(), 2), f, true, opt));
  const uint64_t R = Polyhedron::read_geo_header(f).record_size();

  // A record half written when the process died: count and checksum were never updated.
  fseeko(f, 0, SEEK_END);
  const Bytes garbage(R / 2, 0xab);
  ASSERT_EQ(fwrite(garbage.data(), 1, garbage.size(), f), garbage.size());
  fflush(f);
  EXPECT_EQ(Polyhedron::read_geo_header(f).count, 2u);
  EXPECT_TRUE(Polyhedron::verify_geo(f));
  const Polyhedron second = Polyhedron::from_geo(f, 1);
  expect_points(second.points, scaled(Polyhedron::C20(), 2).points);

  const Polyhedron third = scaled(Polyhedron::C20(), 3);
  ASSERT_TRUE(Polyhedron::to_geo(third, f, true, opt));
  EXPECT_EQ(bytes_of(f).size(), 32 + 3 * R);
  EXPECT_TRUE(Polyhedron::verify_geo(f));
  const Polyhedron read_back = Polyhedron::from_geo(f, 2);
  expect_points(read_back.points, third.points);
  fclose(f);
}

TEST(GeoFormat, ChecksumDetectsCorruption) {
  const geo_options opt = declared({}, 3, 0);
  FILE* f = temp();
  for (int i = 0; i < 3; i++) ASSERT_TRUE(Polyhedron::to_geo(scaled(Polyhedron::C20(), 1 + i), f, true, opt));
  const uint64_t R = Polyhedron::read_geo_header(f).record_size();
  const Bytes good = bytes_of(f);
  fclose(f);

  Bytes bad = good;   // one flipped mantissa bit in record 1
  bad[32 + R + 3] ^= 0x10;
  f = file_with(bad);
  EXPECT_FALSE(Polyhedron::verify_geo(f));
  EXPECT_NO_THROW(Polyhedron::from_geo(f, 1)) << "reading a record does not check the whole-file checksum";
  fclose(f);

  bad = good;         // records 0 and 1 swapped
  std::swap_ranges(bad.begin() + 32, bad.begin() + 32 + R, bad.begin() + 32 + R);
  f = file_with(bad);
  EXPECT_FALSE(Polyhedron::verify_geo(f));
  fclose(f);

  bad = good;         // count 2: the third record no longer counts
  bad[16] = 2;
  f = file_with(bad);
  EXPECT_FALSE(Polyhedron::verify_geo(f));
  fclose(f);

  bad = good;         // a changed fixed header field changes every record's seed
  bad[5] |= 0x04;     // the triangulation flag: same record size, and verify does not decode
  f = file_with(bad);
  EXPECT_FALSE(Polyhedron::verify_geo(f));
  fclose(f);

  f = file_with(good);
  EXPECT_TRUE(Polyhedron::verify_geo(f));
  fclose(f);
}

// ===========================================================================
// Reader rejections
// ===========================================================================

TEST(GeoFormat, HeaderRejections) {
  const auto rejects = [](Code code, Bytes b) {
    FILE* f = file_with(b);
    expect_code(code, [&] { Polyhedron::from_geo(f); });
    fclose(f);
  };
  const auto with = [](size_t i, uint8_t value) {
    Bytes b = worked_example;
    b[i] = value;
    return b;
  };
  rejects(Code::UnsupportedFormat, with(0, 0x30));   // version 1
  rejects(Code::UnsupportedFormat, with(0, 0x18));   // reserved spec bit
  rejects(Code::UnsupportedFormat, with(0, 0x01));   // no coordinates, no graph
  rejects(Code::UnsupportedFormat, with(0, 0x12));   // f32 with a width
  rejects(Code::UnsupportedFormat, with(4, 31));     // width
  rejects(Code::UnsupportedFormat, with(4, 1));
  rejects(Code::UnsupportedFormat, with(5, 0x0d));   // reserved flag bit
  rejects(Code::UnsupportedFormat, with(5, 0x04));   // file scale 0
  rejects(Code::UnsupportedFormat, with(7, 28));     // deg_bits
  rejects(Code::UnsupportedFormat, with(6, 1));      // degrees exactly 1: no surface
  rejects(Code::UnsupportedFormat, with(1, 2));      // N = 2 with a graph
  rejects(Code::UnsupportedFormat, with(9, 1));      // a file scale beside the per-record flag
  rejects(Code::UnsupportedFormat, with(12, 1));     // reserved bytes
  rejects(Code::MalformedFile, with(16, 2));         // two records claimed, one present
  Bytes truncated = worked_example;
  truncated.pop_back();
  rejects(Code::MalformedFile, truncated);
  rejects(Code::MalformedFile, Bytes(20, 0));
  rejects(Code::UnsupportedFormat, Bytes(32, 0));    // a zero header (an interrupted write)

  // The width byte is reported as a number: the library's generic to_string template
  // would print a uint8_t as a raw character, which is not even valid UTF-8 for 136.
  FILE* f = file_with(with(4, 136));
  try {
    Polyhedron::from_geo(f);
    ADD_FAILURE() << "width 136 accepted";
  } catch (const mesh_io_error& e) {
    EXPECT_NE(std::string(e.what()).find("width 136 for"), std::string::npos) << e.what();
  }
  fclose(f);
  try {
    geo::resolve(fixed_options(31), {});
    ADD_FAILURE() << "width 31 accepted";
  } catch (const std::invalid_argument& e) {
    EXPECT_NE(std::string(e.what()).find("width 31 outside"), std::string::npos) << e.what();
  }
}

TEST(GeoFormat, RecordRejections) {
  const auto code_of = [](const Bytes& b) {
    FILE* f = file_with(b);
    try {
      Polyhedron::from_geo(f);
    } catch (const mesh_io_error& e) {
      fclose(f);
      return int(e.code);
    }
    fclose(f);
    return -1;
  };
  const auto with = [](std::initializer_list<std::pair<size_t, uint8_t>> changes) {
    Bytes b = worked_example;
    for (auto [i, v] : changes) b[i] = v;
    return b;
  };
  EXPECT_EQ(code_of(with({{51, 0x01}})), int(Code::MalformedFile));            // padding byte
  EXPECT_EQ(code_of(with({{35, 0x7f}, {34, 0xff}})), int(Code::MalformedFile));// scale NaN / inf
  EXPECT_EQ(code_of(with({{32, 0}, {33, 0}, {34, 0}, {35, 0}})), int(Code::MalformedFile)); // scale 0
  EXPECT_EQ(code_of(with({{48, 0x6f}})), int(Code::InvalidTopology));          // twin 15 >= A = 12
  EXPECT_EQ(code_of(with({{48, 0x60}})), int(Code::InvalidTopology));          // twin 0 for half-edge 0
  EXPECT_EQ(code_of(with({{48, 0x64}})), int(Code::InvalidTopology));          // a torus (GEO-FORMAT.md sec. 11)
  EXPECT_EQ(code_of(with({{48, 0x31}, {49, 0x76}, {50, 0xb9}})), int(Code::InvalidTopology)); // a self-loop face
  EXPECT_EQ(code_of(worked_example), -1);
}

// Every twin block the tetrahedron's degrees admit: the reader accepts exactly the
// spheres with triangular faces, and whatever it accepts re-encodes to the same bytes.
TEST(GeoFormat, EveryTwinBlockOfTheTetrahedron) {
  geo_header H = geo::parse_header(std::span<const uint8_t, 32>(worked_example.data(), 32));
  const Bytes record(worked_example.begin() + 32, worked_example.end());

  int accepted = 0, rejected = 0;
  std::vector<int> entries;
  std::vector<char> matched(12, 0);
  const std::function<void()> enumerate = [&] {
    int h = 0;
    while (h < 12 && matched[h]) h++;
    if (h == 12) {
      Bytes r = record;
      r[16] = r[17] = r[18] = 0;
      for (int e = 0; e < 6; e++) r[16 + e / 2] |= uint8_t(entries[e] << (4 * (e % 2)));
      try {
        const geo_record decoded = geo::decode(H, r);
        EXPECT_EQ(geo::encode(H, decoded), r);
        accepted++;
      } catch (const mesh_io_error& e) {
        EXPECT_EQ(int(e.code), int(Code::InvalidTopology)) << e.what();
        rejected++;
      }
      return;
    }
    matched[h] = 1;
    for (int t = h + 1; t < 12; t++) {
      if (matched[t]) continue;
      matched[t] = 1;
      entries.push_back(t);
      enumerate();
      entries.pop_back();
      matched[t] = 0;
    }
    matched[h] = 0;
  };
  enumerate();
  EXPECT_EQ(accepted + rejected, 10395);   // 11!!
  EXPECT_EQ(accepted, 162);                // counted independently (face walk + Euler, in Python)
}

// ===========================================================================
// File names and delta-complexes
// ===========================================================================

TEST(GeoFormat, FileNamesChooseTheCoordinates) {
  namespace fs = std::filesystem;
  char tmpl[] = "/tmp/geo-io-test-XXXXXX";
  ASSERT_NE(mkdtemp(tmpl), nullptr);
  const fs::path dir(tmpl);
  const Polyhedron P = Polyhedron::C20();

  const auto header_of = [&](const std::string& name) {
    const std::string path = (dir / name).string();
    EXPECT_TRUE(Polyhedron::to_file(P, path));
    FILE* f = fopen(path.c_str(), "rb");
    const geo_header H = Polyhedron::read_geo_header(f);
    fclose(f);
    EXPECT_EQ(rows_of(Polyhedron::from_file(path)), rows_of(P));
    return H;
  };
  EXPECT_EQ(header_of("a.geo").opt.type, geo_type::F64);
  EXPECT_EQ(header_of("a.f32.geo").opt.type, geo_type::F32);
  EXPECT_EQ(header_of("x.y.geo").opt.type, geo_type::F64);
  const geo_header q = header_of("a.q12.geo");
  EXPECT_EQ(q.opt.type, geo_type::FIXED);
  EXPECT_EQ(q.opt.width, 12);
  EXPECT_TRUE(q.record_scale());
  EXPECT_FALSE(q.opt.offset);
  EXPECT_TRUE(q.opt.graph);
  const Polyhedron lossless = Polyhedron::from_file((dir / "a.geo").string());
  expect_points(lossless.points, P.points);
  EXPECT_THROW(Polyhedron::to_file(P, (dir / "a.q31.geo").string()), std::invalid_argument);
  expect_code(Code::NullFile, [&] { Polyhedron::from_file((dir / "missing.geo").string()); });
  fs::remove_all(dir);
}

namespace {

// The file's edge k (GEO-FORMAT.md sec. 7.5) as a half-edge of D, re-derived here from
// the rows D's writer uses: counter-clockwise from v_out.
std::vector<int> file_edges(const DelaunayTriangulation& D) {
  std::vector<int> order, slot(D.nh, -1);
  for (int v = 0; v < D.nv; v++) {
    int h = D.v_out[v];
    do { slot[h] = order.size(); order.push_back(h); h = D.ccw(h); } while (h != D.v_out[v]);
  }
  std::vector<int> edges;
  std::vector<char> matched(order.size(), 0);
  for (size_t i = 0; i < order.size(); i++) {
    if (matched[i]) continue;
    edges.push_back(order[i]);
    matched[i] = matched[slot[D.twin(order[i])]] = 1;
  }
  return edges;
}

}  // namespace

TEST(GeoFormat, DeltaComplexThroughDelaunayTriangulation) {
  // C60 isomer #1264: a 12-cone intrinsic Delaunay delta-complex with a multi-edge
  // (see DCEL.IdtRoundTrip in delaunay-test.cc).
  BuckyGen::buckygen_queue Q = BuckyGen::start(60, false, false);
  Triangulation T;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) { if (idx == 1264) break; idx++; }
  BuckyGen::stop(Q);
  const DelaunayTriangulation D = DelaunayTriangulation::compute(T);
  ASSERT_FALSE(D.is_simplicial());

  const std::vector<int> edges = file_edges(D);
  const auto length = [&](const DelaunayTriangulation&, int h) { return D.he_length[edges[h / 2]]; };
  const auto orig_degree = [&](int v) { return D.v_orig_degree[v]; };

  std::vector<coord3d> x;
  for (int v = 0; v < D.nv; v++) x.push_back({double(v), double(v * v) / 7, 1.0 / (v + 1)});
  for (geo_type type : {geo_type::NONE, geo_type::F64}) {
    geo_options opt;
    opt.type = type;
    const std::span<const coord3d> positions = type == geo_type::NONE ? std::span<const coord3d>{}
                                                                      : std::span<const coord3d>(x);
    FILE* f = temp();
    ASSERT_TRUE(DelaunayTriangulation::to_geo(D, positions, f, false, opt));
    EXPECT_TRUE(Polyhedron::read_geo_header(f).opt.triangulation);
    expect_code(type == geo_type::NONE ? Code::UnsupportedFormat : Code::NonSimplicial,
                [&] { Polyhedron::from_geo(f); });

    std::vector<coord3d> y;
    const DelaunayTriangulation E = DelaunayTriangulation::from_geo(f, 0, length, orig_degree, &y);
    EXPECT_TRUE(E.check_consistency());
    EXPECT_FALSE(E.is_simplicial());
    expect_points(y, positions);
    ASSERT_EQ(E.nv, D.nv);
    ASSERT_EQ(E.nh, int(2 * edges.size()));

    // An isomorphism of the two DCELs that fixes every vertex.
    std::vector<int> m(E.nh);
    for (size_t k = 0; k < edges.size(); k++) { m[2 * k] = edges[k]; m[2 * k + 1] = D.twin(edges[k]); }
    for (int h = 0; h < E.nh; h++) {
      EXPECT_EQ(E.he_origin[h], D.he_origin[m[h]]);
      EXPECT_EQ(m[E.he_next[h]], D.he_next[m[h]]);
      EXPECT_EQ(E.he_length[h], D.he_length[m[h]]);
    }
    for (int v = 0; v < D.nv; v++) EXPECT_NEAR(E.vertex_angle_sum(v), D.vertex_angle_sum(v), 1e-12);

    // Re-encoding reproduces the file.
    FILE* g = temp();
    ASSERT_TRUE(DelaunayTriangulation::to_geo(E, positions, g, false, opt));
    EXPECT_EQ(bytes_of(g), bytes_of(f));
    fclose(g);

    const auto too_long = [&](const DelaunayTriangulation& G, int h) { return h == 0 ? 1e9 : length(G, h); };
    EXPECT_THROW(DelaunayTriangulation::from_geo(f, 0, too_long, orig_degree), std::invalid_argument);
    const auto zero = [](const DelaunayTriangulation&, int) { return 0.0; };
    EXPECT_THROW(DelaunayTriangulation::from_geo(f, 0, zero, orig_degree), std::invalid_argument);
    fclose(f);
  }

  // A record that is not flagged as a triangulation is not read as one.
  FILE* f = temp();
  ASSERT_TRUE(Polyhedron::to_geo(tetrahedron(), f));
  expect_code(Code::NotATriangulation, [&] {
    DelaunayTriangulation::from_geo(f, 0, [](const DelaunayTriangulation&, int) { return 1.0; },
                                    [](int) { return 3; });
  });
  fclose(f);
}
