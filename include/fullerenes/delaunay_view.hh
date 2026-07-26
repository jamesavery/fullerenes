#pragma once

// DelaunayView -- trivially-copyable span view of the delta-complex DCEL.
//
// The owning type is DelaunayTriangulation (delaunay.hh), which inherits this
// view and points the spans into its own vectors (the graphview.hh migration
// pattern).  Everything here is device-includable: no <vector>, no
// <functional>, no exceptions, no I/O -- a GPU translation unit includes this
// header alone.
//
// Twin convention:      twin(h) = h^1  (half-edges 2k, 2k+1 are twins).
// Face orientation:     he_next traverses each face CCW.
// Vertex circulation:   cw(h) rotates CW around origin(h).
//
// CONST DOES NOT PROTECT THE ARRAYS: the members are span<T>, so element
// access through a const view (or a const owner) yields mutable references.
// The vectors' old const-propagation is gone; a read-only method mutating
// the complex is no longer a compile error.  Const-ness of a view is
// documentation, stated here because the compiler no longer states it.
//
// Scope (stage 2 of the span-DCEL promotion, PROMOTION-DESIGN.md): the SoA
// arrays, navigation, the intrinsic-geometry methods, the Diamond
// predicates, the capacity formulas, check_consistency, and the canonical
// field order.  Owner-level today and scheduled to move here in stage 3
// with the bounded workspace: the mutation machinery (allocation over free
// lists, flips, sweeps, vertex removal) AND the remaining allocation-free
// readers that ride with it (is_delaunay_edge, is_delaunay,
// count_non_delaunay, min_live_degree, is_well_formed over caller bits) --
// their staying put is sequencing, not design.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>

// The DCEL's element type is int, bound as the lib's 32-bit index width
// (PROMOTION-DESIGN.md Q4): machine-check the premise where the arrays live.
static_assert(std::is_same_v<int, std::int32_t>,
    "DelaunayView assumes int is the lib's 32-bit index type");

// ============================================================================
// Intrinsic geometry primitives.
// ============================================================================
namespace delaunay_detail {

inline constexpr double two_pi = 2 * std::numbers::pi_v<double>;

// Tolerance bands of the floating-point predicates.  These absorb FP noise
// in the cotangent/Heron evaluations, nothing else: the mathematics is
// is_delaunay <=> cot-sum >= 0 and strict convexity <=> both endpoint forms
// > 0.  Never widen a band to make a case pass.
inline constexpr double delaunay_band   = -1e-10;  // accept near-tight as tight
inline constexpr double convexity_band  =  1e-12;  // strictness margin
inline constexpr double cot_degenerate  =  1e15;   // sentinel for a degenerate triangle

// Heron product: H(a,b,c) = (a+b+c)(-a+b+c)(a-b+c)(a+b-c) = 16*Area^2.
// Returns 0 if the triangle inequality is violated.
inline double heron_product(double a, double b, double c) {
  double s1 = -a + b + c;
  double s2 =  a - b + c;
  double s3 =  a + b - c;
  if (s1 < 0 || s2 < 0 || s3 < 0) return 0;
  return (a + b + c) * s1 * s2 * s3;
}

// The Heron product in SQUARED-length coordinates, exact over the integers:
// for squared sides x = a^2, y = b^2, z = c^2,
//   H(x,y,z) = 2(xy + yz + zx) - (x^2 + y^2 + z^2) = 16*Area^2.
// The integer form the exact cocircularity predicate rationalizes over.
inline long long heron_product_sq(long long x, long long y, long long z) {
  return 2*(x*y + y*z + x*z) - (x*x + y*y + z*z);
}

// Cotangent of the angle opposite side `opp` in a triangle with sides
// (opp, b, c): cot(alpha) = (b^2 + c^2 - opp^2) / sqrt(H).
// Returns +/-cot_degenerate on a degenerate triangle (H <= 0) -- a value no
// sane tolerance ever brackets, so degenerate diamonds classify as
// non-tight rather than tripping the tight test.
inline double cot_opposite(double opp, double b, double c) {
  double H = heron_product(opp, b, c);
  double num = b*b + c*c - opp*opp;
  if (H <= 0) return (num >= 0) ? cot_degenerate : -cot_degenerate;
  return num / std::sqrt(H);
}

// Angle of the corner adjacent to sides `adj1`, `adj2` in a triangle whose
// opposite side is `opp`.  Law of cosines, clamped for floating-point safety
// at triangle-inequality boundaries.
inline double triangle_angle(double adj1, double adj2, double opp) {
  double c = (adj1*adj1 + adj2*adj2 - opp*opp) / (2 * adj1 * adj2);
  return std::acos(std::clamp(c, -1.0, 1.0));
}

}  // namespace delaunay_detail

// ============================================================================
// Diamond: the local geometry around an edge in a metrized triangulation.
//
//      B
//     / \          upper triangle: sides (e, a, b)
//    a   b         a = side adjacent to u, b = side adjacent to v
//   /     \  .
//  u---e---v       e = diagonal being tested/flipped
//   \     /
//    c   d         lower triangle: sides (e, c, d)
//     \ /          c = side adjacent to u, d = side adjacent to v
//      D
//
// All geometric predicates (Delaunay, convexity, flipped length) depend only
// on these five edge lengths.  No vertex IDs or topology needed.
// ============================================================================
struct Diamond {
  double e, a, b, c, d;

  // cot(angle_B) + cot(angle_D) >= 0, within delaunay_band.
  bool is_delaunay() const {
    using namespace delaunay_detail;
    return cot_opposite(e, a, b) + cot_opposite(e, c, d) >= delaunay_band;
  }

  // Angle sum < pi at both u and v, strictly (convexity_band margin).
  // The two endpoints are the SAME half-predicate applied to the two ends of
  // the diagonal (the diamond read from v swaps a<->b, c<->d): sin(angle) is
  // proportional to sqrt(Ha)*Q + P*sqrt(Hd), with (P, Q) the endpoint's two
  // law-of-cosines numerators.  Both tests share one evaluation of the two
  // Heron roots, so the arithmetic per endpoint is identical.
  bool is_convex() const {
    using namespace delaunay_detail;
    double e2 = e*e;
    double Ha = heron_product(e, a, b), Hd = heron_product(e, c, d);
    double sHa = (Ha > 0) ? std::sqrt(Ha) : 0;
    double sHd = (Hd > 0) ? std::sqrt(Hd) : 0;
    auto convex_at = [&](double P, double Q) {
      return sHa * Q + P * sHd > convexity_band;
    };
    return convex_at(e2 + a*a - b*b, e2 + c*c - d*d)      // at u
        && convex_at(e2 + b*b - a*a, e2 + d*d - c*c);     // at v
  }

  // Length of BD, the other diagonal.
  double flipped_length() const {
    using delaunay_detail::heron_product;
    // f^2 = a^2 + c^2 - (PQ - sqrt(Ha*Hd)) / (2e^2)
    double e2 = e*e, a2 = a*a, b2 = b*b, c2 = c*c, d2 = d*d;
    double P = e2 + a2 - b2;
    double Q = e2 + c2 - d2;
    double Ha = heron_product(e, a, b), Hd = heron_product(e, c, d);
    double sqrtHH = (Ha > 0 && Hd > 0) ? std::sqrt(Ha * Hd) : 0;
    double f2 = a2 + c2 - (P * Q - sqrtHH) / (2 * e2);
    return (f2 > 0) ? std::sqrt(f2) : 0;
  }

  // Cocircular ("tight") test: cot(angle_B) + cot(angle_D) == 0 exactly,
  // i.e. the four points u, v, B, D are concyclic on the surface.  In this
  // case both triangulations of the diamond are equally-valid Delaunay
  // refinements of the same cell.  Uses exact integer arithmetic on
  // length-squared (valid when all five lengths square to non-negative
  // integers, e.g. equilateral triangulations and their flips).
  // See CANONICAL-TESSELATION.md for the derivation.
  bool is_cocircular() const {
    // Tight Delaunay: cot(angle_B) + cot(angle_D) == 0 exactly.  Equivalent
    // to s1 * area_2 + s2 * area_1 = 0 where s1 = a^2+b^2-e^2, s2 = c^2+d^2-e^2.
    // Squaring after sign-check: tight iff sign(s1) != sign(s2) and
    // s1^2 * H2 == s2^2 * H1 (with H = heron_product_sq = 16*area^2).  Or:
    // both s1, s2 == 0.  Done in integer length-squared arithmetic so the
    // predicate is exact for equilateral triangulations and any sequence of
    // flips.
    using delaunay_detail::heron_product_sq;
    long long Le = (long long)std::llround(e * e);
    long long La = (long long)std::llround(a * a);
    long long Lb = (long long)std::llround(b * b);
    long long Lc = (long long)std::llround(c * c);
    long long Ld = (long long)std::llround(d * d);
    long long s1 = La + Lb - Le;
    long long s2 = Lc + Ld - Le;
    if (s1 == 0 && s2 == 0) return true;
    if (s1 == 0 || s2 == 0) return false;
    if ((s1 > 0) == (s2 > 0)) return false;          // same sign: not tight
    return s1 * s1 * heron_product_sq(Le, Lc, Ld)
        == s2 * s2 * heron_product_sq(Le, La, Lb);
  }

  // Floating-point cocircular test for general (non-equilateral) metrics,
  // where length-squared is not integer so the exact predicate above does not
  // apply: tight iff |cot(angle_B) + cot(angle_D)| < tol.  Scale-invariant
  // (cotangents are dimensionless), so tol is a pure angle threshold.
  // cot_opposite returns +/-1e15 on a degenerate triangle, which never lands
  // within a sane tol, so degenerate diamonds are correctly reported non-tight.
  bool is_cocircular(double tol) const {
    using delaunay_detail::cot_opposite;
    double cotB = cot_opposite(e, a, b);
    double cotD = cot_opposite(e, c, d);
    return std::abs(cotB + cotD) < tol;
  }
};

// ============================================================================
// dcel_capacities(nv0): the workspace sizes a bounded DCEL over a closed
// genus-0 triangulation with nv0 vertices is built to (nh = 2E = 6*nv0 - 12,
// nf = 2*nv0 - 4; the free lists hold at most one entry per edge / face).
// Reduction never exceeds the build sizes (flips rewire in place; removal is
// net-decreasing), so these are also the lifetime capacities of the whole
// compute() pipeline.  General delta-complexes (arbitrary cone count / genus,
// or post-build vertex insertion) do NOT satisfy these identities and must
// size explicitly.
// ============================================================================
struct DcelCapacities {
  long nh_cap;          // half-edge slots   = 6*nv0 - 12
  long nf_cap;          // face slots        = 2*nv0 - 4
  long free_edges_cap;  // recycled edge ids = nh_cap / 2  (dead until stage 3)
  long free_faces_cap;  // recycled face ids = nf_cap      (dead until stage 3)
};

// @pre nv0 >= 3 (the smallest closed triangulation); smaller inputs -- e.g.
// an unfilled batch slot -- yield all-zero capacities rather than the
// negative values the Euler identity would produce (which a size_t cast
// would turn into near-SIZE_MAX allocations).
constexpr DcelCapacities dcel_capacities(int nv0) {
  if (nv0 < 3) return { 0, 0, 0, 0 };
  const long nh = 6L * nv0 - 12;
  const long nf = 2L * nv0 - 4;
  return { nh, nf, nh / 2, nf };
}

// Pin the Euler identity on the icosahedron (nv0 = 12: E = 30, F = 20) and
// the degenerate guard.
static_assert(dcel_capacities(12).nh_cap == 60 && dcel_capacities(12).nf_cap == 20
              && dcel_capacities(12).free_edges_cap == 30
              && dcel_capacities(12).free_faces_cap == 20,
              "dcel_capacities must satisfy the Euler identity");
static_assert(dcel_capacities(0).nh_cap == 0 && dcel_capacities(2).nh_cap == 0,
              "dcel_capacities must guard nv0 < 3");

// ============================================================================
// DelaunayView: span SoA view of the half-edge DCEL.
// ============================================================================
struct DelaunayView {
  // --- Counts (dead slots included in nh/nf) ---
  int nv = 0;   // live vertices
  int nh = 0;   // allocated half-edges (including dead slots)
  int nf = 0;   // allocated faces (including dead slots)

  // --- Half-edge topology (indexed 0..nh-1) ---
  std::span<int>    he_next;    // next half-edge CCW in same face
  std::span<int>    he_origin;  // origin vertex (-1 = dead)
  std::span<int>    he_face;    // face to the left

  // --- Half-edge metric (indexed 0..nh-1) ---
  std::span<double> he_length;  // edge length (same for h and twin(h))
  std::span<double> he_angle;   // angle at origin(h) in face(h)

  // --- Per-vertex (indexed 0..nv-1) ---
  std::span<int>    v_out;          // one outgoing half-edge (-1 = dead vertex)
  std::span<double> v_cone_angle;   // cone angle = total vertex angle
  std::span<int>    v_orig_degree;  // degree in the original triangulation

  // --- Per-face (indexed 0..nf-1) ---
  std::span<int>    f_he;       // one boundary half-edge (-1 = dead face)

  // -------------------------------------------------------------------------
  // Navigation.
  // -------------------------------------------------------------------------
  int  twin(int h)  const { return h ^ 1; }
  int  edge(int h)  const { return h >> 1; }
  int  prev(int h)  const { return he_next[he_next[h]]; }  // only for triangulations
  int  dest(int h)  const { return he_origin[twin(h)]; }
  bool alive(int h) const { return he_origin[h] >= 0; }

  // Bigon edge: both half-edges of h bound the same face.  Arises in
  // delta-complexes around low-degree cone vertices (an i-j edge of an "iji"
  // face); dihedral quantities across such an edge are undefined.
  bool is_bigon(int h) const { return he_face[h] == he_face[twin(h)]; }

  // CW rotation around origin(h): next outgoing half-edge clockwise.
  int cw(int h) const { return he_next[twin(h)]; }

  // CCW rotation around origin(h): next outgoing half-edge counterclockwise.
  int ccw(int h) const { return twin(prev(h)); }

  // The three half-edges of (triangular) face f, in he_next order starting
  // from its representative.  Pre: f is live (f_he[f] >= 0).
  std::array<int,3> face_halfedges(int f) const {
    const int h = f_he[f];
    return {h, he_next[h], prev(h)};
  }
  // The three corner vertices of face f, CCW (origins of face_halfedges(f)).
  std::array<int,3> face_vertices(int f) const {
    const auto h = face_halfedges(f);
    return {he_origin[h[0]], he_origin[h[1]], he_origin[h[2]]};
  }

  // Count outgoing half-edges from v; 0 for a dead vertex.
  int vertex_degree(int v) const {
    int deg = 0;
    for ([[maybe_unused]] int h : incident(v)) deg++;   // empty when v_out[v] < 0
    return deg;
  }

  // Range over the outgoing half-edges of v (the cw ring from v_out[v]);
  // empty if v has no incident live edge.  The canonical vertex traversal.
  // SELF-CONTAINED VALUE: the range and its iterator carry the one span the
  // walk needs (cw(h) = he_next[h^1]), so a range taken from a temporary
  // view -- e.g. a by-value batch entry -- never dangles.
  struct IncidentHalfEdges {
    std::span<const int> he_next; int h0;
    struct iterator {
      std::span<const int> he_next; int start, h; bool done;
      int       operator*()  const { return h; }
      iterator& operator++()       { h = he_next[h ^ 1]; done = (h == start); return *this; }
      bool      operator!=(const iterator&) const { return !done; }
    };
    iterator begin() const { return {he_next, h0, h0, h0 < 0}; }
    iterator end()   const { return {}; }
  };
  IncidentHalfEdges incident(int v) const { return {he_next, v_out[v]}; }

  // Range over one (even) half-edge per live edge.  The canonical edge
  // traversal, skipping dead slots.  Self-contained value, like incident():
  // alive(h) = he_origin[h] >= 0 is the only fact the walk reads.
  struct LiveEdges {
    std::span<const int> he_origin; int nh;
    struct iterator {
      std::span<const int> he_origin; int h, nh;
      void advance() { while (h < nh && he_origin[h] < 0) h += 2; }
      int       operator*()  const { return h; }
      iterator& operator++()       { h += 2; advance(); return *this; }
      bool      operator!=(const iterator& o) const { return h != o.h; }
    };
    iterator begin() const { iterator it{he_origin, 0, nh}; it.advance(); return it; }
    iterator end()   const { return {he_origin, nh, nh}; }
  };
  LiveEdges edges() const { return {he_origin, nh}; }

  // -------------------------------------------------------------------------
  // Intrinsic geometry.
  // -------------------------------------------------------------------------

  // The five-length local geometry around edge h.
  Diamond diamond(int h) const {
    // h: u->v.  Face left of h has third vertex B = dest(next(h)).
    // Twin face has third vertex D = dest(next(twin(h))).
    int t = twin(h);

    double e_len = he_length[h];
    // a = edge from u to B, b = edge from v to B (in face of h)
    int h_vB = he_next[h];          // v->B
    int h_Bu = he_next[h_vB];       // B->u
    double a = he_length[h_Bu];     // side adjacent to u (B-u)
    double b = he_length[h_vB];     // side adjacent to v (v-B)

    // c = edge from u to D, d = edge from v to D (in face of twin)
    int h_uD = he_next[t];          // u->D
    int h_Dv = he_next[h_uD];       // D->v
    double c = he_length[h_uD];     // side adjacent to u (u-D)
    double d = he_length[h_Dv];     // side adjacent to v (D-v)

    return {e_len, a, b, c, d};
  }

  // Recompute the three corner angles of face f from its edge lengths.
  void recompute_face_angles(int f) {
    if (f_he[f] < 0) return;
    // h_i: u_i -> u_{i+1} with length L_i.  Angle at origin(h_i) is the
    // corner between sides L_i (outgoing) and L_{i-1} (incoming), opposite
    // to L_{i+1}.
    const auto h = face_halfedges(f);
    double L[3] = { he_length[h[0]], he_length[h[1]], he_length[h[2]] };
    for (int i = 0; i < 3; i++)
      he_angle[h[i]] = delaunay_detail::triangle_angle(L[i], L[(i + 2) % 3], L[(i + 1) % 3]);
  }

  // Recompute every corner angle (he_angle) from the current edge lengths,
  // then refresh the per-vertex cone-angle cache (v_cone_angle) that
  // curvature() / is_flat() / remove_flat_vertices() read.  Both are derived
  // from he_length, so this is THE entry point that re-establishes angle
  // coherence after any change to the metric.  Delaunay flips do NOT need
  // it: the cone angle is flip-invariant (the quad's interior angle at each
  // diamond vertex is independent of which diagonal is present -- paper
  // lem:coneflip), so flip_edge keeps v_cone_angle correct on its own.
  void recompute_all_angles() {
    recompute_all_face_angles();
    recompute_cone_angles();
  }

  // Recompute he_angle for every face, WITHOUT refreshing the v_cone_angle
  // cache.  For a hot caller (a line search) that reads only he_angle per
  // trial; pair with recompute_cone_angles once the kept metric needs the
  // cone cache.
  void recompute_all_face_angles() {
    for (int f = 0; f < nf; f++)
      recompute_face_angles(f);
  }

  // Refresh the cone-angle cache at one vertex.  @pre he_angle current.
  void recompute_cone_angle(int v) {
    v_cone_angle[v] = vertex_angle_sum(v);
  }

  // The equilateral metric: every edge length 1, every corner angle pi/3,
  // cone angle Theta_v = deg(v) * pi/3 (so kappa_v = (6 - deg v) * pi/3).
  // Writes exactly the nh live-slot prefix and the nv vertices.
  void set_equilateral_metric() {
    constexpr double pi = std::numbers::pi_v<double>;
    for (int h = 0; h < nh; h++) he_length[h] = 1.0;
    for (int h = 0; h < nh; h++) he_angle[h]  = pi / 3.0;
    for (int v = 0; v < nv; v++)
      v_cone_angle[v] = v_orig_degree[v] * pi / 3.0;   // (deg*pi)/3: the
                                  // reference association, kept to the ULP
  }

  // Refresh the cone-angle cache v_cone_angle[v] = vertex_angle_sum(v) for
  // every live vertex.  @pre he_angle current.  O(nh).
  void recompute_cone_angles() {
    for (int v = 0; v < nv; v++)
      if (v_out[v] >= 0)
        recompute_cone_angle(v);
  }

  // Total angle at vertex v = sum of incident corner angles (the paper's
  // Theta_v).  @pre angles current (recompute_all_angles after any length
  // change).
  double vertex_angle_sum(int v) const {
    // he_angle[h] is the angle at origin(h) in face(h): one corner per
    // outgoing half-edge.  incident(v) is empty when v_out[v] < 0.
    double sum = 0.0;
    for (int h : incident(v)) sum += he_angle[h];
    return sum;
  }

  // Discrete Gaussian curvature (angle defect) at v: kappa_v = 2*pi - Theta_v.
  // @pre v_cone_angle current (recompute_all_angles after any length change).
  double curvature(int v) const {
    return delaunay_detail::two_pi - v_cone_angle[v];
  }

  // Total curvature Sigma_v kappa_v over the live vertices.  Gauss-Bonnet:
  // equal to 4*pi on any closed genus-0 complex -- the cheapest global check
  // of the cone-angle cache, exposed for tests and validators.
  double total_curvature() const {
    double sum = 0.0;
    for (int v = 0; v < nv; v++)
      if (v_out[v] >= 0) sum += curvature(v);
    return sum;
  }

  // A vertex is flat (zero cone curvature, hence removable) iff
  // |kappa_v| < flat_tol.  flat_tol must separate the metric's curvature
  // noise floor from the smallest real cone curvature (pi/3 ~ 1.047): 1e-6
  // suits exact (equilateral) metrics, where any tol in (0, pi/3) agrees
  // exactly with the degree-6 test; a numerically solved metric (e.g. a
  // CEPS conformal solve, whose layout-cut seams leave ~5e-4 spurious
  // curvature at genuinely flat vertices) needs ~1e-2.
  // @pre v_cone_angle current.
  bool is_flat(int v, double flat_tol = 1e-6) const {
    return std::abs(curvature(v)) < flat_tol;
  }

  // -------------------------------------------------------------------------
  // Structural + metric validation: the DCEL class invariant, executable.
  // True iff all nine facts hold on every live element:
  //   @inv 1  twin closure: twin(h) < nh for live h
  //   @inv 2  triangular faces: he_next^3(h) == h
  //   @inv 3  origin chaining: dest(h) == origin(next(h))
  //   @inv 4  face agreement: face(next(h)) == face(h)
  //   @inv 5  v_out witnesses: v_out[v] live and originating at v
  //   @inv 6  f_he witnesses: f_he[f] live and lying in face f
  //   @inv 7  positive lengths on live half-edges
  //   @inv 8  twin length equality: length(h) == length(twin h)
  //   @inv 9  triangle inequality on every live face (1e-9 relative slack)
  // Allocation-free; O(nh).
  // -------------------------------------------------------------------------
  bool check_consistency() const {
    for (int h = 0; h < nh; h++)                                    // @inv 1
      if (alive(h) && twin(h) >= nh) return false;

    for (int h = 0; h < nh; h++)                                    // @inv 2
      if (alive(h) && he_next[he_next[he_next[h]]] != h) return false;

    for (int h = 0; h < nh; h++)                                    // @inv 3
      if (alive(h) && dest(h) != he_origin[he_next[h]]) return false;

    for (int h = 0; h < nh; h++)                                    // @inv 4
      if (alive(h) && he_face[he_next[h]] != he_face[h]) return false;

    for (int v = 0; v < nv; v++)                                    // @inv 5
      if (v_out[v] >= 0 && (!alive(v_out[v]) || he_origin[v_out[v]] != v))
        return false;

    for (int f = 0; f < nf; f++)                                    // @inv 6
      if (f_he[f] >= 0 && (!alive(f_he[f]) || he_face[f_he[f]] != f))
        return false;

    for (int h = 0; h < nh; h++)                                    // @inv 7
      if (alive(h) && he_length[h] <= 0) return false;

    for (int h = 0; h < nh; h += 2)                                 // @inv 8
      if (alive(h) && he_length[h] != he_length[twin(h)]) return false;

    for (int h = 0; h < nh; h++) {                                  // @inv 9
      if (!alive(h)) continue;
      int h1 = he_next[h], h2 = he_next[h1];
      double L0 = he_length[h], L1 = he_length[h1], L2 = he_length[h2];
      double eps = 1e-9 * (L0 + L1 + L2);
      if (L0 > L1 + L2 + eps) return false;
      if (L1 > L0 + L2 + eps) return false;
      if (L2 > L0 + L1 + eps) return false;
    }

    return true;
  }

  // -------------------------------------------------------------------------
  // Canonical field order.
  //
  // to_tuple() lists the 9 SoA arrays in THE canonical order (the scalar
  // counts nv/nh/nf are not fields):
  //   { he_next, he_origin, he_face, he_length, he_angle,
  //     v_out, v_cone_angle, v_orig_degree, f_he }
  // This order is the contract the owner's repoint() fold and any field-wise
  // state comparator build on.  Allocation SIZING is dcel_capacities' job --
  // this view deliberately does not model the parent's batch::batchable_view
  // concept (its per-vertex multipliers cannot express the DCEL's affine
  // sizes); a batch adapter derives its tables from dcel_capacities.
  // -------------------------------------------------------------------------
  static constexpr std::size_t n_fields = 9;

  auto to_tuple() {
    return std::forward_as_tuple(he_next, he_origin, he_face,
                                 he_length, he_angle,
                                 v_out, v_cone_angle, v_orig_degree, f_he);
  }
  auto to_tuple() const {
    return std::forward_as_tuple(he_next, he_origin, he_face,
                                 he_length, he_angle,
                                 v_out, v_cone_angle, v_orig_degree, f_he);
  }
};

static_assert(std::is_trivially_copyable_v<DelaunayView>,
              "DelaunayView must be trivially copyable");
static_assert(std::tuple_size_v<decltype(std::declval<DelaunayView>().to_tuple())>
                  == DelaunayView::n_fields,
              "to_tuple arity must equal n_fields");
