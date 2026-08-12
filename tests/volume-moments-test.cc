// Moments of the solid enclosed by a Polyhedron, and the two MASS MODELS of the
// inertia tensor.
//
// Since 2026-08-07 `PolyhedronView<double>::inertia_matrix(MassModel)` answers
// about the ENCLOSED SOLID by default (exact signed-tetrahedron sums over the
// faces) and about the VERTICES only on explicit opt-in, and BOTH models are
// CENTRAL -- taken about the mass distribution's own centre, so the tensor is
// invariant under translating the polyhedron.  These tests pin:
//
//   1. the moment kernel against an ANALYTIC SOLID: a box, whose volume,
//      centroid and covariance are closed form.  This is the only oracle that
//      pins the kernel's CONSTANTS -- every other assertion here is scale-free
//      or self-referential, and `inertia_from_second_moment` is linear, so a
//      wrong scalar factor would propagate through all of them unchanged.
//      Splitting one face of the same box also pins that the answer does NOT
//      move when a patch is triangulated more finely, and drives the
//      parallel-axis correction at nonzero magnitude (the vertex mean, which is
//      the kernel's apex, then differs from the volume centroid).
//   2. volume_moments().mass == volume_tetra() on a SHIFTED cage -- shifted
//      because on a centred one the two reduce to the same floating-point sum
//      and the comparison is vacuous.  Shifted, it also measures what the
//      kernel's re-centred apex buys over volume_tetra()'s origin apex.
//   3. s^3 / s^5 scaling of volume and second moment.
//   4. the vertex branch == the centred vertex sum written out by hand
//   5. translation invariance of BOTH branches (the point of the change), with
//      the pre-2026-08-07 origin-based sum shown to move under the same shift,
//      and the centroid shown to move by exactly the shift vector
//   6. the pinned molecular call sites are UNCHANGED by the centring: each does
//      move_to_origin() first, and there inertia_matrix(ATOMS) reproduces the old
//      origin-based sum to roundoff
//   7. the NAMED OUTCOMES of a surface that cannot be integrated: an empty
//      triangle list, non-finite coordinates, and -- the one that is silently
//      wrong without a check -- a MIS-ORIENTED face list, whose covariance comes
//      out finite with a negative eigenvalue
//   8. the two models' TENSORS agree to 5.9% (||dI||/||I||, unit-trace) but their
//      FRAMES do not: an axis both models RESOLVE still turns by 3.40 deg (3.70
//      at a marginal gap), so the two frames are NOT interchangeable on a
//      low-symmetry cage.  That measurement is why the existing program call
//      sites are pinned to MassModel::Atoms -- the molecular convention -- rather
//      than migrated to the new default.  Symmetric cages are the exception: the
//      squashed icosahedral C80 agrees to 0.0006 deg (gated at 0.01).  Raw
//      per-axis angles reach 6.27 deg and are printed but deliberately NOT gated
//      -- inside a near-degenerate eigenplane no individual direction is
//      determined, so the angle there measures nothing physical.

#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "fullerenes/polyhedron.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/buckygen-wrapper.hh"

using namespace std;

namespace {

// ---------------------------------------------------------------------------
// Comparison vocabulary.  One relative-difference operator everywhere, so no
// assertion in this file compares a max-element difference against a Frobenius
// norm (they differ by up to a factor 3, which silently eats headroom).
// ---------------------------------------------------------------------------

// Relative difference of two tensors, Frobenius throughout: ||A-B|| / ||B||.
double rel_diff(const matrix3d& A, const matrix3d& B) { return (A - B).norm() / B.norm(); }

// The SHAPE of a mass distribution, scale removed: the tensor normalised to unit
// trace.  Needed only where two DIFFERENT mass models are compared -- they carry
// different total masses, so their tensors differ by a scale factor that says
// nothing about the body.
matrix3d unit_trace(const matrix3d& M) { return M * (1.0 / M.trace()); }

// Angle in degrees between two axes, insensitive to eigenvector sign.
double axis_angle_deg(const coord3d& u, const coord3d& v) {
  const double c = min(1.0, fabs(u.dot(v)) / (u.norm() * v.norm()));
  return acos(c) * 180.0 / M_PI;
}

coord3d row(const matrix3d& M, int i) { return coord3d(M(i, 0), M(i, 1), M(i, 2)); }

// Argmax over a stream of named measurements: the largest value seen, and who
// produced it.  One word for what was seven hand-written `if (x > w) { w = x; }`.
struct Worst {
  string name;
  double value = -1;
  void observe(const string& who, double x) { if (x > value) { value = x; name = who; } }
};

// ---------------------------------------------------------------------------
// Oracles (deliberately written out in full, not calling the library).
// ---------------------------------------------------------------------------

// The pre-2026-08-07 vertex sum: uniform mass at the vertices, about the ORIGIN.
matrix3d vertex_inertia_about_origin(const vector<coord3d>& pts) {
  matrix3d I;
  for (size_t k = 0; k < pts.size(); k++) {
    const coord3d& x(pts[k]);
    const long double xx(x.dot(x));
    for (int i = 0; i < 3; i++) {
      I(i, i) += xx;
      for (int j = 0; j < 3; j++) I(i, j) -= x[i] * x[j];
    }
  }
  return I;
}

// The same sum, centred on the vertex mean -- what inertia_matrix(ATOMS) is now.
matrix3d vertex_inertia_centred(const vector<coord3d>& pts) {
  coord3d c(0, 0, 0);
  for (const coord3d& p : pts) c += p;
  c /= (double)pts.size();

  matrix3d I;
  for (size_t k = 0; k < pts.size(); k++) {
    const coord3d x(pts[k] - c);
    const long double xx(x.dot(x));
    for (int i = 0; i < 3; i++) {
      I(i, i) += xx;
      for (int j = 0; j < 3; j++) I(i, j) -= x[i] * x[j];
    }
  }
  return I;
}

// volume_tetra()'s arithmetic written out: signed tetrahedra from the ORIGIN,
// with no re-centred apex.  The comparison that shows what the kernel's apex
// choice is for.
double tetra_volume_about_origin(const vector<coord3d>& p, const vector<tri_t>& tris) {
  double V = 0;
  for (const tri_t& t : tris) V += p[t[0]].dot(p[t[1]].cross(p[t[2]]));
  return fabs(V / 6.0);
}

// The kernel's arithmetic WITHOUT its positive-semidefiniteness postcondition --
// i.e. what volume_moments() would return if the check were removed.  Used only
// to show that on a mis-oriented surface that value is finite, passes every
// other guard, and has a negative eigenvalue.
matrix3d unguarded_solid_covariance(const vector<coord3d>& p, const vector<tri_t>& tris,
                                    double* volume_out) {
  coord3d o(0, 0, 0);
  for (const coord3d& x : p) o += x;
  o /= (double)p.size();

  double det_sum = 0;
  coord3d m1(0, 0, 0);
  matrix3d m2;
  for (const tri_t& t : tris) {
    const coord3d a = p[t[0]] - o, b = p[t[1]] - o, c = p[t[2]] - o;
    const double det = a.dot(b.cross(c));
    const coord3d s = a + b + c;
    det_sum += det;
    m1 += s * det;
    m2 += (s.outer(s) + a.outer(a) + b.outer(b) + c.outer(c)) * (det / 120.0);
  }
  const double sgn = (det_sum < 0) ? -1.0 : 1.0;
  const double vol = sgn * det_sum / 6.0;
  *volume_out = vol;
  const coord3d cr = m1 * (sgn / (24.0 * vol));
  return m2 * sgn - cr.outer(cr) * vol;
}

// ---------------------------------------------------------------------------
// An analytic solid: the axis-aligned box [0,l] + origin.
// ---------------------------------------------------------------------------
// 8 corners indexed by their bits, corner i = origin + (bit0*lx, bit1*ly, bit2*lz),
// and the 12 outward-oriented triangles of its surface.
struct Box {
  coord3d l, origin;
  vector<coord3d> points;
  vector<tri_t> tris;

  double volume() const { return l[0] * l[1] * l[2]; }
  coord3d centroid() const { return origin + l * 0.5; }
  // ∫ (x-c)(x-c)^T dV over the box = diag(lx^2, ly^2, lz^2) * V/12.
  matrix3d covariance() const {
    matrix3d C;
    for (int i = 0; i < 3; i++) C(i, i) = l[i] * l[i] * volume() / 12.0;
    return C;
  }
};

Box make_box(const coord3d& l, const coord3d& origin = coord3d(0, 0, 0)) {
  Box B{l, origin, {}, {}};
  for (int i = 0; i < 8; i++)
    B.points.push_back(origin + coord3d((i & 1) * l[0], ((i >> 1) & 1) * l[1],
                                        ((i >> 2) & 1) * l[2]));
  B.tris = {{0,2,3},{0,3,1},   // -z
            {4,5,7},{4,7,6},   // +z
            {0,1,5},{0,5,4},   // -y
            {2,6,7},{2,7,3},   // +y
            {0,4,6},{0,6,2},   // -x
            {1,3,7},{1,7,5}};  // +x
  return B;
}

// Split the +z face into a fan of 4 around its centre, leaving the SOLID
// untouched and the surface still closed and outward-oriented.  Two things at
// once: the answer must not move when a patch is triangulated more finely, and
// the added point drags the vertex MEAN -- which is the kernel's tetrahedron
// apex -- off the volume centroid, so the parallel-axis correction
// `- cr (x) cr * vol` is exercised at nonzero magnitude instead of at ~0.
void split_top_face(Box& B) {
  const int c = (int)B.points.size();
  B.points.push_back(B.origin + coord3d(B.l[0] / 2, B.l[1] / 2, B.l[2]));
  B.tris.erase(B.tris.begin() + 2, B.tris.begin() + 4);   // the two +z triangles
  for (const tri_t& t : vector<tri_t>{{4,5,c},{5,7,c},{7,6,c},{6,4,c}})
    B.tris.push_back(t);
}

// A triangle surface is closed and consistently oriented iff every directed edge
// appears exactly once (so its reverse carries the neighbouring triangle).  The
// kernel's stated precondition, checked rather than assumed -- the kernel takes
// a bare triangle LIST, which carries no rotation system for
// GraphView::is_consistently_oriented() to look at.
bool is_closed_oriented(const vector<tri_t>& tris) {
  map<arc_t, int> seen;
  for (const tri_t& t : tris)
    for (int i = 0; i < 3; i++) seen[{t[i], t[(i + 1) % 3]}]++;
  for (const auto& [a, n] : seen)
    if (n != 1 || seen.count({a.second, a.first}) == 0) return false;
  return true;
}

// ---------------------------------------------------------------------------
// Eigenvalue conditioning: which axes are determined, and which only span a
// subspace.  Both are read off the ADJACENT relative gaps of the ascending
// spectrum -- for an ascending triple every non-adjacent gap is a sum of
// adjacent ones, so the two adjacent gaps carry all the information.
// ---------------------------------------------------------------------------
using Gaps = array<double, 2>;

Gaps relative_gaps(const coord3d& lam) {
  const double scale = max(fabs(lam[0]), max(fabs(lam[1]), fabs(lam[2])));
  return {(lam[1] - lam[0]) / scale, (lam[2] - lam[1]) / scale};
}

// The more pessimistic of two spectra, gap by gap: an axis counts as determined
// only if BOTH mass models resolve it.  Deriving the partition from one model
// and applying it to both would treat the two symmetric subjects asymmetrically,
// and in one direction would pass on broken code.
Gaps common_gaps(const Gaps& a, const Gaps& b) { return {min(a[0], b[0]), min(a[1], b[1])}; }

// Per-axis separation: the smaller of the gaps to its neighbours.
coord3d separations(const Gaps& g) { return coord3d(g[0], min(g[0], g[1]), g[1]); }

// Group axes whose gap is below `tol`: inside a group no individual direction is
// determined, so only the SUBSPACE is comparable.  Exactly complementary to
// separations(): axis i is a singleton here iff separations(g)[i] >= tol, which
// is why the two gates below share ONE tolerance -- see kResolvedGap.
// @pre ascending: g[0] >= 0 && g[1] >= 0
vector<vector<int>> clusters(const Gaps& g, double tol) {
  vector<vector<int>> cl{{0}};
  for (int i = 1; i < 3; i++) {
    if (g[i - 1] < tol) cl.back().push_back(i);
    else                cl.push_back({i});
  }
  return cl;
}

// Orthogonal projector onto the span of the eigenvector rows in `idx`.
matrix3d projector(const matrix3d& axes, const vector<int>& idx) {
  matrix3d P;
  for (int i : idx) { const coord3d u = row(axes, i); P += u.outer(u); }
  return P;
}

// The largest principal angle between two subspaces of equal dimension:
// sin(theta_max) = ||(Id - P_U) P_V||_2, obtained from the eigenvalues of M^T M.
// Unlike an eigenvector angle this is well defined when the eigenvalues inside
// the subspace are degenerate -- rotating within the subspace does not change it.
// @pre projectors: (Pu*Pu - Pu).norm() < 1e-9 && (Pv*Pv - Pv).norm() < 1e-9
// @pre same_rank:  Pu.trace() == Pv.trace()
double max_principal_angle_deg(const matrix3d& Pu, const matrix3d& Pv) {
  const matrix3d M = (matrix3d::unit_matrix() - Pu) * Pv;
  const coord3d s2 = (M.transpose() * M).symmetric_part().eigensystem().first;
  const double s = sqrt(max(0.0, max(s2[0], max(s2[1], s2[2]))));
  return asin(min(1.0, s)) * 180.0 / M_PI;
}

// The worst subspace angle between two frames over a shared cluster partition.
double subspace_angle_deg(const matrix3d& A, const matrix3d& B, const vector<vector<int>>& cl) {
  double worst = 0;
  for (const vector<int>& c : cl)
    worst = max(worst, max_principal_angle_deg(projector(A, c), projector(B, c)));
  return worst;
}

// ---------------------------------------------------------------------------
// Fixtures.
// ---------------------------------------------------------------------------

// A fixed rotation (axis (1,2,3)/|.|, angle 0.7 rad), so the test cages are NOT
// already sitting in their own principal frame -- an identity answer would then
// be the guard fallback, not the real frame.
matrix3d fixed_rotation() { return matrix3d::rotation(coord3d(1, 2, 3), 0.7); }

// C60-Ih (leapfrog of the icosahedron) with force-field geometry -- no database
// and no buckygen needed.
Polyhedron make_C60() {
  Triangulation ico(vector<int>(12, 5));
  Triangulation gc = ico.GCtransform(1, 1);
  return Polyhedron::fullerene_polyhedron(FullereneGraph(gc.dual_graph()));
}

// C80 = GC(2,0) of the icosahedron, same route.
Polyhedron make_C80() {
  Triangulation ico(vector<int>(12, 5));
  Triangulation gc = ico.GCtransform(2, 0);
  return Polyhedron::fullerene_polyhedron(FullereneGraph(gc.dual_graph()));
}

// The first `count` general C_N isomers off buckygen -- ordinary, low-symmetry
// cages (the icosahedral ones have isotropic tensors, whose eigenvectors carry
// no information at all, so an axis comparison on them is vacuous).
vector<Polyhedron> make_isomers(int N, int count) {
  vector<Polyhedron> out;
  auto Q = BuckyGen::start(N, false, false);
  Graph G;
  for (int i = 0; i < count && BuckyGen::next_fullerene(Q, G); i++)
    out.push_back(Polyhedron::fullerene_polyhedron(
        FullereneGraph(Triangulation(G).dual_graph())));
  BuckyGen::stop(Q);
  EXPECT_EQ((int)out.size(), count) << "buckygen returned too few C" << N << " isomers";
  return out;
}

// fullerene_polyhedron() leaves the cage sitting in its own principal frame, so
// rotate it out: an identity answer from principal_axes() then means the
// degeneracy guard fired, not that the frame happens to be the coordinate axes.
void rotate_out_of_frame(Polyhedron& P) { P.transform(fixed_rotation()); }

// An icosahedral cage is isotropic -- both models give a multiple of the
// identity, whose eigenvectors mean nothing.  Squash it (distinct factors on all
// three axes) so the eigenvalues separate and the frame becomes comparable.
void squash(Polyhedron& P) { P.scale(coord3d(1.0, 1.15, 1.35)); }

}  // namespace

// ---------------------------------------------------------------------------
// 1. The kernel against a closed-form solid.
// ---------------------------------------------------------------------------

// The ONLY absolute oracle in this file.  Every other assertion on the second
// moment is scale-free or compares the kernel with itself, and
// inertia_from_second_moment is linear, so a wrong constant (det/60 instead of
// det/120, a dropped s (x) s term, swapped diagonal and off-diagonal entries of
// the reference tetrahedron's S) would survive all of them unchanged.  A box
// pins the constant AND the structure: its covariance is diagonal with three
// DIFFERENT entries in a known ratio.
TEST(MassMoments, BoxMatchesClosedForm) {
  const Box B = make_box(coord3d(3, 5, 7), coord3d(2, -3, 5));
  ASSERT_TRUE(is_closed_oriented(B.tris)) << "box fixture is not a closed oriented surface";

  const MassMoments m = volume_moments(std::span<const coord3d>(B.points), B.tris);

  EXPECT_EQ(m.code, MassMoments::Code::Ok);
  EXPECT_NEAR(m.mass, B.volume(), 1e-13 * B.volume());
  EXPECT_LT((m.centroid - B.centroid()).norm(), 1e-13 * B.l.norm());
  EXPECT_LT(rel_diff(m.covariance, B.covariance()), 1e-13)
      << "got " << m.covariance << " want " << B.covariance();
  // A second moment is a sum of outer products x (x) x, whose (i,j) and (j,i)
  // entries are the SAME product -- so exact symmetry, not symmetry to roundoff.
  EXPECT_EQ((m.covariance - m.covariance.transpose()).norm(), 0.0);
  const matrix3d I = matrix3d::unit_matrix() * m.covariance.trace() - m.covariance;
  EXPECT_EQ((I - I.transpose()).norm(), 0.0);
}

// The same solid with one face triangulated four ways instead of two.  Nothing
// about the body changed, so nothing about the moments may change -- this is the
// "a densely triangulated patch does not move the answer" claim, measured.  It
// also moves the vertex mean (the kernel's apex) off the volume centroid, which
// is the only way the parallel-axis correction is exercised at all: on every
// symmetric fixture the two coincide and the correction is identically zero.
TEST(MassMoments, RefiningAFaceMovesNothing) {
  const Box plain = make_box(coord3d(3, 5, 7), coord3d(2, -3, 5));
  Box refined = make_box(coord3d(3, 5, 7), coord3d(2, -3, 5));
  split_top_face(refined);
  ASSERT_TRUE(is_closed_oriented(refined.tris)) << "refined box is not closed/oriented";

  // 8 corners (mean = the box centre) plus one point on the +z face pull the
  // mean up by lz/2 * 1/9, i.e. the apex sits lz/18 above the centroid.  Exact,
  // so the magnitude the parallel-axis term is exercised at is a stated number,
  // not "some nonzero amount".
  coord3d apex(0, 0, 0);
  for (const coord3d& p : refined.points) apex += p;
  apex /= (double)refined.points.size();
  const double offset = (apex - refined.centroid()).norm();
  printf("  refined box: apex-to-centroid |cr| = %.6g (expected lz/18 = %.6g)\n",
         offset, refined.l[2] / 18);
  ASSERT_NEAR(offset, refined.l[2] / 18, 1e-12 * refined.l[2])
      << "the parallel-axis term is not being exercised where this test thinks it is";

  const MassMoments m0 = volume_moments(std::span<const coord3d>(plain.points), plain.tris);
  const MassMoments m1 = volume_moments(std::span<const coord3d>(refined.points), refined.tris);

  EXPECT_EQ(m1.code, MassMoments::Code::Ok);
  EXPECT_NEAR(m1.mass, m0.mass, 1e-13 * m0.mass);
  EXPECT_LT((m1.centroid - m0.centroid).norm(), 1e-13 * plain.l.norm());
  EXPECT_LT(rel_diff(m1.covariance, m0.covariance), 1e-13);
  EXPECT_LT(rel_diff(m1.covariance, refined.covariance()), 1e-13);
}

// An inward-wound surface describes the same solid.  The kernel normalises the
// orientation sign so no caller has to know the winding; before this case that
// branch was never entered, since every fixture in the file is outward-oriented.
TEST(MassMoments, FlippedWindingIsTheSameSolid) {
  const Box B = make_box(coord3d(3, 5, 7), coord3d(2, -3, 5));
  Box F = B;
  for (tri_t& t : F.tris) t.flip();
  ASSERT_TRUE(is_closed_oriented(F.tris)) << "flipping every triangle broke consistency";

  const MassMoments m = volume_moments(std::span<const coord3d>(B.points), B.tris);
  const MassMoments f = volume_moments(std::span<const coord3d>(F.points), F.tris);

  EXPECT_EQ(f.code, MassMoments::Code::Ok);
  EXPECT_EQ(f.mass, m.mass);
  EXPECT_EQ((f.centroid - m.centroid).norm(), 0.0);
  EXPECT_EQ((f.covariance - m.covariance).norm(), 0.0);
}

// What the kernel's re-centred apex buys.  volume_tetra() sums tetrahedra from
// the ORIGIN, so on a body sitting far from it the terms are huge and the answer
// is the small difference between them; the kernel puts the apex at the vertex
// mean, where each term is the size of the body.  The two are the same integral
// in exact arithmetic -- this measures how far apart floating point drives them.
TEST(MassMoments, ApexChoiceSurvivesLargeOffsets) {
  const coord3d l(3, 5, 7);
  const double V = l[0] * l[1] * l[2];
  double worst_kernel = 0;
  for (double d : {0.0, 1e3, 1e6, 1e8}) {
    const Box B = make_box(l, coord3d(d, d, d));
    const MassMoments m = volume_moments(std::span<const coord3d>(B.points), B.tris);
    const double rel_kernel = fabs(m.mass - V) / V;
    const double rel_origin = fabs(tetra_volume_about_origin(B.points, B.tris) - V) / V;
    worst_kernel = max(worst_kernel, rel_kernel);
    printf("  box at offset %-8.0e : kernel rel err %.3g, origin-apex rel err %.3g\n",
           d, rel_kernel, rel_origin);
    EXPECT_LT(rel_kernel, 1e-9) << "kernel lost the volume at offset " << d;
  }
  printf("  worst kernel relative error over the offsets: %.3g\n", worst_kernel);
}

// ---------------------------------------------------------------------------
// 2-3. The method against volume_tetra(), and scaling.
// ---------------------------------------------------------------------------

// Both integrate centroid_surface(), so they are the same integral and must
// agree.  On a CENTRED cage the comparison is nearly vacuous -- the two reduce to
// the same floating-point sum and the difference is exactly 0 -- so the cage is
// shifted first, which is where the two arithmetics actually differ.
TEST(MassMoments, ShiftedCageVolumeMatchesTetra) {
  for (Polyhedron P : {Polyhedron::C20(), make_C60()}) {
    const double centred_diff =
        fabs(P.volume_moments().mass - P.volume_tetra()) / P.volume_tetra();

    P.move(coord3d(17.3, -5.1, 2.9));
    const MassMoments m = P.volume_moments();
    const double rel = fabs(m.mass - P.volume_tetra()) / P.volume_tetra();
    printf("  N=%d: rel diff centred %.3g, shifted %.3g\n", P.N, centred_diff, rel);

    EXPECT_EQ(m.code, MassMoments::Code::Ok);
    EXPECT_GT(m.mass, 0);
    EXPECT_LT(rel, 1e-9);
  }
}

// Scaling the polyhedron by s scales the enclosed volume by s^3 and the second
// moment (an integral of x x^T dV) by s^5.
TEST(MassMoments, ScalesAsCube) {
  Polyhedron P = make_C60();
  const MassMoments m0 = P.volume_moments();

  const double s = 2.5;
  P.scale(coord3d(s, s, s));
  const MassMoments m1 = P.volume_moments();

  EXPECT_NEAR(m1.mass, m0.mass * s * s * s, 1e-12 * m0.mass * s * s * s);
  EXPECT_LT(rel_diff(m1.covariance, m0.covariance * pow(s, 5)), 1e-12);
}

// ---------------------------------------------------------------------------
// 4-5. The two mass models.
// ---------------------------------------------------------------------------

// inertia_matrix(ATOMS) is the centred vertex sum, written out by hand here.
TEST(MassMoments, VertexBranchIsTheCentredVertexSum) {
  for (Polyhedron P : {Polyhedron::C20(), make_C60()}) {
    vector<coord3d> pts(P.points.begin(), P.points.end());
    const matrix3d oracle = vertex_inertia_centred(pts);
    const matrix3d got = P.inertia_matrix(MassModel::Atoms);
    // Same arithmetic in a different accumulation order: equal to roundoff.
    EXPECT_LT(rel_diff(got, oracle), 1e-12)
        << "N=" << P.N << " got=" << got << " oracle=" << oracle;
  }
}

// The point of the 2026-08-07 change: neither model depends on where the
// polyhedron sits.  The pre-change origin-based sum, on the same shift, does.
TEST(MassMoments, BothBranchesAreTranslationInvariant) {
  Polyhedron P = make_C60();
  const coord3d d(17.3, -5.1, 2.9);

  const matrix3d Iv0 = P.inertia_matrix(MassModel::Solid);
  const matrix3d Ix0 = P.inertia_matrix(MassModel::Atoms);
  const coord3d  c0  = P.volume_moments().centroid;
  vector<coord3d> pts0(P.points.begin(), P.points.end());
  const matrix3d Iold0 = vertex_inertia_about_origin(pts0);

  P.move(d);
  const matrix3d Iv1 = P.inertia_matrix(MassModel::Solid);
  const matrix3d Ix1 = P.inertia_matrix(MassModel::Atoms);
  const coord3d  c1  = P.volume_moments().centroid;
  vector<coord3d> pts1(P.points.begin(), P.points.end());
  const matrix3d Iold1 = vertex_inertia_about_origin(pts1);

  const double res_v = rel_diff(Iv1, Iv0), res_x = rel_diff(Ix1, Ix0);
  const double res_c = (c1 - c0 - d).norm() / d.norm();
  printf("  C60 shifted by |d| = %.4g: residuals  volume %.3g  vertex %.3g  centroid %.3g"
         "  (gates 1e-9, 1e-9, 1e-12)\n", d.norm(), res_v, res_x, res_c);

  EXPECT_LT(res_v, 1e-9) << "volume branch moved";
  EXPECT_LT(res_x, 1e-9) << "vertex branch moved";
  // Covariance is translation-INVARIANT; the centroid is translation-COVARIANT.
  // Only the latter proves the centroid is a point on the body rather than a
  // constant: returning {0,0,0} unconditionally would pass every other gate here.
  EXPECT_LT(res_c, 1e-12) << "centroid did not follow the shift";

  // The behaviour that changed: the old origin-based tensor is dominated by the
  // offset (parallel-axis term ~ N|d|^2) and is a different tensor entirely.
  EXPECT_GT(rel_diff(Iold1, Iold0), 1.0)
      << "the pre-change sum was expected to move under translation";
  printf("  old origin-based tensor changed by %.3g relative under the same shift\n",
         rel_diff(Iold1, Iold0));
}

// ---------------------------------------------------------------------------
// 6. The pinned molecular call sites.
// ---------------------------------------------------------------------------

// Pinning the molecular call sites to MassModel::Atoms restores the
// PRE-2026-08-07 answer only if the polyhedron was already centred when they
// called -- the vertex sum used to be taken about the ORIGIN and is now taken
// about the vertex mean.  All of them do `move_to_origin()` immediately before,
// on the same object with no intervening change to the points, so the centring
// subtraction is a no-op there.  "No-op" is not "exact", though: move_to_origin()
// zeroes the mean only to roundoff, and the branch was re-expressed as
// tr(M)Id - M over a second-moment sum rather than the old interleaved long
// double accumulation.  This measures what is left.
//
// NOT ALL of the pinned sites are verified by building them: five are live
// (Polyhedron::fullerene_polyhedron and programs/{fullerene-isomers-polyhedra,
// fullerene-polyhedron,polyhedron-optimize,halma-polyhedron}.cc), while
// programs/density.cc is commented out of programs/CMakeLists.txt and
// programs/{gaudi,gs-ex}.cc appear in no build file at all and no longer compile
// against the current API.  What this test pins is therefore the PATTERN
// (move_to_origin(); then the ATOMS model) that all eight share.
TEST(MassMoments, PinnedMolecularCallSitesUnchanged) {
  // A tensor perturbation of relative size eps rotates an eigenvector with
  // relative eigenvalue gap g by about eps/g (Davis-Kahan).  With eps bounded by
  // kMaxTensorRelDiff below, an axis whose gap is at least kMinGap turns by at
  // most 1e-14/1e-6 = 1e-8 rad = 5.7e-7 deg -- under kMaxAngleDeg.  That is the
  // relation between the three constants; without it they are three unrelated
  // magic numbers.
  const double kMaxTensorRelDiff = 1e-14;
  const double kMinGap           = 1e-6;
  const double kMaxAngleDeg      = 1e-6;

  vector<Polyhedron> cages = make_isomers(60, 18);      // #17 is the worst frame case
  cages.push_back(make_C60());

  Worst w_rel, w_sub, w_raw;
  for (size_t i = 0; i < cages.size(); i++) {
    Polyhedron& P = cages[i];
    const string name = "C60 #" + to_string(i);
    P.move_to_origin();                                 // what every pinned site does first

    const vector<coord3d> pts(P.points.begin(), P.points.end());
    const matrix3d I_old = vertex_inertia_about_origin(pts);   // the pre-change answer
    const matrix3d I_new = P.inertia_matrix(MassModel::Atoms);  // the pinned answer
    w_rel.observe(name, rel_diff(I_new, I_old));

    // The frame, compared the only way that means anything: a cage with an
    // exactly degenerate pair (many isomers have one by symmetry) hands the two
    // eigenvectors back in either order, which shows up as a spurious 90 deg
    // between individual vectors while the SUBSPACE is identical.
    const matrix3d A_old = I_old.eigensystem().second;
    const matrix3d A_new = P.principal_axes(MassModel::Atoms);
    for (int k = 0; k < 3; k++)
      w_raw.observe(name, axis_angle_deg(row(A_old, k), row(A_new, k)));
    w_sub.observe(name, subspace_angle_deg(A_old, A_new,
        clusters(relative_gaps(I_new.eigensystem().first), kMinGap)));
  }
  printf("  pinned pattern over %zu centred cages: worst ||dI||/||I|| = %.3g (%s), worst frame "
         "subspace angle = %.3g deg (%s)\n         worst raw per-vector %.3g deg (%s) -- "
         "degenerate pairs, not gated\n",
         cages.size(), w_rel.value, w_rel.name.c_str(), w_sub.value, w_sub.name.c_str(),
         w_raw.value, w_raw.name.c_str());

  // Roundoff, not a change of answer: the alignment the molecular programs
  // produce is the one they produced before.
  EXPECT_LT(w_rel.value, kMaxTensorRelDiff);
  EXPECT_LT(w_sub.value, kMaxAngleDeg);
}

// ---------------------------------------------------------------------------
// 7. Surfaces that cannot be integrated -- the named outcomes.
// ---------------------------------------------------------------------------

// The default path needs a face structure where the old vertex sum needed none:
// volume_moments() calls faces().  So the ways a surface can fail to exist are
// part of this method's contract, and each one has a NAME.
//
// Until 2026-08-09 one of them had no name and no test: on an EDGELESS view
// (N points, deg == 0) faces() returned an empty list without throwing, and
// centroid_triangulation({}) then took its "already a triangulation" branch
// (vacuously true for zero faces) into orient_triangulation(), which indexed
// tris[0] unconditionally -- SIGSEGV, through all four integrals alike.  Case (a)
// could therefore only be made against the KERNEL.  The invariant work (plan step
// 1 + 4) removed both halves of that: the edgeless view now FAILS
// is_consistently_oriented() instead of passing vacuously, and
// orient_triangulation is gone, so centroid_surface() answers with
// CentroidSurface::Code::NoFaces and case (a) exercises the method.
//
// This was a TRACKED outcome, not folklore: claude-projects/curvature-flow/
// refactor-debt.md, entry 2026-08-07-orientation-invariant-not-explicit, which
// named the cure (make orientation an invariant maintained by construction
// rather than repaired at every use site) and the verification.
TEST(MassMoments, SurfacesThatCannotBeIntegrated) {
  const matrix3d Id = matrix3d::unit_matrix();
  using SurfaceCode = Polyhedron::CentroidSurface::Code;

  // (a) NO FACES: an edgeless view -- N points, not one edge.  Both the kernel's
  // own guard (which still protects direct callers, who pass a triangle list with
  // no view behind it) and the METHOD, which now names the outcome instead of
  // crashing on it.
  {
    vector<coord3d> pts{{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1}};
    const MassMoments m = volume_moments(std::span<const coord3d>(pts), {});
    EXPECT_EQ(m.code, MassMoments::Code::Degenerate);
    EXPECT_EQ(m.mass, 0.0);
    EXPECT_EQ(m.covariance.norm(), 0.0);

    Polyhedron P;
    static_cast<Owned<PolyhedronView<double>>&>(P) = Owned<PolyhedronView<double>>(int(pts.size()));
    P.set_points(pts);
    ASSERT_EQ(P.count_edges(), 0u);

    const Polyhedron::CentroidSurface S = P.centroid_surface();
    EXPECT_EQ(S.code, SurfaceCode::NoFaces);
    EXPECT_TRUE(S.points.empty());
    EXPECT_TRUE(S.tris.empty());

    // The documented value of every scalar integral when there is no surface.
    EXPECT_EQ(P.volume_moments().code, MassMoments::Code::Degenerate);
    EXPECT_EQ(P.surface_area(), 0.0);
    EXPECT_EQ(P.volume(),       0.0);
    EXPECT_EQ(P.volume_tetra(), 0.0);
  }

  // (b) NaN COORDINATES on a valid topology (is_invalid() == true).
  {
    Polyhedron P = Polyhedron::C20();
    P.points[0] = coord3d(NAN, NAN, NAN);
    ASSERT_TRUE(P.is_invalid());
    const MassMoments m = P.volume_moments();
    const matrix3d I = P.inertia_matrix();
    const InertialFrame F = P.principal_frame();
    printf("  NaN points:  volume=%g ||I||=%g ||A-Id||=%g\n",
           m.mass, I.norm(), (F.axes - Id).norm());
    EXPECT_EQ(m.code, MassMoments::Code::Degenerate);   // !(NaN > 0) -> the zeroed guard value
    EXPECT_EQ(m.mass, 0.0);
    EXPECT_EQ(I.norm(), 0.0);
    EXPECT_EQ(F.code, InertialFrame::Code::DegenerateMass);
    EXPECT_LT((F.axes - Id).norm(), 1e-12);
    EXPECT_LT((P.principal_axes() - Id).norm(), 1e-12);   // the matrix-only overload agrees
  }

  // (c) MIS-ORIENTED FACES -- the one that is confidently WRONG without a check.
  // Half the surface wound inwards makes det_sum cancel while the first moment
  // does not, so the enclosed "volume" survives the `volume > 0` guard as a
  // near-zero number, the centroid ratio blows up, and the covariance comes out
  // FINITE with a negative eigenvalue.  A second moment is a Gram matrix and is
  // positive semidefinite by construction, so lambda_min < 0 is an exact
  // postcondition violation -- no tolerance is involved, and it is the only
  // signal this input produces.
  {
    Box B = make_box(coord3d(1, 1, 1), coord3d(-0.5, -0.5, -0.5));
    B.points[7] += coord3d(1e-3, 0, 0);          // break the exact cancellation
    for (int i : {2, 3, 6, 7, 10, 11}) B.tris[i].flip();   // +z, +y and +x faces
    ASSERT_FALSE(is_closed_oriented(B.tris)) << "the fixture is supposed to be mis-oriented";

    double unguarded_volume = 0;
    const matrix3d cov = unguarded_solid_covariance(B.points, B.tris, &unguarded_volume);
    const double lam_min = cov.eigenvalues()[2];
    printf("  mis-oriented faces: unguarded volume = %.3g (passes 'volume > 0'), "
           "||cov|| = %.3g (finite), lambda_min = %.3g\n",
           unguarded_volume, cov.norm(), lam_min);
    ASSERT_GT(unguarded_volume, 0.0) << "the volume guard was expected NOT to fire here";
    ASSERT_TRUE(std::isfinite(cov.norm())) << "the covariance was expected to be finite";
    ASSERT_LT(lam_min, 0.0) << "the covariance was expected to have a negative eigenvalue";

    const MassMoments m = volume_moments(std::span<const coord3d>(B.points), B.tris);
    EXPECT_EQ(m.code, MassMoments::Code::NotPositiveSemidefinite);
    EXPECT_EQ(m.mass, 0.0);
    EXPECT_EQ(m.covariance.norm(), 0.0);
  }

  // (d) NOT ORIENTABLE: a cage with faces, but not a sphere's.  Transposing two
  // entries of one degree-3 row of C20 merges three pentagons into one 15-arc
  // face -- 12 faces become 10, chi = 0 -- so it is a perfectly consistent
  // rotation system OF A TORUS (tests/oriented-surface-test.cc makes that case
  // in full).  faces() refuses it by throwing, which is right for faces() and
  // useless to a scalar integral, so centroid_surface() is where that becomes a
  // named outcome the three doubles can be documented against.
  {
    Polyhedron P = Polyhedron::C20();
    ASSERT_EQ(P.centroid_surface().code, SurfaceCode::Ok);

    std::swap(P[0][0], P[0][1]);
    ASSERT_EQ(P.oriented_surface().code, OrientedSurface::Code::GenusMismatch);
    ASSERT_EQ(P.oriented_surface().genus, 1);
    EXPECT_THROW((void)P.faces(), std::logic_error);

    const Polyhedron::CentroidSurface S = P.centroid_surface();
    EXPECT_EQ(S.code, SurfaceCode::NotOrientable);
    EXPECT_TRUE(S.points.empty());
    EXPECT_TRUE(S.tris.empty());

    EXPECT_EQ(P.volume_moments().code, MassMoments::Code::Degenerate);
    EXPECT_EQ(P.surface_area(), 0.0);
    EXPECT_EQ(P.volume(),       0.0);
    EXPECT_EQ(P.volume_tetra(), 0.0);
  }
}

// The two triangulations preserve the faces' orientation, so a MIXED
// triangle-and-polygon mesh comes out consistent.  It did not before 2026-08-09:
// the polygon branch wound its boundary arc backwards while a triangle face was
// copied through unchanged, so on the square pyramid below the quad's four
// triangles and the four lateral faces traversed their shared edges in the SAME
// direction -- a surface that is not oriented at all, previously reaching the
// integrals unchecked (centroid_triangulation's orientation-repair call was
// commented out, and the repair pass would have abort()ed the process anyway).
//
// The face list is the whole input: centroid_triangulation reads only N and
// `faces`, so the view underneath needs no edges.
TEST(CentroidTriangulation, MixedTriangleAndPolygonFacesStayConsistent) {
  // Square pyramid: quad base (0,3,2,1) and four lateral triangles, every
  // directed edge appearing exactly once -- the outward-oriented surface.
  const vector<face_t> faces{face_t(vector<node_t>{0,3,2,1}), face_t(vector<node_t>{0,1,4}),
                             face_t(vector<node_t>{1,2,4}),   face_t(vector<node_t>{2,3,4}),
                             face_t(vector<node_t>{3,0,4})};

  Polyhedron P;
  static_cast<Owned<PolyhedronView<double>>&>(P) = Owned<PolyhedronView<double>>(5);

  const vector<tri_t> centroid = P.centroid_triangulation(faces);
  EXPECT_EQ(centroid.size(), 8u);          // 4 around the base centroid + 4 laterals
  EXPECT_TRUE(is_closed_oriented(centroid));

  // The fan is the same claim about the other triangulation, and doubles as the
  // check that the fixture face list was oriented to begin with (a fan adds only
  // diagonals, so it is oriented iff the faces are).
  const vector<tri_t> fan = P.triangulation(faces);
  EXPECT_EQ(fan.size(), 6u);               // 2 for the base + 4 laterals
  EXPECT_TRUE(is_closed_oriented(fan));

  // A face too small to be a polygon is refused, not read past its end.
  const vector<face_t> degenerate{face_t(vector<node_t>{0,1})};
  EXPECT_THROW((void)P.centroid_triangulation(degenerate), std::logic_error);
  EXPECT_THROW((void)P.triangulation(degenerate),          std::logic_error);
}

// ---------------------------------------------------------------------------
// 8. How far the default moves the answer for the existing callers.
// ---------------------------------------------------------------------------

// Everything measurable about the two frames of one cage.  Pure: it reads the
// polyhedron, computes, returns -- no captured accumulators, and NO gtest
// assertion inside, because an ASSERT_* in a helper returns from the HELPER, not
// from the TEST, and a "verified" precondition that merely stops one lambda
// early is worse than none.  The test body does all the gating.
namespace {
struct FrameComparison {
  matrix3d Av, Ax;              // the two frames, rows = axes
  InertialFrame::Code code_v, code_x;
  coord3d  lam;                 // ascending eigenvalues of the unit-trace SOLID tensor
  coord3d  sep;                 // per-axis separation, the more pessimistic of the two models
  double   drel = 0;            // ||Ivn - Ixn|| / ||Ivn||, both unit-trace
  double   ang[3] = {0,0,0};    // per-axis angle between the frames, degrees
  double   sub = 0;             // worst subspace angle over the common clusters
  double   fallback_v = 0, fallback_x = 0;    // ||A - Id||: ~0 means the guard fired
  double   unitarity_v = 0, unitarity_x = 0;  // ||A A^T - Id||
  bool     ascending_v = false, ascending_x = false;
};

// `cluster_tol` is the gap at which an axis stops being individually determined;
// it is the SAME number the caller uses to decide which axes to gate individually
// (kResolvedGap), which is what makes the two gates a partition.
FrameComparison compare_frames(const Polyhedron& P, double cluster_tol) {
  FrameComparison c;
  const InertialFrame Fv = P.principal_frame(MassModel::Solid);
  const InertialFrame Fx = P.principal_frame(MassModel::Atoms);
  c.Av = Fv.axes;      c.Ax = Fx.axes;
  c.code_v = Fv.code;  c.code_x = Fx.code;

  const matrix3d Id = matrix3d::unit_matrix();
  c.fallback_v = (c.Av - Id).norm();                    c.fallback_x = (c.Ax - Id).norm();
  c.unitarity_v = (c.Av * c.Av.transpose() - Id).norm();
  c.unitarity_x = (c.Ax * c.Ax.transpose() - Id).norm();

  // Scale-free comparison of the two mass models: both tensors normalised to
  // unit trace, so only the SHAPE of the mass distribution is compared.
  const matrix3d Ivn = unit_trace(P.inertia_matrix(MassModel::Solid));
  const matrix3d Ixn = unit_trace(P.inertia_matrix(MassModel::Atoms));
  c.drel = rel_diff(Ixn, Ivn);   // relative to the SOLID tensor -- the default answer

  // matrix3d::eigensystem sorts by |lambda|, which coincides with ascending
  // lambda only for a positive-semidefinite tensor.  An inertia tensor is one,
  // so index 0 is the smallest moment = the LONGEST geometric axis -- but that
  // is a property of the input, so it is checked rather than assumed.
  const coord3d lam_v = Ivn.eigensystem().first, lam_x = Ixn.eigensystem().first;
  c.ascending_v = lam_v[0] <= lam_v[1] && lam_v[1] <= lam_v[2];
  c.ascending_x = lam_x[0] <= lam_x[1] && lam_x[1] <= lam_x[2];
  c.lam = lam_v;

  const Gaps g = common_gaps(relative_gaps(lam_v), relative_gaps(lam_x));
  c.sep = separations(g);
  // Only the MULTI-axis clusters, so the two gates partition the three axes
  // rather than overlap: a singleton cluster is exactly an axis the caller gates
  // individually, and its "subspace angle" is that same axis angle.  A cage whose
  // three eigenvalues all cluster contributes a 3-dimensional subspace, whose
  // principal angle to itself is 0 -- vacuous, and correctly so: on a near-
  // isotropic tensor no direction is determined and only `drel` says anything.
  vector<vector<int>> multi;
  for (const vector<int>& cl : clusters(g, cluster_tol))
    if (cl.size() > 1) multi.push_back(cl);
  c.sub = subspace_angle_deg(c.Av, c.Ax, multi);
  for (int i = 0; i < 3; i++) c.ang[i] = axis_angle_deg(row(c.Av, i), row(c.Ax, i));
  return c;
}
}  // namespace

// The quantity that is WELL CONDITIONED is the tensor difference: the two mass
// models answer with the same shape to ~6% (Frobenius, unit-trace normalised)
// over the cages below.  An eigenVECTOR angle is not well conditioned -- inside a
// near-degenerate eigenplane no individual direction is determined at all, and
// any rotation there is an equally valid frame -- so this test gates on
//   (i)  the tensor difference, and
//   (ii) angles that survive the conditioning: the per-axis angle for axes
//        separated from BOTH neighbours, and the largest principal angle between
//        eigen-SUBSPACES for the rest.
// The two (ii) gates share ONE constant, kResolvedGap, and that is not a
// coincidence to be maintained by hand: axis i is a singleton cluster exactly
// when its separation reaches the same tolerance, so with one number the two
// gates PARTITION the three axes -- every axis is gated by exactly one of them.
// Two literals that merely happen to be equal would silently stop partitioning
// the moment one of them was tuned.
//
// The threshold sits clear of the knee rather than 10% from it.  MEASURED
// 2026-08-08 over the 41 cages below (the run prints the live version of every
// number here), separation vs worst per-axis angle:
//     sep 0.016  ->  6.27 deg      (C70 #2)
//     sep 0.041  ->  6.17 deg      (C70 #7)
//     sep 0.058  ->  3.70 deg      (C60 #17)
//     sep 0.072  ->  2.52 deg      (C60 #17)
//     sep 0.084  ->  3.40 deg      (C70 #16)   <- worst ADMITTED at kResolvedGap
// The cliff is at separation ~0.04, not at 0.07: the previous gate at 0.05 sat
// 10% from it (C70 #7 had separation 0.045 under the old one-sided derivation and
// per-axis angles 6.17/6.15, ABOVE the 6.0 bound, and passed only because the
// separation test excluded it -- an 11% shift in that eigenvalue gap would have
// failed the test).  At 0.07 the worst admitted angle is 3.40 deg and the worst
// angle among axes JUST BELOW the threshold is 3.70 deg, so the 6.0 deg bound
// would still hold if the threshold drifted all the way down to ~0.05.  The
// earlier "~60% headroom" claim was true of the admitted population and false at
// its boundary; what is true here is stated per gate below.
TEST(MassMoments, VolumeAndVertexFramesAgreeOnCages) {
  const double kMaxTensorRelDiff = 0.10;    // measured worst 0.0593 (C70 #2)
  const double kMaxFrameAngleDeg = 6.0;     // measured worst 3.40 admitted, 3.70 just below
  // ONE gap tolerance, used twice -- see the note above.  Chosen clear of the
  // knee in the separation-vs-angle data, not adjacent to it.
  const double kResolvedGap = 0.07;
  const double kSymmetricAngleDeg = 0.01;   // squashed icosahedral cage, measured 0.0006

  struct Case { string name; Polyhedron P; double angle_bound; };
  vector<Case> cases;
  for (int N : {60, 70}) {
    vector<Polyhedron> cages_N = make_isomers(N, 20);
    for (size_t i = 0; i < cages_N.size(); i++)
      cases.push_back({"C" + to_string(N) + " #" + to_string(i), cages_N[i], kMaxFrameAngleDeg});
  }
  // An icosahedral cage, squashed until its eigenvalues separate.  Symmetry
  // forces the two models onto the same axes, which bounds where the
  // low-symmetry disagreement above comes from: not from the weighting as such,
  // but from the shape asymmetry the two models resolve differently.
  { Polyhedron P = make_C80(); squash(P); cases.push_back({"C80-GC(2,0) squashed", P, kSymmetricAngleDeg}); }

  // Fold 1: MEASURE.
  vector<FrameComparison> measured;
  for (Case& c : cases) {
    rotate_out_of_frame(c.P);   // so an identity frame can only be the guard fallback
    measured.push_back(compare_frames(c.P, kResolvedGap));
  }
  ASSERT_EQ(measured.size(), 41u) << "the sweep did not produce the expected cage set";

  // Fold 2: GATE.
  int wellsep_axes = 0;
  for (size_t k = 0; k < cases.size(); k++) {
    const Case& c = cases[k];
    const FrameComparison& f = measured[k];

    EXPECT_EQ(f.code_v, InertialFrame::Code::Ok) << c.name << " (solid frame)";
    EXPECT_EQ(f.code_x, InertialFrame::Code::Ok) << c.name << " (atoms frame)";
    EXPECT_GT(f.fallback_v, 1e-3) << c.name << ": solid frame is the guard fallback";
    EXPECT_GT(f.fallback_x, 1e-3) << c.name << ": atoms frame is the guard fallback";
    EXPECT_LT(f.unitarity_v, 1e-9) << c.name << ": solid frame not unitary";
    EXPECT_LT(f.unitarity_x, 1e-9) << c.name << ": atoms frame not unitary";
    EXPECT_TRUE(f.ascending_v) << c.name << ": solid eigenvalues not ascending";
    EXPECT_TRUE(f.ascending_x) << c.name << ": atoms eigenvalues not ascending";
    EXPECT_LT(f.drel, kMaxTensorRelDiff) << c.name;

    for (int i = 0; i < 3; i++)
      if (f.sep[i] >= kResolvedGap) {          // determined axes: a real disagreement
        wellsep_axes++;
        EXPECT_LT(f.ang[i], c.angle_bound)
            << c.name << " axis " << i << " (separation " << f.sep[i] << ")";
      }
    EXPECT_LT(f.sub, c.angle_bound) << c.name << " (subspace)";   // the rest
  }
  EXPECT_GE(wellsep_axes, 80) << "too few determined axes to make the gate meaningful";  // 94 measured

  // Fold 3: WORST, and the separation-vs-angle table the threshold is read off.
  Worst w_raw, w_sub, w_wellsep, w_long, w_drel;
  Worst w_below;   // worst angle among axes JUST BELOW the threshold -- the knee
  for (size_t k = 0; k < cases.size(); k++) {
    const string& name = cases[k].name;
    const FrameComparison& f = measured[k];
    w_drel.observe(name, f.drel);
    w_long.observe(name, f.ang[0]);
    w_sub.observe(name, f.sub);
    for (int i = 0; i < 3; i++) {
      w_raw.observe(name, f.ang[i]);
      if (f.sep[i] >= kResolvedGap) w_wellsep.observe(name, f.ang[i]);
      else                          w_below.observe(name, f.ang[i]);
    }
    if (f.ang[0] > 2.0 || f.ang[1] > 2.0 || f.ang[2] > 2.0)
      printf("    %-22s lambda %5.3f %5.3f %5.3f  sep %6.4f %6.4f %6.4f  "
             "angles %6.3f %6.3f %6.3f  subspace %6.3f  dI/I %6.4f\n",
             name.c_str(), f.lam[0], f.lam[1], f.lam[2], f.sep[0], f.sep[1], f.sep[2],
             f.ang[0], f.ang[1], f.ang[2], f.sub, f.drel);
  }
  const FrameComparison& sym = measured.back();   // the squashed icosahedral cage
  printf("    %-22s sep %6.4f %6.4f %6.4f  angles %7.5f %7.5f %7.5f  subspace %7.5f  "
         "dI/I %6.4f   <- symmetry forces agreement\n",
         cases.back().name.c_str(), sym.sep[0], sym.sep[1], sym.sep[2],
         sym.ang[0], sym.ang[1], sym.ang[2], sym.sub, sym.drel);
  printf("  %zu cages, %d axes resolved at gap >= %.2f | worst dI/I %.4f (%s)\n"
         "  RAW %.3f deg (%s, NOT meaningful) | SUBSPACE %.3f deg (%s) | RESOLVED %.3f deg (%s)"
         " | LONGEST %.3f deg (%s)\n"
         "  worst angle among UNRESOLVED axes (gap < %.2f): %.3f deg (%s) -- the knee the "
         "threshold must stay clear of\n",
         cases.size(), wellsep_axes, kResolvedGap, w_drel.value, w_drel.name.c_str(),
         w_raw.value, w_raw.name.c_str(), w_sub.value, w_sub.name.c_str(),
         w_wellsep.value, w_wellsep.name.c_str(), w_long.value, w_long.name.c_str(),
         kResolvedGap, w_below.value, w_below.name.c_str());
}
