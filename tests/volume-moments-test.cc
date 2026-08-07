// Moments of the solid enclosed by a Polyhedron, and the two weightings of the
// inertia tensor.
//
// Since 2026-08-07 `PolyhedronView<double>::inertia_matrix(bool)` answers about
// the ENCLOSED SOLID by default (exact signed-tetrahedron sums over the faces)
// and about the VERTICES only on explicit opt-in, and BOTH weightings are
// CENTRAL -- taken about the mass distribution's own centre, so the tensor is
// invariant under translating the polyhedron.  These tests pin:
//
//   * volume_moments().volume == volume_tetra()   (the same solid, two routes)
//   * the vertex branch == the centred vertex sum written out by hand
//   * translation invariance of BOTH branches (the point of the change), with
//     the pre-2026-08-07 origin-based sum shown to move under the same shift
//   * the pinned molecular call sites are UNCHANGED by the centring: each does
//     move_to_origin() first, and there inertia_matrix(true) reproduces the old
//     origin-based sum to 5.8e-16 relative (frame subspace angle 1.3e-14 deg)
//   * the two weightings' TENSORS agree to ~6% (||dI||/||I|| < 0.10, measured
//     0.0593 over 41 cages) -- the well-conditioned invariant, and the only
//     agreement between the weightings that survives scrutiny
//   * their FRAMES do not: a well-separated axis still turns by up to 3.70 deg,
//     so the two frames are NOT interchangeable on a low-symmetry cage.  That
//     measurement is why the existing program call sites are pinned to
//     inertia/principal_axes/align_with_axes(vertex_weighted = true) -- the
//     molecular convention -- rather than migrated to the new default.
//     Symmetric cages are the exception: the icosahedral C80 agrees to 0.001 deg.
//     Raw per-axis angles are printed but deliberately NOT gated -- inside a
//     near-degenerate eigenplane (C70: two short axes 1.6% apart) no individual
//     direction is determined, so the angle there measures nothing physical.
//   * scaling by s scales the volume by s^3 and the second moment by s^5

#include <gtest/gtest.h>
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

// The same sum, centred on the vertex mean -- what inertia_matrix(true) is now.
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

double max_abs_diff(const matrix3d& A, const matrix3d& B) {
  double m = 0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m = max(m, fabs(A(i, j) - B(i, j)));
  return m;
}

// Angle in degrees between two axes, insensitive to eigenvector sign.
double axis_angle_deg(const coord3d& u, const coord3d& v) {
  const double c = min(1.0, fabs(u.dot(v)) / (u.norm() * v.norm()));
  return acos(c) * 180.0 / M_PI;
}

coord3d row(const matrix3d& M, int i) { return coord3d(M(i, 0), M(i, 1), M(i, 2)); }

matrix3d identity3() { matrix3d Id; Id(0,0) = Id(1,1) = Id(2,2) = 1; return Id; }

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
double max_principal_angle_deg(const matrix3d& Pu, const matrix3d& Pv) {
  const matrix3d M = (identity3() - Pu) * Pv;
  const coord3d s2 = (M.transpose() * M).symmetric_part().eigensystem().first;
  const double s = sqrt(max(0.0, max(s2[0], max(s2[1], s2[2]))));
  return asin(min(1.0, s)) * 180.0 / M_PI;
}

// Group ASCENDING eigenvalues into clusters of near-degenerate ones: consecutive
// eigenvalues join the same cluster when their gap, relative to the largest
// eigenvalue, is below `tol`.  Inside a cluster no individual direction is
// determined, so only the SUBSPACE is comparable.
vector<vector<int>> degenerate_clusters(const coord3d& lam, double tol) {
  const double scale = max(fabs(lam[0]), max(fabs(lam[1]), fabs(lam[2])));
  vector<vector<int>> cl{{0}};
  for (int i = 1; i < 3; i++) {
    if ((lam[i] - lam[i-1]) / scale < tol) cl.back().push_back(i);
    else                                   cl.push_back({i});
  }
  return cl;
}

// Per-axis separation: the smaller of the two relative eigenvalue gaps to the
// other axes.  An axis is only DETERMINED when this is comfortably nonzero.
coord3d axis_separations(const coord3d& lam) {
  const double scale = max(fabs(lam[0]), max(fabs(lam[1]), fabs(lam[2])));
  coord3d s;
  for (int i = 0; i < 3; i++) {
    double g = 1e300;
    for (int j = 0; j < 3; j++) if (j != i) g = min(g, fabs(lam[i] - lam[j]));
    s[i] = g / scale;
  }
  return s;
}

// A fixed rotation (axis (1,2,3)/|.|, angle 0.7 rad) via Rodrigues, so the test
// cages are NOT already sitting in their own principal frame -- an identity
// answer would then be the guard fallback, not the real frame.
matrix3d fixed_rotation() {
  const coord3d k = coord3d(1, 2, 3) / coord3d(1, 2, 3).norm();
  const double th = 0.7, c = cos(th), s = sin(th);
  matrix3d R;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) R(i, j) = (i == j ? c : 0) + (1 - c) * k[i] * k[j];
  R(0, 1) += -s * k[2]; R(0, 2) += s * k[1];
  R(1, 0) += s * k[2];  R(1, 2) += -s * k[0];
  R(2, 0) += -s * k[1]; R(2, 1) += s * k[0];
  return R;
}

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
void rotate_out_of_frame(Polyhedron& P) {
  const matrix3d R = fixed_rotation();
  for (node_t u = 0; u < P.N; u++) P.points[u] = R * P.points[u];
}

// An icosahedral cage is isotropic -- both weightings give a multiple of the
// identity, whose eigenvectors mean nothing.  Squash it (distinct factors on all
// three axes) so the eigenvalues separate and the frame becomes comparable.
void squash(Polyhedron& P) { P.scale(coord3d(1.0, 1.15, 1.35)); }

}  // namespace

// ---------------------------------------------------------------------------
// The volume: the moment kernel and volume_tetra() integrate the same solid.
// ---------------------------------------------------------------------------

// Both integrate the centroid triangulation of the faces, so they describe the
// same solid and must agree to roundoff.
TEST(VolumeMoments, DodecahedronVolumeMatchesTetra) {
  Polyhedron P = Polyhedron::C20();
  const VolumeMoments m = P.volume_moments();

  EXPECT_GT(m.volume, 0);
  EXPECT_NEAR(m.volume, P.volume_tetra(), 1e-9 * P.volume_tetra());
  // Centred at the origin by construction.
  EXPECT_LT((m.centroid - coord3d(0, 0, 0)).norm(), 1e-12 * pow(m.volume, 1.0 / 3));
}

TEST(VolumeMoments, C60VolumeMatchesTetra) {
  Polyhedron P = make_C60();
  const VolumeMoments m = P.volume_moments();

  EXPECT_GT(m.volume, 0);
  const double rel = fabs(m.volume - P.volume_tetra()) / P.volume_tetra();
  cout << "  C60: volume_moments = " << m.volume << ", volume_tetra = "
       << P.volume_tetra() << ", rel diff = " << rel << "\n";
  EXPECT_LT(rel, 1e-9);
}

// Scaling the polyhedron by s scales the enclosed volume by s^3 and the second
// moment (an integral of x x^T dV) by s^5.
TEST(VolumeMoments, ScalesAsCube) {
  Polyhedron P = make_C60();
  const VolumeMoments m0 = P.volume_moments();

  const double s = 2.5;
  P.scale(coord3d(s, s, s));
  const VolumeMoments m1 = P.volume_moments();

  EXPECT_NEAR(m1.volume, m0.volume * s * s * s, 1e-12 * m0.volume * s * s * s);
  EXPECT_LT(max_abs_diff(m1.covariance, m0.covariance * pow(s, 5)),
            1e-12 * m0.covariance.norm() * pow(s, 5));
}

// ---------------------------------------------------------------------------
// The two weightings.
// ---------------------------------------------------------------------------

// inertia_matrix(true) is the centred vertex sum, written out by hand here.
TEST(VolumeMoments, VertexBranchIsTheCentredVertexSum) {
  for (Polyhedron P : {Polyhedron::C20(), make_C60()}) {
    vector<coord3d> pts(P.points.begin(), P.points.end());
    const matrix3d oracle = vertex_inertia_centred(pts);
    const matrix3d got = P.inertia_matrix(true);
    // Same arithmetic in a different accumulation order: equal to roundoff.
    EXPECT_LT(max_abs_diff(got, oracle), 1e-12 * oracle.norm())
        << "N=" << P.N << " got=" << got << " oracle=" << oracle;
  }
}

// The point of the 2026-08-07 change: neither weighting depends on where the
// polyhedron sits.  The pre-change origin-based sum, on the same shift, does.
TEST(VolumeMoments, BothBranchesAreTranslationInvariant) {
  Polyhedron P = make_C60();
  const coord3d d(17.3, -5.1, 2.9);

  const matrix3d Iv0 = P.inertia_matrix(false);   // volume
  const matrix3d Ix0 = P.inertia_matrix(true);    // vertex
  vector<coord3d> pts0(P.points.begin(), P.points.end());
  const matrix3d Iold0 = vertex_inertia_about_origin(pts0);

  P.move(d);
  const matrix3d Iv1 = P.inertia_matrix(false);
  const matrix3d Ix1 = P.inertia_matrix(true);
  vector<coord3d> pts1(P.points.begin(), P.points.end());
  const matrix3d Iold1 = vertex_inertia_about_origin(pts1);

  EXPECT_LT(max_abs_diff(Iv1, Iv0), 1e-9 * Iv0.norm()) << "volume branch moved";
  EXPECT_LT(max_abs_diff(Ix1, Ix0), 1e-9 * Ix0.norm()) << "vertex branch moved";

  // The behaviour that changed: the old origin-based tensor is dominated by the
  // offset (parallel-axis term ~ N|d|^2) and is a different tensor entirely.
  EXPECT_GT(max_abs_diff(Iold1, Iold0), Iold0.norm())
      << "the pre-change sum was expected to move under translation";
  cout << "  C60 shifted by |d| = " << d.norm()
       << ": old-frame tensor changed by " << max_abs_diff(Iold1, Iold0)
       << " (||I_old|| = " << Iold0.norm() << "), new branches unchanged\n";
}

// Pinning the eight molecular call sites to vertex_weighted = true restores the
// PRE-2026-08-07 answer only if the polyhedron was already centred when they
// called -- the vertex sum used to be taken about the ORIGIN and is now taken
// about the vertex mean.  All eight do `move_to_origin()` (or, for
// halma-polyhedron.cc, move_to_origin() then principal_axes(true)) immediately
// before, on the same object with no intervening change to the points, so the
// centring subtraction is a no-op there.  "No-op" is not "exact", though:
// move_to_origin() zeroes the mean only to roundoff, and the branch was
// re-expressed as tr(M)Id - M over a second-moment sum rather than the old
// interleaved long double accumulation.  This measures what is left.
TEST(VolumeMoments, PinnedMolecularCallSitesUnchanged) {
  vector<Polyhedron> cages = make_isomers(60, 18);      // #17 is the worst frame case
  cages.push_back(make_C60());

  double worst_rel = 0, worst_sub = 0, worst_raw = 0;
  for (Polyhedron& P : cages) {
    P.move_to_origin();                                 // what every pinned site does first

    const vector<coord3d> pts(P.points.begin(), P.points.end());
    const matrix3d I_old = vertex_inertia_about_origin(pts);   // the pre-change answer
    const matrix3d I_new = P.inertia_matrix(true);             // the pinned answer

    const double rel = max_abs_diff(I_new, I_old) / I_old.norm();
    worst_rel = max(worst_rel, rel);

    // The frame, compared the only way that means anything: a cage with an
    // exactly degenerate pair (many isomers have one by symmetry) hands the two
    // eigenvectors back in either order, which shows up as a spurious 90 deg
    // between individual vectors while the SUBSPACE is identical.
    const matrix3d A_old = I_old.eigensystem().second;
    const matrix3d A_new = P.principal_axes(true);
    for (int i = 0; i < 3; i++)
      worst_raw = max(worst_raw, axis_angle_deg(row(A_old, i), row(A_new, i)));
    for (const vector<int>& cl : degenerate_clusters(I_new.eigensystem().first, 1e-9))
      worst_sub = max(worst_sub,
                      max_principal_angle_deg(projector(A_old, cl), projector(A_new, cl)));
  }
  printf("  pinned sites over %zu centred cages: worst ||dI||/||I|| = %.3g, worst frame "
         "subspace angle = %.3g deg (worst raw per-vector %.3g deg, degenerate pairs)\n",
         cages.size(), worst_rel, worst_sub, worst_raw);

  // Roundoff, not a change of answer: the alignment the molecular programs
  // produce is the one they produced before.
  EXPECT_LT(worst_rel, 1e-14);
  EXPECT_LT(worst_sub, 1e-6);
}

// The default path now needs a face structure where the old vertex sum needed
// none: volume_moments() calls faces().  This pins what the failure looks like
// when the topology cannot supply one.
//
// MEASURED 2026-08-07, and it is NOT clean for one input: on an EDGELESS view
// (N points, deg == 0) faces() returns an empty list without throwing, and
// centroid_triangulation({}) then takes its "already a triangulation" branch
// (vacuously true for zero faces) into orient_triangulation(), which indexes
// tris[0] unconditionally -- SIGSEGV.  That is PRE-EXISTING and not specific to
// the new path: volume_tetra(), volume_divergence() and surface_area() all reach
// it through exactly the same two calls and crash identically on the same input
// (verified by calling volume_tetra() directly on this view).  What the changed
// default does is make inertia_matrix()/principal_axes()/align_with_axes()
// inherit volume_tetra()'s precondition -- a closed surface with at least one
// face -- where before they were pure point arithmetic needing no topology.  The
// kernel's own tris.empty() guard is therefore unreachable THROUGH THE METHOD
// for that input; it still protects direct kernel callers, which is what the
// first case below checks.
TEST(VolumeMoments, SurfacesThatCannotBeIntegrated) {
  matrix3d Id = identity3();

  // (a) POINTS ONLY: faces() is clean and empty; the kernel's guard holds when
  // it is reached.  The method is NOT called here -- see the note above.
  {
    vector<coord3d> pts{{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1}};
    vector<node_t> nbr(8 * 3, -1);
    vector<uint8_t> deg(8, 0);
    PolyhedronView<double> V(8, 3, nbr, deg, pts);

    EXPECT_EQ(V.count_edges(), 0);
    EXPECT_TRUE(V.is_consistently_oriented());
    EXPECT_TRUE(V.faces().empty());

    const VolumeMoments m = volume_moments(std::span<const coord3d>(pts), {});
    EXPECT_EQ(m.volume, 0.0);
    EXPECT_EQ(m.covariance.norm(), 0.0);
  }

  // (b) NaN COORDINATES on a valid topology (is_invalid() == true).
  {
    Polyhedron P = Polyhedron::C20();
    P.points[0] = coord3d(NAN, NAN, NAN);
    ASSERT_TRUE(P.is_invalid());
    const VolumeMoments m = P.volume_moments();
    const matrix3d I = P.inertia_matrix();
    const matrix3d A = P.principal_axes();
    printf("  NaN points:  volume=%g ||I||=%g ||A-Id||=%g\n",
           m.volume, I.norm(), (A - Id).norm());
    EXPECT_EQ(m.volume, 0.0);            // !(NaN > 0) -> the zeroed guard value
    EXPECT_EQ(I.norm(), 0.0);
    EXPECT_LT((A - Id).norm(), 1e-12);
  }

}

// How far the new default moves the answer for the existing
// align_with_axes()/principal_axes() callers, measured rather than assumed.
//
// The quantity that is WELL CONDITIONED is the tensor difference: the two
// weightings answer with the same shape to 5.9% (Frobenius, unit-trace
// normalised) over the 41 cages below.  An eigenVECTOR angle is not
// well conditioned -- inside a near-degenerate eigenplane no individual
// direction is determined at all, and any rotation there is an equally valid
// frame -- so this test gates on
//   (i)  the tensor difference, and
//   (ii) angles that survive the conditioning: the largest principal angle
//        between eigen-SUBSPACES grouped by near-degenerate eigenvalues, and
//        the per-axis angle restricted to axes separated from BOTH neighbours.
//
// MEASURED 2026-08-07 on the first 20 general C60 + 20 C70 isomers + a squashed
// icosahedral C80 (the raw per-cage table is in the commit message):
//   worst tensor difference           5.93%   (C70 #2)
//   worst RAW per-axis angle          6.27 deg (C70 #2)   <- NOT meaningful
//   worst SUBSPACE angle              3.70 deg (C60 #17)
//   worst WELL-SEPARATED axis angle   3.70 deg (C60 #17), over 109 such axes
//   worst LONGEST-axis angle          3.25 deg (C70 #16)
//   icosahedral C80                   0.001 deg (symmetry forces agreement)
// The two worst RAW numbers are indeed degeneracy artifacts: C70 #2 and #7 are
// rugby balls whose two SHORT axes differ by 1.6% and 4.5%, and their 6.2 deg
// swing inside that plane collapses to 0.95 and 0.50 deg once the plane is
// compared as a subspace.  But the disagreement does NOT vanish with the
// degeneracy: C60 #17 turns a WELL-SEPARATED axis (gaps 7.4% and 9.7%) by
// 3.70 deg, C60 #14 turns one with a 10.3% gap by 3.33 deg, and the longest axis
// -- always the best separated one here (gap 0.10-0.52) -- moves up to 3.25 deg.
// So the frames are close, but "< 1 degree" is FALSE for ordinary cages even
// after conditioning is accounted for; sub-degree agreement holds for symmetric
// cages and for near-spherical shells, not in general.  Bounds carry ~60%
// headroom over the measured worst case.
TEST(VolumeMoments, VolumeAndVertexFramesAgreeOnCages) {
  const double kMaxTensorRelDiff = 0.10;    // measured 0.0593
  const double kMaxFrameAngleDeg = 6.0;     // measured 3.70, subspace and well-separated alike
  const double kWellSeparated = 0.05;       // relative eigenvalue gap to BOTH neighbours
  const double kDegenerateTol = 0.05;       // eigenvalues this close span a subspace, not axes

  struct Worst { string name; double v = -1; };
  Worst w_raw, w_sub, w_wellsep, w_long;
  double worst_drel = 0;
  int cages = 0, wellsep_axes = 0;

  auto examine = [&](const string& name, Polyhedron& P, double angle_bound) {
    rotate_out_of_frame(P);   // so an identity frame can only be the guard fallback

    const matrix3d Iv = P.inertia_matrix(false);
    const matrix3d Ix = P.inertia_matrix(true);
    const matrix3d Av = P.principal_axes(false);
    const matrix3d Ax = P.principal_axes(true);

    const matrix3d Id = identity3();
    ASSERT_GT((Av - Id).norm(), 1e-3) << name << ": volume frame is the guard fallback";
    ASSERT_GT((Ax - Id).norm(), 1e-3) << name << ": vertex frame is the guard fallback";
    ASSERT_LT((Av * Av.transpose() - Id).norm(), 1e-9) << name << ": frame not unitary";
    ASSERT_LT((Ax * Ax.transpose() - Id).norm(), 1e-9) << name << ": frame not unitary";

    // Scale-free comparison of the two weightings: both tensors normalised to
    // unit trace, so only the SHAPE of the mass distribution is compared.
    const matrix3d Ivn = Iv * (1.0 / (Iv(0,0) + Iv(1,1) + Iv(2,2)));
    const matrix3d Ixn = Ix * (1.0 / (Ix(0,0) + Ix(1,1) + Ix(2,2)));
    const double drel = (Ivn - Ixn).norm() / Ivn.norm();
    worst_drel = max(worst_drel, drel);
    EXPECT_LT(drel, kMaxTensorRelDiff) << name;

    // Eigenvalues ASCENDING, so index 0 is the smallest moment = the LONGEST
    // geometric axis, and index 2 the shortest.
    const coord3d lam = Ivn.eigensystem().first;
    const coord3d sep = axis_separations(lam);
    cages++;

    double ang[3], raw = 0;
    for (int i = 0; i < 3; i++) {
      ang[i] = axis_angle_deg(row(Av, i), row(Ax, i));
      raw = max(raw, ang[i]);
    }
    if (raw > w_raw.v) w_raw = {name, raw};
    if (ang[0] > w_long.v) w_long = {name, ang[0]};

    // (ii-a) Axes that ARE determined: separated from both neighbours.  A swing
    // here is a real disagreement about the body's shape, not a reparametrised
    // degenerate plane.
    for (int i = 0; i < 3; i++)
      if (sep[i] >= kWellSeparated) {
        wellsep_axes++;
        if (ang[i] > w_wellsep.v) w_wellsep = {name, ang[i]};
        EXPECT_LT(ang[i], angle_bound)
            << name << " axis " << i << " (separation " << sep[i] << ")";
      }

    // (ii-b) SUBSPACES, grouping eigenvalues too close to determine an individual
    // direction: invariant to rotation inside a degenerate eigenplane, which the
    // per-axis angle is not.
    double sub = 0;
    for (const vector<int>& cl : degenerate_clusters(lam, kDegenerateTol))
      sub = max(sub, max_principal_angle_deg(projector(Av, cl), projector(Ax, cl)));
    if (sub > w_sub.v) w_sub = {name, sub};
    EXPECT_LT(sub, angle_bound) << name << " (subspace)";

    // The cages where the raw angle is large enough to be worth reading: print
    // the conditioning alongside it so the two can never be confused again.
    if (raw > 2.0)
      printf("    %-22s lambda %5.3f %5.3f %5.3f  sep %6.4f %6.4f %6.4f  "
             "angles %6.3f %6.3f %6.3f  subspace %6.3f  dI/I %6.4f\n",
             name.c_str(), lam[0], lam[1], lam[2], sep[0], sep[1], sep[2],
             ang[0], ang[1], ang[2], sub, drel);
  };

  for (int N : {60, 70}) {
    vector<Polyhedron> cages_N = make_isomers(N, 20);
    for (size_t i = 0; i < cages_N.size(); i++)
      examine("C" + to_string(N) + " #" + to_string(i), cages_N[i], kMaxFrameAngleDeg);
  }

  // An icosahedral cage, squashed until its eigenvalues separate.  Symmetry
  // forces the two weightings onto the same axes (0.001 deg), which bounds where
  // the low-symmetry disagreement above comes from: not from the weighting as
  // such, but from the shape asymmetry the two weightings resolve differently.
  { Polyhedron P = make_C80(); squash(P); examine("C80-GC(2,0) squashed", P, 0.01); }

  EXPECT_GE(cages, 40) << "sweep did not run";
  EXPECT_GE(wellsep_axes, 100) << "too few determined axes to make the gate meaningful";
  printf("  %d cages, %d well-separated axes | worst dI/I %.4f | RAW %.3f deg (%s, NOT "
         "meaningful) | SUBSPACE %.3f deg (%s) | WELL-SEPARATED %.3f deg (%s) | LONGEST "
         "%.3f deg (%s)\n",
         cages, wellsep_axes, worst_drel, w_raw.v, w_raw.name.c_str(),
         w_sub.v, w_sub.name.c_str(), w_wellsep.v, w_wellsep.name.c_str(),
         w_long.v, w_long.name.c_str());
}
