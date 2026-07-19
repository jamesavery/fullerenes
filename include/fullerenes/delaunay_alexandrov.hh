#pragma once

// Alexandrov embeddings of fullerene metrics (Bobenko-Izmestiev).
//
// Given an intrinsic Delaunay triangulation T of a polyhedral cone metric
// on S^2 with geodesic edge lengths, find the unique convex polyhedron P
// in R^3 whose boundary is isometric to T.  Two entry points:
//
//   AlexandrovSolver   -- the B-I solver on any cone iDT (n-generic; the
//                         fullerene dual's 12-cone metric in production).
//   AlexandrovIDTCubic -- the CUBIC polyhedral metric of a fullerene (flat
//                         regular unit pentagons/hexagons, 20..60 cones),
//                         a thin wrapper composing the solver.
//
// Solver algorithm (Bobenko-Izmestiev with continuation and endgame):
//   1. Initialize equal radii r = R*1
//   2. Continuation: trace the homotopy F(t,r) = kappa(r) - t*kappa1 = 0
//      from t=1 toward t=0, with adaptive steps and edge flips
//   3. Endgame: extrapolate r to t=0 via Lagrange interpolation
//   4. Polish: trust-region Newton on kappa(r) = 0
//   5. Reconstruct: place vertices in R^3 from converged (T, r)

#include "fullerenes/delaunay.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/triangulation.hh"
#include <array>
#include <cmath>
#include <functional>
#include <vector>

struct AlexandrovSolver {
  // ValidationStatus: outcome of the post-convergence validator.  solve()
  // always returns reconstructed positions (so failed cases can be
  // inspected); callers must check stats_status == OK before using them
  // as valid Alexandrov-polytope coordinates.
  enum class ValidationStatus {
    OK,                          // valid output: simple, well-formed, convex
    FAIL_KAPPA_NOT_CONVERGED,    // |κ| > 0.01: continuation + Newton didn't reach κ = 0
    FAIL_NOT_SIMPLE,             // T̄(0) has multi-edges in cells, or F < 3
    FAIL_RECONSTRUCT,            // Gram-BFS yielded NaN (negative perp²)
    FAIL_VOLUME_DEGENERATE,      // vol_norm < 0.01 (drum-cap or near-flat)
    FAIL_SELF_INTERSECTING,      // two non-adjacent triangles cross in 3D
    FAIL_NOT_CONVEX,             // some vertex on outside of a face plane
  };
  static const char* status_str(ValidationStatus s);

  DelaunayTriangulation D;
  std::vector<double> r;
  // Optional override of the initial radii.  When non-empty and of size
  // D.nv, solve() uses this in place of the default initial_radii(D)
  // (= 2·R_max·1).  Used by warm-start callers.
  std::vector<double> r_init_override;

  // Optional replacement for the stage-4 trust-region Newton polish.
  // When set, solve() invokes it with the exact post-extrapolation state
  // (D, r) in place of the internal polish; the callable must drive
  // max|kappa(D, r)| toward 0, flipping D to weighted-Delaunay as it
  // moves r (AlexandrovSolver::flip_to_weighted_delaunay), and return
  // whether it converged.  The internal trace/diag recorders are not
  // populated on this path.  Incubation seam for the optimize framework
  // (claude-projects/optimize) to run its re-expressed polish inside the
  // full production pipeline for head-to-head validation.
  std::function<bool(DelaunayTriangulation&, std::vector<double>&)> polish_override;

  // The homotopy κ(r)=t·κ₁ is traced by natural t-continuation (BI eq. 38,
  // dr/dt=J⁻¹κ₁): the path is monotone in t (BI Thm 5: J non-degenerate,
  // constant Lorentzian signature for 0<κᵢ<δᵢ), so no arclength is needed,
  // and the tracker is scale-invariant by construction.  The legacy
  // pseudo-arclength path and its experiment knobs (Directions 1–6) are
  // retired to src/c++/attic/delaunay_alexandrov_palc.cc.attic.

  bool verbose = false;

  bool trace_jacobian = false;       // record per-step spectrum of J
  int stats_steps = 0, stats_flips = 0, stats_newton_total = 0;
  double stats_final_kappa = 0;
  double stats_extrap_kappa = 0;   // max|kappa| right after endgame extrapolation
  std::vector<double> r_before_extrap;   // last continuation iterate (for diagnostics)
  ValidationStatus stats_status = ValidationStatus::FAIL_KAPPA_NOT_CONVERGED;
  bool valid() const { return stats_status == ValidationStatus::OK; }

  // Post-convergence verification (per CLAUDE.md invariants I-1, I-2).
  // Populated by solve() after continuation + Newton; meaningful iff
  // stats_final_kappa is small (the κ=0 hypothesis was achieved numerically).
  bool stats_t0_simplicial = false;       // D has no multi-edges/self-loops/bigons
  bool stats_tbar_simple_polygonal = false; // T̄(0) cells all have ≥3 distinct labels
  int  stats_tbar_n_cells = 0;            // number of polygonal cells in T̄(0)
  double stats_volume_norm = 0;           // |V| / ⟨ℓ⟩³ on the reconstructed polytope
                                           // (1.03M-scan median ≈ 1.0; degenerate ≤ 1e-6;
                                           //  healthy lower-bound ≈ 0.12)
  bool stats_polytope_convex = false;     // every non-face vertex on inside of every face plane
  bool stats_polytope_no_self_intersect = false; // no two non-adjacent triangles intersect in 3D

  // One entry per continuation step (if trace_jacobian) or Newton step
  // (phase='N').
  struct TraceEntry {
    char   phase;                    // 'T' for continuation, 'N' for Newton polish
    int    step;
    double t;                        // homotopy parameter (continuation; else 0)
    double ds;                       // dt step (continuation) / trust radius (Newton)
    int    nit;                      // corrector iterations (continuation only)
    double kappa_max;                // max|kappa|
    double kappa_norm;               // ||kappa||_2
    std::vector<double> eigvals;     // eigenvalues of J, ascending
  };
  std::vector<TraceEntry> trace;

  // Per-step trajectory diagnostic (Task #28 characterization run).  Filled
  // when record_diag=true.  Captures the geometric proximity of (T,r) to
  // ∂P(M) so we can compare drum-cap vs non-degenerate trajectories.
  bool record_diag = false;
  struct DiagEntry {
    char   phase;                    // 'T' for continuation step, 'N' for Newton iter
    int    step;
    double t;                        // homotopy parameter
    double ds;                       // dt step (continuation) / trust radius (Newton)
    int    nit;                      // corrector iters (continuation) / 1 (Newton)
    double kappa_max;                // max|κ|
    // Distance to ∂P(M) on the (ConcQuadr) side: min over non-bigon edges
    // of (π − θ_e).  Small values mean the iterate is close to the
    // boundary; 0 means at the boundary; negative would mean past it.
    double theta_min_dist_to_pi;
    // Counts of edges with (π − θ_e) below thresholds {0.1, 0.01, 1e-3}.
    // Tracks the *distribution* shape: drum-cap convergence has many
    // edges piling up near π simultaneously; non-degenerate convergence
    // has only the flat-face-diagonal edges approaching.
    int    n_near_pi_01, n_near_pi_001, n_near_pi_0001;
    // F(T) margin: smallest Cayley-Menger pyramid h_sq.
    double min_h_sq;
    // Radius spread: std(r) / mean(r).  Drum-cap collapses all radii
    // together; non-degenerate keeps spread.
    double r_cv;
    // Hessian/Jacobian degeneracy: sign(det J).  In P(M) interior with
    // Lorentzian Hessian (1, n−1), this should be a fixed sign.  Sign
    // change indicates a fold / degeneracy.
    int    det_J_sign;
    // Cumulative flips at this step.
    int    n_flips_cum;
    // Number of alive non-bigon edges (denominator for n_near_pi_* ratios).
    int    n_non_bigon_alive;
  };
  std::vector<DiagEntry> diag_trace;

  // Per-step homotopy trajectory, for visualizing the B-I continuation
  // (opt-in via record_trajectory; off by default = zero cost).  One entry per
  // accepted-or-rejected continuation/Newton step, capturing the FULL
  // generalized-convex-polytope (GCP) state needed to draw the deforming
  // pyramids: the reconstructed cone positions (apex at the origin, |p_v| = r_v),
  // the per-cone angle defect κ_v at the radial edge a−v (→ 0 as the pyramids
  // close into a polytope), and the triangular faces of T (which change at edge
  // flips, so each entry is self-contained).  positions is empty for that step
  // iff Gram-BFS reconstruction failed (κ too high / degenerate iterate).
  bool record_trajectory = false;
  struct TrajEntry {
    char   phase;                          // 'T' continuation, 'N' Newton polish
    int    step;
    double t;                              // homotopy parameter (0 in Newton)
    double kappa_max;                      // max|κ|
    std::vector<double> kappa;             // per-cone angle defect κ_v
    std::vector<coord3d> positions;        // reconstructed cones (apex at origin)
    std::vector<std::array<int,3>> faces;  // triangular faces of T at this step
    std::vector<std::array<double,3>> face_len;  // base edge lengths of faces[f], in its
                                           // vertex order (ℓ_01, ℓ_12, ℓ_20) — the metric,
                                           // needed to build each rigid pyramid exactly
                                           // (Gram-BFS distorts cotree edges, so |pᵢ−pⱼ|
                                           // is unreliable; the radii rᵥ=|pᵥ| are exact).
  };
  std::vector<TrajEntry> trajectory;

  // Returns the 3D coordinates of the n cone points.  ALWAYS returns
  // positions when reconstruction is possible — even on failed
  // validation — so that failure cases can be visualized and debugged.
  // Callers must check `stats_status == OK` (or `valid()`) before
  // treating the result as a valid Alexandrov polytope.  The final
  // triangulation after edge flips is left in D.  Empty result only
  // when the Gram-BFS reconstruction itself yields NaN, in which case
  // stats_status == FAIL_RECONSTRUCT.
  // @anchor alexandrov-solve
  // @pre  D is a valid cone iDT: well-formed DCEL, every live vertex has
  //       cone angle < 2π (strictly positive curvature), lengths > 0
  // @post stats_status is one of the ValidationStatus outcomes; on OK the
  //       result is the unique (up to rigid motion) convex polytope
  //       isometric to D's metric, with |result[v] − apex| == r[v]
  // @post result.size() == D.nv, or empty iff
  //       stats_status == FAIL_RECONSTRUCT
  std::vector<coord3d> solve();

  // Full polytope output: vertex positions + 1-/2-skeleton.  Bundles
  //   - solve()                                            (12 cone-point positions)
  //   - polytope_tesselation(D, r, vertex_labels)         (T̄(0): polygonal 2-faces)
  // into one struct so callers receive the complete Alexandrov polytope.
  // Empty positions on failure (continuation didn't converge or
  // post-convergence invariants violated; see stats_*).  `vertex_labels[k]`
  // is the external label for DCEL vertex k; defaults to identity (k → k).
  // @anchor alexandrov-solve-polytope
  // @pre  as alexandrov-solve; vertex_labels empty or of size D.nv
  // @post result.status == stats_status; on ok() the tesselation's cells
  //       are the polytope's flat polygonal 2-faces in external labels
  struct AlexandrovPolytope {
    std::vector<coord3d> positions;       // V(P): 12 cone-point positions in R³.
                                           // Always populated when reconstruction
                                           // succeeded; check `status` for validity.
    CanonicalTesselation tesselation;     // T̄(0): polygonal 2-skeleton of P
    ValidationStatus status = ValidationStatus::FAIL_KAPPA_NOT_CONVERGED;
    bool ok() const { return status == ValidationStatus::OK; }
  };
  AlexandrovPolytope solve_polytope(const std::vector<int>& vertex_labels = {});

  // Place vertices in R^3 from a triangulation + per-vertex radii using
  // BFS over the Gram matrix entries pos[u].pos[v] = (r_u^2+r_v^2-L_uv^2)/2.
  // Exposed for debugging: lets callers visualise non-converged radii
  // (e.g. the last continuation iterate of a failed solve).  Returns empty
  // on gross inconsistency (negative perpendicular squared distance).
  static std::vector<coord3d> reconstruct(const DelaunayTriangulation& T,
                                          const std::vector<double>& r);

  // ------ Diagnostic accessors (exposed for verification / unit tests) ------

  // κ(T, r): per-vertex angle deficit at the radial edge a−v.
  // @anchor gcp-kappa
  // @pre  r.size() == T.nv with r > 0 componentwise
  // @post result.size() == T.nv; result[v] is NaN iff some pyramid at v
  //       fails to close (r outside F(T))
  static std::vector<double> kappa(const DelaunayTriangulation& T,
                                    const std::vector<double>& r);

  // BI total scalar curvature H(T, r) = Σ_v r_v κ_v + Σ_e ℓ_e (π − θ_e).
  // ∂H/∂r_v = κ_v(r) (BI 2008, Proposition 5, eq. 13).  Hessian is
  // Lorentzian (signature (1, n−1)) by BI Theorem 4 + Lemma 3.4.
  // @anchor gcp-H
  // @pre  as gcp-kappa
  static double H(const DelaunayTriangulation& T,
                   const std::vector<double>& r);

  // J = ∂κ/∂r at (T, r): the dense symmetric B-I Jacobian (Hessian of H).
  // Exported so external solvers (e.g. the Nv-vertex generalization) can
  // build on the same Jacobian instead of re-deriving it.
  // @anchor gcp-jacobian
  // @pre  as gcp-kappa
  // @post result is T.nv × T.nv and symmetric up to roundoff; entries are
  //       NaN where the pyramid geometry is degenerate
  static matrix<double> jacobian(const DelaunayTriangulation& T,
                                  const std::vector<double>& r);

  // Eigenvalues of J = ∂κ/∂r at (T, r), sorted ascending.
  // Empty on eigensolver failure.
  // @anchor gcp-jacobian-eigvals
  // @pre  as gcp-kappa
  static std::vector<double> jacobian_eigvals(const DelaunayTriangulation& T,
                                               const std::vector<double>& r);

  // sign(det J(T, r)), computed via LU pivot product.  Returns 0 on
  // numerical singularity.  Used as a fold detector.
  // @anchor gcp-jacobian-det-sign
  // @pre  as gcp-kappa
  // @post result ∈ {−1, 0, +1}
  static int jacobian_det_sign(const DelaunayTriangulation& T,
                                const std::vector<double>& r);

  // Feasibility predicate: r ∈ F(T) iff every incident pyramid closes.
  // @anchor gcp-feasible
  // @pre  r.size() == T.nv
  static bool feasible(const DelaunayTriangulation& T,
                        const std::vector<double>& r);

  // Clip step δ so r + δ' ∈ F(T) strictly: δ' = δ when already feasible,
  // else (0.95 · s_max) · δ with s_max the largest feasible fraction
  // (bisected).  The single source of truth for the F(T)-feasible step
  // rule, used by the endgame extrapolation and the Newton polish; exposed
  // so external polish implementations (polish_override) apply the SAME
  // rule instead of transcribing it.  `clipped` (if non-null) reports
  // whether a scale was applied.
  // @anchor gcp-feasible-step
  // @pre  r ∈ F(T); delta.size() == r.size()
  // @post r + result ∈ F(T) strictly
  static std::vector<double> feasible_step(const DelaunayTriangulation& T,
                                            const std::vector<double>& r,
                                            const std::vector<double>& delta,
                                            bool* clipped = nullptr);

  // Flip T to the weighted-Delaunay triangulation of r (B-I bad-edge
  // rule), up to the internal iteration cap.  Returns the number of
  // flips performed.  Exposed for external polish implementations
  // (polish_override) — same single-source-of-truth rationale as
  // feasible_step.
  // @anchor topology-flip-to-weighted-delaunay
  // @pre  as gcp-kappa
  static int flip_to_weighted_delaunay(DelaunayTriangulation& T,
                                        const std::vector<double>& r);

  // ------ B-I tesselation extraction (Theorem Del=Pol, B-I §3.4) ------

  // GCP dihedral θ_e at base edge h: sum of pyramid dihedrals from the two
  // adjacent faces.  At κ=0, equals the polytope's interior dihedral at e.
  // For a non-bigon edge this is well-defined; for bigon edges
  // (DelaunayTriangulation::is_bigon) the result is not geometrically
  // meaningful.
  // @anchor gcp-theta
  // @pre  h is a live half-edge of T; r as in gcp-kappa
  static double theta(const DelaunayTriangulation& T,
                       const std::vector<double>& r, int h);

  // Per-half-edge "inessential" mask: tight[h] iff |θ_e − π| < ε for the
  // GCP dihedral at edge h.  At κ=0 the inessential edges are precisely the
  // diagonals of flat 2-faces of the polytope (B-I lines 798–820).
  // Both halves of an edge are set consistently.  ε defaults to 1e-7 (rad).
  // @anchor bi-inessential-edges
  // @pre  as gcp-kappa
  // @post result.size() == T.nh; result[h] == result[T.twin(h)]
  static std::vector<bool> inessential_edges(const DelaunayTriangulation& T,
                                              const std::vector<double>& r,
                                              double eps = 1e-7);

  // Polytope tesselation T̄ at (T, r): the cell decomposition obtained by
  // collapsing all inessential (θ = π) edges of T.  This is the polygonal
  // 2-skeleton of the Alexandrov polytope in B-I's framework.  At κ=0 with
  // a non-degenerate polytope, T̄ is a simple polygonal tesselation
  // (every cell is a flat polygon of P, no multi-edges, no self-loops).
  //
  // `vertex_labels[k]` maps DCEL vertex k to an external label (typically
  // the cone point's index in the input dual triangulation).
  // @anchor bi-polytope-tesselation
  // @pre  (T, r) at κ ≈ 0 (a converged solve); vertex_labels.size() == T.nv
  static CanonicalTesselation polytope_tesselation(
      const DelaunayTriangulation& T,
      const std::vector<double>& r,
      const std::vector<int>& vertex_labels,
      double inessential_eps = 1e-7);

  // ------ T(0) / T̄(0) verification (per CLAUDE.md invariants I-1, I-2) ------

  // True iff T contains no multi-edges, self-loops, or bigons — the only
  // T-shape compatible with a non-degenerate Alexandrov polytope.  Per
  // invariant I-1, any T(0) for a fullerene metric MUST be simplicial; a
  // failure here indicates misconvergence or incomplete flip mechanics.
  // Delegates to DelaunayTriangulation::is_simplicial.
  // @anchor bi-is-simplicial
  static bool is_simplicial(const DelaunayTriangulation& T);

  // True iff every polygon in `tess` has ≥ 3 distinct vertex labels and no
  // repeated label on its boundary.  Combined with `is_simplicial(T)` this
  // certifies T̄ is a simple polygonal tesselation.
  // @anchor bi-is-simple-polygonal
  static bool is_simple_polygonal(const CanonicalTesselation& tess);

  // ------ Reconstructed-polytope geometric checks ------

  // True iff `pos` describes a convex non-degenerate polytope.  Two checks:
  //   (a) Defensive precondition: signed volume (in the CCW half-edge
  //       convention) is strictly positive.  Rejects flat (drum-cap, vol=0)
  //       and globally-inverted (vol<0) configurations — for either, the
  //       outward-normal direction inferred from CCW order is wrong, and
  //       the vertex test below would give a meaningless answer.
  //   (b) Convexity test: every vertex v ∉ f has signed distance
  //       ≤ `tol·mean_edge_length` from the plane of f, where the outward
  //       normal is `(b−a) × (c−a)` for three consecutive vertices in
  //       `he_next` order.  No spherical-approximation assumption — works
  //       for nanotubes, oblate polytopes, irregular shapes.  For T̄(0)-
  //       collapsed flat faces, all triangles within a face share a plane
  //       so the per-triangle check is correct.
  // O(V·F) + O(F).
  // @anchor bi-is-convex
  // @pre  pos.size() == T.nv
  static bool is_convex(const DelaunayTriangulation& T,
                          const std::vector<coord3d>& pos,
                          double tol = 1e-3);

  // True iff some pair of non-adjacent triangles in T (sharing no vertex)
  // intersect in 3D.  Möller's triangle-triangle test.  O(F²) face pairs;
  // ≤ 400 pair-tests on V=12.
  //
  // A self-intersecting "polytope" is not embedded in R³ — it's not a
  // valid polytope at all.  This check is a core definitional gate of
  // validity, on equal footing with convexity and simplicity.  Convexity
  // does imply non-self-intersection for closed 2-spheres, so the two
  // checks are not independent on healthy outputs; but the gate must be
  // enforced anyway so that any failure mode is reported under its
  // correct label.
  // @anchor bi-has-self-intersection
  // @pre  pos.size() == T.nv
  static bool has_self_intersection(const DelaunayTriangulation& T,
                                      const std::vector<coord3d>& pos,
                                      double tol = 1e-6);
};

// ============================================================================
// AlexandrovIDTCubic -- Alexandrov embedding of the CUBIC polyhedral metric
// of a fullerene ("flat pentagons instead of sharp cones").
//
// The metric M_cubic: every face of the cubic fullerene graph is a flat
// regular unit-edge polygon (12 pentagons, N/2-10 hexagons), glued
// isometrically along the unit edges.  All curvature sits at the cubic
// vertices: a vertex incident to k pentagons (k in {1,2,3}; in a cubic
// graph, faces sharing a vertex share an edge, so k counts edge-fused
// pentagons) has angle sum (3-k)*2pi/3 + k*3pi/5 = 2pi - k*pi/15, i.e.
//
//     kappa(v) = k * pi/15,     sum over v of k_v = 60  =>  total 4pi.
//
// The cone set is the pentagon-incident cubic vertices: 20 (C20, all k=3)
// up to 60 (no two pentagons fused; all IPR isomers), CONSTANT in N.
// Hexagon-only vertices and all face interiors are flat, so the
// hexagon-removed iDT -- and hence the solver's problem -- never sees
// the hexagons.
//
// Exact relation to the dual deltahedron metric M_dual (one unit
// equilateral triangle per cubic vertex): M_cubic = sqrt(3) * M_dual with
// each of the 12 pi/3 cones intrinsically truncated at geodesic radius 1
// (the cut runs through the barycenters of the 5 incident triangles, five
// unit chords) and capped by a flat regular unit pentagon.  The two
// metrics are isometric outside those 12 constant-size caps; "cubic"
// opens each pi/3 cone into up to 5 cones of pi/15 spread around a flat
// pentagon.
//
// By Alexandrov's theorem M_cubic has a unique convex realization
// P_cubic: a polytope on the 20..60 cones, independent of N.  Its edge
// skeleton is an OUTPUT of the solve (the theta_e = pi collapse Tbar(0)):
// pentagons are intrinsically flat but not guaranteed to be faces of
// P_cubic, and hexagon regions generally crease (a nanotube's tube
// becomes long flat strips).  Closed forms: C20 -> regular dodecahedron,
// C60-Ih -> truncated icosahedron.  Empirically (5917 isomers C20-C100),
// the realization is fully simplicial for all but exactly-symmetric
// isomers, and is rounder than the dual realization for all but C20.
//
// Implementation: a thin wrapper over the n-generic AlexandrovSolver.
//   1. kis subdivision of the cubic graph, built directly from the dual
//      triangulation T: kis vertices 0..Nv-1 = T vertices (face centers),
//      Nv+t = dual triangle t (= cubic vertex); one kis triangle per
//      (face, boundary arc) incidence.  Lengths: cubic-cubic edge 1,
//      hexagon spoke 1, pentagon spoke R5 = 1/(2 sin(pi/5)).  Every face
//      center is flat by construction (5 wedges of 2pi/5, resp. 6 of
//      pi/3), and hexagon-only cubic vertices are flat, so the metric's
//      cones are exactly the pentagon-incident cubic vertices.
//   2. DelaunayTriangulation::compute(K, edge_length_fn, FLAT_TOL,
//      &new_to_old) removes the flat vertices and reports the survivors.
//   3. AlexandrovSolver on the resulting 20..60-cone iDT, unmodified.
//
// NOTE: the metric is not Loeschian (pentagon geometry brings sqrt(5));
// none of the integer-exact Eisenstein machinery applies to solver.D --
// only the float predicates (Diamond, cocircular_edges(tol)) may be used.
// ============================================================================

struct AlexandrovIDTCubic {
  // Named guard tolerances (referenced by the contracts below).
  // FLAT_TOL sits far below the smallest real cone curvature pi/15 ~ 0.209;
  // KAPPA_TOL / TOTAL_KAPPA_TOL bound the float error of the kis metric's
  // angle sums (exact values would be 0 in exact arithmetic).
  static constexpr double FLAT_TOL        = 1e-6;
  static constexpr double KAPPA_TOL       = 1e-9;
  static constexpr double TOTAL_KAPPA_TOL = 1e-8;

  // Circumradius of the flat regular unit-edge pentagon: the kis spoke
  // length from a pentagon center to its corners (hexagon spokes are 1).
  static inline const double R5 = 0.5 / std::sin(M_PI / 5);

  // The n-generic B-I solver.  Configure knobs before build()/solve(),
  // read stats/r/D after; build() only replaces solver.D.
  AlexandrovSolver solver;

  // Cone bookkeeping, filled by build().  Cone i is vertex i of solver.D.
  std::vector<tri_t> cone_triangle;  // dual triangle (CCW, T labels) = the cubic vertex
  std::vector<int>   cone_npent;     // k = #pentagon corners (deg-5 T vertices), 1..3

  // The kis complex of the dual triangulation T, plus its metric.
  struct KisMetric {
    Triangulation      K;         // 0..Nv-1 face centers, Nv+t = triangle t
    std::vector<tri_t> triangle;  // t -> CCW corners in T labels (arc convention)
    std::vector<int>   fdeg;      // T-vertex degree = cubic face size (5 or 6)
    int                Nv = 0;    // T.N
    // The prescribed intrinsic metric of the kis complex:
    // cubic-cubic edge 1, hexagon spoke 1, pentagon spoke R5.
    DelaunayTriangulation::EdgeLengthFn edge_length_fn() const;
  };

  // @anchor cubic-kis-metric
  // @pre  oriented:  T's neighbour rings are consistently CCW-oriented
  // @pre  fullerene: all_of(indices(T.N), [&](int u){
  //           return T.degree(u) == 5 || T.degree(u) == 6; })
  // @post result.triangle.size() == 2*T.N - 4
  // @post result.K.N == T.N + 2*T.N - 4
  // @throws std::logic_error when pre(fullerene) is violated
  static KisMetric kis_metric(const Triangulation& T);

  // kis -> flat removal (lib compute) -> cone iDT in solver.D + cone labels.
  // @anchor cubic-build
  // @pre  as cubic-kis-metric (T is an oriented fullerene dual)
  // @post cone_triangle.size() == solver.D.nv &&
  //           cone_npent.size() == solver.D.nv
  // @post kappa_is_k_pi_15: all_of(indices(solver.D.nv), [&](int i){
  //           return fabs((2*M_PI - solver.D.v_cone_angle[i])
  //                       - cone_npent[i]*M_PI/15) <= KAPPA_TOL; })
  // @post gauss_bonnet: |sum_i (2*M_PI - v_cone_angle[i]) - 4*M_PI|
  //           <= TOTAL_KAPPA_TOL
  // @throws std::logic_error when a cone guard trips (a flat vertex
  //         survived removal, kappa != k*pi/15, or total curvature != 4pi
  //         -- all "can't happen" on a correct kis metric)
  void build(const Triangulation& T);

  // build(T) + solver.solve().
  // @anchor cubic-solve
  // @pre  as cubic-build
  // @post result.size() == solver.D.nv; result[i] is cone i's position,
  //       valid as an Alexandrov polytope iff solver.valid()
  std::vector<coord3d> solve(const Triangulation& T);

  // build(T) + solver.solve_polytope(): positions + Tbar(0) (cells labeled
  // by cone index) + validation status, mirroring the base solver's API.
  // @anchor cubic-solve-polytope
  // @pre  as cubic-build
  // @post on ok: result.positions.size() == solver.D.nv
  AlexandrovSolver::AlexandrovPolytope solve_polytope(const Triangulation& T);

  // Flat-face census of Tbar(0) against the cubic graph's face lattice.
  // Face u of the cubic graph (= vertex of T) is realized FLAT iff all its
  // corners are cones and they form one cell of tess.  (A face with a
  // flat-vertex corner cannot be a cell of the cone tesselation.)
  struct FlatFaceCensus {
    int n_cells   = 0;   // cells of Tbar(0)
    int pent_flat = 0;   // pentagons realized flat (0..12)
    int hex_flat  = 0;   // hexagons realized flat
    int n_hex     = 0;   // hexagon count (Nv - 12)
    // The realization IS the face lattice (every face a flat polygon).
    bool face_lattice() const { return pent_flat == 12 && hex_flat == n_hex; }
  };
  // @anchor cubic-flat-face-census
  // @pre  T is the dual triangulation build() consumed
  // @pre  tess is the solved polytope's Tbar(0) in cone labels
  //       (solve_polytope(T).tesselation)
  // @post result.pent_flat <= 12 && result.hex_flat <= result.n_hex &&
  //           result.n_hex == T.N - 12
  FlatFaceCensus flat_face_census(const Triangulation& T,
                                  const CanonicalTesselation& tess) const;
};
