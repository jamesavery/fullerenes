#pragma once

// Bobenko-Izmestiev algorithm for the Alexandrov embedding of a
// polyhedral metric on S^2 with 12 cone points (fullerene iDT).
//
// Given an intrinsic Delaunay triangulation T with geodesic edge lengths,
// find the unique convex polyhedron P in R^3 whose boundary is isometric to T.
//
// Algorithm (Bobenko-Izmestiev with PALC and extrapolation endgame):
//   1. Initialize equal radii r = R*1
//   2. PALC: trace the homotopy F(t,r) = kappa(r) - t*kappa1 = 0
//      from t=1 toward t=0, with adaptive arc-length steps and edge flips
//   3. Endgame: extrapolate r to t=0 via Lagrange interpolation
//   4. Polish: trust-region Newton on kappa(r) = 0
//   5. Reconstruct: place vertices in R^3 from converged (T, r)

#include "fullerenes/delaunay.hh"
#include "fullerenes/geometry.hh"
#include <vector>

struct AlexandrovSolver {
  // ValidationStatus: outcome of the post-convergence validator.  solve()
  // always returns reconstructed positions (so failed cases can be
  // inspected); callers must check stats_status == OK before using them
  // as valid Alexandrov-polytope coordinates.
  enum class ValidationStatus {
    OK,                          // valid output: simple, well-formed, convex
    FAIL_KAPPA_NOT_CONVERGED,    // |κ| > 0.01: PALC + Newton didn't reach κ = 0
    FAIL_NOT_SIMPLE,             // T̄(0) has multi-edges in cells, or F < 3
    FAIL_RECONSTRUCT,            // Gram-BFS yielded NaN (negative perp²)
    FAIL_VOLUME_DEGENERATE,      // vol_norm < 0.01 (drum-cap or near-flat)
    FAIL_SELF_INTERSECTING,      // two non-adjacent triangles cross in 3D
    FAIL_NOT_CONVEX,             // some vertex on outside of a face plane
  };
  static const char* status_str(ValidationStatus s);

  DelaunayTriangulation D;
  std::vector<double> r;
  // Optional override of the initial radii for PALC.  When non-empty
  // and of size D.nv, solve() uses this in place of the default
  // PALC::initial_radii(D) (= 2·R_max·1).  Used by the symmetry-
  // breaking perturbation experiments and by warm-start callers.
  std::vector<double> r_init_override;
  // Coefficient for the strict-interior margin schedule applied during
  // PALC: at homotopy parameter t ∈ (0, 1], iterates are kept in P(M)
  // with safety distance c·t from ∂P(M) on the (ConcQuadr) and
  // (CloGeod) sides.  At t = 0 (Newton polish), margin → 0 reproduces
  // the closure check.  c = 0 disables margin enforcement (pre-#28
  // behaviour).  Tunable for experimentation.
  double palc_interior_margin = 0.0;
  // Stochastic-perturbation experiment (Direction 1): after each accepted
  // PALC step, multiplicatively perturb r per-vertex by uniform random
  // factor in [1−ε, 1+ε] before the next predictor.  Continuously breaks
  // any symmetry the metric induces in the trajectory.  ε = 0 (default)
  // disables the perturbation.  `stochastic_seed` makes runs reproducible.
  double stochastic_perturbation_eps = 0.0;
  uint32_t stochastic_seed = 1;
  // Synchronized batch multi-flip experiment (Direction 2): after each
  // accepted PALC step (post the standard flip_to_delaunay phase),
  // perform a batch flip of all alive non-bigon edges with
  // θ_e > π − batch_multiflip_threshold, all in one pass without
  // re-evaluating θ between flips.  Then a final flip_to_delaunay
  // cleans up.  Threshold = 0 disables.
  double palc_batch_multiflip_threshold = 0.0;
  // Deflated homotopy (Direction 4): if r_deflate_target.size() == D.nv
  // and deflate_strength > 0, augment the PALC residual with a repulsive
  // deflation term that diverges at r_deflate_target — the empirically-
  // known drum-cap fixed point from a prior failed run.
  //   F_def(r, t) = κ(r) − t·κ₁ + (1−t)·α·(r − r*) / ‖r − r*‖²
  // r_deflate_target = r*, deflate_strength = α.  At t = 1 the deflation
  // is off (factor (1−t) = 0); at t = 0 it dominates near r* and forces
  // PALC away from drum-cap.  After the PALC track, an un-deflated
  // Newton polish gives a clean κ = 0 solution.  Inactive when target
  // is empty or strength = 0.
  std::vector<double> r_deflate_target;
  double deflate_strength = 0.0;

  // Continuation method for the BI homotopy κ(r)=t·κ₁.
  //   NATURAL — t-parameterized continuation (BI eq. 38, dr/dt=J⁻¹κ₁).
  //             The path is monotone in t (BI Thm 5: J non-degenerate,
  //             constant Lorentzian signature for 0<κᵢ<δᵢ), so no arclength
  //             is needed; scale-invariant by construction. This is the
  //             default.
  //   PALC    — pseudo-arclength continuation (legacy; arclength ds²=dt²+‖dr‖²
  //             is scale-dependent). Kept for comparison/benchmarking.
  enum class Continuation { NATURAL, PALC };
  Continuation continuation = Continuation::NATURAL;

  bool verbose = false;
  bool trace_jacobian = false;       // record per-step spectrum of J
  int stats_steps = 0, stats_flips = 0, stats_newton_total = 0;
  double stats_final_kappa = 0;
  double stats_extrap_kappa = 0;   // max|kappa| right after endgame extrapolation
  std::vector<double> r_before_extrap;   // last PALC iterate (for diagnostics)
  ValidationStatus stats_status = ValidationStatus::FAIL_KAPPA_NOT_CONVERGED;
  bool valid() const { return stats_status == ValidationStatus::OK; }

  // Post-convergence verification (per CLAUDE.md invariants I-1, I-2).
  // Populated by solve() after PALC + Newton; meaningful iff stats_final_kappa
  // is small (the κ=0 hypothesis was achieved numerically).
  bool stats_t0_simplicial = false;       // D has no multi-edges/self-loops/bigons
  bool stats_tbar_simple_polygonal = false; // T̄(0) cells all have ≥3 distinct labels
  int  stats_tbar_n_cells = 0;            // number of polygonal cells in T̄(0)
  double stats_volume_norm = 0;           // |V| / ⟨ℓ⟩³ on the reconstructed polytope
                                           // (1.03M-scan median ≈ 1.0; degenerate ≤ 1e-6;
                                           //  healthy lower-bound ≈ 0.12)
  bool stats_polytope_convex = false;     // every non-face vertex on inside of every face plane
  bool stats_polytope_no_self_intersect = false; // no two non-adjacent triangles intersect in 3D

  // One entry per PALC step (if trace_jacobian) or Newton step (phase='N').
  struct TraceEntry {
    char   phase;                    // 'P' for PALC, 'N' for Newton polish
    int    step;
    double t;                        // homotopy parameter (PALC only; else 0)
    double ds;                       // PALC arc-length step
    int    nit;                      // corrector iterations (PALC only)
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
    char   phase;                    // 'P' for PALC step, 'N' for Newton iter
    int    step;
    double t;                        // homotopy parameter
    double ds;                       // arc-length step (PALC) / trust radius (Newton)
    int    nit;                      // corrector iters (PALC) / 1 (Newton)
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

  // Returns the 3D coordinates of the n cone points.  ALWAYS returns
  // positions when reconstruction is possible — even on failed
  // validation — so that failure cases can be visualized and debugged.
  // Callers must check `stats_status == OK` (or `valid()`) before
  // treating the result as a valid Alexandrov polytope.  The final
  // triangulation after edge flips is left in D.  Empty result only
  // when the Gram-BFS reconstruction itself yields NaN, in which case
  // stats_status == FAIL_RECONSTRUCT.
  std::vector<coord3d> solve();

  // Full polytope output: vertex positions + 1-/2-skeleton.  Bundles
  //   - solve()                                            (12 cone-point positions)
  //   - polytope_tesselation(D, r, vertex_labels)         (T̄(0): polygonal 2-faces)
  // into one struct so callers receive the complete Alexandrov polytope.
  // Empty positions on failure (PALC didn't converge or post-convergence
  // invariants violated; see stats_*).  `vertex_labels[k]` is the external
  // label for DCEL vertex k; defaults to identity (k → k).
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
  // (e.g. the last PALC iterate of a failed solve).  Returns empty on
  // gross inconsistency (negative perpendicular squared distance).
  static std::vector<coord3d> reconstruct(const DelaunayTriangulation& T,
                                          const std::vector<double>& r);

  // ------ Diagnostic accessors (exposed for verification / unit tests) ------

  // κ(T, r): per-vertex angle deficit at the radial edge a−v.
  static std::vector<double> kappa(const DelaunayTriangulation& T,
                                    const std::vector<double>& r);

  // BI total scalar curvature H(T, r) = Σ_v r_v κ_v + Σ_e ℓ_e (π − θ_e).
  // ∂H/∂r_v = κ_v(r) (BI 2008, Proposition 5, eq. 13).  Hessian is
  // Lorentzian (signature (1, n−1)) by BI Theorem 4 + Lemma 3.4.
  static double H(const DelaunayTriangulation& T,
                   const std::vector<double>& r);

  // Eigenvalues of J = ∂κ/∂r at (T, r), sorted ascending.
  // Empty on LAPACK failure.
  static std::vector<double> jacobian_eigvals(const DelaunayTriangulation& T,
                                               const std::vector<double>& r);

  // sign(det J(T, r)), computed via LU pivot product.  Returns 0 on
  // numerical singularity or LAPACK failure.  Used as a fold detector.
  static int jacobian_det_sign(const DelaunayTriangulation& T,
                                const std::vector<double>& r);

  // Feasibility predicate: r ∈ F(T) iff every incident pyramid closes.
  static bool feasible(const DelaunayTriangulation& T,
                        const std::vector<double>& r);

  // ------ B-I tesselation extraction (Theorem Del=Pol, B-I §3.4) ------

  // GCP dihedral θ_e at base edge h: sum of pyramid dihedrals from the two
  // adjacent faces.  At κ=0, equals the polytope's interior dihedral at e.
  // For a non-bigon edge this is well-defined; for bigon edges (h and h^1
  // in the same face) the result is not geometrically meaningful.
  static double theta(const DelaunayTriangulation& T,
                       const std::vector<double>& r, int h);

  // Per-half-edge "inessential" mask: tight[h] iff |θ_e − π| < ε for the
  // GCP dihedral at edge h.  At κ=0 the inessential edges are precisely the
  // diagonals of flat 2-faces of the polytope (B-I lines 798–820).
  // Both halves of an edge are set consistently.  ε defaults to 1e-7 (rad).
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
  static CanonicalTesselation polytope_tesselation(
      const DelaunayTriangulation& T,
      const std::vector<double>& r,
      const std::vector<int>& vertex_labels,
      double inessential_eps = 1e-7);

  // ------ T(0) / T̄(0) verification (per CLAUDE.md invariants I-1, I-2) ------

  // True iff T contains no multi-edges, self-loops, or bigons — the only
  // T-shape compatible with a non-degenerate Alexandrov polytope.  Per
  // invariant I-1, any T(0) for a fullerene metric MUST be simplicial; a
  // failure here indicates PALC misconvergence or incomplete flip mechanics.
  // O(nh) on the DCEL.
  static bool is_simplicial(const DelaunayTriangulation& T);

  // True iff every polygon in `tess` has ≥ 3 distinct vertex labels and no
  // repeated label on its boundary.  Combined with `is_simplicial(T)` this
  // certifies T̄ is a simple polygonal tesselation.
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
  // O(V·F) + O(F) on V=12, F≤20.
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
  static bool has_self_intersection(const DelaunayTriangulation& T,
                                      const std::vector<coord3d>& pos,
                                      double tol = 1e-6);
};
