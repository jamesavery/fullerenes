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
#include "fullerenes/delaunay_polytope.hh"
#include "fullerenes/eisenstein.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/triangulation.hh"
#include <array>
#include <cmath>
#include <functional>
#include <span>
#include <vector>

struct AlexandrovSolver {
  // ValidationStatus: outcome of the post-convergence validator.  solve()
  // always returns reconstructed positions (so failed cases can be
  // inspected); callers must check stats_status == OK before using them
  // as valid Alexandrov-polytope coordinates.
  enum class ValidationStatus {
    OK,                          // valid output: simple, well-formed, convex
    FAIL_KAPPA_NOT_CONVERGED,    // Newton polish didn't reach its κ < 1e-10 target
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
  // populated on this path, and stats_flips reports only the
  // continuation-stage flips (the override's own flips are not
  // counted -- the callable returns no flip count).  Incubation seam
  // for the optimize framework
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
  //
  // The gluing is recorded per face SLOT, not per vertex pair: T is a delta-complex
  // (multi-edges, self-loops, bigon faces), so "the face across edge (i,j)" is
  // ambiguous from labels alone.  Slot s of face f is its base edge from corner s to
  // corner s+1 (mod 3); face_twin[f][s] = 3*g + s' names the face g and slot s' glued
  // to it (the DCEL twin half-edge's face and cycle slot), and face_theta[f][s] is
  // the GCP dihedral θ at that base edge.  With r and face_len every pyramid is a
  // rigid tetrahedron and the complex is reconstructible exactly, whether or not
  // positions is available -- the seam for drawing the abstract GCP (the
  // delaunay-geometry GUI cuts every radial face and spreads the defect).
  bool record_trajectory = false;
  struct TrajEntry {
    char   phase;                          // 'T' continuation, 'N' Newton polish
    int    step;
    double t;                              // homotopy parameter (0 in Newton)
    double kappa_max;                      // max|κ|
    std::vector<double> kappa;             // per-cone angle defect κ_v
    std::vector<double> r;                 // radii r_v (exact; positions may be empty)
    std::vector<coord3d> positions;        // reconstructed cones (apex at origin)
    std::vector<std::array<int,3>> faces;  // triangular faces of T at this step
    std::vector<std::array<double,3>> face_len;  // base edge lengths of faces[f], in its
                                           // vertex order (ℓ_01, ℓ_12, ℓ_20) — the metric,
                                           // needed to build each rigid pyramid exactly
                                           // (Gram-BFS distorts cotree edges, so |pᵢ−pⱼ|
                                           // is unreliable; the radii rᵥ=|pᵥ| are exact).
    std::vector<std::array<int,3>>    face_twin;   // 3*g + s' glued across slot s (−1: none)
    std::vector<std::array<double,3>> face_theta;  // θ at slot s (NaN where undefined)
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
  // The scale s of that rule alone (feasible_step == s * delta): 1 when
  // r + delta is feasible, else FEAS_SAFETY * s_max.  The optimizer
  // framework's step-clip hook.
  // @pre  as feasible_step
  // @post 0 < result <= 1; *clipped iff result < 1
  static double feasible_fraction(const DelaunayTriangulation& T,
                                  const std::vector<double>& r,
                                  const std::vector<double>& delta,
                                  bool* clipped = nullptr);

  // The flip cap of one weighted-Delaunay repair: a GUARD against a scan
  // that does not terminate, never a knob.  A legal-flip sequence
  // terminates (the Bobenko-Izmestiev piecewise-quadratic extension
  // strictly increases at every flip), so the cap is sized to the complex
  // and generous: one repair legitimately flips every cocircular
  // quadrilateral whose tie a step breaks, dozens on the cubic metric's
  // equal-radius start.  ONE spelling, for this library and for the
  // parallel-primitives port.
  static constexpr int flip_cap(int nh) { return 4 * nh; }

  // The outcome of one repair: T is weighted-Delaunay for r with r
  // feasible (Delaunay); r is infeasible for T before or after a flip
  // (Infeasible); a bad edge's diamond does not support a flip, the B-I
  // obstruction and evidence against the weights (Unflippable); the cap
  // was reached with a bad edge left, an unfinished computation and never
  // a verdict on the weights (Budget: T's status latch trips
  // BudgetExceeded, so the run is refused by name).
  enum class Repair { Delaunay, Infeasible, Unflippable, Budget };

  // Repair T to the weighted-Delaunay triangulation of r by legal flips
  // (θ_e > π on a strictly convex diamond), with r required feasible for T
  // before the first flip and after every flip; `flips` counts the flips
  // applied.  A non-finite dihedral is infeasibility, never a bad edge.
  // The continuation's corrector, the polish entry and the polish trials
  // all run this; the optimize framework's cell-resolved model does too.
  // @anchor topology-repair
  // @pre  as gcp-kappa
  // @post Delaunay: T is weighted-Delaunay for r and r ∈ F(T).  Otherwise
  //       T holds the flips applied so far: a caller discards the copy or
  //       refuses the state, never evaluates on it as if admissible.
  static Repair repair(DelaunayTriangulation& T, const std::vector<double>& r,
                       int& flips);

  // The unconditional flip loop: flip every bad edge, a non-finite
  // dihedral included, until none remains or the cap is reached, with NO
  // feasibility gate.  Returns the number of flips performed.  Not on any
  // solve path (those repair); kept as the elementary operation the
  // optimize sub-project's tests exercise directly.
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

  // ── The doubling witness: the exact certificate that a twelve-cone DUAL
  //    metric is a convex polygon doubled along its boundary, i.e. that its
  //    Alexandrov realization is flat (delaunay-geometry/
  //    validated-alexandrov.tex, "Flat realizations", whose verifier this
  //    transcribes onto the solver's own complex).  A flat metric is what
  //    the continuation cannot deliver: the pyramids flatten as the
  //    curvature is driven down, kappa reaches zero exactly and no placement
  //    closes.  The floor state's complex nevertheless holds the witness:
  //    the polygon's twelve boundary edges are edges of the limiting
  //    polytope (the fold), so they persist in the weighted-Delaunay complex
  //    near it, their dihedral angles tending to 0 while the eighteen
  //    diagonals inside the two sheets tend to pi.
  //
  //    The verdict is the PROPOSAL's (twelve edges), never the state's, and
  //    nothing numerical enters it.  The twelve must form one cycle through
  //    the twelve cones whose complement is two sheets of ten faces.  Each
  //    sheet is developed in the Eisenstein lattice from its faces' integer
  //    squared lengths -- the dual metric's lengths are square roots of
  //    lattice norms, so each rounds to its norm, and place_third_eis_total
  //    places every face's apex across the sheet's interior edges, the root
  //    edge's few unit-orbit representatives (Sector0Reps) tried in turn --
  //    and every occurrence of a cone must land at one position (a cone
  //    inside the sheet would show as a holonomy).  The developed boundary,
  //    read with the sheet on its left, must turn left by exactly 30 degrees
  //    at every corner, an integer identity on consecutive edge vectors u, v:
  //    wedge(u, v) > 0, dot2(u, v) > 0, dot2(u, v)^2 = 3 |u|^2 |v|^2.  A flat
  //    disk without interior cone whose boundary is straight between twelve
  //    corners of interior angle 150 degrees is a convex polygon with those
  //    corners; the other sheet has the complementary 150 degrees at every
  //    cone (the cone angle being 300 degrees) and the same sides in the
  //    same order, hence is the congruent polygon; the gluing along the
  //    cycle is the complex's own.  A proposal that fails proves nothing.
  //
  //    A doubled polygon is a CORRECT realization of its metric (Alexandrov's
  //    theorem includes the degenerate ones) that no consumer can use as a
  //    three-dimensional geometry: callers report it by name, never as a
  //    delivered polytope and never as a solver failure. ──
  enum class DoublingVerdict {
    Witness,               // the twelve edges bound a doubled convex polygon
    NotTwelveCones,        // the complex is not a twelve-cone complex
    NotTwelveEdges,        // a repeated, dead or loop edge in the proposal
    NotOneCycle,           // a cone is not on the cycle exactly twice, or the
                           // cycle does not visit each cone once
    NotTwoSheets,          // the faces off the cycle are not two connected
                           // components of ten faces
    LengthNotANorm,        // a squared edge length rounds to no lattice norm
                           // (not the dual metric)
    NoLatticeDevelopment,  // no root orbit develops a sheet on the lattice
    SheetNotFlat,          // two routes to a face, or two occurrences of a
                           // cone, disagree: a cone inside the sheet
    NotThirtyDegrees,      // a corner does not turn left by exactly 30 degrees
    NoFoldFound,           // the SEARCH below examined every fold the complex
                           // contains and none verified
  };
  static const char* doubling_verdict_str(DoublingVerdict v);

  struct DoubledPolygon {
    DoublingVerdict verdict = DoublingVerdict::NotTwelveEdges;
    bool ok() const { return verdict == DoublingVerdict::Witness; }
    std::array<int, 12> boundary{};        // the fold: half-edges with the first sheet on
                                           // their left, in cycle order
    std::array<int, 12> corner{};          // the cones in that order (corner[k] = origin
                                           // of boundary[k])
    std::array<Eisenstein, 12> polygon{};  // the corners' exact lattice positions in the
                                           // first sheet's development
    std::array<int, 10> sheet_a{}, sheet_b{};   // the faces of the two sheets
    long long area = 0;                    // of the polygon, in unit triangles: half the
                                           // dual's face count
  };
  // FLATNESS IS A PROPERTY OF THE METRIC, and every condition above is a
  // condition on the complex: the cycle and the two sheets are combinatorial,
  // the developments and the 30-degree identity are exact integer arithmetic
  // on the edge lengths.  Nothing in the certificate mentions radii,
  // coordinates or a solve.  So the input complex decides it, by SEARCH
  // rather than by proposal: the folds a twelve-cone complex contains are
  // the twelve-edge cycles through its cones, of which the degree-two and
  // single-cycle conditions leave a few dozen, each verified exactly.  On
  // the five fullerenes known to be flat (C96, C120, C170, C180 and the IPR
  // C384) the fold is already an edge set of the INTRINSIC DELAUNAY complex
  // -- no solve, and no flipping toward it, is needed -- and in each case
  // exactly one of the candidates verifies.
  //
  // Prefer this entry point: its verdict cannot be lost to a solve that
  // refused before reaching a state near the flat limit, and it is the same
  // in every scalar tier.  Exhaustive recognition of flatness is still not
  // claimed -- a fold needing a refinement the complex does not contain
  // would be reported NoFoldFound, never as a wrong answer.
  // @anchor alexandrov-doubled-polygon-search
  // @pre  D is a live twelve-cone complex of the DUAL metric
  // @post result.ok() only if some fold of D verified; then the polygon is
  //       that fold's, and area == (live faces of D) / 2
  static DoubledPolygon doubled_polygon(const DelaunayView& D);

  // The proposal from a state, verified: the twelve live non-loop edges of
  // smallest dihedral angle theta at (D, r), a non-finite theta sorting last.
  // Kept for diagnostics -- it says whether the state a solve reached is the
  // fold -- but the search above is what decides flatness.
  // @anchor alexandrov-doubled-polygon
  // @pre  D is a live complex of the DUAL metric (every live edge length is
  //       the square root of an Eisenstein norm); r.size() >= D.nv
  // @post result.ok() only if the twelve edges verify as in the paragraph
  //       above; then result.area == (number of live faces of D) / 2 and
  //       result.corner is a permutation of the twelve cones
  static DoubledPolygon doubled_polygon(const DelaunayView& D, std::span<const double> r);
  // The verification of a proposal, edge e being the half-edges 2e and 2e+1.
  // @pre  as alexandrov-doubled-polygon (no radii needed)
  static DoubledPolygon verify_doubled_polygon(const DelaunayView& D,
                                               std::span<const int, 12> edges);

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

  // ------ Post-convergence polytope validation (the realize gate) ------

  // Per-check record of the validation ladder; the instance stats_* fields
  // mirror these after solve().
  struct PolytopeValidation {
    bool   t0_simplicial         = false;
    bool   tbar_simple_polygonal = false;
    int    tbar_n_cells          = 0;
    double volume_norm           = 0;
    bool   no_self_intersect     = false;
    bool   convex                = false;
  };

  // The three-property polytope gate on an arbitrary realized (T, r, pos):
  //   SIMPLICITY       T̄(0) is a simple polygonal tesselation with F ≥ 3;
  //   WELL-FORMEDNESS  vol/⟨ℓ⟩³ ≥ 0.01 and no 3D self-intersection;
  //   CONVEXITY        every non-face vertex inside every face plane.
  // This IS the gate solve() applies — its instance path delegates here and
  // copies the per-check record into stats_* — exposed so a caller holding
  // a (T, r, pos) produced elsewhere (e.g. a device batch solve) can apply
  // the identical acceptance rather than re-deriving the ladder.
  // Checks that an early failure skips are left at their defaults in `out`.
  // @pre  κ ≈ 0 (a converged solve)
  // @throws std::invalid_argument unless r.size() >= T.nv and
  //         pos.size() == T.nv (enforced up front; the ladder indexes both)
  static ValidationStatus validate_polytope(const DelaunayTriangulation& T,
                                            const std::vector<double>& r,
                                            const std::vector<coord3d>& pos,
                                            PolytopeValidation* out = nullptr,
                                            bool verbose = false);
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
//   2. Flat-vertex removal, in one of two regimes (a regime = which
//      implementation of the predicates decides):
//        build(T)         THE DEFAULT, exact: remove_and_complete_cyclotomic_kis
//                         -- flatness by the integer cone-excess count, every
//                         other predicate the exact sign of one element of
//                         Z[2cos(pi/15)] (fullerenes/cyclotomic.hh,
//                         delaunay_cyclotomic.hh), then the canonical
//                         completion, so solver.D is the canonical
//                         triangulation of the kis surface: a function of
//                         the labelled graph, not of the order in which flips
//                         are performed (cyclotomic-idt.tex, "Invariance").
//        build_banded(T)  DelaunayTriangulation::compute(K, edge_length_fn,
//                         FLAT_TOL, &new_to_old): floating-point predicates
//                         with a tolerance band and the flatness tolerance.
//                         No canonical completion.  Kept as the reference
//                         the parallel (CPU/GPU) implementation of this
//                         pipeline still mirrors, until that implementation
//                         carries the exact policy too; measured on every
//                         isomer through C100 to reach the same reduction as
//                         the exact regime (claude-projects/delaunay/tools/
//                         bench_cubic_regimes).
//      Both compact to cone ids and label the cones the same way.
//   3. AlexandrovSolver on the resulting 20..60-cone iDT, unmodified.
//
// NOTE: the metric is not Loeschian (pentagon geometry brings sqrt(5)), so
// the integer-exact Eisenstein machinery does not apply to solver.D; the
// exact regime is the cyclotomic one above.  After either build only the
// float predicates (Diamond, cocircular_edges(tol)) can be asked of
// solver.D: the exact algebraic lengths exist only during the build, and
// the double lengths that survive in solver.D do not determine them.
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
  // FLOAT SHADOW of the exact authority: 25*R5^2 == the ring constant
  // Real30::lsq_pentagon_spoke() (fullerenes/cyclotomic.hh); the two are
  // gated at double rounding by test_cyclotomic_algebra's [I] bridge.
  static inline const double R5 = 0.5 / std::sin(M_PI / 5);

  // The n-generic B-I solver.  Configure knobs before build()/solve(),
  // read stats/r/D after; build() only replaces solver.D.
  AlexandrovSolver solver;

  // Configure-before-build knob: when set, the flat-vertex removal inside
  // build() / build_banded() TRACKS every removed kis vertex (DelaunayTriangulation point
  // tracker), so after solve() the kappa=0 solver.D carries each hexagon
  // center, pentagon center and hexagon-only cubic vertex as a
  // (cell, barycentric) location on the cubic polytope's surface --
  // transported through the removal and all homotopy flips.  Tracker
  // labels are kis-complex ids (see KisMetric: id < T.N = dual vertex /
  // cubic face center; id T.N + t = cubic vertex t in T.triangles() order).
  bool track_removed = false;

  // Cone bookkeeping, filled by build() and build_banded().  Cone i is
  // vertex i of solver.D.
  std::vector<tri_t> cone_triangle;  // dual triangle (CCW, T labels) = the cubic vertex
  std::vector<int>   cone_npent;     // k = #pentagon corners (deg-5 T vertices), 1..3
  std::vector<int>   cone_kis_vertex; // kis id of cone i (= T.N + its triangle index)

  // The same three, as one value: what cone_labels computes and both
  // builds install.
  struct ConeLabels {
    std::vector<tri_t> triangle;
    std::vector<int>   npent;
    std::vector<int>   kis_vertex;
  };

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
  static KisMetric kis_metric(const TriangulationView& T);

  // kis -> exact flat removal and canonical completion -> compaction ->
  // cone iDT in solver.D + cone labels.  Returns the completion's counts.
  // @anchor cubic-build
  // @pre  as cubic-kis-metric (T is an oriented fullerene dual)
  // @post cone_triangle.size() == solver.D.nv &&
  //           cone_npent.size() == solver.D.nv
  // @post kappa_is_k_pi_15: all_of(indices(solver.D.nv), [&](int i){
  //           return fabs((2*M_PI - solver.D.v_cone_angle[i])
  //                       - cone_npent[i]*M_PI/15) <= KAPPA_TOL; })
  // @post gauss_bonnet: |sum_i (2*M_PI - v_cone_angle[i]) - 4*M_PI|
  //           <= TOTAL_KAPPA_TOL
  // @post every cocircular cell is triangulated canonically: the fan from
  //       its unique least-rotation corner, or, when exactly two corners tie,
  //       the symmetric split about the diameter joining them (result.fanned
  //       and result.periodic_completed count the two).  Only a word with
  //       three or four least corners, or a failed disk gate, is refused --
  //       result.ambiguous / result.nondisk, neither of which can occur on a
  //       cone surface
  // @throws std::logic_error when a cone guard trips (a flat vertex
  //         survived removal, kappa != k*pi/15, or total curvature != 4pi
  //         -- all "can't happen" on a correct kis metric);
  //         std::runtime_error when the exact predicates refuse a decision
  //         or a completion invariant trips (the throw set of
  //         remove_and_complete_cyclotomic_kis)
  DelaunayView::CompletionStats build(const TriangulationView& T);

  // The tolerance-based regime of the same build (class banner, step 2):
  // DelaunayTriangulation::compute under the float predicates, no
  // completion, then the same compaction and cone labelling.
  // @anchor cubic-build-banded
  // @pre  as cubic-build
  // @post cone_triangle.size() == solver.D.nv && cone_npent.size() ==
  //       solver.D.nv; kappa_is_k_pi_15 and gauss_bonnet as cubic-build
  // @throws std::logic_error as cubic-build's cone guards
  void build_banded(const TriangulationView& T);

  // build(T) + solver.solve().
  // @anchor cubic-solve
  // @pre  as cubic-build
  // @post result.size() == solver.D.nv; result[i] is cone i's position,
  //       valid as an Alexandrov polytope iff solver.valid()
  std::vector<coord3d> solve(const TriangulationView& T);

  // build(T) + solver.solve_polytope(): positions + Tbar(0) (cells labeled
  // by cone index) + validation status, mirroring the base solver's API.
  // @anchor cubic-solve-polytope
  // @pre  as cubic-build
  // @post on ok: result.positions.size() == solver.D.nv
  AlexandrovSolver::AlexandrovPolytope solve_polytope(const TriangulationView& T);

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
  FlatFaceCensus flat_face_census(const TriangulationView& T,
                                  const CanonicalTesselation& tess) const;

 private:
  // The cone labels of the compacted complex D: cone i is kis vertex
  // new_to_old[i], a cubic vertex T.N + t with its dual triangle and its
  // pentagon-corner count k.  Pure.
  // @throws std::logic_error when a face centre survived the removal
  static ConeLabels cone_labels(const KisMetric& M, const DelaunayTriangulation& D,
                                const std::vector<int>& new_to_old);
  // The curvature quantum: every cone's double-precision curvature is
  // k*pi/15 within KAPPA_TOL and they sum to 4*pi within TOTAL_KAPPA_TOL.
  // @throws std::logic_error otherwise
  static void check_curvature_quantum(const DelaunayTriangulation& D,
                                      const ConeLabels& L);
  void set_cones(ConeLabels L);
};
