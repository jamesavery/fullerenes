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
  DelaunayTriangulation D;
  std::vector<double> r;
  bool verbose = false;
  bool trace_jacobian = false;       // record per-step spectrum of J
  int stats_steps = 0, stats_flips = 0, stats_newton_total = 0;
  double stats_final_kappa = 0;
  double stats_extrap_kappa = 0;   // max|kappa| right after endgame extrapolation
  std::vector<double> r_before_extrap;   // last PALC iterate (for diagnostics)

  // Post-convergence verification (per CLAUDE.md invariants I-1, I-2).
  // Populated by solve() after PALC + Newton; meaningful iff stats_final_kappa
  // is small (the κ=0 hypothesis was achieved numerically).
  bool stats_t0_simplicial = false;       // D has no multi-edges/self-loops/bigons
  bool stats_tbar_simple_polygonal = false; // T̄(0) cells all have ≥3 distinct labels
  int  stats_tbar_n_cells = 0;            // number of polygonal cells in T̄(0)

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

  // Returns the 3D coordinates of the n cone points, or an empty vector on
  // failure.  The final triangulation after edge flips is left in D.
  std::vector<coord3d> solve();

  // Full polytope output: vertex positions + 1-/2-skeleton.  Bundles
  //   - solve()                                            (12 cone-point positions)
  //   - polytope_tesselation(D, r, vertex_labels)         (T̄(0): polygonal 2-faces)
  // into one struct so callers receive the complete Alexandrov polytope.
  // Empty positions on failure (PALC didn't converge or post-convergence
  // invariants violated; see stats_*).  `vertex_labels[k]` is the external
  // label for DCEL vertex k; defaults to identity (k → k).
  struct AlexandrovPolytope {
    std::vector<coord3d> positions;       // V(P): 12 cone-point positions in R³
    CanonicalTesselation tesselation;     // T̄(0): polygonal 2-skeleton of P
    bool ok() const { return !positions.empty(); }
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
};
