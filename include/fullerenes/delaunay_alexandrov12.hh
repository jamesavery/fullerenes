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

  // Place vertices in R^3 from a triangulation + per-vertex radii using
  // BFS over the Gram matrix entries pos[u].pos[v] = (r_u^2+r_v^2-L_uv^2)/2.
  // Exposed for debugging: lets callers visualise non-converged radii
  // (e.g. the last PALC iterate of a failed solve).  Returns empty on
  // gross inconsistency (negative perpendicular squared distance).
  static std::vector<coord3d> reconstruct(const DelaunayTriangulation& T,
                                          const std::vector<double>& r);
};
