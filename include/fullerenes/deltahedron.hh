#pragma once
#include "fullerenes/triangulation.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/span_vector.hh"
#include <functional>

class Polyhedron;  // forward declaration

namespace buckinverse { struct ExtensionPath; }

enum class OptMethod { CG, LBFGS, STEIHAUG };

class Deltahedron : public Triangulation {
public:
  // Callback for fromExtensionPathOptimized diagnostics.
  // Called at each sub-phase: (step_index, phase_name, snapshot_deltahedron).
  // Phases: "seed", "placed", "reflected", "patched", "cg", "final".
  // "patched" = after patch optimize, BEFORE full-graph CG (the key one).
  using StepCallback = std::function<void(int step, const char* phase, const Deltahedron& D)>;
  Spanify::SpanVector<coord3d> points;
  int iterations_used = 0;  // Set by optimize()
  double final_gmax_L = 0;  // Set by optimize(): max_i(||g_i||*L) at final iteration
  double final_angle_relerr = 0;  // Set by optimize(): max per-angle |theta-pi/3|/(pi/3)
  int final_n_concave = 0;        // Set by optimize(): count of vertices with h < 0
  int n_energy_evals = 0;   // Set by optimize(): number of energy-only evaluations
  int n_grad_evals = 0;     // Set by optimize(): number of energy+gradient evaluations
  int n_hv_evals = 0;       // Set by optimize(): number of Hessian-vector products
  FILE* opt_log = nullptr;  // If set, optimize() writes periodic diagnostics here
  double opt_k_flat = 2.0;  // E_flat coefficient for optimize(). Set to 0 to skip phase 1.
  double opt_k_conv = 0;   // E_conv coefficient for optimize(). Quadratic one-sided penalty:
                            // E_conv = k * sum_v max(0, -h_v)^2.  Zero at any convex geometry,
                            // pushes concave vertices toward convexity.  Set to ~10 in the
                            // extension path pipeline.  See CONVEX-BFGS-design.md.
  bool opt_convex_constraint = false;  // Steihaug rejects steps that make convex vertices
                                       // concave (hard constraint). Use with k_conv=0.
  bool opt_skip_post_reflect = false;   // Skip post-optimization reflect_concave + CG polish.
                                        // Use when a subsequent phase will handle convexity.
  OptMethod opt_method = OptMethod::CG;  // Optimization method for optimize()

  // Constructors
  Deltahedron() = default;
  Deltahedron(const Triangulation& T, const vector<coord3d>& points);
  Deltahedron(const Triangulation& T, std::span<const coord3d> points);
  Deltahedron(const Polyhedron& P);  // must be a triangulation

  // Build from extension path with incremental geometry.
  // Uses tridiagonal Laplacian to place strip vertices at each step.
  static Deltahedron fromExtensionPath(const buckinverse::ExtensionPath& ep);

  // Build from extension path with per-step optimization.
  // Same as fromExtensionPath, but calls optimize() after each expansion
  // step to relax geometry before the next strip is placed.
  // max_work_per_step: work budget per optimize() call. 0 = default (20*N^3).
  static Deltahedron fromExtensionPathOptimized(const buckinverse::ExtensionPath& ep, FILE* log = nullptr,
                                                  StepCallback diag = nullptr,
                                                  OptMethod method = OptMethod::LBFGS,
                                                  double step_tol = 1e-4,
                                                  double final_tol = 1e-11,
                                                  long long max_work_per_step = 0,
                                                  double step_angle_tol = 0,
                                                  double final_angle_tol = 0,
                                                  OptMethod final_method = OptMethod::STEIHAUG,
                                                  double patch_grad_tol = 1e-10,
                                                  bool global_post_patch_reflect = false);

  // Quality metrics
  double max_angle_relerr() const;  // max over face angles of |theta - pi/3| / (pi/3)
  int count_concave() const;        // count of vertices with h < 0 (fan-normal based)

  // Faces from triangles on the fly (no storage)
  vector<face_t> compute_dual_faces() const;

  // Laplacian smoothing
  void smooth(double q);

  // GC transform with 3D coordinates
  Deltahedron GCtransform(unsigned k, unsigned l) const;
  Deltahedron halma_transform(int m) const;

  // Optimize geometry toward equilateral triangles.
  // Replaces this->points with optimized coordinates.
  // target_L: desired edge length (0 = compute from mean of initial edges).
  // grad_tol: dimensionless convergence tolerance (default 1e-10).
  //           Convergence is declared when max_i(||g_i|| * L) < grad_tol,
  //           i.e. the largest per-vertex force in dimensionless units is
  //           below this threshold.  This is scale-invariant: the same
  //           tolerance gives the same geometric quality regardless of N or L.
  // max_work: work budget = n_energy + N*n_grad + N*n_hv. 0 = default (20*N^3).
  // angle_tol: if > 0, converge when max_angle_relerr() < angle_tol and no concave vertices.
  // Returns true if converged (vs budget exhaustion).
  bool optimize(std::span<const coord3d> initial_geometry, double target_L = 0,
                double grad_tol = 1e-10,
                const vector<bool>& fixed = {},
                long long max_work = 0,
                double angle_tol = 0);

  // Optimize a small local patch using trust-region Newton.
  // Called on a sub-Deltahedron extracted by extractPatch (O(1) vertices).
  // free_mask[i] = true for vertices that may move; false = fixed boundary.
  // Uses E_bond + E_angle + E_conv (softplus convexity bias) + hard
  // convexity check (h >= -0.05*L for all free vertices).
  // Solves the trust-region subproblem via lambda-search on (H+lambda*I),
  // which naturally handles indefinite Hessians.
  // target_L: desired edge length (0 = compute from boundary edges).
  // Returns true if converged.
  bool optimize_patch(std::span<const coord3d> initial_geometry,
                      const vector<bool>& free_mask,
                      const vector<bool>& interior_mask = {},
                      double target_L = 0,
                      int max_iter = 150,
                      double grad_tol = 1e-6,
                      bool convex_constraint = true);

  // Reflect concave vertices through their neighbor centroid plane.
  // Vertices with h < -threshold are reflected (h = signed height above
  // centroid along fan normal; h < 0 = concave).
  // If fixed is non-empty, fixed[v]=true vertices are skipped.
  // Returns number of vertices reflected.
  int reflect_concave(std::span<coord3d> pts, double threshold = 0,
                      const vector<bool>& fixed = {}) const;

  // Repeatedly reflect concave vertices until none remain (or 20-pass limit).
  // Returns total number of vertices reflected across all passes.
  int reflect_all_concave(std::span<coord3d> pts, double threshold = 0,
                          const vector<bool>& fixed = {}) const;

  // Finite-difference gradient check. Returns max relative error.
  // Uses the given geometry (or this->points if empty).
  double gradient_check(std::span<const coord3d> geometry, double target_L = 0, double eps = 1e-6) const;

  // Finite-difference Hessian check for optimize_patch's analytical Hessian.
  // Compares assemble_patch_hessian against FD of the gradient.
  // Returns max relative error. Prints worst entries if verbose.
  double hessian_check(std::span<const coord3d> geometry,
                       const vector<bool>& free_mask,
                       const vector<bool>& interior_mask = {},
                       double target_L = 0,
                       double eps = 1e-5,
                       bool verbose = false) const;
};
