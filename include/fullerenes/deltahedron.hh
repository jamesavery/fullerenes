#pragma once
#include "fullerenes/triangulation.hh"
#include "fullerenes/geometry.hh"
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
  std::span<coord3d> points;           // view -- always valid when N > 0
  std::vector<coord3d> owned_points;   // owned storage (empty when viewing external data)
  int iterations_used = 0;  // Set by optimize()
  double final_gmax_L = 0;  // Set by optimize(): max_i(||g_i||*L) at final iteration
  double final_angle_relerr = 0;  // Set by optimize(): max per-angle |theta-pi/3|/(pi/3)
  int final_n_concave = 0;        // Set by optimize(): count of vertices with h < 0
  int n_energy_evals = 0;   // Set by optimize(): number of energy-only evaluations
  int n_grad_evals = 0;     // Set by optimize(): number of energy+gradient evaluations
  int n_hv_evals = 0;       // Set by optimize(): number of Hessian-vector products
  FILE* opt_log = nullptr;  // If set, optimize() writes periodic diagnostics here
  double opt_k_flat = 2.0;  // E_flat coefficient for optimize(). Set to 0 to skip phase 1.
  OptMethod opt_method = OptMethod::CG;  // Optimization method for optimize()

  void repoint_coords() { points = std::span<coord3d>(owned_points); }

  void set_points(std::vector<coord3d> pts) {
    owned_points = std::move(pts);
    repoint_coords();
  }

  // Constructors
  Deltahedron() = default;
  Deltahedron(const Triangulation& T, std::span<const coord3d> pts);
  Deltahedron(const Polyhedron& P);  // must be a triangulation

  // Rule of 5 (span needs repointing after copy/move)
  Deltahedron(const Deltahedron& D)
    : Triangulation(D), owned_points(D.owned_points),
      iterations_used(D.iterations_used), final_gmax_L(D.final_gmax_L),
      final_angle_relerr(D.final_angle_relerr), final_n_concave(D.final_n_concave),
      n_energy_evals(D.n_energy_evals), n_grad_evals(D.n_grad_evals),
      n_hv_evals(D.n_hv_evals), opt_log(D.opt_log), opt_k_flat(D.opt_k_flat),
      opt_method(D.opt_method) {
    if (!owned_points.empty()) repoint_coords();
    else points = D.points;
  }
  Deltahedron(Deltahedron&& D) noexcept
    : Triangulation(std::move(D)), owned_points(std::move(D.owned_points)),
      iterations_used(D.iterations_used), final_gmax_L(D.final_gmax_L),
      final_angle_relerr(D.final_angle_relerr), final_n_concave(D.final_n_concave),
      n_energy_evals(D.n_energy_evals), n_grad_evals(D.n_grad_evals),
      n_hv_evals(D.n_hv_evals), opt_log(D.opt_log), opt_k_flat(D.opt_k_flat),
      opt_method(D.opt_method) {
    if (!owned_points.empty()) repoint_coords();
    else points = D.points;
    D.points = {};
  }
  Deltahedron& operator=(const Deltahedron& D) {
    if (this != &D) {
      Triangulation::operator=(D);
      owned_points = D.owned_points;
      iterations_used = D.iterations_used;
      final_gmax_L = D.final_gmax_L;
      final_angle_relerr = D.final_angle_relerr;
      final_n_concave = D.final_n_concave;
      n_energy_evals = D.n_energy_evals;
      n_grad_evals = D.n_grad_evals;
      n_hv_evals = D.n_hv_evals;
      opt_log = D.opt_log; opt_k_flat = D.opt_k_flat;
      opt_method = D.opt_method;
      if (!owned_points.empty()) repoint_coords();
      else points = D.points;
    }
    return *this;
  }
  Deltahedron& operator=(Deltahedron&& D) noexcept {
    Triangulation::operator=(std::move(D));
    owned_points = std::move(D.owned_points);
    iterations_used = D.iterations_used;
    final_gmax_L = D.final_gmax_L;
    final_angle_relerr = D.final_angle_relerr;
    final_n_concave = D.final_n_concave;
    n_energy_evals = D.n_energy_evals;
    n_grad_evals = D.n_grad_evals;
    n_hv_evals = D.n_hv_evals;
    opt_log = D.opt_log; opt_k_flat = D.opt_k_flat;
    opt_method = D.opt_method;
    if (!owned_points.empty()) repoint_coords();
    else points = D.points;
    D.points = {};
    return *this;
  }

  // Build from extension path with incremental geometry.
  // Uses tridiagonal Laplacian to place strip vertices at each step.
  static Deltahedron fromExtensionPath(const buckinverse::ExtensionPath& ep);

  // Build from extension path with per-step CG optimization.
  // Same as fromExtensionPath, but calls optimize() after each expansion
  // step to relax geometry before the next strip is placed.
  // max_iter_per_step: CG iterations per expansion step (default 200).
  // max_iter_per_step: CG iterations per expansion step. 0 (default) = use 2*D.N.
  // After all steps, runs a final 12*N optimization pass.
  static Deltahedron fromExtensionPathOptimized(const buckinverse::ExtensionPath& ep, int max_iter_per_step = 0, FILE* log = nullptr,
                                                  StepCallback diag = nullptr,
                                                  OptMethod method = OptMethod::CG,
                                                  double step_tol = 1e-3,
                                                  double final_tol = 1e-5,
                                                  long long max_work_per_step = 0,
                                                  double step_angle_tol = 0,
                                                  double final_angle_tol = 0);

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
  bool optimize(std::span<const coord3d> initial_geometry, double target_L = 0, int max_iter = 3000, double grad_tol = 1e-10,
                const vector<bool>& fixed = {},
                long long max_work = 0,
                double angle_tol = 0);

  // Optimize a small local patch using trust-region Newton.
  bool optimize_patch(std::span<const coord3d> initial_geometry,
                      const vector<bool>& free_mask,
                      const vector<bool>& interior_mask = {},
                      double target_L = 0,
                      int max_iter = 150,
                      double grad_tol = 1e-6);

  // Reflect concave vertices through their neighbor centroid plane.
  int reflect_concave(std::span<coord3d> pts, double threshold = 0,
                      const vector<bool>& fixed = {}) const;

  // Finite-difference gradient check. Returns max relative error.
  double gradient_check(std::span<const coord3d> geometry, double target_L = 0, double eps = 1e-6) const;

  // Finite-difference Hessian check for optimize_patch's analytical Hessian.
  double hessian_check(std::span<const coord3d> geometry,
                       const vector<bool>& free_mask,
                       const vector<bool>& interior_mask = {},
                       double target_L = 0,
                       double eps = 1e-5,
                       bool verbose = false) const;
};
