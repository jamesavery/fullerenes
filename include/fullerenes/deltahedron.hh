#pragma once
#include "fullerenes/owned.hh"
#include "fullerenes/geometry.hh"
#include <functional>

class Polyhedron;  // forward declaration

namespace buckinverse { struct ExtensionPath; }

enum class OptMethod { CG, LBFGS, STEIHAUG };

// Deltahedron: owned triangulation with 3D vertex coordinates.
// Inherits geometry/optimizer methods from DeltahedronView via Owned<DeltahedronView>.
// Adds optimizer state fields and static factory methods.
class Deltahedron : public Owned<DeltahedronView> {
  using base_t = Owned<DeltahedronView>;
public:
  // Callback for fromExtensionPathOptimized diagnostics.
  using StepCallback = std::function<void(int step, const char* phase, const Deltahedron& D)>;

  // Optimizer state (set by optimize())
  int iterations_used = 0;
  double final_gmax_L = 0;
  double final_angle_relerr = 0;
  int final_n_concave = 0;
  int n_energy_evals = 0;
  int n_grad_evals = 0;
  int n_hv_evals = 0;
  FILE* opt_log = nullptr;
  double opt_k_flat = 2.0;
  double opt_k_conv = 0;
  bool opt_convex_constraint = false;
  bool opt_skip_post_reflect = false;
  OptMethod opt_method = OptMethod::CG;

  // Constructors
  Deltahedron() = default;
  Deltahedron(const TriangulationView& T, const vector<coord3d>& points);
  Deltahedron(const TriangulationView& T, std::span<const coord3d> points);
  Deltahedron(const Polyhedron& P);  // must be a triangulation

  // Replace owned coordinate storage and repoint the span.
  void set_points(std::vector<coord3d> pts) {
    owned_points = std::move(pts);
    repoint();
  }

  // Optimize geometry (uses Deltahedron-specific state fields)
  bool optimize(std::span<const coord3d> initial_geometry, double target_L=0,
                double grad_tol=1e-10, const vector<bool>& fixed={},
                long long max_work=0, double angle_tol=0);
  bool optimize_patch(std::span<const coord3d> initial_geometry,
                      const vector<bool>& free_mask, const vector<bool>& interior_mask={},
                      double target_L=0, int max_iter=150, double grad_tol=1e-6,
                      bool convex_constraint=true);

  // Build from extension path with incremental geometry.
  static Deltahedron fromExtensionPath(const buckinverse::ExtensionPath& ep);

  // Build from extension path with per-step optimization.
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
};
