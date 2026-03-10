#pragma once
#include "fullerenes/triangulation.hh"
#include "fullerenes/geometry.hh"
#include <functional>
#include <cstdint>

class Polyhedron;  // forward declaration

namespace buckinverse { struct ExtensionPath; }

enum class OptMethod { CG, LBFGS, STEIHAUG };

// Exit status from Deltahedron::optimize() and fromExtensionPathOptimized.
// Pinpoints why the optimizer stopped.
enum class OptResult {
  CONVERGED,         // Gradient or angle tolerance met
  STAGNATED,         // Energy stopped decreasing for 50+ iterations
  BUDGET_EXHAUSTED,  // Hit work budget without converging
  CONVEXITY_STUCK    // Reflect-optimize loop failed to eliminate concavity
};

// Human-readable name for OptResult (defined in deltahedron.cc)
const char* opt_result_name(OptResult r);

// Pipeline diagnostics bitfield.
// Tracks which mechanisms fired during fromExtensionPathOptimized, for both
// per-isomer debugging and population-level statistics.  Two packed uint32_t:
//   flags    — boolean events (OR'd across steps/phases)
//   counters — small integers (max/count, packed in 4-bit nibbles)
struct PipelineDiag {
    uint32_t flags = 0;

    // --- Per-step events (OR'd across all expansion steps) ---
    static constexpr uint32_t REFLECT_CYCLING_STEP  = 1u << 0;   // reflect not decreasing concavity
    static constexpr uint32_t HULL_USED_STEP        = 1u << 1;   // hull projection escalation
    static constexpr uint32_t HULL_CYCLING_STEP     = 1u << 2;   // hull also cycling (gave up)
    static constexpr uint32_t CONVEXITY_FAIL_STEP   = 1u << 3;   // step ended with concave vertices
    static constexpr uint32_t PATCH_CYCLING         = 1u << 4;   // patch loop broke early
    static constexpr uint32_t STAG_STEP             = 1u << 5;   // optimize() stagnated during a step

    // --- Final phase events ---
    static constexpr uint32_t REFLECT_CYCLING_FINAL = 1u << 8;   // reflect not decreasing at final
    static constexpr uint32_t HULL_USED_FINAL       = 1u << 9;   // hull escalation at final
    static constexpr uint32_t HULL_CYCLING_FINAL    = 1u << 10;  // hull cycling at final
    static constexpr uint32_t STAG_FINAL            = 1u << 11;  // final reflect-opt stagnated
    static constexpr uint32_t STAG_CONSTRAINED      = 1u << 12;  // constrained Steihaug stagnated
    static constexpr uint32_t BUDGET_CONSTRAINED    = 1u << 13;  // constrained Steihaug hit budget

    // --- Optimizer-level events (OR'd across all optimize() calls) ---
    static constexpr uint32_t NEG_CURVATURE         = 1u << 16;  // Steihaug negative curvature
    static constexpr uint32_t TR_BOUNDARY           = 1u << 17;  // Steihaug trust-region boundary
    static constexpr uint32_t STEP_REJECTED         = 1u << 18;  // Steihaug rejected step (rho < 0.1)
    static constexpr uint32_t CVX_REJECTED          = 1u << 19;  // constrained Steihaug convexity rejection
    static constexpr uint32_t LBFGS_RESET           = 1u << 20;  // L-BFGS history reset

    // --- Structural ---
    static constexpr uint32_t HAS_F_RING            = 1u << 24;  // at least one F-ring step

    // --- Small counters packed into second uint32_t ---
    // Each field is 4 bits (0-15, saturating).
    uint32_t counters = 0;
    //  bits  0- 3: max_patch_rounds across all steps
    //  bits  4- 7: max_reflect_rounds across all steps
    //  bits  8-11: final_reflect_rounds
    //  bits 12-15: n_cvx_fail_steps (saturating count)
    //  bits 16-19: n_stag_steps (saturating count)
    //  bits 20-23: n_step_rejected (Steihaug rejections in constrained phase, sat.)
    //  bits 24-25: seed type (0=C20, 1=C28, 2=C30)
    //  bits 26-27: final result (OptResult enum value)

    // --- Generic counter helpers ---
    void set_counter(int bit_offset, int width, int value) {
        uint32_t mask = ((1u << width) - 1) << bit_offset;
        counters = (counters & ~mask) | (((uint32_t)value << bit_offset) & mask);
    }
    int get_counter(int bit_offset, int width) const {
        return (int)((counters >> bit_offset) & ((1u << width) - 1));
    }
    void max_counter(int bit_offset, int width, int value) {
        int cur = get_counter(bit_offset, width);
        int sat = (1 << width) - 1;
        if (value > sat) value = sat;
        if (value > cur) set_counter(bit_offset, width, value);
    }
    void inc_counter(int bit_offset, int width) {
        int cur = get_counter(bit_offset, width);
        int sat = (1 << width) - 1;
        if (cur < sat) set_counter(bit_offset, width, cur + 1);
    }

    // --- Named accessors ---
    int max_patch_rounds()        const { return get_counter(0, 4); }
    int max_reflect_rounds_step() const { return get_counter(4, 4); }
    int final_reflect_rounds()    const { return get_counter(8, 4); }
    int n_cvx_fail_steps()        const { return get_counter(12, 4); }
    int n_stag_steps()            const { return get_counter(16, 4); }
    int n_step_rejected()         const { return get_counter(20, 4); }
    int seed_type_bits()          const { return get_counter(24, 2); }
    int final_result_bits()       const { return get_counter(26, 2); }

    void set_max_patch_rounds(int v)        { max_counter(0, 4, v); }
    void set_max_reflect_rounds_step(int v) { max_counter(4, 4, v); }
    void set_final_reflect_rounds(int v)    { set_counter(8, 4, std::min(v, 15)); }
    void inc_cvx_fail_steps()               { inc_counter(12, 4); }
    void inc_stag_steps()                   { inc_counter(16, 4); }
    void inc_step_rejected()                { inc_counter(20, 4); }
    void set_seed_type(int v)               { set_counter(24, 2, v); }
    void set_final_result(OptResult r)      { set_counter(26, 2, (int)r); }

    bool any(uint32_t mask) const { return (flags & mask) != 0; }

    // Human-readable flag name for a single flag bit (defined in deltahedron.cc)
    static const char* flag_name(uint32_t single_flag);
};

class Deltahedron : public Triangulation {
public:
  // Callback for fromExtensionPathOptimized diagnostics.
  // Called at each sub-phase: (step_index, phase_name, snapshot_deltahedron).
  // Phases: "seed", "placed", "reflected", "patched", "cg", "final".
  // "patched" = after patch optimize, BEFORE full-graph CG (the key one).
  using StepCallback = std::function<void(int step, const char* phase, const Deltahedron& D)>;
  vector<coord3d> points;
  int iterations_used = 0;  // Set by optimize()
  double final_gmax_L = 0;  // Set by optimize(): max_i(||g_i||*L) at final iteration
  double final_angle_relerr = 0;  // Set by optimize(): max per-angle |theta-pi/3|/(pi/3)
  int final_n_concave = 0;        // Set by optimize(): count of vertices with h < 0
  OptResult final_opt_result = OptResult::BUDGET_EXHAUSTED;  // Set by optimize() and fromExtensionPathOptimized
  PipelineDiag diag;        // Set by fromExtensionPathOptimized: pipeline diagnostics
  uint32_t opt_diag_flags = 0;  // Set by optimize(): optimizer-level diagnostic flags (PipelineDiag flag bits)
  int n_energy_evals = 0;   // Set by optimize(): number of energy-only evaluations
  int n_grad_evals = 0;     // Set by optimize(): number of energy+gradient evaluations
  int n_hv_evals = 0;       // Set by optimize(): number of Hessian-vector products
  FILE* opt_log = nullptr;  // If set, optimize() writes periodic diagnostics here
  double opt_k_flat = 2.0;  // E_flat coefficient for optimize(). Set to 0 to skip phase 1.
  double opt_k_conv = 0;   // E_conv coefficient for optimize(). Quadratic one-sided penalty:
                            // E_conv = k * sum_v max(0, -h_v)^2.  Zero at any convex geometry,
                            // pushes concave vertices toward convexity.  Not used in the
                            // extension path pipeline (k_conv=0); reflect-optimize loops
                            // handle convexity instead.  E_conv(k=5, softplus) is used only
                            // in the patch optimizer where the analytical Hessian is tractable.
  bool opt_convex_constraint = false;  // Steihaug rejects steps that make convex vertices
                                       // concave (hard constraint). Use with k_conv=0.
  bool opt_skip_post_reflect = false;   // Skip post-optimization reflect_concave + CG polish.
                                        // Use when a subsequent phase will handle convexity.
  OptMethod opt_method = OptMethod::LBFGS;  // Optimization method for optimize() (m=10)

  // Constructors
  Deltahedron() = default;
  Deltahedron(const Triangulation& T, const vector<coord3d>& points);
  Deltahedron(const Polyhedron& P);  // must be a triangulation

  // Build from extension path with incremental geometry.
  // Uses tridiagonal Laplacian to place strip vertices at each step.
  static Deltahedron fromExtensionPath(const buckinverse::ExtensionPath& ep);

  // Build from extension path with per-step optimization.
  // Each step: place strip → lift → reflect → patch reflect-optimize loop →
  // full-graph reflect-optimize loop (LBFGS). Final: reflect-optimize (STEIHAUG)
  // + constrained Steihaug polish. See PATCH-GEOMETRY.md for details.
  // max_work_per_step: work budget per optimize() call. 0 = default (400*Nv^2).
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
  // Returns OptResult indicating why the optimizer stopped.
  OptResult optimize(const vector<coord3d>& initial_geometry, double target_L = 0,
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
  bool optimize_patch(const vector<coord3d>& initial_geometry,
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
  int reflect_concave(vector<coord3d>& pts, double threshold = 0,
                      const vector<bool>& fixed = {},
                      vector<bool>* reflected_mask = nullptr) const;

  // Repeatedly reflect concave vertices until none remain (or 20-pass limit).
  // Returns total number of vertices reflected across all passes.
  // If reflected_mask is non-null, sets reflected_mask[v]=true for every vertex
  // that was reflected in any pass (must be pre-sized to N).
  int reflect_all_concave(vector<coord3d>& pts, double threshold = 0,
                          const vector<bool>& fixed = {},
                          vector<bool>* reflected_mask = nullptr) const;

  // Project concave vertices onto the convex hull of the point set.
  // Computes the convex hull, identifies vertices interior to it (concave),
  // and moves each to the nearest point on the hull surface.
  // Preserves graph topology — only coordinates change.
  // Returns number of vertices projected.
  int project_onto_convex_hull(vector<coord3d>& pts) const;

  // Finite-difference gradient check. Returns max relative error.
  // Uses the given geometry (or this->points if empty).
  double gradient_check(const vector<coord3d>& geometry, double target_L = 0, double eps = 1e-6) const;

  // Finite-difference Hessian check for optimize_patch's analytical Hessian.
  // Compares assemble_patch_hessian against FD of the gradient.
  // Returns max relative error. Prints worst entries if verbose.
  double hessian_check(const vector<coord3d>& geometry,
                       const vector<bool>& free_mask,
                       const vector<bool>& interior_mask = {},
                       double target_L = 0,
                       double eps = 1e-5,
                       bool verbose = false) const;
};
