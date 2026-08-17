#pragma once
#include "fullerenes/owned.hh"
#include "fullerenes/geometry.hh"
#include <functional>
#include <cstdint>
#include <cstdio>

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

// Deltahedron: owned triangulation with 3D vertex coordinates.
// Inherits geometry/optimizer methods from DeltahedronView<double> via Owned.
// Adds optimizer state fields and static factory methods.
class Deltahedron : public Owned<DeltahedronView<double>> {
  using base_t = Owned<DeltahedronView<double>>;
public:
  // Callback for fromExtensionPathOptimized diagnostics.
  using StepCallback = std::function<void(int step, const char* phase, const Deltahedron& D)>;

  // Optimizer state (set by optimize())
  int iterations_used = 0;
  double final_gmax_L = 0;
  double final_angle_relerr = 0;
  int final_n_concave = 0;
  OptResult final_opt_result = OptResult::BUDGET_EXHAUSTED;  // Set by optimize() and fromExtensionPathOptimized
  PipelineDiag diag;        // Set by fromExtensionPathOptimized: pipeline diagnostics
  uint32_t opt_diag_flags = 0;  // Set by optimize(): optimizer-level diagnostic flags (PipelineDiag flag bits)
  int n_energy_evals = 0;
  int n_grad_evals = 0;
  int n_hv_evals = 0;
  FILE* opt_log = nullptr;
  double opt_k_flat = 2.0;
  double opt_k_conv = 0;   // E_conv coefficient for optimize(). Quadratic one-sided penalty:
                            // E_conv = k * sum_v max(0, -h_v)^2.  Zero at any convex geometry,
                            // pushes concave vertices toward convexity.  Not used in the
                            // extension path pipeline (k_conv=0); reflect-optimize loops
                            // handle convexity instead.  E_conv(k=5, softplus) is used only
                            // in the patch optimizer where the analytical Hessian is tractable.
  bool opt_convex_constraint = false;
  bool opt_skip_post_reflect = false;
  double opt_reflect_threshold = 0.05;  // In-loop reflection: every iteration, vertices with
                            // h < -opt_reflect_threshold*L are mirrored through their
                            // neighbour-centroid plane before the step -- the fold-escape
                            // operator that keeps a cold-start descent out of folded
                            // (convex, wrong-basin) minima.  0 disables.
  OptMethod opt_method = OptMethod::LBFGS;  // Optimization method for optimize() (m=10)

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

  // PLY mesh I/O (see polyhedron-io.cc / deltahedron-io.cc). from_ply builds an
  // oriented Polyhedron and requires it to be a triangulation.
  // @post   result.is_consistently_oriented() && all faces are triangles
  // @throws mesh_io_error  as Polyhedron::from_ply, plus NotATriangulation on a non-triangular face
  static Deltahedron from_ply(FILE *file);
  // @post   result == (no write error occurred)
  // @throws mesh_io_error  on a null file
  static bool to_ply(const Deltahedron &D, FILE *file, bool binary=true);

  // Optimize geometry (uses Deltahedron-specific state fields).
  // Returns OptResult indicating why the optimizer stopped.
  OptResult optimize(std::span<const coord3d> initial_geometry, double target_L=0,
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
