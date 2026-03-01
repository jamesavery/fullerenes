#pragma once
#include "fullerenes/triangulation.hh"
#include "fullerenes/geometry.hh"

class Polyhedron;  // forward declaration

namespace buckinverse { struct ExtensionPath; }

class Deltahedron : public Triangulation {
public:
  vector<coord3d> points;
  int iterations_used = 0;  // Set by optimize()

  // Constructors
  Deltahedron() = default;
  Deltahedron(const Triangulation& T, const vector<coord3d>& points);
  Deltahedron(const Polyhedron& P);  // must be a triangulation

  // Build from extension path with incremental geometry.
  // Uses tridiagonal Laplacian to place strip vertices at each step.
  static Deltahedron fromExtensionPath(const buckinverse::ExtensionPath& ep);

  // Build from extension path with per-step CG optimization.
  // Same as fromExtensionPath, but calls optimize() after each expansion
  // step to relax geometry before the next strip is placed.
  // max_iter_per_step: CG iterations per expansion step (default 200).
  static Deltahedron fromExtensionPathOptimized(const buckinverse::ExtensionPath& ep, int max_iter_per_step = 200);

  // Faces from triangles on the fly (no storage)
  vector<face_t> compute_dual_faces() const;

  // Laplacian smoothing
  void smooth(double q);

  // GC transform with 3D coordinates (l=0 for now)
  Deltahedron GCtransform(unsigned k, unsigned l=0) const;

  // Optimize geometry toward equilateral triangles.
  // Replaces this->points with optimized coordinates.
  // target_L: desired edge length (0 = compute from mean of initial edges).
  // max_iter: maximum CG iterations (default 3000).
  // grad_tol: gradient norm convergence tolerance (default 1e-12).
  //           Phase 2 exits early when gradient drops 100x from phase start.
  //           Phase 3 exits when gradient norm < grad_tol.
  //           Use ~1e-6 for intermediate relaxation steps.
  // Returns true if converged.
  bool optimize(const vector<coord3d>& initial_geometry, double target_L = 0, int max_iter = 3000, double grad_tol = 1e-12);

  // Finite-difference gradient check. Returns max relative error.
  // Uses the given geometry (or this->points if empty).
  double gradient_check(const vector<coord3d>& geometry, double target_L = 0, double eps = 1e-6) const;
};
