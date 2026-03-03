# Deltahedron Geometry Optimization

## Overview

The deltahedron optimizer produces equilateral triangulations of fullerene duals
using conjugate gradient minimization with a two-phase energy function and
periodic vertex reflection for convexity.

## Energy Function

The total energy is:

    E = E_bond + E_angle + E_curv + E_flat

### E_bond (k_bond = 1.0)

Harmonic bond energy. For each edge (i,j):

    E_bond = (k/2) sum_edges (|x_i - x_j| - L)^2

where L is the target edge length (default: mean of initial edge lengths).

### E_angle (k_angle = 1.0)

Angular deviation from 60 degrees. For each triangle with angles theta_0,
theta_1, theta_2:

    E_angle = (k/2) sum_triangles sum_c (theta_c - pi/3)^2

### E_curv (k_curv = 2.0)

Discrete curvature penalty. For each interior edge shared by two triangles,
penalizes deviation of the dihedral angle from 180 degrees (flat). Uses the
cotangent formula for the discrete Laplace-Beltrami operator.

### E_flat (k_flat = 2.0, phase 1 only)

Eigenvalue-based flatness measure from flatness.tex. For each vertex v of
degree d <= 6, compute the centroids of the d surrounding triangles. The
smallest eigenvalue lambda_0 of the 3x3 covariance matrix X^T X (where X is
the centroid-corrected matrix) measures deviation from coplanarity.

    E_flat = (k_flat/2) sum_v (lambda_0(v) / tr(A(v)))

The gradient uses the quotient rule and propagates through the centroid
dependencies. The gradient contribution to vertex v itself cancels by symmetry
(v appears in all centroids equally).

## Convexity: Periodic Vertex Reflection

Convexity is maintained not by an energy term, but by a geometric intervention:
every iteration, check all degree-<=6 vertices for deep concavity
(h < -0.1 * L, where h is the signed height above the neighbor centroid plane).
Any concave vertex is reflected through its neighbor centroid plane:

    x[v] = centroid + |h| * n_hat

After reflection, the CG state is reset (energy/gradient/direction recomputed).

This approach replaced an earlier softplus convexity penalty (E_conv) because:

- It is cheaper: O(V) check with small constant vs O(V) gradient with cross
  products and sigmoid evaluations
- It is simpler: no k_conv parameter to tune, no extra optimization phase
- It is more effective: 0 concave out of 913 isomers vs 1 with E_conv
- It converges faster: the energy landscape without E_conv is simpler, so
  the CG optimizer reaches equilateral solutions in fewer iterations

### Why E_conv failed

An earlier approach used a softplus convexity penalty:

    E_conv = k_conv * sigma * log(1 + exp(-h/sigma)),  sigma = 0.2 * L

with exact gradient including the normal-rotation term. This was tested with
k_conv from 1 to 50. The fundamental problem: for certain topologies, degree-5
vertices invert through the neighbor plane within the first 10 iterations.
The equilateral forces overpower any reasonable k_conv value, and higher k_conv
values degrade triangle quality globally. Periodic reflection avoids this
tradeoff entirely by operating outside the energy function.

The E_conv code is retained in deltahedron_energy_and_gradient() (guarded by
if(k_conv > 0)) for gradient_check() validation and potential future use.

## Two-Phase Optimization

Each phase gets up to max_iter/2 iterations.

### Phase 1: Flatness + equilateral

Active terms: E_bond + E_angle + E_curv + E_flat.

Exits early when the gradient norm drops 100x from its value at the start
of phase 1.

### Phase 2: Pure equilateral (k_flat = 0)

Active terms: E_bond + E_angle + E_curv.

Converges when gradient norm < grad_tol (default 1e-12).

## Test Results (max_iter = 12*N)

| Test                  | Isomers | Time | NC  | Concave | Notes                  |
|-----------------------|---------|------|-----|---------|------------------------|
| SmallFullerenes C20-50| 812     | 12s  |  83 | 0       |                        |
| C60 Hard (50 worst)   | 50      |  6s  |  25 | 0       | All perfect geometry   |
| C60 IPR               | 1       |  0s  |   0 | 0       | r5/r6 = 1.169839 exact |
| C100 IPR sample       | 50      |  2s  |   0 | 0       | All convex             |
| **Total**             | **913** |**20s**|    |         |                        |

NC = not converged (gradient norm > 1e-12). All NC cases have perfect
equilateral geometry (L_cv ~ 0, angles = 60 degrees) -- they just need more
iterations for the gradient to reach the tight tolerance.
