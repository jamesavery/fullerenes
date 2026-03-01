# Deltahedron Geometry Optimization

## Overview

The deltahedron optimizer produces equilateral triangulations of fullerene duals
using conjugate gradient minimization with a multi-phase energy function. This
document describes the energy terms, the three-phase optimization strategy, and
the post-optimization vertex reflection fix for concave local minima.

## Energy Function

The total energy is:

    E = E_bond + E_angle + E_curv + E_flat + E_conv

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

### E_flat (k_flat = 2.0, phases 1-2 only)

Eigenvalue-based flatness measure from flatness.tex. For each vertex v of
degree d <= 6, compute the centroids of the d surrounding triangles. The
smallest eigenvalue lambda_0 of the 3x3 covariance matrix X^T X (where X is
the centroid-corrected matrix) measures deviation from coplanarity.

    E_flat = (k_flat/2) sum_v (lambda_0(v) / tr(A(v)))

The gradient uses the quotient rule and propagates through the centroid
dependencies. The gradient contribution to vertex v itself cancels by symmetry
(v appears in all centroids equally).

### E_conv (k_conv = 2.0, phase 1 only)

Smooth convexity bias via softplus. For each vertex v with degree d <= 6,
compute the signed height h above the neighbor centroid plane:

    h = (x[v] - centroid) . n_hat

where n_hat is the normalized fan normal. The penalty is:

    E_conv = k_conv * sigma * log(1 + exp(-h/sigma)),  sigma = 0.2 * L

This is nearly zero for convex vertices (h > 0), linear in |h| for concave
vertices (h < 0), and smooth everywhere (no kink at h = 0).

The gradient is exact, including the normal-rotation term:

    dh/dx[v]   = n_hat
    dh/dx[n_j] = -n_hat/d + r_perp x (e_{j-1} - e_{j+1}) / |N|

where r_perp = (x[v] - centroid) - h * n_hat, and e_k = x[n_k] - x[v].

## Three-Phase Optimization

Each phase gets up to max_iter/3 iterations. Phase transitions are checked
every 50 iterations.

### Phase 1: Convexity + Flatness

Active terms: E_bond + E_angle + E_curv + E_flat + E_conv.

Exits early when all vertices are convex, or when the concave vertex count
stops improving (stalls for 2 consecutive checks = 100 iterations).

### Phase 2: Flatness only (k_conv = 0)

Active terms: E_bond + E_angle + E_curv + E_flat.

Exits early when the gradient norm drops 100x from its value at the start
of phase 2.

### Phase 3: Pure equilateral (k_flat = 0, k_conv = 0)

Active terms: E_bond + E_angle + E_curv.

Converges when gradient norm < grad_tol (default 1e-12).

## Post-Optimization Vertex Reflection

### Problem: degree-5 vertex inversion

For certain fullerene topologies, a degree-5 vertex can invert through the
neighbor centroid plane within the first 10 iterations, reaching h ~ -1.33
(about 53% of edge length). This happens during phase 1 while E_conv is
active -- the equilateral forces overpower the k_conv = 2 penalty. The
optimizer then converges to a true local minimum with perfect equilateral
triangles but one deeply concave vertex.

This was observed in 4 of 812 isomers tested (C44 #2, C44 #5, C50 #15,
C50 #41). In each case:

- Exactly one degree-5 vertex inverts
- The flip happens within the first 10 iterations
- The optimizer converges to a perfect equilateral solution (L_cv = 0,
  all angles = 60 degrees) with h ~ -1.33 at the inverted vertex
- Increasing k_conv does not help: higher values (e.g. 50) cause global
  degradation without fixing the inversions

### Fix: vertex reflection

After optimization converges, check for deeply concave vertices (h < -0.1 * L).
For each, reflect through the neighbor centroid plane:

    x[v] = centroid + |h| * n_hat

Then re-run the full three-phase optimization with the remaining iteration
budget. The equilateral energy landscape is roughly symmetric around h = 0, so
the reflected position lands in the convex basin.

This fixes all 4 degree-5 inversions with zero impact on runtime or quality.

### Remaining case: C38 isomer 1

One mild concave case (h = -0.087, about 3.5% of edge length) persists. This
involves 6 degree-6 vertices in a symmetric configuration, all at exactly the
same concavity depth. Reflecting all 6 simultaneously does not help because
they are coupled -- the re-optimization falls back into the same symmetric
basin. This appears to be a genuine local minimum of the equilateral energy
for this topology.

## Test Results (max_iter = 12*N)

| Test                  | Isomers | Time | NC  | Concave | Notes                    |
|-----------------------|---------|------|-----|---------|--------------------------|
| SmallFullerenes C20-50| 812     | 16s  | 102 | 1       | C38 #1 only (h = -0.087) |
| C60 Hard (50 worst)   | 50      | 6s   | 22  | 0       | All perfect geometry     |
| C60 IPR               | 1       | 0s   | 0   | 0       | r5/r6 = 1.169839 exact   |
| C100 IPR sample       | 50      | 2s   | 0   | 0       | All convex               |
| **Total**             | **913** |**24s**|     |         |                          |

NC = not converged (gradient norm > 1e-12). All NC cases have perfect
equilateral geometry (L_cv ~ 0, angles = 60 degrees) -- they just need more
iterations for the gradient to reach the tight tolerance.

## k_conv Sweep Results

Tested on the 149 hardest C60 isomers (max_iter = 2000, exact gradient):

| k_conv |  NC | Convex | Concave | Stuck |
|--------|-----|--------|---------|-------|
|      1 |  53 |    148 |       1 |     0 |
|    **2**|**49**|**148**|   **1** | **0** |
|      5 |  52 |    146 |       3 |     0 |
|      7 |  52 |    145 |       2 |     2 |
|     10 |  50 |    141 |       6 |     2 |
|     20 |  59 |    138 |       8 |     3 |
|     50 | 149 |     65 |      63 |    21 |

k_conv = 2 is optimal: fewest NC, fewest concave, zero stuck.
