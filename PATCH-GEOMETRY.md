# Full-Graph Optimizer: CG vs L-BFGS vs Steihaug-Toint

## Background

The Deltahedron full-graph optimizer minimizes a combined energy
(E_bond + E_angle + E_curv + E_flat + E_conv) to produce equilateral
triangulations. The original implementation uses Polak-Ribiere nonlinear
CG. Benchmarks on C60-C100 (41,917 isomers) showed ~1-2% converge slowly,
ending with ang_relerr ~0.5% after the 12*N iteration budget. These are
NOT local minima -- extra CG iterations fix them all. They are
ill-conditioned problems where CG's linear convergence is too slow.

Two additional methods were implemented:
- **L-BFGS**: quasi-Newton with limited-memory inverse Hessian (m=10 history).
  Superlinear convergence, O(mN) per iteration.
- **Steihaug-Toint**: trust-region Newton-CG using matrix-free Hessian-vector
  products. Quadratic convergence near the solution.

All three methods share the same energy/gradient function, two-phase
structure (phase 1 with E_flat, phase 2 without), periodic convexity
reflection, and convergence criterion (max_i(||g_i|| * L) < grad_tol).

## Implementation

Files modified:
- `include/fullerenes/deltahedron.hh` -- OptMethod enum {CG, LBFGS, STEIHAUG},
  opt_method member, evaluation counters (n_energy_evals, n_grad_evals, n_hv_evals)
- `src/c++/deltahedron.cc` -- Hessian-vector product, L-BFGS two-loop recursion,
  Steihaug-Toint trust-region solver, refactored optimize() dispatch

### Hessian-vector product

`deltahedron_hv_product()` computes H*v without assembling H:
- **E_bond**: per-edge, O(|E|). M*dv = k*(1-L/r)*dv + k*(L/r^3)*(d.dv)*d
- **E_angle**: per-triangle-corner, exact arm-space blocks from the same
  formulas as assemble_patch_hessian
- **E_curv**: per-vertex rank-1 term (dK x dK) + curvature correction
  (dev * d^2K/dx^2)
- **E_flat**: finite-difference fallback (phase 1 only)

FD validation: gradient_check ~1e-9, cross-method solution difference < 1e-10 * L.

### L-BFGS

Two-loop recursion with m=10 history depth. Armijo backtracking line search
(same as CG). History cleared on reflection or phase transition. Curvature
safeguard: skip update if y.s < 1e-10 * s.s.

### Steihaug-Toint

Outer trust-region loop with inner CG solving H*z = -g. Handles negative
curvature by stepping to trust-region boundary. Accept/reject based on
rho = actual_reduction / predicted_reduction. Trust region expands on
rho > 0.75, shrinks by 4x on rejection.

## Benchmark: C60 (1812 isomers)

Budget = 24*Nv, tol = 1e-10, Nv = 32.

```
CG      : 1785/1812 converged (98.51%)
  ang_re  : med=2.87e-12  p95=6.50e-12  max=1.29e-04
  edge_re : med=2.44e-12  p95=6.43e-12  max=7.15e-05
  ang_std : med=2.24e-10  p95=6.04e-10  max=3.08e-02  (degrees)
  K_re    : med=1.91e-11  p95=5.63e-11  max=1.77e-03  (Gaussian curv relerr)
  K_std   : med=1.20e-11  p95=4.54e-11  max=1.71e-03  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0832  p05=0.0220  min=0.0001
  concave : med=0  max=0
  iters   : med=105  p95=356  max=768
  work    : med=3620  p95=12225  max=26884
  gmax*L  : med=7.94e-11  p95=9.90e-11  max=2.88e-06
  time    : total=8.2s  mean=4.5ms

L-BFGS  : 1808/1812 converged (99.78%)
  ang_re  : med=3.77e-12  p95=6.31e-12  max=1.29e-04
  edge_re : med=3.12e-12  p95=5.56e-12  max=7.15e-05
  ang_std : med=2.83e-10  p95=4.91e-10  max=3.08e-02  (degrees)
  K_re    : med=2.20e-11  p95=4.37e-11  max=1.77e-03  (Gaussian curv relerr)
  K_std   : med=1.26e-11  p95=2.92e-11  max=1.71e-03  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0832  p05=0.0220  min=0.0000
  concave : med=0  max=0
  iters   : med=63  p95=149  max=768
  work    : med=2084  p95=4924  max=25403
  gmax*L  : med=8.08e-11  p95=9.85e-11  max=1.94e-07
  time    : total=3.4s  mean=1.9ms

Steihaug: 1812/1812 converged (100.00%)
  ang_re  : med=1.49e-13  p95=3.64e-13  max=1.29e-04
  edge_re : med=1.24e-13  p95=3.29e-13  max=7.15e-05
  ang_std : med=1.12e-11  p95=2.92e-11  max=3.08e-02  (degrees)
  K_re    : med=7.73e-13  p95=2.55e-12  max=1.77e-03  (Gaussian curv relerr)
  K_std   : med=4.33e-13  p95=1.79e-12  max=1.71e-03  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0832  p05=0.0220  min=0.0000
  concave : med=0  max=0
  iters   : med=4  p95=4  max=24
  work    : med=2371  p95=3971  max=38967
  gmax*L  : med=3.40e-12  p95=7.30e-12  max=9.43e-11
  time    : total=9.0s  mean=5.0ms
```

## Benchmark: C70 (8149 isomers)

Budget = 24*Nv, tol = 1e-10, Nv = 37.

```
CG      : 7993/8149 converged (98.09%)
  ang_re  : med=3.09e-12  p95=7.43e-12  max=5.43e-05
  edge_re : med=2.68e-12  p95=7.40e-12  max=3.08e-05
  ang_std : med=2.47e-10  p95=7.38e-10  max=1.43e-02  (degrees)
  K_re    : med=2.46e-11  p95=7.16e-11  max=6.88e-04  (Gaussian curv relerr)
  K_std   : med=1.41e-11  p95=5.57e-11  max=7.17e-04  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0651  p05=0.0157  min=0.0002
  concave : med=0  max=0
  iters   : med=131  p95=521  max=888
  work    : med=5157  p95=20501  max=35414
  gmax*L  : med=8.08e-11  p95=9.90e-11  max=3.56e-06
  time    : total=57.5s  mean=7.1ms

L-BFGS  : 8137/8149 converged (99.85%)
  ang_re  : med=4.13e-12  p95=7.07e-12  max=5.43e-05
  edge_re : med=3.50e-12  p95=6.74e-12  max=3.08e-05
  ang_std : med=3.13e-10  p95=5.68e-10  max=1.43e-02  (degrees)
  K_re    : med=2.99e-11  p95=6.05e-11  max=6.88e-04  (Gaussian curv relerr)
  K_std   : med=1.51e-11  p95=3.68e-11  max=7.17e-04  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0651  p05=0.0157  min=0.0000
  concave : med=0  max=0
  iters   : med=82  p95=213  max=888
  work    : med=3120  p95=8101  max=33822
  gmax*L  : med=8.28e-11  p95=9.86e-11  max=2.62e-07
  time    : total=24.1s  mean=3.0ms

Steihaug: 8149/8149 converged (100.00%)
  ang_re  : med=1.73e-13  p95=5.55e-13  max=5.43e-05
  edge_re : med=1.47e-13  p95=5.10e-13  max=3.08e-05
  ang_std : med=1.32e-11  p95=4.63e-11  max=1.43e-02  (degrees)
  K_re    : med=1.13e-12  p95=4.26e-12  max=6.88e-04  (Gaussian curv relerr)
  K_std   : med=5.58e-13  p95=2.63e-12  max=7.17e-04  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0651  p05=0.0157  min=0.0000
  concave : med=0  max=0
  iters   : med=4  p95=4  max=21
  work    : med=3296  p95=6108  max=43828
  gmax*L  : med=3.68e-12  p95=9.79e-12  max=9.98e-11
  time    : total=59.1s  mean=7.3ms
```

## Benchmark: C80 (10642 isomers, stride=3)

Budget = 24*Nv, tol = 1e-10, Nv = 42.

```
CG      : 10447/10642 converged (98.17%)
  ang_re  : med=3.23e-12  p95=8.50e-12  max=6.99e-05
  edge_re : med=2.74e-12  p95=7.70e-12  max=4.23e-05
  ang_std : med=2.58e-10  p95=8.44e-10  max=1.84e-02  (degrees)
  K_re    : med=2.99e-11  p95=9.19e-11  max=9.33e-04  (Gaussian curv relerr)
  K_std   : med=1.48e-11  p95=6.41e-11  max=7.96e-04  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0577  p05=0.0138  min=0.0002
  concave : med=0  max=0
  iters   : med=150  p95=586  max=1008
  work    : med=6619  p95=25907  max=45112
  gmax*L  : med=8.22e-11  p95=9.90e-11  max=1.09e-05
  time    : total=94.2s  mean=8.9ms

L-BFGS  : 10629/10642 converged (99.88%)
  ang_re  : med=4.35e-12  p95=7.55e-12  max=6.99e-05
  edge_re : med=3.63e-12  p95=6.68e-12  max=4.23e-05
  ang_std : med=3.33e-10  p95=6.31e-10  max=1.84e-02  (degrees)
  K_re    : med=3.83e-11  p95=7.66e-11  max=9.33e-04  (Gaussian curv relerr)
  K_std   : med=1.72e-11  p95=4.31e-11  max=7.96e-04  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0577  p05=0.0138  min=0.0000
  concave : med=0  max=0
  iters   : med=104  p95=280  max=1008
  work    : med=4476  p95=12054  max=43427
  gmax*L  : med=8.43e-11  p95=9.87e-11  max=9.58e-06
  time    : total=44.8s  mean=4.2ms

Steihaug: 10642/10642 converged (100.00%)
  ang_re  : med=1.97e-13  p95=8.38e-13  max=6.99e-05
  edge_re : med=1.65e-13  p95=7.21e-13  max=4.23e-05
  ang_std : med=1.50e-11  p95=7.01e-11  max=1.84e-02  (degrees)
  K_re    : med=1.57e-12  p95=7.17e-12  max=9.33e-04  (Gaussian curv relerr)
  K_std   : med=6.93e-13  p95=3.94e-12  max=7.96e-04  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0577  p05=0.0138  min=0.0000
  concave : med=0  max=0
  iters   : med=4  p95=4  max=35
  work    : med=4581  p95=8907  max=176224
  gmax*L  : med=4.13e-12  p95=1.49e-11  max=9.98e-11
  time    : total=107.5s  mean=10.1ms
```

## Benchmark: C100 (9860 isomers, stride=29)

Budget = 24*Nv, tol = 1e-10, Nv = 52.

```
CG      : 9669/9860 converged (98.06%)
  ang_re  : med=3.39e-12  p95=9.16e-12  max=6.24e-05
  edge_re : med=2.85e-12  p95=8.09e-12  max=3.56e-05
  ang_std : med=2.72e-10  p95=9.62e-10  max=1.81e-02  (degrees)
  K_re    : med=4.03e-11  p95=1.18e-10  max=9.80e-04  (Gaussian curv relerr)
  K_std   : med=1.61e-11  p95=7.30e-11  max=8.01e-04  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0409  p05=0.0094  min=0.0002
  concave : med=0  max=0
  iters   : med=191  p95=763  max=1248
  work    : med=10327  p95=41256  max=68295
  gmax*L  : med=8.55e-11  p95=9.93e-11  max=4.04e-05
  time    : total=146.0s  mean=14.8ms

L-BFGS  : 9822/9860 converged (99.61%)
  ang_re  : med=4.94e-12  p95=9.29e-12  max=6.24e-05
  edge_re : med=4.13e-12  p95=8.06e-12  max=3.56e-05
  ang_std : med=3.85e-10  p95=8.42e-10  max=1.81e-02  (degrees)
  K_re    : med=5.82e-11  p95=1.19e-10  max=9.80e-04  (Gaussian curv relerr)
  K_std   : med=2.19e-11  p95=6.05e-11  max=8.01e-04  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0409  p05=0.0094  min=0.0002
  concave : med=0  max=0
  iters   : med=157  p95=461  max=1248
  work    : med=8330  p95=24446  max=66247
  gmax*L  : med=8.58e-11  p95=9.89e-11  max=2.50e-07
  time    : total=86.6s  mean=8.8ms

Steihaug: 9860/9860 converged (100.00%)
  ang_re  : med=2.51e-13  p95=1.55e-12  max=6.24e-05
  edge_re : med=2.10e-13  p95=1.31e-12  max=3.56e-05
  ang_std : med=1.94e-11  p95=1.27e-10  max=1.81e-02  (degrees)
  K_re    : med=2.73e-12  p95=1.72e-11  max=9.80e-04  (Gaussian curv relerr)
  K_std   : med=1.00e-12  p95=7.17e-12  max=8.01e-04  (K deviation std)
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0409  p05=0.0094  min=0.0001
  concave : med=0  max=0
  iters   : med=4  p95=5  max=43
  work    : med=8115  p95=18360  max=324626
  gmax*L  : med=5.04e-12  p95=2.96e-11  max=9.95e-11
  time    : total=203.0s  mean=20.6ms
```

## Scaling test (first isomer per size)

Budget = 24*Nv, tol = 1e-10. Per-element work = n_energy + Nv * n_grad + Nv * n_hv.

```
N      Nv   | CG_it  CG_work  CG_gmax    | LB_it  LB_work  LB_gmax    | ST_it  ST_work  ST_gmax
------+------+-----------------------------+-----------------------------+----------------------------
C60    32   |    92     3165   4.90e-11   |    59     1952   5.27e-11   |     4     2147   2.42e-12
C80    42   |    50     2202   6.33e-11   |    45     1936   7.86e-11   |     4     2817   1.45e-12
C100   52   |   260    14046   8.35e-11   |   154     8169   7.82e-11   |     4     8479   6.84e-12
C120   62   |   194    12414   9.93e-11   |   187    11789   7.91e-11   |     4    11845   5.80e-12
C140   72   |   452    33446   8.63e-11   |   488    35645   9.86e-11   |     4    27795   6.18e-12
C160   82   |  1968   174753   8.05e-02 N |   938    78043   9.45e-11   |    19   106782   8.13e-12
C200  102   |  2304   239614   8.16e-11   |   820    84485   7.98e-11   |    15   220742   8.63e-12
```

CG fails to converge on C160 first isomer even with 24*Nv (gmax*L = 0.08).
L-BFGS and Steihaug both converge on all sizes.

## Metric definitions

- **ang_re**: angle relative error = sum|theta_i - 60 deg| / sum(60 deg), over all
  fan angles at all vertices. Zero means all triangles are perfectly equilateral.
- **edge_re**: edge relative error = sum|L_e - L_mean| / (n_edges * L_mean).
  Zero means all edges have identical length.
- **ang_std**: standard deviation of angles in degrees.
- **K_re**: Gaussian curvature relative error = sum|K_i - K_target_i| / sum|K_target_i|
  where K_i = 2*pi - angle_sum_i (angle deficit) and K_target = pi/3 for deg-5,
  0 for deg-6. Zero means curvature is distributed exactly as prescribed by topology.
- **K_std**: standard deviation of K_i - K_target_i. Measures curvature uniformity.
- **K_total**: total Gaussian curvature = sum K_i. Must equal 4*pi (Gauss-Bonnet).
- **h_min**: minimum signed height of any vertex above its neighbor centroid plane.
  Positive = all convex. Negative = concave vertices exist.
- **concave**: number of vertices with h < 0.
- **work**: per-element operations = n_energy + Nv * n_grad + Nv * n_hv. Counts each
  scalar component of a gradient or Hv evaluation as one unit, making the cost of
  different operation types comparable. One CG iteration costs ~1 gradient + ~2 energy
  evaluations. One Steihaug outer iteration costs ~1 gradient + ~1 energy + k Hv
  products where k is the inner CG iteration count.
- **mean time**: wall-clock time per isomer for the final optimize() call
  (fromExtensionPathOptimized + optimize), not including buckygen enumeration
  or graph reduction. Measured on a single core.

## Analysis

1. **All three methods converge to the same geometry.** The max quality metrics
   (ang_re, edge_re, K_re) are identical across methods for each size. The
   worst-case isomers are the same ~1-2% outlier population regardless of
   optimizer. Cross-validation confirms max point difference < 1e-10 * L.

2. **Steihaug achieves 100% convergence** on C60 and C70 (9961 isomers total)
   within 24*Nv budget, while CG fails on 1.5-2% and L-BFGS on 0.15-0.2%.

3. **Steihaug produces ~20x tighter median quality** than CG/L-BFGS. Median
   ang_re ~1.5e-13 vs ~3e-12, median K_re ~1e-12 vs ~2.5e-11. This is because
   Steihaug converges to gmax*L ~3e-12 while CG/L-BFGS stop at ~8e-11 (just
   under the 1e-10 tolerance).

4. **K_total = 4*pi for every isomer**, confirming Gauss-Bonnet is satisfied
   exactly (to machine precision) regardless of optimizer or quality level.

5. **L-BFGS is fastest in wall-clock** (3.0ms mean on C70 vs 7.1-7.3ms for
   CG and Steihaug). Its per-iteration cost is similar to CG but it needs
   fewer iterations. Steihaug's Hv products make each outer iteration ~10x
   more expensive, which offsets the dramatic iteration count reduction.

6. **Steihaug's inner CG cost grows with Nv.** At C200 (Nv=102), the inner
   CG averages ~143 Hv products per outer iteration, approaching the
   max_inner=200 cap. This makes Steihaug's total work comparable to or
   exceeding CG for larger sizes, despite far fewer outer iterations.

7. **The outlier population** (ang_re ~5e-5, K_re ~7e-4 for C70) appears
   identical across all three methods. These isomers converge to a geometry
   with small but nonzero angle error, suggesting the initial geometry from
   the extension path pipeline puts them in a basin where convergence is
   slow. All three methods reach the same final quality given enough budget.

## Benchmark with 48*Nv budget (30,463 isomers)

With a more generous iteration budget (48*Nv), CG and L-BFGS convergence rates
improve but still fall short of Steihaug's 100%.

### Convergence summary

```
         Isomers       CG          L-BFGS       Steihaug
C60        1,812   99.3% (1800)   99.9% (1810)   100% (1812)
C70        8,149   99.5% (8107)  100.0% (8148)   100% (8149)
C80       10,642   99.5% (10591) 100.0% (10636)  100% (10642)
C100       9,860   99.6% (9822)   99.9% (9854)   100% (9860)
```

### C60 (1812 isomers, 48*Nv budget)

```
CG      : 1800/1812 converged (99.34%)
  ang_re  : med=2.87e-12  p95=6.43e-12  max=1.29e-04
  edge_re : med=2.44e-12  p95=6.39e-12  max=7.15e-05
  ang_std : med=2.24e-10  p95=6.04e-10  max=3.08e-02  (degrees)
  K_re    : med=1.91e-11  p95=5.61e-11  max=1.77e-03
  K_std   : med=1.20e-11  p95=4.49e-11  max=1.71e-03
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0832  p05=0.0220  min=0.0000
  concave : med=0  max=0
  iters   : med=105  p95=356  max=1536
  work    : med=3620  p95=12225  max=53284
  gmax*L  : med=7.93e-11  p95=9.88e-11  max=6.81e-07
  time    : total=9.2s  mean=5.1ms

L-BFGS  : 1810/1812 converged (99.89%)
  ang_re  : med=3.77e-12  p95=6.31e-12  max=1.29e-04
  edge_re : med=3.12e-12  p95=5.56e-12  max=7.15e-05
  ang_std : med=2.83e-10  p95=4.91e-10  max=3.08e-02  (degrees)
  K_re    : med=2.20e-11  p95=4.37e-11  max=1.77e-03
  K_std   : med=1.26e-11  p95=2.92e-11  max=1.71e-03
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0832  p05=0.0220  min=0.0000
  concave : med=0  max=0
  iters   : med=63  p95=149  max=1536
  work    : med=2084  p95=4924  max=50772
  gmax*L  : med=8.08e-11  p95=9.85e-11  max=2.16e-07
  time    : total=3.7s  mean=2.0ms

Steihaug: 1812/1812 converged (100.00%)
  ang_re  : med=1.49e-13  p95=3.64e-13  max=1.29e-04
  edge_re : med=1.24e-13  p95=3.29e-13  max=7.15e-05
  ang_std : med=1.12e-11  p95=2.92e-11  max=3.08e-02  (degrees)
  K_re    : med=7.73e-13  p95=2.55e-12  max=1.77e-03
  K_std   : med=4.33e-13  p95=1.79e-12  max=1.71e-03
  K_total : med=12.566371  min=12.566371  max=12.566371  (target=12.566371)
  h_min   : med=0.0832  p05=0.0220  min=0.0000
  concave : med=0  max=0
  iters   : med=4  p95=4  max=24
  work    : med=2371  p95=3971  max=38967
  gmax*L  : med=3.40e-12  p95=7.30e-12  max=9.43e-11
  time    : total=9.5s  mean=5.3ms
```

### Histogram analysis (48*Nv budget, all sizes)

Histograms and per-isomer CSV data are in
`claude-projects/buckinverse/benchmark-data/`.
Regenerate with: `python3 claude-projects/buckinverse/benchmark-data/plot_quality.py`

**Angle relative error** (log10 scale): Steihaug's distribution peaks at ~1e-13,
cleanly separated from CG/L-BFGS which peak at ~3e-12 -- a 1.5 decade gap at
every size. CG/L-BFGS distributions overlap almost completely. All three methods
share a common tail extending to ~1e-4 from the same ~dozen outlier isomers.

**Gaussian curvature relative error** (log10 scale): Same pattern as angle error.
Steihaug peaks at ~1e-12, CG/L-BFGS at ~2e-11. The curvature metric tracks angle
error because K = 2*pi - angle_sum, so angle precision directly determines
curvature precision.

**Timing**: L-BFGS has the tightest, leftmost distribution at every size (~2ms at
C60, ~5ms at C100). CG and Steihaug have similar medians but CG has a long right
tail of slow outliers (up to 70ms at C60, 800ms at C100) from the ~0.5% that exhaust
their iteration budget. Steihaug's timing distribution is tighter than CG's.

**Per-element work** (log10 scale): All three methods have broadly overlapping work
distributions. L-BFGS is slightly cheaper (median ~2000 at C60). Steihaug is
comparable (~2400) despite only 4 outer iterations, because each outer iteration
involves inner CG solves with multiple Hv products. CG's distribution has the
widest spread, reflecting the high variance in iteration counts.

**Scaling trends** (summary plot, medians across sizes):
- Angle and curvature error: Steihaug maintains ~1.5 decade advantage at all sizes
- Timing: L-BFGS ~2x faster than CG/Steihaug, all scale roughly linearly with N
- Work: All three scale similarly, converging at larger N as Steihaug's inner CG
  costs grow

**Convergence failures**: Even with 48*Nv, CG fails on 0.4-0.7% and L-BFGS on
0.01-0.1% of isomers. These are the ill-conditioned cases where superlinear or
quadratic convergence is needed. Steihaug handles all of them.

## Reproducing

```bash
cd build2
cmake --build . --target bench_quality
./benchmarks/bench_quality 70                         # all C70, 48*Nv budget
./benchmarks/bench_quality 60 1 48                    # all C60, 48*Nv budget
./benchmarks/bench_quality 100 10 48 /tmp/out.csv     # every 10th C100, with CSV output
```

Selecting optimizer method in code:
```cpp
Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);
D.opt_method = OptMethod::STEIHAUG;  // or CG, LBFGS
D.opt_k_flat = 0;                    // extension path pipeline
D.optimize(D.points, 0, 24*D.N, 1e-10);
```
