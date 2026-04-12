# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fullerenes is a C++17/Fortran application for topological analysis and 3D structure generation of fullerene isomers (carbon cage molecules). It performs graph-theoretic analysis, geometry optimization, and supports GPU acceleration via SYCL (Intel DPC++). The codebase combines a legacy Fortran solver with a modern C++ front-end and GPU-accelerated computational kernels.

## Build Commands

```bash
# Configure (from repo root)
mkdir build && cd build
cmake ..

# Build everything
cmake --build build

# Build with SYCL GPU support (requires Intel DPC++ compiler)
# CMake auto-detects sycl-ls and sets SYCL_TARGETS to NVIDIA, AMD, or X86
cmake -DENABLE_SYCL=ON -DSYCL_TARGETS=NVIDIA ..

# Run all tests (CTest)
cd build && ctest

# Run a single SYCL test
./build/tests/sycl-tests/eigen-functor-test
```

## Key Build Options

- `ENABLE_SYCL` — Build GPU code (auto-detected from `sycl-ls`)
- `SYCL_TARGETS` — GPU backend: `NVIDIA`, `AMD`, or `X86` (single target per build)
- `CMAKE_CUDA_ARCHITECTURES` — NVIDIA SM architecture (e.g. 86 for RTX 3090)
- `FORTRAN_NMAX` — Max vertices for Fortran static allocation (default 5000)
- `GPU_MAXREGCOUNT` — GPU register pressure limit (default 80)

## Dependencies

Required: CMake 3.5+, C++17 compiler, Fortran compiler (gfortran), BLAS, LAPACK, GSL, GTest, OpenMP.
Optional (GPU): Intel DPC++ compiler (`icpx`/`clang++`), CUDA toolkit, oneDPL.

## Architecture

### Class Hierarchy (Graph Domain Model)

The core domain model is an inheritance chain in `include/fullerenes/`:

```
Graph                    — Undirected/directed graph with adjacency lists (graph.hh)
  └─ PlanarGraph         — Adds planar embedding + 2D Tutte layout (planargraph.hh)
       └─ CubicGraph     — Specialized for 3-valent (degree-3) graphs (cubicgraph.hh)
            └─ FullereneGraph — Fullerene-specific: spirals, Goldberg-Coxeter, leapfrog (fullerenegraph.hh)

Polyhedron : PlanarGraph — 3D coordinates + volume/area/inertia + force-field optimization (polyhedron.hh)
Triangulation            — Surface triangulation from dual spiral (triangulation.hh)
```

### Key Data Types (geometry.hh, auxiliary.hh)

```cpp
typedef int node_t;
typedef vector<vector<node_t>> neighbours_t;  // Sparse adjacency
typedef pair<int,int> edge_t;                  // Undirected edge {min,max}
typedef pair<int,int> arc_t;                   // Directed edge
typedef vector<int> face_t;                    // Face as node cycle
struct coord2d : public pair<double,double>;   // 2D with operator overloads
struct coord3d;                                // 3D with cross product, polar angles, etc.
```

### Spiral Nomenclature (spiral.hh)

Fullerene topology is compactly encoded as spiral indices and jump sequences. The `general_spiral` and `spiral_nomenclature` types are first-class representations used throughout for isomer identification and generation.

### GPU/SYCL Subsystem

Two header directories support GPU computation:

- `include/fullerenes/sycl-headers/` — Data structures for GPU: `Fullerene<T,K>` template (float/double, uint16_t), span-based views, device queue management, status enums for pipeline stages
- `include/fullerenes/kernel-headers/` — Kernel functors inheriting from `base-functor.hh`: eigen (QR decomposition), forcefield-optimize (molecular dynamics), dualize, tutte, hessian, spherical-projection, geometry-properties
- `src/sycl/` — Implementations compiled with DPC++, each with tuned `--maxrregcount` flags

The GPU pipeline processes batches of fullerene isomers through stages: graph dualization → Tutte layout → spherical projection → eigenvalue decomposition → force-field geometry optimization.

### Library Structure

- `fullerenes` (shared lib, `src/c++/`) — Core C++ graph/geometry/spiral library
- `sycl_fullerene_lib` + component libs (`src/sycl/`) — GPU kernels, each compiled separately with per-module register constraints
- `fortran_opt`, `fortran_lib` (`src/fortran/`) — Legacy Fortran solvers (static allocation, configured via `config.f.in`)
- `programs/` — 20+ CLI tools for conversion, generation, and analysis
- `src/contrib/mgmres.cc` — External GMRES solver contribution

### Fortran Integration

Fortran routines handle ring spiral indexing, coordinate generation, force-field optimization, and Schlegel diagrams. They use static allocation bounded by `FORTRAN_NMAX`. The C++/Fortran interface is in `fortran.hh` and `graph_fortran.cc`.

## Branches

- `master` — main/stable branch
- `development` — active development branch

## Database

`database/` contains fullerene isomer data (ring spiral pentagon indices). The database is complete up to C130. Subdirectories: `All/` (general), `IPR/` (isolated pentagon rule), `Yoshida/`, `HOG/`. Path configured via `FULLERENE_DATABASE_PATH` CMake variable.

**IMPORTANT: Always use the IsomerDB C++ interface** (`include/fullerenes/isomerdb.hh`, `src/c++/isomerdb.cc`) to access isomer data. Never try to read the raw database files directly or guess their format.

Key API:
```cpp
// Get number of isomers (hard-coded in Nisomers_data arrays, no file I/O needed)
int64_t n = IsomerDB::number_isomers(60, "Any", false);  // 1812 general C60 isomers
int64_t n = IsomerDB::number_isomers(60, "Any", true);   // 1 IPR C60 isomer

// Read all isomers for a given N
IsomerDB db = IsomerDB::readPDB(60, false);  // all 1812 C60 isomers
IsomerDB db = IsomerDB::readPDB(80, true);   // 7 IPR C80 isomers

// Get a single isomer and build its graph
IsomerDB::Entry e = IsomerDB::getIsomer(60, idx, false);
FullereneGraph G = IsomerDB::makeIsomer(60, e);
```

Isomer counts (general): C20=1, C24=1, C26=1, C28=2, C30=3, C32=6, C34=6, C36=15, C38=17, C40=40, C42=45, C44=89, C46=116, C48=199, C50=271, C60=1812, C70=8149, C80=31924, C100=285914.
Isomer counts (IPR): C60=1, C70=1, C72=1, C74=1, C76=2, C78=5, C80=7, C82=9, C84=24, C86=19, C88=35, C90=46, C92=86, C94=134, C96=187, C98=259, C100=450.
Fullerenes exist for all even N >= 20 except N=22.

## Design Priorities

1. **Correctness first**: Algorithms must *never* fail. They will run unsupervised on 10^11 graphs and must produce the correct result for every one. Arguments like "X is rare" are not acceptable — a single failure among millions means bugs remain.
2. **Cleanness and conciseness**: Code must be easily understood and reasoned about. The goal is formal proofs of correctness.
3. **Efficiency**: Target O(N) or close, with small coefficients, to handle 10^11 graphs.

## Rules

- **No git commits**: Never commit to git yourself. Stage the changed files with git add and propose a commit message. The user always reviews and commits manually.
- **No backticks in commit messages**: Commit messages get pasted into an editor where backticks cause problems. Use plain text instead.
- **Never stage claude-projects/**: The `claude-projects/` directory has its own separate git repo. Never add, stage, or commit files under `claude-projects/` to the fullerene repository. It is listed in `.gitignore`.
- **NEVER kill background processes without EXPLICIT user approval**: Long-running computations (benchmarks, enumerations) may represent hours of irreplaceable work. ALWAYS ask "Can I stop task X? It has been running for Y time and has Z partial results" and WAIT for the user to confirm. Even if you discover a bug and want to relaunch, do NOT stop running processes — modify code, rebuild separately, and ask the user whether to stop the old run. This applies to TaskStop, Ctrl-C, `kill`, and any other means of termination. NO EXCEPTIONS.
- **NEVER do a git reset without EXPLICIT user approval**: This may throw away days of work. Always ask the user.
- **Long-running experiments must support partial results and resumption**: When writing benchmarks or experiments that may run for a long time: (1) Write partial results incrementally (e.g., flush JSON after each isomer, or write checkpoints periodically) so that killing or crashing doesn't lose all progress. (2) Where possible, support restarting from a partially completed state (e.g., save enumeration state, or accept a "start from" parameter). Never design an experiment that only writes output at the very end. (3) Use un-buffered output.
- **NEVER relax guards or validation thresholds**: When a guard condition, assertion, or validation test fails, the solution is NEVER to relax the constraint to make it pass. Guards exist to detect errors — finding errors is essential for correct code. Instead: investigate the root cause, fix the underlying algorithm, not the check. Only change a threshold if there is a clear mathematical proof that the new value is correct. The goal is not a passing test — the goal is correct code that never produces the wrong result.
- **Never prefix Bash commands with `cd dir &&`**: The permissions system matches the entire command string against glob patterns like `Bash(cmake*)`. A command `cd /some/dir && cmake ...` does NOT match `Bash(cmake*)` because it starts with `cd`. Instead, use absolute paths: `cmake --build /path/to/build` or `cmake -S /path/to/src -B /path/to/build`.

## Invariants

- **Orientation**: All graphs and triangulations must ALWAYS be oriented. This is not optional — orientation must be maintained at all times when inserting or removing edges or vertices. Never call `orient_neighbours()` or `orient_triangulation()` — if you find yourself needing them, the graph was constructed wrong. New code must always produce oriented neighbour lists directly and preserve orientation through every topological operation. If you encounter a situation where a graph or triangulation is not oriented, treat it as a bug.
- **`zero_order_geometry()`**: This function is fragile and often fails for larger graphs. It is OK to use for seed graphs (C20/C28/C30) but should be avoided for larger fullerenes.
- **No radial assumptions**: Never assume fullerene geometry is spherical. Fullerenes can be elongated (nanotubes), oblate, or irregular. Never use "distance from origin" or "radial projection" as a geometric primitive — it only works for spheres. Instead, use local geometric information: outward face/fan normals derived from the CW/CCW orientation of neighbor lists, neighbor centroids, and edge vectors. The orientation convention is CCW-on-outside (fan normal from consecutive edge cross products points outward).
- **Iteration budgets are safeguards, not tuning knobs**: Never try to speed up optimization by reducing max iteration counts. The max iterations should be generous enough to never be reached. Instead, improve convergence speed by fixing convergence criteria (tolerances, scaling), improving the energy landscape (better preconditioners, better initial geometry), or reducing per-iteration cost. If the optimizer is hitting the iteration budget, the problem is that convergence is too slow — solve that rather than forcing a cruder result by cutting the budget.
- **Alexandrov's embedding theorem**: By Alexandrov's theorem, every deltahedron (convex polyhedron with all-equilateral triangular faces) has a unique convex isometric embedding (up to rigid motion) where ALL triangles are perfectly equilateral. If the optimizer produces non-equilateral triangles (angle error > 0), that is an optimizer failure — the true energy minimum has zero angle error. Never claim angle error is "intrinsic to the topology."
- **No coordinate-based symmetry inference**: Symmetry (the assignment of O(3) matrices to combinatorial automorphisms) must be derived PURELY from graph combinatorics, NEVER from 3D coordinates. Never use Procrustes, coordinate cross-correlations, intertwiner averaging, or any method that infers the permutation-to-matrix pairing from 3D positions. The Symmetry class provides site symmetries for every vertex, edge, and face — use these to unambiguously assign O(3) operations to automorphisms via their cycle structure (fixed vertices → rotation axes, etc.).

## Deltahedron (Dual Geometry Optimizer)

The `Deltahedron` class (`include/fullerenes/deltahedron.hh`, `src/c++/deltahedron.cc`) computes equilateral-triangle embeddings of fullerene dual graphs. Key entry point:

```cpp
Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);
```

This reduces a fullerene to a seed (C20/C28/C30), uses precomputed seed geometry, and incrementally expands with per-step patch optimization (trust-region Newton) and reflect-optimize loops. See PATCH-GEOMETRY.md for the full pipeline description.

### Optimizer methods

Three methods available via `OptMethod`: CG (Polak-Ribiere), LBFGS (L-BFGS with m=10), STEIHAUG (trust-region Newton-CG with Hessian-vector products). `fromExtensionPathOptimized` takes both `method` (per-step) and `final_method` (final polish) parameters, allowing mixed configs like LBFGS+STEIHAUG.

### Work budget

The optimizer uses a unified work budget: `n_energy + 2*n_grad + 7*n_hv` (calibrated cost ratios). Default: 400*Nv^2 where Nv = dual vertex count = N/2+2. Phase 1 (E_flat) gets 1/4 of budget in standalone optimize(); E_flat is always OFF (k_flat=0) in the extension path pipeline.

### Convexity: reflect-optimize loops

Convexity is maintained by reflect-optimize loops at every level, not by energy penalties. `reflect_all_concave` mirrors concave vertices through their neighbor centroid plane (basin-switching, up to 20 passes), then the optimizer converges in the convex basin. Both the patch optimizer and full-graph optimizer use this loop: optimize freely (no hard constraint) → reflect concave on full graph → re-optimize if anything was reflected. The final phase uses constrained Steihaug (h>=0 trust region) to lock in convexity. E_conv (softplus, k=5) is used in the patch optimizer for smooth guidance but is NOT used in the full-graph optimizer (k_conv=0) because its Hessian dominates and corrupts angle quality.

### Stagnation detection

When angle_tol is set, the optimizer tracks whether energy decreases meaningfully. After 50 iterations without improvement, it breaks early (stuck at local minimum).

### Benchmark and diagnostic tools (benchmarks/)

- **bench_epopt** — Batch quality benchmark. Generates M random isomers of size N via buckygen + reservoir sampling, runs fromExtensionPathOptimized on each, outputs JSON with per-isomer quality stats (edge_cv, edge_relerr_max, h_min, n_concave, ang_min/max/std, ang_relerr_max, gmax_L, iters, ms) and summary statistics. Supports resumption from partial JSON output.
  Usage: `bench_epopt N M output.json`
- **bench_quality_pipeline** — Multi-method pipeline benchmark. Tests different method/tolerance configs across buckygen-enumerated isomers. Supports `--method` (CG/LBFGS/ST), `--final-method`, `--step-tol`, `--final-tol`, `--work-factor`, `--angle-tol`. Outputs CSV with per-isomer metrics.
  Usage: `bench_quality_pipeline N M output.csv [--method LBFGS --final-method ST ...]`
- **bench_report.py** — Formats bench_epopt JSON output as markdown tables with summary across sizes and worst-5 outliers.
  Usage: `python3 bench_report.py bench_C60.json bench_C70.json ...`
- **patch_diag** — Per-expansion-step diagnostic (GTest). Replays the expansion pipeline manually, dumping mol2 after each sub-phase (lift, patch optimize, full CG). Prints per-vertex h values, edge lengths, signed triangle volumes. Takes N and buckygen_idx as CLI args.
  Usage: `patch_diag N idx` or `patch_diag nanotube N`
- **step_mol2** — Step-by-step quality evolution. Builds partial ExtensionPaths (first k steps) and shows quality stats at each intermediate size.
- **difficult_isomers.json** — 78 difficult isomers (C200/C250/C300) that fail to converge under LBFGS+ST with default budget. Contains spiral representations for instant reconstruction without buckygen enumeration.

Build (from build2/): `cmake --build . --target bench_epopt` (or bench_quality_pipeline, patch_diag, step_mol2).


## Recent Development Notes

Linear algebra operations have been externalized to a separate "BatchLAS" repository. CUDA dependency has been removed in favor of SYCL-only GPU support. Matrix view semantics (`DenseMatView`) are replacing handle-based approaches.
