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

`database/` contains fullerene isomer data (ring spiral pentagon indices) in subdirectories: `All/` (general, up to C150), `IPR/` (isolated pentagon rule, up to C200), `Yoshida/`, `HOG/`. Path configured via `FULLERENE_DATABASE_PATH` CMake variable.

## Rules

- **No git commits**: Never commit to git yourself. Instead, show the proposed commit message and let the user commit manually.
- **No backticks in commit messages**: Commit messages get pasted into an editor where backticks cause problems. Use plain text instead.
- **Never stage claude-projects/**: The `claude-projects/` directory has its own separate git repo. Never add, stage, or commit files under `claude-projects/` to the fullerene repository. It is listed in `.gitignore`.

## Invariants

- **Orientation**: All graphs must be constructed oriented from the start (`is_oriented == true`). Never call `orient_neighbours()` — if you find yourself needing it, the graph was constructed wrong. New code should always produce oriented neighbour lists directly.
- **`zero_order_geometry()`**: Do not use this function in new code. It is fragile and often fails.

## Recent Development Notes

Linear algebra operations have been externalized to a separate "BatchLAS" repository. CUDA dependency has been removed in favor of SYCL-only GPU support. Matrix view semantics (`DenseMatView`) are replacing handle-based approaches.
