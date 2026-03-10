# Span-Based Views in the Fullerenes Library

The graph and polyhedron classes support two modes of operation:

1. **Owning** -- the object allocates and manages its own memory (the default)
2. **View** -- the object references externally-owned memory with zero copy

This enables GPU batch processing, Python/numpy interop, and memory-mapped
data without changing any algorithms.

## Quick Reference

### Creating owned objects (the usual way)

```cpp
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/polyhedron.hh"

// Graph with N vertices, max degree dmax (allocates internally)
Graph G(60, 3);

// From brace-init
Graph G2{{1,2,3}, {0,2,3}, {0,1,3}, {0,1,2}};

// Fill constructor: N vertices, each row initialized to the given vector
Graph G3(32, std::vector<node_t>(6, -1));  // 32 vertices, dmax=6, all -1

// Fullerene from the database
FullereneGraph fg = FullereneGraph::C20();

// Polyhedron with owned points
Polyhedron P = Polyhedron::C20();
```

### Creating views over external memory

```cpp
// External flat buffers (e.g. from GPU, numpy, or batch allocation)
std::vector<node_t> values(N * dmax);
std::vector<uint8_t> deg(N);

// Graph view: wraps the external buffers, no allocation
Graph view(N, dmax,
           std::span<node_t>(values),
           std::span<uint8_t>(deg));

// Polyhedron view: graph + coordinate view
std::vector<coord3d> points(N);
Polyhedron poly(PlanarGraph(view),
                std::span<coord3d>(points), /*face_max=*/6);
```

### How to tell if an object is a view or owned

```cpp
Graph G = FullereneGraph::C20();
G.owned_values.empty();   // false -> owned

Graph view(N, dmax, span_vals, span_deg);
view.owned_values.empty(); // true -> view

Polyhedron P = Polyhedron::C20();
P.points.owned.empty();    // false -> owned points
```

## Design

Each class follows the **SpanVector pattern**: a single type that can be either
a view or owning, discriminated by whether internal owned storage is empty.

```
DenseGraph<K>              view base: span<K> values, span<uint8_t> deg
  Graph : DenseGraph<int>  + optional owned_values, owned_deg vectors
    PlanarGraph : Graph    no extra data (view-capable through Graph)
      CubicGraph           no extra data
        FullereneGraph     no extra data
      Triangulation        no extra data
      Polyhedron           + SpanVector<coord3d> points (view or owned)
    Deltahedron            + SpanVector<coord3d> points (view or owned)
```

Key properties:

- **`DenseGraph<K>::values`** and **`deg`** are `std::span`, not raw pointers
- **Views propagate through the hierarchy**: copying a view Graph into a
  CubicGraph produces a CubicGraph view (no allocation, if dmax is already 3)
- **Mutations write through**: modifying a view's adjacency or coordinates
  modifies the external buffer
- **`faces()`** and **`triangles()`** are computed on demand from the oriented
  adjacency -- they are not stored

## Batch Slicing (GPU Pattern)

The primary motivation for views is processing batches of graphs from one
contiguous allocation:

```cpp
const int B = 1000;  // batch size
const int N = 60;    // vertices per graph
const int dmax = 3;  // cubic graphs

// One big allocation for the batch
std::vector<node_t>  all_values(B * N * dmax);
std::vector<uint8_t> all_deg(B * N);
std::vector<coord3d> all_points(B * N);

// ... fill batch data (e.g. from GPU kernel output) ...

// Process each graph as a view -- zero allocation per graph
for (int b = 0; b < B; b++) {
    std::span<node_t>  vals(all_values.data() + b * N * dmax, N * dmax);
    std::span<uint8_t> degs(all_deg.data() + b * N, N);
    std::span<coord3d> pts(all_points.data() + b * N, N);

    Graph g(N, dmax, vals, degs);
    Polyhedron poly(PlanarGraph(g), pts, 6);

    // All graph/geometry algorithms work on the view:
    double V = poly.volume();
    auto faces = poly.faces();
    bool connected = g.is_connected();
}
```

## Copy and Move Semantics

| Operation | Owned source | View source |
|-----------|-------------|-------------|
| Copy construct | Deep copy (owned) | Shallow copy (view) |
| Move construct | Transfer ownership | Transfer view |
| Copy assign | Deep copy (owned) | Shallow copy (view) |
| Move assign | Transfer ownership | Transfer view |

To force a deep copy of a view into owned storage:

```cpp
Graph view(N, dmax, span_vals, span_deg);

// Option 1: construct from base_t (always copies data)
Graph owned(static_cast<const Graph::base_t&>(view));

// Option 2: manually copy
Graph owned2(view.N, view.dmax);
for (node_t u = 0; u < view.N; u++)
    owned2.assign_row(u, view[u]);
```

## Lifetime Safety

Views do not manage the lifetime of the memory they reference. The caller
must ensure the external memory outlives all views into it.

```cpp
// SAFE: view and buffer have matching lifetimes
std::vector<node_t> buf(N * dmax);
std::vector<uint8_t> deg(N);
Graph view(N, dmax, std::span(buf), std::span(deg));
// ... use view while buf and deg are alive ...

// UNSAFE: dangling view
Graph make_dangling_view() {
    std::vector<node_t> buf(60);
    std::vector<uint8_t> deg(20);
    return Graph(20, 3, std::span(buf), std::span(deg));
    // buf and deg destroyed here -- returned Graph is dangling!
}
```

## Faces and Triangles

`Polyhedron::faces()` and `Triangulation::triangles()` compute their results
on demand from the oriented adjacency. They are O(N) methods, not stored
members.

For code that accesses faces in a hot loop (e.g. an optimization inner loop),
cache the result locally:

```cpp
auto fs = poly.faces();  // compute once
for (int iter = 0; iter < 10000; iter++) {
    for (const auto& f : fs) {
        // ... use f ...
    }
}
```

## Class-Specific Notes

### CubicGraph

The constructor from Graph validates that all vertices have degree 3 and
restrides to dmax=3 if needed. When the input already has dmax=3 (e.g. a
view into a cubic batch buffer), no restride occurs and the view is preserved.

### FullereneGraph

Validates the fullerene property (cubic, planar, exactly 12 pentagons and
the rest hexagons). View is preserved if the input is already dmax=3.

### Polyhedron

Has two independently view-or-owned components:
- Adjacency (inherited from Graph)
- Coordinates (`SpanVector<coord3d> points`)

Both can independently be views or owned. A Polyhedron constructed with
`std::span<coord3d>` gets a coordinate view; one constructed with
`const vector<coord3d>&` gets owned coordinates. The adjacency follows
whichever Graph was passed in.

### Deltahedron

Same pattern as Polyhedron: inherits Triangulation (adjacency) and adds
`SpanVector<coord3d> points` for coordinates.
