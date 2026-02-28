# Standalone Expansion via Extension Paths

## Problem

Given a fullerene dual graph, reduce it to a seed and record the reduction
path. Then reconstruct the original graph from ONLY the seed enum and the
recorded path -- no original graph, no saved bitmap state.

## Representation

### ExtensionPath

```cpp
struct ExtensionPath {
    SeedType seed;
    int full_N;                           // total vertex count of full graph

    struct SeedVertex {
        node_t id;
        node_t nbr[6];
        uint8_t active;                   // bitmask of active positions
    };
    std::vector<SeedVertex> seed_state;   // full vertex data after reduction

    std::vector<ExtensionStep> steps;     // expansion order (seed -> full)
};
```

### ExtensionStep

```cpp
struct ExtensionStep {
    ExpKind kind;                         // L_i, B(i,j), or F
    Dir dir;                              // DRight or DLeft
    std::vector<node_t> strip;            // vertices to create
    std::vector<node_t> path;             // existing vertices (reconnected)
    std::vector<node_t> tp;               // existing vertices (reconnected)
};
```

Vertex IDs are in the full (original) graph's numbering throughout.

### Why seed_state stores full Vertex data

The seed vertices after full reduction have both active AND inactive
neighbour positions. Inactive positions hold "shadow" strip vertex IDs
left behind by prior splice/snip operations during reduction. These
shadows are needed by unsnip during round-trip expansion. For standalone
expansion (which uses insertAfter instead of unsnip), the inactive
positions are not strictly needed, but storing them keeps the seed state
identical for both expansion modes.

## Two Expansion Modes

### Round-trip: expand(ReductionStep)

Uses saved strip vertex states (full nbr[6] + active bitmask captured
before kill). Restores strip vertices via unkill, then undoes splice/snip
operations in reverse. Produces bit-exact output matching the original.

### Standalone: expand(ExtensionStep)

Creates strip vertices from CW adjacency formulas. Does not use any saved
state beyond what is in the ExtensionPath. Uses insertAfter instead of
unsnip for the snipped endpoints. Produces CW-rotation-equivalent output
(same circular neighbour ordering, potentially different physical starting
position in the neighbour array).

## The insertAfter Primitive

### The unsnip problem

During reduction, reconnectStraight/Bent calls snip(path[0], strip[0])
and snip(tp[n], strip[n-1]), which clear the bit for strip[0]/strip[n-1]
in those vertices' active bitmasks without changing the nbr[] value. The
round-trip expand reverses this with unsnip, which finds the strip vertex
ID at an inactive position and re-sets the bit.

Standalone expansion cannot use unsnip because strip vertices created from
CW formulas have all positions active (degree 5 uses positions 0-4, degree
6 uses positions 0-5). There are no inactive positions holding shadow IDs.

### The solution

insertAfter(u, v, after) inserts vertex v into u's CW ordering immediately
after vertex `after`. It requires u to be degree 5 (one unused position)
and promotes u to degree 6 by rebuilding the nbr[] array contiguously:

```cpp
void insertAfter(node_t u, node_t v, node_t after) {
    assert(degree(u) == 5);
    // Extract 5 active neighbors in CW order
    node_t cw[5]; int ci = 0;
    for (int p = 0; p < D_MAX; p++)
        if (V[u].active & (1u << p))
            cw[ci++] = V[u].nbr[p];
    // Find insertion point
    int ins = index of `after` in cw[];
    // Build new 6-element array: cw[0..ins], v, cw[ins+1..4]
    V[u].nbr = new_cw;
    V[u].active = 0x3f;
}
```

### Insertion points by type and direction

For each expansion type, the two snipped endpoints (path[0] and tp[n] for
L-type, path[0] and path[bentLen+4] for B-type) need insertAfter with the
correct CW insertion point. These were derived from face-adjacency analysis
of the expansion topology:

**L-type:**
- path[0]: insert strip[0] after tp[0] (DRight), after path[1] (DLeft)
- tp[n]: insert strip[n-1] after path[n] (DRight), after tp[n-1] (DLeft)

**B(0,0):**
- path[0]: insert strip[0] after tp[0] (DRight), after path[1] (DLeft)
- path[4]: insert strip[2] after path[3] (DRight), after tp[2] (DLeft)

**F-ring:** No insertAfter needed -- reconnectRing uses only splice,
not snip, so expandRing (which reverses splices) works unchanged.

## Standalone Expand Algorithm

```
expand(ExtensionStep step):
  1. Remove path[0] and tp.back() from deg5 set (L/B only)
  2. Create strip vertices with CW formulas (type/direction dependent)
  3. For snipped far endpoint: insertAfter(tp[n], strip[n-1], ...)
  4. Reverse tp splice operations (same as expandStraight/Bent minus unsnip)
  5. Reverse path splice operations
  6. For snipped near endpoint: insertAfter(path[0], strip[0], ...)
```

For F-ring: create strip vertices with CW formulas (chirality detected from
reduced graph via next(ring[0], outer[0])), then call expandRing directly.

## graphFromExtensionPath

```cpp
Graph graphFromExtensionPath(const ExtensionPath& ep) {
    ReducibleDual rd(ep.full_N);           // all-dead capacity array
    for (auto& sv : ep.seed_state)         // place seed vertices
        rd.V[sv.id] = {sv.nbr, sv.active};
    for (auto& step : ep.steps)            // apply expansion steps
        rd.expand(step);
    return rd.toGraph();
}
```

## CW-Rotation Equivalence

The standalone expansion produces graphs that are CW-rotation-equivalent
to the original, not bit-exact. The difference is that insertAfter rebuilds
neighbour arrays contiguously at positions 0-5, while the original graph
has neighbours at positions determined by the reduction history (with gaps
at previously-inactive positions). toGraph() extracts neighbours by
iterating set bits in order, so different physical layouts produce different
linear starting points within the same CW circular sequence.

The round-trip expansion (using saved state + unsnip) preserves the exact
physical layout and produces bit-exact output.

## Verification

Both expansion modes verified on all 1,456,593 fullerene isomers C20-C100:
- Round-trip (saved state): 0 failures, exact match
- Standalone (extension path only): 0 failures, CW-rotation match
- Seed distribution: C20 = 1,328,794; C28 = 49,458; C30 = 78,341
