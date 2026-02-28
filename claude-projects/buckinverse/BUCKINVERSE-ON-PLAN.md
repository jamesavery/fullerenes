# O(N) Reduction to Seed: Design and Implementation

## Problem

The reduce-to-seed loop was O(N^2) per isomer: O(N) reduction steps, each
O(N) due to full adjacency copy, vertex renumbering, degree-5 scan, and
F-ring scan. We want O(N) total.

## Data Structure: Bitmask Over Fixed Array

Following RemainingGraph (triangulation.cc:572), each vertex stores a
fixed-size neighbour array with a bitmask selecting active entries:

```cpp
struct Vertex {
    node_t nbr[6];     // CW-ordered neighbours (max degree 6)
    uint8_t active;     // Bitmask: bit i set iff nbr[i] is present
};
```

The full structure maintains:
- `vector<Vertex> V` indexed by original vertex ID (never resized)
- `set<node_t> deg5` — the 12 degree-5 vertices (maintained incrementally)
- `int n_alive` — live vertex count

## Three Mutations

Every reduction (L, B, F) decomposes into sequences of three O(1) operations:

1. **splice(u, old_v, new_v)**: Replace old_v with new_v at the same CW
   position. The bitmask bit stays set; only the array value changes.

2. **snip(u, v)**: Clear the bit for v in u's bitmask. Degree drops by 1.

3. **kill(v)**: Zero v's bitmask, decrement n_alive, remove from deg5.
   No neighbour updates needed — reconnection already replaced all
   references to v before kill is called.

CW navigation (next/prev) uses circular bit scan over the 6-bit bitmask.
degree(u) = popcount(active).

## Site Finding — O(1) Amortised

All site finders scan `deg5` (12 entries, never the full vertex set):

- **L0**: Adjacent degree-5 pair. For each u in deg5, check if any
  active neighbour v is also degree-5.
- **Li**: Degree-5 vertex adjacent to degree-6 chain. Walk via advanceCW(3)
  through degree-6 vertices until hitting another degree-5.
- **B(0,0)**: Triple (deg-5, deg-6, deg-5) with the two deg-5 vertices
  non-adjacent. CW ordering around all three vertices validates the site.
- **F-ring**: Nanotube pole detection — a deg-5 vertex with all deg-5
  neighbours. O(12 x 5) = O(1). Ring vertices are the deg-6 neighbours
  of pole-neighbours, distinguished from strip vertices by also being
  adjacent to pole-neighbours. 5-cycle traced via advanceCW(3).

`findSite()` returns the first valid site, short-circuiting.

## Validity Filter

C22 does not exist as a fullerene (valid sizes are C20, C24, C26, ...).
Any reduction whose result vertex count equals 13 (C22's dual) or falls
below 12 (C20's dual) is rejected. This prevents structurally valid-looking
sites from producing impossible intermediate graphs.

## Reduction Loop

```cpp
SeedType reduceToSeed(int maxRedLen) {
    while (auto site = findSite(maxRedLen))
        reduce(*site);
    // n_alive is now 12 (C20), 16 (C28-Td), or 17 (C30-D5h)
}
```

O(N) construction + O(N) steps x O(1) per step = O(N) total.

## Bugs Found During Implementation

1. **B(0,0) degree check**: Checked degree of q2 (tp endpoint) instead of
   p4 (path endpoint). The snipped vertices are path[0] and path[4], so
   those are the ones that must have degree 6 before reduction.

2. **B(0,0) DLeft q1 validation**: Used `prev(b, c)` (gives p3) instead
   of `next(b, c)` (gives q1). This made the DLeft direction for B(0,0)
   always fail, missing half of all B(0,0) sites.

3. **F-ring ring_second selection**: Picked "any degree-6 neighbour of
   ring_start" as the second ring vertex, but ring_start has 4 degree-6
   neighbours (2 ring + 2 strip). Fixed by requiring ring_second to be
   adjacent to a pole-neighbour, which distinguishes ring from strip.

4. **Invalid intermediate sizes**: B(0,0) from C28 (16v) produces C22
   (13v), which doesn't exist. The reduction surgery produced degree-4
   vertices. Fixed by filtering reductions whose result size is invalid.

## Empirical Results

Verified on all 1,456,590 isomers C32-C100, zero failures.

| N    | Dual verts | Isomers  | red us/isomer | us/vertex |
|------|-----------|----------|---------------|-----------|
| C60  | 32        | 1,812    | 14.5          | 0.45      |
| C80  | 42        | 31,924   | 20.4          | 0.49      |
| C100 | 52        | 285,914  | 25.4          | 0.49      |

Per-vertex time is constant at ~0.49 us, confirming O(N) per isomer.
Buckygen generation is ~6 us/isomer (amortized constant, independent of N).
Speedup vs O(N^2) reference: 53x at C80.

Seed distribution: C20 = 91.2%, C28 = 3.4%, C30 = 5.4%.
All C28 and C30 results verified as genuine seeds (no further reductions).

## Files

- `buckinverse.hh` — ReducibleDual declaration (struct, inline queries)
- `buckinverse.cc` — Implementation (~570 lines added, old code retained)
- `test_scale.cc` — Comparison test (reference vs fast, with seed validation)
- `test_fast.cc` — Fast-path-only test for scaling beyond C80
