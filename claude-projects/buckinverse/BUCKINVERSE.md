# Buckinverse: Inverting Buckygen Extension Operations

## 1. Mathematical Foundations

### 1.1 Fullerene Dual Graphs as Triangulations

A fullerene C_n is a 3-regular planar graph on n vertices where every face is
a pentagon or hexagon. Its planar dual is a triangulation T on n/2 + 2 vertices,
where each vertex has degree 5 (dual to a pentagon) or degree 6 (dual to a
hexagon). By Euler's formula for polyhedra (V - E + F = 2) combined with
the handshaking lemma, every such triangulation has exactly 12 vertices of
degree 5 and (n/2 + 2 - 12) = n/2 - 10 vertices of degree 6.

We work exclusively in the dual representation. The dual graph T is a
triangulation of the sphere: every face is a triangle, every edge borders
exactly two triangles, and the neighbour list of each vertex is cyclically
ordered (clockwise when viewed from outside the sphere). This cyclic ordering
is the *planar embedding*, and we maintain the invariant that all neighbour
lists are CW-ordered throughout every operation.

### 1.2 Cyclic Navigation on Triangulations

Let N(u) denote the CW-ordered cyclic neighbour list of vertex u. For an edge
u-v, the position of v in N(u) determines a local coordinate system. We define:

- **nextCW(u, v)**: The neighbour of u immediately after v in CW order.
  Equivalently, if v = N(u)[i], then nextCW(u, v) = N(u)[(i+1) mod deg(u)].

- **prevCW(u, v)**: The neighbour of u immediately before v in CW order.
  Equivalently, N(u)[(i-1) mod deg(u)].

- **advanceCW(u, v, k)**: The neighbour k positions after v in CW order:
  N(u)[(i+k) mod deg(u)].

These primitives support four derived navigation operations that encode the
geometry of "walking along a strip" in a triangulation:

**straightAhead(dir, u, from)**: Advance through vertex u when arriving from
vertex `from`, continuing "straight" along the strip. The offset from `from`
depends on degree and direction:
- At degree-6 vertices: always advance 3 positions CW (the vertex directly
  "across" from the entry point in a hexagonal neighbourhood).
- At degree-5 vertices: advance 3 positions for DRight, advance 2 positions
  for DLeft. This asymmetry is the geometric consequence of a pentagonal
  neighbourhood having odd degree.

**turnAhead(dir, u, from)**: Advance through u with a "bent" turn instead of
going straight. The offset is:
- DRight: advance 2 positions CW.
- DLeft: advance deg(u) - 2 positions CW (equivalently, 2 positions CCW).

**sideNbr(dir, a, b)**: The "parallel path" vertex --- the neighbour of a on
the specified side of the directed edge a -> b:
- DRight: prevCW(a, b) (the neighbour just before b in CW order).
- DLeft: nextCW(a, b) (the neighbour just after b in CW order).

**flipDir**: DRight <-> DLeft. Swapping the direction mirrors all operations
to the opposite side of the strip.

### 1.3 The Three Seed Graphs

Buckygen generates all fullerene isomers by recursive expansion starting from
three irreducible seed graphs. These are the unique smallest fullerene duals
that cannot be obtained by expanding a smaller fullerene:

| Seed | Fullerene | Dual vertices | Degree-5 | Degree-6 | Symmetry |
|------|-----------|---------------|----------|----------|----------|
| S1   | C20       | 12            | 12       | 0        | Ih       |
| S2   | C28       | 16            | 12       | 4        | Td       |
| S3   | C30       | 17            | 12       | 5        | D5h      |

The three seeds have distinct dual vertex counts (12, 16, 17), which provides
a trivial identification criterion. C20 is the icosahedron (all faces
pentagonal). C28 and C30 have specific ring spiral pentagon index (RSPI)
sequences that uniquely identify them among their respective isomer sets.

Every fullerene C_n with n >= 32 (equivalently, every triangulation with
>= 18 dual vertices) can be reduced to one of these seeds by a sequence
of reduction operations. The seed to which a fullerene reduces is not unique
in general --- different reduction paths may lead to different seeds.

### 1.4 Extension Operations

Buckygen uses three families of extension operations that grow a triangulation
by inserting a strip of new vertices. Each operation is parametrised by a
starting directed edge (u, v) and a direction d in {DLeft, DRight}.

#### 1.4.1 The Strip Insertion Model

All extension types share a common geometric structure. Given a path
P = [p_0, p_1, ..., p_m] through the triangulation and a parallel path
Q = [q_0, q_1, ..., q_m] on one side of P (determined by the direction d),
the extension inserts a strip S = [s_0, s_1, ..., s_{k-1}] of new vertices
between P and Q. The new vertices are wired to both the path vertices and
the parallel vertices, effectively "widening" the strip by one row.

Before extension, each path vertex p_i is adjacent to the corresponding
parallel vertex q_i. After extension, the strip vertex s_j is interposed:
p_i is now adjacent to s_j instead of q_i, and q_i is adjacent to s_j
instead of p_i. The strip vertices inherit the triangular face structure,
maintaining the triangulation property.

The key degree changes are:
- The two strip endpoints s_0 and s_{k-1} become degree-5 (new pentagons).
- The path endpoints p_0 and p_m (or q_0 and q_m for certain types) gain
  one neighbour, going from degree-5 to degree-6. This preserves the
  invariant of exactly 12 degree-5 vertices.
- All interior strip vertices and path vertices maintain degree 6.

#### 1.4.2 L_i Extensions (Straight)

An L_i extension inserts a straight strip of i+2 new vertices.

**Path**: Computed by repeated straightAhead steps from the starting edge:
- path[0] = u, path[1] = v
- path[k] = straightAhead(d, path[k-1], path[k-2]) for k >= 2
- Total path length: i+3 vertices

**Parallel path**: One sideNbr per path vertex (matching Haskell convention):
- par[k] = sideNbr(d, path[k], direction_to_next) for k = 0, ..., i+2
- For k < i+2: direction_to_next = path[k+1]
- For k = i+2: direction_to_next = straightAhead(d, path[k], path[k-1])
- Total parallel length: i+3 vertices

**Strip**: i+2 new vertices [s_0, ..., s_{i+1}]:
- s_0 and s_{i+1} are degree-5 (the new pentagons)
- s_j for 0 < j < i+1 are degree-6

**Degree changes at existing vertices**:
- path[0] gains neighbour s_0 (degree 5 -> 6)
- par[i+2] gains neighbour s_{i+1} (degree 5 -> 6)
- All other existing vertices maintain their degree (edge replacements only)

**Special case L_0**: The shortest straight extension. Path has 3 vertices,
strip has 2 vertices (both degree-5). The path[1] vertex (middle of path)
need not be degree-5.

#### 1.4.3 B_{i,j} Extensions (Bent)

A B_{i,j} extension inserts a bent strip of i+j+3 new vertices, with a
"turn" at position i+1 in the strip. The bend creates an elbow in the strip.

**Path structure**:
- Pre-bend: i+2 straight steps from the starting edge, producing path[0..i+1]
- Turn: turnAhead at path[i+2] = straightAhead from the last pre-bend edge
- Post-bend: j+1 straight steps from the turn vertex, producing path[i+3..i+j+4]
- Total path length: i+j+5 vertices

**Parallel path**: sideNbr at each edge, but *skipping the edge entering
the turn vertex* (the edge from path[i+1] to path[i+2]). This gives:
- Pre-bend: par[k] = sideNbr(d, path[k], path[k+1]) for k = 0, ..., i
  (i+1 entries)
- Post-bend: par[k] = sideNbr(d, path[k], path[k+1]) for k = i+2, ..., i+j+3
  (j+2 entries)
- Total parallel length: i+j+3 vertices

The edge skip at the bend is a critical detail. The turn vertex path[i+2]
has three consecutive path neighbours in its cyclic list (path[i+1],
path[i+2] in the turn direction, and path[i+3]), which creates the elbow
geometry. The strip vertex at the bend position (s_{i+1}) correspondingly
has degree 6 with three path neighbours.

**Strip**: i+j+3 new vertices:
- s_0 and s_{i+j+2} are degree-5
- s_{i+1} (the bend vertex) is degree-6 with three path neighbours
- All other interior vertices are degree-6

**Special case B_{0,0}**: Path has 5 vertices, strip has 3 vertices [a, b, c]
where a and c are degree-5, b is degree-6 (the bend vertex). Parallel path
has 3 entries.

**Degree changes**:
- path[0] gains neighbour s_0 (degree 5 -> 6)
- path[i+j+4] gains neighbour s_{i+j+2} (degree 5 -> 6)

#### 1.4.4 F Extensions (Nanotube Ring)

The F extension inserts a ring of 5 degree-6 vertices between two parallel
5-cycles of vertices, widening a nanotube equator by one row. Unlike L and B
extensions, F does not change any vertex's degree — it adds 5 degree-6 vertices
and the number of degree-5 vertices stays at 12.

**Structure**: Given a "ring" 5-cycle R = [r_0, ..., r_4] of degree-6 vertices
and a set of 5 "outer" (cap) vertices O = [o_0, ..., o_4], the expansion inserts
5 new degree-6 vertices S = [s_0, ..., s_4] between R and O.

**Adjacency of new vertex s_i** (CW order):
```
s_i neighbours: [r_i, s_{(i-1) mod 5}, o_i, o_{(i+1) mod 5}, s_{(i+1) mod 5}, r_{(i+1) mod 5}]
```

**Edge replacements at existing vertices**:
- r_i: replace o_{(i-1) mod 5} → s_{(i-1) mod 5}, replace o_i → s_i
- o_i: replace r_i → s_{(i-1) mod 5}, replace r_{(i+1) mod 5} → s_i

**Occurrence**: F extensions produce the (5,0) nanotube family. Starting from
seed C30 (D5h), successive F extensions produce C40, C50, C60, ... Each
intermediate graph is a (5,0) nanotube. These are the only fullerenes that
require F reduction; all other fullerenes reduce to seeds using only L and B.

### 1.5 Reductions as Inverse Operations

A reduction is the inverse of an extension: given a fullerene dual graph T
and a valid reduction site (kind, u, v, d), the reduction removes the strip
vertices and reconnects path and parallel vertices to produce a smaller
triangulation T' with fewer vertices.

The reduction is characterised by the same triple (kind, u, v, d) as the
corresponding extension. The extension (kind, u, v, d) applied to T'
would produce T. Thus the sequence of reductions from T to a seed S
encodes the sequence of extensions from S to T (read in reverse).

**Vertex accounting**:
- L_i reduction: removes i+2 vertices (T has n vertices, T' has n-i-2)
- B_{i,j} reduction: removes i+j+3 vertices
- F reduction: removes 5 vertices

**Degree restoration**:
- For L and B: The two endpoints that gained a degree during extension
  (path[0] and par[i+2] for straight, or path[0] and path[i+j+4] for bent)
  lose a neighbour, restoring them to degree-5. The two strip endpoints that
  were degree-5 are removed entirely. Net change: the set of degree-5
  vertices changes by removing two (the strip endpoints) and adding back
  two (the path/parallel endpoints), preserving the count of 12.
- For F: All 5 removed vertices are degree-6, so the degree-5 count
  is unchanged. Ring and outer vertices regain their original adjacencies.


## 2. Algorithmic Analysis

### 2.1 Reduction Enumeration

The first algorithmic task is to enumerate all valid reduction sites on a
given triangulation T. A reduction site is a tuple (kind, u, v, d) where
kind in {L_0, L_1, ..., B_{0,0}, B_{1,0}, ...}, u is a degree-5 vertex,
(u, v) is a directed edge, and d in {DLeft, DRight}.

The enumeration proceeds in four phases, controlled by a maximum reduction
length parameter (default 5).

#### 2.1.1 L_0 Reductions

For each degree-5 vertex u and each neighbour v with deg(v) = 5:
- For each direction d in {DLeft, DRight}:
  - Check the **flanking condition**: the vertex 2 positions from v in u's
    CW list (on the side determined by d) must be degree-6, and likewise
    the vertex 2 positions from u in v's CW list.
  - Formally:
    - DRight: deg(advanceCW(u, v, 2)) = 6 and deg(advanceCW(v, u, 2)) = 6
    - DLeft: deg(advanceCW(u, v, deg(u)-2)) = 6 and deg(advanceCW(v, u, deg(v)-2)) = 6

The flanking condition ensures that the surgery produces a valid triangulation.
It rules out reductions that would create degree-4 vertices or other
topological defects.

#### 2.1.2 L_i Reductions (i >= 1)

For each degree-5 vertex u and each neighbour v with deg(v) != 5 (L_0 is
handled separately):
- **Follow the straight path**: Starting from (u, v), repeatedly apply
  straightAhead until either:
  - A degree-5 vertex w is reached at distance dist (success)
  - A cycle is detected (failure)
  - The maximum distance is exceeded (failure)
- If a degree-5 endpoint w is found at distance dist:
  - The reduction type is L_{dist-1}
  - For each direction d, check the flanking condition at *both endpoints*:
    the start (u, v) and the end (w, prevW), using the same formula as L_0.

The direction is irrelevant during the straight-ahead walk at degree-6
interior vertices (straightAhead gives the same result for both DLeft and
DRight at degree-6), so the walk uses an arbitrary direction. Direction
matters only at degree-5 endpoints for the flanking check.

**Cycle detection**: Uses a visited set to detect if the walk revisits a
vertex. This is necessary because on small triangulations, a straight walk
can wrap around and revisit the starting region.

#### 2.1.3 B_{0,0} Reductions

For each degree-5 vertex u and each neighbour v (any degree):
- For each direction d:
  - Compute the B_{0,0} path: [u, v, straightAhead(v,u), turnAhead, straightAhead]
  - Compute the 3-entry parallel path (skipping the turn edge)
  - Check:
    1. Path has no duplicate vertices (non-self-intersecting)
    2. Far endpoint path[4] is degree-5
    3. path[0] != path[4] (not a loop)
    4. Path and parallel path are vertex-disjoint

Note that B_{0,0} does *not* require the flanking condition that straight
reductions require. It has a weaker validity check.

#### 2.1.4 B_{i,j} Reductions (i+j > 0)

For each degree-5 vertex u, each neighbour v, each direction d, and each
pair (i, j) with 1 <= i+j <= maxRedLen - 2:
- Compute the bent path via computeBentPathSafe (using a uniform parallel
  path that does *not* skip the turn edge, for validation purposes)
- Check:
  1. Path is non-self-intersecting
  2. Far endpoint is degree-5
  3. Not a loop
  4. All intermediate vertices (path[1..n-2]) are degree-6
  5. The flanking vertex at the far endpoint is degree-6 and not on the path

The flanking direction for bent reductions is the *opposite* of straight
reductions. Where a straight reduction checks advanceCW(..., 2) for DRight,
a bent reduction checks advanceCW(..., deg-2) for DRight (and vice versa).
This reversal occurs because the bend in the path flips the "colour side"
of the strip.

The intermediate degree-6 condition (check 4) ensures that we find the
*shortest* reduction at each site. If an intermediate vertex were degree-5,
a shorter reduction would exist starting or ending there.

### 2.2 Reduction Surgery: Strip-Based Inversion

Given a fullerene dual graph, the inversion algorithm identifies strip vertices
directly from local graph structure, then removes them and reconnects the
remaining vertices. This approach works on arbitrary graphs, not just those
recently produced by expansion.

#### 2.2.1 The Expansion-Site vs. Strip-Site Distinction

The Phase 2 `allReductions` function enumerates *expansion sites* — locations
where an expansion COULD be applied to grow the graph. These are characterised
by degree-5 starting vertices and flanking conditions. However, expansion
sites do not directly identify which vertices to REMOVE for reduction.

The key insight: the strip vertices (added by expansion) are identifiable by
their local topology:
- **L_0 strip**: Two adjacent degree-5 vertices
- **L_i strip (i >= 1)**: A chain [s_0(deg-5), s_1(deg-6), ..., s_i(deg-6),
  s_{i+1}(deg-5)] where interior vertices are connected by straight-ahead
  (advanceCW 3, which is direction-independent at degree 6)
- **B_{0,0} strip**: Three vertices [a(deg-5), b(deg-6), c(deg-5)] where
  a adj b, b adj c, but a NOT adj c
- **F strip**: A 5-cycle of degree-6 vertices parallel to another 5-cycle
  of degree-6 vertices (the "ring")

Once the strip is identified, path and true-parallel vertices are extracted
from the CW ordering around each strip vertex.

#### 2.2.2 CW Ordering Around Strip Vertices (L-Type)

The expansion inserts strip vertices between path and true-parallel. In the
expanded graph, each strip vertex has a specific CW neighbourhood pattern
that reveals the path and true-parallel vertices.

**Near endpoint s_0 (degree 5, i+2 entries in strip):**
```
DRight CW: [path[0], tp[0], tp[1], s_1, path[1]]
DLeft  CW: [path[0], path[1], s_1, tp[1], tp[0]]
```
Extraction: path[1] = next(s_0, s_1) for DRight, prev(s_0, s_1) for DLeft.

**Interior s_j (degree 6, j = 1..i):**
```
DRight CW: [s_{j-1}, tp[j], tp[j+1], s_{j+1}, path[j+1], path[j]]
DLeft  CW: [s_{j-1}, path[j], path[j+1], s_{j+1}, tp[j+1], tp[j]]
```
The strip chain walk s_{j+1} = advanceCW(s_j, s_{j-1}, 3) is direction-
independent at degree 6. Path and tp are extracted incrementally: tp[j+1]
and path[j+1] are derived from s_j's CW ordering given known values of
tp[j] and path[j] from the previous step.

**Far endpoint s_{i+1} (degree 5):**
```
DRight CW: [s_i, tp[i+1], tp[i+2], path[i+2], path[i+1]]
DLeft  CW: [s_i, path[i+1], path[i+2], tp[i+2], tp[i+1]]
```

**Consistency checks at each step**: The CW ordering implies constraints
(e.g., g.next(s_j, tp[j+1]) == s_{j+1}) that are verified. If any check
fails, the candidate is rejected.

#### 2.2.3 CW Ordering Around Strip Vertices (B_{0,0} Type)

For B_{0,0}, the strip has 3 vertices: a (deg-5 endpoint), b (deg-6 bend),
c (deg-5 endpoint). The path has 5 vertices [p0..p4] and tp has 3 [q0,q1,q2].

```
DRight:  a CW: [p0, q0, q1, b, p1]
         b CW: [a, q1, c, p3, p2, p1]
         c CW: [b, q1, q2, p4, p3]

DLeft:   a CW: [p0, p1, b, q1, q0]
         b CW: [a, p1, p2, p3, c, q1]
         c CW: [b, p3, p4, q2, q1]
```

The bend vertex b is degree-6 with 3 path neighbours (p1, p2, p3), distinguishing
it from interior L-type strip vertices (which have only 2 path neighbours).

#### 2.2.4 F-Ring Detection and CW Ordering

F-ring detection uses straight-ahead tracing to find a 5-cycle of degree-6
vertices (the "ring"), then locates the parallel 5-cycle (the "strip") and
the outer (cap) vertices.

**Step 1: Trace ring.** For each pair of adjacent degree-6 vertices (u, v),
trace 5 steps of straightAhead (advanceCW 3). If the trace returns to u,
we have a ring R = [r_0, ..., r_4].

**Step 2: Find strip.** For each ring edge (r_i, r_{(i+1) mod 5}), find
the common non-ring degree-6 neighbour. This is strip[i]. Verify the 5
strip vertices form a 5-cycle.

**Step 3: Resolve outers.** Each strip[i] has exactly 2 non-ring non-strip
neighbours: {outer[i], outer[(i+1) mod 5]}. The outer array is resolved by
intersection: outer[i] = common element of strip[(i-1) mod 5]'s extras and
strip[i]'s extras.

**Step 4: Validate.** All 15 vertices (5 ring + 5 strip + 5 outer) must be
distinct. All strip vertex neighbours must be in {ring, strip, outer}
(strip topology clean).

#### 2.2.5 Reconnection Patterns

The reconnection reverses the expansion's edge rewiring. For each vertex
type:

**L-type path vertices**:
- path[0]: remove strip[0] (degree 6 -> 5)
- path[k] (1 <= k <= pathlength): replace strip[k-1] -> tp[k-1],
  and if k < pathlength: replace strip[k] -> tp[k]
- tp[pathlength]: remove strip[pathlength-1] (degree 6 -> 5)

**L-type true parallel vertices**:
- tp[k] (0 <= k < pathlength): replace strip[k] -> path[k+1],
  and if k > 0: replace strip[k-1] -> path[k]

**B-type**: Before the bend, follows the straight pattern. At the bend,
the bend strip vertex (degree 6, 3 path neighbours) requires 3 replacements
at the corresponding true parallel vertex. After the bend, straight pattern
with shifted indices.

**F-ring reconnection**:
```
ring[i]:  replace strip[(i-1) mod 5] -> outer[(i-1) mod 5],
          replace strip[i] -> outer[i]
outer[i]: replace strip[(i-1) mod 5] -> ring[i],
          replace strip[i] -> ring[(i+1) mod 5]
```

After reconnection, strip vertices are deleted and remaining vertices
renumbered to 0..n'-1.

#### 2.2.6 Detailed Reconnection for L and B Types

**L_0 reduction** (strip = [a, b], path = [p0, p1, p2], tp = [q0, q1, q2]):

```
p0: remove a                    (deg 6 -> 5)
p1: replace a -> q0, b -> q1    (deg unchanged)
p2: replace b -> q1             (deg unchanged)
q0: replace a -> p1             (deg unchanged)
q1: replace a -> p1, b -> p2    (deg unchanged)
q2: remove b                    (deg 6 -> 5)
```

**L_i reduction** (i >= 1; strip has i+2, path has i+3, tp has i+3):

```
path[0]:    remove strip[0]                           (deg 6 -> 5)
path[k]:    replace strip[k-1] -> tp[k-1],
            replace strip[k] -> tp[k]                 for k = 1..i+1
path[i+2]:  replace strip[i+1] -> tp[i+1]            (deg unchanged)

tp[0]:      replace strip[0] -> path[1]
tp[k]:      replace strip[k-1] -> path[k],
            replace strip[k] -> path[k+1]             for k = 1..i
tp[i+1]:    replace strip[i+1] -> path[i+2]
tp[i+2]:    remove strip[i+1]                         (deg 6 -> 5)
```

**B_{0,0} reduction** (strip = [a, b, c], path = [p0..p4], tp = [q0, q1, q2]):

```
p0: remove a                    (deg 6 -> 5)
p1: replace a -> q0, b -> q1
p2: replace b -> q1
p3: replace b -> q1, c -> q2
p4: remove c                    (deg 6 -> 5)

q0: replace a -> p1
q1: replace a -> p1, b -> p2, c -> p3
q2: replace c -> p3
```

Note: q1 (the true parallel at the bend) gains three replacements because
the bend strip vertex b was adjacent to three path vertices (p1, p2, p3)
and all three of those adjacencies were through q1 in the original graph.

**B_{i,j} reduction** (strip has i+j+3 vertices, generalising B_{0,0}):

The pattern extends naturally. Before the bend, reconnection follows the
straight pattern. At the bend, the strip's bend vertex (which has three
path neighbours) requires three replacements at the corresponding true
parallel vertex. After the bend, reconnection follows the straight pattern
with shifted indices.

### 2.3 Iterative Reduction to a Seed

The complete inversion algorithm uses `allInvSites` (strip-based detection)
rather than `allReductions` (expansion-site detection):

```
function reduceToSeed(T):
    while T.N > 17:
        sites = allInvSites(T)
        for site in sites:
            T' = invertReduction(T, site)
            if T' is valid triangulation:
                T = T'
                break
    return identifySeed(T)
```

**Termination**: Every inversion strictly decreases the vertex count
(by at least 2). The smallest possible fullerene duals have 12, 16, or 17
vertices (the seeds). Since every non-seed fullerene has at least one valid
inversion site (empirically verified for all 84 isomers C32--C40), the
algorithm terminates.

**Choice of reduction**: The current implementation tries inversion sites
in order and uses the first that produces a valid triangulation. Buckygen's
canonical form uses a specific ordering (the 5-tuple cascade) to select
reductions; this is not yet implemented.

**Empirical results (C32--C100, 1,456,593 isomers)**:
- All 1,456,593 isomers reduce to seeds with 0 failures
- Seed distribution: C20 = 1,328,794; C28 = 49,458; C30 = 78,341
- Only nanotube-family isomers require F-ring reduction;
  all others reduce via L and B operations alone

### 2.4 O(N) In-Place Reduction: ReducibleDual

The O(N^2) per-isomer cost of the graph-copying approach (Section 2.5)
is eliminated by an in-place bitmask-over-fixed-array representation.

#### 2.4.1 Data Structure

```cpp
struct ReducibleDual {
    struct Vertex {
        node_t nbr[6];    // CW-ordered neighbours (fixed positions)
        uint8_t active;    // Bitmask: bit i set iff nbr[i] is present
    };
    vector<Vertex> V;      // Indexed by original vertex ID (never resized)
    set<node_t> deg5;      // The 12 degree-5 vertices (maintained)
    int n_alive;
};
```

The key property: vertex IDs are never renumbered. A killed vertex has
active = 0 but retains its position in V[]. Degree is popcount(active).
The array V is allocated once at full graph size and never resized.

#### 2.4.2 O(1) Mutation Primitives

Four mutation primitives, each O(1):

**splice(u, old_v, new_v)**: Replace one active neighbour value with another.
Finds old_v in u's active positions, writes new_v in its place. Self-inverse:
splice(u, new_v, old_v) undoes it. Used for reconnection during both
reduction and expansion.

**snip(u, v)**: Remove v from u's active neighbours by clearing its bit.
The value v remains in nbr[] at the cleared position (a "shadow"). Reduces
degree by 1. Used at path/tp endpoints during reduction.

**kill(v)**: Set V[v].active = 0, decrement n_alive, remove from deg5.
Used to delete strip vertices after reconnection.

**unsnip(u, v)**: Inverse of snip. Finds v in u's nbr[] at an INACTIVE
position and re-sets its bit. Works because snip only clears the bit
without changing the value, and splice only operates on active positions.
Used during round-trip expansion to restore snipped edges.

**unkill(v, saved)**: Restore a killed vertex's full state from a saved
Vertex struct. Increments n_alive, adds to deg5 if degree 5. Used during
round-trip expansion to restore strip vertices.

**insertAfter(u, v, after)**: Insert v into u's CW ordering immediately
after vertex `after`, promoting u from degree 5 to degree 6. Rebuilds the
nbr[] array contiguously. Used during standalone expansion (see Section 2.6)
where unsnip is not available.

#### 2.4.3 Reduction Loop

```
reduceToSeed():
    while n_alive > 17:
        site = findSite()       // O(1) amortised
        reduce(site)            // O(1): reconnect + kill
    return identifySeed()
```

findSite() tries L0 first (scan 12 deg-5 vertices, O(1)), then L_i
(bounded chain walk, O(1)), then B(0,0) (O(1)), then F-ring (O(N) but
rare). reduce() calls reconnectStraight/Bent/Ring (bounded splice count)
then kills strip vertices. Total: O(N) amortised.

#### 2.4.4 O(N) Complexity

Each vertex is killed at most once and each splice/snip/unsnip touches
O(1) vertices. Total across all reduction steps: O(N) kills, O(N) splices.
Site finding is O(1) amortised (F-ring O(N) but occurs at most O(1) times
per isomer in practice). Grand total: O(N) per isomer.

**Empirical performance (C20--C100, 1.46M isomers):**
Reduction of all 1,456,593 isomers completes in 42 seconds total.

### 2.5 Round-Trip Expansion

The reduction path can be recorded and replayed in reverse to reconstruct
the original graph from the seed state. Because ReducibleDual preserves
original vertex IDs throughout (never renumbers), the reconstructed graph
is bit-exact identical to the original.

#### 2.5.1 ReductionStep

```cpp
struct ReductionStep {
    InvSite site;
    vector<pair<node_t, Vertex>> saved;  // strip vertex states before kill
};
```

Each saved entry stores a strip vertex's ID and full Vertex struct (nbr[6]
+ active bitmask) captured just before kill(). This is needed because prior
reduction steps may have modified strip vertices via splice (a vertex that
was path[k] in step i may become strip[j] in step i+3).

#### 2.5.2 expand(ReductionStep)

Reverses reduce() in three stages:
1. Remove path[0] and tp.back() from deg5 (L/B only)
2. Unkill strip vertices from saved state
3. Undo reconnection via expandStraight/Bent/Ring

expandStraight/Bent/Ring reverse every splice/snip from the corresponding
reconnect function in exact reverse order. splice is self-inverse
(swap old/new arguments), and unsnip reverses snip.

#### 2.5.3 Verification

Round-trip expansion verified on all 1,456,593 isomers C20--C100.
Comparison is element-by-element exact (same physical neighbour array
layout). Zero failures.

### 2.6 Standalone Expansion via Extension Paths

The standalone expansion reconstructs a fullerene dual graph from ONLY a
seed graph enum and an ExtensionPath -- no original graph or saved bitmap
state required.

#### 2.6.1 ExtensionPath Representation

```cpp
struct ExtensionPath {
    SeedType seed;
    int full_N;
    vector<SeedVertex> seed_state;   // full vertex data after reduction
    vector<ExtensionStep> steps;     // expansion order (seed -> full)
};

struct ExtensionStep {
    ExpKind kind;
    Dir dir;
    vector<node_t> strip;   // vertices to create
    vector<node_t> path;    // existing vertices (reconnected)
    vector<node_t> tp;      // existing vertices (reconnected)
};
```

The seed_state stores full Vertex data (nbr[6] + active bitmask) for all
alive vertices after full reduction to seed. Steps are in expansion order
(reversed from reduction order). All vertex IDs are in the full graph's
numbering.

#### 2.6.2 The insertAfter Solution

Standalone expansion cannot use unsnip because strip vertices created from
CW formulas have all nbr[] positions active -- there are no inactive
positions holding shadow strip IDs. The solution is insertAfter(u, v, after),
which inserts v into u's CW ordering after vertex `after` by extracting the
5 active neighbours in CW order, inserting v at the correct position, and
writing back as a contiguous 6-element array.

The CW insertion points were derived from face-adjacency analysis:
- L-type DRight: path[0] insert strip[0] after tp[0];
  tp[n] insert strip[n-1] after path[n]
- L-type DLeft: path[0] insert strip[0] after path[1];
  tp[n] insert strip[n-1] after tp[n-1]
- B(0,0) DRight: path[0] insert strip[0] after tp[0];
  path[4] insert strip[2] after path[3]
- B(0,0) DLeft: path[0] insert strip[0] after path[1];
  path[4] insert strip[2] after tp[2]
- F-ring: no insertAfter needed (expandRing uses only splice)

#### 2.6.3 Strip Vertex CW Formulas

Each strip vertex's CW neighbour list is computed from the expansion type,
direction, and the known path/tp/strip vertex IDs. The formulas are:

**L-type s_0 (degree 5):**
- DRight: [path[0], tp[0], tp[1], strip[1], path[1]]
- DLeft: [path[0], path[1], strip[1], tp[1], tp[0]]

**L-type interior s_j (degree 6):**
- DRight: [strip[j-1], tp[j], tp[j+1], strip[j+1], path[j+1], path[j]]
- DLeft: [strip[j-1], path[j], path[j+1], strip[j+1], tp[j+1], tp[j]]

**L-type s_{n-1} (degree 5):**
- DRight: [strip[n-2], tp[n-1], tp[n], path[n], path[n-1]]
- DLeft: [strip[n-2], path[n-1], path[n], tp[n], tp[n-1]]

**B(0,0) a (degree 5):**
- DRight: [path[0], tp[0], tp[1], b, path[1]]
- DLeft: [path[0], path[1], b, tp[1], tp[0]]

**B(0,0) b (degree 6, bend vertex):**
- DRight: [a, tp[1], c, path[3], path[2], path[1]]
- DLeft: [a, path[1], path[2], path[3], c, tp[1]]

**B(0,0) c (degree 5):**
- DRight: [b, tp[1], tp[2], path[4], path[3]]
- DLeft: [b, path[3], path[4], tp[2], tp[1]]

**F-ring s_i (degree 6), chirality A (next(ring[0], outer[0]) == outer[4]):**
- [ring[i], ring[i+1], strip[i+1], outer[i+1], outer[i], strip[i-1]]

**F-ring s_i (degree 6), chirality B (next(ring[0], outer[0]) == ring[1]):**
- [ring[i], strip[i-1], outer[i], outer[i+1], strip[i+1], ring[i+1]]

F-ring chirality is detected from the reduced graph: if next(ring[0],
outer[0]) == outer[4] then ordering A, else ordering B.

#### 2.6.4 CW-Rotation Equivalence

The standalone expansion produces graphs that are CW-rotation-equivalent
to the original: each vertex has the same circular neighbour ordering, but
the linear starting position in the neighbour array may differ. This is
because insertAfter rebuilds arrays contiguously (positions 0-5), while the
original has neighbours at positions determined by the bitmask history.

The round-trip expansion preserves the exact physical layout and produces
bit-exact output.

#### 2.6.5 Verification

Standalone expansion verified on all 1,456,593 isomers C20--C100 using
CW-rotation-aware comparison. Zero failures. Timing: 49.5s total for
reduce + standalone expand on all isomers.

### 2.7 Three Levels of Validity

An important subtlety: there are three distinct notions of "valid site":

**Expansion-site validity** (`allReductions`, Phase 2):
Checks whether the graph could be the TARGET of an expansion. Uses degree-5
starting vertices and flanking conditions. These are the sites enumerated by
buckygen's canonical algorithm. They answer "where could we expand?" not
"where should we reduce?"

**Inversion-site validity** (`allInvSites`, Phase 3):
Checks whether a set of vertices forms a strip that can be surgically removed.
Uses strip identification (CW ordering), strip topology cleanliness (all strip
neighbours in {path, strip, tp}), and post-surgery degree validation. These
are the sites used for actual graph reduction.

**The subset relationship**: Every inversion site is also a valid expansion
site (reduction site), but not vice versa. `allInvSites` finds a strict
subset of `allReductions`. The additional constraint is `stripTopologyClean`:
an expansion site may have the right degree-5 endpoints and flanking
conditions, but the strip vertices it implies may have neighbours outside
the expected {path, strip, tp} set, making surgery impossible.

**The mismatch on arbitrary graphs**: On an arbitrary graph (not freshly
expanded), the expansion-site enumeration finds sites where expansion would
be valid, but the corresponding reduction may fail because the strip topology
isn't clean. This was discovered empirically: applying `allReductions` results
directly to `applyReduction` yields 0% success on arbitrary graphs.

**Coincidence on expansion-produced graphs**: On a graph that was just
produced by applying an expansion (i.e., the canonical generation use case),
the two sets coincide. The freshly-added strip vertices are guaranteed to have
clean topology because the expansion just created them with exactly the right
connectivity. This is why buckygen's canonical algorithm can use `allReductions`
directly — it only ever reduces graphs that are the result of its own
expansions.

**The resolution**: For inverting *arbitrary* graphs, the strip-based approach
(`allInvSites`) is required. It identifies strip vertices from their CW
neighbourhood structure and validates strip topology before attempting surgery.
This achieves 100% success on all tested isomers.

### 2.8 Computational Complexity

**Graph-copying approach (invertReduction): O(N^2) per isomer.**

The Phase 3 graph-copying implementation has O(N) per step (full adjacency
copy + renumbering) times O(N) steps = O(N^2).

**In-place approach (ReducibleDual): O(N) per isomer.**

The ReducibleDual implementation (Section 2.4) achieves O(N) total: each
vertex is killed at most once (O(1) per kill), each splice/snip touches
O(1) vertices, site finding is O(1) amortised. See Section 2.4.4 for
the full analysis.

### 2.9 Empirical Scaling Results

**Phase 3 (graph-copying, invertReduction):**

Tested on all isomers C32--C80 (131,191 isomers). Every isomer reduces to
a seed. No B(i,j) with i+j > 0 was needed.

| Range | Isomers | Failures | Time (s) |
|-------|---------|----------|----------|
| C32-C40 | 84 | 0 | 0.03 |
| C32-C60 | 5,762 | 0 | 3.8 |
| C32-C80 | 131,191 | 0 | 186 |

**ReducibleDual (in-place, O(N)):**

Tested on all isomers C20--C100 (1,456,593 isomers). Zero failures in
all three modes: reduction, round-trip expansion, standalone expansion.

| Operation | Time (s) | Description |
|-----------|----------|-------------|
| Reduction | 42.1 | Reduce all 1.46M isomers to seeds |
| Round-trip expansion | 4.9 | Expand from saved state, verify exact match |
| Standalone expansion | 49.5 | Reduce + expand from extension path only |

Seed distribution (C20--C100): C20 = 1,328,794 (91.2%); C28 = 49,458
(3.4%); C30 = 78,341 (5.4%).


## 3. Code Analysis

### 3.1 Implementation Structure

The implementation lives in two files:

- **buckinverse.hh**: Type definitions, inline navigation primitives,
  function declarations.
- **buckinverse.cc**: Seed construction, path computation, reduction
  enumeration (Phase 2), and strip-based inversion surgery (Phase 3).

The code operates on the main project's `Graph` class (parent of
`Triangulation`), which provides:
- `neighbours[u][i]`: CW-ordered adjacency list
- `degree(u)`: vertex degree
- `next(u, v)` / `prev(u, v)`: nextCW / prevCW
- `arc_ix(u, v)`: position of v in N(u), or -1 if not adjacent
- `N`: number of vertices

### 3.2 Type Definitions

```cpp
enum class Dir { DRight, DLeft };

struct ExpKind {
    enum Type { L_type, B_type, F_type };
    Type type;
    int i, j;
    int reductionLength() const;  // L_i: i+1, B_{i,j}: i+j+2, F: 0
    int newVertices() const;      // L_i: i+2, B_{i,j}: i+j+3, F: 5
};

struct Expansion {
    ExpKind kind;
    node_t u, v;   // Starting directed edge
    Dir dir;
};
using Reduction = Expansion;  // Same triple, inverse semantics

struct PathInfo {
    vector<node_t> path;      // Main path vertices
    vector<node_t> parallel;  // Parallel / strip vertices
    bool valid = false;        // false if self-intersection detected
};

// Validated inversion site with all vertices identified (Phase 3).
struct InvSite {
    ExpKind kind;
    Dir dir;
    vector<node_t> strip;    // Strip vertices to remove
    vector<node_t> path;     // Path/ring vertices (remain, get reconnected)
    vector<node_t> tp;       // True parallel/outer vertices (remain, get reconnected)
};
// For F-type: path = ring, tp = outer (cap) vertices.
```

Factory functions `Lk(i)`, `Bk(i,j)`, `Fk()` construct `ExpKind` values.

**ReducibleDual types (in-place reduction/expansion):**

```cpp
struct ReducibleDual {
    struct Vertex {
        node_t nbr[6];     // CW-ordered neighbours (fixed positions)
        uint8_t active;      // Bitmask: bit i set iff nbr[i] is present
    };

    struct ReductionStep {
        InvSite site;
        vector<pair<node_t, Vertex>> saved;  // strip states before kill
    };
    // ...
};

struct ExtensionStep {
    ExpKind kind;
    Dir dir;
    vector<node_t> strip;   // vertices to create
    vector<node_t> path;    // existing vertices (reconnected)
    vector<node_t> tp;      // existing vertices (reconnected)
};

struct ExtensionPath {
    SeedType seed;
    int full_N;
    struct SeedVertex { node_t id; node_t nbr[6]; uint8_t active; };
    vector<SeedVertex> seed_state;
    vector<ExtensionStep> steps;   // expansion order
};
```

### 3.3 Navigation Primitives

All navigation functions are inline and map directly to the mathematical
definitions in Section 1.2:

| Function | Formula | Notes |
|----------|---------|-------|
| `advanceCW(g, u, v, k)` | N(u)[(pos(v) + k) mod deg(u)] | Double-modulo for negative k |
| `straightAhead(g, d, u, from)` | advanceCW(u, from, k) | k=3 always at deg-6; k=3 (DRight) or k=2 (DLeft) at deg-5 |
| `turnAhead(g, d, u, from)` | advanceCW(u, from, k) | k=2 (DRight), k=deg-2 (DLeft) |
| `sideNbr(g, d, a, b)` | prev(a,b) or next(a,b) | DRight: prev, DLeft: next |

### 3.4 Path Computation

Three functions compute paths, each returning a `PathInfo`:

**computeStraightPath(g, u, v, dir, numEntries)**:
Iterates straightAhead to build a path of numEntries vertices. The parallel
path has numEntries entries (matching the Haskell convention): one per path
vertex, where the last entry uses the "continuation" direction beyond the
path endpoint. This last entry is the "far parallel" — the true parallel
at the far endpoint, not a strip vertex.

**computeBentZeroPath(g, u, v, dir)**:
Builds the fixed-structure B_{0,0} path:
```
path = [u, v, straightAhead(v,u), turnAhead(path[2],path[1]), straightAhead(path[3],path[2])]
par  = [sideNbr(path[0],path[1]), sideNbr(path[2],path[3]), sideNbr(path[3],path[4])]
```
The parallel path skips the edge entering the turn vertex (edge 1->2).

**computeBentPath(g, u, v, dir, bi, bj)**:
Generalised bent path with bi pre-bend straight steps and bj post-bend
straight steps. Path has bi+bj+5 entries. Parallel has bi+bj+3 entries,
skipping the edge entering the turn vertex.

All three functions check for self-intersection (duplicate vertices in
the combined path + parallel set) and set `valid = false` if found.

### 3.5 Reduction Enumeration

**allReductions(g, maxRedLen)** returns a vector of all valid reduction sites.

The enumeration follows the structure described in Section 2.1, using
helper functions:

| Helper | Purpose | Mapping to Haskell |
|--------|---------|-------------------|
| `isValidL0Direction` | Flanking check for L_0 | `isValidL0Direction` in Canonical.hs |
| `followStraightToFive` | Walk to next degree-5 | `followStraightToFive` in Canonical.hs |
| `isValidStraightReduction` | Flanking check for L_i | `isValidStraightReduction` in Canonical.hs |
| `isValidB00Reduction` | B_{0,0} validity | `isValidB00Reduction` in Canonical.hs |
| `isValidBentReduction` | B_{i,j} validity | `isValidBentReduction` in Canonical.hs |
| `bentEndpointFlank` | Flank vertex for bent | `bentEndpointFlank` in Canonical.hs |
| `computeBentPathSafe` | Uniform parallel path | `computeBentPathSafe` in Expansion.hs |

**Enumeration order and bounds**:
- L_0: enumerated unconditionally
- L_i (i >= 1) and B_{0,0}: require maxRedLen >= 2
- B_{i,j} (i+j > 0): require maxRedLen >= 3; bentLen ranges 1 to maxRedLen-2

### 3.6 Seed Construction

Seeds are constructed via the main project's `FullereneDual(N, rspi)` function,
which builds a Triangulation from the ring spiral pentagon indices:

```cpp
Triangulation makeSeedC20() { return FullereneDual(20, {0,1,2,3,4,5,6,7,8,9,10,11}); }
Triangulation makeSeedC28() { return FullereneDual(28, {0,1,2,4,6,8,9,10,11,12,13,14}); }
Triangulation makeSeedC30() { return FullereneDual(30, {0,1,2,3,4,5,11,12,13,14,15,16}); }
```

The RSPIs are 0-indexed (the isomer database stores 1-indexed values).
`identifySeed(t)` identifies seeds by dual vertex count (12/16/17).

### 3.7 Mapping from Haskell Reference

The implementation is a faithful translation of the Haskell reference in
`buckygen-revival/src/Expansion.hs` and `Canonical.hs`, with the following
adaptations:

1. **Data structure**: Haskell uses an IntMap of lists for adjacency;
   C++ uses the existing `Graph` class with `vector<vector<node_t>>`.

2. **Immutability**: Haskell functions produce new IntMaps at each step;
   C++ will modify adjacency lists in-place during surgery (Phase 3).

3. **Vertex numbering**: Haskell assumes strip vertices are the last k
   vertices (nv-k..nv-1), which is only true for graphs just expanded.
   The C++ inversion code identifies strip vertices from their CW
   neighbourhood structure (see Section 2.2.2), avoiding this assumption.

4. **Path convention**: `computeStraightPath` now produces numEntries
   parallel entries (matching the Haskell convention), including a
   "continuation" entry for the far parallel vertex.

5. **Inversion approach**: The Haskell `reduceL0`/`reduceStraight`/etc.
   operate on the PathInfo from the pre-expansion graph. The C++
   `invertReduction` operates on InvSite structs that contain the strip,
   path, and tp vertices identified directly from the expanded graph's
   CW structure.

### 3.8 Test Results

**Phase 1** (navigation, seeds, paths): **2059 tests passed, 0 failed.**

**Phase 2** (reduction enumeration): **1037 tests passed, 0 failed.**
Validated on buckygen-generated fullerenes C32--C40:

| Fullerene | Isomers | Total reductions | L_0 | L_i | B_{0,0} | B_{i,j} |
|-----------|---------|------------------|-----|-----|---------|---------|
| C32       | 6       | 676              | 5   | 9   | 6       | 1       |
| C34       | 6       | 642              | 6   | 10  | 6       | 1       |
| C36       | 15      | 1654             | 14  | 25  | 15      | 11      |
| C38       | 17      | 1868             | 16  | 43  | 17      | 25      |
| C40       | 40      | 4293             | 38  | 85  | 40      | 71      |

**Phase 3** (inversion surgery): **257 tests passed, 0 failed.**
All 84 isomers (C32--C40) iteratively reduce to seeds.

**Phase 4** (round-trip expansion): **1,456,593 isomers, 0 failures.**
Reduce to seed recording ReductionSteps, expand in reverse, verify exact
adjacency match. Tested on all isomers C20--C100.

**Phase 5** (standalone expansion): **1,456,593 isomers, 0 failures.**
Reduce to ExtensionPath, reconstruct from seed + path only, verify
CW-rotation-equivalent adjacency. Tested on all isomers C20--C100.

Key observations:
- Every non-seed isomer has at least one valid inversion site
- L_0 inversions are the most common, present on almost every isomer
- Only nanotube-family isomers require F-ring reduction
- Seed graphs identified by dual vertex count: 12 (C20), 16 (C28), 17 (C30)
- B(i,j) with i+j > 0 not needed for any isomer up to C100


## 4. Implementation Notes

### 4.1 The Arbitrary-Graph Challenge and its Resolution

The Haskell reduction functions (reduceL0, reduceStraight, etc.) assume
strip vertices are the highest-numbered vertices in the graph. This fails
for arbitrary graphs where strip vertices can be anywhere.

Two approaches were attempted:

**Failed approach (expansion-site based)**: Use `allReductions` to find
expansion sites, compute path/parallel, identify parallel as strip, find
true parallel from strip neighbourhoods, apply reconnection. This achieved
0% success because on arbitrary graphs the parallel path vertices often
have neighbours outside the expected {path, strip, tp} set, causing
reconnection to fail.

**Working approach (strip-based)**: Identify strip vertices directly from
their CW neighbourhood pattern, extract path and tp vertices from the CW
ordering, validate strip topology, apply reconnection. This achieves 100%
success on all tested isomers. The key insight is that the CW ordering
around each strip vertex encodes exactly which vertices are path and which
are true parallel (Section 2.2.2--2.2.4).

### 4.2 Strip Topology Cleanliness

A crucial validation step: before attempting surgery, verify that every
neighbour of every strip vertex is in {path, strip, tp}. If any strip
vertex has an "extra" neighbour, the site doesn't have clean expansion
topology and cannot be surgically inverted. This check eliminates false
positives from the CW-based identification.

### 4.3 Post-Surgery Validation

After reconnection, verify:
1. All non-strip vertices have degree 5 or 6
2. No non-strip vertex has any strip vertex remaining in its neighbour list
3. The result constructs a valid Graph with `is_oriented = true`

If any check fails, `invertReduction` returns Graph() (N == 0) to signal
failure, and the caller tries the next inversion site.

### 4.4 Functions Summary

**Phase 3 (graph-copying inversion):**

| Function | Purpose |
|----------|---------|
| `findL0InvSites(g)` | Find L_0 strips: adjacent degree-5 pairs |
| `findLiInvSites(g, maxRedLen)` | Find L_i strips: deg-5/deg-6.../deg-5 chains |
| `findB00InvSites(g)` | Find B_{0,0} strips: deg-5/deg-6/deg-5 triples |
| `findFRingInvSites(g)` | Find F strips: parallel degree-6 5-cycles |
| `allInvSites(g, maxRedLen)` | Combine all of the above |
| `invertReduction(g, site)` | Apply surgery, return reduced graph |
| `reconnectStraight(adj, path, strip, tp)` | L-type edge rewiring |
| `reconnectBent(adj, path, strip, tp, ...)` | B-type edge rewiring |
| `stripTopologyClean(g, path, strip, tp)` | Validate strip neighbourhood |

**ReducibleDual (in-place O(N) reduction/expansion):**

| Function | Purpose |
|----------|---------|
| `ReducibleDual(const Graph& g)` | Construct from graph, O(N) |
| `ReducibleDual(int capacity)` | Construct all-dead array for standalone expansion |
| `splice(u, old_v, new_v)` | Replace active neighbour, O(1) |
| `snip(u, v)` | Clear neighbour bit, O(1) |
| `unsnip(u, v)` | Restore inactive neighbour bit, O(1) |
| `kill(v)` / `unkill(v, saved)` | Delete / restore vertex, O(1) |
| `insertAfter(u, v, after)` | CW insertion, degree 5 -> 6, O(1) |
| `findSite(maxRedLen)` | Find first valid inversion site, O(1) amortised |
| `reduce(site)` | Apply reconnection + kill, O(1) |
| `expand(ReductionStep)` | Round-trip: unkill + reverse reconnect, O(1) |
| `expand(ExtensionStep)` | Standalone: CW formulas + insertAfter, O(1) |
| `reduceToSeed()` | Full reduction loop, O(N) |
| `reduceToSeed(path)` | With ReductionStep recording, O(N) |
| `reduceToExtensionPath()` | Reduce and produce ExtensionPath, O(N) |
| `toGraph()` | Extract compacted Graph, O(N) |
| `graphFromExtensionPath(ep)` | Reconstruct from seed + path, O(N) |

### 4.5 Not Yet Implemented

- **B_{i,j} with i+j > 0 standalone expansion**: B(0,0) is implemented.
  General B(i,j) inversion works in Phase 3 but has not been needed for any
  isomer up to C100. If needed, the standalone expand formulas would
  generalise from the B(0,0) case with longer pre-bend and post-bend loops.
- **Canonical reduction selection**: The 5-tuple cascade from buckygen's
  canonical algorithm. Currently any valid inversion site is used.


## Appendix A: Notation Summary

| Symbol | Meaning |
|--------|---------|
| C_n | Fullerene with n carbon atoms |
| T | Dual triangulation (n/2 + 2 vertices) |
| N(u) | CW-ordered neighbour list of vertex u |
| deg(u) | Degree of vertex u (5 or 6) |
| L_i | Straight extension/reduction of index i |
| B_{i,j} | Bent extension/reduction with pre-bend index i, post-bend index j |
| d | Direction: DLeft or DRight |
| (u, v, d) | Reduction site: directed edge plus direction |
| S = [s_0, ..., s_{k-1}] | Strip of new vertices inserted by extension |
| P = [p_0, ..., p_m] | Main path through the triangulation |
| Q = [q_0, ..., q_m] | Parallel path (side neighbours of P) |

## Appendix B: Degree Changes Summary

| Operation | Vertices gaining degree | Vertices losing degree | Vertices removed |
|-----------|------------------------|----------------------|-----------------|
| L_i expansion | path[0]: 5->6, par[i+2]: 5->6 | (none) | (none) |
| L_i reduction | (none) | path[0]: 6->5, tp[i+2]: 6->5 | strip[0..i+1] |
| B_{i,j} expansion | path[0]: 5->6, path[i+j+4]: 5->6 | (none) | (none) |
| B_{i,j} reduction | (none) | path[0]: 6->5, path[i+j+4]: 6->5 | strip[0..i+j+2] |
| F expansion | (none) | (none) | (none, adds 5 deg-6) |
| F reduction | (none) | (none) | strip[0..4] (all deg-6) |
