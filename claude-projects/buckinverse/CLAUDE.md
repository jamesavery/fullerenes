# Buckinverse — Fullerene Reduction Path Finder

## Current Status (Standalone expansion complete, C20-C100 verified)

Phases 1-5 complete. All 1,456,593 isomers C20-C100 verified.
- Phase 1: 2059 passed (navigation, seeds)
- Phase 2: 1037 passed (reduction enumeration)
- Phase 3: 257 passed (inversion surgery)
- Phase 4: Round-trip (reduce + expand via saved state): 0 failures on 1.46M isomers
- Phase 5: Standalone expansion (seed + ExtensionPath only): 0 failures on 1.46M isomers
- Seed tally: C20=1,328,794  C28=49,458  C30=78,341

## Implemented Inversion Types

### L0: Adjacent degree-5 pair
Strip = [s0, s1] (both deg-5, adjacent). CW extraction of path[0..2] and tp[0..2].

### L_i (i >= 1): Straight strip chain
Strip = [s0(deg5), s1..s_i(deg6), s_{i+1}(deg5)]. Chain walk via advanceCW(3) at
degree 6 (direction-independent). Incremental path/tp from CW ordering at each vertex.

### B(0,0): Bent strip triple
Strip = [a(deg5), b(deg6), c(deg5)] where a adj b, b adj c, a NOT adj c.
CW extraction from all 3 vertices with direction-dependent validation.

### F (nanotube ring): Equatorial ring of 5
Strip = 5-cycle of degree-6 vertices between ring (another deg-6 5-cycle) and
5 outer (cap) vertices. Detection: trace straight-ahead 5-cycle, find common
non-ring neighbors for strip, resolve outers by intersection of extras.
Reconnection: ring[i] replaces strip refs with outer refs; outer[i] replaces
strip refs with ring refs. C40 #40 (the (5,0) nanotube) uses this.

## Key Insight: Strip-Based Inversion

The original `allReductions` finds EXPANSION sites (degree-5 starting vertices), NOT reduction
sites. For arbitrary graph inversion, start from STRIP vertices (the vertices added by expansion):

- **L0/L_i strip**: Identified by degree-5 endpoints with deg-6 interior chain
- **B(i,j) strip**: Identified by deg-5/deg-6/deg-5 triple with turn
- **F strip**: Identified by two parallel 5-cycles of degree-6 vertices

### CW Ordering Around Strip Vertices (L-type)

**s0 (degree 5, first endpoint):**
- DRight CW: [path[0], tp[0], tp[1], s1, path[1]]
- DLeft CW: [path[0], path[1], s1, tp[1], tp[0]]

**Interior s_j (degree 6):**
- DRight CW: [s_{j-1}, tp[j], tp[j+1], s_{j+1}, path[j+1], path[j]]
- DLeft CW: [s_{j-1}, path[j], path[j+1], s_{j+1}, tp[j+1], tp[j]]

**Far endpoint s_{i+1} (degree 5):**
- DRight CW: [s_i, tp[i+1], tp[i+2], path[i+2], path[i+1]]
- DLeft CW: [s_i, path[i+1], path[i+2], tp[i+2], tp[i+1]]

### F-Ring Reconnection
- ring[i]: replace strip[(i-1)%5] → outer[(i-1)%5], strip[i] → outer[i]
- outer[i]: replace strip[(i-1)%5] → ring[i], strip[i] → ring[(i+1)%5]
- Outer resolution: outer[i] = common element of strip[i-1]'s extras and strip[i]'s extras

## Architecture

### Files
- `buckinverse.hh` — Types (ExpKind, Expansion, InvSite, PathInfo, ExtensionPath), navigation inlines, declarations
- `buckinverse.cc` — All implementation: seeds, paths, reduction, inversion, expansion, extension paths
- `test_phase1.cc` — Navigation and seed tests (2059 pass)
- `test_phase2.cc` — Reduction enumeration tests (1037 pass)
- `test_phase3.cc` — Inversion surgery tests (257 pass)
- `test_roundtrip.cc` — Round-trip + standalone expansion tests (1.46M isomers, 0 failures)
- `Makefile` — Direct build (g++, links against libfullerenes.so from ../../build2)

### Key Types
```
InvSite { kind, dir, strip[], path[], tp[] }    — validated inversion site
  For L/B: strip=removed vertices, path=boundary, tp=parallel (reconnected)
  For F: strip=removed ring, path=kept ring, tp=outer cap vertices
ExpKind { L_type/B_type/F_type, i, j }          — expansion type
Expansion/Reduction { kind, u, v, dir }          — expansion-site formulation (Phase 2)
ExtensionStep { kind, dir, strip[], path[], tp[] } — one expansion step (vertex IDs in full graph)
ExtensionPath { seed, full_N, seed_state[], steps[] } — portable seed-to-fullerene route
ReducibleDual::ReductionStep { site, saved[] }   — reduction step with saved strip vertex states
```

### Key Functions
- `findL0InvSites(g)` — Find L0 strip pairs (adjacent deg-5)
- `findLiInvSites(g, maxRedLen)` — Find L_i strip chains (static, in .cc)
- `findB00InvSites(g)` — Find B(0,0) strip triples (static, in .cc)
- `findFRingInvSites(g)` — Find F ring strip 5-cycles
- `allInvSites(g, maxRedLen)` — All inversion sites (L0 + L_i + B00 + F)
- `invertReduction(g, site)` — Apply surgery, return reduced graph
- `reconnectStraight(adj, path, strip, tp)` — Reconnection for L-type
- `reconnectBent(adj, path, strip, tp, bentPos, bentLen)` — Reconnection for B-type
- `allReductions(g)` — Phase 2 expansion-site enumeration (NOT for surgery)

## Conventions
- All graphs oriented (is_oriented == true). Never call orient_neighbours().
- CW-ordered neighbour lists. next(u,v) advances CW, prev(u,v) goes CCW.
- Dual has n/2+2 vertices for fullerene C_n.
- Degree-5 = pentagon, degree-6 = hexagon. Exactly 12 degree-5 vertices.
- Three seeds: C20 (12 dv), C28 (16 dv), C30 (17 dv).

## Standalone Expansion Architecture

The ExtensionPath representation is self-contained: (seed enum, seed_state, expansion steps).
- seed_state stores full Vertex data (nbr[6] + active bitmask) for seed vertices,
  preserving inactive positions needed by unsnip during expansion.
- Each ExtensionStep stores strip/path/tp vertex IDs in the full graph's numbering.
- graphFromExtensionPath() constructs ReducibleDual from seed state, applies steps, returns Graph.
- Strip vertices are created from CW formulas (verified on all 1.46M isomers).
- insertAfter primitive promotes degree-5 to degree-6 by CW insertion (replaces unsnip for standalone).
- Standalone expansion produces CW-rotation-equivalent output (same circular neighbor ordering,
  potentially different physical starting position in neighbor array).

## Not Yet Implemented
- B(i,j) with i+j > 0 inversion (not needed for C32-C100, may be for larger)
- Canonical reduction selection (deterministic choice among multiple valid reduction sites)

## Rules
- Follow parent project rules from /home/avery/work/fullerenes/CLAUDE.md
- No git commits (show proposed message, let user commit)
- No backticks in commit messages
- Never call orient_neighbours() or use zero_order_geometry()
