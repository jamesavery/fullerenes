# CANONICAL-TESSELATION.md — the exact cocircularity predicate, the canonical Delaunay tesselation, and its canonical completion

Referenced from `src/c++/delaunay.cc` (`canonical_tesselation`) and
`include/fullerenes/delaunay_geometry.hh` (`Diamond::is_cocircular_exact`).
The implementation lives in those files plus
`include/fullerenes/delaunay_view.hh` (`canonical_completion`,
`next_cell_boundary`); the full-space validation record (C20–C100,
1,456,598 isomers, 0 failures on six gates) is
`claude-projects/delaunay/tools/validate_canonical_completion` — see the
parallel-primitives README arc log of 2026-08-23.

## 1. The exact cocircularity predicate (F == 0)

For an edge e with diamond sides (e, a, b) above and (e, c, d) below
(squared lengths; B and D name the two apexes, opposite e), define

    F = s_upper * sqrt(H_lower) + s_lower * sqrt(H_upper),

where s_upper = a + b − e and s_lower = c + d − e are the law-of-cosines
numerators at the two apexes and H = 16·Area² is the Heron product of each
triangle. For non-degenerate triangles (H_upper, H_lower > 0),

    cot(angle_B) + cot(angle_D) = F / sqrt(H_upper · H_lower)

with a positive denominator, so sign(F) is the local Delaunay
classification: F > 0 strictly locally Delaunay, F == 0 cocircular
("tight" — the four diamond points lie on one circle), F < 0 must flip.
A diamond with a degenerate or non-lattice triangle (tau <= 0) is outside
the predicate's domain and is refused, never classified.

On an equilateral-derived metric every face is an Eisenstein-lattice
triangle — the standing hypothesis the exact regime's entry boundary
establishes (`derive_exact_lsq_carry`); the Heron–Eisenstein identity then
gives H = 3·tau² with tau the integer lattice area number, so

    F = sqrt(3) · (s_upper·tau_lower + s_lower·tau_upper).

That integer linear form (it returns F/sqrt(3)) is
`DiamondSq::delaunay_form()`: the classification, cocircularity included,
is the sign of one integer, tolerance-free.
`Diamond::is_cocircular_exact()` enters it through the integrality trust
boundary (`lsq_integrality_band`); `ExactIntegerMetric::cocircular` reads
it from the verified Lsq carry directly.

## 2. The canonical tesselation

By Bobenko–Springborn, a piecewise-flat surface has a unique intrinsic
Delaunay TESSELATION: the cells are the maximal cocircular polygons, and
any intrinsic Delaunay triangulation refines it by triangulating each cell
arbitrarily. Tight edges are exactly the cell-interior edges; every
non-tight half-edge lies on the boundary of exactly one cell.

`canonical_tesselation` materializes this object canonically: each cell's
boundary is walked (crossing tight edges with the cell-boundary step
`next_cell_boundary`), recorded as the cyclic word of
(vertex label, integer squared boundary-edge length) entries, normalized
to its lexicographically minimal rotation, and the cells are sorted. Two
complexes over the same labels that refine the same tesselation have equal
canonical tesselations. The converse is not claimed: the cell multiset
forgets the gluing (see §5), so equality is the operative conformance test
(spec §5.0's A3), not a proof of identical tesselations. `fingerprint()`
is the FNV-seeded hash combine for certificates; gates compare the full
cells.

## 3. The canonical completion

Downstream consumers read a triangulation, not a tesselation, so the
per-cell triangulation freedom is a reproducibility hole: two runs that
flip in different orders hand downstream different (equally valid)
refinements. `DelaunayView::canonical_completion` closes it by
retriangulating every cocircular cell that has interior edges from its
CANONICAL CORNER: the corner whose boundary rotation word is
lexicographically minimal, keyed on (vertex id, exact squared length).
(`canonical_tesselation` compares by the same key shape over the caller's
label map; the orders coincide when that map is monotone, identity
included.) Conversion uses only tight flips — zero-energy Delaunay moves
inside the cell's circle — so the tesselation, the SURFACE metric, and
every vertex cone angle are unchanged (the edge-length field is not: each
flip replaces one diagonal of a cyclic quadrilateral with the other). The
completed triangulation is then a function of the labeled input complex
alone.

How many corners attain the minimum selects the construction:

- **one** — the cell becomes the FAN from that corner. `fanned` counts
  these; it is the only case any fullerene surface has yet produced.
- **two** — §4 shows the two corners are then antipodal and the chord
  joining them is a diameter of the circumcircle. The cell becomes that
  diameter plus the fan of each half from its own end. Read from either
  end the diagonal set is the same, so nothing is chosen between the two
  minima and the result is again a function of the labeled complex.
  `periodic_completed` counts these. The construction, its proof, its
  O(d) cost and its flip bound are in
  `claude-projects/delaunay/CANONICAL-TOTALITY-DESIGN.md`.
- **three or four** — §4 excludes this on every convex polyhedral
  metric, which is every surface this pipeline reads; the known
  realizations are flat tori consisting of a single cell. `ambiguous`
  counts them, they are left untouched, and such a complex compares at
  the tesselation level.

One refusal class remains, counted by name and never guessed at:

- **nondisk** — the component fails the disk Euler count. A triangulated
  disk with d boundary edges and no interior vertices has exactly d−3
  interior edges, each crossed once from each side by the boundary walk:
  X == 2·(d−3).

The Euler count is NECESSARY for the fan conversion's domain, not
sufficient: with n interior vertices X = 2(d−3) + 6n − Σ_j deg(w_j), so
interior vertices of total degree 6n would pass it. Correctness on
accepted cells rests on the geometry, not the count: **a cocircular cell
of a Delaunay complex has no interior vertices** — by the empty-circumdisk
property no vertex lies strictly inside a cell's circumcircle, so every
vertex of a cell is a boundary corner. (Equivalently: if every spoke of a
vertex v were tight, chaining the four-point concyclicity conditions puts
v and its whole link on one circle, so v is ON the circle, not interior.)
Both arguments assume the cell develops injectively; a non-embedded cell
(§4) is outside their scope — which is why the gate refuses rather than
reasons, as the cheap fail-loud backstop against corrupt or exotic
complexes, with the fan loop's step budget behind it. On the production
pipeline the point is doubly moot: a reduced complex has no live flat
vertices at all, and the owner entry checks `is_delaunay()` before
running.

## 4. How many corners can attain the minimum

Let a cocircular cell have d corners (d >= 4: it has an interior edge) and
boundary word of least period p | d, p < d, with k = d/p >= 2 repeats.
The word is fixed by exactly k of its d rotations, so k is also the
number of corners attaining the least rotation. Period p makes the
corner labels p-periodic, so the corners fall into p residue classes of
k corners each; two classes may share a vertex, so let q <= p be the
number of distinct vertices.

By §3 the cell has no interior vertices, so it is a triangulated disk of
exactly d−2 triangles and its corner angles sum to (d−2)π — a triangle
count, independent of whether the cell develops injectively. The corners
at one vertex are disjoint sectors of that vertex's link, so

    (d − 2)π = Σ_i α_i <= Σ_{v on C} Θ_v <= 2πq.

**Sphere bound.** If every Θ_v <= 2π and the surface is a sphere, the
inequality is strict: equality throughout would make every vertex on the
cell flat with its whole link inside this one cell, so no other triangle
could meet those vertices, the cell would be the entire surface, and the
total curvature would be 0 rather than 4π. So (d−2)π < 2πq <= 2πp, that
is p(k−2) < 2. Then k = 3 forces p = 1 and d = 3, excluded; k >= 4 fails
outright. So **k = 2** for p >= 2, and 2p − 2 < 2q with q <= p forces
**q = p**, every vertex occurring at exactly two corners.

**The hypothesis is convexity.** By Alexandrov's realisation theorem a
sphere with every cone angle at most 2π is precisely a convex polyhedral
metric, degenerate doubly covered polygons included. So the bound holds
on every convex polyhedral metric and on nothing larger, and it needs no
further cone hypothesis: it covers the fullerene dual and the flattened
kis alike, where **the minimum is attained once or twice, never more**.
Positive total curvature is what the proof spends, which is why the
bound fails on a flat torus (§6 of the design document: one-cell
hexagonal and square tori attain k = 3 and k = 4).

**Reduced dual.** With Θ_v = 5π/3 the same inequality gives 2p − 2 <=
5p/3, so p <= 6; p = 6 would force equality throughout, hence a single
cell with V = E = p and F = 1, Euler characteristic 1, impossible for a
closed oriented surface. So **p ∈ {2,3,4,5}, d ∈ {4,6,8,10}**.

**Reduced flattened kis.** Its surviving vertices carry curvature
r_v·π/15 with r_v ∈ {1,2,3}; with q = p the same argument gives
Σ r_v <= 29, hence **p <= 29 and d <= 58**. Neither 10 nor 58 is a
hard-coded capacity anywhere in the code.

**The half-turn.** Periodicity includes the squared side lengths, and
equal chords of a circle subtend equal gaps, so corresponding gaps agree
and each period subtends 2π/k. For k = 2 the two minima are antipodal and
the chord between them is a diameter of the circumcircle — a geodesic
loop at that vertex, which the current triangulation may or may not
realize as an edge. That diameter is what the two-minimum construction of
§3 finds or creates. The proofs are in
`claude-projects/delaunay/CANONICAL-TOTALITY-DESIGN.md`; the argument
above corrects an earlier version of this section that asserted strict
inequality vertex by vertex without the sphere hypothesis.

The minimal lattice realization of the k = 2 half-turn shape is an
inscribed rectangle: sides² {1, 3}, diameter² 4 (two 30-60-90 lattice
triangles, tau = 2, F == 0 on the diagonal exactly).
`claude-projects/delaunay/tools/test_canonical_completion_branches` builds
two witnesses, both flat tori, so both exercise the branch rather than
§4's cone hypothesis. The first glues two such rectangles; its cells carry
period-2 boundary words, (u,1),(x,3),(u,1),(x,3) and its u↔x partner, and
their diameters are already the whole diagonal set, so the completion
makes no flip — while flipping both diameters away and re-completing must
return the identical labelled complex. The second is the quotient of the
plane by the lattice generated by (2,0) and (1,√3), whose single
cocircular cell is a regular hexagon of period 3: its two minima are
three steps apart and are two OCCURRENCES of one vertex, so the chord
between them is a self-loop that no lookup by vertex identifier could
find. Every triangulation of that hexagon is enumerated by search over
cocircular flips, and all of them complete to one labelled complex.

Empirically the periodic class is EMPTY over the full C20–C100 reduced
space (0 cells with two or more minima; the validation record above) —
measured fact, not theorem; the validator's independent multiplicity
cross-check, computed from the pre-completion tesselation, keeps
under-detection falsifiable.

## 5. What the gates compare

Per PRIMITIVES-SPEC §5.0 (parallel-primitives, ruling 2026-08-23):
Delaunay-bearing results conform at the canonical-tesselation level (A3);
after canonical completion, the exact-regime triangulation is unique
wherever no cell is refused (three or more minima, or non-disk), and
completed
complexes compare by the canonical DCEL word (lex-min BFS serialization —
a complete invariant of a CONNECTED labeled complex; the face-multiset
form is not, since it forgets the gluing). DCEL bytes remain meaningful
only as same-schedule regression detectors and under the forced-order
debugging mode.
