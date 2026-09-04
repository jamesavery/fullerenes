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
retriangulating every cocircular cell that has interior edges — subject to
the two refusal classes below — as the FAN from its canonical corner: the
corner whose boundary rotation word is lexicographically minimal, keyed on
(vertex id, exact squared length). (`canonical_tesselation` compares by
the same key shape over the caller's label map; the orders coincide when
that map is monotone, identity included.) Fan conversion uses only tight
flips — zero-energy Delaunay moves inside the cell's circle — so the
tesselation, the SURFACE metric, and every vertex cone angle are unchanged
(the edge-length field is not: each flip replaces one diagonal of a cyclic
quadrilateral with the other). The completed triangulation of every fanned
cell is then a function of the labeled input complex alone.

Two refusal classes are counted by name and left untouched, never guessed
at:

- **ambiguous** — the boundary word is periodic (several corners share the
  minimal rotation): no label-determined apex exists; such isomers compare
  at the tesselation level.
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

## 4. When ambiguity can occur at all

Let a cocircular cell have d corners (d >= 4: it has an interior edge) and
boundary word of period p | d, p < d, with k = d/p >= 2 repeats.
Developed injectively, the cell is a polygon inscribed in a circle with
interior angle sum (d−2)π. Each corner is a sector at a cone; no cone is
interior (§3), so the cell's total angle at any one cone is strictly less
than that cone's angle — 5π/3 on a REDUCED FULLERENE DUAL (every live
vertex a 5π/3 cone), the hypothesis this bound needs. Period p makes the
corner labels p-periodic, so the corners fall into p residue classes of k
corners each; two classes may share a cone, so let q <= p be the number of
distinct cones. Summing corners by cone:

    (d − 2)π < q · 5π/3 <= p · 5π/3   ⟹   k < 5/3 + 2/p.

p = 1 gives k = d < 11/3, i.e. d <= 3 — excluded by d >= 4. For p >= 2
this forces **k = 2**, and then p < 6: **p ∈ {2,3,4,5}, d ∈ {4,6,8,10}**.
Sharper: with d = 2p, q <= p−1 would give 2p − 2 < 5(p−1)/3, i.e. p < 1 —
impossible — so q = p and **every cone appears at exactly two corners**,
at antipodal positions i and i + d/2: the cell is invariant under the
half-turn about its circumcenter, and the chord joining each antipodal
pair is a diameter of the circumcircle — a geodesic loop at that cone,
which the current triangulation may or may not realize as an edge. (With
flat 2π vertices instead of 5π/3 cones the same computation still forces
k = 2 but no longer bounds p; and a non-embedded cell — an empty disk
wrapping past the injectivity radius — relaxes the angle sum, so the
class is constrained, not impossible, off the reduced-dual hypothesis.)

The minimal lattice realization of the k = 2 half-turn shape is an
inscribed rectangle: sides² {1, 3}, diameter² 4 (two 30-60-90 lattice
triangles, tau = 2, F == 0 on the diagonal exactly).
`claude-projects/delaunay/tools/test_canonical_completion_branches` builds a flat
torus from two such rectangles — its vertices are flat 2π points, so it
realizes the ambiguity BRANCH rather than §4's reduced-dual hypothesis
class — whose cells carry period-2 boundary words, (u,1),(x,3),(u,1),(x,3)
and its u↔x partner, and asserts the completion refuses both as
ambiguous: the branch is exercised, not merely trusted.

Empirically the class is EMPTY over the full C20–C100 reduced space
(0 ambiguous; the validation record above) — measured fact, not theorem;
the validator's independent period cross-check (computed from the
pre-completion tesselation) keeps under-detection falsifiable.

## 5. What the gates compare

Per PRIMITIVES-SPEC §5.0 (parallel-primitives, ruling 2026-08-23):
Delaunay-bearing results conform at the canonical-tesselation level (A3);
after canonical completion, the exact-regime triangulation is unique
wherever no cell is refused (ambiguous or nondisk), and completed
complexes compare by the canonical DCEL word (lex-min BFS serialization —
a complete invariant of a CONNECTED labeled complex; the face-multiset
form is not, since it forgets the gluing). DCEL bytes remain meaningful
only as same-schedule regression detectors and under the forced-order
debugging mode.
