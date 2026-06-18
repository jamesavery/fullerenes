"""FullereneDual surface-geodesics bindings: surface_distances + simple geodesics
over the intrinsic equilateral metric (Eisenstein straight-line walks)."""
import numpy as np

import fullerenes as fl


def _c20_dual():
    return fl.FullereneGraph.C20().dual()      # icosahedron dual, 12 vertices


def test_surface_distances_is_a_metric_matrix():
    d = _c20_dual()
    D = d.surface_distances()
    assert D.shape == (d.N, d.N) and D.dtype == np.float64
    assert np.allclose(D, D.T)                          # symmetric
    assert np.allclose(np.diag(D), 0.0)                 # zero diagonal
    off = D[~np.eye(d.N, dtype=bool)]
    assert (off > 0).all()                              # distinct vertices are apart


def test_simple_geodesics_from_reaches_every_other_vertex():
    d = _c20_dual()
    g = d.simple_geodesics_from(0)
    assert {v for (v, _ax, _a, _b, _d2) in g} == set(range(1, d.N))   # all others, not self
    for v, axis, a, b, d2 in g:
        assert d2 == a * a + a * b + b * b              # Eisenstein norm
        assert d2 > 0


def test_end_of_the_line_unit_step_is_a_neighbour():
    d = _c20_dual()
    assert d.end_of_the_line(0, 0, 1, 0) in set(d.adjacency()[0])     # one (1,0) step


def test_geodesic_strip_starts_at_source_and_is_consistent():
    d = _c20_dual()
    strip = d.geodesic_strip(0, 0, 1, 1)                # raises if terminus disagrees (fail-loud)
    assert strip and all(len(run) > 0 for run in strip)
    v0, x0, y0 = strip[0][0]
    assert v0 == 0 and (x0, y0) == (0, 0)               # source at the lattice origin


def test_simple_geodesics_trace_invariants():
    d = _c20_dual()
    M, steps = d.simple_geodesics_trace(0)              # (search radius, per-ray steps)
    assert isinstance(M, int) and M >= 1 and len(steps) > 0
    for axis, a, b, v, d2, H_before, improved in steps:
        assert d2 == a * a + a * b + b * b              # Eisenstein norm
        # Non-tautological cross-checks (vs an independent binding + the returned
        # radius), so a field-mapping or wrong-M bug can't pass:
        assert v == d.end_of_the_line(0, axis, a, b)    # the ray's landing vertex
        assert a <= M and d2 <= M * M                   # every probe bounded by the radius
        assert improved == (d2 <= H_before)             # witness (re)assigned iff <= prior
