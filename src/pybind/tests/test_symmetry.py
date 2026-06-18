"""Phase 5: Symmetry (point group, orbits, 3D representation) and transforms."""
import numpy as np

import fullerenes as fl


def test_c20_is_icosahedral():
    sym = fl.Symmetry(fl.FullereneGraph.C20().dual())
    assert sym.point_group() == "Ih"
    assert sym.order() == 120                         # |Ih| = 120


def test_c60_buckminsterfullerene_is_icosahedral():
    # The unique IPR C60 (or buckygen's first C60) -- both are the Ih cage.
    dual = next(iter(fl.buckygen(60, IPR=True)))
    sym = fl.Symmetry(dual)
    assert sym.point_group() == "Ih"
    assert sym.order() == 120


def test_equivalence_classes_partition_vertices():
    dual = fl.FullereneGraph.C20().dual()
    sym = fl.Symmetry(dual)
    classes = sym.equivalence_classes()
    flat = sorted(v for cls in classes for v in cls)
    assert flat == list(range(dual.N))                # a partition of all vertices
    # Ih C20 dual: all 12 vertices in one orbit
    assert len(classes) == 1


def test_representation_3d_is_orthogonal():
    sym = fl.Symmetry(fl.FullereneGraph.C20().dual())
    R = sym.representation_3d()
    assert R.shape == (sym.order(), 3, 3) and R.dtype == np.float64
    for M in R:
        assert np.allclose(M @ M.T, np.eye(3), atol=1e-9)     # orthogonal
        assert abs(abs(np.linalg.det(M)) - 1.0) < 1e-9        # det = +/-1


def test_point_group_subset_of_known_symmetries():
    syms = set(s.strip() for s in fl.symmetries(60))
    for dual in fl.buckygen(60, IPR=True):
        assert fl.Symmetry(dual).point_group().strip() in syms


# --- transforms ---

def test_fullerene_graph_transforms():
    c20 = fl.FullereneGraph.C20()
    lf = c20.leapfrog()                                # C20 -> C60
    assert lf.N == 60 and lf.is_a_fullerene()
    gc = c20.gc_transform(2, 0)                        # C20 -> C80
    assert gc.N == 80 and gc.is_a_fullerene()
    assert lf.name() == c20.gc_transform(1, 1).name()  # leapfrog == GC(1,1)


def test_fullerene_dual_transforms_round_trip_with_graph():
    c20d = fl.FullereneGraph.C20().dual()
    # GC on the dual matches GC on the cubic then dualizing (same canonical name).
    assert c20d.gc_transform(2, 0).name() == fl.FullereneGraph.C20().gc_transform(2, 0).name()
