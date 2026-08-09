"""Phase 4: geometry (Polyhedron, Deltahedron) and zero-copy coordinate write-back."""
import numpy as np
import pytest

import fullerenes as fl


# --- Polyhedron --------------------------------------------------------------

def test_polyhedron_c20_geometry():
    P = fl.Polyhedron.C20()
    assert P.N == 20
    assert P.points.shape == (20, 3) and P.points.dtype == np.float64
    assert P.points.base is P and P.points.flags.writeable
    assert P.volume() > 0 and P.surface_area() > 0 and P.diameter() > 0
    assert len(P.faces()) == 12                       # 12 pentagons
    assert P.inertia_matrix().shape == (3, 3)
    lo, hi = P.bounding_box()
    assert lo.shape == (3,) and np.all(hi >= lo)


def _low_symmetry_C60():
    """buckygen's C60 #17 -- the cage where the two mass models disagree most
    (3.70 deg between axes that are individually determined; measured in
    tests/volume-moments-test.cc). buckygen's FIRST C60 is the icosahedral one,
    whose inertia tensors are isotropic and whose frames therefore agree
    trivially -- exactly the wrong fixture for these two tests."""
    with fl.buckygen(60) as gen:
        dual = [d for _, d in zip(range(18), gen)][-1]
    return fl.Polyhedron.from_fullerene(dual.dual())


def _offdiag(I):
    """How far a tensor is from diagonal, relative -- i.e. how far the coordinate
    axes are from being its principal axes."""
    return np.linalg.norm(I - np.diag(np.diag(I))) / np.linalg.norm(I)


def test_polyhedron_mass_model_is_a_named_enum():
    # The mass model is a NAMED enum argument, not a bare bool: MassModel.SOLID
    # (the default) integrates the ENCLOSED SOLID, MassModel.ATOMS puts uniform
    # mass at the vertices -- the molecular convention. Accepted BY NAME.
    P = _low_symmetry_C60()
    solid = P.inertia_matrix(mass_model=fl.MassModel.SOLID)
    atoms = P.inertia_matrix(mass_model=fl.MassModel.ATOMS)
    assert solid.shape == atoms.shape == (3, 3)
    assert np.allclose(P.inertia_matrix(), solid)          # SOLID is the default
    assert np.allclose(P.principal_axes(), P.principal_axes(mass_model=fl.MassModel.SOLID))

    # Two mass models, two different answers. Compared unit-trace normalised, so
    # only the SHAPE of the mass distribution is compared and the differing total
    # masses drop out: the two agree to percent, not to roundoff.
    shape = lambda I: I / np.trace(I)
    d = np.linalg.norm(shape(atoms) - shape(solid)) / np.linalg.norm(shape(solid))
    assert 1e-3 < d < 0.10, f"the two mass models differ by {d}"


def test_polyhedron_from_fullerene_leaves_the_cage_in_the_ATOMS_frame():
    # from_fullerene() centres the cage and aligns it with MassModel.ATOMS -- the
    # molecular convention. A DEFAULTED align_with_axes() is MassModel.SOLID, so
    # it silently rotates the cage back OUT of that frame (by up to 3.70 deg on a
    # low-symmetry cage). Documented on both Python methods; pinned here.
    P = _low_symmetry_C60()
    assert _offdiag(P.inertia_matrix(mass_model=fl.MassModel.ATOMS)) < 1e-9   # in it
    assert _offdiag(P.inertia_matrix()) > 1e-4                                # not the SOLID one

    P.align_with_axes()                                      # defaulted -> MassModel.SOLID
    assert _offdiag(P.inertia_matrix()) < 1e-9                               # now the SOLID one
    assert _offdiag(P.inertia_matrix(mass_model=fl.MassModel.ATOMS)) > 1e-4  # left the molecular one


def test_polyhedron_move_to_origin_writes_through():
    P = fl.Polyhedron.C20()
    pts = P.points
    P.move_to_origin()
    assert np.allclose(pts.mean(axis=0), 0, atol=1e-9)  # same buffer, mutated


def test_polyhedron_byo_optimize_writes_back():
    P0 = fl.Polyhedron.C20()
    nb = np.array(P0.neighbours, dtype=np.int32)
    pts = np.array(P0.points, dtype=np.float64)
    pts += np.random.default_rng(0).normal(0, 0.05, pts.shape)
    P = fl.Polyhedron.from_arrays(nb, pts)
    ptr0 = pts.ctypes.data
    before = pts.copy()
    P.optimize()
    assert pts.ctypes.data == ptr0                    # not reallocated
    assert not np.allclose(pts, before)               # optimizer wrote into caller's array


def test_polyhedron_from_arrays_infers_nonuniform_degree():
    # A Polyhedron is general: a fullerene dual is deg 5/6. from_arrays must infer
    # those degrees from the -1 padding, NOT assume cubic all-3 (review finding #1).
    dual = fl.FullereneGraph.C20().dual()           # icosahedron: all 12 vertices deg 5
    nb = np.array(dual.neighbours, dtype=np.int32)
    deg_true = np.asarray(dual.deg)
    pts = np.zeros((nb.shape[0], 3), dtype=np.float64)
    P = fl.Polyhedron.from_arrays(nb, pts)          # deg omitted -> inferred from padding
    assert np.array_equal(np.asarray(P.deg), deg_true)   # 5s
    assert not np.all(np.asarray(P.deg) == 3)            # not silently cubic


def test_polyhedron_transforms_return_new_objects():
    P = fl.Polyhedron.C20()
    assert fl.Polyhedron.from_arrays  # sanity
    D = P.dual()
    assert isinstance(D, fl.Polyhedron) and D is not P


# --- Deltahedron -------------------------------------------------------------

def test_deltahedron_from_dual_is_convex_equilateral():
    D = fl.Deltahedron.from_dual(fl.FullereneGraph.C20().dual())
    assert D.N == 12 and D.n_triangles() == 20
    assert D.triangles().shape == (20, 3)
    assert D.points.base is D and D.points.flags.writeable
    # C20 dual is the icosahedron: exact equilateral, fully convex.
    assert D.count_concave() == 0
    assert D.max_angle_relerr() < 1e-6


def test_deltahedron_optimize_returns_result_and_stats():
    D = fl.Deltahedron.from_dual(fl.FullereneGraph.C20().dual())
    r = D.optimize()
    assert r in (fl.OptResult.CONVERGED, fl.OptResult.STAGNATED)
    assert D.iterations_used >= 0
    assert D.final_angle_relerr < 1e-6


def test_deltahedron_stats_raise_before_a_run():
    # Reading an optimizer stat on geometry that has not been optimized (a fresh
    # from_arrays wrapper, or a transform result) raises -- it does not return a
    # fabricated default verdict. After optimize(), the snapshot is populated.
    D0 = fl.Deltahedron.from_dual(fl.FullereneGraph.C20().dual())
    D = fl.Deltahedron.from_arrays(np.array(D0.neighbours, dtype=np.int32),
                                   np.array(D0.points, dtype=np.float64))
    with pytest.raises(RuntimeError):
        _ = D.final_opt_result
    D.optimize()
    assert D.final_opt_result in (fl.OptResult.CONVERGED, fl.OptResult.STAGNATED)


def test_deltahedron_from_arrays_argorder_and_infer():
    # from_arrays(neighbours, points, deg=None) — same positional order as
    # Polyhedron.from_arrays; deg inferred (5/6) from the -1 padding.
    D0 = fl.Deltahedron.from_dual(fl.FullereneGraph.C20().dual())
    nb = np.array(D0.neighbours, dtype=np.int32)
    pts = np.array(D0.points, dtype=np.float64)
    D = fl.Deltahedron.from_arrays(nb, pts)          # deg omitted -> inferred
    assert D.N == D0.N
    assert np.array_equal(np.asarray(D.deg), np.asarray(D0.deg))   # 5s, inferred
    assert np.shares_memory(pts, D.points)           # zero-copy


def test_deltahedron_gc_transform():
    D = fl.Deltahedron.from_dual(fl.FullereneGraph.C20().dual())
    D2 = D.gc_transform(2, 0)
    assert D2.N == 42                                 # C20 -> GC(2,0) = C80 dual
    assert isinstance(D2, fl.Deltahedron)
