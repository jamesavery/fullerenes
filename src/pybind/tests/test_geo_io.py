"""The compact binary .geo format (GEO-FORMAT.md) from Python: Polyhedron.to_geo /
to_geo_batch / from_geo / read_geo_header / verify_geo, GeoOptions and GeoHeader,
file-name options through write/from_file, and the error mapping. The format itself
is tested in tests/geo-io-test.cc; these tests cover the binding."""
import numpy as np
import pytest

import fullerenes as fl


def _fixed(width, **kw):
    return fl.GeoOptions(type=fl.GeoType.FIXED, width=width, **kw)


def _moved(P, factor, shift=(0.0, 0.0, 0.0)):
    P.scale([factor, factor, factor])
    P.move(list(shift))
    return P


def _same_polyhedron(Q, P):
    assert Q.adjacency() == P.adjacency()
    np.testing.assert_array_equal(Q.points, P.points)


# --- Round trips ---------------------------------------------------------------

def test_f64_round_trip_is_exact(tmp_path):
    P = fl.Polyhedron.C20()
    path = str(tmp_path / "c20.geo")
    assert P.to_geo(path)
    _same_polyhedron(fl.Polyhedron.from_geo(path), P)

    H = fl.Polyhedron.read_geo_header(path)
    assert (H.count, H.options.type, H.options.N) == (1, fl.GeoType.F64, 20)
    assert (H.options.deg_min, H.options.deg_bits) == (3, 0)
    assert H.record_offset(0) == 32 and H.record_size() % 8 == 0
    assert H.edge_capacity() == 30
    assert (tmp_path / "c20.geo").stat().st_size == 32 + H.record_size()
    assert fl.Polyhedron.verify_geo(path)
    assert "GeoType.F64" in repr(H)


@pytest.mark.parametrize("width", [2, 5, 12, 16, 30])
@pytest.mark.parametrize("offset", [False, True])
def test_fixed_point_error_is_within_half_a_step(tmp_path, width, offset):
    P = _moved(fl.Polyhedron.C20(), 3.7, (100.25, -50.5, 3.0))
    path = str(tmp_path / "c20.geo")
    assert P.to_geo(path, options=_fixed(width, offset=offset))
    Q = fl.Polyhedron.from_geo(path)
    assert Q.adjacency() == P.adjacency()

    x = P.points
    origin = np.zeros(3)
    if offset:
        lo, hi = x.min(axis=0), x.max(axis=0)
        origin = (lo / 2 + hi / 2).astype(np.float32).astype(np.float64)
    step = np.abs(x - origin).max() / (2 ** (width - 1) - 1)   # the per-record scale, up to f32 rounding
    assert np.abs(Q.points - x).max() <= step / 2 * (1 + 1e-6)


def test_file_names_choose_the_coordinates(tmp_path):
    P = fl.Polyhedron.C20()
    for name, type_, width in [("a.geo", fl.GeoType.F64, 0), ("a.f32.geo", fl.GeoType.F32, 0),
                               ("a.q12.geo", fl.GeoType.FIXED, 12)]:
        path = str(tmp_path / name)
        assert P.write(path)
        H = fl.Polyhedron.read_geo_header(path)
        assert (H.options.type, H.options.width, H.options.graph) == (type_, width, True)
        assert fl.Polyhedron.from_file(path).adjacency() == P.adjacency()
    with pytest.raises(ValueError):
        P.write(str(tmp_path / "a.q31.geo"))


def test_bring_your_own_buffer_polyhedron_writes(tmp_path):
    P = fl.Polyhedron.C20()
    V = fl.Polyhedron.from_arrays(P.neighbours.copy(), P.points.copy(), P.deg.copy())
    path = str(tmp_path / "view.geo")
    assert V.to_geo(path)
    _same_polyhedron(fl.Polyhedron.from_geo(path), P)


# --- Appending, batches, capacities ----------------------------------------------

def test_append_random_access_and_verify(tmp_path):
    path = str(tmp_path / "many.geo")
    opt = fl.GeoOptions(deg_min=3, deg_bits=0)
    Ps = [_moved(fl.Polyhedron.C20(), 1 + i / 10) for i in range(4)]
    for P in Ps:
        assert P.to_geo(path, append=True, options=opt)
    assert fl.Polyhedron.read_geo_header(path).count == 4
    for i in (2, 0, 3, 1):
        _same_polyhedron(fl.Polyhedron.from_geo(path, i), Ps[i])
    with pytest.raises(IndexError):
        fl.Polyhedron.from_geo(path, 4)
    assert fl.Polyhedron.verify_geo(path)

    with open(path, "r+b") as f:        # one flipped mantissa bit in record 2
        H = fl.Polyhedron.read_geo_header(path)
        f.seek(H.record_offset(2) + 3)
        b = f.read(1)
        f.seek(-1, 1)
        f.write(bytes([b[0] ^ 0x10]))
    assert not fl.Polyhedron.verify_geo(path)


def test_append_rules(tmp_path):
    P = fl.Polyhedron.C20()
    path = str(tmp_path / "rules.geo")
    with pytest.raises(ValueError, match="degree range"):
        P.to_geo(path, append=True)                     # a file created by appending must declare it
    assert P.to_geo(path, append=True, options=fl.GeoOptions(deg_min=3, deg_bits=0))
    with pytest.raises(RuntimeError, match="different"):
        P.to_geo(path, append=True, options=_fixed(12, deg_min=3, deg_bits=0))
    with pytest.raises(RuntimeError):                   # 12 vertices in a file of 20
        P.dual().to_geo(path, append=True)
    assert fl.Polyhedron.read_geo_header(path).count == 1


def test_batch_with_a_vertex_count_per_record(tmp_path):
    C20 = fl.Polyhedron.C20()
    ico = C20.dual()                                    # 12 vertices of degree 5
    path = str(tmp_path / "batch.geo")
    assert fl.Polyhedron.to_geo_batch(path, [C20, ico, C20], fl.GeoOptions(record_n=True))
    H = fl.Polyhedron.read_geo_header(path)
    assert (H.count, H.options.N, H.options.record_n) == (3, 20, True)
    assert (H.options.deg_min, H.options.deg_bits) == (3, 2)
    for i, P in enumerate([C20, ico, C20]):
        _same_polyhedron(fl.Polyhedron.from_geo(path, i), P)

    # Without record_n every record must have exactly N vertices; the refused
    # write leaves the existing file as it was.
    before = (tmp_path / "batch.geo").read_bytes()
    with pytest.raises(RuntimeError, match="record_n"):
        fl.Polyhedron.to_geo_batch(path, [C20, ico])
    assert (tmp_path / "batch.geo").read_bytes() == before


def test_coordinates_only_needs_a_graph(tmp_path):
    P = fl.Polyhedron.C20()
    path = str(tmp_path / "coords.geo")
    assert P.to_geo(path, options=_fixed(16, graph=False))
    with pytest.raises(RuntimeError, match="no graph"):
        fl.Polyhedron.from_geo(path)
    for graph in (P, fl.FullereneGraph.C20()):
        Q = fl.Polyhedron.from_geo(path, graph=graph)
        assert Q.adjacency() == P.adjacency()
        assert np.abs(Q.points - P.points).max() < 1e-4
    with pytest.raises(ValueError):                     # 12 vertices, the record has 20
        fl.Polyhedron.from_geo(path, graph=P.dual())
    with pytest.raises(TypeError):
        fl.Polyhedron.from_geo(path, graph="C20")


def test_graph_only_archive_holds_no_coordinates(tmp_path):
    path = str(tmp_path / "graph.geo")
    assert fl.Polyhedron.C20().to_geo(path, options=fl.GeoOptions(type=fl.GeoType.NONE))
    assert fl.Polyhedron.read_geo_header(path).options.type == fl.GeoType.NONE
    with pytest.raises(RuntimeError, match="no coordinates"):
        fl.Polyhedron.from_geo(path)


# --- Options and errors ------------------------------------------------------------

def test_geo_options_fields():
    o = fl.GeoOptions()
    assert (o.type, o.width, o.scale, o.offset, o.graph) == (fl.GeoType.F64, 0, 0.0, False, True)
    assert (o.triangulation, o.record_n, o.N, o.deg_min, o.deg_bits, o.sync) == (False, False, 0, -1, -1, False)
    o.type, o.width, o.scale, o.N = fl.GeoType.FIXED, 10, 0.25, 7
    assert (o.type, o.width, o.scale, o.N) == (fl.GeoType.FIXED, 10, 0.25, 7)
    assert "width=10" in repr(o)
    with pytest.raises(ValueError):
        o.width = 256
    with pytest.raises(ValueError):
        fl.GeoOptions(width=-1)
    with pytest.raises(TypeError):
        fl.GeoOptions(fl.GeoType.F32)                   # keyword-only


def test_errors(tmp_path):
    P = fl.Polyhedron.C20()
    with pytest.raises(FileNotFoundError):
        fl.Polyhedron.from_geo(str(tmp_path / "missing.geo"))
    with pytest.raises(FileNotFoundError):
        P.to_geo(str(tmp_path / "no-such-dir" / "x.geo"))

    fresh = tmp_path / "refused.geo"
    for bad in (_fixed(1), fl.GeoOptions(offset=True), fl.GeoOptions(scale=1.0),
                fl.GeoOptions(type=fl.GeoType.NONE, graph=False)):
        with pytest.raises(ValueError):
            P.to_geo(str(fresh), options=bad)
        assert not fresh.exists(), "a refused write leaves no new file behind"

    with pytest.raises(RuntimeError):                   # C20 has pentagons
        P.to_geo(str(fresh), options=fl.GeoOptions(triangulation=True))

    garbage = tmp_path / "garbage.geo"
    garbage.write_bytes(b"\x30" + bytes(40))            # version 1
    with pytest.raises(RuntimeError, match="version"):
        fl.Polyhedron.from_geo(str(garbage))
    with pytest.raises(RuntimeError):
        fl.Polyhedron.verify_geo(str(garbage))


def test_is_consistently_oriented_takes_a_genus():
    P = fl.Polyhedron.C20()
    assert P.is_consistently_oriented()
    assert P.is_consistently_oriented(0) and not P.is_consistently_oriented(genus=1)
    fg = fl.FullereneGraph.C20()
    assert fg.is_consistently_oriented() and fg.dual().is_consistently_oriented()
    assert fl.Deltahedron.from_dual(fg.dual()).is_consistently_oriented()


# --- DelaunayTriangulation -----------------------------------------------------------

def _nonsimplicial_idt():
    """C60 isomer #1264: its iDT is a 12-cone delta-complex with a parallel edge
    (DCEL.IdtRoundTrip in tests/delaunay-test.cc)."""
    with fl.buckygen(60) as gen:
        dual = [d for _, d in zip(range(1265), gen)][-1]
    return fl.DelaunayTriangulation.compute(dual)


def _file_edges(D):
    """The file's edge k (GEO-FORMAT.md sec. 7.5) as a half-edge of D: rows run
    counter-clockwise from v_out, ccw(h) = twin(prev(h)), and edge k is the k-th
    still-unmatched half-edge of the walk."""
    nxt, vout = D.he_next, D.v_out
    order, slot = [], {}
    for v in range(D.nv):
        h = int(vout[v])
        while True:
            slot[h] = len(order)
            order.append(h)
            h = int(nxt[nxt[h]]) ^ 1
            if h == vout[v]:
                break
    edges, matched = [], [False] * len(order)
    for i, h in enumerate(order):
        if not matched[i]:
            edges.append(h)
            matched[i] = matched[slot[h ^ 1]] = True
    return np.array(edges)


def test_delaunay_triangulation_round_trip(tmp_path):
    D = _nonsimplicial_idt()
    assert D.nv == 12 and D.check_consistency() and not D.is_simplicial()
    edges = _file_edges(D)
    x = np.column_stack([np.arange(12.0), np.arange(12.0) ** 2 / 7, 1 / (np.arange(12.0) + 1)])

    for points, opt in ((None, fl.GeoOptions(type=fl.GeoType.NONE)), (x, fl.GeoOptions())):
        path = str(tmp_path / "idt.geo")
        assert D.to_geo(path, points=points, options=opt)
        assert fl.Polyhedron.read_geo_header(path).options.triangulation

        origin, nxt, stored = fl.DelaunayTriangulation.read_geo_topology(path)
        assert origin.shape == nxt.shape == (2 * len(edges),)
        if points is None:
            assert stored is None
        else:
            np.testing.assert_array_equal(stored, x)

        E = fl.DelaunayTriangulation.from_geo(path, lengths=D.he_length[edges],
                                              orig_degree=D.v_orig_degree)
        assert E.check_consistency() and not E.is_simplicial()
        np.testing.assert_array_equal(E.he_origin, origin)
        np.testing.assert_array_equal(E.he_next, nxt)
        # An isomorphism of the two DCELs that fixes every vertex.
        m = np.empty(E.nh, dtype=int)
        m[0::2], m[1::2] = edges, edges ^ 1
        np.testing.assert_array_equal(E.he_origin, D.he_origin[m])
        np.testing.assert_array_equal(m[E.he_next], D.he_next[m])
        np.testing.assert_array_equal(E.he_length, D.he_length[m])
        np.testing.assert_allclose(E.vertex_angle_sums(), D.vertex_angle_sums(), rtol=0, atol=1e-12)
        np.testing.assert_array_equal(E.v_orig_degree, D.v_orig_degree)

        again = str(tmp_path / "again.geo")
        assert E.to_geo(again, points=points, options=opt)
        assert (tmp_path / "again.geo").read_bytes() == (tmp_path / "idt.geo").read_bytes()

    with pytest.raises(RuntimeError, match="self-loop or a parallel edge"):
        fl.Polyhedron.from_geo(path)                     # the F64 file: not a simple polyhedron


def test_delaunay_triangulation_appends_and_takes_a_constant_degree(tmp_path):
    fg = fl.FullereneGraph.C20()
    D = fl.DelaunayTriangulation.compute(fg.dual())     # the icosahedron: 12 cones of degree 5
    assert D.is_simplicial()
    path = str(tmp_path / "ico.geo")
    opt = fl.GeoOptions(type=fl.GeoType.FIXED, width=10, deg_min=5, deg_bits=0)
    x = np.random.default_rng(1).normal(size=(12, 3))
    for _ in range(3):
        assert D.to_geo(path, points=x, append=True, options=opt)
    assert fl.Polyhedron.read_geo_header(path).count == 3
    assert fl.Polyhedron.verify_geo(path)
    ones = np.ones(len(_file_edges(D)))
    E = fl.DelaunayTriangulation.from_geo(path, 2, lengths=ones, orig_degree=5)
    np.testing.assert_allclose(E.vertex_angle_sums(), 5 * np.pi / 3, atol=1e-12)
    assert list(E.v_orig_degree) == [5] * 12
    with pytest.raises(IndexError):
        fl.DelaunayTriangulation.read_geo_topology(path, 3)


def test_delaunay_triangulation_errors(tmp_path):
    D = fl.DelaunayTriangulation.compute(fl.FullereneGraph.C20().dual())
    path = str(tmp_path / "ico.geo")
    with pytest.raises(ValueError, match="positions"):
        D.to_geo(path)                                   # F64 needs points
    with pytest.raises(ValueError, match="points"):
        D.to_geo(path, points=np.zeros((12, 2)))
    with pytest.raises(ValueError):
        D.to_geo(path, points=np.zeros((11, 3)))
    with pytest.raises(ValueError, match="graph"):
        D.to_geo(path, options=fl.GeoOptions(type=fl.GeoType.NONE, graph=False))
    assert D.to_geo(path, options=fl.GeoOptions(type=fl.GeoType.NONE))

    E = len(_file_edges(D))
    for lengths, degree, message in [(np.ones(E - 1), 5, "entries for"),
                                     (np.ones(E), [5] * 11, "entries for"),
                                     (np.r_[0.0, np.ones(E - 1)], 5, "length"),
                                     (np.r_[1e9, np.ones(E - 1)], 5, "check_consistency"),
                                     (np.ones(E), -1, "original degree"),
                                     (np.ones((E, 1)), 5, "one-dimensional")]:
        with pytest.raises(ValueError, match=message):
            fl.DelaunayTriangulation.from_geo(path, lengths=lengths, orig_degree=degree)
    with pytest.raises(TypeError):
        fl.DelaunayTriangulation.from_geo(path, np.ones(E), 5)   # the metric is keyword-only

    square = str(tmp_path / "c20.geo")                   # pentagons: not flagged as a triangulation
    assert fl.Polyhedron.C20().to_geo(square)
    with pytest.raises(RuntimeError, match="triangulation"):
        fl.DelaunayTriangulation.from_geo(square, lengths=np.ones(30), orig_degree=3)
    with pytest.raises(FileNotFoundError):
        fl.DelaunayTriangulation.read_geo_topology(str(tmp_path / "missing.geo"))
