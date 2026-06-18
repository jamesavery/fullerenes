"""Phase 2: the zero-copy ownership model, validated on FullereneGraph.

Proves both modes:
  Owned — C++ produced the object; .neighbours is a zero-copy numpy view whose
          base keeps the wrapper alive.
  View  — caller's numpy buffer is shared in place (np.shares_memory).
And that real library computations (is_a_fullerene, canonical name) run through
the transient view in either mode.
"""
import gc

import numpy as np
import pytest

import fullerenes as fl


def test_owned_mode_is_zero_copy_view():
    fg = fl.FullereneGraph.C20()
    assert fg.N == 20 and fg.dmax == 3
    nb = fg.neighbours
    assert nb.shape == (20, 3) and nb.dtype == np.int32
    assert nb.base is fg          # keepalive: array holds the wrapper alive
    assert not nb.flags.owndata   # not a copy
    assert not nb.flags.writeable  # topology is read-only


def test_owned_view_survives_owner_gc():
    arr = fl.FullereneGraph.C20().neighbours  # wrapper has no other strong ref
    gc.collect()
    assert int(arr[0, 0]) >= 0    # base keepalive kept storage alive


def test_real_methods_run_through_the_view():
    fg = fl.FullereneGraph.C20()
    assert fg.is_a_fullerene()
    assert fg.name() == "C20-[GS:1,2,3,4,5,6,7,8,9,10,11,12]-fullerene"
    assert fg.adjacency()[0] == [1, 4, 7]


def test_view_mode_shares_caller_buffer():
    fg = fl.FullereneGraph.C20()
    src = np.array(fg.neighbours, dtype=np.int32)        # writable copy
    deg = np.full(src.shape[0], 3, dtype=np.uint8)
    fg2 = fl.FullereneGraph.from_arrays(src, deg)
    out = fg2.neighbours
    assert np.shares_memory(src, out)                    # genuinely zero-copy
    assert out.base is fg2
    assert fg2.is_a_fullerene()
    assert fg2.name() == fg.name()                       # same isomer


def test_view_mode_rejects_readonly_buffer():
    fg = fl.FullereneGraph.C20()
    ro = np.asarray(fg.neighbours)                       # read-only owned view
    deg = np.full(ro.shape[0], 3, dtype=np.uint8)
    with pytest.raises(ValueError):
        fl.FullereneGraph.from_arrays(ro, deg)
