"""Regression tests for the code-review findings."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

import fullerenes as fl


# --- #1: from_rspi validates length and bounds (no OOB heap write) ----------

def test_from_rspi_rejects_wrong_length():
    with pytest.raises(ValueError):
        fl.FullereneDual.from_rspi(60, [1, 2, 3])             # not 12


def test_from_rspi_rejects_out_of_range():
    with pytest.raises(ValueError):
        fl.FullereneDual.from_rspi(60, list(range(1, 12)) + [10_000])
    with pytest.raises(ValueError):
        fl.FullereneDual.from_rspi(60, [0] + list(range(2, 13)))  # 0 < 1-based min


# --- #2: rspi is 1-based, matching name() and IsomerEntry.rspi --------------

def test_rspi_is_one_based():
    idx, _ = fl.FullereneGraph.C20().dual().rspi()
    assert sorted(idx) == list(range(1, 13))                 # C20: all 12 positions


def test_rspi_round_trips_one_based():
    d0 = next(iter(fl.buckygen(60)))
    idx, jumps = d0.rspi()
    assert min(idx) >= 1
    d1 = fl.FullereneDual.from_rspi(60, idx, jumps or None)
    assert d1.name() == d0.name()


# --- #3: from_arrays never silently copies on dtype mismatch ----------------

def test_from_arrays_rejects_wrong_dtype():
    fg = fl.FullereneGraph.C20()
    nb64 = np.array(fg.neighbours, dtype=np.int64)           # default int on Linux
    with pytest.raises(ValueError):
        fl.FullereneGraph.from_arrays(nb64)


def test_from_arrays_zero_copy_with_int32():
    fg = fl.FullereneGraph.C20()
    nb = np.array(fg.neighbours, dtype=np.int32)
    fg2 = fl.FullereneGraph.from_arrays(nb)
    assert np.shares_memory(nb, fg2.neighbours)              # genuinely aliased


# --- #5 / #8: from_arrays validates degree bounds and deg shape -------------

def test_from_arrays_rejects_degree_over_dmax():
    nb = np.zeros((4, 3), dtype=np.int32)
    deg = np.array([3, 3, 3, 7], dtype=np.uint8)             # 7 > dmax 3
    with pytest.raises(ValueError):
        fl.FullereneGraph.from_arrays(nb, deg)


def test_from_arrays_rejects_bad_deg_ndim():
    nb = np.zeros((4, 3), dtype=np.int32)
    deg2d = np.zeros((4, 2), dtype=np.uint8)
    with pytest.raises(ValueError):
        fl.FullereneGraph.from_arrays(nb, deg2d)


# --- #12: generator return-type parser handles tricky signatures ------------

def _load_generator():
    path = Path(fl.__file__).resolve().parent.parent / "tools" / "gen_bindings.py"
    spec = importlib.util.spec_from_file_location("gen_bindings", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod      # needed for dataclass string-annotation resolution
    spec.loader.exec_module(mod)
    return mod


def test_parse_signature_handles_tricky_returns():
    g = _load_generator()
    assert g.parse_signature("int () const") == ("int", True)
    assert g.parse_signature("bool ()") == ("bool", False)
    assert g.parse_signature("coord3d &() const") == ("coord3d &", True)
    assert g.parse_signature("node_t *() const") == ("node_t *", True)
    assert g.parse_signature("std::function<void (int)> () const") == \
        ("std::function<void (int)>", True)
    assert g.parse_signature("matrix<int> (const vector<node_t> &) const")[0] == "matrix<int>"
