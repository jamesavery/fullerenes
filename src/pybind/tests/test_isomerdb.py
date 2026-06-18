"""Phase 3: IsomerDB.

number_isomers / symmetries are file-free and always run. IsomerDB.read /
make_isomer need the database; those tests skip cleanly if it is absent.
"""
import os

import pytest

import fullerenes as fl


def _db_available():
    root = fl.get_database_path()
    return bool(root) and os.path.isdir(os.path.join(root, "All"))


requires_db = pytest.mark.skipif(not _db_available(), reason="isomer database not installed")


def test_number_isomers_table():
    assert fl.number_isomers(20) == 1
    assert fl.number_isomers(60) == 1812
    assert fl.number_isomers(100) == 285914


@requires_db
def test_read_general_c60():
    db = fl.IsomerDB.read(60)
    assert db.N == 60 and not db.IPR
    assert len(db) == 1812 == db.n_isomers


@requires_db
def test_ipr_c60_is_buckminsterfullerene():
    ipr = fl.IsomerDB.read(60, True)
    assert len(ipr) == 1
    e = ipr[0]
    assert e.group.strip() == "Ih"
    assert len(e.rspi) == 12
    fg = ipr.make_isomer(0)
    assert fg.N == 60 and fg.is_a_fullerene()
    # canonical name's pentagon indices echo the database RSPI
    assert fg.name() == "C60-[GS:1,7,9,11,13,15,18,20,22,24,26,32]-fullerene"


@requires_db
def test_indexing_bounds():
    db = fl.IsomerDB.read(60)
    assert db[-1].group  # negative index works
    with pytest.raises(IndexError):
        _ = db[len(db)]


@requires_db
def test_make_isomer_negative_index_parity():
    # make_isomer must accept negative indices, like __getitem__ (finding #11).
    db = fl.IsomerDB.read(60)
    assert db.make_isomer(-1).name() == db.make_isomer(len(db) - 1).name()
    with pytest.raises(IndexError):
        db.make_isomer(len(db))


@requires_db
def test_entry_rspi_feeds_from_rspi():
    # IsomerEntry.rspi (1-based) feeds straight into from_rspi (finding #2).
    ipr = fl.IsomerDB.read(60, True)
    e = ipr[0]
    fd = fl.FullereneDual.from_rspi(60, list(e.rspi))
    assert fd.name() == ipr.make_isomer(0).name()
