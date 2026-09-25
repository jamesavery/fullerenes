"""Phase 3: canonical naming and spiral round-trips.

The canonical name is invariant under representation: graph and dual agree, the
cubic<->dual round-trip preserves it, and rspi() feeds back into from_rspi().
"""
import fullerenes as fl


def test_graph_and_dual_agree():
    fg = fl.FullereneGraph.C20()
    assert fg.name() == fg.dual().name()
    assert fg.dual().dual().name() == fg.name()      # cubic -> dual -> cubic


def test_rspi_round_trips_through_from_rspi():
    for dual in [next(iter(fl.buckygen(60))), next(iter(fl.buckygen(70)))]:
        N = 2 * dual.N - 4                            # carbon count
        idx, jumps = dual.rspi()
        rebuilt = fl.FullereneDual.from_rspi(N, idx, jumps or None)
        assert rebuilt.name() == dual.name()


def test_name_prefix_and_format():
    name = next(iter(fl.buckygen(60))).name()
    # The fullerene (cubic graph) is named: no "T"; a jump-free spiral has no
    # search tag.
    assert name.startswith("C60-[") and name.endswith("]-fullerene")
    assert "T" not in name[4:] and "GS:" not in name


def test_enumerated_names_are_distinct():
    names = [d.name() for d in fl.buckygen(40)]
    assert len(names) == fl.number_isomers(40) == len(set(names))
