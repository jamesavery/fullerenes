"""Phase 3: buckygen enumeration.

The headline cross-check: the number of isomers buckygen yields equals the
hard-coded number_isomers() table. Also covers the generator's iterator and
context-manager semantics.
"""
import fullerenes as fl


def test_counts_match_number_isomers():
    for N in (20, 24, 28, 60):
        assert sum(1 for _ in fl.buckygen(N)) == fl.number_isomers(N)


def test_ipr_counts():
    assert sum(1 for _ in fl.buckygen(70, IPR=True)) == fl.number_isomers(70, "Any", True)
    assert sum(1 for _ in fl.buckygen(80, IPR=True)) == fl.number_isomers(80, "Any", True)


def test_yields_named_duals():
    g = fl.buckygen(60)
    first = next(iter(g))
    assert type(first).__name__ == "FullereneDual"
    assert first.N == 32                      # C60 dual has N/2+2 = 32 vertices
    assert first.name().startswith("C60-[GS:")


def test_context_manager_and_close():
    with fl.buckygen(80) as gen:
        a, b = next(gen), next(gen)
        assert a.N == 42 and b.N == 42
    # early close: a fresh generator can be closed without exhausting it
    g = fl.buckygen(100)
    next(iter(g))
    g.close()
