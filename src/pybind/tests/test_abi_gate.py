"""Phase 1 ABI gate.

Proves a clang++-23 / gnu++23 extension loaded from conda Python 3.13 is
ABI-compatible with libfullerenes.so across the std::string / std::vector
boundary. These checks use only the file-free hard-coded tables, so they pass
with no isomer database installed.
"""
import fullerenes as fl


def test_version_is_a_string():
    v = fl.version()
    assert isinstance(v, str) and v


def test_isomer_counts_match_known_values():
    # General isomer counts (CLAUDE.md reference table).
    assert fl.number_isomers(20) == 1
    assert fl.number_isomers(60) == 1812
    assert fl.number_isomers(80) == 31924
    assert fl.number_isomers(100) == 285914


def test_ipr_isomer_counts():
    assert fl.number_isomers(60, "Any", True) == 1
    assert fl.number_isomers(70, "Any", True) == 1
    assert fl.number_isomers(100, "Any", True) == 450


def test_symmetries_cross_the_so_boundary():
    # vector<string> built inside libfullerenes.so, marshalled to list[str].
    syms = fl.symmetries(60)
    assert isinstance(syms, list)
    assert syms and all(isinstance(s, str) for s in syms)
    # C60 includes the icosahedral buckminsterfullerene.
    assert any("Ih" in s for s in syms)
