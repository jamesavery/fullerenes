"""Phase 2: the AST binding generator.

Verifies the committed generated bindings are not stale w.r.t. the live clang
AST + libfullerenes.so (the fail-closed drift gate), that the methods it bound
are present and callable, and that it never binds a method missing from the
shared library (the dead-declaration guard).
"""
import subprocess
import sys
from pathlib import Path

import pytest

import fullerenes as fl

PROJ = Path(__file__).resolve().parent.parent
GEN = PROJ / "tools" / "gen_bindings.py"


def test_committed_bindings_are_not_stale():
    # Fail-closed drift gate: the committed compiled artifacts (bindings_generated.cc,
    # _fullerenes.pyi) must equal a fresh regeneration from the current clang AST +
    # spec.toml. If they drift (an undeclared C++ API change, or an un-regenerated
    # spec/generator edit), `check` exits non-zero -- this test fails before stale
    # bindings can land. `check` does NOT mutate the files.
    if subprocess.run(["which", "clang++"], capture_output=True).returncode != 0:
        pytest.skip("clang++ not available")
    r = subprocess.run([sys.executable, str(GEN), "check"], capture_output=True, text=True)
    assert r.returncode == 0, (
        "committed bindings are stale; regenerate with: python3 tools/gen_bindings.py\n"
        + r.stderr)
    assert (PROJ / "src" / "bindings_generated.cc").exists()
    assert (PROJ / "fullerenes" / "_fullerenes.pyi").exists()


def test_generated_methods_are_callable():
    fg = fl.FullereneGraph.C20()
    assert fg.is_cubic() is True
    assert fg.is_triangulation() is False
    assert fg.max_degree() == 3
    assert fg.count_edges() == 30          # C20 has 30 edges
    assert fg.is_consistently_oriented() is True
    assert fg.pentagon_distance_mtx().shape == (12, 12)


def test_dead_declaration_is_not_bound():
    # GraphView::hamiltonian_cycle_count() is declared but not defined in the
    # .so (its definition is staged in claude-projects/unfortran until
    # promotion); the generator must drop it, else the module would fail to
    # import.  When the definition lands, this test inverts: the method must
    # then BE bound.
    assert not hasattr(fl.FullereneGraph.C20(), "hamiltonian_cycle_count")
