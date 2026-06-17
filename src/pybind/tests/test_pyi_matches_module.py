"""extra-stub-sync gate: for each class the .pyi declares, its members must equal
the built module's public surface -- so a hand-written `extra_stub` can't drift
from the .cc (a stub for a method the .cc no longer binds, or a bound method with
no stub). The AST drift gate (test_generator.py) can't see this: it copies
extra_stub from spec.toml *verbatim* without reading the hand-written .cc.

Two checks:
- per declared class, its members equal the built module's (the per-class gate);
- the .pyi declares exactly the module's public classes AND module-level functions
  (the module->pyi direction), so a wholly-unstubbed class/function fails too.

Caveats: names only -- pybind lambdas expose no typed runtime signature, so arg
types AND defaults in a stub are not machine-checkable here (existence is the
guard); and it compares the committed .pyi against the *built* module, so a stale
.so (a .cc change not rebuilt) makes the gate vacuously agree with a stale .pyi.
"""
import ast
from pathlib import Path

import pytest

import fullerenes as fl

_PYI = Path(fl.__file__).resolve().parent / "_fullerenes.pyi"
_ENUM_INTERNALS = {"name", "value"}        # pybind adds these to every py::enum_


def _parse_pyi():
    """One pass over the .pyi: ({public class -> public member names}, {public
    module-level function names}). Underscore-prefixed names are dropped."""
    classes, funcs = {}, set()
    for node in ast.parse(_PYI.read_text()).body:
        if isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
            names = set()
            for it in node.body:
                if isinstance(it, ast.FunctionDef):
                    names.add(it.name)                              # def / @property / @staticmethod
                elif isinstance(it, ast.AnnAssign) and isinstance(it.target, ast.Name):
                    names.add(it.target.id)                         # N: int, CONVERGED: OptResult
            classes[node.name] = {n for n in names if not n.startswith("_")}
        elif isinstance(node, ast.FunctionDef) and not node.name.startswith("_"):
            funcs.add(node.name)
    return classes, funcs


_PYI_MEMBERS, _PYI_FUNCS = _parse_pyi()
# Floor: a class-less .pyi would make the parametrize empty -> a vacuous green.
assert _PYI_MEMBERS, "no classes parsed from _fullerenes.pyi -- gate would pass vacuously"


def _module_members(cls):
    names = {n for n in dir(cls) if not n.startswith("_")}
    if hasattr(cls, "__members__"):        # a pybind enum -> drop its auto name/value
        names -= _ENUM_INTERNALS
    return names


@pytest.mark.parametrize("class_name", sorted(_PYI_MEMBERS))
def test_pyi_matches_built_module(class_name):
    cls = getattr(fl, class_name)
    declared = _PYI_MEMBERS[class_name]
    actual = _module_members(cls)
    missing = actual - declared            # bound on the class but no .pyi stub
    stale = declared - actual              # in the .pyi but not on the built class
    assert not missing, f"{class_name}: bound but missing from .pyi (add an extra_stub): {sorted(missing)}"
    assert not stale, f"{class_name}: in .pyi but not on the built class (stale stub): {sorted(stale)}"


def test_pyi_declares_every_bound_class_and_function():
    # Module->pyi direction: the stub must list every public class and module-level
    # function the extension binds -- so a wholly-unstubbed class/function fails.
    from fullerenes import _fullerenes as ext
    pub = {n for n in dir(ext) if not n.startswith("_")}
    ext_classes = {n for n in pub if isinstance(getattr(ext, n), type)}
    ext_funcs = {n for n in pub if callable(getattr(ext, n)) and n not in ext_classes}
    # Every public extension name is a class or a function -- a public non-callable
    # (e.g. a module-level constant) would otherwise be invisible to both checks.
    assert pub == ext_classes | ext_funcs, (
        f"public extension names that are neither class nor function (not gated here): "
        f"{sorted(pub - ext_classes - ext_funcs)}")
    classes = set(_PYI_MEMBERS)
    assert ext_classes == classes, (
        f"class drift -- bound-not-stubbed={sorted(ext_classes - classes)}, "
        f"stubbed-not-bound={sorted(classes - ext_classes)}")
    assert ext_funcs == _PYI_FUNCS, (
        f"module-function drift -- bound-not-stubbed={sorted(ext_funcs - _PYI_FUNCS)}, "
        f"stubbed-not-bound={sorted(_PYI_FUNCS - ext_funcs)}")
