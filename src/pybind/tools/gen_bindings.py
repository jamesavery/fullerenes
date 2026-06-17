#!/usr/bin/env python3
"""gen_bindings.py — generate the mechanical pybind11 surface from the clang AST.

The fullerene C++ API changes often, so the binding *method surface* is generated
rather than hand-typed. This walks the clang JSON AST of tools/binding_surface.hh
(never a text parser -- the AST is the single source of truth), harvests the
public methods of the view chain named in tools/spec.toml, classifies each by its
return/parameter types via a fixed rule table, and emits:

  src/bindings_generated.cc      register_generated_<PyClass>(cls) per class
  fullerenes/_fullerenes.pyi     type stubs (same AST pass)
  unbound_report.md              every method it could NOT bind, with the reason
                                 (fail-loud: nothing is silently dropped)

Follows the in-project precedent buckygen-revival/tools/extract_bucky_c.py.

Usage:
    gen_bindings.py [emit]    [--spec tools/spec.toml] [--clang clang++]  # write the files
    gen_bindings.py check     # fail (exit 1) if the committed .cc/.pyi drifted (no writes)
    gen_bindings.py --inspect FullereneGraphView                          # dump AST methods
"""

from __future__ import annotations

import argparse
import difflib
import json
import os
import re
import subprocess
import sys
try:
    import tomllib                      # Python 3.11+
except ModuleNotFoundError:             # 3.9/3.10: pip install tomli (see the test extra)
    import tomli as tomllib
from dataclasses import dataclass, field

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
FULL_ROOT = os.path.abspath(os.path.join(PROJ, "..", ".."))
INCLUDES = [
    os.path.join(FULL_ROOT, "include"),
    os.path.join(FULL_ROOT, "build", "include"),
]
# The library whose symbols the dead-declaration gate checks against (see
# load_defined_nullary_keys). Defaults to the in-tree build, but is overridable
# via --so: an out-of-tree build tree (or a promoted location whose build dir is
# not <root>/build) needs to point the gate at the actual freshly-built .so.
DEFAULT_SO = os.path.join(FULL_ROOT, "build", "src", "c++", "libfullerenes.so")

# --- Return-type rule table -------------------------------------------------
# Each entry maps a normalized C++ return type to (call_wrapper, py_stub_type).
# call_wrapper is a format string with {expr} = the C++ call expression.
SCALAR_TYPES = {
    "bool": "bool", "int": "int", "unsigned int": "int", "unsigned": "int",
    "long": "int", "long long": "int", "size_t": "int", "unsigned long": "int",
    "int64_t": "int", "uint64_t": "int", "node_t": "int",
    "double": "float", "float": "float",
}
STRING_TYPES = {"std::string", "string", "std::__cxx11::basic_string<char>"}

# Owned C++ types: returning these needs a Py* wrapper mapping (cross-class).
# Deferred for now -> reported, not bound.
OWNED_TYPES = {
    "FullereneGraph", "FullereneDual", "Triangulation", "PlanarGraph",
    "Polyhedron", "Deltahedron", "Graph", "CubicGraph",
}


@dataclass
class Method:
    name: str
    ret: str            # normalized return type
    qual: str           # full function qualType
    nparams: int
    const: bool
    static: bool
    owner: str          # which harvested class declared it
    has_body: bool      # defined inline in a header (no .so symbol needed)
    hidden: bool = False  # name hidden by a more-derived class (C++ name hiding)


@dataclass
class Bound:
    name: str
    cc_expr: str        # the .def body expression
    pystub: str         # python return annotation


@dataclass
class Unbound:
    cls: str
    name: str
    qual: str
    reason: str


# --- AST utilities ----------------------------------------------------------

def dump_ast(clang: str) -> dict:
    surface = os.path.join(HERE, "binding_surface.hh")
    cmd = [clang, "-std=gnu++23", "-Xclang", "-ast-dump=json",
           "-fsyntax-only", "-ferror-limit=0", "-w", "-xc++", surface]
    for inc in INCLUDES:
        cmd += ["-I", inc]
    raw = subprocess.run(cmd, check=True, stdout=subprocess.PIPE).stdout
    return json.loads(raw)


def parse_signature(qual: str) -> tuple[str, bool]:
    """From a clang function qualType 'RET (PARAMS) CV...' return (return_type, is_const).

    Strips a *balanced* trailing parameter list rather than splitting on the first
    ' (' -- the latter breaks for reference/pointer returns ('X &(...)' has no space
    before '(') and for return types that themselves contain ' (' (e.g.
    'std::function<void (int)>'). The cv-qualifier is read from after the param list."""
    s = qual.strip()
    is_const = False
    changed = True
    while changed:                         # peel trailing method cv/ref/exception specs
        changed = False
        for suf in (" const", " volatile", " noexcept", " &&", " &"):
            if s.endswith(suf):
                if suf == " const":
                    is_const = True
                s = s[: -len(suf)].rstrip()
                changed = True
    if not s.endswith(")"):
        return s, is_const                 # no param list (unexpected) -> whole string
    depth = 0
    for k in range(len(s) - 1, -1, -1):
        if s[k] == ")":
            depth += 1
        elif s[k] == "(":
            depth -= 1
            if depth == 0:
                return s[:k].strip(), is_const
    return s, is_const


def load_defined_nullary_keys(so_path: str) -> set[str]:
    """Set of exact 'Class::method' keys for every NULLARY defined symbol in the .so.

    Out-of-line methods (declared in a header, no inline body) must be defined here,
    else binding them yields an undefined-symbol error at import. The key is exact
    (not a substring of the demangled blob), so 'GraphView::f' is never conflated
    with 'CubicGraphView::f', and arity-aware (requires '()'), so a defined multi-arg
    overload does not vouch for an undefined nullary of the same name.

    The .so must be freshly built: a STALE library makes the gate wrong (a newly
    defined method looks dead and is dropped; a removed one looks live and crashes
    import). We can't cheaply prove freshness here, but we fail loudly if it is
    absent or unreadable rather than miskeying silently."""
    if not os.path.exists(so_path):
        sys.exit(f"gen_bindings: {so_path} not found -- build libfullerenes.so first "
                 f"(the dead-declaration gate reads its symbols). It must also be "
                 f"up to date with the headers, or the gate miskeys.")
    r = subprocess.run(["nm", "-DC", "--defined-only", so_path],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if r.returncode != 0:
        sys.exit(f"gen_bindings: nm failed on {so_path}: "
                 f"{r.stderr.decode('utf-8', 'replace').strip()}")
    out = r.stdout.decode("utf-8", "replace")
    keys: set[str] = set()
    for line in out.splitlines():
        parts = line.split(" ", 2)         # ADDR TYPE  demangled-symbol
        if len(parts) < 3:
            continue
        name_part, sep, rest = parts[2].strip().partition("(")
        if not sep or not rest.startswith(")"):    # require '()' -> nullary
            continue
        # Strip ALL template arg lists (innermost-first, handles nesting) BEFORE
        # splitting on '::' -- otherwise a symbol like
        # 'Foo<std::pair<int,int>>::bar' would split inside the template args and
        # mis-key. This makes a templated view's symbol
        # (PolyhedronView<double>::foo) reduce to the AST owner name
        # (PolyhedronView::foo).
        flat = name_part
        while True:
            stripped = re.sub(r"<[^<>]*>", "", flat)
            if stripped == flat:
                break
            flat = stripped
        comps = flat.split("::")
        if len(comps) >= 2:
            keys.add(f"{comps[-2]}::{comps[-1]}")   # 'Class::method'
    return keys


def find_record(root: dict, name: str) -> dict | None:
    """Find the defining CXXRecordDecl for `name` anywhere in the AST."""
    best = None
    stack = [root]
    while stack:
        node = stack.pop()
        if not isinstance(node, dict):
            continue
        if node.get("kind") == "CXXRecordDecl" and node.get("name") == name:
            if node.get("inner"):          # the definition, not a forward decl
                best = node
        for ch in node.get("inner", []) or []:
            stack.append(ch)
    return best


def methods_of(record: dict) -> list[Method]:
    """Public non-special methods declared directly on `record`."""
    out = []
    tag = record.get("tagUsed", "class")
    access = "public" if tag == "struct" else "private"
    for ch in record.get("inner", []) or []:
        if not isinstance(ch, dict):
            continue
        kind = ch.get("kind")
        if kind == "AccessSpecDecl":
            access = ch.get("access", access)
            continue
        cur = ch.get("access", access)        # explicit per-decl access wins
        if kind != "CXXMethodDecl":
            continue
        if cur != "public":
            continue
        name = ch.get("name", "")
        if not name or name.startswith("operator") or name.startswith("~"):
            continue
        qual = ch.get("type", {}).get("qualType", "")
        ret, is_const = parse_signature(qual)
        inner = ch.get("inner", []) or []
        nparams = sum(1 for g in inner
                      if isinstance(g, dict) and g.get("kind") == "ParmVarDecl")
        has_body = any(isinstance(g, dict) and g.get("kind") == "CompoundStmt"
                       for g in inner)
        out.append(Method(
            name=name, ret=ret, qual=qual, nparams=nparams,
            const=is_const,
            static=ch.get("storageClass") == "static",
            owner=record.get("name", ""), has_body=has_body))
    return out


def harvest(root: dict, harvest_classes: list[str]) -> list[Method]:
    """Methods across the chain, leaf-first. A derived class hides a base-class
    method of the same name (C++ name hiding), but ALL overloads declared in the
    hiding class are kept (so a nullary overload isn't lost to a multi-arg one)."""
    claimed, out = set(), []
    for cname in harvest_classes:
        rec = find_record(root, cname)
        if rec is None:
            print(f"  warning: class {cname} not found in AST", file=sys.stderr)
            continue
        ms = methods_of(rec)
        for m in ms:
            if m.name in claimed:        # hidden by a more-derived class
                m.hidden = True          # keep (don't silently drop) so it can be reported
            out.append(m)
        claimed |= {m.name for m in ms}
    return out


# --- Classification ---------------------------------------------------------

def classify(m: Method, skip: set[str], defined_keys: set[str]) -> tuple[Bound | None, str | None]:
    """Return (Bound, None) if bindable, else (None, reason)."""
    if m.name in skip:
        return None, "skip (provided by hand-written core)"
    if m.static:
        return None, "static method (bind by hand if needed)"
    if not m.const:
        return None, "non-const (mutator; bind by hand under the realloc rule)"
    if m.nparams != 0:
        return None, "has parameters (deferred to a later generator pass)"
    # Out-of-line nullary methods must actually be defined in the .so, else the
    # module fails to load with an undefined symbol. Inline (header-body) methods
    # are compiled into the binding TU and are always safe. Exact 'Class::method'
    # match (arity-aware) -- substring matching would conflate GraphView with
    # CubicGraphView and a nullary with a multi-arg overload.
    if not m.has_body and f"{m.owner}::{m.name}" not in defined_keys:
        return None, "declared but not defined in libfullerenes.so (dead declaration)"

    ret = m.ret.removeprefix("const ").strip()
    if ret.endswith("&") or ret.endswith("*"):
        return None, "reference/pointer return (lifetime; bind by hand)"
    call = f"w.view().{m.name}()"

    if ret in SCALAR_TYPES:
        return Bound(m.name, call, SCALAR_TYPES[ret]), None
    if ret in STRING_TYPES:
        return Bound(m.name, call, "str"), None
    if ret in ("vector<tri_t>", "std::vector<tri_t>"):
        return Bound(m.name, f"pyf::tris_copy({call})", "numpy.ndarray"), None
    if ret in ("matrix<int>", "matrix<double>"):
        elem = ret[ret.index("<") + 1: ret.index(">")]
        return Bound(m.name, f"pyf::matrix_copy<{elem}>({call})", "numpy.ndarray"), None
    if ret == "matrix3d":
        return Bound(m.name, f"pyf::matrix3d_copy({call})", "numpy.ndarray"), None
    if ret in ("coord3<T>", "coord3<double>", "coord3d"):
        return Bound(m.name, f"pyf::coord3d_copy({call})", "numpy.ndarray"), None
    if ret in ("vector<face_t>", "std::vector<face_t>"):
        return Bound(m.name, f"pyf::faces_copy({call})", "list[list[int]]"), None
    if ret in OWNED_TYPES:
        return None, f"owned return {ret} (needs Py wrapper mapping; deferred)"
    return None, f"unsupported return type '{m.ret}'"


# --- Emission ---------------------------------------------------------------

CC_HEADER = """// bindings_generated.cc -- GENERATED by tools/gen_bindings.py. DO NOT EDIT.
// Regenerate with: python3 tools/gen_bindings.py
//
// One register_generated_<PyClass>(cls) per spec'd class, adding the mechanical
// method surface harvested from the clang AST. The hand-written per-class .cc
// creates the py::class_, calls bind_graph_common, then calls these.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "common.hh"
#include "np_interop.hh"

#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/deltahedron.hh"

namespace py = pybind11;
"""


def wrapper_template(k: dict) -> str:
    return "PyGeom" if k.get("wrapper") == "geom" else "PyGraph"


def emit_cc(klasses: list[dict], bound: dict[str, list[Bound]]) -> str:
    out = [CC_HEADER]
    for k in klasses:
        py, owned, view = k["pyname"], k["owned"], k["view"]
        wt = wrapper_template(k)
        out.append(f"\nvoid register_generated_{py}(py::class_<pyf::{wt}<{owned}, {view}>>& cls) {{")
        out.append(f"    using W = pyf::{wt}<{owned}, {view}>;")
        for b in bound[py]:
            out.append(f'    cls.def("{b.name}", [](W& w) {{ return {b.cc_expr}; }});')
        out.append("}")
    out.append("")
    return "\n".join(out)


# Accessors shared by every wrapper (mirrors bind_adjacency_accessors), plus the
# per-flavour extra: is_a_fullerene for graph wrappers, a points property for geom.
PYI_ADJ_CORE = [
    "    N: int",
    "    dmax: int",
    "    @property",
    "    def neighbours(self) -> numpy.ndarray: ...",
    "    @property",
    "    def deg(self) -> numpy.ndarray: ...",
    "    def adjacency(self) -> list[list[int]]: ...",
]
PYI_GRAPH_EXTRA = ["    def is_a_fullerene(self) -> bool: ..."]
PYI_GEOM_EXTRA = ["    @property", "    def points(self) -> numpy.ndarray: ..."]


def pyi_core(klass: dict) -> list:
    extra = PYI_GEOM_EXTRA if klass.get("wrapper") == "geom" else PYI_GRAPH_EXTRA
    return PYI_ADJ_CORE + extra


def emit_pyi(klasses: list[dict], bound: dict[str, list[Bound]], spec: dict) -> str:
    out = ["# _fullerenes.pyi -- GENERATED by tools/gen_bindings.py. DO NOT EDIT.",
           "import numpy", ""]
    # Hand-written stubs (spec data) for classes the generator does NOT harvest +
    # the module-level functions -- emitted with the same indent-and-join path as
    # the per-klass extra_stub below, so the hand-written class/function stub text
    # lives in spec.toml (the harvested-class shared core is still PYI_ADJ_CORE etc.).
    for pk in spec.get("preamble_klass", []):
        out.append(f"class {pk['pyname']}:")
        out += [f"    {line}" for line in pk["stub"]]
        out.append("")
    out += spec.get("meta", {}).get("pyi_module_functions", [])
    out.append("")
    for k in klasses:
        out.append(f"class {k['pyname']}:")
        out += pyi_core(k)
        for b in bound[k["pyname"]]:
            out.append(f"    def {b.name}(self) -> {b.pystub}: ...")
        # Hand-written (non-generated) methods, declared in spec.toml so the stub
        # is complete. Each entry is a full def/decorator line (no indent).
        for line in k.get("extra_stub", []):
            out.append(f"    {line}")
        out.append("")
    return "\n".join(out)


def emit_report(unbound: list[Unbound], n_bound: int) -> str:
    out = ["# Unbound report -- GENERATED by tools/gen_bindings.py",
           "",
           f"{n_bound} methods bound; {len(unbound)} not bound. Each line below is a "
           "public method the generator could not classify -- handle by hand "
           "(special-case in the .cc) or extend the type-rule table.", ""]
    by_reason: dict[str, list[Unbound]] = {}
    for u in unbound:
        by_reason.setdefault(u.reason, []).append(u)
    for reason in sorted(by_reason):
        out.append(f"## {reason}  ({len(by_reason[reason])})")
        for u in sorted(by_reason[reason], key=lambda x: (x.cls, x.name)):
            out.append(f"- `{u.cls}::{u.name}` -- `{u.qual}`")
        out.append("")
    return "\n".join(out)


# --- Driver -----------------------------------------------------------------

# The three generated artifacts, by path. The two COMPILED ones are the drift
# gate's surface (GATED); the unbound report is emit-only -- it embeds raw,
# un-normalized clang type spellings, so gating it would turn a clang-version
# difference into false drift, whereas the .cc/.pyi go through the rule table.
GEN_CC     = os.path.join(PROJ, "src", "bindings_generated.cc")
GEN_PYI    = os.path.join(PROJ, "fullerenes", "_fullerenes.pyi")
GEN_REPORT = os.path.join(PROJ, "unbound_report.md")
GATED      = (GEN_CC, GEN_PYI)


# The generated artifacts as {path: text}. A PURE function of the clang AST of
# binding_surface.hh and spec.toml -- the committed compiled files are its
# fixpoint, asserted by `check`. Both `emit` (write) and `check` (diff) call it,
# so there is one source of the generated surface.
def generate(spec_path: str, clang: str, so_path: str = DEFAULT_SO):
    root = dump_ast(clang)
    defined_keys = load_defined_nullary_keys(so_path)

    with open(spec_path, "rb") as f:
        spec = tomllib.load(f)
    skip = set(spec.get("meta", {}).get("skip", []))
    klasses = spec["klass"]

    bound: dict[str, list[Bound]] = {}
    unbound: list[Unbound] = []
    n_bound = 0
    summary: list[str] = []
    for k in klasses:
        methods = harvest(root, k["harvest"])
        kskip = skip | set(k.get("skip", []))      # global + per-class skip
        kbound: list[Bound] = []
        for m in methods:
            b, reason = classify(m, kskip, defined_keys)
            if m.hidden:
                # C++ name hiding: not callable via the derived view, so never
                # bind it -- but if it WOULD otherwise be bindable, report it so
                # nothing bindable is silently dropped.
                if b is not None:
                    unbound.append(Unbound(k["pyname"], m.name, m.qual,
                        "hidden by a more-derived overload (C++ name hiding); not callable via the view"))
                continue
            if b is not None:
                kbound.append(b)
            elif not reason.startswith("skip"):
                unbound.append(Unbound(k["pyname"], m.name, m.qual, reason))
        bound[k["pyname"]] = kbound
        n_bound += len(kbound)
        summary.append(f"  {k['pyname']}: {len(kbound)} bound, "
                       f"{sum(1 for u in unbound if u.cls == k['pyname'])} unbound")

    texts = {
        GEN_CC:     emit_cc(klasses, bound),
        GEN_PYI:    emit_pyi(klasses, bound, spec),
        GEN_REPORT: emit_report(unbound, n_bound),
    }
    return texts, n_bound, len(unbound), summary   # summary is a list; caller formats


# Compare the GATED (compiled) artifacts to the committed files; print a unified
# diff + return 1 on any drift, 0 if in sync. Fail-closed: an undeclared C++ API
# change (or a spec/generator edit) leaves the committed bindings stale and the
# pytest suite / `/build check` breaks before they can land. Does not write.
def check(texts: dict[str, str]) -> int:
    drift = []
    for path in GATED:
        # utf-8 + default (universal) newline: a CRLF checkout reads back as \n,
        # so it matches the '\n'-joined generated text rather than false-drifting.
        try:
            with open(path, encoding="utf-8") as f:
                on_disk, missing = f.read(), False
        except FileNotFoundError:
            on_disk, missing = "", True
        if on_disk != texts[path]:
            drift.append((path, on_disk, texts[path], missing))
    for path, old, new, missing in drift:
        rel = os.path.relpath(path, PROJ)
        note = " (missing on disk)" if missing else ""
        sys.stderr.write(f"\nDRIFT in {rel}{note}:\n")
        sys.stderr.writelines(difflib.unified_diff(
            old.splitlines(keepends=True), new.splitlines(keepends=True),
            fromfile=f"committed/{rel}", tofile=f"regenerated/{rel}"))
    if drift:
        sys.stderr.write(f"\nbindings differ from the current C++ AST / spec "
                         f"({len(drift)} file(s)); regenerate with:\n"
                         f"    python3 tools/gen_bindings.py\n")
        return 1
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="emit", choices=["emit", "check"],
                    help="emit: write the generated files; check: fail if they drifted")
    ap.add_argument("--spec", default=os.path.join(HERE, "spec.toml"))
    ap.add_argument("--clang", default="clang++")
    ap.add_argument("--so", default=DEFAULT_SO,
                    help="libfullerenes.so the dead-declaration gate checks against "
                         "(default: the in-tree build at <root>/build/src/c++)")
    ap.add_argument("--inspect", metavar="CLASS",
                    help="print harvested methods of CLASS and exit")
    args = ap.parse_args()

    if args.inspect:
        rec = find_record(dump_ast(args.clang), args.inspect)
        if rec is None:
            print(f"class {args.inspect} not found", file=sys.stderr)
            return 2
        for m in methods_of(rec):
            flags = "".join(c for c, on in
                            [("s", m.static), ("c", m.const)] if on)
            print(f"  {m.ret:<28} {m.name}()  [{m.nparams} params, {flags}]  :: {m.qual}")
        return 0

    texts, n_bound, n_unbound, summary = generate(args.spec, args.clang, args.so)

    if args.mode == "check":
        rc = check(texts)
        if rc == 0:
            print(f"bindings up to date ({n_bound} bound, {n_unbound} unbound)")
        return rc

    for path, text in texts.items():
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
    print("\n".join(summary), file=sys.stderr)
    print(f"generated: {n_bound} methods bound, {n_unbound} reported unbound")
    return 0


if __name__ == "__main__":
    sys.exit(main())
