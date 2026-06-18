# pybind — Python bindings for the fullerene library

Zero-copy Python bindings (pybind11) over `libfullerenes.so`: topology, geometry,
optimization, enumeration and canonical naming of fullerene isomers. Built on the
library's **view + owned** architecture, so numpy buffers and the C++ algorithms
share one allocation — no marshalling on the hot path.

## Status

Complete and in-tree. The binding builds as an optional CMake target
(`-DENABLE_PYTHON=ON`); the surface is **9 classes + 6 module functions**, covered
by the pytest suite and two fail-closed drift gates (see
`doc/pybind-interface.tex` for the full write-up). GPU/SYCL batch is deferred by
design.

## Install (plain pip, recommended)

Build from source against a prebuilt `libfullerenes.so` — no conda required.
Prerequisites: a C++ compiler matching the library (clang++ by default; see ABI
below), the Python development headers (`python3-dev`), and the library built once
at `<repo>/build`.

```bash
python3 -m venv .venv && . .venv/bin/activate     # from the repo root
pip install -e src/pybind                         # build isolation fetches pybind11 + scikit-build-core
python -c "import fullerenes as f; print(f.version())"
```

`numpy` comes as a dependency; `pybind11` and `scikit-build-core` are pulled into
pip's isolated build env automatically. To run the test suite too:

```bash
pip install -e "src/pybind[test]"
pytest src/pybind/tests
```

For a faster dev loop, pre-install the build deps
(`pip install "pybind11>=2.13" scikit-build-core`) and add `--no-build-isolation`.
The editable install rebuilds the extension on import, so C++ edits take effect on
the next `import`.

> ABI: the module must match the library's compiler/stdlib. The build pins
> `clang++` (what `libfullerenes.so` is built with), and system Python's
> `libstdc++` is the very one the library links — so plain system Python is
> ABI-safe with no conda. If your library was built with gcc, override:
> `pip install -e src/pybind --config-settings=cmake.define.CMAKE_CXX_COMPILER=g++-13`.
> A non-editable `pip install` from an sdist is not supported yet (the wheel does
> not bundle `libfullerenes.so` + its native deps — deferred).

## Build with the library (in-tree CMake)

To build the bindings as part of a full library build instead of via pip:

```bash
cmake -S <repo> -B <repo>/build -DENABLE_PYTHON=ON \
  -DPython_EXECUTABLE=$(which python3) -Dpybind11_DIR=$(python3 -m pybind11 --cmakedir)
cmake --build <repo>/build --target _fullerenes   # drops the .so into src/pybind/fullerenes/
ctest --test-dir <repo>/build -R pybind_pytest
```

This imports from the source tree; don't combine it with an editable install (the
editable finder shadows the source-tree `.so`). Conda works for either path (its
`libstdc++` is a superset of the library's) — just run with the conda interpreter
active.

## Design contract

- **Two ownership modes, zero copies.** Each Python object wraps either a C++
  `Owned<View>` (C++ produced it; `.neighbours`/`.points` are zero-copy numpy views
  with the wrapper as `base`) or caller-supplied numpy arrays (bring-your-own-buffer;
  in-place algorithms like `optimize()` write straight back into them).
- **Reallocation invariant.** A wrapper's backing storage never reallocates in place.
  Operations that change N/dmax return a *new* object, so outstanding numpy views
  never dangle.
- **`.neighbours` is the raw padded `(N, dmax)` int32 view** (columns `>= deg[u]` are
  `-1`); `.adjacency()` returns a clean ragged copy.
- **The method surface is generated** from the clang AST (`tools/gen_bindings.py`),
  not hand-typed, so it tracks the C++ API as it changes; the zero-copy core is
  hand-written and stable. Two fail-closed drift gates keep the committed bindings
  in sync with the C++ AST and the built module (see `doc/pybind-interface.tex`).

## Database

`number_isomers` / `symmetries` / `buckygen` need no files. `IsomerDB.read/get/make`
need the database. At import the package honours the library's compiled-in
`FULLERENE_DATABASE_PATH` (queryable via `get_database_path()`) when that directory
exists; only if it is absent does it fall back to the conventional local copy at
`~/fullerene-data/database`. Override at runtime with
`fullerenes.set_database_path(path)` (a setter — never an environment variable).

## Layout

- `src/`        the extension: `module.cc`, hand-written core (`common.hh`,
  `np_interop.hh`, `naming.hh`), generated `bindings_generated.cc`. (The type-rule
  table lives in `tools/gen_bindings.py`.)
- `tools/`      `gen_bindings.py` (clang-AST → bindings + stubs + report), `spec.toml`,
  `binding_surface.hh`.
- `fullerenes/` the importable package (`__init__.py`, built `_fullerenes*.so`, `.pyi`).
- `tests/`      pytest suites.
- `CMakeLists.txt` (dual-mode: in-tree subdirectory or standalone), `pyproject.toml`.
