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

## Build

> ABI: the module **must** be built with the same compiler/stdlib as the library
> (clang++ 23, `-std=gnu++23`); as an in-tree target it inherits the parent
> toolchain automatically. The system pybind11 (2.9.1) is too old for Python 3.13 —
> use the pip `pybind11` (>=2.13): `-Dpybind11_DIR=$(python3 -m pybind11 --cmakedir)`.

In-tree, with the rest of the library:

```bash
cmake -S <repo> -B <repo>/build -DENABLE_PYTHON=ON \
  -DPython_EXECUTABLE=$(which python3) -Dpybind11_DIR=$(python3 -m pybind11 --cmakedir)
cmake --build <repo>/build --target _fullerenes
ctest --test-dir <repo>/build -R pybind_pytest
```

Standalone dev loop (this directory as its own project; links a prebuilt
`../../build/libfullerenes.so`):

```bash
cmake -S . -B build-py -DPython_EXECUTABLE=$(which python3) \
  -Dpybind11_DIR=$(python3 -m pybind11 --cmakedir)
cmake --build build-py -j4        # drops _fullerenes*.so into ./fullerenes/
python3 -m pytest tests -q
```

Tier-2 editable install (so any program can `import fullerenes`):

```bash
pip install -e . --no-build-isolation     # needs the library prebuilt at ../../build
```

A non-editable `pip install .` (from an sdist) is not supported yet: the wheel
does not bundle `libfullerenes.so` and its native deps (deferred).

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
