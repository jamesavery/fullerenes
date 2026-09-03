// _fullerenes — Python bindings for the fullerene C++ library.
//
// Phase 1 (this file, for now): the ABI gate. A minimal module that links
// libfullerenes.so and exercises the std::vector<std::string> / std::string
// boundary, to prove that a clang++-23 / gnu++23 module loaded from the conda
// Python interpreter is ABI-compatible with the shared library before any of
// the real binding surface is built.
//
// Later phases replace the hand-written defs below with:
//   - a hand-written zero-copy core (common.hh / np_interop.hh), and
//   - a generated method surface (src/bindings_generated.cc, emitted by
//     tools/gen_bindings.py from the clang JSON AST).

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>   // std::vector<std::string> <-> list[str]

#include <string>
#include <vector>
#include <cstdint>

#include "fullerenes/isomerdb.hh"
#include "fullerenes/graphview.hh"   // pentagon_error

namespace py = pybind11;

// Per-class registration functions (one per translation unit).
void register_fullerene_graph(py::module_& m);
void register_fullerene_dual(py::module_& m);
void register_isomerdb(py::module_& m);
void register_enumeration(py::module_& m);
void register_polyhedron(py::module_& m);
void register_deltahedron(py::module_& m);
void register_symmetry(py::module_& m);

PYBIND11_MODULE(_fullerenes, m) {
    m.doc() = "Python bindings for the fullerene library.";

    // pentagon_error is a contract violation on the caller's input (e.g.
    // from_arrays handed a graph that is not a fullerene dual) -- surface it
    // as ValueError like the neighbouring input validation, not the
    // RuntimeError that std::logic_error would default to.
    py::register_exception_translator([](std::exception_ptr p) {
        try { if (p) std::rethrow_exception(p); }
        catch (const pentagon_error& e) {
            PyErr_SetString(PyExc_ValueError, e.what());
        }
    });

    // Pure pybind<->Python string marshalling (no .so boundary crossed).
    m.def("version", []() -> std::string {
        return "pyfullerenes 0.0.1 (abi-gate)";
    });

    // File-free hard-coded isomer counts (no database needed). Proves we can
    // link and call a real libfullerenes.so symbol and return a scalar.
    m.def("number_isomers",
          [](int N, const std::string& sym, bool IPR) -> int64_t {
              return IsomerDB::number_isomers(N, sym, IPR);
          },
          py::arg("N"), py::arg("sym") = "Any", py::arg("IPR") = false,
          "Number of fullerene isomers of size N (hard-coded table; no DB I/O).");

    // File-free hard-coded symmetry lists. THIS is the real ABI test: a
    // std::vector<std::string> is constructed inside libfullerenes.so and
    // returned across the .so boundary, then marshalled to list[str].
    m.def("symmetries",
          [](int N, bool IPR) -> std::vector<std::string> {
              return IsomerDB::symmetries(N, IPR);
          },
          py::arg("N"), py::arg("IPR") = false,
          "Point-group symmetries present among size-N isomers.");

    // --- Domain classes (Phase 2+) ---
    register_fullerene_graph(m);
    register_fullerene_dual(m);
    register_isomerdb(m);
    register_enumeration(m);
    register_polyhedron(m);
    register_deltahedron(m);
    register_symmetry(m);
}
