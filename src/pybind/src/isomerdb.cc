// isomerdb.cc — Python bindings for the isomer database.
//
// number_isomers() / symmetries() (file-free hard-coded tables) live in
// module.cc. Here we bind set/get_database_path and the IsomerDB reader, which
// DO need the database files (default ~/fullerene-data/database, set by
// fullerenes/__init__.py via set_database_path -- a setter, never an env var).

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstring>
#include <string>
#include <vector>

#include "common.hh"

#include "fullerenes/isomerdb.hh"
#include "fullerenes/fullerenegraph.hh"

namespace py = pybind11;

using PyFG = pyf::PyGraph<FullereneGraph, FullereneGraphView>;

static std::string group_str(const IsomerDB::Entry& e) {
    return std::string(e.group, strnlen(e.group, sizeof(e.group)));
}

// Normalize a Python index (negative = from the end) and bounds-check it.
static size_t resolve_index(const IsomerDB& db, py::ssize_t i) {
    const py::ssize_t n = (py::ssize_t)db.entries.size();
    if (i < 0) i += n;
    if (i < 0 || i >= n) throw py::index_error();
    return (size_t)i;
}

void register_isomerdb(py::module_& m) {
    m.def("set_database_path", [](const std::string& p) { IsomerDB::database_path = p; },
          py::arg("path"),
          "Set the IsomerDB root directory (a setter -- never an environment variable).");
    m.def("get_database_path", []() { return IsomerDB::database_path; });

    py::class_<IsomerDB::Entry>(m, "IsomerEntry", "One database isomer record.")
        .def_property_readonly("group", &group_str, "Point group, e.g. 'Ih'.")
        .def_property_readonly("rspi", [](const IsomerDB::Entry& e) {
            return std::vector<int>(e.RSPI, e.RSPI + 12);
        }, "The 12 pentagon ring-spiral positions (1-based, matching .name() and "
           "FullereneDual.from_rspi / .rspi()).")
        .def_readonly("hlgap", &IsomerDB::Entry::HLgap, "HOMO-LUMO gap.")
        .def("__repr__", [](const IsomerDB::Entry& e) {
            return "<IsomerEntry " + group_str(e) + ">";
        });

    py::class_<IsomerDB>(m, "IsomerDB",
        "All isomers of one size, loaded from the database.")
        .def_static("read", [](int N, bool IPR) { return IsomerDB::readPDB(N, IPR); },
            py::arg("N"), py::arg("IPR") = false,
            "Load all (general or IPR) isomers of size N from the database.")
        .def_property_readonly("N", [](const IsomerDB& db) { return db.N; })
        .def_property_readonly("IPR", [](const IsomerDB& db) { return db.IPR; })
        .def_property_readonly("n_isomers",
            [](const IsomerDB& db) { return (int)db.entries.size(); })
        .def("__len__", [](const IsomerDB& db) { return db.entries.size(); })
        .def("__getitem__", [](const IsomerDB& db, py::ssize_t i) {
            return db.entries[resolve_index(db, i)];
        })
        .def("make_isomer", [](const IsomerDB& db, py::ssize_t i) {
            return PyFG::from_owned(IsomerDB::makeIsomer(db.N, db.entries[resolve_index(db, i)]));
        }, py::arg("i"), "Build the FullereneGraph for isomer i (supports negative indexing).")
        .def("__repr__", [](const IsomerDB& db) {
            return "<IsomerDB C" + std::to_string(db.N) +
                   (db.IPR ? " IPR" : "") + ": " +
                   std::to_string(db.entries.size()) + " isomers>";
        });
}
