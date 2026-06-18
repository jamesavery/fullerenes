// symmetry.cc — Python binding for Symmetry: the point group, vertex orbits and
// 3D matrix representation of a fullerene, derived PURELY from the dual graph's
// combinatorics (never from 3D coordinates).

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstring>
#include <string>
#include <vector>

#include "common.hh"

#include "fullerenes/symmetry.hh"
#include "fullerenes/triangulation.hh"

namespace py = pybind11;

using PyFD = pyf::PyGraph<FullereneDual, FullereneDualView>;

namespace {

struct PySymmetry {
    Symmetry sym;
};

// Representation3D.R (one matrix3d per group element) -> (order, 3, 3) float64.
py::array_t<double> rep3d_array(const Representation3D& rep) {
    const py::ssize_t n = (py::ssize_t)rep.R.size();
    py::array_t<double> a({n, (py::ssize_t)3, (py::ssize_t)3});
    double* d = a.mutable_data();
    for (py::ssize_t i = 0; i < n; ++i)
        std::memcpy(d + i * 9, rep.R[i].values, 9 * sizeof(double));
    return a;
}

}  // namespace

void register_symmetry(py::module_& m) {
    py::class_<PySymmetry>(m, "Symmetry",
        "Combinatorial symmetry of a fullerene dual: point group, vertex orbits, "
        "and the 3D matrix representation -- all derived purely from the graph, "
        "never from coordinates.")
        .def(py::init([](PyFD& dual) {
            return PySymmetry{ Symmetry(Triangulation(dual.view())) };
        }), py::arg("dual"), "Compute the symmetry of a FullereneDual (triangulation).")
        .def("point_group", [](PySymmetry& s) { return s.sym.point_group().to_string(); },
             "Point-group symbol, e.g. 'Ih', 'C2v', 'D5d'.")
        .def("order", [](PySymmetry& s) { return (int)s.sym.G.size(); },
             "Order of the automorphism group (number of symmetry operations).")
        .def("equivalence_classes", [](PySymmetry& s) {
            return s.sym.equivalence_classes(s.sym.G);    // vector<vector<node_t>>
        }, "Vertex orbits under the symmetry group, as list[list[int]].")
        .def("representation_3d", [](PySymmetry& s) {
            return rep3d_array(s.sym.representation_3d());
        }, "Orthogonal 3x3 matrix per group element: (order, 3, 3) float64. "
           "det = +1 (rotation) or -1 (improper).")
        .def("__repr__", [](PySymmetry& s) {
            return "<Symmetry " + s.sym.point_group().to_string() +
                   " order=" + std::to_string((int)s.sym.G.size()) + ">";
        });
}
