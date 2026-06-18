// fullerene_graph.cc — Python bindings for FullereneGraph (the cubic molecule
// graph, dmax=3). Hand-written for now; the mechanical accessors come from
// bind_graph_common, the class-specific methods will later be generated.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>

#include "common.hh"
#include "naming.hh"

#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/planargraph.hh"
#include "fullerenes/spiral.hh"

namespace py = pybind11;

using PyFG = pyf::PyGraph<FullereneGraph, FullereneGraphView>;
using PyFD = pyf::PyGraph<FullereneDual, FullereneDualView>;

// Mechanical method surface, emitted by tools/gen_bindings.py.
void register_generated_FullereneGraph(py::class_<PyFG>& cls);

void register_fullerene_graph(py::module_& m) {
    py::class_<PyFG> cls(m, "FullereneGraph",
        "Cubic fullerene graph (3-regular, 12 pentagons). Zero-copy over numpy.");

    pyf::bind_graph_common<FullereneGraph, FullereneGraphView>(cls);
    register_generated_FullereneGraph(cls);   // generated mechanical methods

    cls.def_static("C20", []() {
        return PyFG::from_owned(FullereneGraph::C20());
    }, "The dodecahedral C20 fullerene graph.");

    cls.def_static("from_arrays",
        [](py::array neighbours, py::object deg_obj) {
            // deg inferred from the -1 padding when omitted -- for a cubic graph
            // (full (N,3) rows) that is all-3; same path as Polyhedron/Deltahedron.
            py::array deg = deg_obj.is_none()
                ? py::array(pyf::infer_deg(neighbours))
                : deg_obj.cast<py::array>();
            return PyFG::from_arrays(std::move(neighbours), std::move(deg));
        },
        py::arg("neighbours"), py::arg("deg") = py::none(),
        "Wrap a caller-owned numpy int32 (N,3) array (zero-copy; no implicit "
        "conversion). deg is inferred from the -1 padding when omitted.");

    // Canonical name via the dual triangulation -- one shared naming helper, so
    // FullereneGraph(g).name() and g.dual().name() can never disagree.
    cls.def("name", [](PyFG& w) {
        FullereneDual dual(w.view().dual_graph());
        return pyf::canonical_name(dual);
    }, "Canonical generalized-spiral name, e.g. 'C20-[GS:...]-fullerene'.");

    cls.def("dual", [](PyFG& w) {
        FullereneGraphView v = w.view();
        return PyFD::from_owned(FullereneDual(v.dual_graph()));
    }, "The dual triangulation (FullereneDual).");

    // --- Transforms (return a new FullereneGraph) ---
    cls.def("gc_transform", [](PyFG& w, unsigned k, unsigned l) {
        return PyFG::from_owned(w.view().GCtransform(k, l));
    }, py::arg("k"), py::arg("l") = 0u, "Goldberg-Coxeter transform GC(k,l).");
    cls.def("leapfrog", [](PyFG& w) {
        return PyFG::from_owned(w.view().leapfrog_fullerene());
    }, "Leapfrog transform (GC(1,1)).");
    cls.def("halma", [](PyFG& w, int n) {
        return PyFG::from_owned(w.view().halma_fullerene(n));
    }, py::arg("n"), "Halma transform of order n (GC(n+1,0)).");

    cls.def("__repr__", [](PyFG& w) {
        return "<FullereneGraph C" + std::to_string(w.N()) + ">";
    });
}
