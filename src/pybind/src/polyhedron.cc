// polyhedron.cc — Python bindings for Polyhedron (planar graph + 3D coords).
//
// Geometry queries and the in-place optimizer run through the transient view, so
// optimize()/move/align write straight back into .points (zero-copy in both owned
// and bring-your-own-buffer modes). I/O statics take a const Polyhedron&, so for
// those we materialise one (owned: copy; view: deep copy from the view).

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include "common.hh"

#include "fullerenes/polyhedron.hh"
#include "fullerenes/fullerenegraph.hh"

namespace py = pybind11;

using PyPoly = pyf::PyGeom<Polyhedron, PolyhedronView<double>>;
using PyFG = pyf::PyGraph<FullereneGraph, FullereneGraphView>;

// Mechanical method surface (volume/surface_area/diameter/inertia/... and the
// planar-graph queries), emitted by tools/gen_bindings.py.
void register_generated_Polyhedron(py::class_<PyPoly>& cls);

// Run f on a const Polyhedron& for the const-& I/O statics. Owned mode passes
// w.owned directly (no copy); only View mode materialises a temporary.
template<class F>
static auto with_polyhedron(PyPoly& w, F&& f) {
    if (w.mode == PyPoly::Mode::Owned) return f(w.owned);
    PolyhedronView<double> v = w.view();
    Polyhedron tmp(v, std::vector<coord3d>(v.points.begin(), v.points.end()));
    return f(tmp);
}

void register_polyhedron(py::module_& m) {
    py::class_<PyPoly> cls(m, "Polyhedron",
        "Planar graph with 3D vertex coordinates. .points is a zero-copy, "
        "writeable numpy view.");

    pyf::bind_geom_common<Polyhedron, PolyhedronView<double>>(cls);
    register_generated_Polyhedron(cls);   // volume/surface_area/diameter/inertia/...

    // --- Construction ---
    cls.def_static("from_fullerene", [](PyFG& fg, bool verbose) {
        pyf::FdSilencer hush(!verbose);   // the Fortran force-field optimizer is chatty
        return PyPoly::from_owned(Polyhedron::fullerene_polyhedron(FullereneGraph(fg.view())));
    }, py::arg("graph"), py::arg("verbose") = false,
       "Build a 3D fullerene polyhedron (force-field pipeline). verbose=True shows "
       "the optimizer log.");

    cls.def_static("C20", []() { return PyPoly::from_owned(Polyhedron::C20()); });

    cls.def_static("from_arrays",
        [](py::array neighbours, py::array points, py::object deg_obj) {
            // A Polyhedron is general (duals are deg 5/6), so deg cannot default
            // to cubic all-3 -- when omitted, infer it from the -1 padding.
            py::array deg = deg_obj.is_none()
                ? py::array(pyf::infer_deg(neighbours))
                : deg_obj.cast<py::array>();
            return PyPoly::from_arrays(std::move(neighbours), std::move(deg), std::move(points));
        },
        py::arg("neighbours"), py::arg("points"), py::arg("deg") = py::none(),
        "Wrap caller-owned int32 (N,dmax) + float64 (N,3) arrays (zero-copy). "
        "deg is inferred from the -1 padding when omitted (any degree, not just cubic).");

    cls.def_static("from_file", [](const std::string& path) {
        return PyPoly::from_owned(Polyhedron::from_file(path));
    }, py::arg("path"), "Read a polyhedron from a file (format by extension).");

    // --- Geometry queries (custom returns; the scalar/matrix ones are generated) ---
    cls.def("bounding_box", [](PyPoly& w) {
        auto bb = w.view().bounding_box();
        return py::make_tuple(pyf::coord3d_copy(bb.first), pyf::coord3d_copy(bb.second));
    }, "(min_corner, max_corner) as two (3,) arrays.");
    cls.def("faces", [](PyPoly& w) { return pyf::faces_copy(w.view().faces()); },
            "Polygon faces as list[list[int]].");

    // --- In-place geometry (writes back through .points) ---
    cls.def("optimize", [](PyPoly& w, int method, double ftol, bool verbose) {
        pyf::FdSilencer hush(!verbose);   // silence the Fortran optimizer log by default
        return w.view().optimize(method, ftol);
    }, py::arg("method") = 3, py::arg("ftol") = 1e-10, py::arg("verbose") = false,
       "Force-field optimize in place (writes into .points). Returns success. "
       "verbose=True shows the optimizer log.");
    cls.def("move_to_origin", [](PyPoly& w) { w.view().move_to_origin(); });
    cls.def("align_with_axes", [](PyPoly& w) { w.view().align_with_axes(); });
    cls.def("scale", [](PyPoly& w, py::handle s) { w.view().scale(pyf::as_coord3d(s)); },
            py::arg("s"), "Scale per-axis by a length-3 vector (in place).");
    cls.def("move", [](PyPoly& w, py::handle d) { w.view().move(pyf::as_coord3d(d)); },
            py::arg("d"), "Translate by a length-3 vector (in place).");

    // --- Transforms (return new Polyhedron) ---
    cls.def("convex_hull", [](PyPoly& w) {
        return PyPoly::from_owned(w.view().incremental_convex_hull());
    });
    cls.def("dual", [](PyPoly& w) { return PyPoly::from_owned(w.view().dual()); });
    cls.def("leapfrog_dual", [](PyPoly& w) { return PyPoly::from_owned(w.view().leapfrog_dual()); });

    // --- I/O ---
    cls.def("write", [](PyPoly& w, const std::string& path) {
        return with_polyhedron(w, [&](const Polyhedron& P) { return Polyhedron::to_file(P, path); });
    }, py::arg("path"), "Write to a file; format chosen by extension (.mol2/.xyz/.obj/...).");
    cls.def("to_povray", [](PyPoly& w) {
        return with_polyhedron(w, [](const Polyhedron& P) { return P.to_povray(); });
    });
    cls.def("to_latex", [](PyPoly& w) {
        return with_polyhedron(w, [](const Polyhedron& P) { return P.to_latex(); });
    });

    cls.def("__repr__", [](PyPoly& w) {
        return "<Polyhedron N=" + std::to_string(w.N()) + ">";
    });
}
