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

// Mechanical method surface (the nullary scalar queries volume/surface_area/
// diameter/width_height_depth/is_invalid and the planar-graph queries), emitted
// by tools/gen_bindings.py.  The inertia methods are NOT among them: they take a
// MassModel argument, so they are hand-written below and listed in spec.toml's
// skip list.
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
    py::enum_<MassModel>(m, "MassModel",
        "Which mass distribution an inertia tensor / principal frame describes. "
        "The two are NOT interchangeable: over 40 C60/C70 cages the tensors agree "
        "to 5.9%, but an axis both models resolve turns by up to 3.40 deg between "
        "the two frames -- 3.70 deg at a marginal eigenvalue gap.")
        .value("SOLID", MassModel::Solid,
               "Uniform density over the ENCLOSED SOLID (the default). Exact per "
               "triangle. Requires a closed, consistently oriented surface.")
        .value("ATOMS", MassModel::Atoms,
               "Uniform mass at the VERTICES -- the molecular convention. Requires "
               "no topology at all.");

    py::class_<PyPoly> cls(m, "Polyhedron",
        "Planar graph with 3D vertex coordinates. .points is a zero-copy, "
        "writeable numpy view.");

    pyf::bind_geom_common<Polyhedron, PolyhedronView<double>>(cls);
    register_generated_Polyhedron(cls);   // volume/surface_area/diameter/...

    // --- Construction ---
    cls.def_static("from_fullerene", [](PyFG& fg, bool verbose) {
        pyf::FdSilencer hush(!verbose);   // the Fortran force-field optimizer is chatty
        return PyPoly::from_owned(Polyhedron::fullerene_polyhedron(FullereneGraph(fg.view())));
    }, py::arg("graph"), py::arg("verbose") = false,
       "Build a 3D fullerene polyhedron (force-field pipeline). verbose=True shows "
       "the optimizer log. The result is centred and aligned with MassModel.ATOMS "
       "-- the molecular convention -- so a subsequent DEFAULTED align_with_axes() "
       "(MassModel.SOLID) rotates it back out of that frame, by up to 3.70 deg.");

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

    // --- Geometry queries (custom returns; the nullary scalar ones are generated,
    //     the matrix ones are hand-written below because they take a MassModel) ---
    cls.def("bounding_box", [](PyPoly& w) {
        auto bb = w.view().bounding_box();
        return py::make_tuple(pyf::coord3d_copy(bb.first), pyf::coord3d_copy(bb.second));
    }, "(min_corner, max_corner) as two (3,) arrays.");
    cls.def("faces", [](PyPoly& w) { return pyf::faces_copy(w.view().faces()); },
            "Polygon faces as list[list[int]].");

    // Parameterised, so hand-written rather than generated (the generator binds
    // nullary methods only). The default is MassModel.SOLID -- the moments of the
    // enclosed solid; MassModel.ATOMS is the molecular convention (uniform mass at
    // the atoms), which is what these returned before 2026-08-07.
    cls.def("inertia_matrix", [](PyPoly& w, MassModel mass_model) {
        return pyf::matrix3d_copy(w.view().inertia_matrix(mass_model));
    }, py::arg("mass_model") = MassModel::Solid,
       "CENTRAL inertia tensor -- about the mass distribution's own centre, so it "
       "does not depend on where the polyhedron sits. Default MassModel.SOLID: the "
       "enclosed solid, uniform density, exact per triangle -- which REQUIRES a "
       "closed, consistently oriented surface (it integrates over faces()). "
       "MassModel.ATOMS: uniform mass at the vertices, needing no topology. "
       "A degenerate or non-integrable distribution yields the ZERO matrix.");
    cls.def("principal_axes", [](PyPoly& w, MassModel mass_model) {
        return pyf::matrix3d_copy(w.view().principal_axes(mass_model));
    }, py::arg("mass_model") = MassModel::Solid,
       "Eigenvectors of inertia_matrix(mass_model) as rows, by ascending "
       "eigenvalue -- so row 0 is the longest axis. Same precondition: "
       "MassModel.SOLID needs a closed oriented surface. Returns the IDENTITY when "
       "no frame can be built (degenerate mass, non-PSD second moment, non-finite "
       "tensor, non-unitary eigenvectors) -- a legal rotation, so it is "
       "indistinguishable from a real frame here.");

    // --- In-place geometry (writes back through .points) ---
    cls.def("optimize", [](PyPoly& w, int method, double ftol, bool verbose) {
        pyf::FdSilencer hush(!verbose);   // silence the Fortran optimizer log by default
        return w.view().optimize(method, ftol);
    }, py::arg("method") = 3, py::arg("ftol") = 1e-10, py::arg("verbose") = false,
       "Force-field optimize in place (writes into .points). Returns success. "
       "verbose=True shows the optimizer log.");
    cls.def("move_to_origin", [](PyPoly& w) { w.view().move_to_origin(); });
    cls.def("align_with_axes", [](PyPoly& w, MassModel mass_model) {
        w.view().align_with_axes(mass_model);
    }, py::arg("mass_model") = MassModel::Solid,
       "Rotate into the principal frame of inertia_matrix(mass_model), in place. "
       "Default MassModel.SOLID, which REQUIRES a closed, consistently oriented "
       "surface; MassModel.ATOMS needs no topology. NOTE: from_fullerene() leaves "
       "the cage in its ATOMS frame, so a defaulted align_with_axes() rotates it "
       "OUT of that frame -- by up to 3.70 deg on a low-symmetry cage. Pass "
       "MassModel.ATOMS to keep the molecular convention.");
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
