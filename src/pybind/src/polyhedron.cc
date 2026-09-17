// polyhedron.cc — Python bindings for Polyhedron (planar graph + 3D coords).
//
// Geometry queries and the in-place optimizer run through the transient view, so
// optimize()/move/align write straight back into .points (zero-copy in both owned
// and bring-your-own-buffer modes). I/O statics take a const Polyhedron&, so for
// those we materialise one (owned: copy; view: deep copy from the view).

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdio>
#include <string>
#include <vector>

#include "common.hh"
#include "geo_io.hh"

#include "fullerenes/deltahedron.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/geo-format.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/triangulation.hh"

namespace py = pybind11;

using PyPoly  = pyf::PyGeom<Polyhedron, PolyhedronView<double>>;
using PyDelta = pyf::PyGeom<Deltahedron, DeltahedronView<double>>;
using PyFG    = pyf::PyGraph<FullereneGraph, FullereneGraphView>;
using PyFD    = pyf::PyGraph<FullereneDual, FullereneDualView>;

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

// A Polyhedron to hand to the C++ writers: the owned object itself, or a copy of
// a View-mode wrapper's arrays.
static Polyhedron materialize(PyPoly& w) {
    if (w.mode == PyPoly::Mode::Owned) return w.owned;
    PolyhedronView<double> v = w.view();
    return Polyhedron(v, std::vector<coord3d>(v.points.begin(), v.points.end()));
}

namespace {

using pyf::open_for_reading;
using pyf::unique_file;
using pyf::write_file;

// The adjacency of any bound graph or geometry object, for coordinate-only files.
PlanarGraphView planar_view_of(py::handle g) {
    if (py::isinstance<PyPoly>(g))  return g.cast<PyPoly&>().view();
    if (py::isinstance<PyDelta>(g)) return g.cast<PyDelta&>().view();
    if (py::isinstance<PyFG>(g))    return g.cast<PyFG&>().view();
    if (py::isinstance<PyFD>(g))    return g.cast<PyFD&>().view();
    throw py::type_error("graph must be a Polyhedron, Deltahedron, FullereneGraph or FullereneDual");
}

std::string repr(const geo_options& o) {
    static const char* const types[] = {"FIXED", "NONE", "F32", "F64"};
    std::string s = std::string("GeoOptions(type=GeoType.") + types[int(o.type)];
    if (o.type == geo_type::FIXED)
        s += ", width=" + std::to_string(o.width) + ", scale=" + py::repr(py::float_(o.scale)).cast<std::string>();
    if (o.offset) s += ", offset=True";
    s += std::string(", graph=") + (o.graph ? "True" : "False");
    if (o.triangulation) s += ", triangulation=True";
    if (o.record_n) s += ", record_n=True";
    s += ", N=" + std::to_string(o.N) + ", deg_min=" + std::to_string(o.deg_min)
       + ", deg_bits=" + std::to_string(o.deg_bits);
    if (o.sync) s += ", sync=True";
    return s + ")";
}

uint8_t checked_width(int width) {
    if (width < 0 || width > 255) throw std::invalid_argument("width " + std::to_string(width) + " outside 0..255");
    return uint8_t(width);
}

void register_geo_types(py::module_& m) {
    py::enum_<geo_type>(m, "GeoType", "Coordinate storage of a .geo file (GEO-FORMAT.md sec. 6).")
        .value("FIXED", geo_type::FIXED, "Fixed point of GeoOptions.width bits (2..30), with a scale per record "
                                         "(GeoOptions.scale == 0) or for the file.")
        .value("NONE", geo_type::NONE, "No coordinates: an archive of surface graphs.")
        .value("F32", geo_type::F32, "IEEE float32.")
        .value("F64", geo_type::F64, "IEEE float64 (lossless).");

    py::class_<geo_options>(m, "GeoOptions",
        "What a .geo writer asks for (GEO-FORMAT.md sec. 13). N = 0 and deg_min = deg_bits = -1 "
        "are derived tightly from the polyhedra a file is created with; a file created by "
        "appending must declare its degree range (with a graph) and N (with record_n).")
        .def(py::init([](geo_type type, int width, float scale, bool offset, bool graph, bool triangulation,
                         bool record_n, uint32_t N, int deg_min, int deg_bits, bool sync) {
                 geo_options o;
                 o.type = type; o.width = checked_width(width); o.scale = scale; o.offset = offset;
                 o.graph = graph; o.triangulation = triangulation; o.record_n = record_n; o.N = N;
                 o.deg_min = deg_min; o.deg_bits = deg_bits; o.sync = sync;
                 return o;
             }),
             py::kw_only(), py::arg("type") = geo_type::F64, py::arg("width") = 0, py::arg("scale") = 0.0f,
             py::arg("offset") = false, py::arg("graph") = true, py::arg("triangulation") = false,
             py::arg("record_n") = false, py::arg("N") = 0, py::arg("deg_min") = -1, py::arg("deg_bits") = -1,
             py::arg("sync") = false)
        .def_readwrite("type", &geo_options::type)
        .def_property("width", [](const geo_options& o) { return int(o.width); },
                               [](geo_options& o, int w) { o.width = checked_width(w); },
                      "Fixed-point bits, 2..30 (0 for the float types and NONE).")
        .def_readwrite("scale", &geo_options::scale, "FIXED: 0 = a scale per record, > 0 = the file's scale.")
        .def_readwrite("offset", &geo_options::offset, "FIXED: an offset per record (the bounding-box midpoint).")
        .def_readwrite("graph", &geo_options::graph, "Store the connectivity (a half-edge twin matching).")
        .def_readwrite("triangulation", &geo_options::triangulation, "Every face is a triangle (checked).")
        .def_readwrite("record_n", &geo_options::record_n, "A vertex count per record; N is then a capacity.")
        .def_readwrite("N", &geo_options::N, "Vertex capacity; 0 = derive.")
        .def_readwrite("deg_min", &geo_options::deg_min, "Smallest admissible degree; -1 = derive.")
        .def_readwrite("deg_bits", &geo_options::deg_bits, "Degree field width; -1 = derive.")
        .def_readwrite("sync", &geo_options::sync, "fsync after each write step.")
        .def("__repr__", [](const geo_options& o) { return repr(o); });

    py::class_<geo_header>(m, "GeoHeader", "A .geo file's header (GEO-FORMAT.md sec. 4).")
        .def_property_readonly("options", [](const geo_header& h) { return h.opt; },
                               "The file's resolved options (a copy).")
        .def_readonly("count", &geo_header::count, "Number of records.")
        .def_readonly("checksum", &geo_header::checksum, "The stored checksum (sec. 9).")
        .def("record_size", &geo_header::record_size, "R, the size of every record in bytes.")
        .def("record_offset", &geo_header::record_offset, py::arg("index"), "Where record `index` starts.")
        .def("edge_capacity", &geo_header::edge_capacity, "E_cap, the edge slots of every record (sec. 7.2).")
        .def("__repr__", [](const geo_header& h) {
            return "<GeoHeader count=" + std::to_string(h.count) + " record_size=" + std::to_string(h.record_size())
                 + " " + repr(h.opt) + ">";
        });
}

}  // namespace

void register_polyhedron(py::module_& m) {
    register_geo_types(m);

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

    // --- Compact binary geometry (.geo, GEO-FORMAT.md) ---
    cls.def("to_geo", [](PyPoly& w, const std::string& path, bool append, const geo_options& opt) {
        return with_polyhedron(w, [&](const Polyhedron& P) {
            return write_file(path, append, [&](FILE* f) { return Polyhedron::to_geo(P, f, append, opt); });
        });
    }, py::arg("path"), py::arg("append") = false, py::arg("options") = geo_options{},
       "Write as a single-record .geo file, or append a record (creating the file when "
       "absent). Returns False when a write failed. Raises ValueError for invalid "
       "options and RuntimeError when the polyhedron does not fit the file.");
    cls.def_static("to_geo_batch", [](const std::string& path, std::vector<PyPoly*> polyhedra,
                                      const geo_options& opt) {
        std::vector<Polyhedron> Ps;
        Ps.reserve(polyhedra.size());
        for (PyPoly* w : polyhedra) Ps.push_back(materialize(*w));
        return write_file(path, false, [&](FILE* f) { return Polyhedron::to_geo(Ps, f, opt); });
    }, py::arg("path"), py::arg("polyhedra"), py::arg("options") = geo_options{},
       "Write a fresh .geo file holding every polyhedron, with the options' derived "
       "fields taken tightly from all of them.");
    cls.def_static("from_geo", [](const std::string& path, uint64_t index, py::object graph) {
        const unique_file f = open_for_reading(path);
        if (graph.is_none()) return PyPoly::from_owned(Polyhedron::from_geo(f.get(), index));
        return PyPoly::from_owned(Polyhedron::from_geo(f.get(), planar_view_of(graph), index));
    }, py::arg("path"), py::arg("index") = 0, py::arg("graph") = py::none(),
       "Read record `index`. A file without connectivity needs `graph` (a Polyhedron, "
       "Deltahedron, FullereneGraph or FullereneDual with the record's vertex count). "
       "IndexError past the last record; RuntimeError for a malformed file or a "
       "record with self-loops or parallel edges.");
    cls.def_static("read_geo_header", [](const std::string& path) {
        return Polyhedron::read_geo_header(open_for_reading(path).get());
    }, py::arg("path"), "The header of a .geo file.");
    cls.def_static("verify_geo", [](const std::string& path) {
        return Polyhedron::verify_geo(open_for_reading(path).get());
    }, py::arg("path"), "Whether a .geo file's stored checksum matches its records.");

    cls.def("__repr__", [](PyPoly& w) {
        return "<Polyhedron N=" + std::to_string(w.N()) + ">";
    });
}
