// delaunay_triangulation.cc — Python bindings for DelaunayTriangulation: the
// intrinsic Delaunay delta-complex (a half-edge DCEL) of a fullerene dual's cone
// metric, and its compact binary .geo serialization (GEO-FORMAT.md).
//
// No mutator is bound, so the arrays are returned as copies and the object is
// immutable from Python. A .geo record stores connectivity only, so reading one
// takes two steps: read_geo_topology gives the half-edges in the exact numbering
// from_geo builds (edge k = half-edges 2k, 2k+1), from which the caller computes
// the edge lengths that from_geo then takes as an array.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <span>
#include <string>
#include <vector>

#include "common.hh"
#include "geo_io.hh"

#include "fullerenes/delaunay.hh"
#include "fullerenes/geo-format.hh"
#include "fullerenes/triangulation.hh"

namespace py = pybind11;

using PyFD = pyf::PyGraph<FullereneDual, FullereneDualView>;
using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

namespace {

template<class T>
py::array_t<T> copy_of(std::span<T> s, int n) {
    py::array_t<T> a{py::ssize_t(n)};
    std::copy(s.begin(), s.begin() + n, a.mutable_data());
    return a;
}

// (n, 3) float array or None -> positions (coord3d is 3 packed doubles, asserted in np_interop.hh).
std::vector<coord3d> points_of(py::handle points) {
    std::vector<coord3d> x;
    if (points.is_none()) return x;
    const DoubleArray a = DoubleArray::ensure(points);
    if (!a || a.ndim() != 2 || a.shape(1) != 3) throw std::invalid_argument("points must be an (n, 3) float array");
    x.resize(size_t(a.shape(0)));
    std::memcpy(x.data(), a.data(), sizeof(double) * 3 * x.size());
    return x;
}

py::array_t<double> points_array(const std::vector<coord3d>& x) {
    py::array_t<double> a({py::ssize_t(x.size()), py::ssize_t(3)});
    std::memcpy(a.mutable_data(), x.data(), sizeof(double) * 3 * x.size());
    return a;
}

}  // namespace

void register_delaunay_triangulation(py::module_& m) {
    py::class_<DelaunayTriangulation> cls(m, "DelaunayTriangulation",
        "Intrinsic Delaunay triangulation of a cone metric: a half-edge DCEL that may be a "
        "non-simplicial delta-complex (parallel edges, self-loops). Half-edges 2k and 2k+1 "
        "are twins; he_next traverses each face counter-clockwise. Immutable from Python; "
        "the array properties are copies, with -1 marking dead half-edge slots.");

    cls.def_static("compute", [](PyFD& w) {
        return DelaunayTriangulation::compute(Triangulation(w.view()));
    }, py::arg("dual"),
       "The iDT of a fullerene dual's equilateral cone metric. Its vertices are the dual's "
       "cones (degree != 6).");

    cls.def_property_readonly("nv", [](const DelaunayTriangulation& D) { return D.nv; }, "Vertex count.");
    cls.def_property_readonly("nh", [](const DelaunayTriangulation& D) { return D.nh; },
                              "Half-edge slots, dead ones included.");
    cls.def_property_readonly("nf", [](const DelaunayTriangulation& D) { return D.nf; },
                              "Face slots, dead ones included.");
    cls.def_property_readonly("he_origin", [](const DelaunayTriangulation& D) { return copy_of(D.he_origin, D.nh); },
                              "Origin vertex of each half-edge (-1: dead slot).");
    cls.def_property_readonly("he_next", [](const DelaunayTriangulation& D) { return copy_of(D.he_next, D.nh); },
                              "The next half-edge counter-clockwise in the same face.");
    cls.def_property_readonly("he_face", [](const DelaunayTriangulation& D) { return copy_of(D.he_face, D.nh); },
                              "The face to the left of each half-edge.");
    cls.def_property_readonly("he_length", [](const DelaunayTriangulation& D) { return copy_of(D.he_length, D.nh); },
                              "Intrinsic edge length (equal on twins).");
    cls.def_property_readonly("he_angle", [](const DelaunayTriangulation& D) { return copy_of(D.he_angle, D.nh); },
                              "Corner angle at each half-edge's origin in its face.");
    cls.def_property_readonly("v_out", [](const DelaunayTriangulation& D) { return copy_of(D.v_out, D.nv); },
                              "One outgoing half-edge per vertex (-1: dead vertex).");
    cls.def_property_readonly("v_orig_degree", [](const DelaunayTriangulation& D) {
        return copy_of(D.v_orig_degree, D.nv);
    }, "Each vertex's degree in the original triangulation.");
    cls.def_property_readonly("v_cone_angle", [](const DelaunayTriangulation& D) {
        return copy_of(D.v_cone_angle, D.nv);
    }, "Each vertex's cone angle.");
    cls.def_property_readonly("f_he", [](const DelaunayTriangulation& D) { return copy_of(D.f_he, D.nf); },
                              "One boundary half-edge per face (-1: dead face).");

    cls.def("check_consistency", &DelaunayTriangulation::check_consistency,
            "The DCEL's structural and metric invariants (triangular faces, chaining, positive "
            "twin-equal lengths, triangle inequalities).");
    cls.def("is_simplicial", &DelaunayTriangulation::is_simplicial,
            "False for a delta-complex with parallel edges or self-loops.");
    cls.def("vertex_angle_sums", [](const DelaunayTriangulation& D) {
        py::array_t<double> a{py::ssize_t(D.nv)};
        auto r = a.mutable_unchecked<1>();
        for (int v = 0; v < D.nv; v++) r(v) = D.v_out[v] >= 0 ? D.vertex_angle_sum(v) : NAN;
        return a;
    }, "The total corner angle at each vertex (NaN for a dead vertex).");

    // --- Compact binary geometry (.geo, GEO-FORMAT.md) ---
    cls.def("to_geo", [](const DelaunayTriangulation& D, const std::string& path, py::object points,
                         bool append, const geo_options& opt) {
        const std::vector<coord3d> x = points_of(points);
        return pyf::write_file(path, append, [&](FILE* f) {
            return DelaunayTriangulation::to_geo(D, x, f, append, opt);
        });
    }, py::arg("path"), py::arg("points") = py::none(), py::arg("append") = false,
       py::arg("options") = geo_options{},
       "Write the connectivity -- flagged as an exact triangulation -- and optional (nv, 3) "
       "positions as one .geo record: a fresh file, or appended. points=None needs "
       "options.type == GeoType.NONE. The metric is not stored.");

    cls.def_static("read_geo_topology", [](const std::string& path, uint64_t index) {
        const pyf::unique_file f = pyf::open_for_reading(path);
        const bool has_points = geo::read_header(f.get()).opt.type != geo_type::NONE;
        // An equilateral stand-in metric, which every triangulation satisfies, only
        // to reuse from_geo's exact half-edge numbering.
        std::vector<coord3d> x;
        const DelaunayTriangulation D = DelaunayTriangulation::from_geo(
            f.get(), index, [](const DelaunayTriangulation&, int) { return 1.0; }, [](int) { return 6; }, &x);
        return py::make_tuple(copy_of(D.he_origin, D.nh), copy_of(D.he_next, D.nh),
                              has_points ? py::object(points_array(x)) : py::object(py::none()));
    }, py::arg("path"), py::arg("index") = 0,
       "(he_origin, he_next, points) of record `index`, in the half-edge numbering from_geo "
       "builds: edge k is half-edges 2k and 2k+1. points is None for a file without "
       "coordinates.");

    cls.def_static("from_geo", [](const std::string& path, uint64_t index, const DoubleArray& lengths,
                                  py::object orig_degree) {
        if (lengths.ndim() != 1) throw std::invalid_argument("lengths must be one-dimensional");
        const bool constant = py::isinstance<py::int_>(orig_degree);
        const int degree = constant ? orig_degree.cast<int>() : 0;
        const std::vector<int> degrees = constant ? std::vector<int>{} : orig_degree.cast<std::vector<int>>();
        const auto L = lengths.unchecked<1>();

        const pyf::unique_file f = pyf::open_for_reading(path);
        return DelaunayTriangulation::from_geo(f.get(), index,
            [&](const DelaunayTriangulation& D, int h) {
                if (L.shape(0) != D.nh / 2)
                    throw std::invalid_argument("lengths has " + std::to_string(L.shape(0)) + " entries for "
                                                + std::to_string(D.nh / 2) + " edges");
                if (!constant && degrees.size() != size_t(D.nv))
                    throw std::invalid_argument("orig_degree has " + std::to_string(degrees.size())
                                                + " entries for " + std::to_string(D.nv) + " vertices");
                return L(h / 2);
            },
            [&](int v) { return constant ? degree : degrees[v]; });
    }, py::arg("path"), py::arg("index") = 0, py::kw_only(), py::arg("lengths"), py::arg("orig_degree"),
       "Record `index` with the caller's metric: lengths[k] is the length of edge k "
       "(half-edges 2k, 2k+1; see read_geo_topology), orig_degree an int for every vertex or "
       "one per vertex. ValueError when the metric does not fit or fails check_consistency; "
       "RuntimeError for a record not flagged as a triangulation or a malformed file.");

    cls.def("__repr__", [](const DelaunayTriangulation& D) {
        return "<DelaunayTriangulation nv=" + std::to_string(D.nv) + " nh=" + std::to_string(D.nh)
             + " nf=" + std::to_string(D.nf) + ">";
    });
}
