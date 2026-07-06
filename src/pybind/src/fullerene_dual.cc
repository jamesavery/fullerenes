// fullerene_dual.cc — Python bindings for FullereneDual (the dual triangulation,
// degree 5/6, dmax=6). This is what buckygen and IsomerDB produce and the source
// of the canonical name. Mechanical accessors come from bind_graph_common +
// the generated surface; the class-specific methods are hand-written here.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <string>
#include <utility>
#include <vector>

#include "common.hh"
#include "naming.hh"

#include "fullerenes/triangulation.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/eisenstein.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"

#include <stdexcept>

namespace py = pybind11;

using PyFD = pyf::PyGraph<FullereneDual, FullereneDualView>;
using PyFG = pyf::PyGraph<FullereneGraph, FullereneGraphView>;

void register_generated_FullereneDual(py::class_<PyFD>& cls);

// Carbon count from dual vertex count: Nv = N/2 + 2  =>  N = 2*Nv - 4.
static int carbon_count(int Nv) { return 2 * Nv - 4; }

// Validate a node index / out-arc axis from Python before it reaches an unchecked
// operator[] in the library's geodesics routines (which index by node id and by
// out-arc, with no bounds check -- an out-of-range value is an OOB read).
static void require_node(PyFD& w, int u) {
    if (u < 0 || u >= w.N())
        throw py::index_error("node index " + std::to_string(u) +
                              " out of range [0," + std::to_string(w.N()) + ")");
}
static void require_axis(PyFD& w, int u, int axis) {
    const int deg = w.deg_data()[u];        // caller must require_node(w,u) first
    if (axis < 0 || axis >= deg)
        throw py::index_error("axis " + std::to_string(axis) + " out of range [0," +
                              std::to_string(deg) + ") for node " + std::to_string(u));
}

void register_fullerene_dual(py::module_& m) {
    py::class_<PyFD> cls(m, "FullereneDual",
        "Dual fullerene triangulation (12 degree-5 vertices, rest degree-6). "
        "Zero-copy over numpy.");

    pyf::bind_graph_common<FullereneDual, FullereneDualView>(cls);
    register_generated_FullereneDual(cls);

    cls.def_static("from_rspi",
        [](int N, std::vector<int> rspi, py::object jumps_obj) {
            const int Nv = N / 2 + 2;                 // dual vertex count
            if (rspi.size() != 12)
                throw std::invalid_argument("rspi must list exactly 12 pentagon positions");
            std::vector<int> rspi0;                   // 0-based for the library ctor
            rspi0.reserve(12);
            for (int i : rspi) {
                if (i < 1 || i > Nv)
                    throw std::invalid_argument("rspi position " + std::to_string(i) +
                        " out of range [1, " + std::to_string(Nv) + "] for C" + std::to_string(N));
                rspi0.push_back(i - 1);
            }
            jumplist_t jumps;
            if (!jumps_obj.is_none())
                for (auto item : jumps_obj.cast<py::list>())
                    jumps.push_back(item.cast<std::pair<int, int>>());
            return PyFD::from_owned(FullereneDual(N, rspi0, jumps));
        },
        py::arg("N"), py::arg("rspi"), py::arg("jumps") = py::none(),
        "From carbon count N and the 12 pentagon spiral positions (1-based, "
        "matching .name() and IsomerEntry.rspi; round-trips with .rspi()).");

    cls.def_static("from_arrays",
        [](py::array neighbours, py::array deg) {
            return PyFD::from_arrays(std::move(neighbours), std::move(deg));
        },
        py::arg("neighbours"), py::arg("deg"),
        "Wrap caller-owned numpy int32 (N,6) + uint8 (N,) arrays (zero-copy; no "
        "implicit conversion). deg is required (degrees are 5/6).");

    cls.def("rspi", [](PyFD& w) {
        FullereneDualView v = w.view();
        std::vector<int> r;
        jumplist_t j;
        v.get_rspi(r, j, /*general=*/true, /*pentagon_start=*/true);
        for (int& x : r) x += 1;                      // 0-based (library) -> 1-based (API)
        py::list jumps;
        for (auto& p : j) jumps.append(py::make_tuple(p.first, p.second));
        return py::make_tuple(r, jumps);
    }, "Ring-spiral pentagon positions (1-based, matching .name() and "
       "IsomerEntry.rspi) and jumps: (list[int], list[(int,int)]). Feeds back "
       "into from_rspi().");

    cls.def("name", [](PyFD& w) {
        return pyf::canonical_name(w.view());
    }, "Canonical generalized-spiral name, e.g. 'C60-[GS:...]-fullerene'.");

    cls.def("dual", [](PyFD& w) {
        FullereneDualView v = w.view();
        PlanarGraph cubic = v.dual_graph();
        // FullereneGraph's constructor abort()s on a non-fullerene; from_arrays
        // does not validate topology, so check here and raise instead of killing
        // the interpreter.
        if (!cubic.is_a_fullerene(/*verbose=*/false))
            throw std::invalid_argument(
                "dual(): this triangulation is not a valid fullerene dual "
                "(from_arrays does not validate topology)");
        return PyFG::from_owned(FullereneGraph(cubic));
    }, "The cubic fullerene graph (dual of this triangulation).");

    // --- Transforms (return a new FullereneDual) ---
    cls.def("gc_transform", [](PyFD& w, unsigned k, unsigned l) {
        return PyFD::from_owned(FullereneDual(w.view().GCtransform(k, l)));
    }, py::arg("k"), py::arg("l") = 0u, "Goldberg-Coxeter transform GC(k,l).");
    cls.def("halma", [](PyFD& w, int m) {
        return PyFD::from_owned(FullereneDual(w.view().halma_transform(m)));
    }, py::arg("m"), "Halma transform of order m.");

    // --- Surface distances + simple geodesics (the surface_distances algorithm) ---
    cls.def("surface_distances",
        [](PyFD& w, std::vector<int> only_nodes, bool calc_self) {
            for (int u : only_nodes) require_node(w, u);
            return pyf::matrix_copy<double>(w.view().surface_distances(only_nodes, calc_self));
        },
        py::arg("only_nodes") = std::vector<int>{}, py::arg("calc_self") = false,
        "Pairwise surface (geodesic) distances over the intrinsic equilateral metric; "
        "(n,n) float64. Empty only_nodes = all vertices.");

    cls.def("simple_geodesics_from",
        [](PyFD& w, int u, bool calc_self) {
            require_node(w, u);
            auto V = w.view();
            auto G = V.simple_geodesics({}, calc_self);    // n x n, identity index map
            py::list out;
            const int n = (int)G.n;
            for (int v = 0; v < n; v++) {
                if (v == u && !calc_self) continue;
                const auto& sg = G(u, v);
                int a = sg.g.first, b = sg.g.second, d2 = sg.g.norm2();
                if (d2 == 0) continue;                     // unreached / self
                out.append(py::make_tuple(v, sg.axis, a, b, d2));
            }
            return out;
        },
        py::arg("u"), py::arg("calc_self") = false,
        "For source u: the shortest simple geodesic to every reachable vertex as "
        "(v, axis, a, b, d2 = a^2+ab+b^2) -- the straight Eisenstein rays the "
        "surface_distances algorithm shoots from u.");

    cls.def("simple_geodesics_trace",
        [](PyFD& w, int u) {
            require_node(w, u);
            auto V = w.view();
            std::vector<TriangulationView::geodesic_step> trace;
            int M = 0;
            V.simple_geodesics({}, /*calc_self=*/false, /*trace_u=*/u, &trace, &M);
            py::list steps;
            for (const auto& s : trace)
                steps.append(py::make_tuple(s.axis, s.a, s.b, (int)s.v,
                                            s.d2, s.H_before, s.improved));
            return py::make_tuple(M, steps);
        },
        py::arg("u"),
        "Per-ray trace of the simple-geodesics search from source u, in search "
        "order (axis outer, then a, then b): (M, steps), where M = the search "
        "radius (max graph distance from u) and each step is "
        "(axis, a, b, v, d2 = a^2+ab+b^2, H_before, improved). improved == "
        "(d2 <= H_before) == 'the witness geodesic to v was (re)assigned' "
        "(H[u][v] = min(d2, H_before)). Drives a step-by-step algorithm animation.");

    cls.def("alexandrov_polytope",
        [](PyFD& w) {
            // Bobenko-Izmestiev Alexandrov embedding of the 12-cone iDT metric:
            // the unique convex polytope isometric to the dual's intrinsic
            // Delaunay triangulation. compute() sorts flats last (cones first),
            // so positions[k] is iDT vertex k; cone_spiral[k] is its original
            // (spiral) dual-vertex index, so callers can match it to the
            // deltahedral cones.
            Triangulation T(w.view());
            std::vector<int> pi = T.sort_flat_last();           // pi[u_old] = u_new
            AlexandrovSolver solver;
            solver.D = DelaunayTriangulation::compute(T);
            auto poly = solver.solve_polytope();

            std::vector<int> orig(T.N);
            for (int u = 0; u < T.N; u++) orig[pi[u]] = u;

            py::list pos, cells, cone_spiral;
            for (const auto& p : poly.positions)
                pos.append(py::make_tuple(p[0], p[1], p[2]));
            for (size_t k = 0; k < poly.positions.size(); k++)
                cone_spiral.append(orig[k]);                    // sorted k -> spiral index
            for (const auto& cell : poly.tesselation.cells) {
                py::list c;
                for (const auto& [lbl, len2] : cell) c.append(lbl);
                cells.append(c);
            }
            return py::make_tuple(pos, cells, cone_spiral, (int)poly.status,
                                  std::string(AlexandrovSolver::status_str(poly.status)));
        },
        "The 12-point Alexandrov (Bobenko-Izmestiev) polytope of the iDT metric: "
        "(positions[12] (x,y,z), cells (CCW vertex-label polygons of the T-bar(0) "
        "2-skeleton), cone_spiral[12] (each iDT vertex's original spiral dual index), "
        "status:int, status_str). A degenerate (drum-cap) status with coplanar "
        "positions is the doubly-covered-polygon case.");

    cls.def("alexandrov_trajectory",
        [](PyFD& w) {
            // The Bobenko-Izmestiev homotopy trajectory toward κ=0: one entry per
            // continuation/Newton step, each the FULL generalized-convex-polytope
            // (GCP) state — apex at the origin, the 12 cones at |p_v|=r_v, the
            // per-cone angle defect κ_v, and the triangular faces of T (which
            // change at edge flips).  Drives the GCP-pyramid-deformation animation.
            // Same vertex indexing as alexandrov_polytope (cones first); cone_spiral
            // maps a trajectory vertex index to its original spiral dual index.
            Triangulation T(w.view());
            std::vector<int> pi = T.sort_flat_last();           // pi[u_old] = u_new
            AlexandrovSolver solver;
            solver.D = DelaunayTriangulation::compute(T);
            solver.record_trajectory = true;
            solver.solve();

            std::vector<int> orig(T.N);
            for (int u = 0; u < T.N; u++) orig[pi[u]] = u;
            py::list cone_spiral;
            for (int k = 0; k < T.N; k++) cone_spiral.append(orig[k]);

            py::list steps;
            for (const auto& e : solver.trajectory) {
                py::list pos, kappa, faces;
                for (const auto& p : e.positions) pos.append(py::make_tuple(p[0], p[1], p[2]));
                for (double k : e.kappa) kappa.append(k);
                for (const auto& f : e.faces) faces.append(py::make_tuple(f[0], f[1], f[2]));
                py::dict d;
                d["phase"]     = std::string(1, e.phase);
                d["step"]      = e.step;
                d["t"]         = e.t;
                d["kappa_max"] = e.kappa_max;
                d["kappa"]     = kappa;
                d["positions"] = pos;
                d["faces"]     = faces;
                steps.append(d);
            }
            return py::make_tuple(steps, cone_spiral, (int)solver.stats_status,
                                  std::string(AlexandrovSolver::status_str(solver.stats_status)));
        },
        "The Bobenko-Izmestiev homotopy trajectory toward κ=0: (steps, cone_spiral, "
        "status:int, status_str), where steps is a list of per-step dicts "
        "{phase:'T'|'P'|'N', step:int, t:float, kappa_max:float, kappa:list[12], "
        "positions:list[12] (x,y,z; apex at origin; empty if Gram-BFS failed at high "
        "kappa), faces:list[(i,j,k)] (triangles of T at that step)}. Drives the "
        "GCP-pyramid-deformation animation; t runs 1->0 over continuation, then 0 in "
        "Newton polish.");

    cls.def("end_of_the_line",
        [](PyFD& w, int u0, int axis, int a, int b) {
            require_node(w, u0); require_axis(w, u0, axis);
            return (int)w.view().end_of_the_line(u0, axis, a, b);
        },
        py::arg("u0"), py::arg("axis"), py::arg("a"), py::arg("b"),
        "Landing vertex of the straight Eisenstein line (a,b) from u0 along the "
        "axis-th out-arc. Axis-aligned prefixes (b=0) give the on-axis vertex path.");

    cls.def("geodesic_strip",
        [](PyFD& w, int u0, int axis, int a, int b) {
            require_node(w, u0); require_axis(w, u0, axis);
            auto V = w.view();
            std::vector<std::vector<Eisenstein>> coords;
            auto quads = V.quads_of_the_line(u0, axis, a, b, &coords);
            // Cross-check the library walk against the lattice annotation (fail loud).
            int term = (int)V.end_of_the_line(u0, axis, a, b);
            if (!quads.empty() && !quads.back().empty()) {
                int last_v = (int)quads.back().back();
                const Eisenstein& last_e = coords.back().back();
                if (last_v != term || last_e.first != a || last_e.second != b)
                    throw std::runtime_error("geodesic_strip: strip terminus disagrees with "
                        "end_of_the_line/(a,b) -- valid only for a,b>=1");
            }
            py::list out;
            for (size_t ri = 0; ri < quads.size(); ri++) {
                py::list run;
                for (size_t k = 0; k < quads[ri].size(); k++)
                    run.append(py::make_tuple((int)quads[ri][k],
                                              coords[ri][k].first, coords[ri][k].second));
                out.append(run);
            }
            return out;
        },
        py::arg("u0"), py::arg("axis"), py::arg("a"), py::arg("b"),
        "The rolling-square strip the straight line (a,b) from u0 crosses: a list of "
        "runs, each a list of (vertex_id, eis_x, eis_y). Requires a,b >= 1.");

    cls.def("__repr__", [](PyFD& w) {
        return "<FullereneDual Nv=" + std::to_string(w.N()) +
               " (C" + std::to_string(carbon_count(w.N())) + ")>";
    });
}
