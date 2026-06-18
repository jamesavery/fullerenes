// deltahedron.cc — Python bindings for Deltahedron (dual triangulation with an
// equilateral-triangle 3D embedding).
//
// The usual source is the dual->deltahedron pipeline (from_dual), which returns
// an owned, already-optimized Deltahedron with zero-copy .points. View-level ops
// (smooth/reflect) write through the points span. optimize() writes into
// .points in owned mode (zero-copy); in view mode it optimizes a temporary and
// copies the result back into the caller's buffer (Deltahedron::optimize writes
// to owned storage, not an arbitrary span).

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <span>
#include <vector>

#include "common.hh"

#include "fullerenes/deltahedron.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/graph.hh"
#include "fullerenes/buckinverse.hh"

namespace py = pybind11;

using PyDelta = pyf::PyGeom<Deltahedron, DeltahedronView<double>>;
using PyFD = pyf::PyGraph<FullereneDual, FullereneDualView>;

// Mechanical method surface (max_angle_relerr/count_concave/triangles/... and the
// triangulation/graph queries), emitted by tools/gen_bindings.py.
void register_generated_Deltahedron(py::class_<PyDelta>& cls);

// Snapshot a Deltahedron's optimizer-stat fields into the wrapper's OptStats.
static pyf::OptStats capture(const Deltahedron& D) {
    return { /*has_run=*/true, D.iterations_used, D.final_gmax_L,
             D.final_angle_relerr, D.final_n_concave, (int)D.final_opt_result };
}

// The last-optimize snapshot, or raise if no optimize()/from_dual has run --
// rather than returning a fabricated default verdict.
static const pyf::OptStats& opt_stats(PyDelta& w) {
    if (!w.opt_stats.has_run)
        throw std::runtime_error("no optimize()/from_dual has run on this Deltahedron");
    return w.opt_stats;
}

void register_deltahedron(py::module_& m) {
    py::enum_<OptResult>(m, "OptResult", "Why the deltahedron optimizer stopped.")
        .value("CONVERGED", OptResult::CONVERGED)
        .value("STAGNATED", OptResult::STAGNATED)
        .value("BUDGET_EXHAUSTED", OptResult::BUDGET_EXHAUSTED)
        .value("CONVEXITY_STUCK", OptResult::CONVEXITY_STUCK);

    py::class_<PyDelta> cls(m, "Deltahedron",
        "Dual triangulation with an equilateral-triangle 3D embedding. "
        ".points is a zero-copy, writeable numpy view.");

    pyf::bind_geom_common<Deltahedron, DeltahedronView<double>>(cls);
    register_generated_Deltahedron(cls);  // max_angle_relerr/count_concave/triangles/...

    // --- Construction ---
    cls.def_static("from_dual", [](PyFD& dual) {
        Graph g(dual.view());                              // ReducibleDual takes an owned Graph
        buckinverse::ReducibleDual rd(g);
        buckinverse::ExtensionPath ep = rd.reduceToExtensionPath(20);
        Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);
        // Post-condition: the reduce/expand round-trip must reproduce the input
        // dual exactly -- guard rather than silently return a different-sized
        // embedding if reduction ever fails to reach a seed.
        if ((int)D.N != (int)g.N)
            throw std::runtime_error(
                "from_dual: reduction/expansion did not reproduce the input dual "
                "(got Nv=" + std::to_string((int)D.N) +
                ", expected " + std::to_string((int)g.N) + ")");
        PyDelta w = PyDelta::from_owned(std::move(D));
        w.opt_stats = capture(w.owned);    // the pipeline already optimized
        return w;
    }, py::arg("dual"),
       "Equilateral 3D embedding of a fullerene dual via the reduction/expansion "
       "pipeline (already optimized).");

    cls.def_static("from_arrays",
        [](py::array neighbours, py::array points, py::object deg_obj) {
            py::array deg = deg_obj.is_none()
                ? py::array(pyf::infer_deg(neighbours))    // 5/6 from the -1 padding
                : deg_obj.cast<py::array>();
            return PyDelta::from_arrays(std::move(neighbours), std::move(deg), std::move(points));
        },
        py::arg("neighbours"), py::arg("points"), py::arg("deg") = py::none(),
        "Wrap caller-owned int32 (N,6) + float64 (N,3) arrays (zero-copy). deg is "
        "inferred from the -1 padding when omitted. Same arg order as "
        "Polyhedron.from_arrays. NOTE: the topology is trusted to be a valid "
        "triangulation; it is not structurally validated.");

    // (max_angle_relerr / count_concave / triangles / n_triangles are generated.)

    // --- In-place geometry (write through .points) ---
    cls.def("smooth", [](PyDelta& w, double q) { w.view().smooth(q); }, py::arg("q"),
            "Laplacian-style smoothing step (in place).");
    cls.def("reflect_concave", [](PyDelta& w, double threshold) {
        auto v = w.view();
        return v.reflect_concave(v.points, threshold);
    }, py::arg("threshold") = 0.0, "Reflect concave vertices once; returns the count reflected.");
    cls.def("reflect_all_concave", [](PyDelta& w, double threshold) {
        auto v = w.view();
        return v.reflect_all_concave(v.points, threshold);
    }, py::arg("threshold") = 0.0, "Reflect concave vertices to convergence; returns total reflected.");

    // --- Optimize (writes into .points; view mode copies back) ---
    cls.def("optimize", [](PyDelta& w, double target_L, double grad_tol,
                           long long max_work, double angle_tol) {
        std::vector<coord3d> init(w.points_data(), w.points_data() + w.N());
        if (w.mode == PyDelta::Mode::Owned) {
            OptResult r = w.owned.optimize(init, target_L, grad_tol, {}, max_work, angle_tol);
            w.opt_stats = capture(w.owned);          // zero-copy: owned IS the geometry
            return r;
        }
        // View mode: optimize a temporary, copy the result back into the caller's
        // buffer, and snapshot the temporary's stats onto the wrapper.
        Deltahedron tmp(w.view(), std::span<const coord3d>(w.points_data(), (size_t)w.N()));
        OptResult r = tmp.optimize(init, target_L, grad_tol, {}, max_work, angle_tol);
        std::copy(tmp.points.begin(), tmp.points.end(), w.points_data());
        w.opt_stats = capture(tmp);
        return r;
    }, py::arg("target_L") = 0.0, py::arg("grad_tol") = 1e-10,
       py::arg("max_work") = 0, py::arg("angle_tol") = 0.0,
       "Optimize the equilateral embedding in place; returns an OptResult.");

    // --- Optimizer stats from the LAST optimize() / from_dual pipeline ---
    // A snapshot on the wrapper (NOT refreshed by smooth()/reflect_*() edits; for
    // live quality use count_concave() / max_angle_relerr()). Reading any of these
    // before an optimize()/from_dual has run raises rather than returning a
    // fabricated default verdict.
    cls.def_property_readonly("iterations_used", [](PyDelta& w) { return opt_stats(w).iterations_used; },
        "Iterations of the last optimize()/from_dual (snapshot; not live).");
    cls.def_property_readonly("final_gmax_L", [](PyDelta& w) { return opt_stats(w).final_gmax_L; },
        "max_i ||g_i||*L at the end of the last optimize() (snapshot; not live).");
    cls.def_property_readonly("final_angle_relerr", [](PyDelta& w) { return opt_stats(w).final_angle_relerr; },
        "Angle rel-error at the last optimize() (snapshot; use max_angle_relerr() for live).");
    cls.def_property_readonly("final_n_concave", [](PyDelta& w) { return opt_stats(w).final_n_concave; },
        "Concave count at the last optimize() (snapshot; use count_concave() for live).");
    cls.def_property_readonly("final_opt_result", [](PyDelta& w) { return OptResult(opt_stats(w).final_opt_result_code); },
        "OptResult of the last optimize()/from_dual.");

    // --- Transforms (return new Deltahedron) ---
    cls.def("gc_transform", [](PyDelta& w, unsigned k, unsigned l) {
        return PyDelta::from_owned(w.view().GCtransform(k, l));
    }, py::arg("k"), py::arg("l") = 0u);
    cls.def("halma", [](PyDelta& w, int m) {
        return PyDelta::from_owned(w.view().halma_transform(m));
    }, py::arg("m"));

    cls.def("__repr__", [](PyDelta& w) {
        return "<Deltahedron Nv=" + std::to_string(w.N()) + ">";
    });
}
