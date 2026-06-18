#pragma once
//
// common.hh — hand-written zero-copy core: the Py* wrapper.
//
// PyGraph<OwnedT, ViewT> backs an adjacency-only Python graph class. It is the
// single source of truth and is in one of two modes:
//
//   Owned  — holds an OwnedT (Owned<View>) the library produced. .neighbours /
//            .deg are zero-copy numpy views into its vectors (base = wrapper).
//   View   — holds the caller's numpy arrays (keepalive) and rebuilds a bare
//            ViewT over their data pointers on demand.
//
// Every operation runs by materialising the transient trivially-copyable ViewT
// via view() and calling the real library method — never a reimplementation.
//
// OwnedT must be IS-A ViewT (e.g. FullereneGraph : ... : FullereneGraphView).

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <span>
#include <utility>
#include <cstdint>
#include <cstdio>
#include <dlfcn.h>
#include <fcntl.h>
#include <unistd.h>

#include "np_interop.hh"

namespace py = pybind11;

namespace pyf {

// Flush libgfortran's stdout unit (6) to the underlying fd. The Fortran
// force-field optimizer buffers its log in libgfortran's own userspace buffer,
// which fflush(stdout) does NOT touch and which otherwise flushes at process
// exit -- after any fd redirect is gone. Looked up via dlsym so it degrades
// gracefully if the symbol is absent.
inline void flush_fortran_stdout() {
    using flush_fn = void (*)(int*);
    if (auto f = reinterpret_cast<flush_fn>(::dlsym(RTLD_DEFAULT, "_gfortran_flush_i4"))) {
        int u6 = 6;
        f(&u6);
    }
}

// Scoped silencer for the chatty Fortran force-field optimizer, which writes its
// iteration log to stdout (libgfortran unit 6). Redirects ONLY fd 1 to /dev/null
// for its lifetime; fd 2 (stderr) is left alone so genuine errors -- Fortran
// STOP, aborts, asserts, C++ what() -- still surface even when silenced. On
// teardown it first flushes libgfortran's unit-6 buffer (so the buffered log
// goes to /dev/null, not the restored terminal), then restores fd 1. If any of
// open/dup fails it bails without redirecting (never leaves stdout pointing at
// /dev/null). Only safe while the GIL is held (a bound call that does not
// release it). A no-op when `silence` is false.
struct FdSilencer {
    int saved1 = -1;
    explicit FdSilencer(bool silence) {
        if (!silence) return;
        std::fflush(stdout);
        int devnull = ::open("/dev/null", O_WRONLY);
        if (devnull < 0) return;
        int s = ::dup(1);
        if (s < 0) { ::close(devnull); return; }     // can't save -> don't redirect
        if (::dup2(devnull, 1) < 0) { ::close(devnull); ::close(s); return; }
        ::close(devnull);
        saved1 = s;
    }
    ~FdSilencer() {
        if (saved1 < 0) return;
        flush_fortran_stdout();          // push libgfortran's unit-6 buffer to /dev/null
        std::fflush(stdout);
        ::dup2(saved1, 1);
        ::close(saved1);
    }
    FdSilencer(const FdSilencer&) = delete;
    FdSilencer& operator=(const FdSilencer&) = delete;
};

// Validate a (neighbours, deg) pair for View-mode wrapping and return {N, dmax}.
// neighbours must be int32 (N, dmax), deg uint8 (N,), both C-contiguous +
// writeable (exact dtype, no implicit conversion -- else the spans would alias a
// converted copy, breaking zero-copy), with every degree fitting the row stride
// (else operator[] reads past the row). Out-of-spec input raises.
inline std::pair<int, int> validate_adjacency(const py::array& neighbours,
                                              const py::array& deg) {
    require_array(neighbours, "neighbours", 'i', 4, 2, "int32");
    require_array(deg, "deg", 'u', 1, 1, "uint8");
    const int N = (int)neighbours.shape(0), dmax = (int)neighbours.shape(1);
    if ((int)deg.shape(0) != N)
        throw std::invalid_argument("deg length must equal N (neighbours.shape[0])");
    const uint8_t* dptr = static_cast<const uint8_t*>(deg.data());
    for (int u = 0; u < N; ++u)
        if (dptr[u] > dmax)
            throw std::invalid_argument("deg[" + std::to_string(u) + "]=" +
                std::to_string(dptr[u]) + " exceeds dmax=" + std::to_string(dmax));
    return {N, dmax};
}

template<typename OwnedT, typename ViewT>
struct PyGraph {
    enum class Mode { Owned, View };
    Mode mode = Mode::Owned;

    OwnedT owned;                      // valid in Owned mode

    // View mode: keepalive refs to the caller's numpy arrays (validated to be
    // exactly int32/uint8, C-contiguous, writeable -- so the spans below alias
    // the caller's buffer rather than a converted copy) + dimensions.
    py::array a_neighbours, a_deg;
    int Nv = 0, dmaxv = 0;

    PyGraph() = default;

    static PyGraph from_owned(OwnedT o) {
        PyGraph w;
        w.mode = Mode::Owned;
        w.owned = std::move(o);
        return w;
    }

    // BYO-buffer: wrap caller numpy arrays WITHOUT copying. validate_adjacency
    // enforces dtype/shape/contiguity/degree up front (out-of-spec raises).
    static PyGraph from_arrays(py::array neighbours, py::array deg) {
        auto [Nv, dmaxv] = validate_adjacency(neighbours, deg);
        PyGraph w;
        w.mode = Mode::View;
        w.Nv = Nv;
        w.dmaxv = dmaxv;
        w.a_neighbours = std::move(neighbours);
        w.a_deg = std::move(deg);
        return w;
    }

    int N()    const { return mode == Mode::Owned ? (int)owned.N    : Nv; }
    int dmax() const { return mode == Mode::Owned ? (int)owned.dmax : dmaxv; }

    // In View mode the arrays were validated int32/uint8 + C-contiguous at
    // from_arrays, so .mutable_data() is a plain pointer into the caller's
    // buffer -- no per-call cast or conversion.
    node_t* neighbours_data() {
        return mode == Mode::Owned
            ? owned.neighbours.data()
            : static_cast<node_t*>(a_neighbours.mutable_data());
    }
    uint8_t* deg_data() {
        return mode == Mode::Owned
            ? owned.deg.data()
            : static_cast<uint8_t*>(a_deg.mutable_data());
    }

    // Transient trivially-copyable view over whichever storage is live.
    ViewT view() {
        if (mode == Mode::Owned)
            return static_cast<ViewT&>(owned);   // OwnedT IS-A ViewT (slice to base)
        std::span<node_t>  sn(neighbours_data(), (size_t)Nv * dmaxv);
        std::span<uint8_t> sd(deg_data(), (size_t)Nv);
        return ViewT((node_t)Nv, dmaxv, sn, sd);
    }
};

// Snapshot of the last optimize()/from_dual run on a geometry wrapper. Primitive
// fields only -- the result is kept as an int code -- so this generic core stays
// free of any algorithm's enum (the typed accessor casts the code back).
// has_run is false until a run populates it: reading a stat before then is an
// error, not a fabricated default.
struct OptStats {
    bool   has_run = false;
    int    iterations_used = 0;
    double final_gmax_L = 0.0;
    double final_angle_relerr = 0.0;
    int    final_n_concave = 0;
    int    final_opt_result_code = -1;
};

// ---------------------------------------------------------------------------
// PyGeom<OwnedT, ViewT>: a PyGraph that also carries coord3d `points`
// (Polyhedron, Deltahedron). It adds only the points array, points_data(), a
// view() that includes the points span, and an opt_stats snapshot of the last
// optimization; mode/owned/the adjacency arrays and accessors/validation are
// inherited. `.points` is a zero-copy WRITEABLE view, so in-place optimize()/
// edits write straight back.
// ---------------------------------------------------------------------------
template<typename OwnedT, typename ViewT>
struct PyGeom : PyGraph<OwnedT, ViewT> {
    using Base = PyGraph<OwnedT, ViewT>;
    using typename Base::Mode;          // the {Owned, View} enum

    py::array a_points;                 // View mode (validated (N,3) float64, keepalive)
    OptStats  opt_stats;                // last optimize()/from_dual outcome (see getters)

    static PyGeom from_owned(OwnedT o) {
        PyGeom w;
        w.mode = Mode::Owned;
        w.owned = std::move(o);
        return w;
    }

    static PyGeom from_arrays(py::array neighbours, py::array deg, py::array points) {
        auto [N, dmax] = validate_adjacency(neighbours, deg);
        require_array(points, "points", 'f', 8, 2, "float64");
        if ((int)points.shape(0) != N || (int)points.shape(1) != 3)
            throw std::invalid_argument("points must have shape (N, 3)");
        PyGeom w;
        w.mode = Mode::View;
        w.Nv = N;
        w.dmaxv = dmax;
        w.a_neighbours = std::move(neighbours);
        w.a_deg = std::move(deg);
        w.a_points = std::move(points);
        return w;
    }

    coord3d* points_data() {
        return this->mode == Mode::Owned ? this->owned.points.data()
                                         : static_cast<coord3d*>(a_points.mutable_data());
    }

    // View carrying the points span (hides the adjacency-only Base::view()).
    ViewT view() {
        if (this->mode == Mode::Owned)
            return static_cast<ViewT&>(this->owned);
        std::span<node_t>  sn(this->neighbours_data(), (size_t)this->Nv * this->dmaxv);
        std::span<uint8_t> sd(this->deg_data(), (size_t)this->Nv);
        std::span<coord3d> sp(points_data(), (size_t)this->Nv);
        return ViewT((node_t)this->Nv, this->dmaxv, sn, sd, sp);
    }
};

// Register the accessors shared by every adjacency graph class on an existing
// py::class_. Class-specific methods (dual, transforms, naming) are added by
// the per-class register_* functions (and, later, by the generator).
// The accessors every wrapper (graph or geom) shares: dims, zero-copy read-only
// neighbours/deg views (base = the wrapper), and the ragged adjacency copy.
// Templated on the wrapper W directly, so w.view() resolves to PyGraph's
// adjacency-only view or PyGeom's points-carrying override as appropriate.
template<typename W, typename PyClass>
void bind_adjacency_accessors(PyClass& cls) {
    cls.def_property_readonly("N",    [](W& w) { return w.N(); });
    cls.def_property_readonly("dmax", [](W& w) { return w.dmax(); });
    cls.def_property_readonly("neighbours", [](py::object self) {
        W& w = self.cast<W&>();
        return neighbours_view(self, w.neighbours_data(), w.N(), w.dmax());
    });
    cls.def_property_readonly("deg", [](py::object self) {
        W& w = self.cast<W&>();
        return deg_view(self, w.deg_data(), w.N());
    });
    cls.def("adjacency", [](W& w) { return adjacency_list(w.view()); },
            "Neighbour lists as list[list[int]] (a copy; padding removed).");
}

// Graph classes: the shared adjacency accessors + is_a_fullerene.
template<typename OwnedT, typename ViewT, typename PyClass>
void bind_graph_common(PyClass& cls) {
    using W = PyGraph<OwnedT, ViewT>;
    bind_adjacency_accessors<W>(cls);
    cls.def("is_a_fullerene", [](W& w) { return w.view().is_a_fullerene(false); });
}

// Geometry classes: the shared adjacency accessors + a zero-copy WRITEABLE
// (N,3) float64 `points` view (in-place edits / optimize write back).
template<typename OwnedT, typename ViewT, typename PyClass>
void bind_geom_common(PyClass& cls) {
    using W = PyGeom<OwnedT, ViewT>;
    bind_adjacency_accessors<W>(cls);
    cls.def_property_readonly("points", [](py::object self) {
        W& w = self.cast<W&>();
        return points_view(self, w.points_data(), w.N());
    });
}

}  // namespace pyf
