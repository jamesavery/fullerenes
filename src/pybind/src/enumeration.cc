// enumeration.cc — Python binding for buckygen isomer enumeration.
//
// `buckygen(N)` returns a lazy iterator yielding every size-N fullerene as a
// FullereneDual (buckygen emits dual triangulations). Serial enumeration via a
// buckygen_queue; the GIL is released around the blocking next_fullerene call.
// The child process is stopped on context exit, close(), or GC.
//
// buckygen is statically linked into libfullerenes.so (buckygen_main is a
// defined symbol there), so no extra link libraries are needed.

#include <pybind11/pybind11.h>

#include <memory>
#include <utility>

#include "common.hh"

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"

namespace py = pybind11;

using PyFD = pyf::PyGraph<FullereneDual, FullereneDualView>;

namespace {

struct PyBuckyGen {
    BuckyGen::buckygen_queue Q;
    bool open_ = false;

    PyBuckyGen(int N, bool IPR, bool only_nontrivial) {
        Q = BuckyGen::start(N, IPR, only_nontrivial);
        open_ = true;
    }
    PyBuckyGen(const PyBuckyGen&) = delete;
    PyBuckyGen& operator=(const PyBuckyGen&) = delete;
    ~PyBuckyGen() { close(); }

    void close() {
        if (open_) { BuckyGen::stop(Q); open_ = false; }
    }

    PyFD next() {
        if (!open_) throw py::stop_iteration();
        FullereneDual fd;
        bool ok;
        {
            py::gil_scoped_release nogil;          // block without holding the GIL
            ok = BuckyGen::next_fullerene(Q, fd);  // template overload -> dual
        }
        if (!ok) { close(); throw py::stop_iteration(); }
        return PyFD::from_owned(std::move(fd));
    }
};

}  // namespace

void register_enumeration(py::module_& m) {
    py::class_<PyBuckyGen>(m, "BuckyGenerator",
        "Lazy iterator over the fullerene isomers of one size (buckygen), "
        "yielding FullereneDual. Usable as an iterator or a context manager.")
        .def("__iter__", [](py::object self) { return self; })
        .def("__next__", &PyBuckyGen::next)
        .def("__enter__", [](py::object self) { return self; })
        .def("__exit__", [](PyBuckyGen& w, py::object, py::object, py::object) {
            w.close();
            return false;
        })
        .def("close", &PyBuckyGen::close,
             "Stop the generator early (also done on context exit / GC).");

    m.def("buckygen", [](int N, bool IPR, bool only_nontrivial) {
        return std::make_unique<PyBuckyGen>(N, IPR, only_nontrivial);
    }, py::arg("N"), py::arg("IPR") = false, py::arg("only_nontrivial") = false,
       "Iterator over size-N fullerene isomers as FullereneDual. Serial buckygen "
       "enumeration; the GIL is released during generation.");
}
