#pragma once
//
// np_interop.hh — hand-written zero-copy core: numpy <-> C++ array makers.
//
// Two kinds of maker:
//   *_view(base, ptr, ...)  — zero-copy numpy view into live C++ storage. The
//                             `base` handle (the Py* wrapper) is kept alive by
//                             the returned array, so the storage cannot vanish
//                             while a view exists.
//   *_copy(...)             — fresh numpy array copied from a returned C++ value
//                             (vector<tri_t>, matrix<int>, ...). One copy, at an
//                             allocation boundary anyway.
//
// Topology views (neighbours/deg) are returned READ-ONLY: mutating adjacency in
// place would violate the no-reallocation invariant. Coordinate views (points)
// are writeable — that is the in-place optimize() write-back path.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>   // vector<vector<int>> -> list[list[int]]

#include <array>
#include <cstring>
#include <string>
#include <vector>
#include <cstdint>

#include "fullerenes/geometry.hh"     // coord3d, coord2d, tri_t, matrix3d
#include "fullerenes/matrix.hh"
#include "fullerenes/dense_graph.hh"  // Spanify::to_neighbours

namespace py = pybind11;

namespace pyf {

inline void set_readonly(py::array& a) {
    py::detail::array_proxy(a.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
}

// View mode aliases the caller's buffer in place, so the array must match the
// expected dtype/ndim EXACTLY (no implicit conversion -- a silent dtype copy
// would break the zero-copy contract: writes would never reach the caller's
// array), be C-contiguous, and be writeable. Fails with a clear message naming
// the required numpy dtype instead of silently copying or a cryptic accessor
// error. `kind`/`itemsize` identify the dtype (e.g. 'i',4 = int32; 'u',1 = uint8).
inline void require_array(const py::array& a, const char* name,
                          char kind, int itemsize, int ndim, const char* npdtype) {
    if (a.dtype().kind() != kind || (int)a.itemsize() != itemsize)
        throw std::invalid_argument(std::string(name) + " must have dtype " + npdtype +
            " (got a different dtype; pass an exact numpy." + npdtype +
            " array -- no implicit conversion, to preserve zero-copy)");
    if (a.ndim() != ndim)
        throw std::invalid_argument(std::string(name) + " must be " +
            std::to_string(ndim) + "-dimensional");
    auto flags = py::detail::array_proxy(a.ptr())->flags;
    if (!(flags & py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_))
        throw std::invalid_argument(std::string(name) + " must be C-contiguous");
    if (!(flags & py::detail::npy_api::NPY_ARRAY_WRITEABLE_))
        throw std::invalid_argument(
            std::string(name) + " must be writeable (view mode shares the buffer in place)");
}

// --- Zero-copy views into live C++ storage (base keeps the owner alive) -----

// neighbours: (N, dmax) int32, row-major. Read-only (padding cols are -1).
inline py::array_t<node_t> neighbours_view(py::handle base, node_t* data, int N, int dmax) {
    py::array_t<node_t> a(
        {(py::ssize_t)N, (py::ssize_t)dmax},
        {(py::ssize_t)(sizeof(node_t) * dmax), (py::ssize_t)sizeof(node_t)},
        data, base);
    set_readonly(a);
    return a;
}

// deg: (N,) uint8. Read-only.
inline py::array_t<uint8_t> deg_view(py::handle base, uint8_t* data, int N) {
    py::array_t<uint8_t> a({(py::ssize_t)N}, {(py::ssize_t)sizeof(uint8_t)}, data, base);
    set_readonly(a);
    return a;
}

// points: (N, 3) float64. Writeable — in-place geometry edits write straight back.
inline py::array_t<double> points_view(py::handle base, coord3d* data, int N) {
    return py::array_t<double>(
        {(py::ssize_t)N, (py::ssize_t)3},
        {(py::ssize_t)sizeof(coord3d), (py::ssize_t)sizeof(double)},
        reinterpret_cast<double*>(data), base);
}

// --- Copies from returned C++ values ----------------------------------------

inline py::array_t<int> tris_copy(const std::vector<tri_t>& tris) {
    py::array_t<int> a({(py::ssize_t)tris.size(), (py::ssize_t)3});
    auto r = a.mutable_unchecked<2>();
    for (py::ssize_t i = 0; i < (py::ssize_t)tris.size(); ++i) {
        r(i, 0) = tris[i][0]; r(i, 1) = tris[i][1]; r(i, 2) = tris[i][2];
    }
    return a;
}

// face_t is `struct face_t : vector<node_t>` (not a typedef), so pybind/stl
// can't auto-convert it; copy each face to a plain vector<int> -> list[list[int]].
inline std::vector<std::vector<int>> faces_copy(const std::vector<face_t>& faces) {
    std::vector<std::vector<int>> out;
    out.reserve(faces.size());
    for (const auto& f : faces) out.emplace_back(f.begin(), f.end());
    return out;
}

// The zero-copy points view and the coord/matrix makers below memcpy/stride on
// the assumption that coord3d is exactly 3 tightly-packed doubles and matrix3d
// 9 -- guard it so a future layout change fails loudly instead of silently
// mis-indexing coordinates.
static_assert(sizeof(coord3d) == 3 * sizeof(double),
    "points_view / coord3d_copy assume coord3d is 3 packed doubles");
static_assert(sizeof(matrix3d) == 9 * sizeof(double),
    "matrix3d_copy assumes matrix3d is 9 packed doubles");

// Infer per-vertex degrees from a (N,dmax) int32 neighbour array's -1 padding
// (RSR convention: valid neighbours at the front, node_t(-1) padding after).
// Works for any degree -- unlike a fixed cubic default, which silently corrupts
// a non-cubic graph. Validates neighbours first so .data() is safe to read.
inline py::array_t<uint8_t> infer_deg(const py::array& neighbours) {
    require_array(neighbours, "neighbours", 'i', 4, 2, "int32");
    const py::ssize_t N = neighbours.shape(0), dmax = neighbours.shape(1);
    const node_t* p = static_cast<const node_t*>(neighbours.data());
    py::array_t<uint8_t> d(std::vector<py::ssize_t>{N});
    auto r = d.mutable_unchecked<1>();
    for (py::ssize_t u = 0; u < N; ++u) {
        int cnt = 0;
        for (py::ssize_t j = 0; j < dmax; ++j)
            if (p[u * dmax + j] != node_t(-1)) ++cnt;
        r(u) = (uint8_t)cnt;
    }
    return d;
}

// coord3d -> (3,) float64 ; matrix3d -> (3,3) float64 (both contiguous, memcpy-safe).
inline py::array_t<double> coord3d_copy(const coord3d& c) {
    py::array_t<double> a(3);
    std::memcpy(a.mutable_data(), &c, 3 * sizeof(double));
    return a;
}

inline py::array_t<double> matrix3d_copy(const matrix3d& m) {
    py::array_t<double> a({(py::ssize_t)3, (py::ssize_t)3});
    std::memcpy(a.mutable_data(), m.values, 9 * sizeof(double));
    return a;
}

// Read a length-3 sequence (list/tuple/ndarray) into a coord3d.
inline coord3d as_coord3d(py::handle h) {
    auto s = h.cast<std::array<double, 3>>();
    return coord3d(s[0], s[1], s[2]);
}

template<typename T>
inline py::array_t<T> matrix_copy(const matrix<T>& M) {
    py::array_t<T> a({(py::ssize_t)M.m, (py::ssize_t)M.n});
    auto r = a.template mutable_unchecked<2>();
    for (py::ssize_t i = 0; i < (py::ssize_t)M.m; ++i)
        for (py::ssize_t j = 0; j < (py::ssize_t)M.n; ++j)
            r(i, j) = M(i, j);
    return a;
}

// Ragged adjacency: list[list[int]], one row per vertex, padding stripped.
// Delegates to the library's Spanify::to_neighbours (dense_graph.hh) so the
// RSR->ragged convention lives in exactly one place; pybind11/stl converts the
// returned vector<vector<int>> to list[list[int]].
template<typename View>
inline std::vector<std::vector<int>> adjacency_list(const View& v) {
    return Spanify::to_neighbours(v);
}

}  // namespace pyf
