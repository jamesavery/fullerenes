#include <CL/sycl.hpp>
#include <numeric>
#include <algorithm>
#include <parallel/algorithm>
#include <execution>
#include <functional>
#include <fullerenes/sycl-wrappers.hh>

namespace primitives{
    template <typename T> T reduce(sycl::queue& Q, const sycl::span<T> vec);
    template <typename T> T reduce(sycl::queue& Q, const sycl::buffer<T>& buf);
    template <typename T> T reduce(sycl::queue& Q, const SyclVector<T>& vec);

    template <typename T, typename UnaryOp> void transform(sycl::queue& Q, T* begin, T* end, T* result, UnaryOp f);
    template <typename T, typename UnaryOp> void transform(sycl::queue& Q, SyclVector<T>& vec, SyclVector<T>& result, UnaryOp f);
    template <typename T, typename UnaryOp> void transform(sycl::queue& Q, Span<T> vec, Span<T> result, UnaryOp f);

    template <typename T, typename BinaryOp, typename UnaryOp> void transform_reduce(sycl::queue& Q, SyclVector<T>& vec, BinaryOp op, UnaryOp f);
    template <typename T, typename BinaryOp, typename UnaryOp> void transform_reduce(sycl::queue& Q, sycl::buffer<T>& buf, BinaryOp op, UnaryOp f);
    template <typename T, typename BinaryOp, typename UnaryOp> void transform_reduce(sycl::queue& Q, sycl::span<T> vec, BinaryOp op, UnaryOp f);

    template <typename T, typename BinaryOp> void exclusive_scan(sycl::queue& Q, SyclVector<T>& vec, BinaryOp op);
    template <typename T, typename BinaryOp> void exclusive_scan(sycl::queue& Q, sycl::buffer<T>& buf, BinaryOp op);
    template <typename T, typename BinaryOp> void exclusive_scan(sycl::queue& Q, sycl::span<T> vec, BinaryOp op);

    template <typename T, typename BinaryOp, typename UnaryOp> void transform_exclusive_scan(sycl::queue& Q, T* begin, T* end, T* store, T init, BinaryOp op, UnaryOp f);

    template <typename T, typename BinaryOp> void inclusive_scan(sycl::queue& Q, SyclVector<T>& vec, BinaryOp op);
    template <typename T, typename BinaryOp> void inclusive_scan(sycl::queue& Q, sycl::buffer<T>& buf, BinaryOp op);
    template <typename T, typename BinaryOp> void inclusive_scan(sycl::queue& Q, sycl::span<T> vec, BinaryOp op);

    template <typename T, typename BinaryOp, typename UnaryOp> void transform_inclusive_scan(sycl::queue& Q, T* begin, T* end, T* store, T init, BinaryOp op,  UnaryOp f);
    template <typename T, typename BinaryOp, typename UnaryOp> void transform_inclusive_scan(sycl::queue& Q, SyclVector<T>& vec, SyclVector<T>& store, T init, BinaryOp op, UnaryOp f);
    template <typename T, typename BinaryOp, typename UnaryOp> void transform_inclusive_scan(sycl::queue& Q, Span<T> vec, Span<T> store, T init, BinaryOp op, UnaryOp f);

}