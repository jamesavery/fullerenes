#include "queue-impl.hh"
#include <fullerenes/sycl-headers/execution-compat.hh>
#include <algorithm>
#include <numeric>
#include <fullerenes/kernel-headers/geometry.hh>
#include "forcefield-includes.hh"

template <typename T>
symMat3<T> inertia_matrix(sycl::group<1>& cta, const std::span<std::array<T,3>> X){
    auto tid = cta.get_local_id(0);
    symMat3<T> I;
    T diag = sycl::reduce_over_group(cta, dot(X[tid], X[tid]), sycl::plus<T>());
    I.a = diag;
    I.d = diag;
    I.f = diag;
    I.a -= sycl::reduce_over_group(cta, X[tid][0]*X[tid][0], sycl::plus<T>());
    I.b -= sycl::reduce_over_group(cta, X[tid][0]*X[tid][1], sycl::plus<T>());
    I.c -= sycl::reduce_over_group(cta, X[tid][0]*X[tid][2], sycl::plus<T>());
    I.d -= sycl::reduce_over_group(cta, X[tid][1]*X[tid][1], sycl::plus<T>());
    I.e -= sycl::reduce_over_group(cta, X[tid][1]*X[tid][2], sycl::plus<T>());
    I.f -= sycl::reduce_over_group(cta, X[tid][2]*X[tid][2], sycl::plus<T>());
    return I;
}

template <typename T>
auto inertia_matrix(SyclQueue& Q, const std::span<std::array<T,3>> X){
    symMat3<T> I;
    Q.wait();
    T diag = std::transform_reduce(FULLERENE_PAR_UNSEQ X.begin(), X.end(), T(0), std::plus<T>{}, [](const auto& x) -> T {return dot(x,x);});
    I.a = diag;
    I.d = diag;
    I.f = diag;
    I.a -= std::transform_reduce(FULLERENE_PAR_UNSEQ X.begin(), X.end(), T(0), std::plus<T>{}, [](const auto& x){return x[0]*x[0];});
    I.b -= std::transform_reduce(FULLERENE_PAR_UNSEQ X.begin(), X.end(), T(0), std::plus<T>{}, [](const auto& x){return x[0]*x[1];});
    I.c -= std::transform_reduce(FULLERENE_PAR_UNSEQ X.begin(), X.end(), T(0), std::plus<T>{}, [](const auto& x){return x[0]*x[2];});
    I.d -= std::transform_reduce(FULLERENE_PAR_UNSEQ X.begin(), X.end(), T(0), std::plus<T>{}, [](const auto& x){return x[1]*x[1];});
    I.e -= std::transform_reduce(FULLERENE_PAR_UNSEQ X.begin(), X.end(), T(0), std::plus<T>{}, [](const auto& x){return x[1]*x[2];});
    I.f -= std::transform_reduce(FULLERENE_PAR_UNSEQ X.begin(), X.end(), T(0), std::plus<T>{}, [](const auto& x){return x[2]*x[2];});
    return I;
}

template <typename T>
auto principal_axes(SyclQueue& Q, const std::span<std::array<T,3>> X){
    auto I = inertia_matrix(Q, X);
    auto [V,lambdas] = I.eigensystem();
    return V;
}

template <typename T>
std::array<std::array<T,3>,3> principal_axes(sycl::group<1>& cta, const std::span<std::array<T,3>> X){
    auto I = inertia_matrix(cta,X);
    auto [V,lambdas] = I.eigensystem();
    return V;
}

// ---------------------------------------------------------------------------
// View-based batch implementations (Phase 7)
// ---------------------------------------------------------------------------

// Forward-declare kernel tag structs used by SYCL parallel_for.
template<typename T, typename K> struct EccentricityFunctorView;
template<typename T, typename K> struct InertiaFunctorView;
template<typename T, typename K> struct TransformCoordinatesFunctorView;
template<typename T, typename K> struct SurfaceAreaFunctorView;
template<typename T, typename K> struct VolumeFunctorView;

template <typename T, typename K>
SyclEvent eccentricity(SyclQueue& Q, std::span<std::array<T,3>> xyz, int N, int capacity,
                       batch::BatchStateView state, std::span<T> out_ellipticity,
                       Workspace /*ws*/) {
    auto statuses = state.status;
    // Empty batch: a zero-size nd_range launch is a no-op on the host/OpenMP backend
    // but cuLaunchKernel rejects a zero grid with CUDA_ERROR_INVALID_VALUE on CUDA.
    if (capacity == 0 || N == 0) return SyclEvent();
    SyclEventImpl ret_val = Q->submit([=](sycl::handler& cgh) {
        cgh.parallel_for<struct EccentricityFunctorView<T,K>>(
            sycl::nd_range<1>(sycl::range<1>(capacity*N), sycl::range<1>(N)),
            [=](sycl::nd_item<1> nditem) {
                auto cta = nditem.get_group();
                auto bid = cta.get_group_linear_id();
                if (statuses[bid].is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;
                auto X = xyz.subspan(bid * N, N);
                auto I = inertia_matrix(cta, X);
                auto [V, lambdas] = I.eigensystem();
                auto elipsoid = rsqrt3(d_sort(d_abs(lambdas)));
                out_ellipticity[bid] = elipsoid[0] / elipsoid[2];
            });
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent inertia(SyclQueue& Q, std::span<std::array<T,3>> xyz, int N, int capacity,
                  batch::BatchStateView state, std::span<std::array<T,3>> out_inertia,
                  Workspace /*ws*/) {
    auto statuses = state.status;
    // Empty batch: a zero-size nd_range launch is a no-op on the host/OpenMP backend
    // but cuLaunchKernel rejects a zero grid with CUDA_ERROR_INVALID_VALUE on CUDA.
    if (capacity == 0 || N == 0) return SyclEvent();
    SyclEventImpl ret_val = Q->submit([=](sycl::handler& cgh) {
        cgh.parallel_for<struct InertiaFunctorView<T,K>>(
            sycl::nd_range<1>(sycl::range<1>(capacity*N), sycl::range<1>(N)),
            [=](sycl::nd_item<1> nditem) {
                auto cta = nditem.get_group();
                auto tid = cta.get_local_linear_id();
                auto bid = cta.get_group_linear_id();
                if (statuses[bid].is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;
                auto X = xyz.subspan(bid * N, N);
                auto I = inertia_matrix(cta, X);
                if (tid == 0) out_inertia[bid] = I.eigenvalues();
            });
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent transform_to_principal_axes(SyclQueue& Q, std::span<std::array<T,3>> xyz, int N, int capacity,
                                                     batch::BatchStateView state,
                                                     Workspace /*ws*/) {
    auto statuses = state.status;
    // Empty batch: a zero-size nd_range launch is a no-op on the host/OpenMP backend
    // but cuLaunchKernel rejects a zero grid with CUDA_ERROR_INVALID_VALUE on CUDA.
    if (capacity == 0 || N == 0) return SyclEvent();
    SyclEventImpl ret_val = Q->submit([=](sycl::handler& cgh) {
        cgh.parallel_for<struct TransformCoordinatesFunctorView<T,K>>(
            sycl::nd_range<1>(sycl::range<1>(capacity*N), sycl::range<1>(N)),
            [=](sycl::nd_item<1> nditem) {
                auto cta = nditem.get_group();
                auto tid = cta.get_local_linear_id();
                auto bid = cta.get_group_linear_id();
                if (statuses[bid].is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;
                auto X = xyz.subspan(bid * N, N);
                auto P = principal_axes(cta, X);
                if (sycl::isnan(P[0][0]) || sycl::isnan(P[0][1]) || sycl::isnan(P[0][2]) ||
                    sycl::isnan(P[1][0]) || sycl::isnan(P[1][1]) || sycl::isnan(P[1][2]) ||
                    sycl::isnan(P[2][0]) || sycl::isnan(P[2][1]) || sycl::isnan(P[2][2])) return;
                X[tid] = dot(P, X[tid]);
            });
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent surface_area(SyclQueue& Q, std::span<std::array<T,3>> xyz, int N, int capacity,
                                            std::span<std::array<K,6>> faces, std::span<uint8_t> deg,
                                            batch::BatchStateView state, std::span<T> out_surface_area,
                                            Workspace /*ws*/) {
    FLOAT_TYPEDEFS(T);
    auto statuses = state.status;
    // Empty batch: a zero-size nd_range launch is a no-op on the host/OpenMP backend
    // but cuLaunchKernel rejects a zero grid with CUDA_ERROR_INVALID_VALUE on CUDA.
    if (capacity == 0 || N == 0) return SyclEvent();
    const int Nf = N / 2 + 2;
    SyclEventImpl ret_val = Q->submit([=](sycl::handler& cgh) {
        auto X_smem2 = sycl::local_accessor<std::array<T,3>, 1>(N, cgh);
        cgh.parallel_for<struct SurfaceAreaFunctorView<T,K>>(
            sycl::nd_range<1>(sycl::range<1>(capacity*Nf), sycl::range<1>(Nf)),
            [=](sycl::nd_item<1> nditem) {
                auto cta = nditem.get_group();
                auto tid = cta.get_local_linear_id();
                auto bid = cta.get_group_linear_id();
                if (statuses[bid].is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;
                auto X = xyz.subspan(bid * N, N);
                for (int i = tid; i < N; i += cta.get_local_range(0)) X_smem2[i] = X[i];
                sycl::group_barrier(cta);
                T A = 0;
                coord3d face_center = {0,0,0};
                auto face = faces[bid * Nf + tid];
                auto face_size = deg[bid * Nf + tid];
                for (int i = 0; i < face_size; i++) face_center += X_smem2[face[i]];
                face_center /= T(face_size);
                for (int i = 0; i < face_size; i++) {
                    auto a = X_smem2[face[i]];
                    auto b = X_smem2[face[(i+1) % face_size]];
                    auto c = face_center;
                    auto u = b - a; auto v = c - a;
                    A += norm(cross(u, v));
                }
                auto result = sycl::reduce_over_group(cta, A, sycl::plus<T>()) / T(2);
                if (tid == 0) out_surface_area[bid] = result;
            });
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent volume(SyclQueue& Q, std::span<std::array<T,3>> xyz, int N, int capacity,
                                       std::span<std::array<K,6>> faces, std::span<uint8_t> deg,
                                       batch::BatchStateView state, std::span<T> out_volume,
                                       Workspace /*ws*/) {
    FLOAT_TYPEDEFS(T);
    auto statuses = state.status;
    // Empty batch: a zero-size nd_range launch is a no-op on the host/OpenMP backend
    // but cuLaunchKernel rejects a zero grid with CUDA_ERROR_INVALID_VALUE on CUDA.
    if (capacity == 0 || N == 0) return SyclEvent();
    const int Nf = N / 2 + 2;
    SyclEventImpl ret_val = Q->submit([=](sycl::handler& cgh) {
        auto X_smem = sycl::local_accessor<std::array<T,3>, 1>(N, cgh);
        cgh.parallel_for<struct VolumeFunctorView<T,K>>(
            sycl::nd_range<1>(sycl::range<1>(capacity*Nf), sycl::range<1>(Nf)),
            [=](sycl::nd_item<1> nditem) {
                auto cta = nditem.get_group();
                auto tid = cta.get_local_linear_id();
                auto bid = cta.get_group_linear_id();
                if (statuses[bid].is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;
                auto X = xyz.subspan(bid * N, N);
                for (int i = tid; i < N; i += cta.get_local_range(0)) X_smem[i] = X[i];
                sycl::group_barrier(cta);
                T V = 0;
                coord3d face_center = {0,0,0};
                auto face = faces[bid * Nf + tid];
                auto face_size = deg[bid * Nf + tid];
                for (int i = 0; i < face_size; i++) face_center += X_smem[face[i]];
                face_center /= T(face_size);
                for (int i = 0; i < face_size; i++) {
                    auto a = X_smem[face[i]];
                    auto b = X_smem[face[(i+1) % face_size]];
                    auto c = face_center;
                    auto u = b - a; auto v = c - a;
                    V += dot(a, cross(u, v)) / T(2);
                }
                auto result = sycl::reduce_over_group(cta, V, sycl::plus<T>()) / T(3);
                if (tid == 0) out_volume[bid] = result;
            });
    });
    return ret_val;
}


// Phase 2: see dualize_buffer_size — all return 0 (local_accessor for scratch).
template <typename T, typename K>
size_t eccentricity_buffer_size(const Device&, int, int) { return 0; }
template <typename T, typename K>
size_t inertia_buffer_size(const Device&, int, int) { return 0; }
template <typename T, typename K>
size_t transform_to_principal_axes_buffer_size(const Device&, int, int) { return 0; }
template <typename T, typename K>
size_t surface_area_buffer_size(const Device&, int, int) { return 0; }
template <typename T, typename K>
size_t volume_buffer_size(const Device&, int, int) { return 0; }
template size_t eccentricity_buffer_size<float, uint16_t>(const Device&, int, int);
template size_t inertia_buffer_size<float, uint16_t>(const Device&, int, int);
template size_t transform_to_principal_axes_buffer_size<float, uint16_t>(const Device&, int, int);
template size_t surface_area_buffer_size<float, uint16_t>(const Device&, int, int);
template size_t volume_buffer_size<float, uint16_t>(const Device&, int, int);

template SyclEvent eccentricity<float, uint16_t>(SyclQueue&, std::span<std::array<float,3>>, int, int, batch::BatchStateView, std::span<float>, Workspace);
template SyclEvent inertia<float, uint16_t>(SyclQueue&, std::span<std::array<float,3>>, int, int, batch::BatchStateView, std::span<std::array<float,3>>, Workspace);
template SyclEvent transform_to_principal_axes<float, uint16_t>(SyclQueue&, std::span<std::array<float,3>>, int, int, batch::BatchStateView, Workspace);
template SyclEvent surface_area<float, uint16_t>(SyclQueue&, std::span<std::array<float,3>>, int, int, std::span<std::array<uint16_t,6>>, std::span<uint8_t>, batch::BatchStateView, std::span<float>, Workspace);
template SyclEvent volume<float, uint16_t>(SyclQueue&, std::span<std::array<float,3>>, int, int, std::span<std::array<uint16_t,6>>, std::span<uint8_t>, batch::BatchStateView, std::span<float>, Workspace);