#include "queue-impl.cc"
#include <fullerenes/sycl-headers/execution-compat.hh>
#include <algorithm>
#include <numeric>
#include <fullerenes/kernel-headers/geometry-functors.hh>
#include "forcefield-includes.cc"

template <typename T>
symMat3<T> inertia_matrix(sycl::group<1>& cta, const Span<std::array<T,3>> X){
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
auto inertia_matrix(SyclQueue& Q, const Span<std::array<T,3>> X){
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
auto principal_axes(SyclQueue& Q, const Span<std::array<T,3>> X){
    auto I = inertia_matrix(Q, X);
    auto [V,lambdas] = I.eigensystem();
    return V;
}

template <typename T>
std::array<std::array<T,3>,3> principal_axes(sycl::group<1>& cta, const Span<std::array<T,3>> X){
    auto I = inertia_matrix(cta,X);
    auto [V,lambdas] = I.eigensystem();
    return V;
}

//Returns the best ellipsoid for the given coordinates, lambda0 = a, lambda1 = b, lambda2 = c.
template <typename T>
std::array<T,3> best_ellipsoid (sycl::group<1>& cta,const Span<std::array<T,3>> X){
    auto I = inertia_matrix(X);
    return rsqrt3(d_sort(d_abs(I.eigenvalues()))); 
}

template <typename T, typename K>
SyclEvent EccentricityFunctor<T,K>::compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_ellipticity){
    SyclEventImpl ret_val = Q -> submit([=](sycl::handler& cgh){
        auto N = batch.N_;
        cgh.parallel_for<struct EccentricityFunctor<T,K>>(sycl::nd_range<1>(sycl::range<1>(batch.capacity()*N), sycl::range<1>(N)), [=](sycl::nd_item<1> nditem){
            auto cta = nditem.get_group();
            auto bid = cta.get_group_linear_id();
            if (batch[bid].m_.flags_.get().is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;
            auto X = batch[bid].d_.X_cubic_;
            auto I = inertia_matrix(cta, X);
            auto [V,lambdas] = I.eigensystem();
            auto elipsoid = rsqrt3(d_sort(d_abs(lambdas)));
            out_ellipticity[bid] = elipsoid[0]/elipsoid[2];
        });
        
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent EccentricityFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T, K> fullerene, Span<T> out_ellipticity){
    if (fullerene.m_.flags_.get().is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return SyclEvent();
    auto X = fullerene.d_.X_cubic_;
    auto I = inertia_matrix(Q, X);
    SyclEventImpl ret_val = Q -> single_task ([=](){
        auto [V,lambdas] = I.eigensystem();
        auto elipsoid = rsqrt3(d_sort(d_abs(lambdas)));
        out_ellipticity[0] = elipsoid[0]/elipsoid[2];
    });
    return ret_val;
    //return elipsoid[0]/elipsoid[2];
}

template <typename T, typename K>
SyclEvent InertiaFunctor<T,K>::compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<std::array<T,3>> out_inertia){
    SyclEventImpl ret_val = Q -> submit([=](sycl::handler& cgh){
        auto N = batch.N_;
        cgh.parallel_for<struct InertiaFunctor<T,K>>(sycl::nd_range<1>(sycl::range<1>(batch.capacity()*N), sycl::range<1>(N)), [=](sycl::nd_item<1> nditem){
            auto cta = nditem.get_group();
            auto tid = cta.get_local_linear_id();
            auto bid = cta.get_group_linear_id();
            if (batch[bid].m_.flags_.get().is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;
            auto X = batch[bid].d_.X_cubic_;
            auto I = inertia_matrix(cta, X);
            if (tid == 0) out_inertia[bid] = I.eigenvalues();
        });
        
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent InertiaFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T, K> fullerene, Span<std::array<T,3>> out_inertia){
    if (fullerene.m_.flags_.get().is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return SyclEvent();
    auto X = fullerene.d_.X_cubic_;
    auto I = inertia_matrix(Q, X);
    SyclEventImpl ret_val = Q -> single_task ([=](){
        out_inertia[0] = I.eigenvalues();
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent TransformCoordinatesFunctor<T,K>::compute(SyclQueue& Q, FullereneBatchView<T, K> batch){
    SyclEventImpl ret_val = Q -> submit([=](sycl::handler& cgh){
        auto N = batch.N_;
        cgh.parallel_for<struct TransformCoordinatesFunctor<T,K>>(sycl::nd_range<1>(sycl::range<1>(batch.capacity()*N), sycl::range<1>(N)), [=](sycl::nd_item<1> nditem){
            auto cta = nditem.get_group();
            auto tid = cta.get_local_linear_id();
            auto bid = cta.get_group_linear_id();
            if (batch[bid].m_.flags_.get().is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;
            auto X = batch[bid].d_.X_cubic_;
            auto P = principal_axes(cta, X);
            if (sycl::isnan(P[0][0]) || sycl::isnan(P[0][1]) || sycl::isnan(P[0][2]) || sycl::isnan(P[1][0]) || sycl::isnan(P[1][1]) || sycl::isnan(P[1][2]) || sycl::isnan(P[2][0]) || sycl::isnan(P[2][1]) || sycl::isnan(P[2][2])) {
                return;
            }
            X[tid] = dot(P,X[tid]);
        });
        
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent TransformCoordinatesFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T, K> fullerene){
    auto X = fullerene.d_.X_cubic_;
    auto P = principal_axes(Q, X);
    Q.wait();
    std::transform(FULLERENE_PAR_UNSEQ X.begin(), X.end(), X.begin(), [P](auto x){return dot(P,x);});
    return Q.get_event();
}

template <typename T, typename K>
SyclEvent SurfaceAreaFunctor<T,K>::compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_surface_area){
    SyclEventImpl ret_val = Q -> submit([=](sycl::handler& cgh){
        FLOAT_TYPEDEFS(T);
        auto N = batch.N_;
        auto Nf = batch.Nf_;
        auto X_smem = sycl::local_accessor<std::array<T,3>, 1>(N, cgh);
        auto smem = sycl::local_accessor<K, 1>(N*3 + Nf * 6, cgh);
        cgh.parallel_for<struct SurfaceAreaFunctor<T,K>>(sycl::nd_range<1>(sycl::range<1>(batch.capacity()*Nf), sycl::range<1>(Nf)), [=](sycl::nd_item<1> nditem){
            auto cta = nditem.get_group();
            auto tid = cta.get_local_linear_id();
            auto bid = cta.get_group_linear_id();
            if (batch[bid].m_.flags_.get().is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;
            for (int i = tid; i < N; i += cta.get_local_range(0)) X_smem[i] = batch[bid].d_.X_cubic_[i];
            sycl::group_barrier(cta);
            T A = 0;
            coord3d face_center = {0,0,0};
            auto fullerene = batch[bid];
            auto face_size = fullerene.d_.deg_[tid];
            auto face = fullerene.d_.faces_cubic_[tid];
            for (int i = 0; i < face_size; i++) face_center += X_smem[face[i]];
            face_center /= T(face_size);
            for (int i = 0; i < face_size; i++){
                auto a = X_smem[face[i]];
                auto b = X_smem[face[(i+1)%face_size]];
                auto c = face_center;
                auto u = b - a;
                auto v = c - a;
                auto n = cross(u,v);
                A += norm(n);
            }
            
            auto result = sycl::reduce_over_group(cta, A, sycl::plus<T>()) / T(2);
            if (tid == 0) out_surface_area[bid] = result;
        });
        
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent SurfaceAreaFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T, K> fullerene, Span<T> out_surface_area, Span<K> indices){
    if (fullerene.m_.flags_.get().is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return SyclEvent();
    auto Nf = fullerene.Nf_;
    if (indices.size() != Nf) throw std::runtime_error("Indices size must be equal to the number of faces");

    Q.wait();
    std::iota(indices.begin(), indices.end(), 0);
    auto result = std::transform_reduce(FULLERENE_PAR_UNSEQ indices.begin(), indices.end(), T(0), std::plus<T>{}, [Nf, fullerene](auto tid){
        T A = 0;
        std::array<T,3> face_center = {0,0,0};
        auto face = fullerene.d_.faces_cubic_[tid];
        auto face_size = fullerene.d_.deg_[tid];
        for (int i = 0; i < face_size; i++) face_center += fullerene.d_.X_cubic_[face[i]];
        face_center /= T(face_size);
        for (int i = 0; i < face_size; i++){
            auto a = fullerene.d_.X_cubic_[face[i]];
            auto b = fullerene.d_.X_cubic_[face[(i+1)%face_size]];
            auto c = face_center;
            auto u = b - a;
            auto v = c - a;
            auto n = cross(u,v);
            A += norm(n);
        }
        return A / T(2);
    });
    Q -> single_task([=](){
        out_surface_area[0] = result;
    });
    return Q.get_event();
}

template <typename T, typename K>
SyclEvent VolumeFunctor<T,K>::compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_volume){
    SyclEventImpl ret_val = Q -> submit([=](sycl::handler& cgh){
        FLOAT_TYPEDEFS(T);
        auto N = batch.N_;
        auto Nf = batch.Nf_;
        auto X_smem = sycl::local_accessor<std::array<T,3>, 1>(N, cgh);
        auto smem = sycl::local_accessor<K, 1>(N*3 + Nf * 6, cgh);
        cgh.parallel_for<struct VolumeFunctor<T,K>>(sycl::nd_range<1>(sycl::range<1>(batch.capacity()*Nf), sycl::range<1>(Nf)), [=](sycl::nd_item<1> nditem){
            auto cta = nditem.get_group();
            auto tid = cta.get_local_linear_id();
            auto bid = cta.get_group_linear_id();
            auto fullerene = batch[bid];
            if (fullerene.m_.flags_.get().is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;
            for (int i = tid; i < N; i += cta.get_local_range(0)) X_smem[i] = fullerene.d_.X_cubic_[i];
            //auto node_graph = NodeNeighbours<K>(cta, fullerene.d_.A_cubic_.data(), smem.get_pointer());
            sycl::group_barrier(cta);
            T V = 0;
            coord3d face_center = {0,0,0};
            auto face = fullerene.d_.faces_cubic_[tid];
            auto face_size = fullerene.d_.deg_[tid];
            for (int i = 0; i < face_size; i++) face_center += X_smem[face[i]];
            face_center /= T(face_size);
            for (int i = 0; i < face_size; i++){
                auto a = X_smem[face[i]];
                auto b = X_smem[face[(i+1)%face_size]];
                auto c = face_center;
                auto u = b - a;
                auto v = c - a;
                auto n = cross(u,v);
                V += dot(a, n) / T(2);
            }
            auto result = sycl::reduce_over_group(cta, V, sycl::plus<T>()) / T(3);
            if (tid == 0) out_volume[bid] = result;
        });
        
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent VolumeFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T, K> fullerene, Span<T> out_volume, Span<K> indices){
    if (fullerene.m_.flags_.get().is_not_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return SyclEvent();
    auto Nf = fullerene.Nf_;
    if (indices.size() != Nf) throw std::runtime_error("Indices size must be equal to the number of faces");

    Q.wait();
    std::iota(indices.begin(), indices.end(), 0);
    auto result = std::transform_reduce(FULLERENE_PAR_UNSEQ indices.begin(), indices.end(), T(0), std::plus<T>{}, [Nf, fullerene](auto tid){
        T V = 0;
        std::array<T,3> face_center = {0,0,0};
        auto face = fullerene.d_.faces_cubic_[tid];
        auto face_size = fullerene.d_.deg_[tid];
        for (int i = 0; i < face_size; i++) face_center += fullerene.d_.X_cubic_[face[i]];
        face_center /= T(face_size);
        for (int i = 0; i < face_size; i++){
            auto a = fullerene.d_.X_cubic_[face[i]];
            auto b = fullerene.d_.X_cubic_[face[(i+1)%face_size]];
            auto c = face_center;
            auto u = b - a;
            auto v = c - a;
            auto n = cross(u,v);
            V += dot(a, n) / T(2);
        }
        return V / T(3);
    });
    Q -> single_task([=](){
        out_volume[0] = result;
    });
    return Q.get_event();
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
SyclEvent EccentricityFunctor<T,K>::compute(SyclQueue& Q, Span<std::array<T,3>> xyz, int N, int capacity,
                                             batch::BatchStateView state, Span<T> out_ellipticity) {
    auto statuses = state.status;
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
SyclEvent InertiaFunctor<T,K>::compute(SyclQueue& Q, Span<std::array<T,3>> xyz, int N, int capacity,
                                        batch::BatchStateView state, Span<std::array<T,3>> out_inertia) {
    auto statuses = state.status;
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
SyclEvent TransformCoordinatesFunctor<T,K>::compute(SyclQueue& Q, Span<std::array<T,3>> xyz, int N, int capacity,
                                                     batch::BatchStateView state) {
    auto statuses = state.status;
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
SyclEvent SurfaceAreaFunctor<T,K>::compute(SyclQueue& Q, Span<std::array<T,3>> xyz, int N, int capacity,
                                            Span<std::array<K,6>> faces, Span<uint8_t> deg,
                                            batch::BatchStateView state, Span<T> out_surface_area) {
    FLOAT_TYPEDEFS(T);
    auto statuses = state.status;
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
SyclEvent VolumeFunctor<T,K>::compute(SyclQueue& Q, Span<std::array<T,3>> xyz, int N, int capacity,
                                       Span<std::array<K,6>> faces, Span<uint8_t> deg,
                                       batch::BatchStateView state, Span<T> out_volume) {
    FLOAT_TYPEDEFS(T);
    auto statuses = state.status;
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

template struct EccentricityFunctor<float, uint16_t>;
template struct InertiaFunctor<float, uint16_t>;
template struct TransformCoordinatesFunctor<float, uint16_t>;
template struct SurfaceAreaFunctor<float, uint16_t>;
template struct VolumeFunctor<float, uint16_t>;