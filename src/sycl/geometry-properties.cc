#include <fullerenes/sycl-headers/geometry-kernels.hh>
#include "forcefield-includes.cc"
#include "queue-impl.cc"
#include "primitives.cc"

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
    T diag = primitives::transform_reduce(Q, X, T(0), sycl::plus<T>(), [](const auto& x) -> T {return dot(x,x);});
    I.a = diag;
    I.d = diag;
    I.f = diag;
    I.a -= primitives::transform_reduce(Q, X, T(0), sycl::plus<T>(), [](const auto& x){return x[0]*x[0];});
    I.b -= primitives::transform_reduce(Q, X, T(0), sycl::plus<T>(), [](const auto& x){return x[0]*x[1];});
    I.c -= primitives::transform_reduce(Q, X, T(0), sycl::plus<T>(), [](const auto& x){return x[0]*x[2];});
    I.d -= primitives::transform_reduce(Q, X, T(0), sycl::plus<T>(), [](const auto& x){return x[1]*x[1];});
    I.e -= primitives::transform_reduce(Q, X, T(0), sycl::plus<T>(), [](const auto& x){return x[1]*x[2];});
    I.f -= primitives::transform_reduce(Q, X, T(0), sycl::plus<T>(), [](const auto& x){return x[2]*x[2];});
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
        cgh.parallel_for<struct EccentricityFunctor<T,K>>(sycl::nd_range<1>(sycl::range<1>(batch.capacity()*N), sycl::range<1>(N)), [=](sycl::group<1> cta){
            auto tid = cta.get_local_linear_id();
            auto bid = cta.get_group_linear_id();
            auto X = batch[bid].d_.X_cubic_.template as_span<std::array<T,3>>();
            auto I = inertia_matrix(cta, X);
            auto [V,lambdas] = I.eigensystem();
            auto elipsoid = rsqrt3(d_sort(d_abs(lambdas)));
            out_ellipticity[bid] = elipsoid[0]/elipsoid[2];
        });
        
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent EccentricityFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T, K> fullerene, Span<T> out_ellipticity, Span<K> indices){
    auto N = fullerene.N_;
    auto X = fullerene.d_.X_cubic_.template as_span<std::array<T,3>>();
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
        cgh.parallel_for<struct InertiaFunctor<T,K>>(sycl::nd_range<1>(sycl::range<1>(batch.capacity()*N), sycl::range<1>(N)), [=](sycl::group<1> cta){
            auto tid = cta.get_local_linear_id();
            auto bid = cta.get_group_linear_id();
            auto X = batch[bid].d_.X_cubic_.template as_span<std::array<T,3>>();
            auto I = inertia_matrix(cta, X);
            if (tid == 0) out_inertia[bid] = I.eigenvalues();
        });
        
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent InertiaFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T, K> fullerene, Span<std::array<T,3>> out_inertia, Span<K> indices){
    auto N = fullerene.N_;
    auto X = fullerene.d_.X_cubic_.template as_span<std::array<T,3>>();
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
        cgh.parallel_for<struct TransformCoordinatesFunctor<T,K>>(sycl::nd_range<1>(sycl::range<1>(batch.capacity()*N), sycl::range<1>(N)), [=](sycl::group<1> cta){
            auto tid = cta.get_local_linear_id();
            auto bid = cta.get_group_linear_id();
            auto X = batch[bid].d_.X_cubic_.template as_span<std::array<T,3>>();
            auto P = principal_axes(cta, X);
            if (isnan(P[0][0]) || isnan(P[0][1]) || isnan(P[0][2]) || isnan(P[1][0]) || isnan(P[1][1]) || isnan(P[1][2]) || isnan(P[2][0]) || isnan(P[2][1]) || isnan(P[2][2])) {
                return;
            }
            X[tid] = dot(P,X[tid]);
        });
        
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent TransformCoordinatesFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T, K> fullerene, Span<K> indices){
    auto N = fullerene.N_;
    auto X = fullerene.d_.X_cubic_.template as_span<std::array<T,3>>();
    auto P = principal_axes(Q, X);
    primitives::transform(Q, X, X, [P](auto x){return dot(P,x);});
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
        cgh.parallel_for<struct SurfaceAreaFunctor<T,K>>(sycl::nd_range<1>(sycl::range<1>(batch.capacity()*N), sycl::range<1>(N)), [=](sycl::group<1> cta){
            auto tid = cta.get_local_linear_id();
            auto bid = cta.get_group_linear_id();
            X_smem[tid] = batch[bid].d_.X_cubic_.template as_span<std::array<T,3>>()[tid];
            auto node_graph = NodeNeighbours<K>(cta, batch[bid].d_.A_cubic_.data(), smem.get_pointer());
            sycl::group_barrier(cta);
            T A = 0;
            if (tid < Nf){
                coord3d face_center = {0,0,0};
                for (int i = 0; i < node_graph.face_size; i++) face_center += X_smem[node_graph.face_nodes[i]];
                face_center /= T(node_graph.face_size);
                for (int i = 0; i < node_graph.face_size; i++){
                    auto a = X_smem[node_graph.face_nodes[i]];
                    auto b = X_smem[node_graph.face_nodes[(i+1)%node_graph.face_size]];
                    auto c = face_center;
                    auto u = b - a;
                    auto v = c - a;
                    auto n = cross(u,v);
                    A += norm(n);
                }
            }
            auto result = sycl::reduce_over_group(cta, A, sycl::plus<T>()) / T(2);
            if (tid == 0) out_surface_area[bid] = result;
        });
        
    });
    return ret_val;
}

template <typename T, typename K>
SyclEvent SurfaceAreaFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T, K> fullerene, Span<T> out_surface_area, Span<K> indices){
    auto N = fullerene.N_;
    auto Nf = fullerene.Nf_;
    //TODO: Implement this

    /* auto X_smem = sycl::local_accessor<std::array<T,3>, 1>(N);
    auto smem = sycl::local_accessor<K, 1>(N*3 + Nf * 6);
    auto X = fullerene.d_.X_cubic_.template as_span<std::array<T,3>>();
    auto node_graph = NodeNeighbours<T,K>(fullerene.d_.A_cubic_.data());
    T A = 0;
    for (int i = 0; i < Nf; i++){
        std::array<T,3> face_center = {0,0,0};
        for (int i = 0; i < node_graph.face_size; i++) face_center += X[node_graph.face_nodes[i]];
        face_center /= node_graph.face_size;
        for (int i = 0; i < node_graph.face_size; i++){
            auto a = X[node_graph.face_nodes[i]];
            auto b = X[node_graph.face_nodes[(i+1)%node_graph.face_size]];
            auto c = face_center;
            auto u = b - a;
            auto v = c - a;
            auto n = cross(u,v);
            A += norm(n);
        }
    }
    out_surface_area[0] = A / T(2);
    return Q.get_event(); */
}

template <typename T, typename K>
SyclEvent VolumeFunctor<T,K>::compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_volume){
    SyclEventImpl ret_val = Q -> submit([=](sycl::handler& cgh){
        FLOAT_TYPEDEFS(T);
        auto N = batch.N_;
        auto Nf = batch.Nf_;
        auto X_smem = sycl::local_accessor<std::array<T,3>, 1>(N, cgh);
        auto smem = sycl::local_accessor<K, 1>(N*3 + Nf * 6, cgh);
        cgh.parallel_for<struct VolumeFunctor<T,K>>(sycl::nd_range<1>(sycl::range<1>(batch.capacity()*N), sycl::range<1>(N)), [=](sycl::group<1> cta){
            auto tid = cta.get_local_linear_id();
            auto bid = cta.get_group_linear_id();
            X_smem[tid] = batch[bid].d_.X_cubic_.template as_span<std::array<T,3>>()[tid];
            auto node_graph = NodeNeighbours<K>(cta, batch[bid].d_.A_cubic_.data(), smem.get_pointer());
            sycl::group_barrier(cta);
            T V = 0;
            if (tid < Nf){
                coord3d face_center = {0,0,0};
                for (int i = 0; i < node_graph.face_size; i++) face_center += X_smem[node_graph.face_nodes[i]];
                face_center /= T(node_graph.face_size);
                for (int i = 0; i < node_graph.face_size; i++){
                    auto a = X_smem[node_graph.face_nodes[i]];
                    auto b = X_smem[node_graph.face_nodes[(i+1)%node_graph.face_size]];
                    auto c = face_center;
                    auto u = b - a;
                    auto v = c - a;
                    auto n = cross(u,v);
                    V += dot(a, n) / T(2);
                }
            }
            auto result = sycl::reduce_over_group(cta, V, sycl::plus<T>()) / T(3);
            if (tid == 0) out_volume[bid] = result;
        });
        
    });
    return ret_val;
}

template struct EccentricityFunctor<float, int>;
template struct InertiaFunctor<float, int>;
template struct TransformCoordinatesFunctor<float, int>;
template struct SurfaceAreaFunctor<float, int>;
template struct VolumeFunctor<float, int>;