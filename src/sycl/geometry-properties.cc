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
T EccentricityFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T, K> batch, Span<K> indices_){
    auto N = batch.N_;
    auto X = batch.d_.X_cubic_.template as_span<std::array<T,3>>();
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
    auto [V,lambdas] = I.eigensystem();
    auto elipsoid = rsqrt3(d_sort(d_abs(lambdas)));
    return elipsoid[0]/elipsoid[2];
}

template struct EccentricityFunctor<float, int>;

