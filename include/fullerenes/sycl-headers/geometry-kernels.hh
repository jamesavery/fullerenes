#pragma once
#include <fullerenes/sycl-headers/base-kernel.hh>

template<typename T, typename K>
struct EccentricityFunctor : public KernelFunctor<EccentricityFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_ellipticity);
    T compute(SyclQueue& Q, Fullerene<T, K> batch, Span<K> indices_);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_ellipticity) const {       return std::make_tuple(std::make_pair(std::ref(indices_), N));}
    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_ellipticity) const { return std::make_tuple();}
};

template<typename T, typename K>
struct InertiaFunctor : public KernelFunctor<InertiaFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_inertia);
    T compute(SyclQueue& Q, Fullerene<T, K> batch, Span<K> indices_);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_inertia) const {       return std::make_tuple(std::make_pair(std::ref(indices_), N));}
    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_inertia) const { return std::make_tuple();}
};

template<typename T, typename K>
struct TransformCoordinatesFunctor : public KernelFunctor<TransformCoordinatesFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch);
    T compute(SyclQueue& Q, Fullerene<T, K> batch, Span<K> indices_);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N ) const {       return std::make_tuple(std::make_pair(std::ref(indices_), N));}
    inline constexpr auto to_tuple_batch(size_t N) const { return std::make_tuple();}
};

template<typename T, typename K>
struct SurfaceAreaFunctor : public KernelFunctor<SurfaceAreaFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_surface_area);
    T compute(SyclQueue& Q, Fullerene<T, K> batch, Span<K> indices_);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_surface_area) const {       return std::make_tuple(std::make_pair(std::ref(indices_), N));}
    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_surface_area) const { return std::make_tuple();}
};

template<typename T, typename K>
struct VolumeFunctor : public KernelFunctor<VolumeFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_volume);
    T compute(SyclQueue& Q, Fullerene<T, K> batch, Span<K> indices_);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_volume) const {       return std::make_tuple(std::make_pair(std::ref(indices_), N));}
    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_volume) const { return std::make_tuple();}
};
