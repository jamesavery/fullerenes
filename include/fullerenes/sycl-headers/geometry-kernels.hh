#pragma once
#include <fullerenes/sycl-headers/base-kernel.hh>

template<typename T, typename K>
struct EccentricityFunctor : public KernelFunctor<EccentricityFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_ellipticity);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<T> out_ellipticity, Span<K> indices_);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_ellipticity) const {       return std::make_tuple(std::make_pair(std::ref(indices_), N));}
    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_ellipticity) const { return std::make_tuple();}
};

template<typename T, typename K>
struct InertiaFunctor : public KernelFunctor<InertiaFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<std::array<T,3>> out_inertia);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<std::array<T,3>> out_inertia, Span<K> indices_);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<std::array<T,3>> out_inertia) const {       return std::make_tuple(std::make_pair(std::ref(indices_), N));}
    inline constexpr auto to_tuple_batch(size_t N, Span<std::array<T,3>> out_inertia) const { return std::make_tuple();}
};

template<typename T, typename K>
struct TransformCoordinatesFunctor : public KernelFunctor<TransformCoordinatesFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<K> indices_);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N ) const {       return std::make_tuple(std::make_pair(std::ref(indices_), N));}
    inline constexpr auto to_tuple_batch(size_t N) const { return std::make_tuple();}
};

template<typename T, typename K>
struct SurfaceAreaFunctor : public KernelFunctor<SurfaceAreaFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_surface_area);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<T> out_surface_area, Span<K> indices_);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_surface_area) const {       return std::make_tuple(std::make_pair(std::ref(indices_), N));}
    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_surface_area) const { return std::make_tuple();}
};

template<typename T, typename K>
struct VolumeFunctor : public KernelFunctor<VolumeFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_volume);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<T> out_volume, Span<K> indices_);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_volume) const {       return std::make_tuple(std::make_pair(std::ref(indices_), N));}
    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_volume) const { return std::make_tuple();}
};
