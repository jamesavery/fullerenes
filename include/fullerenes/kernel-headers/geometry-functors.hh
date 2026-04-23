#pragma once
#include <fullerenes/kernel-headers/base-functor.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>

template<typename T, typename K>
struct EccentricityFunctor : public KernelFunctor<EccentricityFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_ellipticity);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<T> out_ellipticity);
    // View-based overload (Phase 7). xyz: capacity*N 3D coords.
    SyclEvent compute(SyclQueue& Q, Span<std::array<T,3>> xyz, int N, int capacity,
                      batch::BatchStateView state, Span<T> out_ellipticity);

    inline constexpr auto to_tuple(size_t N, Span<T> out_ellipticity) const {       return std::make_tuple();}
    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_ellipticity) const { return std::make_tuple();}
};

template<typename T, typename K>
struct InertiaFunctor : public KernelFunctor<InertiaFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<std::array<T,3>> out_inertia);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<std::array<T,3>> out_inertia);
    // View-based overload (Phase 7). xyz: capacity*N 3D coords.
    SyclEvent compute(SyclQueue& Q, Span<std::array<T,3>> xyz, int N, int capacity,
                      batch::BatchStateView state, Span<std::array<T,3>> out_inertia);

    inline constexpr auto to_tuple(size_t N, Span<std::array<T,3>> out_inertia) const {       return std::make_tuple();}
    inline constexpr auto to_tuple_batch(size_t N, Span<std::array<T,3>> out_inertia) const { return std::make_tuple();}
};

template<typename T, typename K>
struct TransformCoordinatesFunctor : public KernelFunctor<TransformCoordinatesFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch);
    // View-based overload (Phase 7). xyz: capacity*N 3D coords, updated in-place.
    SyclEvent compute(SyclQueue& Q, Span<std::array<T,3>> xyz, int N, int capacity,
                      batch::BatchStateView state);

    inline constexpr auto to_tuple(size_t N ) const {       return std::make_tuple();}
    inline constexpr auto to_tuple_batch(size_t N) const { return std::make_tuple();}
};

template<typename T, typename K>
struct SurfaceAreaFunctor : public KernelFunctor<SurfaceAreaFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_surface_area);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<T> out_surface_area, Span<K> indices_);
    // View-based overload (Phase 7).
    // xyz: capacity*N 3D coords; faces: capacity*Nf arrays of up to 6 node indices;
    // deg: capacity*Nf face degree (5 or 6).
    SyclEvent compute(SyclQueue& Q, Span<std::array<T,3>> xyz, int N, int capacity,
                      Span<std::array<K,6>> faces, Span<uint8_t> deg,
                      batch::BatchStateView state, Span<T> out_surface_area);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_surface_area) const {
        auto Nf = N/2 + 2;       
        return std::make_tuple(std::make_pair(std::ref(indices_), Nf));}
    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_surface_area) const { return std::make_tuple();}
};

template<typename T, typename K>
struct VolumeFunctor : public KernelFunctor<VolumeFunctor<T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_volume);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<T> out_volume, Span<K> indices_);
    // View-based overload (Phase 7). Same parameters as SurfaceAreaFunctor view overload.
    SyclEvent compute(SyclQueue& Q, Span<std::array<T,3>> xyz, int N, int capacity,
                      Span<std::array<K,6>> faces, Span<uint8_t> deg,
                      batch::BatchStateView state, Span<T> out_volume);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_volume) const {
        auto Nf = N/2 + 2;
        return std::make_tuple(std::make_pair(std::ref(indices_), Nf));}
    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_volume) const { return std::make_tuple();}
};
