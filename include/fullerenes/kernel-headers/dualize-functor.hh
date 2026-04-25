#pragma once
#include <fullerenes/kernel-headers/base-functor.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

template <typename T, typename K>
struct DualizeFunctor : public KernelFunctor<DualizeFunctor<T,K>> {
    // View-based batch overload (Phase 7).
    // Uses Spanify::RSRAdjacencyView<K> directly so the index width matches
    // the functor's K. (The named graph views in graphview.hh are fixed to
    // node_t=int32 and will be templated on K in a later commit; until then
    // kernel entry-points address the sized adjacency view directly.)
    //
    //   src: dual-graph batch (Nv=Nf, dmax<=6). Reads neighbours/deg.
    //   dst: cubic-graph batch (Nv=N,  dmax==3). Writes neighbours.
    //   state: per-entry status; honours DUAL_INITIALIZED and sets
    //          FULLERENEGRAPH_PREPARED on success.
    //   faces_cubic: capacity*Nf output of triangle indices per dual vertex.
    //   faces_dual : capacity*N  output of dual-vertex triples per triangle.
    SyclEvent compute(SyclQueue& Q,
                      batch::BatchView<Spanify::RSRAdjacencyView<K>> src,
                      batch::BatchView<Spanify::RSRAdjacencyView<K>> dst,
                      batch::BatchStateView                          state,
                      Span<std::array<K,6>>                          faces_cubic,
                      Span<std::array<K,3>>                          faces_dual);

    mutable FunctorArrays<K> cannon_ixs_;
    mutable FunctorArrays<K> rep_count_;
    mutable FunctorArrays<K> scan_array_;
    mutable FunctorArrays<K> triangle_numbers_;
    mutable FunctorArrays<K> arc_list_;
    
    inline constexpr auto to_tuple(size_t N) const {
        size_t Nin = (N/2) + 2;
        size_t Nout = N;
        size_t MaxDegree = 6;
        return  std::make_tuple(
                std::make_pair(std::ref(cannon_ixs_),       Nin * MaxDegree), 
                std::make_pair(std::ref(rep_count_),        Nin), 
                std::make_pair(std::ref(scan_array_),       Nin), 
                std::make_pair(std::ref(triangle_numbers_), Nin*MaxDegree),
                std::make_pair(std::ref(arc_list_),         Nout*2 ));
    }

    inline constexpr auto to_tuple_batch(size_t N) const {return std::make_tuple();}
};