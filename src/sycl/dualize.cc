#include <CL/sycl.hpp>
#include "numeric"
#include <vector>
#include <iostream>
#include <fstream>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/sycl-headers/dualize-kernel.hh>
#include "primitives.cc"
#include "queue-impl.cc"
#include "forcefield-includes.cc"
#include "kernel.cc"

template<int MaxDegree, typename K>
struct DeviceDualGraph{
    //Check that K is integral
    INT_TYPEDEFS(K);
    
    const K* dual_neighbours;                          //(Nf x MaxDegree)
    const K* face_degrees;                            //(Nf x 1)
    
    DeviceDualGraph(const K* dual_neighbours, const K* face_degrees) : dual_neighbours(dual_neighbours), face_degrees(face_degrees) {}

    K arc_ix(const K u, const K v) const{
        for (uint8_t j = 0; j < face_degrees[u]; j++){
            if (dual_neighbours[u*MaxDegree + j] == v) return j;
        }

        assert(false);
	return -1;		// Make compiler happy
    }
    K arc_ix(const node2& e){ return arc_ix(e[0], e[1]); }

    K next(const K u, const K v) const{
        K j = arc_ix(u,v);
        return dual_neighbours[u*MaxDegree + ((j+1)%face_degrees[u])];
    }
    
    K prev(const K u, const K v) const{
        K j = arc_ix(u,v);
        return dual_neighbours[u*MaxDegree + ((j-1+face_degrees[u])%face_degrees[u])];
    }

    K next_on_face(const K u, const K v) const{
        return prev(v,u);
    }

    K prev_on_face(const K u, const K v) const{
        return next(v,u);
    }

    node2 get_canonical_triangle_arc(const K u, const K v) const{
        //In a triangle u, v, w there are only 3 possible representative arcs, the canonical arc is chosen as the one with the smalles source node.
        node2 min_edge = {u,v};
        K w = next(u,v);
        if (v < u && v < w) min_edge = {v, w};
        if (w < u && w < v) min_edge = {w, u};
        return min_edge;
    }

    node2 get_canonical_arc(K u, K v) const{
        auto start_node = u;
        node2 min_edge = {u,v};
        while (v!= start_node){
            K w = next_on_face(u, v);
            u = v; v = w;
            if(u < min_edge[0]) min_edge = {u,v};
        }
        //assert(next_on_face(u,v) == start_node);
        return min_edge;
    }
};

int roundUp(int numToRound, int multiple) 
{
    assert(multiple);
    return ((numToRound + multiple - 1) / multiple) * multiple;
}


template<typename T, typename K, int MaxDegIn, int MaxDegOut>
SyclEvent dualize_general_impl(  SyclQueue& Q, 
                            Span<K> G_in, 
                            Span<K> Deg_in, 
                            Span<K> G_out, 
                            Span<K> Deg_out,
                            Span<K> cannon_ixs_acc,
                            Span<K> rep_count_acc,
                            Span<K> scan_array_acc,
                            Span<K> triangle_numbers_acc,
                            Span<K> arc_list_acc, 
                            int Nin, 
                            int Nout){
    INT_TYPEDEFS(K);
    sycl::device d = Q->get_device();

    //Find maximum workgroup size
    auto max_workgroup_size = d.get_info<sycl::info::device::max_work_group_size>();

    size_t workgroup_size1 = std::min((int)max_workgroup_size, Nin);
    size_t workgroup_size2 = std::min((int)max_workgroup_size, Nout);
    size_t grid_size1 = roundUp(Nin, workgroup_size1);
    size_t grid_size2 = roundUp(Nout, workgroup_size2);
    
    auto work_distribution = Q->submit([&](handler &h) {
        h.parallel_for(nd_range(range{grid_size1}, range{workgroup_size1}), [=](nd_item<1> nditem) {
            auto thid = nditem.get_global_linear_id();
            DeviceDualGraph<MaxDegIn, K> FD(G_in.data(), Deg_in.data());
            K rep_count = 0;
            if (thid < Nin){
                for (int i = 0; i < FD.face_degrees[thid]; i++){
                    auto canon_arc = FD.get_canonical_arc(thid, FD.dual_neighbours[thid*MaxDegIn + i]);
                    if (canon_arc[0] == thid){
                        cannon_ixs_acc[thid*MaxDegIn + rep_count] = i;
                        ++rep_count;
                    }
                }
                rep_count_acc[thid] = rep_count;
            }
        });
    });
    
    work_distribution.wait();
    primitives::exclusive_scan(Q, rep_count_acc, scan_array_acc, K(0), Plus{});

    auto arc_list_event = Q->submit([&](handler &h) {
        h.parallel_for(nd_range(range{grid_size1}, range{workgroup_size1}), [=](nd_item<1> nditem) {
            auto idx = nditem.get_global_linear_id();
            DeviceDualGraph<MaxDegIn, K> FD(G_in.data(), Deg_in.data());
            if (idx < Nin){
            K rep_count = rep_count_acc[idx];
            K scan_result = scan_array_acc[idx];
            for (int ii = 0; ii < rep_count; ii++){
                K i = cannon_ixs_acc[idx*MaxDegIn + ii];
                K triangle_id = scan_result + ii;
                triangle_numbers_acc[idx*MaxDegIn + i] = triangle_id;
                arc_list_acc[triangle_id*2 + 0] = idx;
                arc_list_acc[triangle_id*2 + 1] = FD.dual_neighbours[idx*MaxDegIn + i];
            }
            }
        });
    });
    
    SyclEventImpl cubic_graph_event = Q->submit([&](handler &h) {
        h.depends_on(arc_list_event);
        h.parallel_for(nd_range(range{grid_size2}, range{workgroup_size2}), [=](nd_item<1> nditem) {
            auto tidx = nditem.get_global_linear_id();
            DeviceDualGraph<MaxDegIn, K> FD(G_in.data(), Deg_in.data());
            if (tidx < Nout){
            K u = arc_list_acc[tidx*2 + 0]; K v = arc_list_acc[tidx*2 + 1];
            auto n_count = 0;
            auto u0 = u;
            node2 edge = FD.get_canonical_arc(v, u); G_out[tidx*MaxDegOut] = triangle_numbers_acc[edge[0]*MaxDegIn + FD.arc_ix(edge)];
            while(v != u0 && n_count < MaxDegOut){
                ++n_count;
                auto w = FD.next_on_face(u,v);
                u = v; v = w;
                edge = FD.get_canonical_arc(v, u); G_out[tidx*MaxDegOut + n_count] = triangle_numbers_acc[edge[0]*MaxDegIn + FD.arc_ix(edge)];
            }
            if (Deg_out.data()) Deg_out[tidx] = n_count+1;
            }
        });
    });
    return SyclEvent(std::move(cubic_graph_event));
}

template <typename T, typename K>
SyclEvent dualize_batch_impl(SyclQueue& Q, FullereneBatchView<T,K> batch){
    INT_TYPEDEFS(K);
    constexpr int     MaxDegree = 6;
    constexpr node_t EMPTY_NODE = std::numeric_limits<node_t>::max();

    auto N = batch.N_;
    auto Nf = batch.Nf_;
    auto capacity =  batch.capacity();
    auto A_dual =    batch.d_.A_dual_;
    auto deg =       batch.d_.deg_;
    auto A_cubic =   batch.d_.A_cubic_;
    auto statuses =  batch.m_.flags_;

    SyclEventImpl cubic_graph_event = Q->submit([&](handler &h) {
        //Create local accessors to shared memory
        local_accessor<node_t, 1>   triangle_numbers(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_neighbours(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_degrees(Nf, h);
        local_accessor<node2, 1>    arc_list(N, h);

        h.parallel_for(sycl::nd_range(sycl::range{N*capacity}, sycl::range{N}), [=](nd_item<1> nditem) {
            auto cta = nditem.get_group();
            node_t f = nditem.get_local_linear_id();    // Face-node index
            auto isomer = nditem.get_group_linear_id(); // Isomer    index
            if(all_set(statuses[isomer], (int)StatusFlag::CUBIC_INITIALIZED)) return;
            sycl::decorated_local_ptr<node_t> cached_neighbours_ptr(cached_neighbours);
            cta.async_work_group_copy(local_ptr<K>(cached_neighbours),    global_ptr<K>(A_dual.begin() + isomer*Nf*MaxDegree), Nf*MaxDegree);
            cta.async_work_group_copy(local_ptr<K>(cached_degrees),       global_ptr<K>(deg.begin()    + isomer*Nf), Nf);

            DeviceDualGraph<MaxDegree, node_t> FD(cached_neighbours.get_pointer(), cached_degrees.get_pointer());

            node_t canon_arcs[MaxDegree]; for(size_t i=0;i<MaxDegree;i++) canon_arcs[i] = EMPTY_NODE; // JA: memset was bug: it writes byte-values, but node_t is 16/32bit.

            node_t rep_count  = 0;
            sycl::group_barrier(cta);

            if (f < Nf){
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    auto canon_arc = FD.get_canonical_triangle_arc(f, FD.dual_neighbours[f*MaxDegree + i]);
                    if (canon_arc[0] == f){
                        canon_arcs[i] = canon_arc[1];
                        rep_count++;
                    }
                }
            }
            sycl::group_barrier(cta);

            node_t scan_result = exclusive_scan_over_group(cta, rep_count, sycl::plus<node_t>{});

            if (f < Nf){
                node_t arc_count = 0;
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    if(canon_arcs[i] != std::numeric_limits<node_t>::max()){
                        triangle_numbers[f*MaxDegree + i] = scan_result + arc_count;
                        ++arc_count;
                    }    
                }
            }
            sycl::group_barrier(cta);

            if (f < Nf){
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    if(canon_arcs[i] != std::numeric_limits<node_t>::max()){
                        node_t u = triangle_numbers[f*MaxDegree + i];
                        arc_list[u] = {f, canon_arcs[i]};
                    }
                }
            }
            sycl::group_barrier(cta);
            auto [u, v] = arc_list[f];
            auto w = FD.next(u,v);
            node2 edge_b = FD.get_canonical_triangle_arc(v, u); A_cubic[isomer*N*3 + f*3 + 0] = triangle_numbers[edge_b[0]*MaxDegree + FD.arc_ix(edge_b)];
            node2 edge_c = FD.get_canonical_triangle_arc(w, v); A_cubic[isomer*N*3 + f*3 + 1] = triangle_numbers[edge_c[0]*MaxDegree + FD.arc_ix(edge_c)];
            node2 edge_d = FD.get_canonical_triangle_arc(u, w); A_cubic[isomer*N*3 + f*3 + 2] = triangle_numbers[edge_d[0]*MaxDegree + FD.arc_ix(edge_d)];

            if (f == 0) statuses[isomer] |= (int)StatusFlag::CUBIC_INITIALIZED;
        });
    });
    return SyclEvent(std::move(cubic_graph_event));
}


/* template <typename T, typename K>
void DualizeFunctor<T,K>::operator()(SyclQueue& Q, FullereneBatchView<T,K> batch, LaunchPolicy policy) {
    
    dispatch_kernel(Q, batch, policy,    
                    dualize_batch_impl<T,K>,
                    [this](SyclQueue& Q, Fullerene<T,K> fullerene){
                        (*this)(Q, fullerene, LaunchPolicy::ASYNC);
                    });
} */

template <typename T, typename K>
template <typename... Data>
SyclEvent DualizeFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T,K> fullerene, Data&&... data){
    if (fullerene.m_.flags_.get() & (int)StatusFlag::CUBIC_INITIALIZED) return SyclEvent(); //Job already done
    if (! (fullerene.m_.flags_.get() & (int)StatusFlag::DUAL_INITIALIZED)) return SyclEvent(); //If the dual graph is not initialized, we cannot proceed.
    auto done_event =  [&](auto&&... data) -> SyclEvent {
        return dualize_general_impl<T,K,6,3>(  Q,
                                        fullerene.d_.A_dual_,
                                        fullerene.d_.deg_,
                                        fullerene.d_.A_cubic_,
                                        Span<K>(),
                                        std::forward<Data>(data)...,
                                        fullerene.Nf_,
                                        fullerene.N_);
    }(std::forward<Data>(data)...);
    fullerene.m_.flags_.get() |= (int)StatusFlag::CUBIC_INITIALIZED;
    return done_event;
}

template <typename T, typename K>
template <typename... Data>
SyclEvent DualizeFunctor<T,K>::compute(SyclQueue& Q, FullereneBatchView<T,K> batch, Data&&... data){
    return dualize_batch_impl<T,K>(Q, batch);
}

template SyclEvent DualizeFunctor<float,uint16_t>::compute(SyclQueue&, FullereneBatchView<float,uint16_t>, Span<uint16_t>&, Span<uint16_t>&, Span<uint16_t>&, Span<uint16_t>&, Span<uint16_t>&);
template SyclEvent DualizeFunctor<float,uint16_t>::compute(SyclQueue&, Fullerene<float, uint16_t>, Span<uint16_t>&, Span<uint16_t>&, Span<uint16_t>&, Span<uint16_t>&, Span<uint16_t>&);
template SyclEvent DualizeFunctor<float,uint32_t>::compute(SyclQueue&, FullereneBatchView<float,uint32_t>, Span<uint32_t>&, Span<uint32_t>&, Span<uint32_t>&, Span<uint32_t>&, Span<uint32_t>&);
template SyclEvent DualizeFunctor<float,uint32_t>::compute(SyclQueue&, Fullerene<float, uint32_t>, Span<uint32_t>&, Span<uint32_t>&, Span<uint32_t>&, Span<uint32_t>&, Span<uint32_t>&);

template struct DualizeFunctor<float,uint16_t>;
template struct DualizeFunctor<float,uint32_t>;
//template struct DualizeFunctor<double,uint16_t>;
//template struct DualizeFunctor<double,uint32_t>;
