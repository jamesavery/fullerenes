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
    
    const Span<std::array<K,MaxDegree>> dual_neighbours;                          //(Nf x MaxDegree)
    const Span<K> face_degrees;                            //(Nf x 1)
    
    DeviceDualGraph(const Span<std::array<K,MaxDegree>> dual_neighbours, const Span<K> face_degrees) : dual_neighbours(dual_neighbours), face_degrees(face_degrees) {}

    K arc_ix(const K u, const K v) const{
        for (uint8_t j = 0; j < face_degrees[u]; j++){
            if (dual_neighbours[u][j] == v) return j;
        }

        assert(false);
	return -1;		// Make compiler happy
    }
    K arc_ix(const node2& e){ return arc_ix(e[0], e[1]); }

    K next(const K u, const K v) const{
        K j = arc_ix(u,v);
        return dual_neighbours[u][((j+1)%face_degrees[u])];
    }
    
    K prev(const K u, const K v) const{
        K j = arc_ix(u,v);
        return dual_neighbours[u][((j-1+face_degrees[u])%face_degrees[u])];
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
    auto faces_cubic = batch.d_.faces_cubic_;
    auto faces_dual = batch.d_.faces_dual_;

    SyclEventImpl cubic_graph_event = Q->submit([&](handler &h) {
        //Create local accessors to shared memory
        local_accessor<std::array<K,MaxDegree>, 1>   triangle_numbers(Nf, h);
        local_accessor<std::array<K,MaxDegree>, 1>   cached_neighbours(Nf, h);
        local_accessor<node_t, 1>   cached_degrees(Nf, h);
        local_accessor<node2, 1>    arc_list(N, h);

        h.parallel_for(sycl::nd_range(sycl::range{N*capacity}, sycl::range{N}), [=](nd_item<1> nditem) {
            auto cta = nditem.get_group();
            node_t f = nditem.get_local_linear_id();    // Face-node index
            auto isomer = nditem.get_group_linear_id(); // Isomer    index
            if(all_set(statuses[isomer], (int)StatusEnum::FULLERENEGRAPH_PREPARED)) return;
            //cta.async_work_group_copy(local_ptr<K>(cached_neighbours),    global_ptr<K>(A_dual.begin() + isomer*Nf*MaxDegree), Nf*MaxDegree);
            //cta.async_work_group_copy(local_ptr<K>(cached_degrees),       global_ptr<K>(deg.begin()    + isomer*Nf), Nf);
            if ( f < Nf){
                cached_neighbours[f] = A_dual[isomer*Nf + f];
                cached_degrees[f] = deg[isomer*Nf + f];
            }

            DeviceDualGraph<MaxDegree, node_t> FD(Span<std::array<K,MaxDegree>>(cached_neighbours.get_pointer(),Nf), Span<K>(cached_degrees.get_pointer(),Nf));

            node_t canon_arcs[MaxDegree]; for(size_t i=0;i<MaxDegree;i++) canon_arcs[i] = EMPTY_NODE; // JA: memset was bug: it writes byte-values, but node_t is 16/32bit.

            node_t rep_count  = 0;
            sycl::group_barrier(cta);


            if (f < Nf){
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    auto canon_arc = FD.get_canonical_triangle_arc(f, FD.dual_neighbours[f][i]);
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
                        triangle_numbers[f][i] = scan_result + arc_count;
                        ++arc_count;
                    }    
                }
            }
            sycl::group_barrier(cta);

            if (f < Nf){
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    if(canon_arcs[i] != std::numeric_limits<node_t>::max()){
                        node_t u = triangle_numbers[f][i];
                        arc_list[u] = {f, canon_arcs[i]};
                    }
                }
            }
            sycl::group_barrier(cta);
            auto [u, v] = arc_list[f];
            auto w = FD.next(u,v);
            node2 edge_b = FD.get_canonical_triangle_arc(v, u); A_cubic[isomer*N + f][0] = triangle_numbers[edge_b[0]][FD.arc_ix(edge_b)];
            node2 edge_c = FD.get_canonical_triangle_arc(w, v); A_cubic[isomer*N + f][1] = triangle_numbers[edge_c[0]][FD.arc_ix(edge_c)];
            node2 edge_d = FD.get_canonical_triangle_arc(u, w); A_cubic[isomer*N + f][2] = triangle_numbers[edge_d[0]][FD.arc_ix(edge_d)];

            //Fill faces_cubic
            if (f < Nf){
                for (int i = 0; i < FD.face_degrees[f]; i++){
                    auto arc = FD.get_canonical_triangle_arc(f, FD.dual_neighbours[f][i]);
                    faces_cubic[isomer*Nf + f][i] = triangle_numbers[arc[0]][FD.arc_ix(arc)];
                }
            }

            //Fill faces_dual
            faces_dual[isomer*N + f] = {u,v,w};

            if (f == 0) statuses[isomer] |= (int)StatusEnum::FULLERENEGRAPH_PREPARED;
        });
    });
    return SyclEvent(std::move(cubic_graph_event));
}

template <typename T, typename K>
SyclEvent prepare_fullerene_graph(SyclQueue& Q, Fullerene<T,K> fullerene, Span<K> cannon_ixs, Span<K> rep_count, Span<K> scan_array, Span<K> triangle_numbers, Span<K> arc_list){
    auto max_workgroup_size = Q->get_device().get_info<sycl::info::device::max_work_group_size>();
    int Nin = fullerene.Nf_;
    int Nout = fullerene.N_;
    size_t workgroup_size1 = std::min((int)max_workgroup_size, Nin);
    size_t workgroup_size2 = std::min((int)max_workgroup_size, Nout);
    size_t grid_size1 = roundUp(Nin, workgroup_size1);
    size_t grid_size2 = roundUp(Nout, workgroup_size2);
    if (cannon_ixs.size() < Nin*6 ||
        rep_count.size() < Nin ||
        scan_array.size() < Nin ||
        triangle_numbers.size() < Nin*6 ||
        arc_list.size() < Nout*2){
        throw std::runtime_error("Insufficient memory for dualization");
    }
    //primitives::transform_exclusive_scan(Q, 

    auto work_distribution = Q -> parallel_for(nd_range(range{grid_size1}, range{workgroup_size1}), [=](nd_item<1> nditem) {
        auto thid = nditem.get_global_linear_id();
        DeviceDualGraph<6, K> FD(fullerene.d_.A_dual_, fullerene.d_.deg_);
        K local_rep_count = 0;
        if (thid < Nin){
            for (int i = 0; i < FD.face_degrees[thid]; i++){
                auto canon_arc = FD.get_canonical_arc(thid, FD.dual_neighbours[thid][i]);
                if (canon_arc[0] == thid){
                    cannon_ixs[thid*6 + local_rep_count] = i;
                    ++local_rep_count;
                }
            }
            rep_count[thid] = local_rep_count;
        }
    });

    work_distribution.wait();
    primitives::exclusive_scan(Q, rep_count, scan_array, K(0), Plus{});

    auto arc_list_event = Q -> parallel_for(nd_range(range{grid_size1}, range{workgroup_size1}), [=](nd_item<1> nditem) {
        auto idx = nditem.get_global_linear_id();
        DeviceDualGraph<6, K> FD(fullerene.d_.A_dual_, fullerene.d_.deg_);
        if (idx < Nin){
        K rep_count_local = rep_count[idx];
        K scan_result = scan_array[idx];
        for (int ii = 0; ii < rep_count_local; ii++){
            auto i = cannon_ixs[idx*6 + ii];
            auto triangle_id = scan_result + ii;
            triangle_numbers[idx*6 + i] = triangle_id;
            arc_list.template as_span<std::array<K,2>>()[triangle_id] = {(K)idx, FD.dual_neighbours[idx][i]};
        }
        }
    });
    
    SyclEventImpl fill_graph_event = Q -> parallel_for(nd_range(range{grid_size2}, range{workgroup_size2}), arc_list_event, [=](nd_item<1> nditem) {
        auto tidx = nditem.get_global_linear_id();
        DeviceDualGraph<6, K> FD(fullerene.d_.A_dual_, fullerene.d_.deg_);
        if (tidx < Nout){
            auto [u, v] = arc_list.template as_span<std::array<K,2>>()[tidx];
            auto w = FD.next(u,v);
            auto edge_b = FD.get_canonical_triangle_arc(v, u); fullerene.d_.A_cubic_[tidx][0] = triangle_numbers[edge_b[0]*6 + FD.arc_ix(edge_b)];
            auto edge_c = FD.get_canonical_triangle_arc(w, v); fullerene.d_.A_cubic_[tidx][1] = triangle_numbers[edge_c[0]*6 + FD.arc_ix(edge_c)];
            auto edge_d = FD.get_canonical_triangle_arc(u, w); fullerene.d_.A_cubic_[tidx][2] = triangle_numbers[edge_d[0]*6 + FD.arc_ix(edge_d)];

            fullerene.d_.faces_dual_.template as_span<std::array<K,3>>() [tidx] = {u,v,w};
        }

        if (tidx < Nin){
            for (int i = 0; i < FD.face_degrees[tidx]; i++){
                auto arc = FD.get_canonical_triangle_arc(tidx, FD.dual_neighbours[tidx][i]);
                fullerene.d_.faces_cubic_[tidx][i] = triangle_numbers[arc[0]*6 + FD.arc_ix(arc)];
            }
        }
        
    });
    return SyclEvent(std::move(fill_graph_event));

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
SyclEvent DualizeFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T,K> fullerene, Span<K> cannon_ixs, Span<K> rep_count, Span<K> scan_array, Span<K> triangle_numbers, Span<K> arc_list){
    if (fullerene.m_.flags_.get() & (int)StatusEnum::FULLERENEGRAPH_PREPARED) return SyclEvent(); //Job already done
    if (! (fullerene.m_.flags_.get() & (int)StatusEnum::DUAL_INITIALIZED)) return SyclEvent(); //If the dual graph is not initialized, we cannot proceed.
    auto done_event = prepare_fullerene_graph(  Q,
                                                fullerene,
                                                cannon_ixs,
                                                rep_count,
                                                scan_array,
                                                triangle_numbers,
                                                arc_list);
    fullerene.m_.flags_.get() |= (int)StatusEnum::FULLERENEGRAPH_PREPARED;
    return done_event;
}

template <typename T, typename K>
SyclEvent DualizeFunctor<T,K>::compute(SyclQueue& Q, FullereneBatchView<T,K> batch){
    return dualize_batch_impl<T,K>(Q, batch);

}

template struct DualizeFunctor<float,uint16_t>;
template struct DualizeFunctor<float,uint32_t>;
//template struct DualizeFunctor<double,uint16_t>;
//template struct DualizeFunctor<double,uint32_t>;
