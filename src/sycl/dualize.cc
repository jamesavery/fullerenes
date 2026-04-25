#include <sycl/sycl.hpp>
#include "numeric"
#include <vector>
#include <iostream>
#include <fstream>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/sycl-headers/execution-compat.hh>
#include <fullerenes/kernel-headers/dualize-functor.hh>
#include "queue-impl.cc"
#include "forcefield-includes.cc"
#include "kernel.cc"

template<int MaxDegree, typename K, typename DegT = K>
struct DeviceDualGraph{
    //Check that K is integral
    INT_TYPEDEFS(K);
    
    const Span<std::array<K,MaxDegree>> dual_neighbours;                          //(Nf x MaxDegree)
    const Span<DegT> face_degrees;                            //(Nf x 1)
    
    DeviceDualGraph(const Span<std::array<K,MaxDegree>> dual_neighbours, const Span<DegT> face_degrees) : dual_neighbours(dual_neighbours), face_degrees(face_degrees) {}

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
    Q.wait();
    std::exclusive_scan(FULLERENE_PAR_UNSEQ rep_count_acc.begin(), rep_count_acc.end(), scan_array_acc.begin(), K(0), std::plus<K>{});

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


// ---------------------------------------------------------------------------
// View-based batch dualize (Phase 7).
// Same algorithm as dualize_batch_impl, but reads the dual graph from a
// BatchView<TriangulationView>, writes cubic adjacency to a
// BatchView<CubicGraphView>, and carries status through BatchStateView.
// Face-arrays (which have no representation in the views) are passed as
// external scratch/output spans.
// ---------------------------------------------------------------------------
template <typename T, typename K>
static SyclEvent dualize_view_batch_impl(SyclQueue& Q,
                                         batch::BatchView<Spanify::RSRAdjacencyView<K>> src,
                                         batch::BatchView<Spanify::RSRAdjacencyView<K>> dst,
                                         batch::BatchStateView                          state,
                                         Span<std::array<K,6>> faces_cubic_scratch,
                                         Span<std::array<K,3>> faces_dual_scratch)
{
    INT_TYPEDEFS(K);
    constexpr int     MaxDegree  = 6;
    constexpr node_t  EMPTY_NODE = std::numeric_limits<node_t>::max();

    const int Nf       = src.N();
    const int N        = dst.N();
    const int capacity = src.size();
    assert(dst.size()        == capacity);
    assert(int(state.size()) == capacity);

    // Raw underlying spans from each view.
    auto [src_adj_std, src_deg_std, src_twin_std] = src.spans();
    auto [dst_adj_std, dst_deg_std, dst_twin_std] = dst.spans();
    (void)src_twin_std; (void)dst_twin_std; (void)dst_deg_std;

    // Wrap as project Spans (pointer,size) with the typed-row shapes the
    // kernel wants.
    Span<std::array<K,MaxDegree>> A_dual(
        reinterpret_cast<std::array<K,MaxDegree>*>(src_adj_std.data()),
        src_adj_std.size() / MaxDegree);
    Span<uint8_t> deg(src_deg_std.data(), src_deg_std.size());
    Span<std::array<K,3>> A_cubic(
        reinterpret_cast<std::array<K,3>*>(dst_adj_std.data()),
        dst_adj_std.size() / 3);
    Span<std::array<K,MaxDegree>> faces_cubic = faces_cubic_scratch;
    Span<std::array<K,3>>         faces_dual  = faces_dual_scratch;
    Span<StatusFlag> statuses(state.status.data(), state.status.size());

    SyclEventImpl cubic_graph_event = Q->submit([&](handler &h) {
        local_accessor<std::array<K,MaxDegree>, 1> triangle_numbers(Nf, h);
        local_accessor<std::array<K,MaxDegree>, 1> cached_neighbours(Nf, h);
        local_accessor<node_t, 1>                  cached_degrees(Nf, h);
        local_accessor<node2, 1>                   arc_list(N, h);

        h.parallel_for(sycl::nd_range(sycl::range{size_t(N)*size_t(capacity)}, sycl::range{size_t(N)}), [=](nd_item<1> nditem) {
            auto cta     = nditem.get_group();
            node_t f     = nditem.get_local_linear_id();
            auto isomer  = nditem.get_group_linear_id();

            if (all_set(statuses[isomer], (int)StatusEnum::FULLERENEGRAPH_PREPARED)) return;
            if (statuses[isomer].is_not_set(StatusEnum::DUAL_INITIALIZED))           return;

            if (f < Nf){
                cached_neighbours[f] = A_dual[isomer*Nf + f];
                cached_degrees[f]    = node_t(deg[isomer*Nf + f]);
            }

            DeviceDualGraph<MaxDegree, node_t, node_t> FD(
                Span<std::array<K,MaxDegree>>(cached_neighbours.get_pointer(), Nf),
                Span<node_t>(cached_degrees.get_pointer(), Nf));

            node_t canon_arcs[MaxDegree];
            for (size_t i = 0; i < MaxDegree; i++) canon_arcs[i] = EMPTY_NODE;

            node_t rep_count = 0;
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
                    if (canon_arcs[i] != EMPTY_NODE){
                        triangle_numbers[f][i] = scan_result + arc_count;
                        ++arc_count;
                    }
                }
            }
            sycl::group_barrier(cta);

            if (f < Nf){
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    if (canon_arcs[i] != EMPTY_NODE){
                        node_t u = triangle_numbers[f][i];
                        arc_list[u] = {f, canon_arcs[i]};
                    }
                }
            }
            sycl::group_barrier(cta);

            auto [u, v] = arc_list[f];
            auto w      = FD.next(u,v);
            node2 edge_b = FD.get_canonical_triangle_arc(v, u);
            A_cubic[isomer*N + f][0] = triangle_numbers[edge_b[0]][FD.arc_ix(edge_b)];
            node2 edge_c = FD.get_canonical_triangle_arc(w, v);
            A_cubic[isomer*N + f][1] = triangle_numbers[edge_c[0]][FD.arc_ix(edge_c)];
            node2 edge_d = FD.get_canonical_triangle_arc(u, w);
            A_cubic[isomer*N + f][2] = triangle_numbers[edge_d[0]][FD.arc_ix(edge_d)];

            if (f < Nf){
                for (int i = 0; i < FD.face_degrees[f]; i++){
                    auto arc = FD.get_canonical_triangle_arc(f, FD.dual_neighbours[f][i]);
                    faces_cubic[isomer*Nf + f][i] = triangle_numbers[arc[0]][FD.arc_ix(arc)];
                }
            }

            faces_dual[isomer*N + f] = {u, v, w};

            if (f == 0) statuses[isomer].set(StatusEnum::FULLERENEGRAPH_PREPARED);
        });
    });
    return SyclEvent(std::move(cubic_graph_event));
}

template <typename T, typename K>
SyclEvent DualizeFunctor<T,K>::compute(SyclQueue& Q,
                                       batch::BatchView<Spanify::RSRAdjacencyView<K>> src,
                                       batch::BatchView<Spanify::RSRAdjacencyView<K>> dst,
                                       batch::BatchStateView                          state,
                                       Span<std::array<K,6>>                          faces_cubic,
                                       Span<std::array<K,3>>                          faces_dual)
{
    return dualize_view_batch_impl<T,K>(Q, src, dst, state, faces_cubic, faces_dual);
}

template struct DualizeFunctor<float,uint16_t>;
template struct DualizeFunctor<float,uint32_t>;
template struct DualizeFunctor<double,uint16_t>;
template struct DualizeFunctor<double,uint32_t>;
