#include <sycl/sycl.hpp>
#include "numeric"
#include <vector>
#include <iostream>
#include <fstream>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/sycl-headers/execution-compat.hh>
#include <fullerenes/kernel-headers/dualize.hh>
#include "queue-impl.hh"
#include "forcefield-includes.hh"
#include "kernel.hh"

template<int MaxDegree, typename K, typename DegT = K>
struct DeviceDualGraph{
    //Check that K is integral
    INT_TYPEDEFS(K);
    
    const std::span<std::array<K,MaxDegree>> dual_neighbours;                          //(Nf x MaxDegree)
    const std::span<DegT> face_degrees;                            //(Nf x 1)
    
    DeviceDualGraph(const std::span<std::array<K,MaxDegree>> dual_neighbours, const std::span<DegT> face_degrees) : dual_neighbours(dual_neighbours), face_degrees(face_degrees) {}

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
                                         std::span<std::array<K,6>> faces_cubic_scratch,
                                         std::span<std::array<K,3>> faces_dual_scratch)
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
    std::span<std::array<K,MaxDegree>> A_dual(
        reinterpret_cast<std::array<K,MaxDegree>*>(src_adj_std.data()),
        src_adj_std.size() / MaxDegree);
    std::span<uint8_t> deg(src_deg_std.data(), src_deg_std.size());
    std::span<std::array<K,3>> A_cubic(
        reinterpret_cast<std::array<K,3>*>(dst_adj_std.data()),
        dst_adj_std.size() / 3);
    std::span<std::array<K,MaxDegree>> faces_cubic = faces_cubic_scratch;
    std::span<std::array<K,3>>         faces_dual  = faces_dual_scratch;
    std::span<StatusFlag> statuses(state.status.data(), state.status.size());

    return launch_per_isomer(Q, N, capacity, [&](handler &h, sycl::nd_range<1> ndr) {
        local_accessor<std::array<K,MaxDegree>, 1> triangle_numbers(Nf, h);
        local_accessor<std::array<K,MaxDegree>, 1> cached_neighbours(Nf, h);
        local_accessor<node_t, 1>                  cached_degrees(Nf, h);
        local_accessor<node2, 1>                   arc_list(N, h);

        h.parallel_for(ndr, [=](nd_item<1> nditem) {
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
                std::span<std::array<K,MaxDegree>>(static_cast<std::array<K,MaxDegree>*>(cached_neighbours.get_pointer()), Nf),
                std::span<node_t>(static_cast<node_t*>(cached_degrees.get_pointer()), Nf));

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
}

template <typename T, typename K>
SyclEvent dualize(SyclQueue& Q,
                  batch::BatchView<Spanify::RSRAdjacencyView<K>> src,
                  batch::BatchView<Spanify::RSRAdjacencyView<K>> dst,
                  batch::BatchStateView                          state,
                  std::span<std::array<K,6>>                          faces_cubic,
                  std::span<std::array<K,3>>                          faces_dual,
                  Workspace                                      /*ws*/)
{
    return dualize_view_batch_impl<T,K>(Q, src, dst, state, faces_cubic, faces_dual);
}

template SyclEvent dualize<float, uint16_t>(SyclQueue&, batch::BatchView<Spanify::RSRAdjacencyView<uint16_t>>, batch::BatchView<Spanify::RSRAdjacencyView<uint16_t>>, batch::BatchStateView, std::span<std::array<uint16_t,6>>, std::span<std::array<uint16_t,3>>, Workspace);
template SyclEvent dualize<float, uint32_t>(SyclQueue&, batch::BatchView<Spanify::RSRAdjacencyView<uint32_t>>, batch::BatchView<Spanify::RSRAdjacencyView<uint32_t>>, batch::BatchStateView, std::span<std::array<uint32_t,6>>, std::span<std::array<uint32_t,3>>, Workspace);
template SyclEvent dualize<double,uint16_t>(SyclQueue&, batch::BatchView<Spanify::RSRAdjacencyView<uint16_t>>, batch::BatchView<Spanify::RSRAdjacencyView<uint16_t>>, batch::BatchStateView, std::span<std::array<uint16_t,6>>, std::span<std::array<uint16_t,3>>, Workspace);
template SyclEvent dualize<double,uint32_t>(SyclQueue&, batch::BatchView<Spanify::RSRAdjacencyView<uint32_t>>, batch::BatchView<Spanify::RSRAdjacencyView<uint32_t>>, batch::BatchStateView, std::span<std::array<uint32_t,6>>, std::span<std::array<uint32_t,3>>, Workspace);

// Phase 2: scratch-size query. Current view-based path uses local_accessor
// for everything so this returns 0; the API is in place so callers can size
// a shared workspace via std::max() over every kernel.
template <typename T, typename K>
size_t dualize_buffer_size(const Device& /*device*/, int /*N*/, int /*capacity*/) {
    return 0;
}
template size_t dualize_buffer_size<float, uint16_t>(const Device&, int, int);
template size_t dualize_buffer_size<float, uint32_t>(const Device&, int, int);
template size_t dualize_buffer_size<double,uint16_t>(const Device&, int, int);
template size_t dualize_buffer_size<double,uint32_t>(const Device&, int, int);
