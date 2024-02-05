#include <CL/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/sycl-isomer-batch.hh>
#include "forcefield-includes.cc"

//Template specialisation for dualize


#define UINT_TYPE_MAX std::numeric_limits<UINT_TYPE>::max()

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

    /**
     * @brief returns the next node in the clockwise order around u
     * @param v the current node around u
     * @param u the node around which the search is performed
     * @return the next node in the clockwise order around u
     */
    K next(const K u, const K v) const{
        K j = arc_ix(u,v);
        return dual_neighbours[u*MaxDegree + ((j+1)%face_degrees[u])];
    }
    
    /**
     * @brief returns the prev node in the clockwise order around u
     * @param v the current node around u
     * @param u the node around which the search is performed
     * @return the previous node in the clockwise order around u
     */
    K prev(const K u, const K v) const{
        K j = arc_ix(u,v);
        return dual_neighbours[u*MaxDegree + ((j-1+face_degrees[u])%face_degrees[u])];
    }

    /**
     * @brief Find the node that comes next on the face. given by the edge (u,v)
     * @param u Source of the edge.
     * @param v Destination node.
     * @return The node that comes next on the face.
     */
    K next_on_face(const K u, const K v) const{
        return prev(v,u);
    }

    /**
     * @brief Find the node that comes next on the face. given by the edge (u,v)
     * @param u Source of the edge.
     * @param v Destination node.
     * @return The node that comes next on the face.
     */
    K prev_on_face(const K u, const K v) const{
        return next(v,u);
    }

    /**
     * @brief Finds the canonical triangle arc of the triangle (u,v,w)
     * 
     * @param u source node
     * @param v target node
     * @return canonical triangle arc 
     */
    node2 get_canonical_triangle_arc(const K u, const K v) const{
        //In a triangle u, v, w there are only 3 possible representative arcs, the canonical arc is chosen as the one with the smalles source node.
        node2 min_edge = {u,v};
        K w = next(u,v);
        if (v < u && v < w) min_edge = {v, w};
        if (w < u && w < v) min_edge = {w, u};
        return min_edge;
    }
};
template <typename T, typename K>
void dualize(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy){
    INT_TYPEDEFS(K);
    constexpr int     MaxDegree = 6;
    constexpr node_t EMPTY_NODE = std::numeric_limits<node_t>::max();
    
    if(policy == LaunchPolicy::SYNC) Q.wait();
    Q.submit([&](handler &h) {
        auto N = batch.N();
        auto Nf = batch.Nf();
        auto capacity = batch.capacity();
	//        std::cout << "Entered dualize" << std::endl;

	// auto d = Q.get_device();
	// std::cout << "Running on " << d.get_info<info::device::name>() << "\n";
	
        //Create local accessors to shared memory
        local_accessor<node_t, 1>   triangle_numbers(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_neighbours(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_degrees(Nf, h);
        local_accessor<node2, 1>    arc_list(N, h);

        //Create device accessors
        accessor     cubic_neighbours_dev(batch.cubic_neighbours, h, write_only);
        accessor     face_degrees_dev    (batch.face_degrees, h, read_only);
        accessor     dual_neighbours_dev (batch.dual_neighbours, h, read_only);
        accessor     statuses_dev        (batch.statuses, h, read_only);
        /* 
        std::cout << N * capacity << std::endl; */
        h.parallel_for(sycl::nd_range(sycl::range{N*capacity}, sycl::range{N}), [=](nd_item<1> nditem) {
            auto cta = nditem.get_group();
            node_t f = nditem.get_local_linear_id();    // Face-node index
            auto isomer = nditem.get_group_linear_id(); // Isomer    index
            if(statuses_dev[isomer] != IsomerStatus::NOT_CONVERGED) return;
	    
            cta.async_work_group_copy(cached_neighbours.get_pointer(), dual_neighbours_dev.get_pointer() + isomer*Nf*MaxDegree, Nf*MaxDegree);
            cta.async_work_group_copy(cached_degrees.get_pointer(),    face_degrees_dev.get_pointer()    + isomer*Nf, Nf);

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

            node_t scan_result = exclusive_scan_over_group(cta, rep_count, plus<node_t>{});

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
//
            auto [u, v] = arc_list[f];
            auto w = FD.next(u,v);
//
            node2 edge_b = FD.get_canonical_triangle_arc(v, u); cubic_neighbours_dev[isomer*N*3 + f*3 + 0] = triangle_numbers[edge_b[0]*MaxDegree + FD.arc_ix(edge_b)];
            node2 edge_c = FD.get_canonical_triangle_arc(w, v); cubic_neighbours_dev[isomer*N*3 + f*3 + 1] = triangle_numbers[edge_c[0]*MaxDegree + FD.arc_ix(edge_c)];
            node2 edge_d = FD.get_canonical_triangle_arc(u, w); cubic_neighbours_dev[isomer*N*3 + f*3 + 2] = triangle_numbers[edge_d[0]*MaxDegree + FD.arc_ix(edge_d)];
        });
    });

    if(policy == LaunchPolicy::SYNC) Q.wait();
}

int roundUp(int numToRound, int multiple) 
{
    assert(multiple);
    return ((numToRound + multiple - 1) / multiple) * multiple;
}

template <typename T, typename K>
void dualize_V1(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy){
    INT_TYPEDEFS(K);
    constexpr int     MaxDegree = 6;
    constexpr node_t EMPTY_NODE = std::numeric_limits<node_t>::max();
    
    if(policy == LaunchPolicy::SYNC) Q.wait();
    auto subgroup_size = Q.get_device().get_info<info::device::sub_group_sizes>()[0];
    size_t lcm = roundUp(batch.Nf(), subgroup_size);
    Q.submit([&](handler &h) {
        auto N = batch.N();
        auto Nf = batch.Nf();
        auto capacity = batch.capacity();
	//        std::cout << "Entered dualize" << std::endl;

	// auto d = Q.get_device();
	// std::cout << "Running on " << d.get_info<info::device::name>() << "\n";
	
        //Create local accessors to shared memory
        local_accessor<node_t, 1>   triangle_numbers(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_neighbours(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_degrees(Nf, h);
        local_accessor<node2, 1>    arc_list(N, h);

        //Create device accessors
        accessor     cubic_neighbours_dev(batch.cubic_neighbours, h, write_only);
        accessor     face_degrees_dev    (batch.face_degrees, h, read_only);
        accessor     dual_neighbours_dev (batch.dual_neighbours, h, read_only);
        accessor     statuses_dev        (batch.statuses, h, read_only);
        /* 
        std::cout << N * capacity << std::endl; */
        h.parallel_for(sycl::nd_range(sycl::range{lcm*capacity}, sycl::range{lcm}), [=](nd_item<1> nditem) {
            auto cta = nditem.get_group();
            node_t f = nditem.get_local_linear_id();    // Face-node index
            auto isomer = nditem.get_group_linear_id(); // Isomer    index
            if(statuses_dev[isomer] != IsomerStatus::NOT_CONVERGED) return;
	    
            cta.async_work_group_copy(cached_neighbours.get_pointer(), dual_neighbours_dev.get_pointer() + isomer*Nf*MaxDegree, Nf*MaxDegree);
            cta.async_work_group_copy(cached_degrees.get_pointer(),    face_degrees_dev.get_pointer()    + isomer*Nf, Nf);

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

            node_t scan_result = exclusive_scan_over_group(cta, rep_count, plus<node_t>{});

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
//          
            for(auto tix = f; tix < N; tix += lcm){
                auto [u, v] = arc_list[tix];
                auto w = FD.next(u,v);
//
                node2 edge_b = FD.get_canonical_triangle_arc(v, u); cubic_neighbours_dev[isomer*N*3 + tix*3 + 0] = triangle_numbers[edge_b[0]*MaxDegree + FD.arc_ix(edge_b)];
                node2 edge_c = FD.get_canonical_triangle_arc(w, v); cubic_neighbours_dev[isomer*N*3 + tix*3 + 1] = triangle_numbers[edge_c[0]*MaxDegree + FD.arc_ix(edge_c)];
                node2 edge_d = FD.get_canonical_triangle_arc(u, w); cubic_neighbours_dev[isomer*N*3 + tix*3 + 2] = triangle_numbers[edge_d[0]*MaxDegree + FD.arc_ix(edge_d)];
            }
        });
    });

    if(policy == LaunchPolicy::SYNC) Q.wait();
}

template void dualize<float,uint16_t>(sycl::queue&Q, IsomerBatch<float,uint16_t>& batch, const LaunchPolicy policy);
template void dualize_V1<float,uint16_t>(sycl::queue&Q, IsomerBatch<float,uint16_t>& batch, const LaunchPolicy policy);
template void dualize<double,uint16_t>(sycl::queue&Q, IsomerBatch<double,uint16_t>& batch, const LaunchPolicy policy);
template void dualize_V1<double,uint16_t>(sycl::queue&Q, IsomerBatch<double,uint16_t>& batch, const LaunchPolicy policy);
