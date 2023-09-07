#include <CL/sycl.hpp>
#include "forcefield_includes.cpp"

#define UINT_TYPE_MAX std::numeric_limits<UINT_TYPE>::max()

template<int MaxDegree, typename K>
struct DeviceDualGraph{
    //Check that K is integral
    static_assert(std::is_integral<K>::value, "K must be integral");

    const K* dual_neighbours;                          //(Nf x MaxDegree)
    const uint8_t* face_degrees;                            //(Nf x 1)
    
    DeviceDualGraph(const K* dual_neighbours, const uint8_t* face_degrees) : dual_neighbours(dual_neighbours), face_degrees(face_degrees) {}

    K dedge_ix(const K u, const K v) const{
        for (uint8_t j = 0; j < face_degrees[u]; j++){
            if (dual_neighbours[u*MaxDegree + j] == v) return j;
        }

        assert(false);
	    return 0;		// Make compiler happy
    }

    /**
     * @brief returns the next node in the clockwise order around u
     * @param v the current node around u
     * @param u the node around which the search is performed
     * @return the next node in the clockwise order around u
     */
    K next(const K u, const K v) const{
        K j = dedge_ix(u,v);
        return dual_neighbours[u*MaxDegree + ((j+1)%face_degrees[u])];
    }
    
    /**
     * @brief returns the prev node in the clockwise order around u
     * @param v the current node around u
     * @param u the node around which the search is performed
     * @return the previous node in the clockwise order around u
     */
    K prev(const K u, const K v) const{
        K j = dedge_ix(u,v);
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
     * @brief Finds the cannonical triangle arc of the triangle (u,v,w)
     * 
     * @param u source node
     * @param v target node
     * @return cannonical triangle arc 
     */
    std::array<K,2> get_cannonical_triangle_arc(const K u, const K v) const{
        //In a triangle u, v, w there are only 3 possible representative arcs, the cannonical arc is chosen as the one with the smalles source node.
        std::array<K,2> min_edge = {u,v};
        K w = next(u,v);
        if (v < u && v < w) min_edge = {v, w};
        if (w < u && w < v) min_edge = {w, u};
        return min_edge;
    }
};

template <typename T, typename K>
void dualise(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy pol){
    if(pol == LaunchPolicy::SYNC) Q.wait();
    Q.submit([&](handler &h) {
        // Create a command group to issue GPU work.
        local_accessor<UINT_TYPE, 1>    triangle_numbers(Nf*MaxDegree, h);
        local_accessor<UINT_TYPE, 1>    cached_neighbours(Nf*MaxDegree, h);
        local_accessor<uint8_t, 1>      cached_degrees(Nf, h);
        local_accessor<std::array<UINT_TYPE,2>, 1> arc_list(N, h);

        h.parallel_for<class dualise>(nd_range(range{N*batch_size}, range{N}), [=](nd_item<1> nditem) {
            auto cta = nditem.get_group();
            auto result = reduce_over_group(cta, 1, plus<int>{}); // Should be size of work-group (N)
            auto thid = nditem.get_local_linear_id();
            auto bid = nditem.get_group_linear_id();
            
            cta.async_work_group_copy(cached_neighbours.get_pointer(), dual_neighbours_dev_ptr + bid*Nf*MaxDegree, Nf*MaxDegree);
            cta.async_work_group_copy(cached_degrees.get_pointer(), face_degrees_dev_ptr + bid*Nf, Nf);
            DeviceDualGraph<MaxDegree, UINT_TYPE> FD(cached_neighbours.get_pointer(), cached_degrees.get_pointer());
            UINT_TYPE cannon_arcs[MaxDegree]; memset(cannon_arcs, UINT_TYPE_MAX, MaxDegree*sizeof(UINT_TYPE));
            UINT_TYPE rep_count  = 0;
            nditem.barrier(access::fence_space::local_space);
            if (thid < Nf){
                for (UINT_TYPE i = 0; i < FD.face_degrees[thid]; i++){
                    auto cannon_arc = FD.get_cannonical_triangle_arc(thid, FD.dual_neighbours[thid*MaxDegree + i]);
                    if (cannon_arc[0] == thid){
                        cannon_arcs[i] = cannon_arc[1];
                        rep_count++;
                    }
                }
            }
            nditem.barrier(access::fence_space::local_space);

            UINT_TYPE scan_result = exclusive_scan_over_group(cta, rep_count, plus<UINT_TYPE>{});

            if (thid < Nf){
                UINT_TYPE arc_count = 0;
                for (UINT_TYPE i = 0; i < FD.face_degrees[thid]; i++){
                    if(cannon_arcs[i] != UINT_TYPE_MAX){
                        triangle_numbers[thid*MaxDegree + i] = scan_result + arc_count;
                        ++arc_count;
                    }    
                }
            }
            nditem.barrier(access::fence_space::local_space);

            if (thid < Nf){
                for (UINT_TYPE i = 0; i < FD.face_degrees[thid]; i++){
                    if(cannon_arcs[i] != UINT_TYPE_MAX){
                        auto idx = triangle_numbers[thid*MaxDegree + i];
                        arc_list[idx] = {UINT_TYPE(thid), cannon_arcs[i]};
                    }
                }
            }
            nditem.barrier(access::fence_space::local_space);
//
            auto [u, v] = arc_list[thid];
           /*  if(bid == 0){
                //printf("ThreadID: %d,\n ", thid);
                sequential_print(cta, rep_count);
            } */
            auto w = FD.next(u,v);
//
            auto edge_b = FD.get_cannonical_triangle_arc(v, u); cubic_neighbours_dev[bid*N*3 + thid*3 + 0] = triangle_numbers[edge_b[0]*MaxDegree + FD.dedge_ix(edge_b[0], edge_b[1])];
            auto edge_c = FD.get_cannonical_triangle_arc(w, v); cubic_neighbours_dev[bid*N*3 + thid*3 + 1] = triangle_numbers[edge_c[0]*MaxDegree + FD.dedge_ix(edge_c[0], edge_c[1])];
            auto edge_d = FD.get_cannonical_triangle_arc(u, w); cubic_neighbours_dev[bid*N*3 + thid*3 + 2] = triangle_numbers[edge_d[0]*MaxDegree + FD.dedge_ix(edge_d[0], edge_d[1])];

        });
    });
}