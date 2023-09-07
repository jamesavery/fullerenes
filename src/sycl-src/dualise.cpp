#include <CL/sycl.hpp>
#include "forcefield_includes.cpp"

#define UINT_TYPE_MAX std::numeric_limits<UINT_TYPE>::max()

template<int MaxDegree, typename K>
struct DeviceDualGraph{
    //Check that K is integral
    INT_TYPEDEFS(K);

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
    node2 get_cannonical_triangle_arc(const K u, const K v) const{
        //In a triangle u, v, w there are only 3 possible representative arcs, the cannonical arc is chosen as the one with the smalles source node.
        node2 min_edge = {u,v};
        K w = next(u,v);
        if (v < u && v < w) min_edge = {v, w};
        if (w < u && w < v) min_edge = {w, u};
        return min_edge;
    }
};

template <typename T, typename K>
void dualise(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy){
    constexpr int MaxDegree = 6;
    INT_TYPEDEFS(K);
    if(policy == LaunchPolicy::SYNC) Q.wait();
    Q.submit([&](handler &h) {
        auto N = batch.N();
        auto Nf = batch.Nf();
        
        //Create local accessors to shared memory
        local_accessor<node_t, 1>   triangle_numbers(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_neighbours(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_degrees(Nf, h);
        local_accessor<node2, 1>    arc_list(N, h);

        //Create device accessors
        accessor<node_t, 1>     cubic_neighbours_dev(batch.cubic_neighbours, h, write_only);
        accessor<node_t, 1>     face_degrees_dev(batch.face_degrees, h, read_only);
        accessor<node_t, 1>     dual_neighbours_dev(batch.dual_neighbours, h, read_only);

        h.parallel_for<class dualise>(sycl::nd_range(sycl::range{N*batch.capacity()}, sycl::range{N}), [=](nd_item<1> nditem) {
            auto cta = nditem.get_group();
            auto thid = nditem.get_local_linear_id();
            auto bid = nditem.get_group_linear_id();
            cta.async_work_group_copy(cached_neighbours.get_pointer(), cubic_neighbours_dev.get_pointer() + bid*Nf*MaxDegree, Nf*MaxDegree);
            cta.async_work_group_copy(cached_degrees.get_pointer(), face_degrees_dev.get_pointer() + bid*Nf, Nf);
            DeviceDualGraph<MaxDegree, node_t> FD(cached_neighbours.get_pointer(), cached_degrees.get_pointer());
            node_t cannon_arcs[MaxDegree]; memset(cannon_arcs, std::numeric_limits<node_t>::max(), MaxDegree*sizeof(node_t));
            node_t rep_count  = 0;
            sycl::group_barrier(cta);
            if (thid < Nf){
                for (node_t i = 0; i < FD.face_degrees[thid]; i++){
                    auto cannon_arc = FD.get_cannonical_triangle_arc(thid, FD.dual_neighbours[thid*MaxDegree + i]);
                    if (cannon_arc[0] == thid){
                        cannon_arcs[i] = cannon_arc[1];
                        rep_count++;
                    }
                }
            }
            sycl::group_barrier(cta);

            node_t scan_result = exclusive_scan_over_group(cta, rep_count, plus<node_t>{});

            if (thid < Nf){
                node_t arc_count = 0;
                for (node_t i = 0; i < FD.face_degrees[thid]; i++){
                    if(cannon_arcs[i] != std::numeric_limits<node_t>::max()){
                        triangle_numbers[thid*MaxDegree + i] = scan_result + arc_count;
                        ++arc_count;
                    }    
                }
            }
            sycl::group_barrier(cta);

            if (thid < Nf){
                for (node_t i = 0; i < FD.face_degrees[thid]; i++){
                    if(cannon_arcs[i] != std::numeric_limits<node_t>::max()){
                        auto idx = triangle_numbers[thid*MaxDegree + i];
                        arc_list[idx] = {node_t(thid), cannon_arcs[i]};
                    }
                }
            }
            sycl::group_barrier(cta);
//
            auto [u, v] = arc_list[thid];
            auto w = FD.next(u,v);
//
            auto edge_b = FD.get_cannonical_triangle_arc(v, u); cubic_neighbours_dev[bid*N*3 + thid*3 + 0] = triangle_numbers[edge_b[0]*MaxDegree + FD.dedge_ix(edge_b[0], edge_b[1])];
            auto edge_c = FD.get_cannonical_triangle_arc(w, v); cubic_neighbours_dev[bid*N*3 + thid*3 + 1] = triangle_numbers[edge_c[0]*MaxDegree + FD.dedge_ix(edge_c[0], edge_c[1])];
            auto edge_d = FD.get_cannonical_triangle_arc(u, w); cubic_neighbours_dev[bid*N*3 + thid*3 + 2] = triangle_numbers[edge_d[0]*MaxDegree + FD.dedge_ix(edge_d[0], edge_d[1])];

        });
    });
    if(policy == LaunchPolicy::SYNC) Q.wait();
}