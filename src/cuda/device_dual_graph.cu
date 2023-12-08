struct DeviceDualGraph{
    const device_node_t* dual_neighbours;                   //Nf x 6
    const uint8_t* face_degrees;                            //Nf x 1

    __device__ DeviceDualGraph(const device_node_t* dual_neighbours, const uint8_t* face_degrees) : dual_neighbours(dual_neighbours), face_degrees(face_degrees) {}

    __device__ device_node_t arc_ix(const device_node_t u, const device_node_t v) const{
        for (uint8_t j = 0; j < face_degrees[u]; j++){
            if (dual_neighbours[u*6 + j] == v) return j;
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
    __device__ device_node_t next(const device_node_t u, const device_node_t v) const{
        device_node_t j = arc_ix(u,v);
        return dual_neighbours[u*6 + ((j+1)%face_degrees[u])];
    }
    
    /**
     * @brief returns the prev node in the clockwise order around u
     * @param v the current node around u
     * @param u the node around which the search is performed
     * @return the previous node in the clockwise order around u
     */
    __device__ device_node_t prev(const device_node_t u, const device_node_t v) const{
        device_node_t j = arc_ix(u,v);
        return dual_neighbours[u*6 + ((j-1+face_degrees[u])%face_degrees[u])];
    }

    /**
     * @brief Find the node that comes next on the face. given by the edge (u,v)
     * @param u Source of the edge.
     * @param v Destination node.
     * @return The node that comes next on the face.
     */
    __device__ device_node_t next_on_face(const device_node_t u, const device_node_t v) const{
        return prev(v,u);
    }

    /**
     * @brief Find the node that comes next on the face. given by the edge (u,v)
     * @param u Source of the edge.
     * @param v Destination node.
     * @return The node that comes next on the face.
     */
    __device__ device_node_t prev_on_face(const device_node_t u, const device_node_t v) const{
        return next(v,u);
    }

    /**
     * @brief Finds the cannonical triangle arc of the triangle (u,v,w)
     * 
     * @param u source node
     * @param v target node
     * @return cannonical triangle arc 
     */
    __device__ device_node2 get_cannonical_triangle_arc(const device_node_t u, const device_node_t v) const{
        //In a triangle u, v, w there are only 3 possible representative arcs, the cannonical arc is chosen as the one with the smalles source node.
        device_node2 min_edge = {u,v};
        device_node_t w = next(u,v);
        if (v < u && v < w) min_edge = {v, w};
        if (w < u && w < v) min_edge = {w, u};
        return min_edge;
    }

    /** @brief Computes the number of triangles each node is a part of. 
     *
     * @param u The node for which to compute how many triangles it is representative of.
     * @param smem A pointer to shared memory.
     */
    __device__ void compute_triangle_numbers(const device_node_t u, device_node_t* triangle_numbers, device_node_t* smem){
        int represent_count = 0;
        if(u < (blockDim.x / 2 +2) ){
        for (int i = 0; i < face_degrees[u]; ++i){
            device_node2 cannon_arc = get_cannonical_triangle_arc(u,dual_neighbours[u*6 + i]);
            if (cannon_arc[0] == u) {represent_count++;}
        }

        smem[u] = represent_count;
        }
    }

    __device__ void compute_cubic_layout(const device_node_t u, const device_node_t* triangle_numbers, device_node_t* cubic_neighbours){
        if (threadIdx.x < blockDim.x / 2 +2){
        for (int i = 0; i < face_degrees[u]; ++i){
            device_node2 cannon_arc = get_cannonical_triangle_arc(u,dual_neighbours[u*6 + i]);
            if (cannon_arc[0] == u) {
                device_node_t v(dual_neighbours[u*6 + i]);
                device_node_t w = prev(u,v);
                device_node2 edge_b = get_cannonical_triangle_arc(v, u); cubic_neighbours[triangle_numbers[u*6 + i]*3 + 0] = triangle_numbers[edge_b[0] * 6 + arc_ix(edge_b[0], edge_b[1])];
                device_node2 edge_c = get_cannonical_triangle_arc(w, v); cubic_neighbours[triangle_numbers[u*6 + i]*3 + 1] = triangle_numbers[edge_c[0] * 6 + arc_ix(edge_c[0], edge_c[1])];
                device_node2 edge_d = get_cannonical_triangle_arc(u, w); cubic_neighbours[triangle_numbers[u*6 + i]*3 + 2] = triangle_numbers[edge_d[0] * 6 + arc_ix(edge_d[0], edge_d[1])];
            };
        }
        }
    }

};