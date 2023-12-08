template <typename node_t>
struct HostCubicGraph{
    const node_t* cubic_neighbours;
    HostCubicGraph(const node_t* cubic_neighbours) : cubic_neighbours(cubic_neighbours) {}

    /** @brief Find the index of the neighbour v in the list of neighbours of u
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: index of v in the list of neighbours of u
    */
    node_t arc_ix(const node_t u, const node_t v) const{
        for (uint8_t j = 0; j < 3; j++)
            if (cubic_neighbours[u*3 + j] == v) return j;

        assert(false);
	return 0;		// Make compiler happy
    }

    /** @brief Find the next neighbour in the clockwise order around u
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: Next neighbour of u after v in the clockwise order
    */
    node_t next(const node_t u, const node_t v) const{
        node_t j = arc_ix(u,v);
        return cubic_neighbours[u*3 + ((j+1)%3)];
    }
    
    /** @brief Find the previous neighbour in the clockwise order around u
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: Previous neighbour of u before v in the clockwise order
    */
    node_t prev(const node_t u, const node_t v) const{
        node_t j = arc_ix(u,v);
        return cubic_neighbours[u*3 + ((j+2)%3)];
    }
    
    /** @brief Find the next node in the face represented by the arc (u,v)
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: The next node in the face
    */
    node_t next_on_face(const node_t u, const node_t v) const{
        return prev(v,u);
    }

    /** @brief Find the previous node in the face represented by the arc (u,v)
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: The previous node in the face
    */
    node_t prev_on_face(const node_t u, const node_t v) const{
        return next(v,u);
    }

    /** @brief Find the size of the face represented by the arc (u,v)
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: size of the face
    */
    node_t face_size(node_t u, node_t v) const{
        node_t d = 1;
        node_t u0 = u;
        while (v != u0)
        {
            node_t w = v;
            v = next_on_face(u,v);
            u = w;
            d++;
        }
        return d;
    }

    uint8_t get_face_oriented(node_t u, node_t v, node_t *f) const{
        constexpr int f_max = 6;
        int i = 0;
	    f[0] = u;
        while (v!=f[0] && i<f_max)
        {   
            i++;
            node_t w = next_on_face(u,v);
            f[i] = v;
            u = v;
            v = w;
        }
        if(i>=f_max) {assert(false); return 0;} //Compiler wants a return statement
        else return i + 1;
    }   

    device_node2 get_face_representation(node_t u, node_t v) const{
        constexpr int f_max =6;
        int i = 0;
        auto start_node = u;
        device_node2 min_edge = {u,v};
        while (v!= start_node && i < f_max){
            node_t w = next_on_face(u, v);
            u = v; v = w;
            if(u < min_edge.x) min_edge = {u,v};
            ++i;
        }
        //assert(next_on_face(u,v) == start_node);
        return min_edge;
    }
};