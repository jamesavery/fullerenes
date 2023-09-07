#include <array>

template <typename K>
struct DeviceCubicGraph{
    static_assert(std::is_integral<K>::value, "K must be integral");
    const accessor<K, 1, access::mode::read> cubic_neighbours;
    const size_t offset;

    inline K operator[](const K i) const{
        return cubic_neighbours[i + offset];
    }

    DeviceCubicGraph(const accessor<K, 1, access::mode::read> cubic_neighbours, size_t offset) : cubic_neighbours(cubic_neighbours), offset(offset) {}

    /** @brief Find the index of the neighbour v in the list of neighbours of u
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: index of v in the list of neighbours of u
    */
    K dedge_ix(const K u, const K v) const{
        for (uint8_t j = 0; j < 3; j++)
            if ((*this)[u*3 + j] == v) return j;

        assert(false);
	return 0;		// Make compiler happy
    }

    /** @brief Find the next neighbour in the clockwise order around u
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: Next neighbour of u after v in the clockwise order
    */
    K next(const K u, const K v) const{
        K j = dedge_ix(u,v);
        return (*this)[u*3 + ((j+1)%3)];
    }
    
    /** @brief Find the previous neighbour in the clockwise order around u
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: Previous neighbour of u before v in the clockwise order
    */
    K prev(const K u, const K v) const{
        K j = dedge_ix(u,v);
        return (*this)[u*3 + ((j+2)%3)];
    }
    
    /** @brief Find the next node in the face represented by the arc (u,v)
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: The next node in the face
    */
    K next_on_face(const K u, const K v) const{
        return prev(v,u);
    }

    /** @brief Find the previous node in the face represented by the arc (u,v)
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: The previous node in the face
    */
    K prev_on_face(const K u, const K v) const{
        return next(v,u);
    }

    /** @brief Find the size of the face represented by the arc (u,v)
    // @param u: source node in the arc (u,v)
    // @param v: target node in the arc (u,v)
    // @return: size of the face
    */
    K face_size(K u, K v) const{
        K d = 1;
        K u0 = u;
        while (v != u0)
        {
            K w = v;
            v = next_on_face(u,v);
            u = w;
            d++;
        }
        return d;
    }

    uint8_t get_face_oriented(K u, K v, K *f) const{
        constexpr int f_max = 6;
        int i = 0;
	    f[0] = u;
        while (v!=f[0] && i<f_max)
        {   
            i++;
            K w = next_on_face(u,v);
            f[i] = v;
            u = v;
            v = w;
        }
        if(i>=f_max) {assert(false); return 0;} //Compiler wants a return statement
        else return i + 1;
    }   

    std::array<K,2> get_face_representation(K u, K v) const{
        constexpr int f_max =6;
        int i = 0;
        auto start_node = u;
        std::array<K,2> min_edge = {u,v};
        while (v!= start_node && i < f_max){
            K w = next_on_face(u, v);
            u = v; v = w;
            if(u < min_edge[0]) min_edge = {u,v};
            ++i;
        }
        //assert(next_on_face(u,v) == start_node);
        return min_edge;
    }
};