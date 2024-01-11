template <typename K>
struct NodeNeighbours{
    std::array<K,3> cubic_neighbours;
    std::array<K,3> next_on_face;
    std::array<K,3> prev_on_face;
    K face_nodes[6] = {UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX}; //Shape 1 x dmax , face associated with the arc a -> b , face associated with the arc a -> c, face associated with the arc a -> d
    std::array<K,3> face_neighbours = {UINT16_MAX, UINT16_MAX, UINT16_MAX};
    unsigned char face_size = __UINT8_MAX__;
    
    /**
     * @brief  This constructor computes the neighbours, outer neighbours, face neighbours for the first Nf threads it stores the nodes that are part of the threadIdx.x^th face.
     * @param  G: The IsomerBatch object that contains the graph data
     * @param  isomer_idx: The index of the isomer that the thread is a part of.
     * @param  sdata: Pointer to shared memory.
     * @return NodeNeighbours object.
     */
    template <Device U, typename T>
    __device__ NodeNeighbours(const IsomerBatch<U>& G, const size_t isomer_idx, T* sdata){
        TEMPLATE_TYPEDEFS(T,K);
        clear_cache(sdata,Block_Size_Pow_2);
        real_t* base_ptr = sdata + Block_Size_Pow_2;
        node_t* L = reinterpret_cast<node_t*>(base_ptr ); //N x 3 list of potential face IDs.
        node_t* A = reinterpret_cast<node_t*>(base_ptr) + blockDim.x * 3; //Uses cache temporarily to store face neighbours. //Nf x 6 
        const DeviceCubicGraph FG(&G.cubic_neighbours[isomer_idx*blockDim.x*3]);
        this->cubic_neighbours   = {FG.cubic_neighbours[threadIdx.x*3], FG.cubic_neighbours[threadIdx.x*3 + 1], FG.cubic_neighbours[threadIdx.x*3 + 2]};
        this->next_on_face = {FG.next_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3]), FG.next_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3 + 1]), FG.next_on_face(threadIdx.x ,FG.cubic_neighbours[threadIdx.x*3 + 2])};
        this->prev_on_face = {FG.prev_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3]), FG.prev_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3 + 1]), FG.prev_on_face(threadIdx.x ,FG.cubic_neighbours[threadIdx.x*3 + 2])};
        int represent_count = 0;
        node2 rep_edges[3] = {{UINT16_MAX, UINT16_MAX}, {UINT16_MAX, UINT16_MAX}, {UINT16_MAX, UINT16_MAX}}; 
        //If a node is representative of the j'th face then the face representation edge will be threadIdx.x -> cubic_neighbour[j]
        bool is_rep[3] = {false, false, false}; 
        node_t edge_idx[3] = {UINT16_MAX, UINT16_MAX, UINT16_MAX};
        for (int j = 0; j  < 3; j++){
            rep_edges[j] = FG.get_face_representation(threadIdx.x, cubic_neighbours[j]);
            edge_idx[j] = FG.arc_ix(rep_edges[j][0], rep_edges[j][1]);
            if(rep_edges[j][0] == threadIdx.x) {++represent_count; is_rep[j] = true;}
        }
        ex_scan<node_t>(reinterpret_cast<node_t*>(sdata), represent_count, blockDim.x);
        auto offset  = reinterpret_cast<node_t*>(sdata)[threadIdx.x];
        int k = 0;
        for(int j = 0; j < 3; j++){
            if(is_rep[j]){
                L[threadIdx.x*3 + j] = offset + k; 
                FG.get_face_oriented(threadIdx.x,cubic_neighbours[j], &A[offset*6 + 6*k]);
                //If the face is a pentagon assign the dummy value UINT16_MAX to the last element.
                if(FG.face_size(threadIdx.x, cubic_neighbours[j])== 5) A[offset*6 + 6*k + 5] = UINT16_MAX;
                ++k;
            }
        }
        BLOCK_SYNC
        face_neighbours = {L[rep_edges[0][0]*3 + edge_idx[0]], L[rep_edges[1][0]*3 + edge_idx[1]], L[rep_edges[2][0]*3 + edge_idx[2]]};
        if(threadIdx.x < (blockDim.x/2) + 2){
            memcpy(&face_nodes[0], &A[threadIdx.x*6], sizeof(node_t)*6);
            face_size = face_nodes[5] == UINT16_MAX ? 5 : 6;
        }
        BLOCK_SYNC
    }
/**
* @brief Constructor for a NodeNeighbours object, which contains the neighbours of a node in the graph and outer neighbours.
* @param G All isomer graphs in the batch.
* @param isomer_idx The index of the isomer to initialize based on.
*/
template <Device U>
__device__ NodeNeighbours(const IsomerBatch<U>& G, const size_t isomer_idx){
        const DeviceCubicGraph FG(&G.cubic_neighbours[isomer_idx*blockDim.x*3]);
        this->cubic_neighbours   = {FG.cubic_neighbours[threadIdx.x*3], FG.cubic_neighbours[threadIdx.x*3 + 1], FG.cubic_neighbours[threadIdx.x*3 + 2]};
        this->next_on_face = {FG.next_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3]), FG.next_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3 + 1]), FG.next_on_face(threadIdx.x ,FG.cubic_neighbours[threadIdx.x*3 + 2])};
        this->prev_on_face = {FG.prev_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3]), FG.prev_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3 + 1]), FG.prev_on_face(threadIdx.x ,FG.cubic_neighbours[threadIdx.x*3 + 2])};
    }
};