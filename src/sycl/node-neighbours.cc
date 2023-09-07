#include <array>
#include <limits>
#include <stdint.h>
template <typename K>
struct NodeNeighbours{
    std::array<K,3> cubic_neighbours;
    std::array<K,3> next_on_face;
    std::array<K,3> prev_on_face;
    std::array<K,6> face_nodes; //= {0, 0, 0, 0, 0, 0}; //Shape 1 x dmax , face associated with the arc a -> b , face associated with the arc a -> c, face associated with the arc a -> d
    std::array<K,3> face_neighbours;
    K face_size = UINT16_MAX;
    
    /**
     * @brief  This constructor computes the neighbours, outer neighbours, face neighbours for the first Nf threads it stores the nodes that are part of the threadIdx.x^th face.
     * @param  G: The IsomerBatch object that contains the graph data
     * @param  isomer_idx: The index of the isomer that the thread is a part of.
     * @param  sdata: Pointer to shared memory.
     * @return NodeNeighbours object.
     */
    NodeNeighbours(cl::sycl::group<1> cta, const sycl::accessor<K, 1, access::mode::read>& cubic_neighbours_acc, K* sdata){
        INT_TYPEDEFS(K);
        face_nodes.fill(UINT16_MAX);
        face_neighbours.fill(UINT16_MAX);
        auto tid = cta.get_local_linear_id();
        auto isomer_idx = cta.get_group_linear_id();
        auto bdim = cta.get_local_linear_range();
        node_t* L = reinterpret_cast<node_t*>(sdata); //N x 3 list of potential face IDs.
        node_t* A = reinterpret_cast<node_t*>(sdata) + bdim * 3; //Uses cache temporarily to store face neighbours. //Nf x 6 
        const DeviceCubicGraph FG(cubic_neighbours_acc, isomer_idx*bdim*3);
        this->cubic_neighbours   = {FG[tid*3], FG[tid*3 + 1], FG[tid*3 + 2]};
        this->next_on_face = {FG.next_on_face(tid, cubic_neighbours[tid*3]), FG.next_on_face(tid, cubic_neighbours[tid*3 + 1]), FG.next_on_face(tid ,cubic_neighbours[tid*3 + 2])};
        this->prev_on_face = {FG.prev_on_face(tid, cubic_neighbours[tid*3]), FG.prev_on_face(tid, cubic_neighbours[tid*3 + 1]), FG.prev_on_face(tid ,cubic_neighbours[tid*3 + 2])};
        int represent_count = 0;
        node2 rep_edges[3] = {{UINT16_MAX, UINT16_MAX}, {UINT16_MAX, UINT16_MAX}, {UINT16_MAX, UINT16_MAX}}; 
        //If a node is representative of the j'th face then the face representation edge will be tid -> cubic_neighbour[j]
        bool is_rep[3] = {false, false, false}; 
        node_t edge_idx[3] = {UINT16_MAX, UINT16_MAX, UINT16_MAX};
        for (int j = 0; j  < 3; j++){
            rep_edges[j] = FG.get_face_representation(tid, cubic_neighbours[j]);
            edge_idx[j] = FG.dedge_ix(rep_edges[j][0], rep_edges[j][1]);
            if(rep_edges[j][0] == tid) {++represent_count; is_rep[j] = true;}
        }
        //ex_scan<node_t>(reinterpret_cast<node_t*>(sdata), represent_count, bdim);
        auto offset  = sycl::exclusive_scan_over_group(cta, represent_count, sycl::plus<node_t>{});  //reinterpret_cast<node_t*>(sdata)[tid];
        int k = 0;
        for(int j = 0; j < 3; j++){
            if(is_rep[j]){
                L[tid*3 + j] = offset + k; 
                FG.get_face_oriented(tid,cubic_neighbours[j], &A[offset*6 + 6*k]);
                //If the face is a pentagon assign the dummy value UINT16_MAX to the last element.
                if(FG.face_size(tid, cubic_neighbours[j])== 5) A[offset*6 + 6*k + 5] = UINT16_MAX;
                ++k;
            }
        }
        sycl::group_barrier(cta);
        face_neighbours = {L[rep_edges[0][0]*3 + edge_idx[0]], L[rep_edges[1][0]*3 + edge_idx[1]], L[rep_edges[2][0]*3 + edge_idx[2]]};
        if(tid < (bdim/2) + 2){
            memcpy(&face_nodes[0], &A[tid*6], sizeof(node_t)*6);
            face_size = face_nodes[5] == UINT16_MAX ? 5 : 6;
        }
        sycl::group_barrier(cta);
    }
/**
* @brief Constructor for a NodeNeighbours object, which contains the neighbours of a node in the graph and outer neighbours.
* @param G All isomer graphs in the batch.
* @param isomer_idx The index of the isomer to initialize based on.
*/

NodeNeighbours(const sycl::accessor<K, 1, access::mode::read>& cubic_neighbours_acc, sycl::group<1>& cta){
        int tid = cta.get_local_linear_id();
        int isomer_idx = cta.get_group_linear_id();
        int blockDim = cta.get_local_linear_range();
        const DeviceCubicGraph FG(cubic_neighbours_acc, isomer_idx*blockDim*3);
        this->cubic_neighbours   = {FG[tid*3], FG[tid*3 + 1], FG[tid*3 + 2]};
        this->next_on_face = {FG.next_on_face(tid, FG[tid*3]), FG.next_on_face(tid, FG[tid*3 + 1]), FG.next_on_face(tid ,FG[tid*3 + 2])};
        this->prev_on_face = {FG.prev_on_face(tid, FG[tid*3]), FG.prev_on_face(tid, FG[tid*3 + 1]), FG.prev_on_face(tid ,FG[tid*3 + 2])};
    }
};