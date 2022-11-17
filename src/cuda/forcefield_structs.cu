#pragma once
#include "fullerenes/gpu/isomerspace_kernel.hh"
#include "fullerenes/gpu/reductions.cuh"
#include <exception>

//Pentagons = 0
//Hexagons = 1
constexpr __constant__ device_real_t optimal_corner_cos_angles[2] = {-0.30901699437494734, -0.5}; 
constexpr __constant__ device_real_t optimal_bond_lengths[3] = {1.479, 1.458, 1.401}; 
constexpr __constant__ device_real_t optimal_dih_cos_angles[8] = {0.7946545571495363, 0.872903607049519, 0.872903607049519, 0.9410338472965512, 0.8162879359966257, 0.9139497166300941, 0.9139497166300941, 1.}; 

#if SEMINARIO_FORCE_CONSTANTS==1
constexpr __constant__ device_real_t angle_forces[2] = {207.924,216.787}; 
constexpr __constant__ device_real_t bond_forces[3] = {260.0, 353.377, 518.992}; 
constexpr __constant__ device_real_t dih_forces[4] = {35.0,65.0,3.772,270.0}; 
constexpr __constant__ device_real_t flat_forces[3] = {0., 0., 0.};
#else
constexpr __constant__ device_real_t angle_forces[2] = {100.0,100.0}; 
constexpr __constant__ device_real_t bond_forces[3] = {260.0,390.0,450.0}; 
constexpr __constant__ device_real_t dih_forces[4] = {35.0,65.0,85.0,270.0}; 
constexpr __constant__ device_real_t flat_forces[3] = {0., 0., 0.};
#endif




template <typename T>
struct CuDeque
{
private:

    device_node_t front, back, q_size, capacity;
    T* array;

public:
    __device__ CuDeque(T* memory, const device_node_t capacity): array(memory), front(0), back(0), q_size(0), capacity(capacity) {}
    
    __device__ device_node_t size(){return q_size;}

    __device__ bool empty(){
        return q_size == 0;
    }

    __device__ bool full(){
        return q_size == capacity;
    }
    
    __device__ T pop_front(){
        if (!empty()){
            T return_val = array[front];
            front = (front + 1) % capacity ;
            q_size--;
            return return_val;
        }
        assert(false);
        return T(); //Compiler wants a return statement
    }

    __device__ T pop_back(){
        if (!empty())
        {
            T return_val = array[back];
            back = back > 0 ? back-1 : capacity-1;
            q_size--;
            return return_val;
        }
        assert(false);
        return T(); //Compiler wants a return statement
    }

    __device__ void push_back(T val){
        assert(!full());
        back = (back + 1) % capacity;
        array[back] = val;
        q_size++;
    }

    __device__ void push_front(T val){
        assert(!full());
        front = front > 0 ? front-1 : capacity-1;
        array[front] = val;
        q_size++;
    }
};



struct DeviceFullereneGraph{
    const device_node_t* cubic_neighbours;
    __device__ DeviceFullereneGraph(const device_node_t* cubic_neighbours) : cubic_neighbours(cubic_neighbours) {}

    __device__ device_node_t dedge_ix(const device_node_t u, const device_node_t v) const{
        for (uint8_t j = 0; j < 3; j++)
            if (cubic_neighbours[u*3 + j] == v) return j;

        assert(false);
	return 0;		// Make compiler happy
    }

    __device__ device_node_t next(const device_node_t u, const device_node_t v) const{
        device_node_t j = dedge_ix(u,v);
        return cubic_neighbours[u*3 + ((j+1)%3)];
    }
    
    __device__ device_node_t prev(const device_node_t u, const device_node_t v) const{
        device_node_t j = dedge_ix(u,v);
        return cubic_neighbours[u*3 + ((j+2)%3)];
    }
    
    __device__ device_node_t next_on_face(const device_node_t u, const device_node_t v) const{
        return prev(v,u);
    }

    __device__ device_node_t prev_on_face(const device_node_t u, const device_node_t v) const{
        return next(v,u);
    }

    __device__ device_node_t face_size(device_node_t u, device_node_t v) const{
        device_node_t d = 1;
        device_node_t u0 = u;
        while (v != u0)
        {
            device_node_t w = v;
            v = next_on_face(u,v);
            u = w;
            d++;
        }
        return d;
    }

    __device__ uint8_t get_face_oriented(device_node_t u, device_node_t v, device_node_t *f) const{
        constexpr int f_max = 6;
        int i = 0;
	    f[0] = u;
        while (v!=f[0] && i<f_max)
        {   
            i++;
            device_node_t w = next_on_face(u,v);
            f[i] = v;
            u = v;
            v = w;
        }
        if(i>=f_max) {assert(false); return 0;} //Compiler wants a return statement
        else return i + 1;
    }   

    __device__ device_node2 get_face_representation(device_node_t u, device_node_t v) const{
        constexpr int f_max =6;
        int i = 0;
        auto start_node = u;
        device_node2 min_edge = {u,v};
        while (v!= start_node && i < f_max){
            device_node_t w = next_on_face(u, v);
            u = v; v = w;
            if(u < min_edge.x) min_edge = {u,v};
            ++i;
        }
        //assert(next_on_face(u,v) == start_node);
        return min_edge;
    }
};

struct DeviceFullereneDual{
    const device_node_t* dual_neighbours;                   //Nf x 6
    const uint8_t* face_degrees;                            //Nf x 1

    __device__ DeviceFullereneDual(const device_node_t* dual_neighbours, const uint8_t* face_degrees) : dual_neighbours(dual_neighbours), face_degrees(face_degrees) {}

    __device__ device_node_t dedge_ix(const device_node_t u, const device_node_t v) const{
        for (uint8_t j = 0; j < face_degrees[u]; j++){
            if (dual_neighbours[u*6 + j] == v) return j;
        }

        assert(false);
	    return 0;		// Make compiler happy
    }

    __device__ device_node_t next(const device_node_t u, const device_node_t v) const{
        device_node_t j = dedge_ix(u,v);
        return dual_neighbours[u*6 + ((j+1)%face_degrees[u])];
    }
    
    __device__ device_node_t prev(const device_node_t u, const device_node_t v) const{
        device_node_t j = dedge_ix(u,v);
        return dual_neighbours[u*6 + ((j-1+face_degrees[u])%face_degrees[u])];
    }

    __device__ device_node_t next_on_face(const device_node_t u, const device_node_t v) const{
        return prev(v,u);
    }

    __device__ device_node_t prev_on_face(const device_node_t u, const device_node_t v) const{
        return next(v,u);
    }

    __device__ device_node2 get_cannonical_triangle_arc(const device_node_t u, const device_node_t v) const{
        //In a triangle u, v, w there are only 3 possible representative arcs, the cannonical arc is chosen as the one with the smalles source node.
        device_node2 min_edge = {u,v};
        device_node_t w = next(u,v);
        if (v < u && v < w) min_edge = {v, w};
        if (w < u && w < v) min_edge = {w, u};
        return min_edge;
    }

    __device__ void compute_triangle_numbers(const device_node_t u, device_node_t* triangle_numbers, device_node_t* smem){
        int represent_count = 0;
        if(u < (blockDim.x / 2 +2) ){
        for (int i = 0; i < face_degrees[u]; ++i){
            device_node2 cannon_arc = get_cannonical_triangle_arc(u,dual_neighbours[u*6 + i]);
            if (cannon_arc.x == u) {represent_count++;}
        }

        smem[u] = represent_count;
        }
        //exclusive_scan(smem, represent_count);
        //if (u < blockDim.x/2 + 2){
        //    for (size_t i = 0; i < represent_count; i++){
        //        triangle_numbers[u*6 + i] = smem[u] + i ;
        //    }
        //}
    }

    __device__ void compute_cubic_layout(const device_node_t u, const device_node_t* triangle_numbers, device_node_t* cubic_neighbours){
        if (threadIdx.x < blockDim.x / 2 +2){
        for (int i = 0; i < face_degrees[u]; ++i){
            device_node2 cannon_arc = get_cannonical_triangle_arc(u,dual_neighbours[u*6 + i]);
            if (cannon_arc.x == u) {
                device_node_t v(dual_neighbours[u*6 + i]);
                device_node_t w = prev(u,v);
                device_node2 edge_b = get_cannonical_triangle_arc(v, u); cubic_neighbours[triangle_numbers[u*6 + i]*3 + 0] = triangle_numbers[edge_b.x * 6 + dedge_ix(edge_b.x, edge_b.y)];
                device_node2 edge_c = get_cannonical_triangle_arc(w, v); cubic_neighbours[triangle_numbers[u*6 + i]*3 + 1] = triangle_numbers[edge_c.x * 6 + dedge_ix(edge_c.x, edge_c.y)];
                device_node2 edge_d = get_cannonical_triangle_arc(u, w); cubic_neighbours[triangle_numbers[u*6 + i]*3 + 2] = triangle_numbers[edge_d.x * 6 + dedge_ix(edge_d.x, edge_d.y)];
            };
        }
        }
    }

};

struct NodeGraph{
    device_node3 cubic_neighbours;
    device_node3 next_on_face;
    device_node3 prev_on_face;
    device_node_t face_nodes[6] = {UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX}; //Shape 1 x dmax , face associated with the arc a -> b , face associated with the arc a -> c, face associated with the arc a -> d
    device_node3 face_neighbours = {UINT16_MAX, UINT16_MAX, UINT16_MAX};
    unsigned char face_size = __UINT8_MAX__;
    __device__ NodeGraph(const IsomerBatch& G, const size_t isomer_idx, device_real_t* sdata){
        clear_cache(sdata,Block_Size_Pow_2);
        device_real_t* base_ptr = sdata + Block_Size_Pow_2;
        device_node_t* L = reinterpret_cast<device_node_t*>(base_ptr ); //N x 3 list of potential face IDs.
        device_node_t* A = reinterpret_cast<device_node_t*>(base_ptr) + blockDim.x * 3; //Uses cache temporarily to store face neighbours. //Nf x 6 
        const DeviceFullereneGraph FG(&G.cubic_neighbours[isomer_idx*blockDim.x*3]);
        this->cubic_neighbours   = {FG.cubic_neighbours[threadIdx.x*3], FG.cubic_neighbours[threadIdx.x*3 + 1], FG.cubic_neighbours[threadIdx.x*3 + 2]};
        this->next_on_face = {FG.next_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3]), FG.next_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3 + 1]), FG.next_on_face(threadIdx.x ,FG.cubic_neighbours[threadIdx.x*3 + 2])};
        this->prev_on_face = {FG.prev_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3]), FG.prev_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3 + 1]), FG.prev_on_face(threadIdx.x ,FG.cubic_neighbours[threadIdx.x*3 + 2])};
        int represent_count = 0;
        device_node2 rep_edges[3] = {UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX, UINT16_MAX}; 
        //If a node is representative of the j'th face then the face representation edge will be threadIdx.x -> cubic_neighbour[j]
        bool is_rep[3] = {false, false, false}; 
        device_node_t edge_idx[3] = {UINT16_MAX, UINT16_MAX, UINT16_MAX};
        for (int j = 0; j  < 3; j++){
            rep_edges[j] = FG.get_face_representation(threadIdx.x, d_get(cubic_neighbours,j));
            edge_idx[j] = FG.dedge_ix(rep_edges[j].x, rep_edges[j].y);
            if(rep_edges[j].x == threadIdx.x) {++represent_count; is_rep[j] = true;}
        }
        ex_scan<device_node_t>(reinterpret_cast<device_node_t*>(sdata), represent_count, blockDim.x);
        auto offset  = reinterpret_cast<device_node_t*>(sdata)[threadIdx.x];
        int k = 0;
        for(int j = 0; j < 3; j++){
            if(is_rep[j]){
                L[threadIdx.x*3 + j] = offset + k; 
                FG.get_face_oriented(threadIdx.x,d_get(cubic_neighbours,j), &A[offset*6 + 6*k]);
                //If the face is a pentagon assign the dummy value UINT16_MAX to the last element.
                if(FG.face_size(threadIdx.x, d_get(cubic_neighbours,j))== 5) A[offset*6 + 6*k + 5] = UINT16_MAX;
                ++k;
            }
        }
        BLOCK_SYNC
        face_neighbours = {L[rep_edges[0].x*3 + edge_idx[0]], L[rep_edges[1].x*3 + edge_idx[1]], L[rep_edges[2].x*3 + edge_idx[2]]};
        if(threadIdx.x < (blockDim.x/2) + 2){
            memcpy(&face_nodes[0], &A[threadIdx.x*6], sizeof(device_node_t)*6);
            face_size = face_nodes[5] == UINT16_MAX ? 5 : 6;
        }
        BLOCK_SYNC
    }
    __device__ NodeGraph(const IsomerBatch& G, const size_t isomer_idx){
        const DeviceFullereneGraph FG(&G.cubic_neighbours[isomer_idx*blockDim.x*3]);
        this->cubic_neighbours   = {FG.cubic_neighbours[threadIdx.x*3], FG.cubic_neighbours[threadIdx.x*3 + 1], FG.cubic_neighbours[threadIdx.x*3 + 2]};
        this->next_on_face = {FG.next_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3]), FG.next_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3 + 1]), FG.next_on_face(threadIdx.x ,FG.cubic_neighbours[threadIdx.x*3 + 2])};
        this->prev_on_face = {FG.prev_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3]), FG.prev_on_face(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3 + 1]), FG.prev_on_face(threadIdx.x ,FG.cubic_neighbours[threadIdx.x*3 + 2])};
    }
};

struct Constants{
    #if USE_CONSTANT_INDICES
    uchar4 i_f_bond;
    uchar4 i_f_inner_angle;
    uchar4 i_f_inner_dihedral;
    uchar4 i_f_outer_angle_m;
    uchar4 i_f_outer_angle_p;
    uchar4 i_f_outer_dihedral;
    uchar4 i_r0;
    uchar4 i_angle0;
    uchar4 i_outer_angle_m0;
    uchar4 i_outer_angle_p0;
    uchar4 i_inner_dih0;
    uchar4 i_outer_dih0_a;
    uchar4 i_outer_dih0_m;
    uchar4 i_outer_dih0_p;

    //Load force constants from neighbouring face information.
    constexpr INLINE device_real_t r0(const uint8_t j) const {return  optimal_bond_lengths[d_get(i_f_bond, j)];}
    constexpr INLINE device_real_t angle0(const uint8_t j)  const {return  optimal_corner_cos_angles[d_get(i_f_inner_angle, j)];}
    constexpr INLINE device_real_t inner_dih0(const uint8_t j) const {return  optimal_dih_cos_angles[d_get(i_f_inner_dihedral, j)];}
    constexpr INLINE device_real_t outer_angle_m0(const uint8_t j) const {return  optimal_corner_cos_angles[d_get(i_f_outer_angle_m, j)];}
    constexpr INLINE device_real_t outer_angle_p0(const uint8_t j) const {return  optimal_corner_cos_angles[d_get(i_f_outer_angle_p, j)];}
    constexpr INLINE device_real_t outer_dih0_a(const uint8_t j) const {return  optimal_dih_cos_angles[d_get(i_f_outer_dihedral, j)];}
    constexpr INLINE device_real_t outer_dih0_m(const uint8_t j) const {return  optimal_dih_cos_angles[d_get(i_r0, j)];}
    constexpr INLINE device_real_t outer_dih0_p(const uint8_t j) const {return  optimal_dih_cos_angles[d_get(i_angle0, j)];}
    constexpr INLINE device_real_t f_bond(const uint8_t j) const {return  bond_forces[d_get(i_outer_angle_m0, j)];}
    constexpr INLINE device_real_t f_inner_angle(const uint8_t j) const {return  angle_forces[d_get(i_outer_angle_p0, j)];}
    constexpr INLINE device_real_t f_inner_dihedral(const uint8_t j) const {return  dih_forces[d_get(i_inner_dih0, j)];}
    constexpr INLINE device_real_t f_outer_angle_m(const uint8_t j) const {return  angle_forces[d_get(i_outer_dih0_a, j)];}
    constexpr INLINE device_real_t f_outer_angle_p(const uint8_t j) const {return  angle_forces[d_get(i_outer_dih0_m, j)];}
    constexpr INLINE device_real_t f_outer_dihedral(const uint8_t j) const {return  dih_forces[d_get(i_outer_dih0_p, j)];}
    constexpr INLINE device_real_t f_flat() const {return 5e2;}
    #else
    device_coord3d f_bond;
    device_coord3d f_inner_angle;
    device_coord3d f_inner_dihedral;
    device_coord3d f_outer_angle_m;
    device_coord3d f_outer_angle_p;
    device_coord3d f_outer_dihedral;
    device_real_t f_flat = 5e2;
    
    device_coord3d r0;
    device_coord3d angle0;
    device_coord3d outer_angle_m0;
    device_coord3d outer_angle_p0;
    device_coord3d inner_dih0;
    device_coord3d outer_dih0_a;
    device_coord3d outer_dih0_m;
    device_coord3d outer_dih0_p;
    #endif

    __device__ __host__ __forceinline__ uint8_t face_index(uint8_t f1, uint8_t f2, uint8_t f3){
        return f1*4 + f2*2 + f3;
    }

    __device__ Constants(const IsomerBatch& G, const size_t isomer_idx){
        //Set pointers to start of fullerene.
        const DeviceFullereneGraph FG(&G.cubic_neighbours[isomer_idx*blockDim.x*3]);
        device_node3 cubic_neighbours = {FG.cubic_neighbours[threadIdx.x*3], FG.cubic_neighbours[threadIdx.x*3 + 1], FG.cubic_neighbours[threadIdx.x*3 + 2]};
        //       m    p
        //    f5_|   |_f4
        //   p   c    b  m
        //       \f1/
        //     f2 a f3
        //        |
        //        d
        //      m/\p
        //       f6
        
        for (uint8_t j = 0; j < 3; j++) {
            //Faces to the right of arcs ab, ac and ad.
            
            uint8_t F1 = FG.face_size(threadIdx.x, d_get(cubic_neighbours, j)) - 5;
            uint8_t F2 = FG.face_size(threadIdx.x, d_get(cubic_neighbours, (j+1)%3)) -5;
            uint8_t F3 = FG.face_size(threadIdx.x, d_get(cubic_neighbours, (j+2)%3)) -5;
            
            //The faces to the right of the arcs ab, bm and bp in no particular order, from this we can deduce F4.
            uint8_t neighbour_F1 = FG.face_size(d_get(cubic_neighbours, j), FG.cubic_neighbours[d_get(cubic_neighbours, j)*3] ) -5;
            uint8_t neighbour_F2 = FG.face_size(d_get(cubic_neighbours, j), FG.cubic_neighbours[d_get(cubic_neighbours, j)*3 + 1] ) -5;
            uint8_t neighbour_F3 = FG.face_size(d_get(cubic_neighbours, j), FG.cubic_neighbours[d_get(cubic_neighbours, j)*3 + 2] ) -5;

            uint8_t F4 = neighbour_F1 + neighbour_F2 + neighbour_F3 - F1 - F3 ;
            
            //Load equillibirium distance, angles and dihedral angles from face information.
            #if USE_CONSTANT_INDICES
            d_set(i_r0,               j,  F3 + F1);
            d_set(i_angle0,           j,  F1);
            d_set(i_inner_dih0,       j,  face_index(F1, F2 , F3));
            d_set(i_outer_angle_m0,   j,  F3);
            d_set(i_outer_angle_p0,   j,  F1);
            d_set(i_outer_dih0_a,     j,  face_index(F3, F4, F1));
            d_set(i_outer_dih0_m,     j,  face_index(F4, F1, F3));
            d_set(i_outer_dih0_p,     j,  face_index(F1, F3, F4));
            
            //Load force constants from neighbouring face information.
            d_set(i_f_bond,           j,  bond_forces[F3 + F1]);
            d_set(i_f_inner_angle,    j,  angle_forces[F1]);
            d_set(i_f_inner_dihedral, j,  dih_forces[F1 + F2 + F3]);
            d_set(i_f_outer_angle_m,  j,  angle_forces[F3]);
            d_set(i_f_outer_angle_p,  j,  angle_forces[F1]);
            d_set(i_f_outer_dihedral, j,  dih_forces[F1 + F3 + F4]);
            #else
            d_set(r0,               j,  optimal_bond_lengths[F3 + F1]);
            d_set(angle0,           j,  optimal_corner_cos_angles[F1]);
            d_set(inner_dih0,       j,  optimal_dih_cos_angles[face_index(F1, F2 , F3)]);
            d_set(outer_angle_m0,   j,  optimal_corner_cos_angles[F3]);
            d_set(outer_angle_p0,   j,  optimal_corner_cos_angles[F1]);
            d_set(outer_dih0_a,     j,  optimal_dih_cos_angles[face_index(F3, F4, F1)]);
            d_set(outer_dih0_m,     j,  optimal_dih_cos_angles[face_index(F4, F1, F3)]);
            d_set(outer_dih0_p,     j,  optimal_dih_cos_angles[face_index(F1, F3, F4)]);
            
            //Load force constants from neighbouring face information.
            d_set(f_bond,           j,  bond_forces[F3 + F1]);
            d_set(f_inner_angle,    j,  angle_forces[F1]);
            d_set(f_inner_dihedral, j,  dih_forces[F1 + F2 + F3]);
            d_set(f_outer_angle_m,  j,  angle_forces[F3]);
            d_set(f_outer_angle_p,  j,  angle_forces[F1]);
            d_set(f_outer_dihedral, j,  dih_forces[F1 + F3 + F4]);
            #endif
        }
    }   
};


