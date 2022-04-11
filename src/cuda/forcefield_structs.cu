#pragma once
#include "fullerenes/gpu/isomerspace_kernel.hh"
#include <exception>

typedef IsomerspaceKernel::device_real_t device_real_t;
typedef IsomerspaceKernel::device_node_t device_node_t;
typedef GPU_REAL3 device_coord3d;
typedef GPU_NODE3 device_node3;

//Pentagons = 0
//Hexagons = 1
__constant__ device_real_t optimal_corner_cos_angles[2] = {-0.30901699437494734, -0.5}; 
__constant__ device_real_t optimal_bond_lengths[3] = {1.479, 1.458, 1.401}; 
__constant__ device_real_t optimal_dih_cos_angles[8] = {0.7946545571495363, 0.872903607049519, 0.872903607049519, 0.9410338472965512, 0.8162879359966257, 0.9139497166300941, 0.9139497166300941, 1.}; 

#if SEMINARIO_FORCE_CONSTANTS==1
__constant__ device_real_t angle_forces[2] = {207.924,216.787}; 
__constant__ device_real_t bond_forces[3] = {260.0, 353.377, 518.992}; 
__constant__ device_real_t dih_forces[4] = {35.0,65.0,3.772,270.0}; 
#else
__constant__ device_real_t angle_forces[2] = {100.0,100.0}; 
__constant__ device_real_t bond_forces[3] = {260.0,390.0,450.0}; 
__constant__ device_real_t dih_forces[4] = {35.0,65.0,85.0,270.0}; 
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
    const device_node_t* neighbours;
    __device__ DeviceFullereneGraph(const device_node_t* neighbours) : neighbours(neighbours) {}

    __device__ device_node_t dedge_ix(const device_node_t u, const device_node_t v) const{
        for (uint8_t j = 0; j < 3; j++)
            if (neighbours[u*3 + j] == v) return j;

        assert(false);
	return 0;		// Make compiler happy
    }

    __device__ device_node_t next(const device_node_t u, const device_node_t v) const{
        device_node_t j = dedge_ix(u,v);
        return neighbours[u*3 + ((j+1)%3)];
    }
    
    __device__ device_node_t prev(const device_node_t u, const device_node_t v) const{
        device_node_t j = dedge_ix(u,v);
        return neighbours[u*3 + ((j+2)%3)];
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
    

};


struct NodeGraph{
    device_node3 neighbours;
    device_node3 next_on_face;
    device_node3 prev_on_face;

    __device__ NodeGraph(const IsomerBatch& G){
        const DeviceFullereneGraph FG(&G.neighbours[blockIdx.x*blockDim.x*3]);
        this->neighbours   = {FG.neighbours[threadIdx.x*3], FG.neighbours[threadIdx.x*3 + 1], FG.neighbours[threadIdx.x*3 + 2]};
        this->next_on_face = {FG.next_on_face(threadIdx.x, FG.neighbours[threadIdx.x*3]), FG.next_on_face(threadIdx.x, FG.neighbours[threadIdx.x*3 + 1]), FG.next_on_face(threadIdx.x ,FG.neighbours[threadIdx.x*3 + 2])};
        this->prev_on_face = {FG.prev_on_face(threadIdx.x, FG.neighbours[threadIdx.x*3]), FG.prev_on_face(threadIdx.x, FG.neighbours[threadIdx.x*3 + 1]), FG.prev_on_face(threadIdx.x ,FG.neighbours[threadIdx.x*3 + 2])};
    }
};

struct Constants{
    device_coord3d f_bond;
    device_coord3d f_inner_angle;
    device_coord3d f_inner_dihedral;
    device_coord3d f_outer_angle_m;
    device_coord3d f_outer_angle_p;
    device_coord3d f_outer_dihedral;

    device_coord3d r0;
    device_coord3d angle0;
    device_coord3d outer_angle_m0;
    device_coord3d outer_angle_p0;
    device_coord3d inner_dih0;
    device_coord3d outer_dih0_a;
    device_coord3d outer_dih0_m;
    device_coord3d outer_dih0_p;
    
    __device__ __host__ uint8_t face_index(uint8_t f1, uint8_t f2, uint8_t f3){
        return f1*4 + f2*2 + f3;
    }

    __device__ Constants(const IsomerBatch& G){
        //Set pointers to start of fullerene.
        const DeviceFullereneGraph FG(&G.neighbours[blockIdx.x*blockDim.x*3]);
        device_node3 neighbours = {FG.neighbours[threadIdx.x*3], FG.neighbours[threadIdx.x*3 + 1], FG.neighbours[threadIdx.x*3 + 2]};
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
            
            uint8_t F1 = FG.face_size(threadIdx.x, d_get(neighbours, j)) - 5;
            uint8_t F2 = FG.face_size(threadIdx.x, d_get(neighbours, (j+1)%3)) -5;
            uint8_t F3 = FG.face_size(threadIdx.x, d_get(neighbours, (j+2)%3)) -5;
            
            //The faces to the right of the arcs ab, bm and bp in no particular order, from this we can deduce F4.
            uint8_t neighbour_F1 = FG.face_size(d_get(neighbours, j), FG.neighbours[d_get(neighbours, j)*3] ) -5;
            uint8_t neighbour_F2 = FG.face_size(d_get(neighbours, j), FG.neighbours[d_get(neighbours, j)*3 + 1] ) -5;
            uint8_t neighbour_F3 = FG.face_size(d_get(neighbours, j), FG.neighbours[d_get(neighbours, j)*3 + 2] ) -5;

            uint8_t F4 = neighbour_F1 + neighbour_F2 + neighbour_F3 - F1 - F3 ;
            
            //Load equillibirium distance, angles and dihedral angles from face information.
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
        }
    }   
};


