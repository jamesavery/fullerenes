#pragma once
#include "coord3d.cu"
#include "fullerenes/gpu/isomerspace_forcefield.hh"
#include <exception>


typedef IsomerspaceForcefield::device_real_t device_real_t;
typedef IsomerspaceForcefield::device_node_t device_node_t;
typedef GPU_REAL3 device_coord3d;
typedef GPU_NODE3 device_node3;

void IsomerspaceForcefield::GenericStruct::allocate(IsomerspaceForcefield::GenericStruct& G, size_t N, const size_t batch_size, const BufferType buffer_type){
    if((!G.allocated)){
        G.buffer_type = buffer_type;
        G.batch_size  = batch_size; 
        G.N           = N; 
        size_t num_elements = N*batch_size;
        if (buffer_type == DEVICE_BUFFER){
            for (size_t i = 0; i < G.pointers.size(); i++) {
                cudaMalloc(get<1>(G.pointers[i]), num_elements* get<2>(G.pointers[i])); 
            }
            printLastCudaError("Failed to allocate device struct");
        }else{
            for (size_t i = 0; i < G.pointers.size(); i++) {
                *get<1>(G.pointers[i])= malloc(num_elements* get<2>(G.pointers[i])); 
            }
        }        
        G.allocated = true;
    }
}

void IsomerspaceForcefield::GenericStruct::free(IsomerspaceForcefield::GenericStruct& G){
    if(G.allocated){
        if (G.buffer_type == DEVICE_BUFFER){    
            for (size_t i = 0; i < G.pointers.size(); i++) {
                cudaFree(*get<1>(G.pointers[i]));
            }
            printLastCudaError("Failed to free device struct"); 
        } else{
            for (size_t i = 0; i < G.pointers.size(); i++) {
                std::free(*get<1>(G.pointers[i])); 
            }
        }
        G.allocated = false;
    }
}
void IsomerspaceForcefield::GenericStruct::copy(IsomerspaceForcefield::GenericStruct& destination, const IsomerspaceForcefield::GenericStruct& source, const size_t num_isomers){
    if(num_isomers > 0){
    for (size_t i = 0; i < destination.pointers.size(); i++)
    {
        cudaMemcpy(*(get<1>(destination.pointers[i])) , *(get<1>(source.pointers[i])), get<2>(source.pointers[i])*source.N*num_isomers, cudaMemcpyKind((source.buffer_type+1) + (source.buffer_type + destination.buffer_type)/2));
    }
    }
    else{
        std::cout << "WARNING: Call to copy made for 0 isomers.";
    }
    printLastCudaError("Failed to copy struct");
}


//Pentagons = 0
//Hexagons = 1
//PPP = 0, {HPP, PHP, PPH} = 1, {PHH, HPH, HHP} = 2, {HHH} = 3
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


struct NodeGraph{
    const device_node3 neighbours;
    const device_node3 next_on_face;
    const device_node3 prev_on_face;
    __device__ NodeGraph(const device_node3& neighbours, const device_node3& next_on_face, const device_node3& prev_on_face) : 
        neighbours(neighbours), next_on_face(next_on_face), prev_on_face(prev_on_face) {}


    __device__ NodeGraph(const IsomerspaceForcefield::IsomerspaceGraph& G):  neighbours(MAKE_NODE3(G.neighbours[(threadIdx.x + blockDim.x*blockIdx.x)*3],G.neighbours[(threadIdx.x + blockDim.x*blockIdx.x)*3 + 1],G.neighbours[(threadIdx.x + blockDim.x*blockIdx.x)*3 + 2])),
                                                                        next_on_face(MAKE_NODE3(G.next_on_face[(threadIdx.x + blockDim.x*blockIdx.x)*3],G.next_on_face[(threadIdx.x + blockDim.x*blockIdx.x)*3 + 1],G.next_on_face[(threadIdx.x + blockDim.x*blockIdx.x)*3 + 2])),
                                                                        prev_on_face(MAKE_NODE3(G.prev_on_face[(threadIdx.x + blockDim.x*blockIdx.x)*3],G.prev_on_face[(threadIdx.x + blockDim.x*blockIdx.x)*3 + 1],G.prev_on_face[(threadIdx.x + blockDim.x*blockIdx.x)*3 + 2])){

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

    __device__ Constants(const IsomerspaceForcefield::IsomerspaceGraph& G){
        //Set pointers to start of fullerene.
        size_t offset = blockDim.x*blockIdx.x;
        device_node_t* neighbours = G.neighbours + offset*3;
        uint8_t* face_right = G.face_right + offset*3;

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
            uint8_t F1 = face_right[threadIdx.x * 3 + j] - 5;
            uint8_t F2 = face_right[threadIdx.x * 3 + (1+ j)%3] - 5;
            uint8_t F3 = face_right[threadIdx.x * 3 + (2 + j)%3] - 5;
            
            //The faces to the right of the arcs ab, bm and bp in no particular order, from this we can deduce F4.
            uint8_t neighbour_F1 = face_right[neighbours[threadIdx.x * 3 + j]*3 ] - 5;
            uint8_t neighbour_F2 = face_right[neighbours[threadIdx.x * 3 + j]*3 + 1 ] - 5;
            uint8_t neighbour_F3 = face_right[neighbours[threadIdx.x * 3 + j]*3 + 2] - 5;

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



