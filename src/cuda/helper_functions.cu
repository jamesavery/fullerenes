#include "coord3d.cu"
#include "coord3d_aligned.cu"
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include "cuda_runtime.h"
#include <assert.h>
#include<iostream>
#include <fstream>

#define __HD__ __device__ __host__ 
namespace cg = cooperative_groups;



template <typename T>
void copy_and_append(T* memory, const T* fullerene, size_t N){
    for (size_t i = 0; i < N; i++)
    {
        memory[i] = fullerene[i];
    }
}

template <typename T>
T* synthetic_array(size_t N, const size_t num_molecules, const T* fullerene){
    size_t array_size = N;
    if (sizeof(T) != sizeof(device_coord3d))
    {
        array_size *= 3;
    }
    T* storage_array = new T[array_size*num_molecules];
    for (size_t i = 0; i < num_molecules; i++)
    {
        copy_and_append(&storage_array[array_size*i],fullerene,array_size);
    }
    return storage_array;
}


__device__ void align16(device_coord3d* input, coord3d_a* output, size_t N){
    cg::sync(cg::this_grid());
    output[threadIdx.x] = {input[threadIdx.x].x, input[threadIdx.x].y, input[threadIdx.x].z, 0};
    cg::sync(cg::this_grid());
}

template <typename T>
__device__ void pointerswap(T **r, T **s)
{
    T *pSwap = *r;
    *r = *s;
    *s = pSwap;
    return;
}



//Pentagons = 0
//Hexagons = 1
//PPP = 0, {HPP, PHP, PPH} = 1, {PHH, HPH, HHP} = 2, {HHH} = 3
__constant__ device_real_t optimal_corner_cos_angles[2] = {-0.30901699437494734, -0.5}; 
__constant__ device_real_t optimal_bond_lengths[3] = {1.479, 1.458, 1.401}; 
__constant__ device_real_t optimal_dih_cos_angles[8] = {0.7946545571495363, 0.872903607049519, 0.872903607049519, 0.9410338472965512, 0.8162879359966257, 0.9139497166300941, 0.9139497166300941, 1.}; 

__constant__ device_real_t angle_forces[2] = {207.924,216.787}; 
__constant__ device_real_t bond_forces[3] = {260.0, 353.377, 518.992}; 
__constant__ device_real_t dih_forces[4] = {35.0,65.0,3.772,270.0}; 

struct BookkeepingData{
    const device_node_t* neighbours;
    const uint8_t* face_right;
    const device_node_t* next_on_face;
    const device_node_t* prev_on_face;
    __device__ BookkeepingData(const device_node_t* neighbours, const uint8_t* face_right, const device_node_t* next_on_face, const device_node_t* prev_on_face) : 
        neighbours(neighbours), face_right(face_right), next_on_face(next_on_face), prev_on_face(prev_on_face) {}

    __device__ BookkeepingData(const BookkeepingData &b1): neighbours(b1.neighbours), face_right(b1.face_right), next_on_face(b1.next_on_face), prev_on_face(b1.prev_on_face){}
};




template <typename T>
struct Constants{
    const T f_bond;
    const T f_inner_angle;
    const T f_inner_dihedral;
    const T f_outer_angle_m;
    const T f_outer_angle_p;
    const T f_outer_dihedral;

    const T r0;
    const T angle0;
    const T outer_angle_m0;
    const T outer_angle_p0;
    const T inner_dih0;
    const T outer_dih0_a;
    const T outer_dih0_m;
    const T outer_dih0_p;
    
    __device__ Constants(const T f_bond, const T f_inner_angle, const T f_inner_dihedral, const T f_outer_angle_m, const T f_outer_angle_p, const T f_outer_dihedral,
                            const T r0, const T angle0, const T outer_angle_m0, const T outer_angle_p0, const T inner_dih0, const T outer_dih0_a, const T outer_dih0_m, const T outer_dih0_p): f_bond(f_bond), f_inner_angle(f_inner_angle),
                            f_inner_dihedral(f_inner_dihedral), f_outer_angle_m(f_outer_angle_m), f_outer_angle_p(f_outer_angle_p), f_outer_dihedral(f_outer_dihedral), r0(r0), angle0(angle0), outer_angle_m0(outer_angle_m0), outer_angle_p0(outer_angle_p0),
                            inner_dih0(inner_dih0), outer_dih0_a(outer_dih0_a), outer_dih0_m(outer_dih0_m), outer_dih0_p(outer_dih0_p) {}


    __device__ Constants(const Constants<T> &c1):   f_bond(c1.f_bond), f_inner_angle(c1.f_inner_angle), f_inner_dihedral(c1.f_inner_dihedral), f_outer_angle_m(c1.f_outer_angle_m), f_outer_angle_p(c1.f_outer_angle_p), f_outer_dihedral(c1.f_outer_dihedral),
                                                    r0(c1.r0), angle0(c1.angle0), outer_angle_m0(c1.outer_angle_m0), outer_angle_p0(c1.outer_angle_p0), inner_dih0(c1.inner_dih0), outer_dih0_a(c1.outer_dih0_a), outer_dih0_m(c1.outer_dih0_m), outer_dih0_p(c1.outer_dih0_p){
        
    }
};

struct EnergyConstants{
    const device_coord3d f_bond;
    const device_coord3d f_inner_angle;
    const device_coord3d f_inner_dihedral;

    const device_coord3d r0;
    const device_coord3d angle0;
    const device_coord3d inner_dih0;
    __device__ EnergyConstants(const device_coord3d f_bond, const device_coord3d f_inner_angle, const device_coord3d f_inner_dihedral, const device_coord3d r0, const device_coord3d angle0, const device_coord3d inner_dih0): f_bond(f_bond), f_inner_angle(f_inner_angle),
                            f_inner_dihedral(f_inner_dihedral),r0(r0), angle0(angle0), inner_dih0(inner_dih0) {}
};

struct ArcConstants{
    const device_real_t f_bond;
    const device_real_t f_inner_angle;
    const device_real_t f_inner_dihedral;
    const device_real_t f_outer_angle_m;
    const device_real_t f_outer_angle_p;
    const device_real_t f_outer_dihedral;

    const device_real_t r0;
    const device_real_t angle0;
    const device_real_t outer_angle_m0;
    const device_real_t outer_angle_p0;
    const device_real_t inner_dih0;
    const device_real_t outer_dih0;
    
    __device__ ArcConstants(const device_real_t f_bond, const device_real_t f_inner_angle, const device_real_t f_inner_dihedral, const device_real_t f_outer_angle_m, const device_real_t f_outer_angle_p, const device_real_t f_outer_dihedral,
                            const device_real_t r0, const device_real_t angle0, const device_real_t outer_angle_m0, const device_real_t outer_angle_p0, const device_real_t inner_dih0, const device_real_t outer_dih0): f_bond(f_bond), f_inner_angle(f_inner_angle),
                            f_inner_dihedral(f_inner_dihedral), f_outer_angle_m(f_outer_angle_m), f_outer_angle_p(f_outer_angle_p), f_outer_dihedral(f_outer_dihedral), r0(r0), angle0(angle0), outer_angle_m0(outer_angle_m0), outer_angle_p0(outer_angle_p0),
                            inner_dih0(inner_dih0), outer_dih0(outer_dih0) {}

};

__device__ __host__ uint8_t face_index(uint8_t f1, uint8_t f2, uint8_t f3){
    return f1*4 + f2*2 + f3;
}

template <typename T>
__device__ Constants<T> compute_constants(BookkeepingData& dat, device_node_t node_id){
    T r0, angle0, inner_dih0, outer_angle_m0, outer_angle_p0, outer_dih0_a, outer_dih0_m, outer_dih0_p;
    T f_bond, f_inner_angle, f_inner_dihedral, f_outer_angle_m, f_outer_angle_p, f_outer_dihedral ;
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
        uint8_t F1 = dat.face_right[node_id * 3 + j] - 5;
        uint8_t F3 = dat.face_right[node_id * 3 + (2 + j)%3] - 5;

        uint8_t f_r_neighbour = dat.face_right[dat.neighbours[node_id * 3 + j]*3 ] - 5;
        uint8_t f_l_neighbour = dat.face_right[dat.neighbours[node_id * 3 + j]*3 + 1 ] - 5;
        uint8_t f_m_neighbour = dat.face_right[dat.neighbours[node_id * 3 + j]*3 + 2] - 5;
        uint8_t F4 = f_r_neighbour + f_l_neighbour + f_m_neighbour - F1 - F3 ;


        uint8_t face_sum = dat.face_right[node_id * 3] - 5 + dat.face_right[node_id * 3 + 1] - 5 + dat.face_right[node_id * 3 + 2] - 5;
        uint8_t dihedral_face_sum = dat.face_right[dat.neighbours[node_id * 3 + j] * 3]-5 + dat.face_right[dat.neighbours[node_id * 3 + j] * 3 + 1]-5 +  dat.face_right[dat.neighbours[node_id * 3 + j] * 3 + 2]-5;

        //Load equillibirium distance, angles and dihedral angles from face information.
        d_set(r0,j,optimal_bond_lengths[ F3 + F1 ]);
        d_set(angle0,j,optimal_corner_cos_angles[ F1 ]);
        d_set(inner_dih0,j,optimal_dih_cos_angles[ face_index(dat.face_right[node_id * 3 + j] - 5, dat.face_right[node_id * 3 + (1+j)%3] - 5 , dat.face_right[node_id * 3 + (2+j)%3] - 5) ]);
        d_set(outer_angle_m0,j,optimal_corner_cos_angles[ F3 ]);
        d_set(outer_angle_p0,j,optimal_corner_cos_angles[ F1 ]);
        
        uint8_t dihedral_index_a = face_index(F3,F4,F1);
        uint8_t dihedral_index_m =  face_index(F4, F1, F3);
        uint8_t dihedral_index_p = face_index(F1,F3, F4);

        d_set(outer_dih0_a,j,optimal_dih_cos_angles[dihedral_index_a]  );
        d_set(outer_dih0_m,j,optimal_dih_cos_angles[dihedral_index_m]  );
        d_set(outer_dih0_p,j,optimal_dih_cos_angles[dihedral_index_p]  );

        //Load force constants from neighbouring face information.
        d_set(f_bond,j,bond_forces[ F3 + F1 ]);
        d_set(f_inner_angle,j,angle_forces[ F1 ]);
        d_set(f_inner_dihedral,j,dih_forces[ face_sum]);
        d_set(f_outer_angle_m,j,angle_forces[ F3 ]);
        d_set(f_outer_angle_p,j,angle_forces[ F1 ]);
        d_set(f_outer_dihedral,j,dih_forces[ dihedral_face_sum]);
    }
    return Constants<T>(f_bond,f_inner_angle,f_inner_dihedral, f_outer_angle_m, f_outer_angle_p, f_outer_dihedral, r0, angle0, outer_angle_m0, outer_angle_p0, inner_dih0, outer_dih0_a, outer_dih0_m, outer_dih0_p);
}

__device__ ArcConstants compute_arc_constants(BookkeepingData &dat, device_node_t node_id, uint8_t j){
    device_real_t r0 ; device_real_t angle0 ; device_real_t inner_dih0 ; device_real_t outer_angle_m0 ; device_real_t outer_angle_p0 ; device_real_t outer_dih0 ;
    device_real_t f_bond ; device_real_t f_inner_angle ; device_real_t f_inner_dihedral ; device_real_t f_outer_angle_m ; device_real_t f_outer_angle_p ; device_real_t f_outer_dihedral ;
    
    

    uint8_t f_r = dat.face_right[node_id * 3 + j] - 5;
    uint8_t f_l = dat.face_right[node_id * 3 + (2 + j)%3] - 5;

    uint8_t face_sum = dat.face_right[node_id * 3] - 5 + dat.face_right[node_id * 3 + 1] - 5 + dat.face_right[node_id * 3 + 2] - 5;
    uint8_t dihedral_face_sum = dat.face_right[dat.neighbours[node_id * 3 + j] * 3]-5 + dat.face_right[dat.neighbours[node_id * 3 + j] * 3 + 1]-5 +  dat.face_right[dat.neighbours[node_id * 3 + j] * 3 + 2]-5;

    //Load equillibirium distance, angles and dihedral angles from face information.
    r0 = optimal_bond_lengths[ f_l + f_r ];
    angle0 = optimal_corner_cos_angles[ f_r ];
    inner_dih0 = optimal_dih_cos_angles[ face_sum ];
    outer_angle_m0 = optimal_corner_cos_angles[ f_l ];
    outer_angle_p0 = optimal_corner_cos_angles[ f_r ];
    outer_dih0 = optimal_dih_cos_angles[ dihedral_face_sum ];

    //Load force constants from neighbouring face information.
    f_bond = bond_forces[ f_l + f_r ];
    f_inner_angle = angle_forces[ f_l ];
    f_inner_dihedral = dih_forces[ face_sum];
    f_outer_angle_m = angle_forces[ f_r ];
    f_outer_angle_p = angle_forces[ f_l ];
    f_outer_dihedral = dih_forces[ dihedral_face_sum];
        

    return ArcConstants(f_bond,f_inner_angle,f_inner_dihedral, f_outer_angle_m, f_outer_angle_p, f_outer_dihedral, r0, angle0, outer_angle_m0, outer_angle_p0, inner_dih0, outer_dih0);
}


//Reduction method for single block fullerenes.
/*
__device__ device_real_t reduction(device_real_t *sdata){

    cg::thread_block block = cg::this_thread_block();
    cg::sync(block);
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(block);
    sdata[threadIdx.x] = cg::reduce(tile32, sdata[threadIdx.x], cg::plus<device_real_t>());
    cg::sync(block);

    device_real_t beta = 0.0;
    if (block.thread_rank() == 0) {
        beta  = 0;
        for (uint16_t i = 0; i < block.size(); i += tile32.size()) {
            beta  += sdata[i];
        }
        sdata[0] = beta;
    }
    cg::sync(block);
    return sdata[0];
}*/


__device__ device_real_t reduction(device_real_t* sdata, const device_real_t data){
    sdata[threadIdx.x] = data;
    cg::thread_block block = cg::this_thread_block();
    cg::sync(block);
    
    if((Block_Size_Pow_2 > 512)){if (threadIdx.x < 512){sdata[threadIdx.x] += sdata[threadIdx.x + 512];} cg::sync(block);}
    if((Block_Size_Pow_2 > 256)){if (threadIdx.x < 256){sdata[threadIdx.x] += sdata[threadIdx.x + 256];} cg::sync(block);}
    if((Block_Size_Pow_2 > 128)){if (threadIdx.x < 128){sdata[threadIdx.x] += sdata[threadIdx.x + 128];} cg::sync(block);}
    if((Block_Size_Pow_2 > 64)){if (threadIdx.x < 64){sdata[threadIdx.x] += sdata[threadIdx.x + 64];} cg::sync(block);}
    if(threadIdx.x < 32){
    if((Block_Size_Pow_2 > 32)){if (threadIdx.x < 32){sdata[threadIdx.x] += sdata[threadIdx.x + 32];} __syncwarp();}
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(block);
    sdata[threadIdx.x] = cg::reduce(tile32, sdata[threadIdx.x], cg::plus<device_real_t>());
    }
    cg::sync(block);
    device_real_t sum = sdata[0];
    cg::sync(block);
    return sum;
}


__device__ device_real_t reduction_max(device_real_t* sdata, const device_real_t data){
    sdata[threadIdx.x] = data;
    cg::thread_block block = cg::this_thread_block();
    cg::sync(block);
    
    if((Block_Size_Pow_2 > 512)){if (threadIdx.x < 512){sdata[threadIdx.x] = max(sdata[threadIdx.x + 512],sdata[threadIdx.x]);} cg::sync(block);}
    if((Block_Size_Pow_2 > 256)){if (threadIdx.x < 256){sdata[threadIdx.x] = max(sdata[threadIdx.x + 256],sdata[threadIdx.x]);} cg::sync(block);}
    if((Block_Size_Pow_2 > 128)){if (threadIdx.x < 128){sdata[threadIdx.x] = max(sdata[threadIdx.x + 128],sdata[threadIdx.x]);} cg::sync(block);}
    if((Block_Size_Pow_2 > 64)){if (threadIdx.x < 64){sdata[threadIdx.x] = max(sdata[threadIdx.x + 64],sdata[threadIdx.x]);} cg::sync(block);}
    if(threadIdx.x < 32){
    if((Block_Size_Pow_2 > 32)){if (threadIdx.x < 32){sdata[threadIdx.x] = max(sdata[threadIdx.x + 32],sdata[threadIdx.x]);} __syncwarp();}
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(block);
    sdata[threadIdx.x] = cg::reduce(tile32, sdata[threadIdx.x], cg::greater<device_real_t>()); 
    }
    cg::sync(block);
    device_real_t max = sdata[0];
    cg::sync(block);
    return max;
}

__device__ half reduction(half *sdata){

    cg::thread_block block = cg::this_thread_block();
    cg::sync(block);
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(block);
    sdata[threadIdx.x] = cg::reduce(tile32, sdata[threadIdx.x], cg::plus<half>());
    cg::sync(block);

    half beta = 0.0;
    if (block.thread_rank() == 0) {
        beta  = 0;
        for (uint16_t i = 0; i < block.size(); i += tile32.size()) {
            beta  += sdata[i];
        }
        sdata[0] = beta;
    }
    cg::sync(block);
    return sdata[0];
}

//Multi purpose reduction algorithm (Small or Large fullerenes).
__device__ void reduction(device_real_t *sdata, device_real_t *gdata, const device_node_t N, const bool single_block_fullerenes){
    cg::thread_block block = cg::this_thread_block();

    cg::sync(block);
    if (((threadIdx.x + blockIdx.x * blockDim.x) >= N) && !single_block_fullerenes)
    {
        sdata[threadIdx.x] = 0;
    }
    cg::sync(block);
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(block);
    sdata[threadIdx.x] = cg::reduce(tile32, sdata[threadIdx.x], cg::plus<device_real_t>());
    cg::sync(block);
    
    device_real_t beta = 0.0;
    if (single_block_fullerenes)
    {
        if (block.thread_rank() == 0) {
            for (uint16_t i = 0; i < block.size(); i += tile32.size()) {
                beta  += sdata[i];
            }
            sdata[0] = beta;
        }
        cg::sync(block);
    }
    else 
    {   
        auto grid = cg::this_grid();
        if (block.thread_rank() == 0) 
        {
            for (uint16_t i = 0; i < block.size(); i += tile32.size()) 
            {
                beta  += sdata[i];
            }
            gdata[blockIdx.x] = beta;
        }
        cg::sync(grid);
        beta = 0.0;
        if (grid.thread_rank() == 0)
        {
            for (uint16_t i = 0; i < gridDim.x; i++) 
            {
                beta  += gdata[i];
            }
            gdata[0] = beta;
        }
        cg::sync(grid);
        if (block.thread_rank() == 0) {sdata[0] = gdata[0];}
        cg::sync(grid);
    }
}

template < class T >
size_t optimize_block_size(size_t N, cudaDeviceProp prop, T kernel){
    int maxActiveBlocks;
    size_t best_size = prop.warpSize;
    size_t min_waste = prop.maxThreadsPerMultiProcessor;
    for (size_t blocksize = prop.warpSize*2; blocksize < prop.maxThreadsPerBlock; blocksize +=prop.warpSize)
    {
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&maxActiveBlocks, kernel, (int)blocksize, (size_t)(sizeof(device_real_t)*blocksize));
        size_t wasted_threads = N % blocksize  +  blocksize*maxActiveBlocks % prop.maxThreadsPerMultiProcessor;

        if (wasted_threads < min_waste)
        {
            min_waste = wasted_threads;
            best_size = blocksize;
        }
    }
    return best_size;
}


__device__ device_real_t AverageBondLength(device_real_t* smem,const device_coord3d* X, const device_node_t* neighbours,const device_node_t node_id, const size_t N){
    device_real_t node_average_bond_length = 0.0;
    for (size_t i = 0; i < 3; i++)
    {
        node_average_bond_length += non_resciprocal_bond_length(X[node_id] - X[neighbours[i]]);
    }
    node_average_bond_length /= (device_real_t)3.0;

    
    return reduction(smem, node_average_bond_length)/(device_real_t)N;
}

__HD__ void print(const device_coord3d& ab){
    printf("[%.8e, %.8e, %.8e]\n",ab.x,ab.y,ab.z);
}
__device__ void print(const half4& ab){
    print_coord(ab);
}

__device__ void print(const half2& ab){
    printf("[%.16e, %.16e] \n", __half2float(ab.x), __half2float(ab.y));
}

__HD__ void print(device_real_t a){
    printf("[%.16e]\n", a);
}

__HD__ void print(bool b){
    printf("[%d]\n",int(b));
}

__HD__ void print(int a){
    printf("[%d]\n",a);
}

__device__ void print(const ushort3& a){
    printf("[%d, %d, %d]\n",a.x,a.y,a.z);
}

__device__ void print(const uchar3& a){
    printf("[%d, %d, %d]\n",a.x,a.y,a.z);
}

__device__ void print(const uint3& a){
    printf("[%d, %d, %d]\n",a.x,a.y,a.z);
}

template <typename T>
__device__ void print_single(T data){
    if (threadIdx.x + blockIdx.x == 0) {
        print(data);
    }
}

template <typename T>
__device__ void sequential_print(T* data){
    for (size_t i = 0; i < blockDim.x; i++)
    {
        if (threadIdx.x == i)
        {
            print(data[i]);
        }
        cg::sync(cg::this_thread_block());
    }
}

template <typename T>
__device__ void sequential_print(T data, size_t fullerene_id){
    if (blockIdx.x == fullerene_id)
    {
    for (size_t i = 0; i < blockDim.x; i++)
    {
        if (threadIdx.x == i)
        {
            print(data);
        }
        cg::sync(cg::this_thread_block());
    }
    }
}

template <typename T>
__host__ void print_array(T* data, size_t N, size_t fullerene_id){
    for (size_t i = 0; i < N; i++)
    {
        print(data[fullerene_id + i]);
    }
}

template <typename T>
__host__ void toBinary(std::string filename,T* data, size_t N, size_t fullerene_id){
    T* pointer =  data + N * fullerene_id;
    std::fstream myFile (filename, std::fstream::out | std::fstream::in | std::fstream::trunc | std::fstream::binary );

    myFile.write(reinterpret_cast<const char*>(pointer), sizeof(T)*N);
    if(!myFile)
      std::cout<<"error";
    myFile.close();
}

template <typename T>
__device__ void sequential_print(T* data, size_t fullerene_id){
    if (blockIdx.x == fullerene_id)
    {
    for (size_t i = 0; i < blockDim.x; i++)
    {
        if (threadIdx.x == i)
        {
            print(data[i]);
        }
        cg::sync(cg::this_thread_block());
    }
    }
}
/*
template <typename T>
__HD__ void swap(T& a, T& b){
    T temp = a;
    a = b;
    b = temp;
}*/

void printLastCudaError(std::string message = ""){
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess){
        std::cout << "\n" << message << " :\t";
        std::cout << cudaGetErrorString(error);
        printf("\n");
    }
}

__device__ void clear_cache(device_real_t* sdata, size_t N){
    cg::sync(cg::this_thread_block());
    for (size_t index = threadIdx.x; index < N; index+=blockDim.x)
    {
        sdata[index] = 0;
    }
    cg::sync(cg::this_thread_block());
}
