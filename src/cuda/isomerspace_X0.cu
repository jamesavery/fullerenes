#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "fullerenes/gpu/isomerspace_kernel.hh"
#include "auxiliary_cuda_functions.cu"
#include "forcefield_structs.cu"
#include "fullerenes/gpu/isomerspace_X0.hh"

typedef IsomerspaceKernel::device_real_t device_real_t;
typedef IsomerspaceKernel::device_node_t device_node_t;
#define BLOCK_SYNC cg::sync(cg::this_thread_block());

__device__
device_node_t multiple_source_shortest_paths(const IsomerBatch& G, device_node_t* distances){
    typedef device_node_t node_t;
    DeviceFullereneGraph FG = DeviceFullereneGraph(&G.neighbours[blockIdx.x*blockDim.x*3]);
    uint8_t Nface = FG.face_size(0,FG.neighbours[0]);
    node_t* outer_face =  FG.get_face_oriented(0, FG.neighbours[0]);
    distances[threadIdx.x] = node_t(NODE_MAX);    
    BLOCK_SYNC
    if (threadIdx.x < Nface)  distances[outer_face[threadIdx.x]] = 0;
    BLOCK_SYNC
    if (threadIdx.x == 0){
        CuDeque<node_t> queue = CuDeque<node_t>(distances + blockDim.x, blockDim.x);
        for (size_t i = 0; i < Nface; i++) queue.push_back(outer_face[i]);
        while (!queue.empty())
        {   
            node_t v = queue.pop_front();
            for (size_t i = 0; i < 3; i++)
            {   
                node_t w = FG.neighbours[v*3 + i];
                if(distances[w] == NODE_MAX) {
                distances[w] = distances[v]+1;
                queue.push_back(w);
                }
            }
        }
    }
    BLOCK_SYNC
    device_node_t distance = distances[threadIdx.x];
    BLOCK_SYNC
    return distance;
}


__device__
device_coord2d spherical_projection(const IsomerBatch& G, device_node_t* sdata){
    typedef device_node_t node_t;
    typedef device_real_t real_t;
    typedef device_coord2d coord2d;

    node_t distance =  multiple_source_shortest_paths(G,reinterpret_cast<node_t*>(sdata));
    BLOCK_SYNC
    clear_cache(reinterpret_cast<real_t*>(sdata), Block_Size_Pow_2); BLOCK_SYNC
    node_t d_max = reduction_max(sdata, distance); BLOCK_SYNC

    clear_cache(reinterpret_cast<real_t*>(sdata), Block_Size_Pow_2); BLOCK_SYNC
    atomicAdd_block(&reinterpret_cast<real_t*>(sdata)[distance],real_t(1.0)); BLOCK_SYNC
    node_t num_of_same_dist = node_t(reinterpret_cast<real_t*>(sdata)[distance]); 
    BLOCK_SYNC
    reinterpret_cast<coord2d*>(sdata)[distance] = {real_t(0.0),real_t(0.0)};
    BLOCK_SYNC
    coord2d xys = reinterpret_cast<coord2d*>(G.xys)[blockIdx.x*blockDim.x + threadIdx.x]; BLOCK_SYNC
    atomicAdd_block(&reinterpret_cast<real_t*>(sdata)[distance*2], xys.x); BLOCK_SYNC
    atomicAdd_block(&reinterpret_cast<real_t*>(sdata)[distance*2+1], xys.y); BLOCK_SYNC
    coord2d centroid = reinterpret_cast<coord2d*>(sdata)[distance] / num_of_same_dist; BLOCK_SYNC   
    coord2d xy = xys - centroid; BLOCK_SYNC
    real_t dtheta = real_t(M_PI)/real_t(d_max+1); BLOCK_SYNC
    real_t phi = dtheta*(distance + real_t(0.5)); BLOCK_SYNC
    real_t theta = atan2f(xy.x, xy.y); BLOCK_SYNC
    coord2d spherical_layout = {theta, phi}; BLOCK_SYNC
    /*
    */
    return spherical_layout;
}

__global__
void kernel_zero_order_geometry(IsomerBatch G, device_real_t scalerad){
    typedef device_real_t real_t;
    typedef device_coord2d coord2d;
    typedef device_coord3d coord3d;

    extern __shared__  device_real_t sdata[];
    clear_cache(sdata, Block_Size_Pow_2);BLOCK_SYNC
    if (G.statuses[blockIdx.x] == NOT_CONVERGED)
    {
    NodeGraph node_graph = NodeGraph(G); BLOCK_SYNC
    coord2d angles = spherical_projection(G,reinterpret_cast<device_node_t*>(sdata));BLOCK_SYNC
    real_t theta = angles.x; real_t phi = angles.y;BLOCK_SYNC
    real_t x = cosf(theta)*sinf(phi), y = sinf(theta)*sinf(phi), z = cosf(phi);BLOCK_SYNC
    coord3d coordinate = {x, y ,z};BLOCK_SYNC

    clear_cache(sdata, Block_Size_Pow_2);BLOCK_SYNC
    x = reduction(sdata, coordinate.x); y = reduction(sdata, coordinate.y); z = reduction(sdata,coordinate.z);BLOCK_SYNC
    coord3d cm = {x, y, z};
    cm /= blockDim.x;BLOCK_SYNC
    coordinate -= cm;BLOCK_SYNC

    real_t Ravg = real_t(0.0);BLOCK_SYNC
    clear_cache(sdata, Block_Size_Pow_2);BLOCK_SYNC
    real_t* base_pointer = sdata + Block_Size_Pow_2; BLOCK_SYNC
    coord3d* X = reinterpret_cast<coord3d*>(base_pointer);BLOCK_SYNC
    X[threadIdx.x] = coordinate;
    BLOCK_SYNC
    real_t local_Ravg = real_t(0.0);BLOCK_SYNC
    for (uint8_t i = 0; i < 3; i++) local_Ravg += norm(X[threadIdx.x] - X[d_get(node_graph.neighbours,i)]);
    BLOCK_SYNC
    Ravg = reduction(sdata, local_Ravg);
    Ravg /= real_t(3*blockDim.x);

    BLOCK_SYNC
    coordinate *= scalerad*1.5/Ravg;
    reinterpret_cast<coord3d*>(G.X)[blockDim.x*blockIdx.x + threadIdx.x] = coordinate;BLOCK_SYNC
    BLOCK_SYNC
    G.statuses[blockIdx.x] = CONVERGED;
    }
    
}
void IsomerspaceX0::check_batch(){
    batch_size = 0;
    for (size_t i = 0; i < device_count; i++)
    {
        cudaSetDevice(i);
        int num_of_not_converged_isomers = 0;
        h_batch[i] <<= d_batch[i];
        
        IsomerStatus* statuses = h_batch[i].statuses;
        for (int j = 0; j < device_capacities[i]; j++)
        {   
            num_of_not_converged_isomers += (int)(statuses[j] == NOT_CONVERGED);
            if (statuses[j] != NOT_CONVERGED)
            {   
                index_queue[i].push(j);
            }
        }
        batch_sizes[i] = num_of_not_converged_isomers;
        batch_size += num_of_not_converged_isomers;
    }
}

void IsomerspaceX0::zero_order_geometry(device_real_t scalerad){
    printLastCudaError("Memcpy Failed! \n");
    auto start = std::chrono::system_clock::now();
    //for (size_t i = 0; i < device_count; i++) {cudaSetDevice(i); d_batch[i] <<= h_batch[i];}
    cudaDeviceSynchronize();
    for (size_t i = 0; i < device_count; i++)
    {
        cudaSetDevice(i);
        void* kernelArgs[] = {(void*)&d_batch[i], (void*)&scalerad};
        safeCudaKernelCall((void*)kernel_zero_order_geometry, dim3(device_capacities[i], 1, 1), dim3(N, 1, 1), kernelArgs, shared_memory_bytes);
    }
    //for (size_t i = 0; i < device_count; i++) {cudaSetDevice(i); h_batch[i] <<= d_batch[i];}
        
    cudaDeviceSynchronize();
    
    auto end = std::chrono::system_clock::now();
    printLastCudaError("IsomerspaceX0 kernel launch failed: ");
}

IsomerspaceX0::IsomerspaceX0(const size_t N) : IsomerspaceKernel::IsomerspaceKernel(N, (void*)kernel_zero_order_geometry){
    this->shared_memory_bytes = sizeof(device_coord3d)*N + sizeof(device_real_t)*Block_Size_Pow_2;
    std::cout << "X0 Kernel Capacity: " << this->batch_capacity << "\n";
    printLastCudaError("X0 kernel class instantiation failed!");
}

IsomerspaceX0::~IsomerspaceX0(){}
