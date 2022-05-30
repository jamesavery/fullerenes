#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "fullerenes/gpu/isomerspace_kernel.hh"
#include "fullerenes/gpu/isomerspace_X0.hh"

__device__
device_node_t multiple_source_shortest_paths(const IsomerBatch& G, device_node_t* distances){
    DEVICE_TYPEDEFS
    
    DeviceFullereneGraph FG = DeviceFullereneGraph(&G.cubic_neighbours[blockIdx.x*blockDim.x*3]);
    node_t outer_face[6];
    uint8_t Nface = FG.get_face_oriented(0, FG.cubic_neighbours[0],outer_face);
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
                node_t w = FG.cubic_neighbours[v*3 + i];
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
    DEVICE_TYPEDEFS

    node_t distance =  multiple_source_shortest_paths(G,reinterpret_cast<node_t*>(sdata));
    BLOCK_SYNC
    clear_cache(reinterpret_cast<real_t*>(sdata), Block_Size_Pow_2); 
    node_t d_max = reduction_max(sdata, distance);

    clear_cache(reinterpret_cast<real_t*>(sdata), Block_Size_Pow_2); 
    ordered_atomic_add(&reinterpret_cast<real_t*>(sdata)[distance],real_t(1.0)); 
    BLOCK_SYNC
    node_t num_of_same_dist = node_t(reinterpret_cast<real_t*>(sdata)[distance]); 
    BLOCK_SYNC
    clear_cache(reinterpret_cast<real_t*>(sdata), Block_Size_Pow_2);
    BLOCK_SYNC
    coord2d xys = reinterpret_cast<coord2d*>(G.xys)[blockIdx.x*blockDim.x + threadIdx.x]; BLOCK_SYNC
    ordered_atomic_add(&reinterpret_cast<real_t*>(sdata)[distance*2], xys.x); 
    ordered_atomic_add(&reinterpret_cast<real_t*>(sdata)[distance*2+1], xys.y); BLOCK_SYNC
    coord2d centroid = reinterpret_cast<coord2d*>(sdata)[distance] / num_of_same_dist; BLOCK_SYNC    
    coord2d xy = xys - centroid;
    real_t dtheta = real_t(M_PI)/real_t(d_max+1); 
    real_t phi = dtheta*(distance + real_t(0.5)); 
    real_t theta = atan2(xy.x, xy.y); 
    coord2d spherical_layout = {theta, phi};
    

    return spherical_layout;
}

__global__
void kernel_zero_order_geometry(IsomerBatch G, device_real_t scalerad){
    typedef device_real_t real_t;
    typedef device_coord2d coord2d;
    typedef device_coord3d coord3d;

    extern __shared__  device_real_t sdata[];
    clear_cache(sdata, Block_Size_Pow_2);
    if (G.statuses[blockIdx.x] == NOT_CONVERGED)
    {
    NodeGraph node_graph = NodeGraph(G, blockIdx.x); 
    coord2d angles = spherical_projection(G,reinterpret_cast<device_node_t*>(sdata));
    real_t theta = angles.x; real_t phi = angles.y;
    real_t x = cos(theta)*sin(phi), y = sin(theta)*sin(phi), z = cos(phi);
    coord3d coordinate = {x, y ,z};

    clear_cache(sdata, Block_Size_Pow_2);
    x = reduction(sdata, coordinate.x); y = reduction(sdata, coordinate.y); z = reduction(sdata,coordinate.z);
    coord3d cm = {x, y, z};
    cm /= blockDim.x;
    coordinate -= cm;

    real_t Ravg = real_t(0.0);
    clear_cache(sdata, Block_Size_Pow_2);
    real_t* base_pointer = sdata + Block_Size_Pow_2; 
    coord3d* X = reinterpret_cast<coord3d*>(base_pointer);
    X[threadIdx.x] = coordinate;
    BLOCK_SYNC
    real_t local_Ravg = real_t(0.0);
    for (uint8_t i = 0; i < 3; i++) local_Ravg += norm(X[threadIdx.x] - X[d_get(node_graph.cubic_neighbours,i)]);
    
    Ravg = reduction(sdata, local_Ravg);
    Ravg /= real_t(3*blockDim.x);
    coordinate *= scalerad*1.5/Ravg;
    reinterpret_cast<coord3d*>(G.X)[blockDim.x*blockIdx.x + threadIdx.x] = coordinate;
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
    synchronize();
    for (size_t i = 0; i < device_count; i++)
    {
        cudaSetDevice(i);
        void* kernelArgs[] = {(void*)&d_batch[i], (void*)&scalerad};
        safeCudaKernelCall((void*)kernel_zero_order_geometry, dim3(device_capacities[i], 1, 1), dim3(N, 1, 1), kernelArgs, shared_memory_bytes, main_stream[i]);
    }
    synchronize();
    
    auto end = std::chrono::system_clock::now();
    printLastCudaError("IsomerspaceX0 kernel launch failed: ");
}

IsomerspaceX0::IsomerspaceX0(const size_t N) : IsomerspaceKernel::IsomerspaceKernel(N, (void*)kernel_zero_order_geometry){
    this->shared_memory_bytes = sizeof(device_coord3d)*N + sizeof(device_real_t)*Block_Size_Pow_2;
    std::cout << "X0 Kernel Capacity: " << this->batch_capacity << "\n";
    printLastCudaError("X0 kernel class instantiation failed!");
}

IsomerspaceX0::~IsomerspaceX0(){}
