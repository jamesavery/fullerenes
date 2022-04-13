#include "fullerenes/gpu/batch_queue.hh"

namespace batch_queue_launch_params{
    bool refill_flag = false;
    size_t block_size;
    int device = 0;
    dim3 block_dims = dim3(1, 1, 1);
    dim3 grid_dims  = dim3(1, 1, 1);
}
__global__ void __refill_batch(IsomerBatch batch, BatchQueue& queue){

}



cudaError_t BatchQueue::refill_batch(IsomerBatch& batch, const cudaStream_t stream){
    using namespace batch_queue_launch_params;
    int current_device; cudaGetDevice(&current_device);
    if (!refill_flag || block_size != batch.n_atoms || device != current_device) 
    {
        refill_flag = true;
        int num_blocks = 0;
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&num_blocks,(void*)__refill_batch, batch.n_atoms, 0);
        block_dims.x = batch.n_atoms;
        grid_dims.x = num_blocks;
        device = current_device;
    }

    void* kargs[] = {(void*)&batch, (void*)&d_props};
    return cudaLaunchCooperativeKernel((void*)__refill_batch, grid_dims, block_dims, kargs, 0, stream);
}

cudaError_t BatchQueue::insert(const Polyhedron& in){
    if (h_props.capacity == h_props.size)
    {
        //resize(h_props.batch, h_props.capacity*2);
        //resize(d_props.batch, d_props.capacity*2);
    }
    for(int i = 0; i < in.neighbours.size(); i++){
        h_props.batch.neighbours[i];
    }
    
}

cudaError_t BatchQueue::insert(const IsomerBatch& input_batch){

}