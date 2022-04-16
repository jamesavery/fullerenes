#include "fullerenes/gpu/reductions.cuh"

__device__ device_node_t max(const device_node_t a, const device_node_t b){
    if (a > b){
        return a;
    }else {
        return b;
    }
}

#if REDUCTION_METHOD==0
    __device__ device_real_t reduction(device_real_t* sdata, const device_real_t data){
        sdata[threadIdx.x] = data;
        BLOCK_SYNC
        if((Block_Size_Pow_2 > 512)){if (threadIdx.x < 512){sdata[threadIdx.x] += sdata[threadIdx.x + 512];} BLOCK_SYNC}
        if((Block_Size_Pow_2 > 256)){if (threadIdx.x < 256){sdata[threadIdx.x] += sdata[threadIdx.x + 256];} BLOCK_SYNC}
        if((Block_Size_Pow_2 > 128)){if (threadIdx.x < 128){sdata[threadIdx.x] += sdata[threadIdx.x + 128];} BLOCK_SYNC}
        if((Block_Size_Pow_2 > 64)){if (threadIdx.x < 64){sdata[threadIdx.x] += sdata[threadIdx.x + 64];} BLOCK_SYNC}
        if(threadIdx.x < 32){
        if((Block_Size_Pow_2 > 32)){if (threadIdx.x < 32){sdata[threadIdx.x] += sdata[threadIdx.x + 32];} __syncwarp();}
        cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
        sdata[threadIdx.x] = cg::reduce(tile32, sdata[threadIdx.x], cg::plus<device_real_t>());
        }
        BLOCK_SYNC
        device_real_t sum = sdata[0];
        BLOCK_SYNC
        return sum;
    }
#elif REDUCTION_METHOD==1
    __device__ device_real_t reduction(device_real_t* sdata, const device_real_t data){
        sdata[threadIdx.x] = data;
        BLOCK_SYNC
        
        if (threadIdx.x < 512){sdata[threadIdx.x] += sdata[threadIdx.x + 512];} BLOCK_SYNC
        if (threadIdx.x < 256){sdata[threadIdx.x] += sdata[threadIdx.x + 256];} BLOCK_SYNC
        if (threadIdx.x < 128){sdata[threadIdx.x] += sdata[threadIdx.x + 128];} BLOCK_SYNC
        if (threadIdx.x < 64){sdata[threadIdx.x] += sdata[threadIdx.x + 64];} BLOCK_SYNC
        if(threadIdx.x < 32){
        if (threadIdx.x < 32){sdata[threadIdx.x] += sdata[threadIdx.x + 32];} __syncwarp();
        cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
        sdata[threadIdx.x] = cg::reduce(tile32, sdata[threadIdx.x], cg::plus<device_real_t>());
        }
        BLOCK_SYNC
        device_real_t sum = sdata[0];
        BLOCK_SYNC
        return sum;
    }
#elif REDUCTION_METHOD==2
    __device__ device_real_t reduction(device_real_t *sdata, const device_real_t data){
        sdata[threadIdx.x] = data;
        cg::thread_block block = cg::this_thread_block();
        BLOCK_SYNC
        cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(block);
        sdata[threadIdx.x] = cg::reduce(tile32, sdata[threadIdx.x], cg::plus<device_real_t>());
        BLOCK_SYNC

        device_real_t beta = 0.0;
        if (block.thread_rank() == 0) {
            beta  = 0;
            for (uint16_t i = 0; i < block.size(); i += tile32.size()) {
                beta  += sdata[i];
            }
            sdata[0] = beta;
        }
        BLOCK_SYNC
        device_real_t sum = sdata[0];
        BLOCK_SYNC
        return sum;
    }
#endif

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

__device__ device_node_t reduction_max(device_node_t* sdata, const device_node_t data){
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
    sdata[threadIdx.x] = cg::reduce(tile32, sdata[threadIdx.x], cg::greater<device_node_t>()); 
    }
    cg::sync(block);
    device_node_t max = sdata[0];
    cg::sync(block);
    return max;
}

__device__ device_real_t global_reduction(device_real_t *sdata, device_real_t *gdata, device_real_t data, bool mask){
    GRID_SYNC
    if(!mask){data = (device_real_t)0.0;}
    device_real_t block_sum    = reduction(sdata,data);
    if(threadIdx.x == 0){gdata[blockIdx.x] = block_sum;}
    GRID_SYNC

        if (gridDim.x > 1024 && threadIdx.x == 0 && ((blockIdx.x + 1024) < gridDim.x))   {if (blockIdx.x < 1024) {gdata[blockIdx.x]  += gdata[blockIdx.x + 1024];}} GRID_SYNC
        if (gridDim.x > 512 && threadIdx.x == 0 && ((blockIdx.x + 512) < gridDim.x))    {if (blockIdx.x < 512)  {gdata[blockIdx.x]  += gdata[blockIdx.x + 512];}} GRID_SYNC
        if (gridDim.x > 256 && threadIdx.x == 0 && ((blockIdx.x + 256) < gridDim.x))    {if (blockIdx.x < 256)  {gdata[blockIdx.x]  += gdata[blockIdx.x + 256];}} GRID_SYNC
        if (gridDim.x > 128 && threadIdx.x == 0 && ((blockIdx.x + 128) < gridDim.x))    {if (blockIdx.x < 128)  {gdata[blockIdx.x]  += gdata[blockIdx.x + 128];}} GRID_SYNC
        if (gridDim.x > 64 && threadIdx.x == 0 && ((blockIdx.x + 64) < gridDim.x))     {if (blockIdx.x < 64)   {gdata[blockIdx.x]  += gdata[blockIdx.x + 64];}} GRID_SYNC
        if (gridDim.x > 32 && threadIdx.x == 0 && ((blockIdx.x + 32) < gridDim.x))     {if (blockIdx.x < 32)   {gdata[blockIdx.x]  += gdata[blockIdx.x + 32];}} GRID_SYNC
        if (threadIdx.x < 32 && blockIdx.x == 0)
        {
            cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
            gdata[threadIdx.x] = cg::reduce(tile32, gdata[threadIdx.x], cg::plus<device_real_t>()); 
        }
    GRID_SYNC
    device_real_t sum = (device_real_t)0.0;
    sum = gdata[0];
    GRID_SYNC
    return sum;
}
