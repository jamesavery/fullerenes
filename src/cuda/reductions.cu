#include "fullerenes/gpu/reductions.cuh"
#include "print_functions.cu"
#define NUM_BANKS 32 
#define LOG_NUM_BANKS 5
#define CONFLICT_FREE_OFFSET(n) ((n) >> NUM_BANKS + (n) >> (2 * LOG_NUM_BANKS)) 
namespace cg = cooperative_groups;

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
#elif REDUCTION_METHOD==3
    __device__ device_real_t reduction(device_real_t* sdata, const device_real_t data = (device_real_t)0){
        auto num_warps = (blockDim.x >> 5) + 1;
        cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
        auto warpid = threadIdx.x >> 5;
        auto laneid = threadIdx.x & 31;
        device_real_t temp = cg::reduce(tile32, data, cg::plus<device_real_t>());
        if (num_warps > 1){
        //sdata[warpid + blockDim.x] = temp;
        BLOCK_SYNC
        if (laneid == 0) sdata[warpid] = temp;
        BLOCK_SYNC
        if (warpid == 0) {
            temp = cg::reduce(tile32, sdata[laneid], cg::plus<device_real_t>());
        }
        }
        if (threadIdx.x == 0) sdata[0] = temp;
        BLOCK_SYNC
        auto result = sdata[0];
        BLOCK_SYNC
        return result;
    }
#elif REDUCTION_METHOD==4
    __device__ device_real_t reduction(device_real_t* sdata, const device_real_t data = (device_real_t)0){
        typedef cub::BlockReduce<device_real_t, 256> BlockReduce;
        auto result = BlockReduce(reinterpret_cast<BlockReduce::TempStorage*>(sdata)).Sum(data)
        if (threadIdx.x == 0) sdata[0] = result;
        BLOCK_SYNC
        result = sdata[0];
        BLOCK_SYNC
        return result;
    }

#endif

//This code finds the maximum value of a set of data using a reduction algorithm.
//The function takes as input a pointer to a shared memory array (sdata) and threadIdx.x's data (data).
//The function returns the maximum value of the data set.
__device__ device_real_t reduction_max(device_real_t* sdata, const device_real_t data){
     auto num_warps = (blockDim.x >> 5) + 1;
        cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
        auto warpid = threadIdx.x >> 5;
        auto laneid = threadIdx.x & 31;
        device_real_t temp = cg::reduce(tile32, data, cg::greater<device_real_t>());
        if (num_warps > 1){
        //sdata[warpid + blockDim.x] = temp;
        BLOCK_SYNC
        if (laneid == 0) sdata[warpid] = temp;
        BLOCK_SYNC
        if (warpid == 0) {
            temp = cg::reduce(tile32, sdata[laneid], cg::greater<device_real_t>());
        }
        }
        if (threadIdx.x == 0) sdata[0] = temp;
        BLOCK_SYNC
        auto result = sdata[0];
        BLOCK_SYNC
        return result;
}

__device__ device_node_t reduction_max(device_node_t* sdata, const device_node_t data){
    sdata[threadIdx.x] = data;
    cg::thread_block block = cg::this_thread_block();
    cg::sync(block);
    
    if((Block_Size_Pow_2 > 512)){if (threadIdx.x < 512){sdata[threadIdx.x] = d_max(sdata[threadIdx.x + 512],sdata[threadIdx.x]);} cg::sync(block);}
    if((Block_Size_Pow_2 > 256)){if (threadIdx.x < 256){sdata[threadIdx.x] = d_max(sdata[threadIdx.x + 256],sdata[threadIdx.x]);} cg::sync(block);}
    if((Block_Size_Pow_2 > 128)){if (threadIdx.x < 128){sdata[threadIdx.x] = d_max(sdata[threadIdx.x + 128],sdata[threadIdx.x]);} cg::sync(block);}
    if((Block_Size_Pow_2 > 64)){if (threadIdx.x < 64){sdata[threadIdx.x] = d_max(sdata[threadIdx.x + 64],sdata[threadIdx.x]);} cg::sync(block);}
    if(threadIdx.x < 32){
    if((Block_Size_Pow_2 > 32)){if (threadIdx.x < 32){sdata[threadIdx.x] = d_max(sdata[threadIdx.x + 32],sdata[threadIdx.x]);} __syncwarp();}
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
template <typename T>
__inline__ __device__ T warpReduceSum(T val) {
#pragma unroll
  for (int offset = warpSize/2; offset > 0; offset >>= 1) 
    val += __shfl_down_sync(-1, val, offset);
  return val;
}

template <typename T>
__device__ void reduction_V0(T* sdata, const T data){
    sdata[threadIdx.x] = data;
    BLOCK_SYNC
    if (threadIdx.x){
        T result = (device_real_t)0.0;
        for(int i = 0; i < blockDim.x; ++i){
            result += sdata[i];
        }
        sdata[0] = result;
    }
    BLOCK_SYNC
}

template <typename T>
__device__ void reduction_V1(T* sdata, const T data){
    sdata[threadIdx.x] = data;
    BLOCK_SYNC
    for (unsigned int s=1; s < blockDim.x; s *= 2) { 
        int index = 2 * s * threadIdx.x;
        if (index < blockDim.x) { sdata[index] += sdata[index + s];}
        BLOCK_SYNC
    }
}

template <typename T>
__device__ void reduction_V2(T* sdata, const T data){
    sdata[threadIdx.x] = data;
    BLOCK_SYNC
    if(blockDim.x > 512){if (threadIdx.x < 512){sdata[threadIdx.x] += sdata[threadIdx.x + 512];} BLOCK_SYNC}
    if(blockDim.x > 256){if (threadIdx.x < 256){sdata[threadIdx.x] += sdata[threadIdx.x + 256];} BLOCK_SYNC}
    if(blockDim.x > 128){if (threadIdx.x < 128){sdata[threadIdx.x] += sdata[threadIdx.x + 128];} BLOCK_SYNC}
    if(blockDim.x > 64){if (threadIdx.x < 64){sdata[threadIdx.x] += sdata[threadIdx.x + 64];} BLOCK_SYNC}
    if(threadIdx.x < 32){
    if((Block_Size_Pow_2 > 32)){if (threadIdx.x < 32){sdata[threadIdx.x] += sdata[threadIdx.x + 32];} __syncwarp();}
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
    sdata[threadIdx.x] = cg::reduce(tile32, sdata[threadIdx.x], cg::plus<T>());
    }
    BLOCK_SYNC
}

template <typename T>
__device__ void reduction_V3(T* sdata, const T data = (T)0){
    auto num_warps = (blockDim.x >> 5) + 1;
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
    auto warpid = threadIdx.x >> 5;
    auto laneid = threadIdx.x & 31;
    T temp = cg::reduce(tile32, data, cg::plus<T>());
    if (num_warps > 1){
    if (laneid == 0) sdata[warpid + blockDim.x] = temp;
    BLOCK_SYNC
    if (warpid == 0) {
        temp = cg::reduce(cg::tiled_partition<32>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
    }

    }
    if (threadIdx.x == 0) sdata[0] = temp;
    BLOCK_SYNC
}


template <typename T>
__device__ void reduction_V4(T* sdata, const T data, unsigned int n_warps, int warp_switch){
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
    auto warpid = threadIdx.x >> 5;
    auto laneid = threadIdx.x & 31;
    T temp = cg::reduce(tile32, data, cg::plus<T>());
    if (n_warps > 1){
    if (laneid == 0) sdata[warpid + blockDim.x] = temp;
    BLOCK_SYNC
    if (warpid == 0) {
        switch (warp_switch)
        {
        case 1:
            temp = cg::reduce(cg::tiled_partition<2>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
            break;
        case 2:
            temp = cg::reduce(cg::tiled_partition<4>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
            break;
        case 3:
            temp = cg::reduce(cg::tiled_partition<8>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
            break;
        case 4:
            temp = cg::reduce(cg::tiled_partition<16>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
            break;
        case 5:
            temp = cg::reduce(cg::tiled_partition<32>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
            break;
        
        default:
            break;
        }
        
    }

    }
    if (threadIdx.x == 0) sdata[0] = temp;
    BLOCK_SYNC
}

template <typename T>
__device__ void reduction_V5(T* sdata, const T data = (T)0, T* result = nullptr){
    *result = (T)0;
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
    auto laneid = threadIdx.x & 31;
    T temp = cg::reduce(tile32, data, cg::plus<T>());
    if (laneid == 0) atomicAdd_block((float*)result, (float)temp);
    BLOCK_SYNC
}

template <typename T>
__device__ void reduction_V6(T* sdata, const T data = (T)0, T* result = nullptr){
    *result = (T)0;
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
    auto warpid = threadIdx.x >> 5;
    T temp = cg::reduce(tile32, data, cg::plus<T>());
    sdata[warpid + blockDim.x] = temp;
    BLOCK_SYNC
    if (warpid == 0) *result = cg::reduce(tile32, sdata[threadIdx.x + blockDim.x], cg::plus<T>());
    BLOCK_SYNC
}

template <typename T>
__device__ void reduction_V7(T* sdata, const T data = (T)0, T* result = nullptr, unsigned int n_warps = 1, int warp_switch = 0){
    *result = (T)0;
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());
    auto warpid = threadIdx.x >> 5;
    auto laneid = threadIdx.x & 31;
    T temp = cg::reduce(tile32, data, cg::plus<T>());
    sdata[warpid + blockDim.x] = temp;
    BLOCK_SYNC
    if (warpid == 0) {
        switch (warp_switch)
        {
        case 1:
            *result = cg::reduce(cg::tiled_partition<2>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
            break;
        case 2:
            *result = cg::reduce(cg::tiled_partition<4>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
            break;
        case 3:
            *result = cg::reduce(cg::tiled_partition<8>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
            break;
        case 4:
            *result = cg::reduce(cg::tiled_partition<16>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
            break;
        case 5:
            *result = cg::reduce(cg::tiled_partition<32>(cg::this_thread_block()), sdata[laneid + blockDim.x], cg::plus<T>());
            break;
        
        default:
            break;
        }
        
    }
    BLOCK_SYNC
}


__device__ void exclusive_scan(device_node_t* sdata,const device_node_t data, const int size = blockDim.x){
    sdata[threadIdx.x] = data;
    BLOCK_SYNC
    if (threadIdx.x == 0){
        auto temp = 0;
        auto temp2 = 0;
        for (size_t i = 1; i < size; i++)
        {   
            temp2 = sdata[i];
            sdata[i] = sdata[i - 1] + temp;
            temp = temp2;
        }
        sdata[0] = 0;
    }
    BLOCK_SYNC
}

template <typename T>
__device__ void prescan(T *sdata, const T data, int n)
{
  int thid = threadIdx.x;
  int offset = 1;
  if(threadIdx.x < n) sdata[thid] = data;
  #pragma unroll
  for (int d = n>>1; d > 0; d >>= 1)                    // build sum in place up the tree
  {
    __syncthreads();
    if (thid < d)
    {
      int ai = offset*(2*thid+1)-1;
      int bi = offset*(2*thid+2)-1;
      sdata[bi] += sdata[ai];
    }
    offset *= 2;
  }
     if (thid==0) { sdata[n - 1] = 0;}
  
  #pragma unroll
  for (int d = 1; d < n; d *= 2) // traverse down tree & build scan
    {
      offset >>= 1;
      __syncthreads();
      if (thid < d)
      {
         int ai = offset*(2*thid+1)-1;
         int bi = offset*(2*thid+2)-1;
         T t = sdata[ai];
         sdata[ai] = sdata[bi];
         sdata[bi] += t;
      }
    }
  __syncthreads();
}

template <typename T>
__device__ void ex_scan(T* sdata, const T data, int n){
    auto warpid = threadIdx.x >> 5;
    auto lane = threadIdx.x & 31;

    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());

    auto result = cg::inclusive_scan(tile32, data);

    if (lane == 31){
        sdata[n+1 + warpid] = result;
    }

    __syncthreads();
    if (warpid == 0){
        auto val = cg::inclusive_scan(tile32, sdata[n+1 + lane]);
        sdata[n+1 + lane] = val;
    }
    __syncthreads();
    if (warpid == 0)
    {
        sdata[threadIdx.x + 1] = result;
    } else{
        if (threadIdx.x < n) {
        sdata[threadIdx.x + 1] =  sdata[n+1 + warpid-1] + result;}
        
    }
    if (threadIdx.x == 0){
        sdata[0] = (device_node_t)0;
    }
    __syncthreads();
}

template <typename T>
__device__ void in_scan(T* sdata, const T data, int n){
    auto warpid = threadIdx.x >> 5;
    auto lane = threadIdx.x & 31;

    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cg::this_thread_block());

    auto result = cg::inclusive_scan(tile32, data);

    if (lane == 31){
        sdata[n+1 + warpid] = result;
    }

    __syncthreads();
    if (warpid == 0){
        auto val = cg::inclusive_scan(tile32, sdata[n+1 + lane]);
        sdata[n+1 + lane] = val;
    }
    __syncthreads();
    if (warpid == 0)
    {
        sdata[threadIdx.x] = result;
    } else{
        if (threadIdx.x < n) {
        sdata[threadIdx.x] =  sdata[n+1 + warpid-1] + result;}
        
    }
    __syncthreads();
}


//WARNING: this function is intended to solve a specific problem efficiently using block scans.
//It does not work if the number of blocks is larger than blocksize^2.
template <typename T>
__device__ __forceinline__ void grid_ex_scan(T* g_out, T* sdata, T data, int n){
    auto gid = blockDim.x * blockIdx.x + threadIdx.x;
    auto n_blocks_needed = (n + blockDim.x - 1) / blockDim.x;  //Fast ceiling integer division.
    g_out[blockIdx.x] = data;
    T result = 0;
    GRID_SYNC
    assert(n_blocks_needed <= gridDim.x);
    if (gid < n){
        sdata[threadIdx.x] = g_out[gid];
    } else{ sdata[threadIdx.x] = 0;}
    if(blockIdx.x < n_blocks_needed) {
        in_scan<T>(sdata, sdata[threadIdx.x], blockDim.x); //First compute individual block inclusive scans.
        result = sdata[threadIdx.x]; //Store this result for later
    }
    GRID_SYNC
    if (threadIdx.x == 0 && blockIdx.x < n_blocks_needed) g_out[blockIdx.x] = sdata[blockDim.x - 1]; 
    GRID_SYNC
    if (blockIdx.x == 0){
        in_scan<T>(sdata, g_out[threadIdx.x],n_blocks_needed); //Now reduce the ends of each scan in a second inclusive scan.
        g_out[threadIdx.x + 1] = sdata[threadIdx.x];  
    }

    GRID_SYNC
    T pre_sum = 0;
    if(blockIdx.x < n_blocks_needed) pre_sum = g_out[blockIdx.x];
    GRID_SYNC
    if(blockIdx.x == 0){
        if (threadIdx.x==0) g_out[0] = 0;
        g_out[gid + 1] = result; //Adding the initial scan results to the scan of the block-ends gives us our final scanned array.
    }
    else if(blockIdx.x < n_blocks_needed){
        g_out[gid + 1] = result + pre_sum; //Adding the initial scan results to the scan of the block-ends gives us our final scanned array.
    }
    GRID_SYNC
}



template <typename T>
__device__ void prescan_noconflict(T *sdata, const T data, int n)
{
    auto thid = threadIdx.x;
    int ai = thid;
    int bi = thid + (n/2);
    int bankOffsetA = CONFLICT_FREE_OFFSET(ai);
    int bankOffsetB = CONFLICT_FREE_OFFSET(bi);
    auto temp_a = sdata[ai];
    auto temp_b = sdata[bi];

    sdata[ai + bankOffsetA] = temp_a;
    sdata[bi + bankOffsetB] = temp_b;
    int offset = 1;
    if(thid < n) sdata[thid] = data;

    for (int d = n>>1; d > 0; d >>= 1)                    // build sum in place up the tree
    {
        __syncthreads();
        if (thid < d)
        {
        int ai = offset*(2*thid+1)-1;
        int bi = offset*(2*thid+2)-1;
        ai += CONFLICT_FREE_OFFSET(ai);
        bi += CONFLICT_FREE_OFFSET(bi);
        sdata[bi] += sdata[ai];
        }
        offset *= 2;
    }
    if (thid==0) { sdata[n - 1 + CONFLICT_FREE_OFFSET(n - 1)] = 0;}
     for (int d = 1; d < n; d *= 2) // traverse down tree & build scan
        {
        offset >>= 1;
        __syncthreads();
        if (thid < d)
        {
            int ai = offset*(2*thid+1)-1;
            int bi = offset*(2*thid+2)-1;
            ai += CONFLICT_FREE_OFFSET(ai);
            bi += CONFLICT_FREE_OFFSET(bi);
            T t = sdata[ai];
            sdata[ai] = sdata[bi];
            sdata[bi] += t;
        }
    }

    __syncthreads();
    temp_a = sdata[ai + bankOffsetA];
    temp_b = sdata[bi + bankOffsetB];
    __syncthreads();
    sdata[ai] = temp_a;
    sdata[bi] = temp_b;
    __syncthreads();


}
