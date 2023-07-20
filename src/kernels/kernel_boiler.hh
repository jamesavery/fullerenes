#if   DEVICE==CPU_BATCH
# define LOCAL_DEFINE(name,T,N) T name##_array[N]
# define LOCAL_USE(ix,name)     auto &name = name##_array[ix]
# define LOOP(ix,N)             for(size_t ix=0;ix<N;ix++)
# define KERNEL_SIGNATURE(T,U)  template <typename T, BufferType U> 

# define BLOCK_SYNC
# define GRID_SYNC       

#elif DEVICE==CPU_ISOMER
# define LOCAL_DEFINE(name,T,N) T name##_array[N]
# define LOCAL_USE(ix,name)     auto &name = name##_array[ix]
# define LOOP(ix,N)   	    \
  #pragma omp parallel for  \
  for(size_t ix=0;ix<N;i++)
# define KERNEL_SIGNATURE(T,U)  template <typename T, BufferType U> 

# define BLOCK_SYNC
# define GRID_SYNC       

#elif DEVICE==GPU_BATCH
# define LOCAL_DEFINE(name,T,N) T name
# define LOCAL_USE(ix,name)     
# define LOOP(ix,N)             BLOCK_SYNC; size_t ix = threadIdx.x; if(ix<N)
# define KERNEL_SIGNATURE(T,U)  template <typename T, BufferType U> __global__

# define BLOCK_SYNC __syncthreads()
# define GRID_SYNC cg::sync(cg::this_grid())

#elif DEVICE==GPU_ISOMER
# define LOCAL_DEFINE(name,T,N) T name
# define LOCAL_USE(ix,name)     
# define LOOP(ix,N)             \
  GRID_SYNC; size_t ix = blockDim.x*blockIdx.x+threadIdx.x; if(ix<N)
# define KERNEL_SIGNATURE(T,U)  template <typename T, BufferType U> __global__

# define BLOCK_SYNC __syncthreads()
# define GRID_SYNC cg::sync(cg::this_grid())

#endif

