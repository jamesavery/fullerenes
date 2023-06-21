#include <cub/cub.cuh>
#include "fullerenes/gpu/cuda_io.hh"


namespace gpu_kernels{
  namespace cub_wrappers {

// Declare, allocate, and initialize device-accessible pointers for input samples and
// output histogram
template <typename T, typename ctrT> 
void histogram(size_t num_samples,    
	       const T* d_samples,
	       size_t num_levels,
	       T lower_level,
	       T upper_level,
	       ctrT*  d_histogram,
	       const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC
	       )
{
  // Determine temporary device storage requirements
  void*    d_temp_storage = NULL;
  size_t   temp_storage_bytes = 0;
  cub::DeviceHistogram::HistogramEven(d_temp_storage, temp_storage_bytes,
				      d_samples, d_histogram, num_levels, lower_level, upper_level, num_samples,ctx.stream);
  // Allocate temporary storage
  cudaMalloc(&d_temp_storage, temp_storage_bytes);
  // Compute histograms
  cub::DeviceHistogram::HistogramEven(d_temp_storage, temp_storage_bytes,
				      d_samples, d_histogram, num_levels, lower_level, upper_level, num_samples,ctx.stream);
  cudaFree(d_temp_storage);
  if(policy == LaunchPolicy::SYNC) ctx.wait();
}

  };
};
