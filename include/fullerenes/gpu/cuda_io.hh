#include <cuda_runtime.h>
#include "isomerspace_kernel.hh"

namespace cuda_io{
    cudaError_t output_to_queue(std::queue<std::pair<Polyhedron, size_t>>& queue, IsomerBatch& batch, const bool copy_2d_layout = true);
    cudaError_t copy(IsomerBatch& target, const IsomerBatch& source);
    cudaError_t resize(IsomerBatch& batch, size_t new_capacity, const cudaStream_t stream = NULL);
}
