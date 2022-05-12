#include <cuda_runtime.h>
#include "isomerspace_kernel.hh"
#include "cuda_execution.hh"

namespace cuda_io{
    cudaError_t output_to_queue(std::queue<std::tuple<Polyhedron, size_t, IsomerStatus>>& queue, IsomerBatch& batch, const bool copy_2d_layout = true);
    cudaError_t copy(IsomerBatch& target, const IsomerBatch& source, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC, const std::pair<int,int>& lhs_range = {-1,-1}, const std::pair<int,int>& rhs_range = {-1,-1});
    cudaError_t resize(IsomerBatch& batch, size_t new_capacity, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC, int front = -1, int back = -1);
    cudaError_t reset_convergence_statuses(IsomerBatch& batch, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    void is_close(const IsomerBatch& a, const IsomerBatch& b, device_real_t tol = 0.0, bool verbose = false);
    int num_converged(const IsomerBatch& batch);
    int num_failed(const IsomerBatch& batch);
    void sort(IsomerBatch& batch);
}
