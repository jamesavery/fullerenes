#ifndef CUDA_IO_H
#define CUDA_IO_H
#include <cuda_runtime.h>
#include "launch_ctx.hh"
#include <chrono>
#include <queue>
#include "fullerenes/gpu/isomer_batch.hh"

namespace cuda_io{
    template <Device T>
    cudaError_t output_to_queue(std::queue<std::tuple<Polyhedron, size_t, IsomerStatus>>& queue, IsomerBatch<T>& batch, const bool copy_2d_layout = true);
    
    template <Device T, Device U>
    cudaError_t copy(IsomerBatch<T>& target, const IsomerBatch<U>& source, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC, const std::pair<int,int>& lhs_range = {-1,-1}, const std::pair<int,int>& rhs_range = {-1,-1});
    
    template <Device T>
    cudaError_t resize(IsomerBatch<T>& batch, size_t new_capacity, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC, int front = -1, int back = -1);
    
    template <Device T>
    cudaError_t reset_convergence_statuses(IsomerBatch<T>& batch, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    
    std::tuple<int, float, float> compare_isomer_arrays(float* a, float* b, int n_isomers, int n_elements_per_isomer, float rtol = 0.0, bool verbose = false, float zero_threshold = 1e-8);
    
    template <Device T>
    void is_close(const IsomerBatch<T>& a, const IsomerBatch<T>& b, float tol = 0.0, bool verbose = false);
    
    template <Device T>
    int count_batch_status(const IsomerBatch<T>& input, const IsomerStatus status);
    
    //Returns the average number of iterations of isomers in the batch that are not empty.
    template <Device T>
    double average_iterations(const IsomerBatch<T>& input);
    
    template <Device T>
    void sort(IsomerBatch<T>& batch, const BatchMember key = IDS, const SortOrder order = ASCENDING);
    
    template <typename T>  
    T mean(std::vector<T>& input);
    std::chrono::nanoseconds sdev(std::vector<std::chrono::nanoseconds>& input);
    std::chrono::microseconds sdev(std::vector<std::chrono::microseconds>& input);
    double sdev(std::vector<double>& input);

}
#endif