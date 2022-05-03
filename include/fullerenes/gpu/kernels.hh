#include "isomerspace_kernel.hh"
#include "cuda_execution.hh"
namespace gpu_kernels{
    namespace isomerspace_forcefield{
        int optimal_batch_size(const int N, const int device_id = 0);
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        float time_spent(); 
        //Uses forcefield optimization to relax the positions of atoms in all isomers of a batch.
        cudaError_t optimize_batch(IsomerBatch& B, const size_t iterations, const size_t max_iterations, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_X0{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        float time_spent(); 
        //Uses spherical projection and 'wrapping a sphere' technique to generate starting coordinates in 3D space for a batch of fullerenes.
        cudaError_t zero_order_geometry(IsomerBatch& B, const device_real_t scalerad, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_tutte{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        float time_spent(); 
        //Computes tutte embedding for an entire batch of fullernes.
        cudaError_t tutte_layout(IsomerBatch& B, const size_t max_iterations = 1000, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }
}
