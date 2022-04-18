#include "isomerspace_kernel.hh"
namespace gpu_kernels{
    namespace isomerspace_forcefield{
        cudaError_t optimize_batch(IsomerBatch& B, const size_t iterations, const size_t max_iterations, const cudaStream_t stream = NULL);
    }

    namespace isomerspace_X0{
        cudaError_t zero_order_geometry(IsomerBatch& B, const device_real_t scalerad, const cudaStream_t stream = NULL);
    }

    namespace isomerspace_tutte{
        cudaError_t tutte_layout(IsomerBatch& B, const size_t max_iterations = 1000, const cudaStream_t = NULL);
    }
}