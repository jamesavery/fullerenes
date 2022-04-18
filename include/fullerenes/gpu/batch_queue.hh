#include "isomerspace_kernel.hh"

namespace cuda_io{
struct BatchQueue
{    
public:
    BatchQueue(const size_t N);
    ~BatchQueue();

    bool is_empty() const {return *props.size == 0;};
    size_t get_size() const {return *props.size;};
    cudaError_t refill_batch(IsomerBatch& target, const cudaStream_t stream = NULL);
    cudaError_t insert(IsomerBatch& input_batch, const cudaStream_t stream = NULL, const bool insert_2d = true);
    cudaError_t insert(const Polyhedron& isomer, const size_t ID, const cudaStream_t stream = NULL, const bool insert_2d = true);
    
    struct QueueProperties
    {
      unsigned int* front{};  
      unsigned int* back{};  
      unsigned int* size{};  
      unsigned int* capacity{};  
      unsigned int* requests{};
    };
    

private:
    const size_t N;
    bool is_mirrored_on_device = false;
    IsomerBatch host_batch   = IsomerBatch(N,1,HOST_BUFFER);
    IsomerBatch device_batch = IsomerBatch(N,1,DEVICE_BUFFER);
    QueueProperties props;

    bool is_host_updated = false, is_device_updated = false;
    cudaError_t to_device(const cudaStream_t stream = NULL);
    cudaError_t to_host(const cudaStream_t stream = NULL);
};
}
