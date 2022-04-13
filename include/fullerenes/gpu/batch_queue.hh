#include "isomerspace_kernel.hh"

struct BatchQueue
{    
public:
    BatchQueue();

    bool is_empty() const {return h_props.size == 0;};
    size_t get_size() const {return h_props.size == 0;};
    cudaError_t refill_batch(IsomerBatch& target, const cudaStream_t stream = NULL);
    cudaError_t insert(const IsomerBatch& input_batch);
    cudaError_t insert(const Polyhedron& isomer);
    cudaError_t insert(const FullereneGraph& isomer);
    cudaError_t insert(const FullereneDual& isomer);

    struct Properties{
        IsomerBatch batch;
        size_t front{};
        size_t back{};
        size_t capacity{};
        size_t size{};
    };

private:
    Properties h_props;
    Properties d_props;
};