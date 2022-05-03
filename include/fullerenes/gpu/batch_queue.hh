#include "isomerspace_kernel.hh"
#include "cuda_execution.hh"

namespace cuda_io{
struct BatchQueue
{    
public:
    BatchQueue(const size_t N);
    ~BatchQueue();

    //Queries about whether the queue is empty.
    bool is_empty() const {return *props.size == 0;};

    //Primitive getter functions.
    int get_size() const {return *props.size;};
    int get_front() const {return *props.front;};
    int get_back() const {return *props.back;};
    int get_capacity() const {return *props.capacity;};
    int get_requests() const {return *props.requests;};

    /** Refill batch 
     * @param target Batch to inspect and fill new isomers into in-place.
     * Checks the status of each isomer in the target batch, if the status is either of: CONVERGED / FAILED or EMPTY,
     * a new isomer will be copied from the queue if available, this happens in parallel.
     * **/
    cudaError_t refill_batch(IsomerBatch& target, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    
    //Insert an entire IsomerBatch into the queue.
    cudaError_t insert(IsomerBatch& input_batch, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC ,const bool insert_2d = true);
    
    //Insert a single Polyhedron.
    cudaError_t insert(const Polyhedron& isomer, const size_t ID, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC, const bool insert_2d = true);
    
    //Insert a single PlanarGraph. This function exists because Polyhedrons are expensive objects to convert to, which occurs implicitly otherwise.
    cudaError_t insert(const PlanarGraph& isomer, const size_t ID, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC, const bool insert_2d = true);
  
    //Resize function, called primarily internally when more space is required. Does what you think it does.
    cudaError_t resize(const size_t new_capacity,const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy = LaunchPolicy::SYNC);
    
    //Simple data structure to contain all the necessary queue counters, makes for smaller function signatures!
    struct QueueProperties
    {
      //Uses pointers because we need to use cuda managed memory such that these counters can be manipulated and accessed on host side and device side.
      int* front{};  
      int* back{};  
      int* size{};  
      int* capacity{};  
      int* requests{};
    };
    
    const size_t N;
    //The underlying data containers are buffers on the host and device side respecitvely, the queue properties essentially wrap these containers to give the functionality of a circular queue.
    IsomerBatch device_batch = IsomerBatch(N,1,DEVICE_BUFFER);
    IsomerBatch host_batch   = IsomerBatch(N,1,HOST_BUFFER);

private:
    bool is_mirrored_on_device = false;
    QueueProperties props;

    //Need to keep track of where the updated version of the batch exists, if we manipulate the host batch, the device is no longer up to date and vice versa.
    bool is_host_updated = false, is_device_updated = false;

    //Internal method for transfering all data to the device if the device is not up to date.
    cudaError_t to_device(const LaunchCtx& ctx = LaunchCtx());
    //Internal method for transfering all data to the host if the host is not up to date.
    cudaError_t to_host(const LaunchCtx& ctx = LaunchCtx());
};
}
