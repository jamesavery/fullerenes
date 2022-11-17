#include "isomerspace_kernel.hh"
#include "cuda_execution.hh"
#include "fullerenes/gpu/cu_array.hh"
namespace cuda_io{
struct IsomerQueue
{    
public:
    IsomerQueue(const size_t N, int device = 0);
    ~IsomerQueue();
    int m_device = 0;
    //Queries about whether the queue is empty.
    bool is_empty() const {return *props.size == 0;};

    //Primitive getter functions.
    int get_size() const;
    int get_front() const;
    int get_back() const;
    int get_capacity() const;
    int get_requests() const;

    /** Refill batch 
     * @param target Batch to inspect and fill new isomers into in-place.
     * Checks the status of each isomer in the target batch, if the status is either of: CONVERGED / FAILED or EMPTY,
     * a new isomer will be copied from the queue if available, this happens in parallel.
     * **/
    cudaError_t refill_batch(IsomerBatch& target, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    
    //Insert an entire IsomerBatch into the queue.
    cudaError_t insert(IsomerBatch& input_batch, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC ,const bool insert_2d = true);
    
    //Push all the finished isomers from a batch into the queue.
    cudaError_t push(IsomerBatch& input_batch, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC );

    Polyhedron pop(const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

    //Insert a single Polyhedron.
    cudaError_t insert(const Polyhedron& isomer, const size_t ID, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC, const bool insert_2d = true);
    
    //Insert a single PlanarGraph. This function exists because Polyhedrons are expensive objects to convert to, which occurs implicitly otherwise.
    cudaError_t insert(const PlanarGraph& isomer, const size_t ID, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC, const bool insert_2d = true);

    //Insert a single Graph. Intended to insert the graph objects (containing a dual-neighbour list) generated by BuckyGen
    cudaError_t insert(const Graph& isomer, const size_t ID, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

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

    friend std::ostream& operator<<(std::ostream& os, const IsomerQueue& a);

private:
    bool is_mirrored_on_device = false;
    QueueProperties props;

    //A pointer to a pointer to a piece of memory reserved for performing the scan operation in the refill function.
    int* g_scan_array{};
    
    //While we wish to allow the user to enqueue queue operations on any cuda stream, we simultaneously wish to ensure ordering of all queue operations.
    //Therefore we place an event in every stream after each kernel launch or memcopy and wait for any previous event launched on the queue.
    cudaEvent_t latest_operation;
    
    //Need to keep track of where the updated version of the batch exists, if we manipulate the host batch, the device is no longer up to date and vice versa.
    bool is_host_updated = false, is_device_updated = false;

    //Internal method for transfering all data to the device if the device is not up to date.
    cudaError_t to_device(const LaunchCtx& ctx = LaunchCtx());
    //Internal method for transfering all data to the host if the host is not up to date.
    cudaError_t to_host(const LaunchCtx& ctx = LaunchCtx());
};
}
