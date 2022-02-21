#pragma once
#include <inttypes.h>
#include <string>
#include <queue>
#include <map>
#include <fullerenes/fullerenegraph.hh> 
#include <fullerenes/polyhedron.hh>

#define GPU_REAL float       
#define GPU_REAL3 float3
#define GPU_REAL2 float2
#define GPU_NODE2 ushort2
#define GPU_NODE3 ushort3
#define MAKE_NODE3 make_ushort3
#define Block_Size_Pow_2 256
#define SEMINARIO_FORCE_CONSTANTS 1
#define USE_MAX_NORM 0
#define REDUCTION_METHOD 0
#define LINESEARCH_METHOD GSS
  

template <typename T>
class IsomerspaceKernel {
public:
  typedef GPU_REAL device_real_t;
  typedef uint16_t device_node_t;
  
  enum IsomerStatus {CONVERGED, FAILED, NOT_CONVERGED, EMPTY};

  size_t get_batch_size()const{return batch_size;}
  void insert_isomer(const T& P, const size_t ID) {insert_queue.push({ID,P});};  //Pushes P to insert_queue. 
  void output_isomer(const size_t idx);                      //Pops the idx'th fullerene from the IsomerBatch to the output_queue.
  
  size_t get_queue_size()const{return insert_queue.size();}
    
  void clear_batch(){batch_size=0;} //Clears batch, this is required after every batch is finished, effectively resets the position of pointer to GPU memory

  size_t get_batch_capacity();          //Uses Cuda API calls to determine the amount of fullerenes of a given size N, that can be optimized simultaneously.
  virtual void check_batch() {};                   //Checks convergence properties of current batch, calculates mean and std of relative bond, angle and dihedral errors of the current batch.
  virtual void update_batch() {};
  virtual void eject_isomer(size_t i, size_t idx) {};

  IsomerspaceKernel(const size_t N, void* kernel); //Simple constructor allocates memory for structs on device and host call this once at the beginning of program, 
                                          //also serves to initialize the cuda default context, which would otherwise be destroyed and recreated every time memory is freed and reallocated.
  ~IsomerspaceKernel();              //Destructor, calls free and delete on GPU and CPU buffers respectively.

  std::queue<std::pair<device_node_t, T>> output_queue;
  std::queue<std::pair<device_node_t, T>> insert_queue;

protected:
  size_t N = 0;                           //Number of verices in each isomer, #of carbon atoms.
  size_t batch_capacity = 0;              //Number of fullerenes that will fit on the GPU concurrently.
  size_t batch_size = 0;                  //Current number of fullerenes copied to GPU.
  size_t shared_memory_bytes = 0;         //Amount of L1 cache to allocate per block.
  
  void* kernel_pointer;
  void* cuda_streams;
  int device_count;

  //The reason why these are vectors is that we might have mutliple devices (GPUs).
  std::vector<int> device_capacities;
  std::vector<int> batch_sizes;
  std::vector<device_real_t*> global_reduction_arrays;  //Array used to communicate across blocks.
  std::vector<std::queue<int>> index_queue;             //Contains indices for valid parts of GPU memory to push new fullerenes to.
};


