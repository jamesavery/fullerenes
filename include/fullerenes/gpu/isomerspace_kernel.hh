#pragma once
#ifndef ISOMERKERNEL_HEADER
#define ISOMERKERNEL_HEADER

#include <inttypes.h>
#include <string>
#include <queue>
#include <map>
#include <fullerenes/fullerenegraph.hh> 
#include <fullerenes/polyhedron.hh>
#include "cu_array.hh"
#include "gpudatastruct.hh"
#include <cuda_runtime.h>
#include <tuple>

#define GPU_REAL float       
#define GPU_REAL3 float3
#define GPU_REAL2 float2
#define GPU_NODE2 ushort2
#define GPU_NODE3 ushort3
#define sin sinf		// change to sin for double precision, sinl for long double
#define cos cosf 		// change to cos for double precision, cosl for long double
#define NODE_MAX UINT16_MAX
#define MAKE_NODE3 make_ushort3
#define Block_Size_Pow_2 256
#define SEMINARIO_FORCE_CONSTANTS 1
#define USE_MAX_NORM 0
#define REDUCTION_METHOD 0
#define LINESEARCH_METHOD GSS

#include "isomer_batch.hh"

class IsomerspaceKernel {
public:
  enum QueueMode {HOST_QUEUE, DEVICE_QUEUE};

  typedef GPU_REAL device_real_t;
  typedef uint16_t device_node_t;
  
  size_t get_batch_size()const{return batch_size;}
  void insert_isomer(const Polyhedron& P, const size_t ID) {insert_queue.push({ID,P});}  //Pushes P to insert_queue. 
  void output_isomer(const size_t idx);                      //Pops the idx'th fullerene from the IsomerBatch to the output_queue.
  
  size_t get_queue_size()const{return insert_queue.size();}
  size_t get_device_queue_size()const{return device_queue_size;}
  size_t get_batch_capacity();          //Uses Cuda API calls to determine the amount of fullerenes of a given size N, that can be optimized simultaneously.
  size_t get_device_capacity(size_t i)const {return device_capacities[i];}
  size_t get_isomer_size() {return N;}
  size_t get_device_count() {return device_count;}
  
  void clear_batch() {
    batch_size = 0; 
    for (size_t i = 0; i < device_count; i++)
    { 
      batch_sizes[i] = 0;
      h_batch[i].n_isomers = 0;
      for (size_t j = 0; j < device_capacities[i]; j++) h_batch[i].statuses[j] = EMPTY;
    }    
  }
  virtual void check_batch(size_t steps) {}                   //Checks convergence properties of current batch, calculates mean and std of relative bond, angle and dihedral errors of the current batch.
  virtual void update_batch();
  void synchronize();
  void output_isomer(size_t i, size_t idx);
  void insert_isomer(size_t i, size_t idx);

  void output_batch_to_queue();
  void insert_queued_isomers();

  void copy_metadata();
  static void kernel_to_kernel_copy(IsomerspaceKernel& source, IsomerspaceKernel& destination);


  IsomerspaceKernel(const size_t N, void* kernel); //Simple constructor allocates memory for structs on device and host call this once at the beginning of program, 
                                          //also serves to initialize the cuda default context, which would otherwise be destroyed and recreated every time memory is freed and reallocated.
  ~IsomerspaceKernel();                   //Destructor, calls free and delete on GPU and CPU buffers respectively.

  std::queue<std::pair<size_t, Polyhedron>> output_queue;
  std::queue<std::pair<size_t, Polyhedron>> insert_queue;
  std::vector<IsomerBatch> h_batch, d_batch, d_output_batch;

protected:
  size_t N = 0;                           //Number of verices in each isomer, #of carbon atoms.
  size_t batch_capacity = 0;              //Number of fullerenes that will fit on the GPU concurrently.
  size_t batch_size = 0;                  //Current number of fullerenes copied to GPU.
  size_t shared_memory_bytes = 0;         //Amount of L1 cache to allocate per block.
  size_t converged_count = 0;             //Total number of converged isomers optimized by this object.
  size_t failed_count = 0;                //Total number of failed isomers optimized by this object.   
  size_t device_queue_size = 0;
  QueueMode queue_mode = DEVICE_QUEUE;
  
  void* kernel_pointer;
  std::vector<cudaStream_t> main_stream;
  std::vector<cudaStream_t> copy_stream;

  int device_count;

  //h_buffer and d_buffer are mirrors and reflect what will be computed, d_input_batch and d_output_batch exist for linking kernels together in a pipeline. 

  //The reason why these are vectors is that we might have mutliple devices (GPUs).
  std::vector<int> device_capacities;
  std::vector<int> batch_sizes;
  std::vector<device_real_t*> global_reduction_arrays;  //Array used to communicate across blocks.
  std::vector<std::queue<int>> index_queue;             //Contains indices for valid parts of GPU memory to push new fullerenes to.
};

#endif
 
