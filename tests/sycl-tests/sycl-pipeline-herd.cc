#include <fullerenes/triangulation.hh>
#include <fullerenes/sycl-kernels.hh>
#include <fullerenes/isomerdb.hh>
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <string>
#include <algorithm>
#include <random>
#define PRINT_CHECK 1

using namespace sycl;
using namespace std;

typedef float device_real_t;		// This seems pretty hacky. Why is this set locally?
typedef uint16_t device_node_t;
constexpr device_node_t EMPTY_NODE = std::numeric_limits<device_node_t>::max(); /* Do this centrally */

// SYCL (at least on LUMI) doesn't follow the one-platform-multiple-devices structure, but assigns one platform for each GPU,
// and gpu_selector_v only returns a single device despite multiple present. 
// So use more robust, but also more verbose, method to get all available GPUs.
// NB: Wait with multi-GPU for now, as rest of program uses single queue.
vector<queue> gpu_queues()
{
  vector<queue> qs;
  for(auto P: platform::get_platforms())
    for(auto d: P.get_devices())
      if(d.is_gpu()) qs.push_back(queue(d,property::queue::in_order{}));

  return qs;
}

// TODO: This should be rolled into buckyherd_queue
vector<size_t> loadbalanced_chunks(size_t N_chunks, size_t n_chunks, size_t my_node_idx=0, int seed=42)
{
  vector<size_t> all_chunks(N_chunks), my_chunks(n_chunks);
  for(size_t i=0;i<N_chunks;i++) all_chunks[i] = i;  
  auto rng = std::default_random_engine(seed); // Deterministic randomness - we need all compute nodes to agree on shuffle.
  shuffle(all_chunks.begin(),all_chunks.end(),rng);
  for(size_t i=0;i<n_chunks;i++) my_chunks[i] = all_chunks[my_node_idx*n_chunks + i];
  
  return my_chunks;
}

  

int main(int argc, char** argv) {
  
  if(argc<2 || (argv[1] != string("cpu") && argv[1] != string("gpu"))){
    fprintf(stderr,"Syntax: %s <device_type=cpu|gpu> [N:60] [batch_size:8192] [N_nodes:1] [my_node_idx:0] [workers_per_node:3] [chunks_per_worker:1] [IPR:0] [only_symmetric:0]\n",argv[0]);
    return -1;
  }
  
  size_t N  = argc>2? strtol(argv[2],0,0) : 60;
  size_t Nf = N/2 + 2;  
  size_t Nisomers = IsomerDB::number_isomers(N);
  assert(N != 22 && N>=20 && (N%2 == 0));

  // TODO: Implement small command line parser for named arguments instead of
  //       1M awkward positional command line arguments. Nice for all programs
  //       that use the library. But no need to add extra dependency.
  // Command line configuration
  size_t batch_size        = argc>3 ? strtol(argv[3],0,0) : std::min<size_t>(8192,Nisomers);
  size_t N_nodes           = argc>4 ? strtol(argv[4],0,0) : 1;
  size_t my_node_idx       = argc>5 ? strtol(argv[5],0,0) : 0;
  size_t workers_per_node  = argc>6 ? strtol(argv[6],0,0) : 3;
  size_t chunks_per_worker = argc>7 ? strtol(argv[7],0,0) : 1;
  bool   IPR               = argc>8 ? strtol(argv[8],0,0) : 0;
  bool   only_symmetric    = argc>9 ? strtol(argv[9],0,0) : 0;
  
  size_t n_chunks          = workers_per_node*chunks_per_worker; /* Number of chunks per compute node / program instance */
  size_t N_chunks          = N_nodes*n_chunks;                   /* Total number of work chunks */

  auto my_chunks = loadbalanced_chunks(N_chunks,n_chunks,my_node_idx);
  BuckyGen::buckyherd_queue BuckyQ(N,N_chunks,workers_per_node,
				   IPR,only_symmetric,my_chunks);
  
  // Set up SYCL queue
  string device_type = argv[1];
  auto selector = device_type == "gpu"? gpu_selector_v : cpu_selector_v;
  queue Q(selector, property::queue::in_order());  
    
  printf("Running on device: %s with %d compute units.\n",
          Q.get_device().get_info<info::device::name>().c_str(),
          Q.get_device().get_info<info::device::max_compute_units>());           

  IsomerBatch<device_real_t,device_node_t> batch(N, batch_size);
  Triangulation G(N);

  size_t ii = 0;      
  {
    printf("Gathering %ld C%ld isomers\n",batch_size,N);  
    sycl::host_accessor acc_dual(batch.dual_neighbours, sycl::write_only);
    sycl::host_accessor acc_degs(batch.face_degrees, sycl::write_only);
    sycl::host_accessor acc_status (batch.statuses, sycl::write_only);
    for (; ii < batch_size; ii++)
      {
        auto more = BuckyQ.next_fullerene(G);
        if(!more) break;

        for (size_t j = 0; j < Nf; j++)
	  {
            for(size_t k = 0; k < G.neighbours[j].size(); k++)
	      {
                acc_dual[ii*Nf*6 + j*6 + k] = G.neighbours[j][k];
	      } 
            if(G.neighbours[j].size() == 5){
	      acc_dual[ii*Nf*6 + j*6 + 5] = std::numeric_limits<device_node_t>::max();
	      acc_degs[ii*Nf + j] = 5;
            } else {
	      acc_degs[ii*Nf + j] = 6;
            }   

	  }
        acc_status[ii] = IsomerStatus::NOT_CONVERGED;
      }
    BuckyQ.stop_all();   
  }
  printf("Processing %ld C%ld isomers\n",ii,N);          
  printf("\tdualize\n");              dualise(Q, batch, LaunchPolicy::SYNC);
  printf("\ttutte_layout\n");         tutte_layout(Q, batch, LaunchPolicy::SYNC);
  printf("\tspherical_projection\n"); spherical_projection(Q, batch, LaunchPolicy::SYNC);
  printf("\tforcefield_optimize\n");  forcefield_optimise(Q, batch, 5*N, 5*N, LaunchPolicy::SYNC);
    
  sycl::buffer<device_real_t, 1> hessians(range<1>(N*90*batch_size));
  sycl::buffer<device_node_t, 1> cols(range<1>(N*90*batch_size));
  printf("\tcompute_hessians\n"); compute_hessians(Q, batch, hessians, cols, LaunchPolicy::SYNC);
  
  return 0;
  
}
    
