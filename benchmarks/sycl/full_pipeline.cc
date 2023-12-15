#include <fullerenes/graph.hh>
#include <fullerenes/sycl-kernels.hh>
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/isomerdb.hh>
#include <string>
#include <algorithm>
#include <random>
using namespace sycl;

template<typename T>
T stddev(const std::vector<T>& data)
{
  if(data.size() <= 1) return 0;
  
    // Calculate the mean
    T mn = mean(data);

    // Calculate the sum of squared differences from the mean
    T sum_of_squares = 0.0;
    for (const T& value : data)
    {
        T diff = value - mn;
        sum_of_squares += diff * diff;
    }

    // Calculate the variance and return the square root
    T variance = sum_of_squares / (data.size() - 1);
    return std::sqrt(variance);
}

vector<size_t> loadbalanced_chunks(size_t N_chunks, size_t n_chunks, size_t my_node_idx=0, int seed=42)
{
  vector<size_t> all_chunks(N_chunks), my_chunks(n_chunks);
  for(size_t i=0;i<N_chunks;i++) all_chunks[i] = i;  
  auto rng = std::default_random_engine(42); // Deterministic randomness - we need all compute nodes to agree on shuffle.
  shuffle(all_chunks.begin(),all_chunks.end(),rng);
  for(size_t i=0;i<n_chunks;i++) my_chunks[i] = all_chunks[my_node_idx*n_chunks + i];
  
  return my_chunks;
}


int main(int argc, char** argv) {
    typedef float real_t;
    typedef uint16_t node_t;

    size_t N  = argc>1 ? strtol(argv[1],0,0) : 200;
    size_t NumNodes = argc>2 ? strtol(argv[2],0,0) : 20000000;
    size_t NisomersInIsomerspace = IsomerDB::number_isomers(N);
    auto BatchSize = NisomersInIsomerspace/std::stoi(getenv("N_TASKS"));

    size_t Nf = N/2 + 2;
    std::string filename = "output/full_pipeline_" + std::string(getenv("SLURM_JOB_ID")) + "_" + std::string(getenv("MY_TASK_ID")) + "_" + std::string(getenv("N_TASKS")) + ".csv";
    ofstream myfile(filename); 

    sycl::queue Q = sycl::queue(sycl::gpu_selector_v, sycl::property::queue::in_order{});

    size_t n_chunks          = 6; /* Number of chunks per compute node / program instance */ //////3
    size_t N_chunks          = std::stoi(getenv("N_TASKS"))*n_chunks;                   /* Total number of work chunks */ //

    auto my_chunks = loadbalanced_chunks(N_chunks,n_chunks,std::stoi(getenv("MY_TASK_ID")));
    BuckyGen::buckyherd_queue BuckyQ(N,N_chunks,3,
				   false,false,my_chunks);

    size_t isomers_in_queue = 0;
    std::vector<node_t> dual_neighbours(BatchSize*1.5*Nf*6);
    std::vector<node_t> face_degrees(BatchSize*1.5*Nf); //We do not know a priori how many isomers our buckygen queue will generate so we reserve space for 1.5 times the number of isomers we want to generate.
    
    
    Graph G(N);

    bool more = true;

    double times_generate; //Times in nanoseconds.
    double times_memcpy; //Times in nanoseconds.
    double times_dual; //Times in nanoseconds.
    double times_tutte; //Times in nanoseconds.
    double times_project; //Times in nanoseconds.
    double times_opt; //Times in nanoseconds.


    auto start = std::chrono::steady_clock::now();
    while (more)
    {
        more &= BuckyQ.next_fullerene(G);
        if(!more) break;
        for (size_t j = 0; j < Nf; j++)
        {
            for(size_t k = 0; k < G.neighbours[j].size(); k++)
            {
                dual_neighbours[isomers_in_queue*Nf*6 + j*6 + k] = G.neighbours[j][k];
            } 
            if(G.neighbours[j].size() == 5){
                dual_neighbours[isomers_in_queue*Nf*6 + j*6 + 5] = std::numeric_limits<node_t>::max();
                face_degrees[isomers_in_queue*Nf + j] = 5;
            } else {
                face_degrees[isomers_in_queue*Nf + j] = 6;
            }   

        }
        //acc_status[ii] = IsomerStatus::NOT_CONVERGED;
        isomers_in_queue++;
    }

    auto end = std::chrono::steady_clock::now(); times_generate = std::chrono::duration<double, std::nano>(end - start).count();

    IsomerBatch<real_t,node_t> batch(N, isomers_in_queue);

    auto fill = [&](IsomerBatch<real_t,node_t>& batch)
    {
        sycl::host_accessor acc_dual(batch.dual_neighbours, sycl::write_only);
        sycl::host_accessor acc_degs(batch.face_degrees, sycl::write_only);
        sycl::host_accessor acc_status (batch.statuses, sycl::write_only);
        bool more = true;
        for(size_t ii = 0; ii < isomers_in_queue; ii++){
            for (size_t j = 0; j < Nf; j++)
            {
                for(size_t k = 0; k < 6; k++)
                {
                    acc_dual[ii*Nf*6 + j*6 + k] = dual_neighbours[ii*Nf*6 + j*6 + k];
                    acc_degs[ii*Nf + j] = face_degrees[ii*Nf + j];
                } 
            }
            acc_status[ii] = IsomerStatus::NOT_CONVERGED;
            isomers_in_queue++;
        }
    };


    auto T1 = std::chrono::steady_clock::now();
    fill(batch);
    auto T2 = std::chrono::steady_clock::now(); times_generate = std::chrono::duration<double, std::nano>(T2 - T1).count();
    nop_kernel(Q, batch, LaunchPolicy::SYNC);
    auto T3 = std::chrono::steady_clock::now(); times_memcpy = std::chrono::duration<double, std::nano>(T3 - T2).count();
    dualize_V1(Q, batch, LaunchPolicy::SYNC);
    auto T4 = std::chrono::steady_clock::now(); times_dual = std::chrono::duration<double, std::nano>(T4 - T3).count();
    tutte_layout(Q, batch, LaunchPolicy::SYNC);
    auto T5 = std::chrono::steady_clock::now(); times_tutte = std::chrono::duration<double, std::nano>(T5 - T4).count();
    spherical_projection(Q, batch, LaunchPolicy::SYNC);
    auto T6 = std::chrono::steady_clock::now(); times_project = std::chrono::duration<double, std::nano>(T6 - T5).count();
    forcefield_optimize(Q, batch, 5*N, 5*N, LaunchPolicy::SYNC);
    auto T7 = std::chrono::steady_clock::now(); times_opt = std::chrono::duration<double, std::nano>(T7 - T6).count();

    myfile << "N, Nf, BatchSize, JOBID, NTASKS, TASK_ID, FILL_ME_UP_SCOTTY, MEMCPY, DUAL, TUTTE, PROJECT, OPT\n" << 
    N << ", " << Nf << ", " << isomers_in_queue << ", " << getenv("SLURM_JOB_ID") << ", " << getenv("N_TASKS") << ", " << getenv("MY_TASK_ID") << ", " << times_generate/isomers_in_queue << ", " << times_memcpy/isomers_in_queue << ", " << times_dual/isomers_in_queue << ", " << times_tutte/isomers_in_queue << ", " << times_project/isomers_in_queue << ", " << times_opt/isomers_in_queue << "\n";
    return 0;
}
