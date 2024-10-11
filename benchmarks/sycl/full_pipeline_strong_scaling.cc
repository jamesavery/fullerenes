#include <fullerenes/graph.hh>
#include <fullerenes/sycl-headers/all-kernels.hh>
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/isomerdb.hh>
#include <string>
#include <algorithm>
#include <random>

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

    size_t N                 = argc>1? stoi(argv[1]) : 200;  // Full isomerspace to process
    size_t N_TASKS_MAX       = argc>2? stoi(argv[2]) : 8; // Total number of work chunks. Should be final/max N_TASKS times 3 to not be dominated by buckygen overhead.
    size_t workers_per_task  = argc>3? stoi(argv[3]) : 3; // How many parallel buckygens are there in the buckyherd_queue?
    size_t chunks_per_worker = argc>4? stoi(argv[4]) : 1; // How many work chunks should each worker get on average?

    size_t NisomersInIsomerspace = IsomerDB::number_isomers(N);

    // We get N_TASKS from environment. Why? Because.
    const char *N_TASKS_str    = getenv("N_TASKS");
    const char *MY_TASK_ID_str = getenv("MY_TASK_ID");

    size_t N_TASKS            = N_TASKS_str?    std::stoi(N_TASKS_str)    : 1;
    size_t MY_TASK_ID         = MY_TASK_ID_str? std::stoi(MY_TASK_ID_str) : 0;
							 
    auto IsomerPerNodeEstimate = NisomersInIsomerspace/N_TASKS;
    auto BatchSize = std::min<int>(IsomerPerNodeEstimate, (1<<17));

    size_t Nf = N/2 + 2;
    std::string filename = "output/full_pipeline_" + std::string(getenv("SLURM_JOB_ID")) + "_" + to_string(MY_TASK_ID) + "_" + to_string(N_TASKS) + ".csv";
    ofstream myfile(filename); 

    //sycl::queue Q = sycl::queue(sycl::gpu_selector_v, sycl::property::queue::in_order{});
    SyclQueue Q("gpu");

    // Now we make the total N_chunks a fixed parameter and vary n_chunks, the number of chunks per task (GCD).
    size_t N_chunks = N_TASKS_MAX * workers_per_task * chunks_per_worker; 
    size_t n_chunks = N_chunks/N_TASKS;    /* Number of chunks per compute node / program instance */ //////3
    
    auto my_chunks = loadbalanced_chunks(N_chunks,n_chunks,MY_TASK_ID);
    BuckyGen::buckyherd_queue BuckyQ(N,N_chunks,workers_per_task,
				   false,false,my_chunks);

    size_t isomers_in_queue = 0;
    Graph G(N);

    bool more = true;

    double times_generate; //Times in nanoseconds.
    double times_memcpy; //Times in nanoseconds.
    double times_dual; //Times in nanoseconds.
    double times_tutte; //Times in nanoseconds.
    double times_project; //Times in nanoseconds.
    double times_opt; //Times in nanoseconds.
    FullereneBatch<real_t,node_t> batch(N, BatchSize);


    auto generate_and_fill = [&](FullereneBatch<real_t, node_t>& batch){
        auto isomer_idx = 0;
        auto acc_dual = batch.d_.A_dual_;
        auto acc_degs = batch.d_.deg_;
        auto acc_status  = batch.m_.flags_;
        while (more && isomer_idx < BatchSize)
        {
            more &= BuckyQ.next_fullerene(G);
            if(!more) break;
            for (size_t j = 0; j < Nf; j++)
            {
                for(size_t k = 0; k < G.neighbours[j].size(); k++)
                {
                    acc_dual[isomer_idx*Nf + j][k] = G.neighbours[j][k];
                } 
                if(G.neighbours[j].size() == 5){
                    acc_dual[isomer_idx*Nf + j][5] = std::numeric_limits<node_t>::max();
                    acc_degs[isomer_idx*Nf + j] = 5;
                } else {
                    acc_degs[isomer_idx*Nf + j] = 6;
                }   

            }
            acc_status[isomer_idx] = StatusFlag::DUAL_INITIALIZED;
            isomer_idx++;
            isomers_in_queue++;
        }
    };

    DualizeFunctor<real_t, node_t> dualize_V1;
    TutteFunctor<real_t, node_t> tutte_layout;
    SphericalProjectionFunctor<real_t, node_t> spherical_projection;
    ForcefieldOptimizeFunctor<PEDERSEN, real_t, node_t> forcefield_optimize;



    while(more){
        auto T1 = std::chrono::steady_clock::now();
        generate_and_fill(batch);
        auto T2 = std::chrono::steady_clock::now(); times_generate = std::chrono::duration<double, std::nano>(T2 - T1).count();
        //nop_kernel(Q, batch, LaunchPolicy::SYNC);
        auto T3 = std::chrono::steady_clock::now(); times_memcpy = std::chrono::duration<double, std::nano>(T3 - T2).count();
        dualize_V1(Q, batch, LaunchPolicy::SYNC);
        auto T4 = std::chrono::steady_clock::now(); times_dual = std::chrono::duration<double, std::nano>(T4 - T3).count();
        tutte_layout(Q, batch, LaunchPolicy::SYNC);
        auto T5 = std::chrono::steady_clock::now(); times_tutte = std::chrono::duration<double, std::nano>(T5 - T4).count();
        spherical_projection(Q, batch, LaunchPolicy::SYNC);
        auto T6 = std::chrono::steady_clock::now(); times_project = std::chrono::duration<double, std::nano>(T6 - T5).count();
        forcefield_optimize(Q, batch, LaunchPolicy::SYNC, 5*N, 5*N);
        auto T7 = std::chrono::steady_clock::now(); times_opt = std::chrono::duration<double, std::nano>(T7 - T6).count();
    }

    myfile << "N, Nf, BatchSize, JOBID, NTASKS, TASK_ID, FILL_ME_UP_SCOTTY, MEMCPY, DUAL, TUTTE, PROJECT, OPT\n" << 
    N << ", " << Nf << ", " << isomers_in_queue << ", " << getenv("SLURM_JOB_ID") << ", " << N_TASKS << ", " << MY_TASK_ID << ", " << times_generate/isomers_in_queue << ", " << times_memcpy/isomers_in_queue << ", " << times_dual/isomers_in_queue << ", " << times_tutte/isomers_in_queue << ", " << times_project/isomers_in_queue << ", " << times_opt/isomers_in_queue << "\n";
    return 0;
}
