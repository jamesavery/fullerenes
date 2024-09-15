#include <fullerenes/graph.hh>
#include <fullerenes/sycl-headers/all-kernels.hh>
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/isomerdb.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>
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
    //auto IsomerPerNodeEstimate = NisomersInIsomerspace/std::stoi(getenv("N_TASKS"));
    auto BatchSize =  std::min((int)NisomersInIsomerspace, (int)1280*10);

    size_t Nf = N/2 + 2;
    //std::string filename = "output/full_pipeline_" + std::string(getenv("SLURM_JOB_ID")) + "_" + std::string(getenv("MY_TASK_ID")) + "_" + std::string(getenv("N_TASKS")) + ".csv";
    //ofstream myfile(filename); 

    SyclQueue Q("gpu");

    //size_t n_chunks          = 6; /* Number of chunks per compute node / program instance */ //////3
    //size_t N_chunks          = std::stoi(getenv("N_TASKS"))*n_chunks;                   /* Total number of work chunks */ //

    //auto my_chunks = loadbalanced_chunks(N_chunks,n_chunks,std::stoi(getenv("MY_TASK_ID")));
    //BuckyGen::buckyherd_queue BuckyQ(N,N_chunks,3,
	//			   false,false,my_chunks);

    size_t isomers_in_queue = 0;
    //std::vector<node_t> dual_neighbours(IsomerPerNodeEstimate*1.5*Nf*6);
    //std::vector<node_t> face_degrees(IsomerPerNodeEstimate*1.5*Nf); //We do not know a priori how many isomers our buckygen queue will generate so we reserve space for 1.5 times the number of isomers we want to generate.
    
    
    Graph G(N);

    bool more = true;

    double times_generate = 0.; //Times in nanoseconds.
    double times_memcpy = 0.; //Times in nanoseconds.
    double times_dual = 0.; //Times in nanoseconds.
    double times_tutte = 0.; //Times in nanoseconds.
    double times_project = 0.; //Times in nanoseconds.
    double times_opt = 0.; //Times in nanoseconds.
    double times_hessian = 0.; //Times in nanoseconds.
    double times_spectral_ends = 0.; //Times in nanoseconds.
    double times_spectral = 0.; //Times in nanoseconds.
    double times_spectral_vectors = 0.; //Times in nanoseconds.

    FullereneBatch<real_t,node_t> batch(N, BatchSize);
    SyclVector<float> hessian_buffer(BatchSize*N*90);
    SyclVector<uint16_t> cols_buffer(BatchSize*N*90);
    SyclVector<float> spectral_ends_buffer(BatchSize*2);
    SyclVector<float> spectral_buffer(BatchSize*N*3);
    SyclVector<float> spectral_vectors_buffer(BatchSize*N*3*N*3);


  /*   auto generate_and_fill = [&](IsomerBatch<real_t, node_t>& batch){
        auto isomer_idx = 0;
        sycl::host_accessor acc_dual(batch.dual_neighbours, sycl::write_only);
        sycl::host_accessor acc_degs(batch.face_degrees, sycl::write_only);
        sycl::host_accessor acc_status (batch.statuses, sycl::write_only);
        while (more && isomer_idx < BatchSize)
        {
            more &= BuckyQ.next_fullerene(G);
            if(!more) break;
            for (size_t j = 0; j < Nf; j++)
            {
                for(size_t k = 0; k < G.neighbours[j].size(); k++)
                {
                    acc_dual[isomer_idx*Nf*6 + j*6 + k] = G.neighbours[j][k];
                } 
                if(G.neighbours[j].size() == 5){
                    acc_dual[isomer_idx*Nf*6 + j*6 + 5] = std::numeric_limits<node_t>::max();
                    acc_degs[isomer_idx*Nf + j] = 5;
                } else {
                    acc_degs[isomer_idx*Nf + j] = 6;
                }   

            }
            acc_status[isomer_idx] = IsomerStatus::NOT_CONVERGED;
            isomer_idx++;
            isomers_in_queue++;
        }
    }; */

    DualizeFunctor<real_t,node_t> dualize_V1;
    TutteFunctor<real_t,node_t> tutte_layout;
    SphericalProjectionFunctor<real_t,node_t> spherical_projection;
    ForcefieldOptimizeFunctor<PEDERSEN,real_t,node_t> forcefield_optimize;
    HessianFunctor<PEDERSEN,real_t,node_t> compute_hessians;
    EigenFunctor<EigensolveMode::ENDS, real_t, node_t> eigensolve_ends;
    EigenFunctor<EigensolveMode::FULL_SPECTRUM, real_t, node_t> eigensolve_full;

    auto Nruns = 10;
    fill(batch);
    for(size_t i = 0; i < Nruns; i++){
        auto T1 = std::chrono::steady_clock::now(); times_generate += std::chrono::duration<double, std::nano>(T1 - T1).count();
        auto T2 = std::chrono::steady_clock::now(); times_generate += std::chrono::duration<double, std::nano>(T2 - T1).count();
        //nop_kernel(Q, batch, LaunchPolicy::SYNC);
        auto T3 = std::chrono::steady_clock::now(); times_memcpy += std::chrono::duration<double, std::nano>(T3 - T2).count();
        dualize_V1(Q, batch, LaunchPolicy::SYNC);
        auto T4 = std::chrono::steady_clock::now(); times_dual += std::chrono::duration<double, std::nano>(T4 - T3).count();
        tutte_layout(Q, batch, LaunchPolicy::SYNC);
        auto T5 = std::chrono::steady_clock::now(); times_tutte += std::chrono::duration<double, std::nano>(T5 - T4).count();
        spherical_projection(Q, batch, LaunchPolicy::SYNC);
        auto T6 = std::chrono::steady_clock::now(); times_project += std::chrono::duration<double, std::nano>(T6 - T5).count();
        forcefield_optimize(Q, batch, LaunchPolicy::SYNC, 4*N, 4*N);
        auto T7 = std::chrono::steady_clock::now(); times_opt += std::chrono::duration<double, std::nano>(T7 - T6).count();
        compute_hessians(Q, batch,LaunchPolicy::SYNC, hessian_buffer, cols_buffer);
        auto T8 = std::chrono::steady_clock::now(); times_hessian += std::chrono::duration<double, std::nano>(T8 - T7).count();
        eigensolve_ends(Q, batch, LaunchPolicy::SYNC, hessian_buffer, cols_buffer, 40, spectral_ends_buffer, spectral_buffer);
        auto T9 = std::chrono::steady_clock::now(); times_spectral_ends += std::chrono::duration<double, std::nano>(T9 - T8).count();
        eigensolve_full(Q, batch, LaunchPolicy::SYNC, hessian_buffer, cols_buffer, 0, spectral_ends_buffer, spectral_buffer);
        auto T10 = std::chrono::steady_clock::now(); times_spectral += std::chrono::duration<double, std::nano>(T10 - T9).count();
        //eigensolve<EigensolveMode::FULL_SPECTRUM_VECTORS>(Q, batch, hessian_buffer, cols_buffer, spectral_buffer, LaunchPolicy::SYNC, 40, spectral_vectors_buffer);
        auto T11 = std::chrono::steady_clock::now(); times_spectral_vectors += std::chrono::duration<double, std::nano>(T11 - T10).count();

    }

    std::cout << "N, Nf, BatchSize, JOBID, NTASKS, TASK_ID, FILL_ME_UP_SCOTTY, MEMCPY, DUAL, TUTTE, PROJECT, OPT\n" << N << ", " << Nf << ", " << (BatchSize*Nruns) << ", " << times_generate/(BatchSize*Nruns) << ", " << times_memcpy/(BatchSize*Nruns) << ", " << times_dual/(BatchSize*Nruns) << ", " << times_tutte/(BatchSize*Nruns) << ", " << times_project/(BatchSize*Nruns) << ", " << times_opt/(BatchSize*Nruns) << ", " << times_hessian/(BatchSize*Nruns) << ", " << times_spectral_ends/(BatchSize*Nruns) << ", " << times_spectral/(BatchSize*Nruns) << ", " << times_spectral_vectors/(BatchSize*Nruns) << std::endl;
    std::cout << std::endl;
    return 0;
}
