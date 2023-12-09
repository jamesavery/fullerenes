#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include <chrono>
#include <fstream>
#include <iomanip> 		// Needed for std::setprecision

using namespace chrono;
using namespace chrono_literals;

#include "fullerenes/isomer_queue.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "numeric"
#include "random"
using namespace gpu_kernels;

int main(int argc, char** argv){

    size_t N_start                = argc > 1 ? strtol(argv[1],0,0) : (size_t)20;        // Argument 1: Start of range of N
    size_t N_limit                = argc > 2 ? strtol(argv[2],0,0) : N_start;           // Argument 2: End of range of N
    size_t N_runs                 = argc > 3 ? strtol(argv[3],0,0) : 3;                 // Argument 3: Number of times to run experiment
    size_t warmup                 = argc > 4 ? strtol(argv[4],0,0) : 0;                 // Argument 4: Seconds to warmup GPU

    ofstream out_file   ("ParBenchmark_" + to_string(N_limit) + ".txt");
    ofstream out_std    ("ParBenchmark_STD_" + to_string(N_limit) + ".txt");
    out_file << "Generate, Samples, Update, Dual, Tutte, X0, Optimize \n";

    cuda_benchmark::warmup_kernel(warmup*1s);
    LaunchCtx ctx(0);
    for (size_t N = N_start; N < N_limit+1; N+=2)
    {   
        if(N == 22) continue;
        Graph G;

        auto batch_size = isomerspace_forcefield::optimal_batch_size(N);
        auto n_fullerenes = IsomerDB::number_fullerenes(N);	
        auto sample_size = min(batch_size*1,n_fullerenes);
        if (n_fullerenes < batch_size){
            sample_size = n_fullerenes;
        } else if (n_fullerenes < batch_size*2){
            sample_size = batch_size;
        }else if (n_fullerenes < batch_size*3){
            sample_size = batch_size*2;
        }else if (n_fullerenes < batch_size*4){
            sample_size = batch_size*3;
        }else{
            sample_size = batch_size*4;
        }
        
        if (N == 22) continue;
        std::cout << sample_size << endl;
        
        bool more_to_generate = true;
        
        std::vector<std::chrono::nanoseconds> 
            T_gens(N_runs),
            T_duals(N_runs),
            T_tuttes(N_runs),
            T_X0s(N_runs),
            T_opts(N_runs),
            T_flat(N_runs),
            T_queue(N_runs),
            T_io(N_runs),
            T_hess(N_runs),
            T_spectrum(N_runs),
            T_eccentricity(N_runs),
            T_volume(N_runs);


        auto Nf = N/2 + 2;
        G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
        G.N = Nf;

        auto path = "isomerspace_samples/dual_layout_" + to_string(N) + "_seed_42";
        ifstream isomer_sample(path,std::ios::binary);
        auto fsize = file_size(path);
        std::vector<device_node_t> input_buffer(fsize/sizeof(device_node_t));
        auto available_samples = fsize / (Nf*6*sizeof(device_node_t));
        isomer_sample.read(reinterpret_cast<char*>(input_buffer.data()), Nf*6*sizeof(device_node_t)*available_samples);

        std::vector<int> random_IDs(available_samples);
        std::iota(random_IDs.begin(), random_IDs.end(), 0);
        std::shuffle(random_IDs.begin(), random_IDs.end(), std::mt19937{42});
        std::vector<int> id_subset(random_IDs.begin(), random_IDs.begin()+sample_size);

        
        auto finished_isomers = 0;
        for (size_t l = 0; l < N_runs; l++)
        {   
            finished_isomers = 0;
            IsomerBatch<GPU> batch0(N,sample_size,0);
            IsomerBatch<GPU> batch1(N,sample_size,0);
            IsomerBatch<GPU> batch2(N,sample_size,0);
            IsomerBatch<CPU> h_batch(N,sample_size);
            LaunchCtx ctx(0);
            device_io::IsomerQueue isomer_q(N,0);
            device_io::IsomerQueue isomer_q_cubic(N,0);
            device_io::IsomerQueue OutQueue(N,0);
            OutQueue.resize(sample_size);
            isomer_q_cubic.resize(min(n_fullerenes,10000 + sample_size));
            isomer_q.resize(min(n_fullerenes,sample_size));
            for (int i = 0; i < sample_size; i++){
                for (size_t j = 0; j < Nf; j++){
                    G.neighbours[j].clear();
                    for (size_t k = 0; k < 6; k++) {
                        auto u = input_buffer[id_subset[i]*Nf*6 + j*6 +k];
                        if(u != UINT16_MAX) G.neighbours[j].push_back(u);
                    }
                }
                isomer_q.insert(G, i);
            }
            batch0.clear();
            isomer_q.refill_batch(batch0);
            auto TDual = isomerspace_dual::time_spent();
            isomerspace_dual::dualize(batch0);
            auto TTutte = isomerspace_tutte::time_spent(); T_duals[l] += isomerspace_dual::time_spent() - TDual;
            isomerspace_tutte::tutte_layout(batch0);
            auto TX0 = isomerspace_X0::time_spent(); T_tuttes[l] += isomerspace_tutte::time_spent() - TTutte;
            isomerspace_X0::zero_order_geometry(batch0, 4.0);
            T_X0s[l] += isomerspace_X0::time_spent() - TX0;

            device_io::copy(batch1, batch0);
            while(isomer_q_cubic.get_size() < min(n_fullerenes,10000)){
                isomer_q_cubic.insert(batch1);
                device_io::copy(batch1, batch0);
            }
            auto TFF = high_resolution_clock::now();
            isomerspace_forcefield::optimize<PEDERSEN>(batch0,N*5,N*5);
            auto TFlat = high_resolution_clock::now(); T_opts[l] += TFlat - TFF;
            isomerspace_forcefield::optimize<FLATNESS_ENABLED>(batch1,N*5,N*5);
            T_flat[l] += high_resolution_clock::now() - TFlat;

            CuArray<device_real_t> hessian(batch1.capacity()*N*90);
            CuArray<device_node_t> cols(batch1.capacity()*N*90);
            CuArray<device_real_t> min_eigs(batch1.capacity());
            CuArray<device_real_t> max_eigs(batch1.capacity());
            CuArray<device_real_t> eccentricities(batch1.capacity());
            CuArray<device_real_t> volumes(batch1.capacity());

            auto THess = high_resolution_clock::now();
            isomerspace_hessian::compute_hessians<PEDERSEN>(batch1, hessian, cols);
            T_hess[l] += high_resolution_clock::now() - THess;

            auto TSpectrum = high_resolution_clock::now();
            isomerspace_eigen::spectrum_ends(batch1, hessian, cols, min_eigs, max_eigs);
            T_spectrum[l] += high_resolution_clock::now() - TSpectrum;

            auto TEcc = high_resolution_clock::now();
            isomerspace_properties::eccentricities(batch1, eccentricities);
            T_eccentricity[l] += high_resolution_clock::now() - TEcc;

            auto TVol = high_resolution_clock::now();
            isomerspace_properties::volume_divergences(batch1, volumes);
            T_volume[l] += high_resolution_clock::now() - TVol;

            OutQueue.push_done(batch2,ctx, LaunchPolicy::SYNC);
            auto T0 = high_resolution_clock::now();
            auto j = isomer_q_cubic.get_size(); 
            if(n_fullerenes >= 10000){
                while(j > sample_size){
                    auto T1 = high_resolution_clock::now();
                    isomer_q_cubic.refill_batch(batch2,ctx, LaunchPolicy::SYNC);
                    auto T2 = high_resolution_clock::now(); T_io[l] += T2 - T1;
                    isomerspace_forcefield::optimize<PEDERSEN>(batch2,N*0.5,N*5,ctx, LaunchPolicy::SYNC);
                    auto T3 = high_resolution_clock::now();
                    OutQueue.push_done(batch2,ctx, LaunchPolicy::SYNC);
                    finished_isomers += OutQueue.get_size();
                    j = isomer_q_cubic.get_size();
                    OutQueue.clear(ctx, LaunchPolicy::SYNC);
                    T_io[l] += high_resolution_clock::now() - T3;
                }
            }else{
                while(finished_isomers < n_fullerenes){
                    auto T1 = high_resolution_clock::now();
                    isomer_q_cubic.refill_batch(batch2);
                    auto T2 = high_resolution_clock::now(); T_io[l] += T2 - T1;
                    isomerspace_forcefield::optimize<PEDERSEN>(batch2,N*0.5,N*5);
                    auto T3 = high_resolution_clock::now();
                    OutQueue.push_done(batch2);
                    j = OutQueue.get_size();
                    finished_isomers += j;
                    OutQueue.clear();
                    T_io[l] += high_resolution_clock::now() - T3;
                }
            }
            T_queue[l] += high_resolution_clock::now() - T0;
            ctx.wait();
            if(OutQueue.get_capacity()> 10000 + sample_size) std::cout << "Warning: OutQueue initial capacity exceeded" << std::endl;
            //std::cout << finished_isomers << ": " <<((high_resolution_clock::now() - Tio)/1us) << "us  " << (T_queue[l]/1us) << "us  " << (T_io[l]/1us) << "us  " << (Trefill/1us) << "us  " << (Tpush/1us) << "us  " << (Tget/1us) << "us  " << (Tclear/1us)/(float)finished_isomers << "us  " << std::endl;

        }
        using namespace device_io;
        //Print out runtimes in us per isomer:

        //Remove the maximum and minimum values from the runtime arrays:
        T_gens.erase(T_gens.begin()+std::distance(T_gens.begin(),std::max_element(T_gens.begin(),T_gens.end())));
        T_duals.erase(T_duals.begin()+std::distance(T_duals.begin(),std::max_element(T_duals.begin(),T_duals.end())));
        T_X0s.erase(T_X0s.begin()+std::distance(T_X0s.begin(),std::max_element(T_X0s.begin(),T_X0s.end())));
        T_tuttes.erase(T_tuttes.begin()+std::distance(T_tuttes.begin(),std::max_element(T_tuttes.begin(),T_tuttes.end())));
        T_opts.erase(T_opts.begin()+std::distance(T_opts.begin(),std::max_element(T_opts.begin(),T_opts.end())));
        T_flat.erase(T_flat.begin()+std::distance(T_flat.begin(),std::max_element(T_flat.begin(),T_flat.end())));
        T_queue.erase(T_queue.begin()+std::distance(T_queue.begin(),std::max_element(T_queue.begin(),T_queue.end())));
        T_io.erase(T_io.begin()+std::distance(T_io.begin(),std::max_element(T_io.begin(),T_io.end())));
        T_hess.erase(T_hess.begin()+std::distance(T_hess.begin(),std::max_element(T_hess.begin(),T_hess.end())));
        T_spectrum.erase(T_spectrum.begin()+std::distance(T_spectrum.begin(),std::max_element(T_spectrum.begin(),T_spectrum.end())));
        T_eccentricity.erase(T_eccentricity.begin()+std::distance(T_eccentricity.begin(),std::max_element(T_eccentricity.begin(),T_eccentricity.end())));
        T_volume.erase(T_volume.begin()+std::distance(T_volume.begin(),std::max_element(T_volume.begin(),T_volume.end())));


        //Print out what fraction of the runtime that each component took:
        out_file << N << ", "<< sample_size << ", "<< finished_isomers << ", " << mean(T_gens)/1ns << ", " << mean(T_duals)/1ns <<", " <<  mean(T_X0s)/1ns <<", " << mean(T_tuttes)/1ns<< ", " << mean(T_opts)/1ns <<  ", " << mean(T_flat)/1ns <<  ", " << mean(T_queue)/1ns << ", " << mean(T_io)/1ns << "," << mean(T_hess)/1ns << "," << mean(T_spectrum)/1ns << "," << mean(T_eccentricity)/1ns << "," << mean(T_volume)/1ns << "\n";
        out_std << N << ", "<< sample_size << ", " << finished_isomers << ", " << sdev(T_gens)/1ns << ", " << sdev(T_duals)/1ns <<", " <<  sdev(T_X0s)/1ns <<", " << sdev(T_tuttes)/1ns<< ", " << sdev(T_opts)/1ns << ", " << sdev(T_flat)/1ns << ", " << sdev(T_queue)/1ns << ", " << sdev(T_io)/1ns << "," << sdev(T_hess)/1ns << "," << sdev(T_spectrum)/1ns << "," << sdev(T_eccentricity)/1ns << "," << sdev(T_volume)/1ns << "\n";
        std::cout << N << "  Dual: " << std::fixed << std::setprecision(2) << (float)(mean(T_duals)/1us)/sample_size << "us  Tutte: " << std::fixed << std::setprecision(2) << (float)(mean(T_tuttes)/1us)/sample_size << "us  X0: " << std::fixed << std::setprecision(2) << (float)(mean(T_X0s)/1us)/sample_size << "us  Opt: " << std::fixed << std::setprecision(2) << (float)(mean(T_opts)/1us)/sample_size << "us  Flat: " << std::fixed << std::setprecision(2) << (float)(mean(T_flat)/1us)/sample_size << "us  Queue: " << std::fixed << std::setprecision(2) << (float)(mean(T_queue)/1us)/finished_isomers << "us" << std::fixed << std::setprecision(2) << "  IO: " << (float)(mean(T_io)/1us)/finished_isomers << "us" << std::fixed << std::setprecision(2) << "  Hess: " << (float)(mean(T_hess)/1us)/sample_size << "us" << std::fixed << std::setprecision(2) << "  Spectrum: " << (float)(mean(T_spectrum)/1us)/sample_size << "us" << std::fixed << std::setprecision(2) << "  Eccentricity: " << (float)(mean(T_eccentricity)/1us)/sample_size << "us" << std::fixed << std::setprecision(2) << "  Volume: " << (float)(mean(T_volume)/1us)/sample_size << "us" << std::endl;
     }
    
}
