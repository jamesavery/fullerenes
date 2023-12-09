#include <csignal>
#include <sys/stat.h>
#include <limits.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <thread>
#include <numeric>
#include <future>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/progress_bar.hh"
#include "fullerenes/gpu/benchmark_functions.hh"

using namespace std;
using namespace std::chrono;
#include "fullerenes/isomer_queue.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/isomerdb.hh"

using namespace gpu_kernels;

int main(int ac, char **argv)
{
    const size_t N                = strtol(argv[1],0,0);     // Argument 1: Number of vertices N
    string output_dir     = ac>=3? argv[2] : "output";    // Argument 2: directory to output files to
    int IPR               = ac>=4? strtol(argv[3],0,0):0; // Argument 3: Only generate IPR fullerenes?
    int only_nontrivial   = ac>=5? strtol(argv[4],0,0):0; // Argument 4: Only generate fullerenes with nontrivial symmetry group?
    int n_best_candidates = ac>=6? strtol(argv[5],0,0):100; // Argument 5: How many best fullerne candidates do you want to store? 

    // Make sure output directory exists
    mkdir(output_dir.c_str(),0777);
    int n_fullerenes = IsomerDB::number_isomers(N,"Any",IPR); 
    int Nd = 1;
    auto batch_size = min(isomerspace_forcefield::optimal_batch_size(N,0)*16, n_fullerenes/Nd);

    IsomerBatch<GPU> B0s[Nd] = {IsomerBatch<GPU>(N,batch_size,0)};
    std::vector<CuArray<device_real_t>> eccentricity(Nd); for (int i = 0; i < Nd; i++) eccentricity[i] = CuArray<device_real_t>(batch_size);
    std::vector<CuArray<device_real_t>> volumes(Nd); for (int i = 0; i < Nd; i++) volumes[i] = CuArray<device_real_t>(batch_size);
    device_io::IsomerQueue Q0s[Nd] = {device_io::IsomerQueue(N,0)}; for (int i = 0; i < Nd; i++) Q0s[i].resize(batch_size*4);
    std::vector<LaunchCtx> gen_ctxs(Nd); for (int i = 0; i < Nd; i++) gen_ctxs[i] = LaunchCtx(i);
    auto policy = LaunchPolicy::ASYNC;
    BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,IPR,only_nontrivial);
    Graph G;
    auto Nf = N/2 + 2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
    G.N = Nf;
    int I=0;			// Global isomer number at start of batch
    auto generate_isomers = [&](int M){
    bool more_isomers = true;
    if(I == n_fullerenes) return false;
    for (int i = 0; i < M; i++){
        if (more_isomers){
            //auto ID = cuda_benchmark::random_isomer("isomerspace_samples/dual_layout_"+to_string(N)+"_seed_42", G);
            more_isomers = BuckyGen::next_fullerene(BuckyQ, G);
            if (!more_isomers) break;
            Q0s[I%Nd].insert(G,I, gen_ctxs[I%Nd], policy);
            I++;
        }
    }
    gen_ctxs[0].wait();
    return true;
    };
    IsomerBatch<CPU> start_batch(N,batch_size,0);

    //Start of program
    generate_isomers(batch_size*4);
    for (int i = 0; i < Nd; i++) Q0s[i].refill_batch(B0s[i]);
    for (int i = 0; i < Nd; i++) gen_ctxs[i].wait();
    for (int i = 0; i < Nd; i++) {
        isomerspace_dual::dualize(B0s[i], gen_ctxs[i], policy);
        isomerspace_tutte::tutte_layout(B0s[i],  50*N,gen_ctxs[i], policy);
        isomerspace_X0::zero_order_geometry(B0s[i], 4.0, gen_ctxs[i], policy);
    }
    for (int i = 0; i < Nd; i++){
        gen_ctxs[i].wait();
        device_io::copy(start_batch, B0s[i]);
    }
    for (int i = 0; i < Nd; i++) {
        isomerspace_forcefield::optimize<PEDERSEN>(B0s[i], 5*N, 5*N, gen_ctxs[i], policy);
        //isomerspace_properties::transform_coordinates(B0s[i], gen_ctxs[i], policy);
        //isomerspace_properties::eccentricities(B0s[i], eccentricity[i], gen_ctxs[i], policy);
        //isomerspace_properties::volume_divergences(B0s[i], volumes[i], gen_ctxs[i], policy);
    }
    for (int i = 0; i < Nd; i++) gen_ctxs[i].wait();
    IsomerBatch<CPU> host_batch(N,batch_size);
    device_io::copy(host_batch, B0s[0]);

    //std::vector<int> nan_positions_0;
    //for (int i = 0; i < eccentricity[0].size(); i++) {
    //    if (std::isnan(eccentricity[0][i])) {
    //        nan_positions_0.push_back(i);
    //    }
    //}

    std::vector<int> failed_positions_0;
    for (int i = 0; i < host_batch.isomer_capacity; i++) {
        if (host_batch.statuses[i] == IsomerStatus::FAILED) {
            failed_positions_0.push_back(i);
        }
    }

    for (int i = 0; i  < failed_positions_0.size(); i++){
        int ID = host_batch.IDs[failed_positions_0[i]];
        Polyhedron P = host_batch.get_isomer(failed_positions_0[i]).value();
        Polyhedron Pstart = start_batch.get_isomer(failed_positions_0[i]).value();

        Polyhedron::to_file(P, "FailedGeometry_"+to_string(ID) + ".mol2");
        Polyhedron::to_file(Pstart, "StartGeometry_"+to_string(ID) + ".mol2");
        Polyhedron Pref(Pstart);
        Pref.optimize();
        Polyhedron::to_file(Pref, "RefGeometry_"+to_string(ID) + ".mol2");
    }
    //std::cout << "Nan positions: " << nan_positions_0 << std::endl;
    std::cout << "Failed positions: " << failed_positions_0 << std::endl;


}
