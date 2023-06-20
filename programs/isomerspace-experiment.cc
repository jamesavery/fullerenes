#include <csignal>
#include <sys/stat.h>
#include <limits.h>
#include <chrono>
#include <filesystem>
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
#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/isomerdb.hh"


using namespace gpu_kernels;

void signal_callback_handler(int signum)
{
  exit(signum);
}


int main(int ac, char **argv)
{
  signal(SIGINT, signal_callback_handler);
  if(ac<2){
    fprintf(stderr,"Syntax: %s <N:int> [output_dir] [IPR:0|1] [only_nontrivial:0|1]\n",argv[0]);
    return -1;
  }
  const size_t N                = strtol(argv[1],0,0);     // Argument 1: Number of vertices N

  string output_dir     = ac>=3? argv[2] : "output";    // Argument 2: directory to output files to
  int IPR               = ac>=4? strtol(argv[3],0,0):0; // Argument 3: Only generate IPR fullerenes?
  int only_nontrivial   = ac>=5? strtol(argv[4],0,0):0; // Argument 4: Only generate fullerenes with nontrivial symmetry group?
  int n_best_candidates = ac>=6? strtol(argv[5],0,0):100; // Argument 5: How many best fullerne candidates do you want to store? 

  // Make sure output directory exists
  mkdir(output_dir.c_str(),0777);
  int n_fullerenes = IsomerDB::number_isomers(N,only_nontrivial?"Nontrivial":"Any",IPR);

  int Nd = LaunchCtx::get_device_count();
  auto batch_size = min(isomerspace_forcefield::optimal_batch_size(N,0)*16, n_fullerenes/Nd);

  IsomerBatch B0s[Nd] = {IsomerBatch(N,batch_size,DEVICE_BUFFER,0), IsomerBatch(N, batch_size, DEVICE_BUFFER,1)};
  IsomerBatch B1s[Nd] = {IsomerBatch(N,batch_size,DEVICE_BUFFER,0), IsomerBatch(N, batch_size, DEVICE_BUFFER,1)};
  IsomerBatch B2s[Nd] = {IsomerBatch(N,batch_size,DEVICE_BUFFER,0), IsomerBatch(N, batch_size, DEVICE_BUFFER,1)};
  IsomerBatch H2s[Nd] = {IsomerBatch(N,batch_size,HOST_BUFFER,0), IsomerBatch(N, batch_size, HOST_BUFFER,1)};
  //CuArray<device_real_t> Qs[Nd] = {CuArray<device_real_t>(N*3*N*3*batch_size,0), CuArray<device_real_t>(N*3*N*3*batch_size,0)};
  //CuArray<device_real_t> Hs[Nd] = {CuArray<device_real_t>(N*3* 10*3*batch_size,0), CuArray<device_real_t>(N*3* 10*3*batch_size,0)};
  CuArray<device_node_t> Cols[Nd] = {CuArray<device_node_t>(N*3*3*10 * batch_size,0), CuArray<device_node_t>(N*3*3*10 * batch_size,0)};
  CuArray<device_real_t> Eigs[Nd] = {CuArray<device_real_t>(N*3 * batch_size,0), CuArray<device_real_t>(N*3 * batch_size,0)};
  std::vector<CuArray<device_real_t>> volume_results(Nd); for (int i = 0; i < Nd; i++) volume_results[i] = CuArray<device_real_t>(batch_size,0);
  std::vector<CuArray<device_real_t>> eccentricity_results(Nd); for (int i = 0; i < Nd; i++) eccentricity_results[i] = CuArray<device_real_t>(batch_size,0);

  

  cuda_io::IsomerQueue Q0s[Nd] = {cuda_io::IsomerQueue(N,0), cuda_io::IsomerQueue(N,1)};
  cuda_io::IsomerQueue Q1s[Nd] = {cuda_io::IsomerQueue(N,0), cuda_io::IsomerQueue(N,1)};
  cuda_io::IsomerQueue Q2s[Nd] = {cuda_io::IsomerQueue(N,0), cuda_io::IsomerQueue(N,1)};

  for (int i = 0; i < Nd; i++) {Q0s[i].resize(batch_size*4); Q1s[i].resize(batch_size*4); Q2s[i].resize(batch_size*4);}

  LaunchCtx gen_ctxs[Nd] = {LaunchCtx(0),LaunchCtx(1)};
  LaunchCtx opt_ctxs[Nd] = {LaunchCtx(0),LaunchCtx(1)};

  BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,IPR,only_nontrivial);  
  ProgressBar progress_bar = ProgressBar('#',30);
  
  
  int I=0;			// Global isomer number at start of batch
  enum {GEN, GEO, OPT, PROP, NUM_STAGES} stage;
  std::array<bool, NUM_STAGES> stage_finished = {false,false,false,false};
  std::array<int, NUM_STAGES> num_finished = {0,0,0,0};
    int num_converged = 0,
        num_failed =0;
  Graph G;
  auto Nf = N/2 + 2;
  G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
  G.N = Nf;

  auto T0 = steady_clock::now();
  auto
    Tgen    = steady_clock::now()-T0,
    Tupdate = steady_clock::now()-T0,
    Tqueue  = steady_clock::now()-T0,
    Tdual   = steady_clock::now()-T0,    
    Ttutte  = steady_clock::now()-T0,
    TX0     = steady_clock::now()-T0,
    Tcopy   = steady_clock::now()-T0,
    Topt    = steady_clock::now()-T0,
    Tfile   = steady_clock::now()-T0,
    Toutq   = steady_clock::now()-T0,
    Tinq    = steady_clock::now()-T0,
    Tinit_geom = steady_clock::now()-T0;



  auto policy = LaunchPolicy::ASYNC;

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
    gen_ctxs[0].wait(); gen_ctxs[1].wait(); 
    num_finished[GEN] = I;
    stage_finished[GEN] = num_finished[GEN] == n_fullerenes;
    return true;
  };

  auto opt_routine = [&](){
    std::vector<int> qsize = {Q2s[0].get_size(), Q2s[1].get_size()};
    for (int i = 0; i < Nd; i++){
      Q1s[i].refill_batch(B1s[i], opt_ctxs[i], policy);
      isomerspace_forcefield::optimise<PEDERSEN>(B1s[i], ceil(0.5*N), 5*N, opt_ctxs[i], policy);
      Q2s[i].push_done(B1s[i], opt_ctxs[i], policy);
    }
    for (int i = 0; i < Nd; i++){
    opt_ctxs[i].wait();
    num_finished[OPT] += Q2s[i].get_size() - qsize[i];
    }
    stage_finished[OPT] = num_finished[OPT] == n_fullerenes;
  };

  auto analyse_routine = [&](){
    std::vector<int> qsize = {Q2s[0].get_size(), Q2s[1].get_size()};
    for (int i = 0; i < Nd; i++){
        Q2s[i].refill_batch(B2s[i], opt_ctxs[i], policy);
        isomerspace_properties::transform_coordinates(B2s[i], opt_ctxs[i], policy);
        isomerspace_properties::eccentricities(B2s[i], eccentricity_results[i], opt_ctxs[i], policy);
        isomerspace_properties::volume_divergences(B2s[i], volume_results[i], opt_ctxs[i], policy);
    }
    for (int i = 0; i < Nd; i++){
        opt_ctxs[i].wait();
        num_finished[PROP] +=  qsize[i] - Q2s[i].get_size();
    }
    stage_finished[PROP] = num_finished[PROP] == n_fullerenes;
  };

  auto loop_iters = 0;
  auto Tstart = steady_clock::now();
  auto n0 = 0;
  while(!stage_finished[PROP] && loop_iters < 100){
      std::cout << "Gen: " << num_finished[GEN] << "  Geometry: " << num_finished[GEO] << "  Opt: " << num_finished[OPT] << "  Prop: " << num_finished[PROP] << std::endl;
      std::cout << "Stage statuses: " << stage_finished[GEN] << ", " << stage_finished[GEO] << ", " << stage_finished[OPT] << ", " << stage_finished[PROP] << std::endl;
      auto T0 = steady_clock::now();
      std::vector<int> qsize_geo = {Q1s[0].get_size(), Q1s[1].get_size()};
      for (int i = 0; i < Nd; i++){
          if(Q0s[i].get_size() > 0){
            Q0s[i].refill_batch(B0s[i], gen_ctxs[i], policy);
            isomerspace_dual::dualise(B0s[i], gen_ctxs[i], policy);
            isomerspace_tutte::tutte_layout(B0s[i], N*10, gen_ctxs[i], policy);
            isomerspace_X0::zero_order_geometry(B0s[i], 4.0f, gen_ctxs[i], policy);
            Q1s[i].push_all(B0s[i], gen_ctxs[i], policy);
          }
      }
      for(int j = 0; j < Nd; j++){
        gen_ctxs[j].wait();
        num_finished[GEO] += Q1s[j].get_size() - qsize_geo[j];
      }
      stage_finished[GEO] = num_finished[GEO] == n_fullerenes;

      auto T1 = steady_clock::now(); Tinit_geom += T1-T0;
      auto handle = std::async(std::launch::async, generate_isomers, 2*batch_size);
      auto T2 = steady_clock::now();
      while(Q1s[0].get_size() > B1s[0].capacity() && Q1s[1].get_size() > B1s[1].capacity()){
          opt_routine();
      }
      while(Q2s[0].get_size() > B2s[0].capacity() && Q2s[1].get_size() > B2s[1].capacity()){
          analyse_routine();
      }
      auto T3 = steady_clock::now(); Topt += T3-T2;
      handle.wait();
      auto T4 = steady_clock::now(); Tgen += T4-T3;
      while(stage_finished[GEN] && stage_finished[GEO] && !stage_finished[OPT]){
          opt_routine();
      }

      while(stage_finished[GEN] && stage_finished[GEO] && stage_finished[OPT] && !stage_finished[PROP]){
          analyse_routine();
      }

      auto T5 = steady_clock::now(); Topt += T5-T4;
      if(loop_iters % 3 == 0 || stage_finished[PROP]){
        auto Titer = steady_clock::now() - Tstart;
        Tstart = steady_clock::now();
        auto Tff   = isomerspace_forcefield::time_spent()/Nd;
        Tqueue = Topt - Tff;
        Ttutte = isomerspace_tutte::time_spent()/Nd;
        TX0    = isomerspace_X0::time_spent()/Nd;
        Tdual  = isomerspace_dual::time_spent()/Nd;
        auto TOverhead = Tinit_geom - Tdual - Ttutte - TX0;
        //progress_bar.update_progress((float)num_finished[PROP]/(float)num_fullerenes.find(N)->second, "Pace: " + to_string((Titer/1us) / max((num_finished[PROP] - n0),1)) + " us/isomer       ", {{"Gen          ", Tgen},{"Init Overhead", TOverhead},{"Opt          ", Topt},{"Dual         ", Tdual},{"Tutte        ", Ttutte},{"X0           ", TX0},{"FF Overhead  ", Tqueue}});
        n0 = num_finished[PROP];
      }
      loop_iters++;
  }
    std::vector<int> volume_sorted_ids_0(volume_results[0].size());
    std::vector<int> volume_sorted_ids_1(volume_results[1].size());
    std::vector<int> eccentricity_sorted_ids_0(eccentricity_results[0].size());
    std::vector<int> eccentricity_sorted_ids_1(eccentricity_results[1].size());
    std::iota(volume_sorted_ids_0.begin(), volume_sorted_ids_0.end(), 0);
    std::iota(volume_sorted_ids_1.begin(), volume_sorted_ids_1.end(), 0);
    std::iota(eccentricity_sorted_ids_0.begin(), eccentricity_sorted_ids_0.end(), 0);
    std::iota(eccentricity_sorted_ids_1.begin(), eccentricity_sorted_ids_1.end(), 0);

    std::sort(volume_sorted_ids_0.begin(), volume_sorted_ids_0.end(), [&](int i1, int i2) {return volume_results[0][i1] < volume_results[0][i2];});
    std::sort(volume_sorted_ids_1.begin(), volume_sorted_ids_1.end(), [&](int i1, int i2) {return volume_results[1][i1] < volume_results[1][i2];});
    std::sort(eccentricity_sorted_ids_0.begin(), eccentricity_sorted_ids_0.end(), [&](int i1, int i2) {return eccentricity_results[0][i1] < eccentricity_results[0][i2];});
    std::sort(eccentricity_sorted_ids_1.begin(), eccentricity_sorted_ids_1.end(), [&](int i1, int i2) {return eccentricity_results[1][i1] < eccentricity_results[1][i2];});

    std::sort(volume_results[0].data, volume_results[0].data + volume_results[0].size());
    std::sort(volume_results[1].data, volume_results[1].data + volume_results[1].size());
    std::sort(eccentricity_results[0].data, eccentricity_results[0].data + eccentricity_results[0].size());
    std::sort(eccentricity_results[1].data, eccentricity_results[1].data + eccentricity_results[1].size());

    std::cout << "Volumes 0: " << volume_results[0] << std::endl;
    std::cout << "Volumes 1: " << volume_results[1] << std::endl;
    std::cout << "Eccentricity 0: " << eccentricity_results[0] << std::endl;
    std::cout << "Eccentricity 1: " << eccentricity_results[1] << std::endl;

    cuda_io::copy(H2s[1], B2s[1]);
    FILE* f = fopen("OMG_SO_SPHERICAL.mol2", "w");
    Polyhedron::to_mol2(H2s[1].get_isomer(eccentricity_sorted_ids_1[0]).value(), f);

    //progress_bar.update_progress((float)num_finished/(float)num_fullerenes.find(N)->second);
    

  //cout << endl << "Finished: " << num_finished << "Failed: , " << num_failed << ", " << num_converged << endl;
  auto Ttot = steady_clock::now() - T0;

  
  
  auto TgeomOverhead = Tinit_geom - Tdual - Ttutte - TX0;
  auto Tsum = Tgen + Tupdate + Tdual + Ttutte + TX0 + Tcopy + Topt + Tfile + Toutq + Tinq + TgeomOverhead;
  //std::cout << std::endl;

  cout << "\n\n\n\n\n\n\n\n";
  /* cout << "Time spent on non:\n"
    "\tTotal Time                     = " << (Ttot/1ms)       << " ms\n"
    "\tTime Unaccounted For           = " << (Ttot-Tsum)/1ms  << " ms\n"
    "\tGeometry Overhead              = " << (TgeomOverhead)/1ms << " ms\n"
    "\tDevice queue                   = " << (Tqueue)/1ms     << " ms\n"
    "\tGenerating graphs              = " << (Tgen/1ms)       << " ms\n"
    "\tDualizing                      = " << (Tdual/1ms)      << " ms\n"
    "\tTutte embedding                = " << (Ttutte/1ms)     << " ms\n"
    "\tInitial geometry               = " << (TX0/1ms)        << " ms\n"
    "\tFF Optimization                = " << (Topt/1ms)       << " ms\n";
   *///

  return 0;
}
