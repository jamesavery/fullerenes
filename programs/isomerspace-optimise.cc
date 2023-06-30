#include <csignal>
#include <sys/stat.h>
#include <limits.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <thread>
#include <future>
#include "fullerenes/isomerdb.hh"
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
  int n_fullerenes = IsomerDB::number_isomers(N,only_nontrivial?"Nontrivial":"Any", IPR);
  int Nd = LaunchCtx::get_device_count();
  auto batch_size = isomerspace_forcefield::optimal_batch_size(N,0)*16;

  IsomerBatch B0s[Nd] = {IsomerBatch(N,batch_size,DEVICE_BUFFER,0), IsomerBatch(N, batch_size, DEVICE_BUFFER,1)};
  IsomerBatch B1s[Nd] = {IsomerBatch(N,batch_size,DEVICE_BUFFER,0), IsomerBatch(N, batch_size, DEVICE_BUFFER,1)};
  CuArray<device_real_t> Qs[Nd] = {CuArray<device_real_t>(N*3*N*3*batch_size,0), CuArray<device_real_t>(N*3*N*3*batch_size,0)};
  CuArray<device_real_t> Hs[Nd] = {CuArray<device_real_t>(N*3* 10*3*batch_size,0), CuArray<device_real_t>(N*3* 10*3*batch_size,0)};
  CuArray<device_node_t> Cols[Nd] = {CuArray<device_node_t>(N*3*3*10 * batch_size,0), CuArray<device_node_t>(N*3*3*10 * batch_size,0)};
  CuArray<device_real_t> Eigs[Nd] = {CuArray<device_real_t>(N*3 * batch_size,0), CuArray<device_real_t>(N*3 * batch_size,0)};

  

  cuda_io::IsomerQueue Q0s[Nd] = {cuda_io::IsomerQueue(N,0), cuda_io::IsomerQueue(N,1)};
  cuda_io::IsomerQueue Q1s[Nd] = {cuda_io::IsomerQueue(N,0), cuda_io::IsomerQueue(N,1)};
  cuda_io::IsomerQueue Q2s[Nd] = {cuda_io::IsomerQueue(N,0), cuda_io::IsomerQueue(N,1)};

  for (int i = 0; i < Nd; i++) {Q0s[i].resize(batch_size*4); Q1s[i].resize(batch_size*4); Q2s[i].resize(batch_size*4);}

  LaunchCtx gen_ctxs[Nd] = {LaunchCtx(0),LaunchCtx(1)};
  LaunchCtx opt_ctxs[Nd] = {LaunchCtx(0),LaunchCtx(1)};

  BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,IPR,only_nontrivial);  
  ProgressBar progress_bar = ProgressBar('#',30);
  
  
  int I=0;			// Global isomer number at start of batch
  int num_finished = 0,
      num_converged = 0,
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

  auto generate_isomers = [&](int M){
    if(I == n_fullerenes) return false;
    for (int i = 0; i < M; i++){
        if (I < n_fullerenes){
              //auto ID = cuda_benchmark::random_isomer("isomerspace_samples/dual_layout_"+to_string(N)+"_seed_42", G);
              BuckyGen::next_fullerene(BuckyQ, G);
              Q0s[I%Nd].insert(G,I, gen_ctxs[I%Nd], LaunchPolicy::ASYNC);
            I++;
        }
    }
    gen_ctxs[0].wait(); gen_ctxs[1].wait(); 
    return true;
  };

  auto opt_routine = [&](){
    for (int i = 0; i < Nd; i++){
      Q1s[i].refill_batch(B1s[i], opt_ctxs[i], LaunchPolicy::ASYNC);
      isomerspace_forcefield::optimise<PEDERSEN>(B1s[i], ceil(0.5*N), 5*N, opt_ctxs[i], LaunchPolicy::ASYNC);
      Q2s[i].push_done(B1s[i], opt_ctxs[i], LaunchPolicy::ASYNC);}
      for (int i = 0; i < Nd; i++){
        opt_ctxs[i].wait();
        num_finished += Q2s[i].get_size();
        Q2s[i].clear();
      }
  };

  auto loop_iters = 0;
  auto Tstart = steady_clock::now();
  auto n0 = 0;
  while(num_finished < n_fullerenes){
      auto T0 = steady_clock::now();
      for (int i = 0; i < Nd; i++){
          if(Q0s[i].get_size() > 0){
            Q0s[i].refill_batch(B0s[i], gen_ctxs[i], LaunchPolicy::ASYNC);
            isomerspace_dual::dualise(B0s[i], gen_ctxs[i], LaunchPolicy::ASYNC);
            isomerspace_tutte::tutte_layout(B0s[i], N*10, gen_ctxs[i], LaunchPolicy::ASYNC);
            isomerspace_X0::zero_order_geometry(B0s[i], 4.0f, gen_ctxs[i], LaunchPolicy::ASYNC);
            Q1s[i].insert(B0s[i], gen_ctxs[i], LaunchPolicy::ASYNC);
          }
      }
      for(int j = 0; j < Nd; j++){
        gen_ctxs[j].wait();
      }
      auto T1 = steady_clock::now(); Tinit_geom += T1-T0;
      auto handle = std::async(std::launch::async, generate_isomers, 2*batch_size);
      auto T2 = steady_clock::now();
      while(Q1s[0].get_size() > B1s[0].capacity() && Q1s[1].get_size() > B1s[1].capacity()){
          opt_routine();
      }
      auto T3 = steady_clock::now(); Topt += T3-T2;
      handle.wait();
      auto T4 = steady_clock::now(); Tgen += T4-T3;
      while(I == n_fullerenes && num_finished < n_fullerenes && Q0s[0].get_size() == 0 && Q0s[1].get_size() == 0){
          opt_routine();
      }
      auto T5 = steady_clock::now(); Topt += T5-T4;
      if(loop_iters % 3 == 0 || num_finished == n_fullerenes){
        auto Titer = steady_clock::now() - Tstart;
        Tstart = steady_clock::now();
        auto Tff   = isomerspace_forcefield::time_spent()/Nd;
        Tqueue = Topt - Tff;
        Ttutte = isomerspace_tutte::time_spent()/Nd;
        TX0    = isomerspace_X0::time_spent()/Nd;
        Tdual  = isomerspace_dual::time_spent()/Nd;
        auto TOverhead = Tinit_geom - Tdual - Ttutte - TX0;
        progress_bar.update_progress((float)num_finished/(float)num_fullerenes.find(N)->second, "Pace: " + to_string((Titer/1us) / max((num_finished - n0),1)) + " us/isomer       ", {{"Gen          ", Tgen},{"Init Overhead", TOverhead},{"Opt          ", Topt},{"Dual         ", Tdual},{"Tutte        ", Ttutte},{"X0           ", TX0},{"FF Overhead  ", Tqueue}});
        n0 = num_finished;
      }
      loop_iters++;
  }

  
    
  
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
