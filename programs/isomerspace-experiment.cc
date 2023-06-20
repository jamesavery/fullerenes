#include <csignal>
#include <sys/stat.h>
#include <limits.h>
#include <cmath>
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

#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/isomerdb.hh"

using namespace std;
using namespace std::chrono;
using namespace cuda_io;

using namespace gpu_kernels;

void signal_callback_handler(int signum)
{
  exit(signum);
}

struct candidate {
  double error;
  vector<int> rspi;

  bool operator<(const candidate &b) const { return error < b.error; }
};


template <typename T> struct k_smallest : public priority_queue<T> {
public:
  size_t k;
  
  k_smallest(size_t k) : k(k) {}

  bool insert(const T& x){
    if(this->size() < k || x < this->top()){
      if(this->size() == k) this->pop();
      this->push(x);
      return true;
    }
    return false;
  }

  vector<T>& as_vector() {
    return (*this).*&k_smallest::c;
  }
  
};



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

  int    Nd = LaunchCtx::get_device_count();
  size_t batch_size = min(isomerspace_forcefield::optimal_batch_size(N,0)*16, n_fullerenes/Nd);

  enum {GEN, GEO, OPT, PROP, NUM_STAGES} stage;  

  // Each stage has one on-device batch per device. TODO: Dynamic number of devices != 2.
  IsomerBatch Bs[4][2] = {
    {IsomerBatch(N,batch_size,DEVICE_BUFFER,0),IsomerBatch(N,batch_size,DEVICE_BUFFER,1)},
    {IsomerBatch(N,batch_size,DEVICE_BUFFER,0),IsomerBatch(N,batch_size,DEVICE_BUFFER,1)},
    {IsomerBatch(N,batch_size,DEVICE_BUFFER,0),IsomerBatch(N,batch_size,DEVICE_BUFFER,1)},
    {IsomerBatch(N,batch_size,DEVICE_BUFFER,0),IsomerBatch(N,batch_size,DEVICE_BUFFER,1)}    
  };
  
  // Final IsomerBatch on host
  IsomerBatch H2s[Nd] = {IsomerBatch(N,batch_size,HOST_BUFFER,0), IsomerBatch(N, batch_size, HOST_BUFFER,1)};
  IsomerBatch host_opt_batch[Nd] = {IsomerBatch(N,batch_size,HOST_BUFFER,0), IsomerBatch(N, batch_size, HOST_BUFFER,1)};
    IsomerBatch host_geo_batch[Nd] = {IsomerBatch(N,batch_size,HOST_BUFFER,0), IsomerBatch(N, batch_size, HOST_BUFFER,1)};
    IsomerBatch host_prop_batch[Nd] = {IsomerBatch(N,batch_size,HOST_BUFFER,0), IsomerBatch(N, batch_size, HOST_BUFFER,1)};  
  
  //CuArray<device_real_t> Qs[Nd] = {CuArray<device_real_t>(N*3*N*3*batch_size,0), CuArray<device_real_t>(N*3*N*3*batch_size,0)};
  //CuArray<device_real_t> Hs[Nd] = {CuArray<device_real_t>(N*3* 10*3*batch_size,0), CuArray<device_real_t>(N*3* 10*3*batch_size,0)};
  CuArray<device_node_t> Cols[Nd] = {CuArray<device_node_t>(N*3*3*10 * batch_size,0), CuArray<device_node_t>(N*3*3*10 * batch_size,0)};
  CuArray<device_real_t> Eigs[Nd] = {CuArray<device_real_t>(N*3 * batch_size,0), CuArray<device_real_t>(N*3 * batch_size,0)};

  vector<CuArray<device_real_t>> volume_results(Nd);       for(int i=0;i<Nd;i++) volume_results[i]       = CuArray<device_real_t>(batch_size,0);
  vector<CuArray<device_real_t>> eccentricity_results(Nd); for(int i=0;i<Nd;i++) eccentricity_results[i] = CuArray<device_real_t>(batch_size,0);

  IsomerQueue Q0s[Nd] = {cuda_io::IsomerQueue(N,0), cuda_io::IsomerQueue(N,1)};
  IsomerQueue Q1s[Nd] = {cuda_io::IsomerQueue(N,0), cuda_io::IsomerQueue(N,1)};
  IsomerQueue Q2s[Nd] = {cuda_io::IsomerQueue(N,0), cuda_io::IsomerQueue(N,1)};

  for (int i = 0; i < Nd; i++) {Q0s[i].resize(batch_size*4); Q1s[i].resize(batch_size*4); Q2s[i].resize(batch_size*4);}

  LaunchCtx gen_ctxs[Nd] = {LaunchCtx(0),LaunchCtx(1)};
  LaunchCtx opt_ctxs[Nd] = {LaunchCtx(0),LaunchCtx(1)};

  BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,IPR,only_nontrivial);  
  ProgressBar progress_bar = ProgressBar('#',30);
  
  int I=0;			// Global isomer number at start of batch

  array<bool, NUM_STAGES> stage_finished = {false,false,false,false};
  array<int,  NUM_STAGES> num_finished = {0,0,0,0};
  vector<vector<int>> num_finished_this_round({vector<int>(Nd),vector<int>(Nd),vector<int>(Nd),vector<int>(Nd)});

  cout << "num_finished_this_round initialized to " << num_finished_this_round[2] << "\n";
  
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

  auto loop_iters = 0;
  auto Tstart = steady_clock::now();
  auto n0 = 0;

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
    //TODO: num_finished_this_round[GEN][d]
    stage_finished[GEN] = num_finished[GEN] == n_fullerenes;
    return true;
  };

  auto geo_routine = [&](){
    std::vector<int> qsize_geo = {Q1s[0].get_size(), Q1s[1].get_size()};
    for (int d = 0; d < Nd; d++){
      if(Q0s[d].get_size() > 0){
	Q0s[d].refill_batch(Bs[GEO][d], gen_ctxs[d], policy);
	isomerspace_dual::dualise(Bs[GEO][d], gen_ctxs[d], policy);
	isomerspace_tutte::tutte_layout(Bs[GEO][d], N*10, gen_ctxs[d], policy);
	isomerspace_X0::zero_order_geometry(Bs[GEO][d], 4.0f, gen_ctxs[d], policy);
	Q1s[d].push_all(Bs[GEO][d], gen_ctxs[d], policy);
	cuda_io::copy(host_geo_batch[d], Bs[GEO][d], gen_ctxs[d], policy);	
      }
    }
    for(int d=0; d<Nd; d++){
      gen_ctxs[d].wait();
      num_finished_this_round[GEO][d] = Q1s[d].get_size() - qsize_geo[d];
      num_finished[GEO]              += num_finished_this_round[GEO][d];
    }
    stage_finished[GEO] = num_finished[GEO] == n_fullerenes;    
  };
  
  auto opt_routine = [&](){
    std::vector<int> qsize = {Q2s[0].get_size(), Q2s[1].get_size()};
    for (int d = 0; d < Nd; d++){
      Q1s[d].refill_batch(Bs[OPT][d], opt_ctxs[d], policy);
      isomerspace_forcefield::optimise<PEDERSEN>(Bs[OPT][d], ceil(0.5*N), 5*N, opt_ctxs[d], policy);
      cuda_io::copy(host_opt_batch[d], Bs[OPT][d], opt_ctxs[d], policy);      
      Q2s[d].push_done(Bs[OPT][d], opt_ctxs[d], policy);
    }
    for (int d = 0; d < Nd; d++){
      opt_ctxs[d].wait();
      num_finished_this_round[OPT][d] = Q2s[d].get_size() - qsize[d];
      num_finished[OPT]              += num_finished_this_round[OPT][d];

      // if(num_finished[OPT]>0){
      // 	int d0 = 0, i0 = 18;
      // 	Polyhedron P1 = host_opt_batch[d0].get_isomer(i0).value();
      // 	Polyhedron::to_file(P1,string("isomer_X-step")+to_string(loop_iters)+"-"
      // 			    +to_string(d0)+"_"+to_string(i0)+".mol2");
      // }
    }
    stage_finished[OPT] = num_finished[OPT] == n_fullerenes;
  };

  auto prop_routine = [&](){
    std::vector<int> qsize = {Q2s[0].get_size(), Q2s[1].get_size()};
    for (int d = 0; d < Nd; d++){
      Q2s[d].refill_batch(Bs[PROP][d], opt_ctxs[d], policy);
      //      isomerspace_properties::transform_coordinates(Bs[PROP][d], opt_ctxs[d], policy);
      //      isomerspace_properties::eccentricities    (Bs[PROP][d], eccentricity_results[d], opt_ctxs[d], policy);
      //isomerspace_properties::volume_divergences(Bs[PROP][d], volume_results[d], opt_ctxs[d], policy);
      cuda_io::copy(host_prop_batch[d], Bs[PROP][d], opt_ctxs[d], policy);      
      Bs[PROP][d].clear(opt_ctxs[d],policy);      
    }
    for (int d = 0; d < Nd; d++){
      opt_ctxs[d].wait();
      num_finished_this_round[PROP][d] = qsize[d] - Q2s[d].get_size();
      num_finished[PROP]              += num_finished_this_round[PROP][d];
    }
    stage_finished[PROP] = num_finished[PROP] == n_fullerenes;
  };

  while(!stage_finished[PROP] && loop_iters < 1000){
      std::cout << "Gen: " << num_finished[GEN] << "  Geometry: " << num_finished[GEO] << "  Opt: " << num_finished[OPT] << "  Prop: " << num_finished[PROP] << std::endl;
      std::cout << "Stage statuses: " << stage_finished[GEN] << ", " << stage_finished[GEO] << ", " << stage_finished[OPT] << ", " << stage_finished[PROP] << std::endl;
      auto T0 = steady_clock::now();

      geo_routine();

      auto T1 = steady_clock::now(); Tinit_geom += T1-T0;
      auto handle = std::async(std::launch::async, generate_isomers, 2*batch_size);
      auto T2 = steady_clock::now();
      
      while(Q1s[0].get_size() > Bs[GEO][0].capacity() && Q1s[1].get_size() > Bs[GEO][1].capacity()){
          opt_routine();
      }
      while(Q2s[0].get_size() > Bs[OPT][0].capacity() && Q2s[1].get_size() > Bs[OPT][1].capacity()){
          prop_routine();
      }
      auto T3 = steady_clock::now(); Topt += T3-T2;
      handle.wait();
      auto T4 = steady_clock::now(); Tgen += T4-T3;
      while(stage_finished[GEN] && stage_finished[GEO] && !stage_finished[OPT]){
          opt_routine();
      }

      while(stage_finished[GEN] && stage_finished[GEO] && stage_finished[OPT] && !stage_finished[PROP]){
          prop_routine();
      }

      auto T5 = steady_clock::now(); Topt += T5-T4;
      // Do stuff on CPU. TODO: Nicer Nd
      vector<device_real_t> volumes_merged(sum(num_finished_this_round[PROP]));
      vector<device_real_t> eccentricity_merged(sum(num_finished_this_round[PROP]));
      vector<pair<int,int>> vol_nanids, ecc_nanids, opt_failed;
    
      int i=0;
      for(int d=0;d<Nd;d++)
	for(int di=0;di<num_finished_this_round[PROP][d];di++,i++){
	  auto v = volumes_merged[i]      = volume_results[d][di];
	  auto e = eccentricity_merged[i] = eccentricity_results[d][di];
	  
	  if(!isfinite(v)) vol_nanids.push_back({d,di});
	  if(!isfinite(e) || di==0){
	    ecc_nanids.push_back({d,di});
	    Polyhedron P0 = host_geo_batch[d].get_isomer(di).value();
	    Polyhedron::to_file(P0,string("failed_X0-")+to_string(d)
				+"_"+to_string(di)+".mol2");

	    Polyhedron P1 = host_opt_batch[d].get_isomer(di).value();
	    Polyhedron P1ref = P0;
	    P1ref.optimise();

	    Polyhedron::to_file(P1,string("failed_X-")+to_string(d)+"_"+to_string(di)+".mol2");
	    Polyhedron::to_file(P1ref,string("failed_Xref-")+to_string(d)+"_"+to_string(di)+".mol2");

	    Polyhedron P2 = host_prop_batch[d].get_isomer(di).value();
	    Polyhedron::to_file(P2,string("failed_X2-")+to_string(d)+"_"+to_string(di)+".mol2");	    
	    }
	  if(host_opt_batch[d].statuses[di] == IsomerStatus::FAILED){
	    opt_failed.push_back({d,di});
	    //	    uint16_t *g = host_geo_batch[d].cubic_neighbours;//+(di*3*N);
	    
	    //	    cout << vector<uint16_t>(g,g+100*3*N)<<"\n";
	  }
	}

      cout
	<< "opt_failed = " << opt_failed << ";\n"
	<< "vol_nanids = " << vol_nanids << ";\n"
	<< "ecc_nanids = " << ecc_nanids << ";\n\n";
      
      
      
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
