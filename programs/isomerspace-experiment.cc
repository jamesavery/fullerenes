//TODO: Nd -> compile time from CMake, assert == device_count()
#include <csignal>
#include <sys/stat.h>
#include <limits.h>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <thread>
#include <numeric>
#include <future>
#include <queue>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/progress_bar.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include <fullerenes/sycl-headers/all-kernels.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>
#include "fullerenes/isomerdb.hh"

#define STORE_HESSIAN 1
#define STORE_EIGENSYSTEM 1


using namespace std;
using namespace std::chrono;


typedef float real_t;	// We'll use the device real type everywhere in this program
typedef uint16_t device_node_t;

  typedef enum {GEN, GEO, OPT, PROP, STAT, NUM_STAGES} stage_t;
  typedef enum {VOLUME,ECCENTRICITY,MIN_FREQ,MAX_FREQ, FREQ_WIDTH, INERTIA, NUM_RESULTS} results_t; // Results on which to do statistics
  string result_names[NUM_RESULTS] = {"volume","eccentricity","minfreq","maxfreq","freqwidth","inertia"};
  string stage_names [NUM_STAGES]  = {"generate","start_geometry","optimized_geometry","properties","statistics"};
  constexpr int result_sizes[NUM_RESULTS] = {1,1,1,1,1,3}; // How many scalars per result?


void signal_callback_handler(int signum)
{
  exit(signum);
}

#include <sys/resource.h>
void ensure_minimal_stacksize(rlim_t minimal_stacksize)
{
  struct rlimit rl;

  int error = getrlimit(RLIMIT_STACK, &rl);
  rl.rlim_cur = std::max(rl.rlim_cur,minimal_stacksize);
  error += setrlimit(RLIMIT_STACK, &rl);

  if(error<0){
    perror("ensure_minimal_stacksize failed: ");
    abort();
  }
}


struct isomer_candidate {
  float value;
  int id;
  array<float,NUM_RESULTS> results;
  
  vector<uint16_t> cubic_neighbours;
  vector<real_t> X;

  isomer_candidate(double value, int id, array<float,NUM_RESULTS> &results, Fullerene<real_t, device_node_t> fullerene):
    value(value), id(id), results(results), cubic_neighbours(3*fullerene.N_),X(3*fullerene.N_)/*, H(90*N), Hcol(90*N)*/ { // Det tager 50% ekstra tid at gemme Hessian'en
    auto X = fullerene.d_.X_cubic_.template as_span<real_t>();
    auto A = fullerene.d_.A_cubic_.template as_span<device_node_t>();
    auto N = fullerene.N_;
    memcpy(&cubic_neighbours[0], A.data(), A.size_bytes());
    memcpy(&X[0],                X.data(), X.size_bytes());
  }

  bool operator<(const isomer_candidate &b) const { return value < b.value; }
  friend ostream &operator<<(ostream &s, const isomer_candidate& c){ s<<"("<<c.value<<","<<c.id<<")"; return s; }
};


template <typename T> struct k_smallest : public std::priority_queue<T> {
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
  const size_t N                = strtol(argv[1],0,0);   // Argument 1: Number of vertices N

  string output_dir     = ac>=3? argv[2] : string("output/C")+to_string(N);     // Argument 2: directory to output files to
  bool IPR               = ac>=4? strtol(argv[3],0,0):0;  // Argument 3: Only generate IPR fullerenes?
  bool only_nontrivial   = ac>=5? strtol(argv[4],0,0):0;  // Argument 4: Only generate fullerenes with nontrivial symmetry group?
  size_t n_best_candidates = ac>=6? strtol(argv[5],0,0):30; // Argument 5: How many best fullerne candidates do you want to store? 

  ensure_minimal_stacksize(N*10000000);
  
  
  // Make sure output directory exists
  mkdir(output_dir.c_str(),0777);
  size_t n_fullerenes = IsomerDB::number_isomers(N,only_nontrivial?"Nontrivial":"Any",IPR);
  n_best_candidates = std::max(size_t(1),std::min(n_best_candidates, n_fullerenes/4));

  
  int    Nd = Device::get_devices(DeviceType::GPU).size(); // Number of devices
  vector<Device> devices = Device::get_devices(DeviceType::GPU);
  size_t batch_size = std::min((size_t)(devices[0].get_property(DeviceProperty::MAX_COMPUTE_UNITS)*10), (size_t)(n_fullerenes/Nd));
  size_t final_batch_size = n_best_candidates*NUM_RESULTS*2/Nd+(Nd-1);

  cout << "Analyzing C"<<N<< " isomer space with " << n_fullerenes << " isomers.\n";
  cout << "n_best_candidates = " << n_best_candidates << "; batch_size = " << batch_size << "; final_batch_size = " << final_batch_size <<"\n";
  // Organize all of our of computation pipeline stages and our scoresl of different results

  array<array<double,2>,NUM_RESULTS> result_bounds = {{
    {0.24*pow(N,3/2.)-50,0.4*pow(N,3/2.)+50}, // Volume bounds
    {1,N/10.},                                // Eccentricity bounds
    {0,10*33.356},                            // Smallest frequency bounds (in cm^{-1} = teraherz*33.356)
    //    {40*33.356,55*33.356},                    // Largest frequency bounds (in cm^{-1})
    {(39+2*log(N-20))*33.356,(43+2*log(N-20))*33.356},                    // Largest frequency bounds (in cm^{-1})    
    {40*33.356,55*33.356},                    // Bandwidth bounds (in cm^{-1})
    {0,0} // What are the inertia bounds?
    }};
  
  array<vector<SyclVector<real_t>>,NUM_RESULTS> results; // TODO: vector<Cuarray> gives linkin error -- needs fixing to allow general Nd.
  array<long double,NUM_RESULTS> means, stdevs;
  array<vector<size_t>,NUM_RESULTS> terrible_outliers;
  vector<SyclVector<real_t>> hessians(Nd);
  vector<SyclVector<device_node_t>> hessian_cols(Nd);
  vector<SyclVector<real_t>> Q(Nd),lams(Nd); 
  vector<SyclVector<real_t>> intermediate_spectrum_ends(Nd);
  array<long double,NUM_RESULTS> Ev, Ev2, K, current_mean, current_stddev;
  size_t n_converged=0;
  for(int r=0; r< NUM_RESULTS; r++) results[r].resize(Nd);
  // Initialize unified memory result CuArray's
  for(int d=0;d<Nd;d++){
    hessians[d]      = SyclVector<real_t>(N*3* 10*3*batch_size,0);
    hessian_cols[d]  = SyclVector<device_node_t>(N*3* 10*3*batch_size,0);
    intermediate_spectrum_ends[d] = SyclVector<real_t>(2*batch_size,0);

    for(int r=0; r< NUM_RESULTS; r++) results[r][d] = SyclVector<real_t>(result_sizes[r]*batch_size,0);
  }

  // We keep track of the n_best_candidates isomers with smallest resp. largest values for each result
  vector<k_smallest<isomer_candidate>>
    result_min(NUM_RESULTS, k_smallest<isomer_candidate>(n_best_candidates)),
    result_max(NUM_RESULTS, k_smallest<isomer_candidate>(n_best_candidates))
    ;
  k_smallest<isomer_candidate> bandwidth_min(n_best_candidates), bandwidth_max(n_best_candidates);

  array<set<pair<real_t,int>>,NUM_RESULTS> result_reference;
  
  // Each stage has one on-device batch per device. TODO: NUM_STAGES vector<vector<IsomerBatch>>, Nd. Dynamic number of devices != 2.
  //std::array<std::vector<FullereneBatch<real_t, uint16_t,2>,NUM_STAGES> Bs;
  std::array<std::vector<FullereneBatch<real_t, device_node_t>>,4> Bs;
  for (int i = 0; i < 4; i++) {
    Bs[i].resize(Nd);
    for (int j = 0; j < Nd; j++) {
      Bs[i][j] = FullereneBatch<real_t, device_node_t>(N,batch_size);
    }
  }
  // Final IsomerBatch on host

  vector<FullereneBatch<real_t, device_node_t>> host_batch(Nd, FullereneBatch<real_t, device_node_t>(N,batch_size));
  vector<FullereneBatch<real_t, device_node_t>> final_batch(Nd, FullereneBatch<real_t, device_node_t>(N,final_batch_size));
  vector<FullereneBatch<real_t, device_node_t>> final_host_batch(Nd, FullereneBatch<real_t, device_node_t>(N,final_batch_size));
  
  // TODO: Organize Qs by stages together with batches.Structure nicely.
  vector<FullereneQueue<real_t, device_node_t>> Q0s(Nd, FullereneQueue(N,batch_size*4)); // Graph-generate to X0-generate
  vector<FullereneQueue<real_t, device_node_t>> Q1s(Nd, FullereneQueue(N,batch_size*4)); // X0-generate to optimization
  vector<FullereneQueue<real_t, device_node_t>> Q2s(Nd, FullereneQueue(N,batch_size*4)); // optimization to properties + stat
  vector<FullereneQueue<real_t, device_node_t>> Q3s(Nd, FullereneQueue(N,batch_size*4)); // k_best result final computations

  DualizeFunctor<real_t, device_node_t> dualize;
  TutteFunctor<real_t, device_node_t> tutte;
  SphericalProjectionFunctor<real_t, device_node_t> spherical_projection;
  ForcefieldOptimizeFunctor<PEDERSEN, real_t, device_node_t> forcefield_optimize;
  HessianFunctor<PEDERSEN, real_t, device_node_t> compute_hessian;
  EigenFunctor<EigensolveMode::ENDS, real_t, device_node_t> spectrum_ends;
  EccentricityFunctor<real_t, device_node_t> eccentricity;
  VolumeFunctor<real_t, device_node_t> volume;
  SurfaceAreaFunctor<real_t, device_node_t> surface_area;
  InertiaFunctor<real_t, device_node_t> inertia;
  TransformCoordinatesFunctor<real_t, device_node_t> transform_coordinates;
  
  
  //for (int i = 0; i < Nd; i++) {Q0s[i].resize(batch_size*4); Q1s[i].resize(batch_size*4); Q2s[i].resize(batch_size*4);}

  // TODO: Should we really have two different streams, or is one enough?
  //std::cout << sizeof(SyclQueue) << std::endl;
  vector<SyclQueue> gen_ctxs; std::generate_n(std::back_inserter(gen_ctxs),Nd,[&,i = 0]() mutable {return SyclQueue(devices[ i++ % devices.size() ],true);});
  vector<SyclQueue> opt_ctxs; std::generate_n(std::back_inserter(opt_ctxs),Nd,[&,i = 0]() mutable {return SyclQueue(devices[ i++ % devices.size() ],true);});

  // TODO: Multi-CPU buckygen
  BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,IPR,only_nontrivial);  
  ProgressBar progress_bar = ProgressBar('#',30);
  
  int I=0;			// Global isomer number at start of batch

  vector<bool> stage_finished(NUM_STAGES,0);
  vector<int>  num_finished  (NUM_STAGES,0);
  vector<vector<int>> num_finished_this_round(NUM_STAGES,vector<int>(Nd));
  
  Graph G;
  auto Nf = N/2 + 2;
  G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
  G.N = Nf;

  constexpr real_t carbon_mass = 1.9944733e-26/*kg*/, aangstrom_length = 1e-10/*m*/;
  
  constexpr size_t num_bins = 1000;
  array<array<size_t,num_bins>,NUM_RESULTS> result_histograms({0});
  vector<array<size_t,num_bins*num_bins>> result_histograms2D(NUM_RESULTS*NUM_RESULTS,{0});  

  // Hack using knowledge about range for easy binning (for now)
  //... histograms
  
  // TODO: Implement proper timer class with tree structure, pass as argument to kernels instead of littering program code
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

  auto policy = LaunchPolicy::SYNC;

  auto generate_isomers = [&](int M){
    bool more_isomers = true;
    if(I == n_fullerenes) return false;
    for (int i = 0; i < M; i++){
        if (more_isomers){
            more_isomers = BuckyGen::next_fullerene(BuckyQ, G);
            if (!more_isomers) break;
            Q0s[I%Nd].push_back(G,I);
            I++;
        }
    }
    num_finished[GEN] = I;
    //TODO: num_finished_this_round[GEN][d]
    stage_finished[GEN] = num_finished[GEN] >= n_fullerenes; // TODO: 
    return true;
  };

  auto geo_routine = [&](){
    std::vector<int> qsize_geo; std::for_each(Q1s.begin(),Q1s.end(),[&](auto &Q){qsize_geo.push_back(Q.size());});
    for (int d = 0; d < Nd; d++){
      if(Q0s[d].size() > 0){
	auto &B   = Bs[GEO][d];
	auto &ctx = gen_ctxs[d];
	real_t scalerad = 4;
  QueueUtil::push(ctx, B, Q0s[d], ConditionFunctor(0, 0, StatusEnum::EMPTY | StatusEnum::CONVERGED_2D | StatusEnum::FAILED_2D) );
	//Q0s[d].refill_batch(B, ctx, policy);
	dualize(ctx, B, policy);
	tutte(ctx, B, policy);
	spherical_projection(ctx, B, policy);
  QueueUtil::push(ctx, Q1s[d], B, ConditionFunctor(0, 0, StatusEnum::CONVERGED_2D | StatusEnum::FAILED_2D | StatusEnum::NOT_CONVERGED) ); //Push all to next stage
	//Q1s[d].push_all(B, ctx, policy);
      }
    }
    for(int d=0; d<Nd; d++){
      gen_ctxs[d].wait();
      num_finished_this_round[GEO][d] = Q1s[d].size() - qsize_geo[d];
      num_finished[GEO]              += num_finished_this_round[GEO][d];
    }
    stage_finished[GEO] = num_finished[GEO] >= n_fullerenes;    // TODO: Debug over-counting
  };
  
  auto opt_routine = [&](bool final=false){
    std::vector<int> qsize_opt; std::for_each(Q2s.begin(),Q2s.end(),[&](auto &Q){qsize_opt.push_back(Q.size());});
    //std::vector<int> qsize_opt = {Q2s[0].get_size(), Q2s[1].get_size()};
    for (int d = 0; d < Nd; d++){
      auto &B   = Bs[OPT][d];
      auto &ctx = opt_ctxs[d];
      //Q1s[d].refill_batch(B, opt_ctxs[d], policy);
      
      QueueUtil::push(ctx, B, Q1s[d], ConditionFunctor(0,0,StatusEnum::CONVERGED_3D | StatusEnum::FAILED_3D | StatusEnum::EMPTY) );
      forcefield_optimize(ctx, B, policy, 3*N, 20*N);      
      QueueUtil::push(ctx, Q2s[d], B, ConditionFunctor(0,0,StatusEnum::CONVERGED_3D | StatusEnum::FAILED_3D), StatusEnum::EMPTY ); //Consume from batch (by setting consumed isomers to EMPTY)
      std::cout << "Q2s[" << d << "].size() = " << Q2s[d].size() << std::endl;
      //Q2s[d].push_done(B, ctx, policy);
      //if(final) device_io::copy(Bs[PROP][d],B,ctx,policy); 
    }
    for (int d = 0; d < Nd; d++){
      opt_ctxs[d].wait();
			
      num_finished_this_round[OPT][d] = Q2s[d].size() - qsize_opt[d];
      num_finished[OPT]              += num_finished_this_round[OPT][d];
    }
    stage_finished[OPT] = num_finished[OPT] >= n_fullerenes; // TODO: Debug over-counting

    if(final){
      for(int d=0;d<Nd;d++){
	size_t NSTAT = int(IsomerStatus::NOT_CONVERGED)+1;
	vector<size_t> n_states(NSTAT);
	
	for(size_t di=0;di<Bs[PROP][d].capacity();di++) n_states[size_t(Bs[PROP][d].m_.flags_[di])]++;

	cout << "Batch " << d << " has statuses: " << n_states << "; num_finished = "<<num_finished[OPT]<<"\n";
      }
    }
    return true;
  };

  auto prop_routine = [&](){
    std::vector<int> qsize(Nd); std::for_each(Q2s.begin(),Q2s.end(),[&](auto &Q){qsize.push_back(Q.size());});
    for (int d = 0; d < Nd; d++){
      auto &B   = Bs[PROP][d];
      auto &ctx = opt_ctxs[d];

      //Q2s[d].refill_batch(B, ctx, policy);
      QueueUtil::push(ctx, B, Q2s[d], ConditionFunctor()); //Empty condition means all indices in the destination batch are valid
      transform_coordinates(ctx, B, policy);
      // Hessian units: kg/s^2, divide eigenvalues by carbon_mass to get 1/s^2, frequency = sqrt(lambda/carbon_mass)
      compute_hessian(ctx, B, policy, hessians[d], hessian_cols[d]);
      SyclVector<real_t> eigenvectorsss;
      spectrum_ends(ctx, B, policy, hessians[d], hessian_cols[d], size_t(40), intermediate_spectrum_ends[d], eigenvectorsss);

      //      isomerspace_eigen     ::lambda_max(B, hessians[d], hessian_cols[d], results[MAX_FREQ][d], 40, ctx, policy);
      // Inertia is in Ångstrøm^2, multiply by carbon_mass and aangstrom_length*aangstrom_length to get SI units kg*m^2
      inertia   (ctx, B, policy, ((Span<real_t>)results[INERTIA][d]).template as_span<std::array<real_t,3>>());      
      eccentricity       (ctx, B, policy, results[ECCENTRICITY][d]);
      volume   (ctx, B, policy, results[VOLUME][d]);
    }
    for (int d = 0; d < Nd; d++){
      opt_ctxs[d].wait();
      num_finished_this_round[PROP][d] = qsize[d] - Q2s[d].size();
      num_finished[PROP]              += num_finished_this_round[PROP][d];
    }
    stage_finished[PROP] = num_finished[PROP] >= n_fullerenes; // TODO: Debug over-counting
  };

  
  auto update_mean_var = [&](real_t v, long double &Ev, long double &Ev2,  double K){
    Ev += v - K;  Ev2 += (v - K)*(v-K);
  };

  auto terrible_outlier = [&](real_t v, int r){
    return v<result_bounds[r][0] || v>result_bounds[r][1];
  };

  auto stat_routine = [&](){
    //    fprintf(stderr,"stat start\n");
    //      vector<real_t> volumes_merged(sum(num_finished_this_round[PROP]));
    //      vector<real_t> eccentricity_merged(sum(num_finished_this_round[PROP]));
    vector<int> opt_failed, opt_not_converged;

    int i=0;
    for(int d=0;d<Nd;d++){
      size_t n = num_finished_this_round[PROP][d];

      // First pass of batch: Update mean and variance for all results.
      for(int di=0;di<n;di++){
	auto &status = Bs[PROP][d].m_.flags_[di];	
	if(status == StatusEnum::CONVERGED_3D){
	  n_converged++;
		  
	  for(int r=0;r<NUM_RESULTS;r++){
	    auto v = results[r][d][di];
	    if(isfinite(v)){
	      if(n_converged==1) K[r] = v;
	      
	      update_mean_var(v,Ev[r],Ev2[r],K[r]);
	    }
	  }	
	}
      }
      for(int r=0;r<NUM_RESULTS;r++){
	current_mean[r]   = K[r] + Ev[r]/n_converged;
	current_stddev[r] = sqrt((Ev2[r] - Ev[r]*Ev[r] / n_converged) / (n_converged - 1.0L));	
      }

      // Now process the results 
      for(int di=0;di<num_finished_this_round[PROP][d];di++,i++){
	int id = Bs[PROP][d].m_.ID_[di];
	auto &status = Bs[PROP][d].m_.flags_[di];
	
	if(status != StatusEnum::EMPTY){
	  num_finished_this_round[STAT][d]++;
	  num_finished[STAT]++;
	}
	  
	if(status == StatusEnum::CONVERGED_3D){
	  array<real_t,NUM_RESULTS> R;

	  // Convert from Hessian eigenvalues to normal mode frequencies in teraherz	    
	  real_t min_freq = sqrt(results[MIN_FREQ][d][di]/carbon_mass)/(2*M_PI)*1e-12, 
	         max_freq = sqrt(results[MAX_FREQ][d][di]/carbon_mass)/(2*M_PI)*1e-12;
	  
	  // Convert frequencies from teraherz to cm^{-1} to compare with experiment
	  min_freq *= 33.356;
	  max_freq *= 33.356;
		 
	  results[MIN_FREQ][d][di]   = min_freq;
	  results[MAX_FREQ][d][di]   = max_freq;
	  results[FREQ_WIDTH][d][di] = max_freq-min_freq;


	  for(int r=0;r<NUM_RESULTS;r++)
	    if(result_sizes[r]==1) R[r] = results[r][d][di];
	  
	  isomer_candidate C(0, id, R, Bs[PROP][d][di]);
	  
	  for(int r=0;r<NUM_RESULTS;r++) if(result_sizes[r]==1){
	      auto v = results[r][d][di];
	      if(terrible_outlier(v,r)){
		v = results[r][d][di] = nan("");
		terrible_outliers[r].push_back(id);
		continue;
	      } 
	      if(isfinite(v)){
		C.value =  v; result_min[r].insert(C);
		C.value = -v; result_max[r].insert(C);
		means[r] += v;       
	      }

	      auto vmin = result_bounds[r][0], vmax = result_bounds[r][1];

	      size_t v_bin = std::min(num_bins-1,std::max(0LU,(size_t)floor(num_bins*((v-vmin)/(vmax-vmin)))));
	      result_histograms[r][v_bin]++;

	      for(int s=0;s<r;s++) if(result_sizes[s]==1){
		  auto w = results[s][d][di];		  
		  if(isfinite(w) && !terrible_outlier(w,s)){
		    auto wmin = result_bounds[s][0], wmax = result_bounds[s][1];
		    size_t w_bin = std::min(num_bins-1,std::max((size_t)0,(size_t)floor(num_bins*((w-wmin)/(wmax-wmin)))));

		    result_histograms2D[(r*(r-1))/2+s][v_bin*num_bins+w_bin]++;
		  }
	      }
	    }

	  
	  //		result_reference[r].insert({v,id+1}); // To check that result_min and result_max actually gets k smallest & largest.
	} else {
	  if(status == StatusEnum::FAILED_3D) opt_failed.push_back(id);
	  if(status == StatusEnum::NOT_CONVERGED) opt_not_converged.push_back(id);	  
	}
      }
    }
  
    cout << "num_finished_this_round = " << num_finished_this_round  << "\n";
    cout << "num_finished            = " << num_finished  << "\n";      
    cout << "opt_failed        = " << opt_failed << "\n"
	 << "opt_not_converged = " << opt_not_converged << "\n";
    
    
    
    if(num_finished[STAT] >= n_fullerenes) stage_finished[STAT] = true; // TODO: We're overcounting somehow, getting num_finished > n_fullerenes. DEBUG!

    //    fprintf(stderr,"stat end\n");    
  };

  auto final_round = [&](int s){
    bool final = !stage_finished[s];
    for(int t=GEN;t<s;t++) final &= stage_finished[t];
    return final;
  };
  
  while(!stage_finished[PROP]){
      cout << loop_iters << " start: Gen: " << num_finished[GEN] << "  Geometry: " << num_finished[GEO] << "  Opt: " << num_finished[OPT] << "  Prop: " << num_finished[PROP] << std::endl;
      cout << loop_iters << " start: Stage sm_.flags_: " << stage_finished[GEN] << ", " << stage_finished[GEO] << ", " << stage_finished[OPT] << ", " << stage_finished[PROP] << std::endl;
      auto T0 = steady_clock::now();

      geo_routine();

      
      auto T1 = steady_clock::now(); Tinit_geom += T1-T0;
      auto handle = std::async(std::launch::async, generate_isomers, 2*batch_size);
      auto T2 = steady_clock::now();
      
      while(Q1s[0].size() > Bs[OPT][0].capacity() && Q1s[1].size() > Bs[OPT][1].capacity()){
          opt_routine();
      }
      if(Q2s[0].size() > Bs[PROP][0].capacity() && Q2s[1].size() > Bs[PROP][1].capacity()){
          prop_routine();
      }
      auto T3 = steady_clock::now(); Topt += T3-T2;
      stat_routine();
      
      handle.wait();
      auto T4 = steady_clock::now(); Tgen += T4-T3;
      ssize_t final_step_max = 100, num_orphans = 0;
      if(final_round(OPT)){
	for(ssize_t final_step=0;final_step < final_step_max  && final_round(OPT); final_step++){
	  fprintf(stderr,"final opt start\n");
	  opt_routine(true);
	  fprintf(stderr,"final opt end\n");
	}
	if(!stage_finished[OPT]){
	  stage_finished[OPT] = true; // Avoid endless loop if miscounting. But also TODO: Find out why miscounting sometimes happens!
	  num_orphans = n_fullerenes - num_finished[OPT];
	  fprintf(stderr,"Finishing final round of optimizing, %ld isomers orphaned.\n",num_orphans);
	}
      }

      if   (final_round(PROP)){ fprintf(stderr,"final prop start\n"); prop_routine(); fprintf(stderr,"final prop end\n"); }
      if   (final_round(STAT)){ fprintf(stderr,"final stat start\n"); stat_routine();  fprintf(stderr,"final stat end\n"); }


      
      if(loop_iters % 3 == 0 || stage_finished[STAT]){
        auto Titer = steady_clock::now() - Tstart;
        Tstart = steady_clock::now();
        //auto Tff   = isomerspace_forcefield::time_spent()/Nd;
        //Tqueue = Topt - Tff;
        //Ttutte = isomerspace_tutte::time_spent()/Nd;
        //TX0    = isomerspace_X0::time_spent()/Nd;
        //Tdual  = isomerspace_dual::time_spent()/Nd;
        //auto TOverhead = Tinit_geom - Tdual - Ttutte - TX0;
        //progress_bar.update_progress((float)num_finished[PROP]/(float)num_fullerenes.find(N)->second, "Pace: " + to_string((Titer/1us) / max((num_finished[PROP] - n0),1)) + " us/isomer       ", {{"Gen          ", Tgen},{"Init Overhead", TOverhead},{"Opt          ", Topt},{"Dual         ", Tdual},{"Tutte        ", Ttutte},{"X0           ", TX0},{"FF Overhead  ", Tqueue}});
        n0 = num_finished[STAT];
      }
      loop_iters++;

      
      cout << loop_iters << " end: Gen: " << num_finished[GEN] << "  Geometry: " << num_finished[GEO] << "  Opt: " << num_finished[OPT] << "  Prop: " << num_finished[PROP] << std::endl;
      cout << loop_iters << " end: Stage statuses: " << stage_finished[GEN] << ", " << stage_finished[GEO] << ", " << stage_finished[OPT] << ", " << stage_finished[PROP] << std::endl;

      auto T5 = steady_clock::now(); Topt += T5-T4;


  }
  cout << "Exited loop, waiting on op_ctxs.\n";
  for(int d=0;d<Nd;d++) opt_ctxs[d].wait();

  array<vector<isomer_candidate>,NUM_RESULTS> smallest_candidates, biggest_candidates;  
  
  for(int r=0;r<NUM_RESULTS;r++){
    means[r] /= num_finished[PROP];
    smallest_candidates[r] = sorted(result_min[r].as_vector());
    biggest_candidates[r]  = sorted(result_max[r].as_vector());
  }

  
  
  cout << "COMPLETED IN " << loop_iters << " rounds.\n";
  cout << "num_finished            = " << num_finished  << "\n";
  for(int r=0;r<NUM_RESULTS;r++) if(result_sizes[r]==1){
      auto &smallest = smallest_candidates[r], &biggest = biggest_candidates[r];      

    cout << result_names[r] << "_min = " << smallest << "\n"
         << "   (removed " << terrible_outliers[r].size() << " outliers, please inspect and fix)\n"
	 << result_names[r] << "_max = " << biggest << "\n\n";
  }
  cout.flush();

  // cout << "all_volumes = " << result_reference[VOLUME] << "\n";
  // cout << "all_eccentricity = " << result_reference[ECCENTRICITY] << "\n";  
  //  cout << "all_maxfreq = " << result_reference[MAX_FREQ] << "\n";
  //  cout << "all_minfreq = " << result_reference[MIN_FREQ] << "\n";


  // GATHER SELECTED CANDIDATES TO RECOMPUTE AND OUTPUT
  auto host_graph = [](const vector<device_node_t>& device_graph){
    const size_t N = device_graph.size()/3;
    Graph G(N,true);
    for(node_t u=0;u<N;u++) G.neighbours[u] = {device_graph[3*u],device_graph[3*u+1],device_graph[3*u+2]};
    return G;
  };

  auto host_points = [](const vector<real_t> &X){
    const size_t N = X.size()/3;
    vector<coord3d> points(N);
    for(node_t u=0;u<N;u++) points[u] = {X[3*u],X[3*u+1],X[3*u+2]};
    return points;
  };

  // Postprocess results on GPU
  {


  
  for(int r=0;r<NUM_RESULTS;r++) if(result_sizes[r]==1){
      // Output candidate isomers
      auto &smallest = smallest_candidates[r], &biggest = biggest_candidates[r];            
      for(int i=0;i<n_best_candidates;i++){
	isomer_candidate C = smallest[i];
	Polyhedron P(host_graph(C.cubic_neighbours), host_points(C.X)); // TODO: Hurtigere!
	Q3s[0].push_back(P,C.id, StatusEnum::CONVERGED_3D); 
      }

      for(int i=0;i<n_best_candidates;i++){
	isomer_candidate C = biggest[i];
	Polyhedron P(host_graph(C.cubic_neighbours), host_points(C.X)); // TODO: Hurtigere!
	Q3s[1].push_back(P,C.id, StatusEnum::CONVERGED_3D);
      }
    }

  if(STORE_HESSIAN){
    for(int d=0;d<Nd;d++){    
    auto &B   = final_batch[d];
    auto &ctx = opt_ctxs[d];
    ctx.wait();
    QueueUtil::push(ctx, B, Q3s[d], ConditionFunctor(0,StatusEnum::CONVERGED_3D) );
    //Q3s[d].refill_batch(B, ctx, policy);
    final_host_batch[d].clear();	

    compute_hessian(ctx, B, policy, hessians[d], hessian_cols[d]);

    if(STORE_EIGENSYSTEM){
      EigenFunctor<EigensolveMode::FULL_SPECTRUM_VECTORS, real_t, device_node_t> eigen;
      Q[d]    = SyclVector<real_t>(3*N*3*N*final_batch_size,0);
      lams[d] = SyclVector<real_t>(3*N*final_batch_size,0);    
      eigen(ctx, B, policy,hessians[d],hessian_cols[d],0,lams[d], Q[d]);
    }
    final_host_batch[d] = B; //device_io::copy(final_host_batch[d],B,ctx,policy);
    ctx.wait();
    }}
  
  }
  
  // GEM RESULTATER
  mkdir(output_dir.c_str(), 0777);
  FILE *f = fopen((output_dir+"/result_bounds.float64").c_str(),"wb");
  fwrite(&result_bounds[0],sizeof(double),2*NUM_RESULTS,f);
  fclose(f);

  f = fopen((output_dir+"/result_names.string").c_str(),"w");
  for(int r=0;r<NUM_RESULTS;r++) fprintf(f,"%s\n",result_names[r].c_str());
  fclose(f);
  
  FILE *f_min    = fopen((output_dir+"/result_mins.float32").c_str(),"wb");
  FILE *f_max    = fopen((output_dir+"/result_maxs.float32").c_str(),"wb");
  FILE *f_hist   = fopen((output_dir+"/result_hists.uint64").c_str(),"wb");
  FILE *f_hist2D = fopen((output_dir+"/result_hists2D.uint64").c_str(),"wb");

  // TODO: Automatiser min/max isf. at al koden er dobbelt
  // TODO: Faktoriser ud i output_routine, så man kan lave små batch-størrelser uden at crashe (Q3s[d].is_empty loop)
  for(int r=0;r<NUM_RESULTS;r++) if(result_sizes[r]==1) {
    string result_dir = output_dir+"/"+result_names[r];
    string min_dir    = result_dir + "/min", max_dir = result_dir+"/max", outlier_dir = result_dir + "/terrible_outliers";
    auto &smallest = smallest_candidates[r], &biggest = biggest_candidates[r];          
    

    
    // Output candidate isomers
    auto result_index = [&](bool big, size_t r, size_t i) 
    {
      return array<size_t,2>{{big,r*n_best_candidates+i}};
    };
    
    //TODO: Redundant due to writing out all results for k_smallest
    vector<real_t> smallest_values(n_best_candidates), biggest_values(n_best_candidates);
    for(size_t i=0;i<n_best_candidates;i++){
      smallest_values[i] = smallest[i].value;
      biggest_values [i] = biggest[i].value;
    }
    fwrite(&smallest_values[0],sizeof(real_t),n_best_candidates,f_min);
    fwrite(&biggest_values[0], sizeof(real_t),n_best_candidates,f_max);    
    
    mkdir(result_dir.c_str(), 0777);
    mkdir(min_dir.c_str(), 0777);

    mkdir(max_dir.c_str(), 0777);
    mkdir(outlier_dir.c_str(), 0777);    
    
    for(size_t i=0;i<n_best_candidates;i++){      
      assert(i<smallest.size());
      assert(i<biggest.size());

      string min_basename = min_dir+"/Pmin"+to_string(i),
	     max_basename = max_dir+"/Pmax"+to_string(i);
      
      isomer_candidate Cmin = smallest[i], Cmax = biggest[i];
      auto [d_min,di_min] = result_index(0,r,i);
      auto [d_max,di_max] = result_index(1,r,i);      

      auto Pmin = (Polyhedron)final_host_batch[d_min][di_min];//.get_isomer(di_min).value();
      auto Pmax = (Polyhedron)final_host_batch[d_max][di_max];//.get_isomer(di_max).value();

      Polyhedron::to_file(Pmin,min_basename+".mol2");
      Polyhedron::to_file(Pmin,min_basename+".spiral");      
      Polyhedron::to_file(Pmax,max_basename+".mol2");
      Polyhedron::to_file(Pmax,max_basename+".spiral");

      {// TODO: binary_write til funktion, som ogsaa tjekker for fejl + perror()'er	
	FILE *f = fopen((min_basename+"-graph.uint16").c_str(),"wb");
	//	fwrite(&final_host_batch[d_min].cubic_neighbours+di_min*3*N, sizeof(device_node_t), 3*N,f);
	fwrite(&Cmin.cubic_neighbours[0], sizeof(device_node_t), 3*N,f);
	fclose(f);

	f = fopen((max_basename+"-graph.uint16").c_str(),"wb");
	fwrite(&Cmax.cubic_neighbours[0], sizeof(device_node_t), 3*N,f);
	fclose(f);

	f = fopen((min_basename+"-PX.float32").c_str(),"wb");
	fwrite(&Cmin.X[0], sizeof(real_t), 3*N,f);
	fclose(f);

	f = fopen((max_basename+"-PX.float32").c_str(),"wb");
	fwrite(&Cmax.X[0], sizeof(real_t), 3*N,f);
	fclose(f);	

	f = fopen((min_basename+"-X.float32").c_str(),"wb");
	fwrite(final_host_batch[d_min].d_.X_cubic_.data()+di_min*3*N, sizeof(real_t), 3*N,f);
	fclose(f);

	
	f = fopen((max_basename+"-X.float32").c_str(),"wb");
	fwrite(final_host_batch[d_max].d_.X_cubic_.data()+di_max*3*N, sizeof(real_t), 3*N,f);
	fclose(f);

	
	if(STORE_HESSIAN){
	f = fopen((min_basename+"-hessians.float32").c_str(),"wb");
	fwrite(&hessians[d_min].data()[di_max*3*3*10*N], sizeof(real_t), 3*3*10*N,f); 
	fclose(f);

	f = fopen((max_basename+"-hessians.float32").c_str(),"wb");
	fwrite(&hessians[d_max].data()[di_max*3*3*10*N], sizeof(real_t), 3*3*10*N,f); 
	fclose(f);	

	f = fopen((min_basename+"-hessian_cols.uint16").c_str(),"wb");
	fwrite(&hessian_cols[d_min].data()[di_max*3*3*10*N], sizeof(uint16_t), 3*3*10*N,f); 
	fclose(f);

	f = fopen((max_basename+"-hessian_cols.uint16").c_str(),"wb");
	fwrite(&hessian_cols[d_max].data()[di_max*3*3*10*N], sizeof(uint16_t), 3*3*10*N,f); 
	fclose(f);	

	if(STORE_EIGENSYSTEM){
	  f = fopen((min_basename+"-hessian_Q.float32").c_str(),"wb");
	  fwrite(&Q[d_min].data()[di_max*3*N*3*N], sizeof(real_t), 3*N*3*N,f); 
	  fclose(f);

	  f = fopen((min_basename+"-hessian_lams.float32").c_str(),"wb");
	  fwrite(&lams[d_min].data()[di_max*3*N], sizeof(real_t), 3*N,f); 
	  fclose(f);			

	  f = fopen((max_basename+"-hessian_Q.float32").c_str(),"wb");
	  fwrite(&Q[d_max].data()[di_max*3*N*3*N], sizeof(real_t), 3*N*3*N,f); 
	  fclose(f);

	  f = fopen((max_basename+"-hessian_lams.float32").c_str(),"wb");
	  fwrite(&lams[d_max].data()[di_max*3*N], sizeof(real_t), 3*N,f); 
	  fclose(f);			
	}
	}
	
	f = fopen((min_basename+"-results.float32").c_str(),"wb");
	fwrite(&Cmin.results[0],sizeof(real_t),NUM_RESULTS-1,f); 
	fclose(f);

	f = fopen((max_basename+"-results.float32").c_str(),"wb");
	fwrite(&Cmax.results[0],sizeof(real_t),NUM_RESULTS-1,f); 
	fclose(f);
      }
    }

    // Output histograms:
    fwrite(&result_histograms[r][0],sizeof(uint64_t),num_bins,f_hist);

    for(int s=0;s<r;s++){
      fwrite(&result_histograms2D[(r*(r-1))/2+s][0],sizeof(uint64_t),num_bins*num_bins,f_hist2D);
    }
  }

  fclose(f_min);
  fclose(f_max);
  fclose(f_hist);
  fclose(f_hist2D);

  

  // REPORTING TIMINGS
  auto Ttot = steady_clock::now() - T0;

  
  
  auto TgeomOverhead = Tinit_geom - Tdual - Ttutte - TX0;
  auto Tsum = Tgen + Tupdate + Tdual + Ttutte + TX0 + Tcopy + Topt + Tfile + Toutq + Tinq + TgeomOverhead;
  //cout << std::endl;

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
