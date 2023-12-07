#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/triangulation.hh>
#include <fullerenes/symmetry.hh>
#include <fullerenes/isomerdb.hh>
#include <random>
#include <chrono>
using namespace std::chrono;

vector<size_t> m_samples(size_t M, size_t m, int seed=42)
{
  vector<size_t> full_range(M), sample_indices(m);
  for(size_t i=0;i<M;i++) full_range[i] = i;
  std::mt19937 rng(seed);

  shuffle(full_range.begin(),full_range.end(),rng);
  for(size_t i=0;i<m;i++) sample_indices[i] = full_range[i];
  sort(sample_indices.begin(), sample_indices.end());

  return sample_indices;
}

int main(int ac, char **av)
{
  size_t default_M = 1000000, default_m = 10000;
  if(ac<2) {
    fprintf(stderr,"Syntax: %s <N> [M:%ld] [m:%ld] [IPR:0] [only_symmetric:0]\n",
	    av[0], default_M, default_m);
    return -1;
  }
  bool
    IPR = ac>4? strtol(av[4],0,0):0,
    only_symmetric = ac>5? strtol(av[5],0,0):0;

  int     N = strtol(av[1],0,0);  
  int64_t Nisomers = IsomerDB::number_isomers(N,only_symmetric?"Nontrivial":"Any",IPR);
  size_t
    M = min<size_t>(ac>2? strtol(av[2],0,0):default_M,Nisomers),
    m = min<size_t>(ac>3? strtol(av[3],0,0):default_m,Nisomers);
  
  printf("Processing %ld samples out of %ld from the C%d isomerspace (IPR=%d, only_symmetric=%d, Nisomers=%ld).\n",
	 m,M,N,IPR,only_symmetric,Nisomers);
  

  // State for buckygen-timing:
 vector<size_t> sample_ix = m_samples(M,m);
 size_t ii=0;
 steady_clock::time_point before_time, after_time;
 using msec = duration<double,milliseconds::period>;
 msec buckygen_time{0};

 // Processer de grafer, som tilhoerer vores sample-liste
 auto process_graph = [&](const Triangulation &g, size_t cnt){
   // Her kan man smide en masse andet processering ind
 }; 
  
  BuckyGen::buckygen_queue BQ = BuckyGen::start(N,IPR,only_symmetric);
  Triangulation g;
  
  size_t cnt = 0;
  after_time = steady_clock::now();
  while(BuckyGen::next_fullerene(BQ,g) && cnt < M){
    before_time = steady_clock::now();

    if(cnt == sample_ix[ii]){
      buckygen_time += before_time-after_time; // Tid som det har taget at vente at buckygen har genereret g
      process_graph(g, cnt);
      ii++;
    }
    cnt++;
    after_time = steady_clock::now();
  }
  
  printf("Processed %ld out of %ld C%d-graphs in %gms (%gus/graph, %gus/vertex)\n",
	 m,cnt,N,
	 buckygen_time.count(),1000*buckygen_time.count()/m,1000*buckygen_time.count()/(m*N));

  
  return 0;
}
