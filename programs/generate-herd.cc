#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/triangulation.hh>
#include <fullerenes/symmetry.hh>

#define TEST_BUCKYHERD 1

int main(int ac, char **av)
{
  if(ac<2) {
    fprintf(stderr,"Syntax: %s <N> [Nchunks:1] [Nworkers:1] [IPR:0] [only_nontrivial:0]\n",
	    av[0]);
    return -1;
  }
  int N = strtol(av[1],0,0),
    Nchunks = ac>2?strtol(av[2],0,0):1, Nworkers = ac>3?strtol(av[3],0,0):1,
    IPR = ac>4? strtol(av[4],0,0):0, only_nontrivial = ac>5?strtol(av[5],0,0):0;

  printf("Generating C%d isomerspace (IPR=%d, only_nontrivial=%d) in %d chunks with %d parallel workers.\n",
	 N,IPR,only_nontrivial,Nchunks,Nworkers);
  
#if TEST_BUCKYHERD==0
  BuckyGen::buckygen_queue BQ = BuckyGen::start(N,IPR,only_nontrivial);
  Triangulation g;

  size_t cnt = 0;
  while(BuckyGen::next_fullerene(BQ,g)){
    assert(g.N == N/2+2);
    cnt++;
  }
#else
  // TODO: Add CPU_SET etc. and sched_setaffinity to buckyherd_queue
  BuckyGen::buckyherd_queue HQ(N,Nchunks,Nworkers,IPR,only_nontrivial);
  Triangulation g;

  size_t cnt = 0;
  while(HQ.next_fullerene(g)){
    assert(g.N == N/2+2);
    cnt++;
  }  
#endif  
  printf("Generated %ld C%d graphs\n",cnt,N);
  
  return 0;
}
