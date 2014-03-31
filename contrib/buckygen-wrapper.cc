#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <signal.h>
#include <assert.h>
#include <unistd.h>

#include "../libgraph/graph.hh"
#include "../libgraph/auxiliary.hh"
#include "buckygen-wrapper.hh"

namespace BuckyGen {


buckygen_queue QGlobal;	// Only access this from child process
bool push_graph(const buckygen_queue& Q);

#define MAIN buckygen_main
#define DONT_OUTPUT 1
#define FILTER(nbtot,nbop,doflip) push_graph(QGlobal)

extern "C" {
  // TODO: Find pragma to turn all the compiler warnings caused by the buckygen code.
#include "buckygen-1.0/buckygen.c"
}

void signal_finished(const buckygen_queue& Q)
{
  struct { long mtype; char text[1]; } msg;
  msg.mtype = 2;
  msgsnd(Q.qid,(void*)&msg,sizeof(msg),0);
}


void stop(const buckygen_queue& Q)
{
  kill(Q.pid,9);
  msgctl(Q.qid,IPC_RMID,0);
}

  
buckygen_queue start(int N, int IPR, bool only_nontrivial, 
		     size_t chunk_index, size_t chunk_number)
{
  buckygen_queue Q;
  Q.qid = msgget(IPC_PRIVATE,IPC_CREAT | 0666);
  Q.Nvertices = N/2+2;
  Q.IPR     = IPR;
  Q.chunk_index  = chunk_index;
  Q.chunk_number = chunk_number;

  assert(Q.qid >= 0);

  if(!(Q.pid = fork())){	// Child?
    QGlobal = Q;
    
    int  npar = 2;
    char *av[5]  = {strdup("buckygen"), strdup((to_string(N)+"d").c_str())};
    char *ipr    = strdup(IPR? "-I" : "");
    char *chunks = strdup((to_string(chunk_index)+"/"+to_string(chunk_number)).c_str());
    
    if(IPR) av[npar++] = ipr;
    if(chunk_number != 1) av[npar++] = chunks;
    if(only_nontrivial)   av[npar++] = strdup("-V");

    buckygen_main(npar, av);
    signal_finished(Q);
    exit(0);
  } else {			// Parent?
    assert(Q.pid >= 0);
    return Q;
  }
}

bool push_graph(const buckygen_queue& Q)
{
  // Uses global buckygen variables
  int nv = Q.Nvertices;
  struct {
    long mtype;
    int neighbours[6*MAXN];
  } msg;

  msg.mtype = 1;
  memset(msg.neighbours,-1,6*nv*sizeof(int));

  for(int u=0;u<nv;u++){
    EDGE *e(firstedge[u]);
    for(int i=0;i<degree[u];i++,e=e->next)
      msg.neighbours[6*u+i] = e->end;
  }
  
  ssize_t length = msgsnd(Q.qid,(void*)&msg,sizeof(long)+6*Q.Nvertices*sizeof(int),0); // Blocking send
  if(length>=0) return true;
  else {
    fprintf(stderr,"In BuckyGen::push_graph: %s\n",strerror(errno));
    return false;
  }
}

bool next_fullerene(const buckygen_queue& Q, Graph& G)
{
  struct {
    long mtype;
    int neighbours[6*MAXN];
  } msg;
  ssize_t length = msgrcv(Q.qid, (void*)&msg, sizeof(long)+6*Q.Nvertices*sizeof(int), -2, 0);

  if(length < 0){
    fprintf(stderr,"In BuckyGen::next_fullerene: %s\n",strerror(errno));
    return false;
  } else if(msg.mtype == 1) {	// Completed graph
    G.N = Q.Nvertices;
    G.neighbours.resize(G.N);

    for(int u=0;u<Q.Nvertices;u++){
      G.neighbours[u].clear();
      for(int i=0; 6>i && (msg.neighbours[u*6+i] != -1); i++)
	G.neighbours[u].push_back(msg.neighbours[u*6+i]);
    }

    return true;
  } else if(msg.mtype == 2) {	// No more graphs to generate
    //    stop(Q);
    return false;
  }
  abort();
}

#undef FILTER
#undef DONT_OUTPUT


}
