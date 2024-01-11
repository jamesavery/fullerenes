#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <unistd.h>
#include <signal.h>
#include <assert.h>
#include <unistd.h>

#include "fullerenes/graph.hh"
#include "fullerenes/auxiliary.hh"
#include "fullerenes/buckygen-wrapper.hh"

namespace BuckyGen {


buckygen_queue QGlobal;	// Only access this from child process
bool push_graph(const buckygen_queue& Q);

#define MAIN buckygen_main
#define DONT_OUTPUT 1
#define FILTER(nbtot,nbop,doflip) push_graph(QGlobal)

extern "C" {
  // TODO: Find pragma to turn all the compiler warnings caused by the buckygen code.
  // NB: Currently must be run in sub-process, as buckygen uses a gazillion global variables.
#include "buckygen-1.0/buckygen.c"
}

enum { INVALID,GRAPH_READY,WORKER_FINISHED } msg_type;
  
void signal_finished(const buckygen_queue& Q)
{
  struct { long mtype; int worker_index[1];  } msg;
  msg.mtype = WORKER_FINISHED;
  msg.worker_index[0] = Q.worker_index;

  msgsnd(Q.qid,(void*)&msg,sizeof(msg),0);
}


void stop(const buckygen_queue& Q)
{
  pid_t gid = getpid();
  sighandler_t old_handler = signal(SIGTERM,SIG_IGN); // Protect ourselves while we kill our children
  killpg(gid,SIGTERM);
  signal(SIGTERM,old_handler);                        // Restore normalcy.
  msgctl(Q.qid,IPC_RMID,0);			      // Kill the Sys-V IPC queue
}

  
buckygen_queue start(int N, bool IPR, bool only_nontrivial, 
		     size_t chunk_index, size_t chunk_number)
{
  buckygen_queue Q;
  Q.qid = msgget(IPC_PRIVATE,IPC_CREAT | 0666);
  Q.Nvertices = N/2+2;
  Q.chunk_index  = chunk_index;
  Q.chunk_number = chunk_number;

  assert(Q.qid >= 0);

  if(!(Q.pid = fork())){	// Child
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
  } else {			// Parent
    pid_t gid = getpid();	// Keep track of children with group ID
    setpgid(Q.pid,gid);	// to kill them all when parent leaves.
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

  msg.mtype = GRAPH_READY;
  memset(msg.neighbours,-1,6*nv*sizeof(int));

  for(int u=0;u<nv;u++){
    EDGE *e(firstedge[u]);
    for(int i=0;i<degree[u];i++,e=e->next)
      msg.neighbours[6*u+i] = e->end;
  }
  
  int snd_result = msgsnd(Q.qid,(void*)&msg,sizeof(long)+6*Q.Nvertices*sizeof(int),0); // Blocking send
  if(snd_result>=0) return true;
  else if(errno == EIDRM){
    //    fprintf(stderr,"BuckyGen::push_graph %d,%d(%s): Queue was removed by parent process - BuckyGen::stop() was called on this queue.\n",snd_result,errno,strerror(errno));
      exit(0);
  } else {
    //    fprintf(stderr,"BuckyGen::push_graph %d,%d(%s): Queue is invalid. Why am I still alive?\n",
    //	    snd_result,errno,
    //	    strerror(errno));
    exit(-1);
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
  } else if(msg.mtype == GRAPH_READY) {	// Completed graph
    G.N = Q.Nvertices;
    G.neighbours.resize(G.N);

    for(int u=0;u<Q.Nvertices;u++){
      G.neighbours[u].clear();
      for(int i=0; 6>i && (msg.neighbours[u*6+i] != -1); i++)
	G.neighbours[u].push_back(msg.neighbours[u*6+i]);
    }

    return true;
  } else if(msg.mtype == WORKER_FINISHED) {	// No more graphs to generate
    stop(Q);
    return false;
  }
  abort();
}


  /************************* HERDING THE BUCKS ***********************/  
  buckygen_queue buckyherd_queue::new_worker(int worker_index) {
    buckygen_queue Q;
    Q.qid          = qid;
    Q.Nvertices    = Nvertices;
    Q.worker_index = worker_index;
    Q.chunk_number = Nchunks;	// TODO: Pick one name.

    // Individual stuff.
    assert(!chunks_todo.empty()); // Don't call on an empty work stack
    Q.chunk_index = chunks_todo.back(); chunks_todo.pop_back();

    if(!(Q.pid = fork())){	// Child process
      QGlobal = Q;		// Different queues for different children: OK due to separate memory spaces.
	
      int  npar = 2;
      char *av[5]  = {strdup("buckygen"), strdup((to_string(N)+"d").c_str())};
      char *ipr    = strdup(IPR? "-I" : "");
      char *chunks = strdup((to_string(Q.chunk_index)+"/"+to_string(Q.chunk_number)).c_str());
  
      if(IPR) av[npar++] = ipr;
      if(Q.chunk_number != 1) av[npar++] = chunks;
      if(only_nontrivial)   av[npar++] = strdup("-V");

      buckygen_main(npar, av);
      signal_finished(Q);
      exit(0);
    } else {			// Parent process
      pid_t gid = getpid();	// Keep track of children with group ID
      setpgid(Q.pid,gid);	// to kill them all when parent leaves.
      assert(Q.pid >= 0);
      active_workers++;
      return Q;
    }
  }
  

  buckyherd_queue::buckyherd_queue(size_t N, size_t Nchunks, size_t Nworkers, 
				   bool IPR, bool only_nontrivial,
				   vector<size_t> _chunks_todo)
    : N(N), Nvertices(N/2+2), Nchunks(Nchunks), Nworkers(Nworkers), active_workers(0),
      IPR(IPR), only_nontrivial(only_nontrivial), chunks_todo(_chunks_todo.begin(), _chunks_todo.end())
  {
    qid = msgget(IPC_PRIVATE,IPC_CREAT|0666); // Create a Sys-V IPC queue
    assert(qid >= 0);

    // Which CPU cores to use?
    free_cpu_cores=0;
    for(int core=1;core<Nworkers+1;core++) free_cpu_cores |= (1<<core);
    // Assign herd leader to CPU core 0, workers to subsequent cores. It's the user's responsibility to assure Nworkers+1 <= number of physical CPU cores
    
    if(chunks_todo.empty())
      for(size_t i=0;i<Nchunks;i++) chunks_todo.push_back(i);
    
    for(size_t i=0;i<Nworkers && !chunks_todo.empty();i++)
      worker_processes.push_back(new_worker(i));
  }
  


  bool buckyherd_queue::next_fullerene(Graph& G) 
  {
    buckyherd_queue &H(*this);
    
    struct {
	long mtype;
	int neighbours[6*MAXN];
    } msg;

    while(true){
      ssize_t length = msgrcv(H.qid, (void*)&msg, sizeof(long)+6*H.Nvertices*sizeof(int), -2, 0);
    
      if(length < 0){
	fprintf(stderr,"In BuckyHerd::next_fullerene: %s\n",strerror(errno));
	return false;
      } else if(msg.mtype == GRAPH_READY) {	// Completed graph
	G.N = H.Nvertices;
	G.neighbours.resize(G.N);
	
	for(int u=0;u<H.Nvertices;u++){
	  G.neighbours[u].clear();
	  for(int i=0; 6>i && (msg.neighbours[u*6+i] != -1); i++)
	    G.neighbours[u].push_back(msg.neighbours[u*6+i]);
	}
	
	return true;
      } else if(msg.mtype == WORKER_FINISHED) {	 // A worker finished!
	H.active_workers--;
	if(H.chunks_todo.empty()){               // If there are no more tasks
	  if(H.active_workers<=0) return false;  // and no active workers, then we're all done!
	} else {
	  int worker_index = msg.neighbours[0];	 // If there are more tasks, give the worker
	  H.worker_processes[worker_index] = H.new_worker(worker_index); // some more work to do.
	  // Since this message didn't contain a graph: don't return, but read the next message.
	}
      }
    }
    abort();
  }

  void buckyherd_queue::stop_all() const {
    pid_t gid = getpid();
    sighandler_t old_handler = signal(SIGTERM,SIG_IGN); // Protect ourselves while we kill our children
    killpg(gid,SIGTERM);
    signal(SIGTERM,old_handler);                        // Restore normalcy.
    msgctl(qid,IPC_RMID,0);			      // Kill the Sys-V IPC queue    
  }

  

  
#undef FILTER
#undef DONT_OUTPUT


}
