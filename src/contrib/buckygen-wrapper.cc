#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <unistd.h>
#include <signal.h>
#include <assert.h>
#include <unistd.h>
#include <signal.h>

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

// A graph is shipped worker->parent as one or more chunks of this message.
// mtext (everything after mtype) = {worker_index, seq, nchunks, ndata, data...}.
// A chunk's mtext size is (4 + ndata)*sizeof(int) and must not exceed the
// queue's msg_qbytes, or msgsnd fails with EINVAL (macOS) / blocks (Linux).
struct BGMsg {
  long mtype;
  int  worker_index;            // sender, for per-worker reassembly in the herd
  int  seq;                     // 0-based chunk index within this graph
  int  nchunks;                 // total chunks for this graph
  int  ndata;                   // number of ints carried in data[] this chunk
  int  data[6*MAXN];            // slice of the flattened (6 per vertex) neighbour array
};
static const int BGMSG_HEADER_INTS = 4;  // worker_index, seq, nchunks, ndata

// Read the queue's byte cap. Default-fallback to the macOS-safe 2048 if the
// stat fails for any reason -- never trust an un-probed limit.
static long queue_qbytes(int qid)
{
  struct msqid_ds ds;
  if(msgctl(qid, IPC_STAT, &ds) == 0 && ds.msg_qbytes > 0) return (long)ds.msg_qbytes;
  return 2048;
}

// Largest number of payload ints that fit in one message on this queue.
// Leave 8 bytes of slack for any per-message accounting the kernel folds in.
static int chunk_capacity_ints(long msg_qbytes)
{
  long usable = msg_qbytes - 8 - BGMSG_HEADER_INTS*(long)sizeof(int);
  long n = usable / (long)sizeof(int);
  return n < 1 ? 1 : (int)n;
}

void signal_finished(const buckygen_queue& Q)
{
  BGMsg msg;
  msg.mtype = WORKER_FINISHED;
  msg.worker_index = Q.worker_index;
  // mtext = just the worker_index int (its offset is 0 within mtext).
  msgsnd(Q.qid, (void*)&msg, sizeof(int), 0);
}


void stop(const buckygen_queue& Q)
{
  pid_t gid = getpid();
  auto old_handler = signal(SIGTERM,SIG_IGN); // Protect ourselves while we kill our children
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
  Q.msg_qbytes = queue_qbytes(Q.qid);  // inherited by the child via QGlobal

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
  // Uses global buckygen variables. Flatten the neighbour lists (6 per vertex,
  // -1 padded) then ship them as msg_qbytes-sized chunks. A single msgsnd of
  // the whole graph overflows the queue cap for large fullerenes (macOS: EINVAL,
  // worker dies, parent hangs in msgrcv); chunking keeps every send within cap.
  const int nv    = Q.Nvertices;
  const int total = 6*nv;
  int flat[6*MAXN];
  memset(flat,-1, total*sizeof(int));
  for(int u=0;u<nv;u++){
    EDGE *e(firstedge[u]);
    for(int i=0;i<degree[u];i++,e=e->next)
      flat[6*u+i] = e->end;
  }

  const int cap     = chunk_capacity_ints(Q.msg_qbytes);
  const int nchunks = (total + cap - 1) / cap;
  BGMsg msg;
  msg.mtype        = GRAPH_READY;
  msg.worker_index = Q.worker_index;
  msg.nchunks      = nchunks;

  for(int seq=0; seq<nchunks; seq++){
    const int off = seq*cap;
    const int nd  = (total-off < cap) ? (total-off) : cap;
    msg.seq   = seq;
    msg.ndata = nd;
    memcpy(msg.data, flat+off, nd*sizeof(int));
    const size_t msgsz = (BGMSG_HEADER_INTS + nd)*sizeof(int);
    if(msgsnd(Q.qid,(void*)&msg,msgsz,0) < 0){      // Blocking send (backpressure ok)
      if(errno == EIDRM) exit(0);                   // parent called BuckyGen::stop()
      else               exit(-1);                  // queue invalid -- should not happen now
    }
  }
  return true;
}

bool next_fullerene(const buckygen_queue& Q, Graph& G)
{
  // Single worker: a graph's chunks arrive contiguously and in seq order, so
  // we just append them until the graph is complete.
  int flat[6*MAXN];
  int fill = 0, got = 0, nchunks = -1;
  BGMsg msg;

  while(true){
    ssize_t length = msgrcv(Q.qid, (void*)&msg, sizeof(msg)-sizeof(long), -2, 0);
    if(length < 0){
      fprintf(stderr,"In BuckyGen::next_fullerene: %s\n",strerror(errno));
      return false;
    } else if(msg.mtype == WORKER_FINISHED) {	// No more graphs to generate
      stop(Q);
      return false;
    } else if(msg.mtype == GRAPH_READY) {
      if(nchunks < 0) nchunks = msg.nchunks;
      memcpy(flat+fill, msg.data, msg.ndata*sizeof(int));
      fill += msg.ndata;
      if(++got == nchunks){			// Completed graph
        Graph adj(Q.Nvertices, GRAPH_DMAX);
        for(int u=0;u<Q.Nvertices;u++)
          for(int i=0; 6>i && (flat[u*6+i] != -1); i++)
            adj.push_back(u, flat[u*6+i]);
        // Buckygen's e->next traversal preserves cyclic planar order.
        G = Graph(adj);
        return true;
      }
    } else abort();
  }
}


  /************************* HERDING THE BUCKS ***********************/  
  buckygen_queue buckyherd_queue::new_worker(int worker_index) {
    buckygen_queue Q;
    Q.qid          = qid;
    Q.Nvertices    = Nvertices;
    Q.worker_index = worker_index;
    Q.chunk_number = Nchunks;	// TODO: Pick one name.
    Q.msg_qbytes   = msg_qbytes; // inherited by the child via QGlobal

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
    msg_qbytes = queue_qbytes(qid);

    // Per-worker graph reassembly state (one slot per worker).
    reasm_buf.assign(Nworkers, vector<int>(6*Nvertices, -1));
    reasm_fill.assign(Nworkers, 0);
    reasm_got.assign(Nworkers, 0);
    reasm_nchunks.assign(Nworkers, -1);

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
    BGMsg msg;

    while(true){
      ssize_t length = msgrcv(H.qid, (void*)&msg, sizeof(msg)-sizeof(long), -2, 0);

      if(length < 0){
	fprintf(stderr,"In BuckyHerd::next_fullerene: %s\n",strerror(errno));
	return false;
      } else if(msg.mtype == GRAPH_READY) {	// One chunk of a graph
	const int w = msg.worker_index;
	if(H.reasm_got[w] == 0) H.reasm_nchunks[w] = msg.nchunks;
	memcpy(H.reasm_buf[w].data()+H.reasm_fill[w], msg.data, msg.ndata*sizeof(int));
	H.reasm_fill[w] += msg.ndata;
	if(++H.reasm_got[w] == H.reasm_nchunks[w]){	// This worker completed a graph
	  const int* flat = H.reasm_buf[w].data();
	  Graph adj(H.Nvertices, GRAPH_DMAX);
	  for(int u=0;u<H.Nvertices;u++)
	    for(int i=0; 6>i && (flat[u*6+i] != -1); i++)
	      adj.push_back(u, flat[u*6+i]);
	  // Buckygen's e->next traversal preserves cyclic planar order.
	  G = Graph(adj);
	  H.reasm_fill[w] = 0; H.reasm_got[w] = 0; H.reasm_nchunks[w] = -1;
	  return true;
	}
      } else if(msg.mtype == WORKER_FINISHED) {	 // A worker finished!
	H.active_workers--;
	if(H.chunks_todo.empty()){               // If there are no more tasks
	  if(H.active_workers<=0) return false;  // and no active workers, then we're all done!
	} else {
	  int worker_index = msg.worker_index;	 // If there are more tasks, give the worker
	  H.worker_processes[worker_index] = H.new_worker(worker_index); // some more work to do.
	  // Since this message didn't contain a graph: don't return, but read the next message.
	}
      }
    }
    abort();
  }

  void buckyherd_queue::stop_all() const {
    pid_t gid = getpid();
    auto old_handler = signal(SIGTERM,SIG_IGN); // Protect ourselves while we kill our children
    killpg(gid,SIGTERM);
    signal(SIGTERM,old_handler);                        // Restore normalcy.
    msgctl(qid,IPC_RMID,0);			      // Kill the Sys-V IPC queue    
  }

  

  
#undef FILTER
#undef DONT_OUTPUT


}
