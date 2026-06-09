#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/wait.h>
#include <unistd.h>
#include <signal.h>
#include <time.h>
#include <assert.h>
#include <stdexcept>

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
  struct { long mtype; int worker_index;  } msg;
  memset(&msg,0,sizeof(msg));
  msg.mtype = WORKER_FINISHED;
  msg.worker_index = Q.worker_index;

  msgsnd(Q.qid,
        (void*)&msg,
        sizeof(msg),0);
}

// Reap a child, retrying across EINTR so a delivered signal can't leave a zombie.
static void reap(pid_t pid)
{
  if(pid <= 0) return;
  while(waitpid(pid, nullptr, 0) < 0 && errno == EINTR) { }
}

// The forked buckygen child. Runs the enumerator in-process (buckygen uses a
// gazillion globals, hence a separate process) and never returns.
//
// Post-fork hygiene, in order:
//  - reset SIGTERM to SIG_DFL: the child inherited the host's dispositions, and
//    a host SIGTERM handler would otherwise run here instead of terminating us,
//    hanging stop()'s waitpid (the host is a long-lived process, e.g. a Python
//    interpreter).
//  - become our own process-group leader so stop()/stop_all() can killpg just
//    this queue's subtree, never the caller's group.
//  - bail if already orphaned (getppid() != the captured parent pid). Combined
//    with the getppid poll in push_graph, this is the portable (Linux + macOS)
//    parent-death guard -- no PR_SET_PDEATHSIG, which would also misfire when
//    only the forking thread (not the process) dies.
[[noreturn]] static void run_buckygen_child(const buckygen_queue& Q, int N, bool IPR,
                                            bool only_nontrivial,
                                            size_t chunk_index, size_t chunk_number)
{
  signal(SIGTERM, SIG_DFL);
  setpgid(0, 0);                         // own group (race-free with the parent's setpgid)
  if(getppid() != Q.parent_pid) _exit(0);

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
  _exit(0);                              // _exit, not exit: skip the parent's atexit/stdio/dtors
}


void stop(buckygen_queue& Q)
{
  // The child leads its OWN process group (set in start()), so killpg(Q.pid)
  // tears down only this queue's subtree -- never the caller's group (which,
  // in a long-lived host, would hit sibling queues and unrelated subprocesses).
  // reap() then collects the leader so no zombie lingers. Nulling pid/qid makes
  // stop() idempotent: a second call (the caller's, after next_fullerene already
  // stopped on WORKER_FINISHED) is a no-op instead of killpg'ing a reused pid.
  if(Q.pid > 0)  killpg(Q.pid, SIGTERM);
  if(Q.qid >= 0) msgctl(Q.qid, IPC_RMID, 0);   // remove the queue (also unblocks a child in msgsnd)
  if(Q.pid > 0)  reap(Q.pid);
  Q.pid = -1;
  Q.qid = -1;
}


buckygen_queue start(int N, bool IPR, bool only_nontrivial,
		     size_t chunk_index, size_t chunk_number)
{
  buckygen_queue Q;
  Q.qid          = msgget(IPC_PRIVATE, IPC_CREAT | 0666);
  Q.Nvertices    = N/2+2;
  Q.chunk_index  = chunk_index;
  Q.chunk_number = chunk_number;
  Q.parent_pid   = getpid();

  if(Q.qid < 0) throw std::runtime_error("BuckyGen::start: msgget failed");

  Q.pid = fork();
  if(Q.pid < 0){                              // handled explicitly (survives -DNDEBUG)
    msgctl(Q.qid, IPC_RMID, 0);
    throw std::runtime_error("BuckyGen::start: fork failed");
  }
  if(Q.pid == 0)
    run_buckygen_child(Q, N, IPR, only_nontrivial, chunk_index, chunk_number);

  setpgid(Q.pid, Q.pid);                      // parent half (child's setpgid(0,0) wins the race)
  return Q;
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
  memset(msg.neighbours,-1, 6*MAXN*sizeof(int));

  for(int u=0;u<nv;u++){
    EDGE *e(firstedge[u]);
    for(int i=0;i<degree[u];i++,e=e->next)
      msg.neighbours[6*u+i] = e->end;
  }

  // Non-blocking send + poll: a blocking msgsnd can't notice the parent dying
  // while it waits on a full queue. On EAGAIN (queue full -> consumer slow OR
  // parent gone) we check getppid() against the captured parent and bail if
  // orphaned, else back off briefly. Normal (drained) path sends first try, so
  // this costs nothing when the consumer keeps up.
  const size_t len = sizeof(long) + 6*Q.Nvertices*sizeof(int);
  for(;;){
    if(msgsnd(Q.qid, (void*)&msg, len, IPC_NOWAIT) >= 0) return true;
    if(errno == EINTR) continue;
    if(errno == EIDRM) _exit(0);                       // queue removed by stop()
    if(errno == EAGAIN){
      if(getppid() != Q.parent_pid) _exit(0);          // parent gone -> bail (Linux + macOS)
      struct timespec ts{0, 2*1000*1000};              // 2 ms
      nanosleep(&ts, nullptr);
      continue;
    }
    _exit(-1);                                          // unexpected send error
  }
}

bool next_fullerene(buckygen_queue& Q, Graph& G)
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
    Graph adj(Q.Nvertices, GRAPH_DMAX);
    for(int u=0;u<Q.Nvertices;u++)
      for(int i=0; 6>i && (msg.neighbours[u*6+i] != -1); i++)
	adj.push_back(u, msg.neighbours[u*6+i]);
    // Buckygen's e->next traversal preserves cyclic planar order.
    G = adj;

    return true;
  } else if(msg.mtype == WORKER_FINISHED) {	// No more graphs to generate
    stop(Q);
    return false;
  }
  throw std::runtime_error("BuckyGen::next_fullerene: unexpected IPC message type "
                           + std::to_string(msg.mtype));
}


  /************************* HERDING THE BUCKS ***********************/
  buckygen_queue buckyherd_queue::new_worker(int worker_index) {
    buckygen_queue Q;
    Q.qid          = qid;
    Q.Nvertices    = Nvertices;
    Q.worker_index = worker_index;
    Q.chunk_number = Nchunks;	// TODO: Pick one name.
    Q.parent_pid   = getpid();

    // Individual stuff.
    assert(!chunks_todo.empty()); // Don't call on an empty work stack
    Q.chunk_index = chunks_todo.back(); chunks_todo.pop_back();

    Q.pid = fork();
    if(Q.pid < 0) throw std::runtime_error("BuckyGen herd: fork failed");
    if(Q.pid == 0)
      run_buckygen_child(Q, (int)N, IPR, only_nontrivial, Q.chunk_index, Q.chunk_number);

    setpgid(Q.pid, Q.pid);	// worker leads its OWN group, so stop_all() kills only its subtree
    active_workers++;
    return Q;
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
	Graph adj(H.Nvertices, GRAPH_DMAX);
	for(int u=0;u<H.Nvertices;u++)
	  for(int i=0; 6>i && (msg.neighbours[u*6+i] != -1); i++)
	    adj.push_back(u, msg.neighbours[u*6+i]);
	// Buckygen's e->next traversal preserves cyclic planar order.
	G = Graph(adj);

	return true;
      } else if(msg.mtype == WORKER_FINISHED) {	 // A worker finished!
	H.active_workers--;
	int worker_index = msg.neighbours[0];
	reap(H.worker_processes[worker_index].pid);   // finished worker exited; reap it (no zombie)
	H.worker_processes[worker_index].pid = -1;     // mark the slot empty so stop_all() skips it
	if(H.chunks_todo.empty()){               // If there are no more tasks
	  if(H.active_workers<=0) return false;  // and no active workers, then we're all done!
	} else {
	  H.worker_processes[worker_index] = H.new_worker(worker_index); // give the worker more work
	  // Since this message didn't contain a graph: don't return, but read the next message.
	}
      } else {
	throw std::runtime_error("BuckyHerd::next_fullerene: unexpected IPC message type "
				 + std::to_string(msg.mtype));
      }
    }
  }

  void buckyherd_queue::stop_all() {
    // Each worker leads its own process group (set in new_worker()): kill only
    // those subtrees, never the caller's group, then reap so none zombie.
    // Clearing worker_processes/qid makes this idempotent (the dtor may re-run it).
    for(auto& W : worker_processes) if(W.pid > 0) killpg(W.pid, SIGTERM);
    if(qid >= 0){ msgctl(qid, IPC_RMID, 0); qid = -1; }
    for(auto& W : worker_processes) reap(W.pid);
    worker_processes.clear();
  }




#undef FILTER
#undef DONT_OUTPUT


}
