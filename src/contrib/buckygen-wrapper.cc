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
  Q.msg_qbytes = queue_qbytes(Q.qid);  // inherited by the child via QGlobal

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
    // Non-blocking send + poll (per chunk): a blocking msgsnd can't notice the
    // parent dying while it waits on a full queue. On EAGAIN (queue full ->
    // consumer slow OR parent gone) check getppid() against the captured parent
    // and bail if orphaned, else back off briefly. Drained path sends first try.
    for(;;){
      if(msgsnd(Q.qid, (void*)&msg, msgsz, IPC_NOWAIT) >= 0) break;
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
  return true;
}

// The worker's payload is the dual's adjacency in the library's storage
// layout: rows of BUCKYGEN_DMAX with -1 padding, in cyclic planar order
// (buckygen's e->next traversal preserves it).  A dual's stride is
// BUCKYGEN_DMAX, so a matching destination takes each row wholesale and only
// the degrees are counted; a wider one is padded.
//
// dst may view any memory the caller owns, device (shared-USM) slots included.
static void fill_dual(const int* flat, node_t Nv, FullereneDualView dst)
{
  if(dst.N != Nv)
    throw std::runtime_error("BuckyGen::next_fullerene: destination holds "
                             + std::to_string((long)dst.N) + " vertices, the dual has "
                             + std::to_string((long)Nv));
  if(dst.dmax < BUCKYGEN_DMAX)
    throw std::runtime_error("BuckyGen::next_fullerene: destination stride "
                             + std::to_string(dst.dmax) + " cannot hold degree "
                             + std::to_string(BUCKYGEN_DMAX));
  if((long)dst.neighbours.size() < (long)Nv*dst.dmax || (long)dst.deg.size() < Nv)
    throw std::runtime_error("BuckyGen::next_fullerene: destination spans are "
                             "smaller than its own N x dmax");

  for(node_t u=0; u<Nv; u++){
    const int*  src = flat + (size_t)u*BUCKYGEN_DMAX;
    node_t*     row = dst.neighbours.data() + (size_t)u*dst.dmax;
    uint8_t     d   = 0;
    while(d < BUCKYGEN_DMAX && src[d] != -1){ row[d] = src[d]; d++; }
    for(int i=d; i<dst.dmax; i++) row[i] = -1;   // pad a wider destination
    dst.deg[u] = d;
  }
}

bool next_fullerene(buckygen_queue& Q, FullereneDualView dst)
{
  // Single worker: a graph's chunks arrive contiguously and in seq order, so
  // we just append them until the graph is complete.
  int flat[BUCKYGEN_DMAX*MAXN];
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
        fill_dual(flat, Q.Nvertices, dst);
        return true;
      }
    } else throw std::runtime_error("BuckyGen::next_fullerene: unexpected IPC message type "
                                    + std::to_string(msg.mtype));
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


  bool buckyherd_queue::next_fullerene(FullereneDualView dst)
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
	  fill_dual(H.reasm_buf[w].data(), (node_t)H.Nvertices, dst);
	  H.reasm_fill[w] = 0; H.reasm_got[w] = 0; H.reasm_nchunks[w] = -1;
	  return true;
	}
      } else if(msg.mtype == WORKER_FINISHED) {	 // A worker finished!
	H.active_workers--;
	int worker_index = msg.worker_index;
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
