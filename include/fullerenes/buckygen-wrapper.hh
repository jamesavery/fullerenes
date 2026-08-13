#pragma once
#include <limits.h>
#include <sys/types.h>
#include <type_traits>
#include "fullerenes/graph.hh"
#include "fullerenes/owned.hh"


namespace BuckyGen {

  // buckygen's payload stride: it emits a fullerene dual, whose maximum degree
  // is 6, as rows of this width with -1 padding.  Taken from the view type so
  // the enumerator and the storage layout state the width once.
  inline constexpr int BUCKYGEN_DMAX = FullereneDualView::default_dmax;

  // Size an owning graph/dual type to hold one buckygen dual -- Nv vertices at
  // BUCKYGEN_DMAX, the payload's own width, whatever the destination type's
  // default stride is -- and return the view the enumerator fills.
  //
  // Reuse is the normal case (one destination filled in a loop): storage is
  // reallocated only when the shape changes, and any twin table is dropped
  // because it describes the previous graph.
  template<class G_t> requires owning_graph<G_t>
  FullereneDualView dual_slot(G_t& out, node_t Nv) {
    if((node_t)out.N != Nv || out.dmax != BUCKYGEN_DMAX ||
       out.owned_neighbours.size() != (size_t)Nv*BUCKYGEN_DMAX){
      out.owned_neighbours.assign((size_t)Nv*BUCKYGEN_DMAX, node_t(-1));
      out.owned_deg.assign((size_t)Nv, 0);
      out.N    = Nv;
      out.dmax = BUCKYGEN_DMAX;
    }
    if constexpr (requires { out.owned_twin; }) out.owned_twin.clear();
    out.repoint();
    return FullereneDualView(out.N, out.dmax, out.neighbours, out.deg);
  }

  typedef struct {
    pid_t pid{};
    pid_t parent_pid{};        // captured pre-fork; the child polls getppid() against it
    int   qid{};
    int   Nvertices{};
    int   chunk_index{}, worker_index{}, chunk_number{};
    long  msg_qbytes{};  // per-queue byte cap (IPC_STAT msg_qbytes); graphs are
                         // chunked to fit it. macOS defaults to 2048, far below
                         // a single large-fullerene message; Linux to 16384.
  } buckygen_queue;

  struct buckyherd_queue {
    int qid{};
    size_t N{}, Nvertices{}, Nchunks{}, Nworkers{}, active_workers{};
    bool IPR{}, only_nontrivial{};
    vector<size_t> chunks_todo{}, chunks_done{};
    vector<buckygen_queue>  worker_processes{};
    size_t free_cpu_cores{};
    long   msg_qbytes{};  // per-queue byte cap; see buckygen_queue::msg_qbytes.

    // Per-worker reassembly buffers: graphs arrive split across messages that
    // interleave across workers on the shared queue, but a single worker's
    // chunks stay FIFO-ordered, so each worker reassembles independently.
    vector<vector<int>> reasm_buf{};
    vector<int>         reasm_fill{}, reasm_got{}, reasm_nchunks{};

    buckygen_queue new_worker(int worker_index);

    // Fill a caller-provided dual view (see the free next_fullerene below).
    bool next_fullerene(FullereneDualView dst);

    // ... or any owning graph/dual type, sized here to the dual's own stride.
    template<class G_t> requires owning_graph<G_t>
    bool next_fullerene(G_t& out) {
      return next_fullerene(dual_slot(out, (node_t)Nvertices));
    }

    void stop_all();           // idempotent: reaps + clears worker_processes

    buckyherd_queue(size_t N, size_t Nchunks, size_t Nworkers,
		    bool IPR, bool only_nontrivial,
		    vector<size_t> chunks_todo = {});

    // Owns forked process groups + an IPC queue; copying would double-teardown
    // (each copy's dtor would killpg/reap the same -- possibly reused -- pids).
    buckyherd_queue(const buckyherd_queue&) = delete;
    buckyherd_queue& operator=(const buckyherd_queue&) = delete;

    ~buckyherd_queue(){ stop_all(); }
  };
  
  buckygen_queue start(int N, bool IPR, bool only_nontrivial = false, 
		       size_t chunk_index=0, size_t chunk_number=1);

  
  // Non-const Q: stop() nulls Q.pid/Q.qid after teardown so it is idempotent
  // (a second stop -- e.g. the caller's, after next_fullerene already stopped on
  // WORKER_FINISHED -- becomes a no-op rather than killpg'ing a reaped/reused pid).
  void stop(buckygen_queue& Q);

  // Fill a caller-provided dual view -- the body every spelling reaches.  The
  // caller owns the memory; it may be anything addressable, a device
  // (shared-USM) slot included, so a batch can have each dual written straight
  // into the arena its kernels read.
  bool next_fullerene(buckygen_queue& Q, FullereneDualView dst);

  // ... or an owning graph/dual type (Graph, Triangulation, FullereneDual,
  // ...), sized here and then filled through the same body.
  template<class G_t> requires owning_graph<G_t>
  bool next_fullerene(buckygen_queue& Q, G_t& out) {
    return next_fullerene(Q, dual_slot(out, (node_t)Q.Nvertices));
  }
}


