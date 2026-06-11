#pragma once
#include <limits.h>
#include <sys/types.h>
#include <type_traits>
#include "fullerenes/graph.hh"


namespace BuckyGen {
  
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
    bool next_fullerene(Graph& G);

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

  bool next_fullerene(buckygen_queue& Q, Graph& G);

  // Backward-compatibility shim: fill any owned graph/triangulation type
  // (Triangulation, FullereneDual, ...) directly. buckygen emits a dual
  // triangulation; we fill a Graph, then deep-copy its adjacency into `out`.
  // In the view-based hierarchy these owned types no longer derive from Graph,
  // so the Graph& overload above cannot bind them.
  template<class G_t>
    requires (!std::is_same_v<std::remove_cv_t<G_t>, Graph>)
  bool next_fullerene(buckygen_queue& Q, G_t& out) {
    Graph G;
    bool r = next_fullerene(Q, G);
    if (r) out = G;
    return r;
  }
}


