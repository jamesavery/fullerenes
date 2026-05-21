#pragma once
#include <limits.h>
#include <sys/types.h>
#include "fullerenes/graph.hh"


namespace BuckyGen {
  
  typedef struct {
    pid_t pid{};
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

    void stop_all() const;
    
    buckyherd_queue(size_t N, size_t Nchunks, size_t Nworkers,
		    bool IPR, bool only_nontrivial,
		    vector<size_t> chunks_todo = {});

    ~buckyherd_queue(){ stop_all(); }
  };   
  
  buckygen_queue start(int N, bool IPR, bool only_nontrivial = false, 
		       size_t chunk_index=0, size_t chunk_number=1);

  
  void stop(const buckygen_queue& Q);

  bool next_fullerene(const buckygen_queue& Q, Graph& G);
}


