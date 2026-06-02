#pragma once
#include <limits.h>
#include <sys/types.h>
#include <type_traits>
#include "fullerenes/graph.hh"


namespace BuckyGen {
  
  typedef struct {
    pid_t pid{};
    int   qid{};
    int   Nvertices{};
    int   chunk_index{}, worker_index{}, chunk_number{};
  } buckygen_queue;

  struct buckyherd_queue {
    int qid{};
    size_t N{}, Nvertices{}, Nchunks{}, Nworkers{}, active_workers{};
    bool IPR{}, only_nontrivial{};
    vector<size_t> chunks_todo{}, chunks_done{};
    vector<buckygen_queue>  worker_processes{};
    size_t free_cpu_cores{};

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

  // Backward-compatibility shim: fill any owned graph/triangulation type
  // (Triangulation, FullereneDual, ...) directly. buckygen emits a dual
  // triangulation; we fill a Graph, then deep-copy its adjacency into `out`.
  // In the view-based hierarchy these owned types no longer derive from Graph,
  // so the Graph& overload above cannot bind them.
  template<class G_t>
    requires (!std::is_same_v<std::remove_cv_t<G_t>, Graph>)
  bool next_fullerene(const buckygen_queue& Q, G_t& out) {
    Graph G;
    bool r = next_fullerene(Q, G);
    if (r) out = G;
    return r;
  }
}


