#pragma once
#include <limits.h>
#include <sys/types.h>
#include "fullerenes/graph.hh"


namespace BuckyGen {
  
  typedef struct {
    pid_t pid;
    int   qid;
    int   Nvertices;
    int   IPR;
    size_t chunk_index, chunk_number;
  } buckygen_queue;

  buckygen_queue start(int N, int IPR, bool only_nontrivial = false, 
		       size_t chunk_index=0, size_t chunk_number=1);
  void stop(const buckygen_queue& Q);

  bool next_fullerene(const buckygen_queue& Q, Graph& G);
}


