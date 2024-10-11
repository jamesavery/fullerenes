#pragma once
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>
#include <fullerenes/buckygen-wrapper.hh>

template <typename T, typename K>
void fill(FullereneBatch<T,K>& B, int mytask_id = 0, int ntasks = 1) {
  int N = B.N_;
  int Nf = B.Nf_;
  int N_graphs = B.capacity();
  auto face_degrees_acc = B.d_.deg_;
  auto dual_neighbours_acc = B.d_.A_dual_;
  auto statuses_acc = B.m_.flags_;
  BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,false, false, mytask_id, ntasks);
  Graph G;

  G.neighbours = neighbours_t(Nf, std::vector<node_t>(6,-1));
  G.N = Nf;
  int num_generated = 0;
  for (int i = 0; i < N_graphs; ++i) {
    bool more_isomers = BuckyGen::next_fullerene(BuckyQ, G);
    if (!more_isomers) break;
    num_generated++;
    statuses_acc[i] = StatusFlag::DUAL_INITIALIZED;
    for(int j = 0; j < Nf; j++) {
      face_degrees_acc[i*Nf + j] = G.neighbours[j].size();
      for(int k = 0; k < G.neighbours[j].size(); k++) 
        dual_neighbours_acc[i*Nf + j][k] = G.neighbours[j][k];

    }
  }
  BuckyGen::stop(BuckyQ);
  if (num_generated < N_graphs) 
    for (int i = num_generated; i < N_graphs; ++i) {
      statuses_acc[i] = StatusFlag::DUAL_INITIALIZED;
      //Repeat the same graphs as already generated.
      for(int j = 0; j < Nf; j++) {
        face_degrees_acc[i*Nf + j] = face_degrees_acc[(i%num_generated)*Nf + j];
        for(int k = 0; k < 6; k++) 
          dual_neighbours_acc[i*Nf + j][k] = dual_neighbours_acc[(i%num_generated)*Nf + j][k];

      }
    }
}