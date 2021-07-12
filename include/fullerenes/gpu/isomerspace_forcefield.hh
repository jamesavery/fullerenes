#pragma once
#include <inttypes.h>

namespace IsomerspaceForcefield {
  typedef float    device_real_t;
  typedef uint16_t device_node_t;
  
  size_t computeBatchSize(size_t N);
  
  void   OptimizeBatch(device_real_t* h_X,
		       device_node_t* h_cubic_neighbours,
		       device_node_t* h_next_on_face,
		       device_node_t* h_prev_on_face,
		       uint8_t* h_face_right,
		       const size_t N,
		       const size_t batch_size);

};
