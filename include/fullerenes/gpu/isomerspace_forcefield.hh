#pragma once
#include <inttypes.h>
#include <string>
#include <fullerenes/fullerenegraph.hh> 

#define GPU_REAL float       
#define GPU_REAL3 float3
#define Block_Size_Pow_2 256
#define SEMINARIO_FORCE_CONSTANTS 1
#define USE_MAX_NORM 1

  


class IsomerspaceForcefield {
public:
  typedef GPU_REAL device_real_t;
  typedef uint16_t device_node_t;
  
  
struct DeviceGraph{
  bool allocated = false;
  size_t batch_size = 0;
  size_t N = 0;

  device_real_t* X;
  device_node_t* neighbours;
  device_node_t* next_on_face;
  device_node_t* prev_on_face;
  uint8_t* face_right;

  DeviceGraph(){}
  DeviceGraph(size_t N, size_t batch_size ,device_real_t* X, device_node_t* neighbours, device_node_t* next_on_face, device_node_t* prev_on_face, uint8_t* face_right): N(N), batch_size(batch_size), X(X), neighbours(neighbours), next_on_face(next_on_face), prev_on_face(prev_on_face), face_right(face_right){}
  DeviceGraph(const DeviceGraph& G):N(G.N), batch_size(G.batch_size), X(G.X), neighbours(G.neighbours), next_on_face(G.next_on_face), prev_on_face(G.prev_on_face), face_right(G.face_right){}
  
  void copy_to_gpu(const DeviceGraph& G);

  void allocate(const size_t N, const size_t batch_size);
  void allocate_host(const size_t N, const size_t batch_size);
  void free();
  void free_host();
}; 



struct InternalCoordinates{
  bool allocated = false;
  size_t N = 0;
  size_t batch_size = 0;

  device_real_t* bonds;
  device_real_t* angles;
  device_real_t* dihedrals;
  device_real_t* outer_angles_m;
  device_real_t* outer_angles_p;
  device_real_t* outer_dihedrals_a;
  device_real_t* outer_dihedrals_m;
  device_real_t* outer_dihedrals_p;

  InternalCoordinates(){}

  void allocate(const size_t N,const size_t batch_size);

  void free();
  static void to_file(const InternalCoordinates& coords, size_t ID_in_batch, std::string fullerene_name);

};




  static size_t get_batch_capacity(const size_t N);
  
  void insert_isomer(const FullereneGraph& G,  const vector<coord3d> &X0);
  void insert_isomer_batch(const DeviceGraph& G);
  
  void optimize_batch(size_t maxIter);
  void check_batch();
  void get_cartesian_coordinates(device_real_t* X);
  void get_internal_coordinates(device_real_t* bonds, device_real_t* angles, device_real_t* dihedrals);
  void clear_batch(){batch_size = 0;}
  IsomerspaceForcefield(const size_t N);
  ~IsomerspaceForcefield();


private:
  size_t N = 0;
  size_t batch_capacity = 0;
  size_t batch_size = 0;
  device_real_t* global_reduction_array;

  DeviceGraph d_graph;
  DeviceGraph h_graph;
  InternalCoordinates d_coords;
  
};


