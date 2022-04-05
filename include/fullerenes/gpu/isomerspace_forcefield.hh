#pragma once
#include <inttypes.h>
#include <string>
#include <queue>
#include <map>
#include "isomerspace_kernel.hh"
#include <fullerenes/fullerenegraph.hh> 
#include <fullerenes/polyhedron.hh>
#include "gpudatastruct.hh"

class IsomerspaceForcefield: public IsomerspaceKernel {
public:
  struct IsomerBatchStats : GPUDataStruct{
    device_real_t* bond_rms;
    device_real_t* angle_rms;
    device_real_t* dihedral_rms;

    device_real_t* bond_mean;
    device_real_t* angle_mean;
    device_real_t* dihedral_mean;

    device_real_t* bond_max;
    device_real_t* angle_max;
    device_real_t* dihedral_max;

    device_real_t* energy;
    device_real_t* grad_norm;
    size_t s = sizeof(device_real_t);
    
    IsomerBatchStats(){
      pointers =   {{"bond_rms",(void**)&bond_rms,s,false}, {"angle_rms",(void**)&angle_rms,s,false}, {"dihedral_rms",(void**)&dihedral_rms,s,false}, 
                    {"bond_mean", (void**)&bond_mean,s,false}, {"angle_mean", (void**)&angle_mean,s,false}, {"dihedral_mean",(void**)&dihedral_mean,s,false}, 
                    {"bond_max", (void**)&bond_max,s,false}, {"angle_max", (void**)&angle_max,s,false}, {"dihedral_max", (void**)&dihedral_max,s,false},
                    {"energy", (void**)&energy,s,false}, {"grad_norm", (void**)&grad_norm,s,false}};
    }
  };

  struct InternalCoordinates : GPUDataStruct{
    device_real_t* bonds;             
    device_real_t* angles;
    device_real_t* dihedrals;
    device_real_t* outer_angles_m;
    device_real_t* outer_angles_p;
    device_real_t* outer_dihedrals_a;
    device_real_t* outer_dihedrals_m;
    device_real_t* outer_dihedrals_p;

    size_t s = sizeof(device_real_t)*3;

    InternalCoordinates(){
      pointers = {{"bonds", (void**)&bonds, s,true}, {"angles", (void**)&angles, s,true}, {"dihedrals", (void**)&dihedrals, s,true}, {"outer_angles_m", (void**)&outer_angles_m, s,true}, 
                  {"outer_angles_p", (void**)&outer_angles_p, s,true}, {"outer_dihedrals_a", (void**)&outer_dihedrals_a, s,true}, {"outer_dihedrals_m", (void**)&outer_dihedrals_m, s,true}, 
                  {"outer_dihedrals_p", (void**)&outer_dihedrals_p, s,true}};
    }
  };

  struct SharedQueue
  {
    IsomerBatch data;
    size_t q_size, capacity, front, back;

    SharedQueue(const size_t N, const size_t capacity): q_size(0), capacity(capacity), front(0), back(0) {GPUDataStruct::allocate(this->data,N,capacity,DEVICE_BUFFER);}
  };
  
  void push_batch_from_kernel(const IsomerspaceKernel& input_kernel);
  void update_device_batches();
  void check_batch(size_t max_iterations);                   //Checks convergence properties of current batch, calculates mean and std of relative bond, angle and dihedral errors of the current batch.
  void optimize_batch(const size_t iterations);

  void get_cartesian_coordinates(device_real_t* X) const;                                                     //Populate target buffer (CPU) with cartesian coordiantes from isomers on GPU.
  void get_internal_coordinates(device_real_t* bonds, device_real_t* angles, device_real_t* dihedrals); //Populate target buffers (CPU) with internal coordinates from isomers on GPU.
  
  void to_file(size_t ID_in_batch);   //Reads and dumps all graph information, cartesian coordinates and harmonic constants to files.
  void batch_statistics_to_file();

  size_t get_converged_count()const{return converged_count;};
  size_t get_failed_count()const{return failed_count;}

  IsomerspaceForcefield(const size_t N);
  ~IsomerspaceForcefield();


protected:
  std::vector<IsomerBatch> d_queues;
  std::vector<InternalCoordinates> d_coords;     //Provided for diagnostic purposes.                           Dimensions: N x 1 x 3
  std::vector<InternalCoordinates> h_coords;     //Provided for diagnostic purposes.                           Dimensions: N x 1 x 3
  std::vector<InternalCoordinates> d_harmonics;  //Provided for diagnostic purposes.                           Dimensions: N x 1 x 3
  std::vector<InternalCoordinates> h_harmonics;  //Provided for diagnostic purposes.                           Dimensions: N x 1 x 3
};


