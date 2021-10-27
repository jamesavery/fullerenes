#pragma once
#include <inttypes.h>
#include <string>
#include <fullerenes/fullerenegraph.hh> 

#define GPU_REAL float       
#define GPU_REAL3 float3
#define GPU_NODE3 ushort3
#define MAKE_NODE3 make_ushort3
#define Block_Size_Pow_2 512
#define SEMINARIO_FORCE_CONSTANTS 1
#define USE_MAX_NORM 0
#define REDUCTION_METHOD 0

  


class IsomerspaceForcefield {
public:
  typedef GPU_REAL device_real_t;
  typedef uint16_t device_node_t;
  enum BufferType {host_buffer, device_buffer};

  struct GenericStruct{
    bool allocated = false;
    size_t N = 0;
    size_t batch_size = 0;
    BufferType buffer_type;
    GenericStruct(){}
    std::vector<std::tuple<std::string,void**,size_t>> pointers;
    static void allocate(GenericStruct& G,const size_t N,const size_t batch_size, const BufferType buffer_type);
    static void free(GenericStruct& G);
    void set_pointers(std::vector<std::tuple<std::string,void**,size_t>> &pointers){this->pointers = pointers;}
  };

  struct IsomerspaceStats : GenericStruct{
    device_real_t* bond_rms;
    device_real_t* angle_rms;
    device_real_t* dihedral_rms;

    device_real_t* bond_mean;
    device_real_t* angle_mean;
    device_real_t* dihedral_mean;

    device_real_t* bond_max;
    device_real_t* angle_max;
    device_real_t* dihedral_max;

    size_t s = sizeof(device_real_t);
    std::vector<std::tuple<std::string,void**,size_t>> pointers =   {{"bond_rms",(void**)&bond_rms, s}, {"angle_rms",(void**)&angle_rms, s}, {"dihedral_rms",(void**)&dihedral_rms ,s}, 
                                                                    {"bond_mean", (void**)&bond_mean, s}, {"angle_mean", (void**)&angle_mean, s}, {"dihedral_mean",(void**)&dihedral_mean, s}, 
                                                                    {"bond_max", (void**)&bond_max, s}, {"angle_max", (void**)&angle_max, s}, {"dihedral_max", (void**)&dihedral_max, s}};

    IsomerspaceStats(){set_pointers(pointers);}
  };
  
  
  struct IsomerspaceGraph : GenericStruct{
    device_real_t* X;             //Cartesian coordinates.
    device_node_t* neighbours;    //Cubic neighbour list.
    device_node_t* next_on_face;  //Next node on face / counter-clockwise direction
    device_node_t* prev_on_face;  //Prev node on face / clockwise direction
    uint8_t* face_right;          //Face size list.

    //Since we are dealing with data of non-uniform types we must encapsulate this information if we are to iterate over the pointers for allocation, deallocation and copying.
    std::vector<std::tuple<std::string,void**,size_t>> pointers = {{"X",(void**)&X,sizeof(device_real_t)}, {"neighbours",(void**)&neighbours, sizeof(device_node_t)}, 
                                                                  {"next_on_face", (void**)&next_on_face, sizeof(device_node_t)}, {"prev_on_face", (void**)&prev_on_face, sizeof(device_node_t)},  
                                                                  {"face_right", (void**)&face_right, sizeof(uint8_t)}};

    
    IsomerspaceGraph(){set_pointers(pointers);}
    IsomerspaceGraph(device_real_t* X, device_node_t* neighbours, device_node_t* next_on_face, device_node_t* prev_on_face, uint8_t* face_right): X(X), neighbours(neighbours), next_on_face(next_on_face), prev_on_face(prev_on_face), face_right(face_right){}
    void copy_to_gpu(const IsomerspaceGraph& G);
  }; 



  struct InternalCoordinates : GenericStruct{
    device_real_t* bonds;             
    device_real_t* angles;
    device_real_t* dihedrals;
    device_real_t* outer_angles_m;
    device_real_t* outer_angles_p;
    device_real_t* outer_dihedrals_a;
    device_real_t* outer_dihedrals_m;
    device_real_t* outer_dihedrals_p;

    size_t s = sizeof(device_real_t);
    std::vector<std::tuple<std::string,void**,size_t>> pointers = {{"bonds", (void**)&bonds, s}, {"angles", (void**)&angles, s}, {"dihedrals", (void**)&dihedrals, s}, {"outer_angles_m", (void**)&outer_angles_m, s}, 
                                                                    {"outer_angles_p", (void**)&outer_angles_p, s}, {"outer_dihedrals_a", (void**)&outer_dihedrals_a, s}, {"outer_dihedrals_m", (void**)&outer_dihedrals_m, s}, 
                                                                    {"outer_dihedrals_p", (void**)&outer_dihedrals_p, s}};
    InternalCoordinates(){set_pointers(pointers);}
  };

  static size_t get_batch_capacity(const size_t N); //Uses Cuda API calls to determine the amount of fullerenes of a given size N, that can be optimized simultaneously.
  size_t get_batch_size()const{return batch_size;}
  void insert_isomer(const FullereneGraph& G,  const vector<coord3d> &X0);  //Essentially adapter pattern converting FullereneGraph objects into 1D arrays and inserting them into an IsomerspaceGraph. 
  void insert_isomer_batch(const IsomerspaceGraph& G);                      //Inserts an entire batch, used at the moment for inserting synthetic loads, 
                                                                            //could be used in the future for more efficient transfer from CPU to GPU.
  
  void optimize_batch(size_t maxIter);  //Performs Conjugate Gradient Forcefield optimization on a fullerene isomer batch.
  void check_batch();                   //Checks convergence properties of current batch, calculates mean and std of relative bond, angle and dihedral errors of the current batch.
  
  void get_cartesian_coordinates(device_real_t* X);                                                     //Populate target buffer (CPU) with cartesian coordiantes from isomers on GPU.
  void get_internal_coordinates(device_real_t* bonds, device_real_t* angles, device_real_t* dihedrals); //Populate target buffers (CPU) with internal coordinates from isomers on GPU.
  
  void clear_batch(){batch_size = 0;} //Clears batch, this is required after every batch is finished, effectively resets the position of pointer to GPU memory
  void to_file(size_t ID_in_batch);   //Reads and dumps all graph information, cartesian coordinates and harmonic constants to files.
  void batch_statistics_to_file();

  IsomerspaceForcefield(const size_t N); //Simple constructor allocates memory for structs on device and host call this once at the beginning of program, 
                                          //also serves to initialize the cuda default context, which would otherwise be destroyed and recreated every time memory is freed and reallocated.
  ~IsomerspaceForcefield();              //Destructor, calls free and delete on GPU and CPU buffers respectively.


protected:
  size_t N = 0;                           //Size of each isomer, #of carbon atoms.
  size_t batch_capacity = 0;              //Number of fullerenes that will fit on the GPU concurrently.
  size_t batch_size = 0;                  //Current number of fullerenes copied to GPU.
  size_t shared_memory_bytes = 0;         //Amount of L1 cache to allocate per block.
  size_t isomer_number = 0;               //Isomer number of the first fullerene in batch.
  device_real_t* global_reduction_array;  //Array used to communicate across blocks.


  IsomerspaceGraph d_graph;         //GPU container for graph information and X0.                 Dimensions: N x M x 3
  IsomerspaceGraph h_graph;         //Host buffer for graph information and X0.                   Dimensions: N x M x 3

  InternalCoordinates d_coords;     //Provided for diagnostic purposes.                           Dimensions: N x 1 x 3
  InternalCoordinates h_coords;     //Provided for diagnostic purposes.                           Dimensions: N x 1 x 3
  InternalCoordinates d_harmonics;  //Provided for diagnostic purposes.                           Dimensions: N x 1 x 3
  InternalCoordinates h_harmonics;  //Provided for diagnostic purposes.                           Dimensions: N x 1 x 3

  IsomerspaceStats d_stats;         //GPU buffers containing statistical batch information.       Dimensions: 1 x M x 1
  IsomerspaceStats h_stats;         //Host buffers containing statistical batch information.      Dimensions: 1 x M x 1      
};


