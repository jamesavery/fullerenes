#include "fullerenes/gpu/isomer_batch.hh"

template<typename T>
IsomerBatch<T,CPU>::IsomerBatch(size_t n_atoms, size_t n_isomers, int device) :
  n_atoms(n_atoms), isomer_capacity(n_isomers), n_faces(n_atoms/2+2), device(device) {
  typedef device_node_t node_t;
  
  pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(node_t)*n_atoms*3,true},
		  {"dual_neighbours", (void**)&dual_neighbours,
		                                     sizeof(node_t)*(n_atoms/2 +2) * 6, true},
		  {"face_degrees", (void**)&face_degrees,sizeof(uint8_t)*(n_atoms/2 +2),true},
		  {"X", (void**)&X, sizeof(device_real_t)*n_atoms*3, true},
		  {"xys", (void**)&xys, sizeof(device_hpreal_t)*n_atoms*2, true},
		  {"statuses", (void**)&statuses, sizeof(IsomerStatus), false},
		  {"IDs", (void**)&IDs, sizeof(size_t), false},
		  {"iterations", (void**)&iterations, sizeof(size_t), false}};

    for (size_t i = 0; i < pointers.size(); i++) {
      //For asynchronous memory transfers host memory must be pinned. 
      cudaMallocHost(get<1>(pointers[i]), n_isomers * get<2>(pointers[i]));
      memset(*get<1>(pointers[i]),0, n_isomers*get<2>(pointers[i]));
    }
    allocated = true;
}

template<typename T>
void IsomerBatch<T,CPU>::operator=(const IsomerBatch<T,CPU>& input){
  typedef device_node_t node_t;  
  cudaSetDevice(m_device);

  pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(node_t)*n_atoms*3,true},
		{"dual_neighbours", (void**)&dual_neighbours,
		                                     sizeof(node_t) * (n_atoms/2 +2)*6,true},
		{"face_degrees", (void**)&face_degrees, sizeof(uint8_t)*(n_atoms/2 +2),true},
		{"X", (void**)&X, sizeof(device_real_t)*n_atoms*3, true},
		{"xys", (void**)&xys, sizeof(device_hpreal_t)*n_atoms*2, true},
		{"statuses", (void**)&statuses, sizeof(IsomerStatus), false},
		{"IDs", (void**)&IDs, sizeof(size_t), false},
		{"iterations", (void**)&iterations, sizeof(size_t), false}};
  
    if (allocated == true){
        for (size_t i = 0; i < pointers.size(); i++) {
            cudaFree(*get<1>(pointers[i]));
        }
        allocated = false;
    }
    //Construct a tempory batch: allocates the needed amount of memory.
    this->m_device = input.m_device;
    this->isomer_capacity = input.isomer_capacity;
    this->n_atoms = input.n_atoms;
    this->n_faces = input.n_faces;
        
    //Copy contents of old batch into newly allocated memory.
        for (size_t i = 0; i < pointers.size(); i++) {
            //For asynchronous memory transfers host memory must be pinned.
#ifdef __CUDACC__	  
            cudaMallocHost(get<1>(pointers[i]), isomer_capacity * get<2>(pointers[i]));
#endif
	}
}

template<typename T>
IsomerBatch<T,CPU>::~IsomerBatch(){
    if (allocated == true);
    {
      for (size_t i = 0; i < pointers.size(); i++) 
	cudaFreeHost(*get<1>(pointers[i])); 
    }
    allocated = false;
}


