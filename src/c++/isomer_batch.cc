#include "fullerenes/isomer_batch.hh"

template<Device T>
IsomerBatch<T>::IsomerBatch(size_t n_atoms, size_t n_isomers, int device) :
  n_atoms(n_atoms), isomer_capacity(n_isomers), n_faces(n_atoms/2+2), m_device(device) {
  typedef device_node_t node_t;
  
  pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(node_t)*n_atoms*3,true},
		  {"dual_neighbours", (void**)&dual_neighbours,
		                                     sizeof(node_t)*(n_atoms/2 +2) * 6, true},
		  {"face_degrees", (void**)&face_degrees,sizeof(uint8_t)*(n_atoms/2 +2),true},
		  {"X", (void**)&X, sizeof(device_real_t)*n_atoms*3, true},
		  {"xys", (void**)&xys, sizeof(device_real_t)*n_atoms*2, true},
		  {"statuses", (void**)&statuses, sizeof(IsomerStatus), false},
		  {"IDs", (void**)&IDs, sizeof(size_t), false},
		  {"iterations", (void**)&iterations, sizeof(size_t), false}};

  for (auto [name, ptr,isomer_size]: pointers) {      
      *ptr = calloc(n_isomers, isomer_size);
      memset(*ptr,0, n_isomers*isomer_size);
    }
    allocated = true;
}

template<Device T>
void IsomerBatch<T>::operator=(const IsomerBatch<T>& input){
  typedef device_node_t node_t;  
  pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(node_t)*n_atoms*3,true},
		{"dual_neighbours", (void**)&dual_neighbours,
		                                     sizeof(node_t) * (n_atoms/2 +2)*6,true},
		{"face_degrees", (void**)&face_degrees, sizeof(uint8_t)*(n_atoms/2 +2),true},
		{"X", (void**)&X, sizeof(device_real_t)*n_atoms*3, true},
		{"xys", (void**)&xys, sizeof(device_real_t)*n_atoms*2, true},
		{"statuses", (void**)&statuses, sizeof(IsomerStatus), false},
		{"IDs", (void**)&IDs, sizeof(size_t), false},
		{"iterations", (void**)&iterations, sizeof(size_t), false}};
  
    if (allocated == true){	// Do we need to destroy an old batch in the LHS?
      for (auto [name,ptr,isomer_size,whatisthis] : pointers)
	free(*ptr);
      allocated = false;
    }
    //Construct a tempory batch: allocates the needed amount of memory.
    m_device = input.m_device;
    isomer_capacity = input.isomer_capacity;
    n_atoms = input.n_atoms;
    n_faces = input.n_faces;
        
    //Copy contents of old batch into newly allocated memory.
    for (size_t i=0;i<pointers.size();i++){
      auto [name,ptr,isomer_size,whatisthis] = pointers[i];
      auto from_ptr = get<1>(input.pointers[i]);
      *ptr = calloc(input.isomer_capacity,isomer_size);
	    
      memcpy(*ptr,*from_ptr,isomer_capacity*isomer_size);
    }
}

template<Device T>
IsomerBatch<T>::~IsomerBatch(){
    if (allocated == true);
    {
      for (size_t i = 0; i < pointers.size(); i++) 
	free(*get<1>(pointers[i])); 
    }
    allocated = false;
}


