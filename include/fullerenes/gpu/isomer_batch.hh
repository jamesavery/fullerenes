#ifndef ISOMERBATCH_STRUCT
#define ISOMERBATCH_STRUCT
#include "gpudatastruct.hh"

enum IsomerStatus {EMPTY, CONVERGED, FAILED, NOT_CONVERGED};
enum class LaunchPolicy {SYNC, ASYNC};

struct IsomerBatch : GPUDataStruct
{

    device_node_t* neighbours;
    device_real_t* X;
    device_real_t* xys;

    IsomerStatus* statuses;
    size_t* IDs;
    size_t* iterations;

    IsomerBatch(){
      pointers =   {{"neighbours",(void**)&neighbours, sizeof(device_node_t)*3, true}, {"X", (void**)&X, sizeof(device_real_t)*3, true}, {"xys", (void**)&xys, sizeof(device_real_t)*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
    }

    IsomerBatch(size_t n_atoms, size_t n_isomers, BufferType buffer_type) {
      pointers =   {{"neighbours",(void**)&neighbours, sizeof(device_node_t)*3, true}, {"X", (void**)&X, sizeof(device_real_t)*3, true}, {"xys", (void**)&xys, sizeof(device_real_t)*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
      allocate(*this,n_atoms, n_isomers, buffer_type);
      initialize(*this);
    }

    bool operator==(const IsomerBatch& b);
    bool operator!=(const IsomerBatch& b) {return !(*this == b);}
    friend std::ostream& operator<<(std::ostream& os, const IsomerBatch& a);
};

#endif
