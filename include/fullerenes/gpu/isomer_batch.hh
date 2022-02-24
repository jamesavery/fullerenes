#include "gpudatastruct.hh"

enum IsomerStatus {CONVERGED, FAILED, NOT_CONVERGED, EMPTY};

struct IsomerBatch : GPUDataStruct
{
    typedef GPU_REAL device_real_t;
    typedef uint16_t device_node_t;

    device_node_t* neighbours;
    device_real_t* X;
    device_real_t* xys;

    IsomerStatus* statuses;
    size_t* IDs;
    size_t* iterations;

    IsomerBatch(){
      pointers =   {{"neighbours",(void**)&neighbours, sizeof(device_node_t)*3, true}, {"X", (void**)&X, sizeof(device_real_t)*3, true}, {"xys", (void**)&xys, sizeof(device_real_t)*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
    }
};