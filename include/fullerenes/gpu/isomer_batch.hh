#ifndef ISOMERBATCH_STRUCT
#define ISOMERBATCH_STRUCT
#include "gpudatastruct.hh"
#include <stdint.h>

enum IsomerStatus {EMPTY, CONVERGED, FAILED, NOT_CONVERGED};
enum class LaunchPolicy {SYNC, ASYNC};

struct IsomerBatch : GPUDataStruct
{

    device_node_t* cubic_neighbours;
    device_node_t* dual_neighbours;
    unsigned char* face_degrees;
    device_real_t* X;
    device_real_t* xys;

    IsomerStatus* statuses;
    size_t* IDs;
    size_t* iterations;

    IsomerBatch(){
      pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(device_node_t)*3, true}, {"dual_neighbours", (void**)&dual_neighbours, sizeof(device_node_t)*4, true}, {"face_degrees", (void**)&face_degrees, sizeof(uint8_t)*1, true}, {"X", (void**)&X, sizeof(device_real_t)*3, true}, {"xys", (void**)&xys, sizeof(device_real_t)*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
    }

    IsomerBatch(size_t n_atoms, size_t n_isomers, BufferType buffer_type);

    bool operator==(const IsomerBatch& b);
    bool operator!=(const IsomerBatch& b) {return !(*this == b);}
    friend std::ostream& operator<<(std::ostream& os, const IsomerBatch& a);
};

#endif
