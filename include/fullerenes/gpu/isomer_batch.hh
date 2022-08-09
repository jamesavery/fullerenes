#ifndef ISOMERBATCH_STRUCT
#define ISOMERBATCH_STRUCT
#include "gpudatastruct.hh"
#include <stdint.h>

enum class IsomerStatus {EMPTY, CONVERGED, FAILED, NOT_CONVERGED};
enum BatchMember {COORDS3D, COORDS2D, CUBIC_NEIGHBOURS, DUAL_NEIGHBOURS, FACE_DEGREES, IDS, ITERATIONS, STATUSES};
enum SortOrder {ASCENDING, DESCENDING};
enum class LaunchPolicy {SYNC, ASYNC};

struct IsomerBatch : GPUDataStruct
{

    device_real_t* X;
    device_real_t* xys;

    device_node_t* cubic_neighbours;
    device_node_t* dual_neighbours;
    unsigned char* face_degrees;

    size_t* IDs;
    size_t* iterations;
    IsomerStatus* statuses;

    IsomerBatch(){
      pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(device_node_t)*3, true}, {"dual_neighbours", (void**)&dual_neighbours, sizeof(device_node_t)*4, true}, {"face_degrees", (void**)&face_degrees, sizeof(uint8_t)*1, true}, {"X", (void**)&X, sizeof(device_real_t)*3, true}, {"xys", (void**)&xys, sizeof(device_real_t)*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
    }

    void operator=(const IsomerBatch &);

    ~IsomerBatch() override;
    IsomerBatch(size_t n_atoms, size_t n_isomers, BufferType buffer_type, int device  = 0);
    void set_print_simple() {verbose = false;}
    void set_print_verbose() {verbose = true;}
    //Prints a specific parameter from the batch
    void print(const BatchMember param, const std::pair<int,int>& range = {-1,-1});
    int get_device_id() const {return m_device;}

    std::optional<Polyhedron> get_isomer(const size_t index) const;
    std::optional<Polyhedron> get_isomer_by_id(const size_t ID) const;
    std::vector<size_t> find_ids(const IsomerStatus status);
    void shrink_to_fit();

    bool operator==(const IsomerBatch& b);
    bool operator!=(const IsomerBatch& b) {return !(*this == b);}
    friend std::ostream& operator<<(std::ostream& os, const IsomerBatch& a);

  private:
    int m_device = 0;
    bool verbose = false;
};

#endif
