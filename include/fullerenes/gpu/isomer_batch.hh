#ifndef ISOMERBATCH_STRUCT
#define ISOMERBATCH_STRUCT
#include "fullerenes/gpu/cuda_definitions.h"
#include "fullerenes/gpu/launch_ctx.hh"
#include "gpudatastruct.hh"
#include <stdint.h>
#include <optional>
#include <vector>
#include "fullerenes/polyhedron.hh"
#include "fullerenes/graph.hh"

enum class IsomerStatus {EMPTY, CONVERGED, PLZ_CHECK, FAILED, NOT_CONVERGED};
enum BatchMember {COORDS3D, COORDS2D, CUBIC_NEIGHBOURS, DUAL_NEIGHBOURS, FACE_DEGREES, IDS, ITERATIONS, STATUSES};
enum SortOrder {ASCENDING, DESCENDING};
enum class LaunchPolicy {SYNC, ASYNC};
enum BufferType   {HOST_BUFFER, DEVICE_BUFFER};
template <BufferType T>
struct IsomerBatch
{

    float* X;
    device_hpreal_t* xys;

    device_node_t* cubic_neighbours;
    device_node_t* dual_neighbours;
    unsigned char* face_degrees;

    size_t* IDs;
    size_t* iterations;
    IsomerStatus* statuses;

    IsomerBatch(){
      pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(device_node_t)*3, true}, {"dual_neighbours", (void**)&dual_neighbours, sizeof(device_node_t)*4, true}, {"face_degrees", (void**)&face_degrees, sizeof(uint8_t)*1, true}, {"X", (void**)&X, sizeof(device_real_t)*3, true}, {"xys", (void**)&xys, sizeof(device_hpreal_t)*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
    }

    void operator=(const IsomerBatch &);

    ~IsomerBatch() override;
    IsomerBatch(size_t n_atoms, size_t n_isomers, int device  = 0);
    void set_print_simple() {verbose = false;} 
    void set_print_verbose() {verbose = true;} 
    bool get_print_mode() const {return verbose;}
    //Prints a specific parameter from the batch
    void print(const BatchMember param, const std::pair<int,int>& range = {-1,-1}); 
    int get_device_id() const {return m_device;} 
    int size() const {return m_size;}
    int capacity() const {return isomer_capacity;}

    std::optional<Polyhedron> get_isomer(const size_t index) const; //Returns the isomer at the given index
    std::optional<Polyhedron> get_isomer_by_id(const size_t ID) const; //Returns the isomer with the given ID
    std::vector<size_t> find_ids(const IsomerStatus status); //Returns a vector of IDs with a given status
    void shrink_to_fit();        
    void append(const Graph& G, const size_t id);  //Appends a graph to the batch and increments the size
    void append(const PlanarGraph& G, const size_t id, const bool copy_2d_layout = true); //Appends a planar graph to the batch and increments the size
    void append(const Polyhedron& P, const size_t id); //Appends a polyhedron to the batch and increments the size
    void clear(const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy = LaunchPolicy::SYNC);                 //Clears the batch and resets the size to 0
    bool operator==(const IsomerBatch& b); //Returns true if the two batches are equal
    bool operator!=(const IsomerBatch& b) {return !(*this == b);}
    friend std::ostream& operator<<(std::ostream& os, const IsomerBatch& a); //Prints the batch to the given stream

  private:
    int m_size = 0;
    int m_device = 0;
    bool verbose = false;

};

#endif
