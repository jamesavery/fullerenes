#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include <chrono>
#include <fstream>


using namespace chrono;
using namespace chrono_literals;

#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/isomer_batch.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "numeric"
#include "random"


int main(int argc, char** argv){
    size_t N                = argc > 1 ? strtol(argv[1],0,0) : (size_t)20;     // Argument 1: Number of vertices N
    size_t ID               = argc > 2 ? strtol(argv[2],0,0) : (size_t)0;       // Argument 2: ID

    //Set up the batch.
    Graph G(N);
    auto queue = BuckyGen::start(N,false,false);
    for(size_t i = 0; i < ID; i++){
        BuckyGen::next_fullerene(queue,G);
    }

    IsomerBatch batch(N,1,DEVICE_BUFFER);
    IsomerBatch h_batch(N,1,HOST_BUFFER);
    cuda_io::IsomerQueue Q(N);
    Q.insert(G,ID);
    Q.refill_batch(batch);
    gpu_kernels::isomerspace_dual::dualise(batch);
    gpu_kernels::isomerspace_tutte::tutte_layout(batch);
    gpu_kernels::isomerspace_X0::zero_order_geometry(batch,4.);
    cuda_io::copy(h_batch,batch);

    //Polyhedron::to_mol2(h_batch.get_isomer_by_id(ID).value(),fopen("BEAUTIFUL_X0.mol2","w+")); 
    cuda_io::reset_convergence_statuses(batch);
    gpu_kernels::isomerspace_forcefield::optimise<PEDERSEN>(batch,N*1,N*5);
    batch.print(BatchMember::COORDS3D);
    Q.insert(batch);
    Polyhedron P = Q.pop();
    std::cout << P.points;
    //Polyhedron::to_mol2(batch.get_isomer_by_id(ID).value(),fopen("BEAUTIFUL_N1.mol2","w+")); 
    //gpu_kernels::isomerspace_forcefield::optimise_batch<PEDERSEN>(batch,N*2,N*5);
    //Polyhedron::to_mol2(batch.get_isomer_by_id(ID).value(),fopen("BEAUTIFUL_N3.mol2","w+")); 
    
    //batch.print(BatchMember::COORDS3D);
    //batch.print(BatchMember::DUAL_NEIGHBOURS);
    //batch.print(BatchMember::CUBIC_NEIGHBOURS);
}
