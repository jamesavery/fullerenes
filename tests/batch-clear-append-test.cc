#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/config.h"
#include <chrono>
#include <fstream>
#include "random"
#include "numeric"
#include "fullerenes/isomer_queue.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/isomerdb.hh"

using namespace chrono;
using namespace chrono_literals;

int main(int argc, char** argv){
    using namespace gpu_kernels;

    //Initialize a BuckyGen queue with N vertices
    const size_t N = strtol(argv[1],0,0);     // Argument 1: Number of vertices N
    BuckyGen::buckygen_queue Q = BuckyGen::start(N,false,false);  

    //Optimal batch size for N vertices
    auto M_b = min(gpu_kernels::isomerspace_forcefield::optimal_batch_size(N,0),(int)IsomerDB::number_isomers(N));

    //Number of faces for N vertices
    auto Nf = N/2 + 2;

    bool more_to_do = true;
    size_t I = 0;

    FullereneDual G(Nf);                                          //Graph for BuckyGen to fill
    IsomerBatch<CPU> B0(N,M_b);                //Host batch     
    IsomerBatch<GPU> B1(N,M_b);              //Device batch
    while (more_to_do)
    {
        while (B0.size() < B0.capacity())
        {
            more_to_do &= BuckyGen::next_fullerene(Q,G);            //Generate next fullerene
            if (!more_to_do) break;                                     
            G.update();                                             
            PlanarGraph pG = G.dual_graph();                        //Compute Cubic Graph
            pG.layout2d    = pG.tutte_layout();                     //Compute 2D Embedding
            Polyhedron P(pG);                                       //Create polyhedron            
            P.points       = P.zero_order_geometry();               //Compute 3D Embedding
            B0.append(P, I);                                        //Append to batch
            I++;
        }
        if (B0.size() == 0) break;
        device_io::copy(B1, B0);                                      //Copy to device B1 <- B0
        isomerspace_forcefield::optimise<PEDERSEN>(B1, N*5, N*5);     //Forcefield Optimization
        B0.clear();
        //Do something with results from B1 next... (Future work)
    }
    LaunchCtx::clear_allocations();
}
