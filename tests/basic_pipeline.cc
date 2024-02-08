#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/cu_array.hh"
#include <iostream>

using namespace gpu_kernels;

int main(int argc, char** argv) {

    size_t N  = argc>1 ? strtol(argv[1],0,0) : 20;
    size_t Nf = N/2 + 2;
    size_t BatchSize = argc>2 ? strtol(argv[2],0,0) : 1;

    BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, 0, 0);
    IsomerBatch<Device::GPU> batch(N, BatchSize);
    IsomerBatch<Device::CPU> batch_cpu(N, BatchSize);
    Graph G(N);
    for (size_t ii = 0; ii < BatchSize; ii++)
    {
        auto more = BuckyGen::next_fullerene(BuckyQ, G);
        if(!more) break;
        
        for (size_t j = 0; j < Nf; j++)
        {
            for(size_t k = 0; k < G.neighbours[j].size(); k++)
            {
                batch_cpu.dual_neighbours[ii*Nf*6 + j*6 + k] = G.neighbours[j][k];
            } 
            if(G.neighbours[j].size() == 5){
                batch_cpu.dual_neighbours[ii*Nf*6 + j*6 + 5] = std::numeric_limits<uint16_t>::max();
                batch_cpu.face_degrees[ii*Nf + j] = 5;
            } else {
                batch_cpu.face_degrees[ii*Nf + j] = 6;
            }   
        }
        batch_cpu.statuses[ii] = IsomerStatus::NOT_CONVERGED;
    }

    device_io::copy(batch, batch_cpu);


    CuArray<float> hessians(N*90*BatchSize);
    CuArray<uint16_t> cols(N*90*BatchSize);
    CuArray<float> eigenvalues(N*3*BatchSize);
    CuArray<float> Qmat (N*3*N*3*BatchSize);
    isomerspace_dual::dualize(batch);
    isomerspace_tutte::tutte_layout(batch);
    isomerspace_X0::zero_order_geometry(batch,4.0);
    isomerspace_forcefield::optimize(batch, 5*N, 5*N);
    isomerspace_hessian::compute_hessians(batch, hessians, cols);
    isomerspace_eigen::eigensolve(batch, Qmat, hessians, cols, eigenvalues);
    
    /* std::cout << "Eigenvalues: " << std::endl;
    for (size_t i = 0; i < BatchSize; i++)
    {
        for (size_t j = 0; j < N*3; j++)
        {
            std::cout << eigenvalues[i*N*3 + j] << ", ";
        }
        std::cout << std::endl;
    } */

}   