#include <fullerenes/graph.hh>
#include <fullerenes/sycl-kernels.hh>
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <string>
#include <algorithm>
#define PRINT_CHECK 1
using namespace sycl;
int main(int argc, char** argv) {
    typedef float real_t;
    typedef uint16_t node_t;

    size_t N  = argc>1 ? strtol(argv[1],0,0) : 20;
    size_t Nf = N/2 + 2;
    size_t BatchSize = argc>2 ? strtol(argv[2],0,0) : 1;
    std::string device_type = argc>3 ? argv[3] : "gpu";

    auto selector =  device_type == "cpu" ? sycl::cpu_selector_v : sycl::gpu_selector_v;

    sycl::queue Q = sycl::queue(selector, sycl::property::queue::in_order{});
    
    BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, 0, 0);
    IsomerBatch<real_t,node_t> batch(N, BatchSize);
    Graph G(N);
    //Fill the batch
    {

    sycl::host_accessor acc_dual(batch.dual_neighbours, sycl::write_only);
    sycl::host_accessor acc_degs(batch.face_degrees, sycl::write_only);
    sycl::host_accessor acc_status (batch.statuses, sycl::write_only);
    for (size_t ii = 0; ii < BatchSize; ii++)
    {
        auto more = BuckyGen::next_fullerene(BuckyQ, G);
        if(!more) break;

        for (size_t j = 0; j < Nf; j++)
        {
            for(size_t k = 0; k < G.neighbours[j].size(); k++)
            {
                acc_dual[ii*Nf*6 + j*6 + k] = G.neighbours[j][k];
            } 
            if(G.neighbours[j].size() == 5){
                acc_dual[ii*Nf*6 + j*6 + 5] = std::numeric_limits<node_t>::max();
                acc_degs[ii*Nf + j] = 5;
            } else {
                acc_degs[ii*Nf + j] = 6;
            }   

        }
        acc_status[ii] = IsomerStatus::NOT_CONVERGED;
    }
    }
    
    dualize(Q, batch, LaunchPolicy::SYNC);
    tutte_layout(Q, batch, LaunchPolicy::SYNC);
    spherical_projection(Q, batch, LaunchPolicy::SYNC);
    forcefield_optimize(Q, batch, 5*N, 5*N, LaunchPolicy::SYNC);

    sycl::buffer<real_t, 1> hessians(range<1>(N*90*BatchSize));
    sycl::buffer<node_t, 1> cols(range<1>(N*90*BatchSize));
    sycl::buffer<real_t, 1> eigenvalues(range<1>(N*3*BatchSize));
    compute_hessians(Q, batch, hessians, cols, LaunchPolicy::SYNC);
    eigensolve<EigensolveMode::FULL_SPECTRUM>(Q, batch, hessians, cols, eigenvalues);

    #if PRINT_CHECK
    std::cout << "Output Cubic Graph:" << std::endl;
    for (size_t ii = 0; ii < BatchSize; ii++)
    {
        sycl::host_accessor acc_cubic(batch.cubic_neighbours, sycl::read_only);
        for (size_t i = 0; i < N; i++)
        {
            std::cout << i << ": ";
            for (size_t j = 0; j < 3; j++)
            {
                std::cout << acc_cubic[ii*N*3 + i*3 + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << "Output 2D Coordinates:" << std::endl;
    for (size_t ii = 0; ii < BatchSize; ii++)
    {
        sycl::host_accessor acc_xys(batch.xys, sycl::read_only);
        for (size_t i = 0; i < N; i++)
        {
            std::cout << i << ": ";
            std::cout << acc_xys[ii*N + i][0] << " ";
            std::cout << acc_xys[ii*N + i][1] << " ";
            std::cout << std::endl;
        }
    }

    std::cout << "Initial 3D Coordinates:" << std::endl;
    for (size_t ii = 0; ii < BatchSize; ii++)
    {
        sycl::host_accessor acc_X(batch.X, sycl::read_only);
        for (size_t i = 0; i < N; i++)
        {
            std::cout << i << ": ";
            std::cout << acc_X[ii*N + i][0] << " ";
            std::cout << acc_X[ii*N + i][1] << " ";
            std::cout << acc_X[ii*N + i][2] << " ";
            std::cout << std::endl;
        }
    }

    {   
        std::cout << "Output Hessians:" << std::endl;
        sycl::host_accessor acc_hessians(hessians, sycl::read_only);
        for (size_t ii = 0; ii < BatchSize; ii++)
        {
            for (size_t i = 0; i < N*3; i++)
            {
                std::cout << i << ": ";
                for (size_t j = 0; j < 30; j++)
                {
                    std::cout << acc_hessians[ii*N*90 + i*30 + j] << ", ";
                }
                std::cout << std::endl;
            }
        }
    }
    
    {   
        std::cout << "Output Eigenvalues:" << std::endl;
        sycl::host_accessor acc_eigenvalues(eigenvalues, sycl::read_only);
        for (size_t ii = 0; ii < BatchSize; ii++)
        {
            for (size_t i = 0; i < N*3; i++)
            {
                    std::cout << acc_eigenvalues[ii*N*3 + i] << ", ";
            }
            std::cout << std::endl;
        }
    }
    #endif

    

    std::cout << "Hello, world!" << std::endl;
    return 0;
}
