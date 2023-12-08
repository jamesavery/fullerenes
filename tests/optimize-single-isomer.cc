#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/isomer_queue.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/isomerdb.hh"

using namespace gpu_kernels;
int main(int argc, char** argv){
    const size_t N                  = argc>1 ? strtol(argv[1],0,0) : 20;     // Argument 1: Number of vertices
    const bool output_geometry      = argc>2 ? strtol(argv[2],0,0) : 1;     // Argument 3: Whether to output geometry
    const size_t max_iter           = argc>3 ? strtol(argv[3],0,0) : N*5;     // Argument 2: Maximum number of iterations

    bool more_to_do = true;
    auto n_isomers = IsomerDB::number_isomers(N);
    Graph G;
    auto Nf = N/2 +2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
    G.neighbours.resize(Nf);
    G.N = Nf;
    PlanarGraph Pg;
    auto batch_size = 1;
    //BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,0,0);  
    IsomerBatch<CPU> Bhost(N,batch_size);
    IsomerBatch<GPU> Bdev(N,batch_size);
    BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,0,0);
    
    for (size_t i = 0; i < batch_size; i++)
    {
        //auto ID = cuda_benchmark::random_isomer("isomerspace_samples/dual_layout_"+to_string(N)+"_seed_42", G);
        auto more = BuckyGen::next_fullerene(BuckyQ, G);
        auto ID = i;
        Bhost.append(G,ID);
    }
    
    device_io::copy(Bdev, Bhost);
    isomerspace_dual::dualise(Bdev);
    isomerspace_tutte::tutte_layout<GPU,float>(Bdev, (int)10*N);
    device_io::copy(Bhost, Bdev);
    std::cout << "Initial 2D geometry:\n";
    for (size_t i = 0; i < N; i++)
    {
        std::cout << i << ": ";
        std::cout << Bhost.xys[i*2 + 0] << " ";
        std::cout << Bhost.xys[i*2 + 1] << " ";
        std::cout << std::endl;
    }

    isomerspace_X0::zero_order_geometry(Bdev, 4.0);
    device_io::copy(Bhost, Bdev);
    std::cout << "Initial 3D geometry:\n";
    for (size_t i = 0; i < N; i++)
    {
        std::cout << i << ": ";
        std::cout << Bhost.X[i*3 + 0] << " ";
        std::cout << Bhost.X[i*3 + 1] << " ";
        std::cout << Bhost.X[i*3 + 2] << " ";
        std::cout << std::endl;
    }


    device_io::copy(Bhost, Bdev);
    
    //Write starting geometry to file
    
    std::ofstream out_geom("starting_geometry.float32");
    std::ofstream out_graph("cubic_graphs.uint16");
    out_geom.write((char*)Bhost.X, 3*N*batch_size*sizeof(float));
    out_graph.write((char*)Bhost.cubic_neighbours, 3*N*batch_size*sizeof(uint16_t));
    out_geom.close();

    //CuArray<float> gradients(N*3);
    //isomerspace_forcefield::get_gradients<PEDERSEN>(Bdev, gradients);
//
    //std::cout << "Gradients:\n" << gradients << std::endl;

    isomerspace_forcefield::optimize<PEDERSEN>(Bdev, 5*N, 5*N);

    CuArray<float> hessians(N*90);
    CuArray<uint16_t> cols(N*90);
    LaunchCtx ctx(0);
    isomerspace_hessian::compute_hessians<PEDERSEN>(Bdev, hessians, cols, ctx, LaunchPolicy::SYNC);

    std::cout << "Hessians:\n";
    for (size_t i = 0; i < N*3; i++)
    {   
        std::cout << i << ": ";
        for (size_t j = 0; j < 30; j++)
        {
            std::cout << hessians[i*30 + j] << ", ";
        }
        std::cout << std::endl;
    }


    device_io::copy(Bhost, Bdev);

    std::cout << "Final geometry:\n";
    for (size_t i = 0; i < N; i++)
    {
        std::cout << i << ": ";
        std::cout << Bhost.X[i*3 + 0] << " ";
        std::cout << Bhost.X[i*3 + 1] << " ";
        std::cout << Bhost.X[i*3 + 2] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Cubic graph:\n";
    for (size_t i = 0; i < N*3; i++)
    {
        std::cout << Bhost.cubic_neighbours[i] << ", ";
    }
    std::cout << std::endl;
}