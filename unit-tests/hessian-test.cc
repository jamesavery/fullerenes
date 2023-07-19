#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/isomerdb.hh"

using namespace gpu_kernels;
int main(int argc, char** argv){
    const size_t N                = argc>1 ? strtol(argv[1],0,0) : 60;     // Argument 1: Number of vertices 
    float reldelta                = argc>2 ? strtof(argv[2],0) : 1e-5;    // Argument 3: Relative delta
    int isomer_num                = argc>3 ? strtol(argv[3],0,0) : 0;     // Argument 4: Isomer number
    std::string spiral_           = argc>4 ? argv[4] : "C60-[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene";          // Argument 2: Spiral
    std::string name_             = argc>5 ? argv[5] : "C60ih";          // Argument 2: Spiral
    std::string filename          = argc>6 ? argv[6] : "hessian_validation";        // Argument 2: Filename

    std::ofstream hess_analytical(filename + "_analytical.csv");
    std::ofstream hess_numerical(filename + "_numerical.csv");
    std::ofstream hess_neighbours(filename + "_neighbours.csv");
    std::ofstream cubic_graph;//("C" + to_string(N) + "_CubicGraph_" + to_string(isomer_num) + ".bin");
    std::ofstream geometry;//("C" + to_string(N) + "_Geometry_" + to_string(isomer_num) + ".bin");
    if (isomer_num == -1){
        cubic_graph = std::ofstream(name_ + "_CubicGraph.bin");
        geometry = std::ofstream(name_ + "_Geometry.bin");
    }

    bool more_to_do = true;
    auto n_isomers = IsomerDB::number_isomers(N);
    Graph G;
    auto Nf = N/2 +2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
    G.neighbours.resize(Nf);
    G.N = Nf;
    PlanarGraph Pg;
    auto batch_size = min(1, (int)n_isomers);
    //BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,0,0);  
    IsomerBatch Bhost(N,batch_size,CPU);
    IsomerBatch Bdev(N,batch_size,GPU);
    
    if (isomer_num == -1) {
        spiral_nomenclature C60name(spiral_);    
        Triangulation C60dual(C60name);
        auto C60cubic = C60dual.dual_graph();
        for(int i = 0;  i < batch_size; i++) Bhost.append(C60cubic,0, false);
    } else {
        for (size_t i = 0; i < batch_size; i++)
        {
            auto ID = cuda_benchmark::random_isomer("isomerspace_samples/dual_layout_"+to_string(N)+"_seed_42", G);
            Bhost.append(G,ID);
        }
    }
    
    cuda_io::copy(Bdev, Bhost);
    if(isomer_num != -1) isomerspace_dual::dualise(Bdev);
    isomerspace_tutte::tutte_layout(Bdev, (int)10*N);
    isomerspace_X0::zero_order_geometry(Bdev, 4.0);
    isomerspace_forcefield::optimise<PEDERSEN>(Bdev, 3*N, 3*N);

    CuArray<device_real_t> hessians(N*3*3*10 * batch_size);
    CuArray<device_node_t> cols(N*3*3*10 * batch_size);
    CuArray<device_real_t> eigenvalues(N*3 * batch_size);

    isomerspace_hessian::compute_hessians<PEDERSEN>(Bdev, hessians, cols);


    isomerspace_eigen::eigensolve_cusolver(Bdev, hessians, cols, eigenvalues);
    auto start = std::chrono::steady_clock::now();
    isomerspace_eigen::eigensolve_cusolver(Bdev, hessians, cols, eigenvalues);
    auto time = std::chrono::steady_clock::now() - start;
    std::cout << "cuSOLVE Eigensolver took: " << (time/1us)/(float)batch_size << " us / graph" << std::endl;
    CuArray<device_real_t> Q(N*3*N*3*batch_size);
    start = std::chrono::steady_clock::now();
    isomerspace_eigen::eigensolve(Bdev, Q, hessians, cols, eigenvalues);
    time = std::chrono::steady_clock::now() - start;
    std::cout << "Custom Eigensolver took: " << (time/1us)/(float)batch_size << " us / graph" << std::endl;

    for(int i = N*2; i < N*3; i++) {
        //std::cout << "Q[" << i  << "] : " << std::vector<device_real_t>(Q.data + i*N*3, Q.data + (i+1)*N*3) << std::endl;
    }
    //std::cout << "Eigenvalues: " << eigenvalues << std::endl;
    //FILE* f = fopen((name_ + "_Geometry.mol2").c_str(), "w");
    //FG.to_mol2(FG,f);

//    std::cout << Bdev << std::endl;
}
