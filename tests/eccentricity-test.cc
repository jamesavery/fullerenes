#include "string"
#include "algorithm"
#include "iostream"
#include "numeric"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/isomer_queue.hh"
#include "fullerenes/isomerdb.hh"

using namespace gpu_kernels;
using namespace isomerspace_dual;
using namespace isomerspace_eigen;
using namespace isomerspace_forcefield;
using namespace isomerspace_hessian;
using namespace isomerspace_X0;
using namespace isomerspace_tutte;

int main(int argc, char** argv){
    size_t N = argc > 1 ? atoi(argv[1]) : 60;
    auto n_isomers = IsomerDB::number_isomers(N);
    Graph G;
    auto Nf = N/2 +2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
    G.neighbours.resize(Nf);
    G.N = Nf;
    PlanarGraph Pg;
    auto batch_size = min(isomerspace_forcefield::optimal_batch_size(N)*8, (int)n_isomers);
    //BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,0,0);  
    IsomerBatch<CPU> Bhost(N,batch_size);
    IsomerBatch<GPU> Bdev(N,batch_size);
    for (size_t i = 0; i < batch_size; i++)
    {
        auto ID = cuda_benchmark::random_isomer("isomerspace_samples/dual_layout_"+to_string(N)+"_seed_42", G);
        Bhost.append(G,ID);
    }
    device_io::copy(Bdev, Bhost);
    dualize(Bdev);
    tutte_layout(Bdev, N*10);
    zero_order_geometry(Bdev, 4.0);
    optimize<PEDERSEN>(Bdev, N*5, N*5);
    device_io::copy(Bhost, Bdev);
    
    CuArray<device_real_t> ECC(batch_size);
    isomerspace_properties::eccentricities(Bdev, ECC);
    auto tol = std::numeric_limits<device_real_t>::epsilon()*1e4;
    for (size_t i = 0; i < batch_size; i++)
    {
        auto P = Bhost.get_isomer(i).value();
        auto matrix = P.inertia_matrix();
        if (isnan(matrix.norm())) continue;
        auto lambda = matrix.eigenvalues();
        if(fabs(lambda[0]) > fabs(lambda[1])) std::swap(lambda[0],lambda[1]);
        if(fabs(lambda[1]) > fabs(lambda[2])) std::swap(lambda[1],lambda[2]);
        if(fabs(lambda[0]) > fabs(lambda[1])) std::swap(lambda[0],lambda[1]);
        auto ref = (1./sqrt(lambda[0]))  / (1./sqrt(lambda[2]));
        auto val = ECC[i];
        if (abs(ref-val)/val > tol){
            std::cout << "Eccentricity test failed for isomer " << i << " relative error = " << abs(ref-val)/val << std::endl;
            return 1;
        }
    }
    std::cout << "Eccentricity test passed for C" << N << std::endl;
}
