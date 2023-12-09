#include "numeric"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/isomer_queue.hh"
#include "fullerenes/progress_bar.hh"
#include "fullerenes/isomerdb.hh"

using namespace gpu_kernels;
using namespace isomerspace_dual;
using namespace isomerspace_eigen;
using namespace isomerspace_forcefield;
using namespace isomerspace_hessian;
using namespace isomerspace_X0;
using namespace isomerspace_tutte;

int main(int argc, char** argv){
    const size_t N                  = argc>1 ? strtol(argv[1],0,0) : 60;     // Argument 1: Number of vertices 
    const size_t batch_size         = min(argc>2 ? strtol(argv[2],0,0) : 100, IsomerDB::number_isomers(N));      // Argument 2: Number of isomers to generate
    
    std::string type = "float32";
    
    bool more_to_do = true;

    Graph G;
    auto Nf = N/2 +2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
    G.neighbours.resize(Nf);
    G.N = Nf;

    CuArray<device_real_t> hessians(N*3*3*10 * batch_size); 
    CuArray<device_node_t> cols(N*3*3*10 * batch_size);
    CuArray<device_real_t> Q(N*3*N*3*batch_size); CuArray<device_real_t> Q_validation(N*3*N*3*batch_size);
    CuArray<device_real_t> eigs(N*3*batch_size); CuArray<device_real_t> eigs_validation(N*3*batch_size);

    PlanarGraph Pg;
    //BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N,0,0);  
    //IsomerQueue Q0(N);
    IsomerBatch<CPU> Bhost(N,batch_size);

    int Nd = LaunchCtx::get_device_count();
    IsomerBatch<GPU> Bdev(N,batch_size,0);
    BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, false, false);  

    for (size_t i = 0; i < batch_size; i++)
    {   
        bool success = BuckyGen::next_fullerene(BuckyQ, G);
        if(success) Bhost.append(G,i);
        else break;
    }

    
    device_io::copy(Bdev, Bhost);
    
    LaunchCtx ctx(0);
    LaunchPolicy policy = LaunchPolicy::SYNC;
    
    dualize(Bdev);
    tutte_layout(Bdev, (int)20*N);
    zero_order_geometry(Bdev, 4.0);
    optimize<PEDERSEN>(Bdev, 5*N, 6*N);
    optimize<PEDERSEN>(Bdev, 1*N, 6*N);
    isomerspace_properties::transform_coordinates(Bdev);
    compute_hessians<PEDERSEN>(Bdev, hessians, cols);

    eigensolve(Bdev, Q, hessians, cols, eigs);    
    eigensolve_special(Bdev, Q_validation, hessians, cols, eigs_validation);
    CuArray<device_real_t> max_eigvects(N*3*batch_size);
    CuArray<device_real_t> max_eigvals(batch_size);
    CuArray<device_real_t> min_eigvals(batch_size);
    CuArray<device_real_t> min_eigvects(N*3*batch_size);

    spectrum_ends(Bdev, hessians, cols, min_eigvals, max_eigvals, min_eigvects, max_eigvects, N*3-6);
    device_io::copy(Bhost, Bdev);
    std::vector<int> indices(N*3);
    
    std::vector<int> indices_validation(N*3);
    std::vector<device_real_t> rel_errs(N*3);
    std::vector<device_real_t> max_rel_errs(batch_size);
    std::cout << " Status" << int(Bhost.statuses[0]) << std::endl;
    
    for (size_t I = 0; I < batch_size; I++){
        auto offset = I*N*3;
        device_real_t* eigs_ = eigs.data + offset;
        device_real_t* eigs_validation_ = eigs_validation.data + offset;
        std::iota(indices.begin(), indices.end(), 0);
        std::iota(indices_validation.begin(), indices_validation.end(), 0); 
        std::sort(indices.begin(), indices.end(), [&](int i, int j){return eigs_[i] < eigs_[j];});
        std::sort(indices_validation.begin(), indices_validation.end(), [&](int i, int j){return eigs_validation_[i] < eigs_validation_[j];});
        std::transform(indices.begin() + 6, indices.end(), indices_validation.begin() + 6, rel_errs.begin() + 6, [&](int i, int j){return std::abs(eigs_[i]  - eigs_validation_[j])/ ( std::abs(eigs_[i])  + float(std::abs(eigs_[i]) < 1e-5) );});
        max_rel_errs[I] = *std::max_element(rel_errs.begin() + 6, rel_errs.end());
    }

    std::cout << "Max rel err: " << max_rel_errs << std::endl;

    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int i, int j){return eigs[i] < eigs[j];});
    std::sort(indices_validation.begin(), indices_validation.end(), [&](int i, int j){return eigs_validation[i] < eigs_validation[j];});
    int off0 = indices[6] * N*3; std::cout << eigs[indices[6]] << std::endl;
    int off1 = indices_validation[6] * N*3; std::cout << eigs_validation[indices_validation[6]] << std::endl;
    std::cout << min_eigvals[0] << std::endl;
    std::transform(indices.begin(), indices.end(), rel_errs.begin(), [&](int i){return std::abs(std::abs(Q[i+ off0]) - std::abs(min_eigvects[i]))/ ( std::abs(Q[i+ off0])  + float(std::abs(Q[i+ off0]) < 1e-5) );});
    std::cout << rel_errs << std::endl;
    std::transform(indices.begin(), indices.end(), rel_errs.begin(), [&](int i){return std::abs(std::abs(Q[i+ off0]) - std::abs(Q_validation[i + off1]))/ ( std::abs(Q[i+ off0])  + float(std::abs(Q[i+ off0]) < 1e-5) );});
    std::cout << rel_errs << std::endl;
    std::transform(indices.begin(), indices.end(), rel_errs.begin(), [&](int i){return std::abs(std::abs(Q_validation[i + off1]) - std::abs(min_eigvects[i]))/ ( std::abs(Q_validation[i + off1])  + float(std::abs(Q_validation[i + off1]) < 1e-5) );});
    std::cout << rel_errs << std::endl;


}
