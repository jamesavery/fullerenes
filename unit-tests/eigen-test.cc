#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/geometry.hh"
using namespace gpu_kernels;
int main(int argc, char** argv){
    auto N = 60;
    auto M = 1000;
    IsomerBatch batch1(N,M,DEVICE_BUFFER,0);
    IsomerBatch batch2(N,M,DEVICE_BUFFER,0);
    IsomerBatch h_batch(N,M,HOST_BUFFER,0);
    cuda_io::IsomerQueue queue1(N,0);
    BuckyGen::buckygen_queue  bucky_Q  = BuckyGen::start(N, false, false);
    Graph G;
    G.N = N;
    G.neighbours = neighbours_t(N);
    CuArray<float> output( M );

    for(int i  = 0; i < M; i++){
        BuckyGen::next_fullerene(bucky_Q, G);
        queue1.insert(G,0);
    }
    queue1.refill_batch(batch1);
    isomerspace_dual::dualise(batch1);
    isomerspace_tutte::tutte_layout(batch1);
    isomerspace_X0::zero_order_geometry(batch1, 4.0);
    //isomerspace_forcefield::test_fun(batch, output);
    cuda_io::reset_convergence_statuses(batch1);
    cuda_io::copy(batch2, batch1);
    isomerspace_forcefield::optimise<PEDERSEN>(batch1, N*5, N*5);
    isomerspace_forcefield::get_bond_rmse<PEDERSEN>(batch1, output);
    std::cout << output << endl;
    isomerspace_forcefield::get_angle_rmse<PEDERSEN>(batch1, output);
    std::cout << output << endl;

    isomerspace_forcefield::optimise<FLATNESS_ENABLED>(batch2, N*5, N*5);
    isomerspace_forcefield::get_bond_rmse<PEDERSEN>(batch2, output);
    std::cout << output << endl;
    isomerspace_forcefield::get_angle_rmse<PEDERSEN>(batch2, output);
    std::cout << output << endl;
    cuda_io::copy(h_batch, batch1);
    h_batch.print(STATUSES);
    //isomerspace_forcefield::test_fun(batch, output);
    auto P  = h_batch.get_isomer_by_id(0);
    Polyhedron::to_file(P.value(), "TestMol.mol2") ;
    LaunchCtx::clear_allocations();

}