#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/batch_queue.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/geometry.hh"
using namespace gpu_kernels;
int main(int argc, char** argv){
    IsomerBatch batch(20,1,DEVICE_BUFFER,0);
    cuda_io::IsomerQueue queue(20,0);
    BuckyGen::buckygen_queue  bucky_Q  = BuckyGen::start(20, false, false);
    CuArray<device_real_t> output(9);
    isomerspace_forcefield::test_fun(batch,output);
    Graph G;
    G.N = 20;
    G.neighbours = neighbours_t(20);
    BuckyGen::next_fullerene(bucky_Q, G);
    queue.insert(G,0);
    queue.refill_batch(batch);
    isomerspace_dual::cubic_layout(batch);
    isomerspace_forcefield::test_fun(batch, output);
    matrix3d test_mat(1., 2., 3., 2., 4., 5., 3., 5., 6.);
    coord3d lambdas = test_mat.eigenvalues();
    std::cout << test_mat.eigenvector(lambdas[0]);
    std::cout << test_mat.eigenvector(lambdas[1]);
    std::cout << test_mat.eigenvector(lambdas[2]);
    std::cout << output;


}