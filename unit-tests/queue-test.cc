#include "numeric"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/gpu/benchmark_functions.hh"


int main(int argc, char** argv) {
    auto iters = 10;
    const size_t N                =   argc > 1 ? strtol(argv[1],0,0) : 20;     // Argument 1: Number of vertices 
    using namespace gpu_kernels;
    auto batch_size = 200;
    std::cout << "Batch size: " << batch_size << std::endl;
    cuda_io::IsomerQueue Q0(N, 0); Q0.resize(5*batch_size);
    cuda_io::IsomerQueue Q1(N, 0); 
    cuda_io::IsomerQueue Q2(N, 0); Q2.resize(5*batch_size);
    auto Nf = N/2 + 2;
    auto ID = 0;
    Graph G(Nf);
    IsomerBatch<GPU> batch(N, batch_size);
    IsomerBatch<GPU> batch1(N, batch_size);
    auto generate_samples = [&](int samples) {
        for(auto i = 0; i < samples; i++) {
            cuda_benchmark::random_isomer("isomerspace_samples/dual_layout_" + std::to_string(N) + "_seed_42", G);
            Q0.insert(G,ID);
            ID++;
        }
        Q0.refill_batch(batch);
        isomerspace_dual::dualise(batch);
        isomerspace_tutte::tutte_layout(batch);
        isomerspace_X0::zero_order_geometry(batch, 4.0);
        Q1.insert(batch);
    };
    while(Q2.get_size() < 1000){
        for(int i = 0; i < iters; i++) {
            generate_samples(batch_size*0.5);
            Q1.refill_batch(batch1);
            isomerspace_forcefield::optimise<PEDERSEN>(batch1, N*1,N*6);
            Q2.push_done(batch1);
        }
    }
    while(Q1.get_size()>0) {
        Q1.refill_batch(batch1);
        isomerspace_forcefield::optimise<PEDERSEN>(batch1, N*6,N*6);
        Q2.push_done(batch1);
    }

    IsomerBatch<CPU> host_batch(N, Q2.get_capacity());
    cuda_io::copy(host_batch,Q2.device_batch);
    host_batch.shrink_to_fit();
    host_batch.print(BatchMember::ITERATIONS);
    cuda_io::sort(host_batch, BatchMember::IDS, SortOrder::ASCENDING);
    std::vector<size_t> ids(host_batch.IDs, host_batch.IDs + host_batch.capacity());
    std::vector<size_t> id_valid(host_batch.capacity(),0);
    std::iota(id_valid.begin(), id_valid.end(), 0);
    std::cout << "IDs: " << ids << std::endl;
    assert(std::equal(ids.begin(), ids.end(), id_valid.begin()));

}
