#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/isomer_queue.hh"

#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/gpu/cu_array.hh"
#include <numeric>
#include <random>

using namespace gpu_kernels;
using namespace device_io;
#define SYNC LaunchPolicy::SYNC
#define ASYNC LaunchPolicy::ASYNC

int main(int ac, char **argv){
    
    int N                = strtol(argv[1],0,0);     // Argument 1: Number of vertices N
    auto n_isomers = IsomerDB::number_isomers(N);
    Graph G;
    auto Nf = N/2 +2;
    G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));

    std::string path = "isomerspace_samples/dual_layout_" + to_string(N) + "_seed_42";
    auto fsize = file_size(path);
    auto n_samples = fsize / (Nf * 6 * sizeof(device_node_t));
    ifstream in_file(path,std::ios::binary);
    std::vector<device_node_t> dual_neighbours(n_samples * Nf * 6);
    in_file.read((char*)dual_neighbours.data(),n_samples * Nf * 6 * sizeof(device_node_t));
    
    int sample_size = min((int)n_samples,300);

    std::vector<int> random_IDs(n_samples);
    std::iota(random_IDs.begin(), random_IDs.end(), 0);
    std::shuffle(random_IDs.begin(), random_IDs.end(), std::mt19937{42});
    auto id_range_end = min((int)n_samples, sample_size);
    std::vector<int> id_subset(random_IDs.begin(), random_IDs.begin()+id_range_end);
    
    auto SEQUENTIAL_INSERT_AND_FILL_TEST = [&](){
        IsomerQueue TestQ(N,0);
        IsomerQueue ControlQ(N,0);
        IsomerBatch<GPU> d_TestB(N,sample_size, 0);
        IsomerBatch<GPU> d_ControlB(N,sample_size, 0);

        IsomerBatch<CPU> h_TestB(N,sample_size);
        IsomerBatch<CPU> h_ControlB(N,sample_size);

        LaunchCtx test_ctx(0);
        LaunchCtx control_ctx(0);

        G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
        G.N = Nf;
        
        for (int i = 0; i < sample_size; ++i){
            for (size_t j = 0; j < Nf; j++){
                G.neighbours[j].clear();
                for (size_t k = 0; k < 6; k++) {
                    auto u = dual_neighbours[id_subset[i]*Nf*6 + j*6 +k];
                    if(u != UINT16_MAX) G.neighbours[j].push_back(u);
                }
            }
            //Insert isomers in both the test and control queue.
            TestQ.insert(G,i, test_ctx, ASYNC);
            ControlQ.insert(G,i, control_ctx, ASYNC);
        }
        //Refill batches 
        TestQ.refill_batch(d_TestB, test_ctx, ASYNC);
        ControlQ.refill_batch(d_ControlB, control_ctx, ASYNC);

        //Copy to host to test equality
        copy(h_TestB, d_TestB, test_ctx, ASYNC);
        copy(h_ControlB, d_ControlB, control_ctx, ASYNC);
        

        //Wait for operations on streams to complete.
        test_ctx.wait(); control_ctx.wait();

        //Test equality.
        bool pass = h_TestB == h_ControlB;
        if (!pass) {
            //std::cout << h_TestB;
            //std::cout << h_ControlB;
        }
        return pass;
    };

    auto BATCH_INSERT_AND_FILL_TEST = [&](){
        IsomerQueue TestQ(N,0);
        IsomerQueue TestQ2(N,0);
        IsomerQueue ControlQ(N,0);
        IsomerQueue ControlQ2(N,0);
        IsomerBatch<GPU> d_TestB(N,sample_size, 0);
        IsomerBatch<GPU> d_TestB2(N,sample_size, 0);
        IsomerBatch<GPU> d_ControlB(N,sample_size, 0);
        IsomerBatch<GPU> d_ControlB2(N,sample_size, 0);

        IsomerBatch<CPU> h_TestB(N,sample_size, 0);
        IsomerBatch<CPU> h_ControlB(N,sample_size, 0);

        LaunchCtx test_ctx(0);
        LaunchCtx control_ctx(0);

        G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
        G.N = Nf;
        
        for (int i = 0; i < sample_size; ++i){
            for (size_t j = 0; j < Nf; j++){
                G.neighbours[j].clear();
                for (size_t k = 0; k < 6; k++) {
                    auto u = dual_neighbours[id_subset[i]*Nf*6 + j*6 +k];
                    if(u != UINT16_MAX) G.neighbours[j].push_back(u);
                }
            }
            //Insert isomers in both the test and control queue.
            TestQ.insert(G,i, test_ctx, ASYNC);
            ControlQ.insert(G,i, control_ctx, ASYNC);
        }
        //Refill batches 
        TestQ.refill_batch(d_TestB, test_ctx, ASYNC);
        ControlQ.refill_batch(d_ControlB, control_ctx, ASYNC);

        TestQ2.insert(d_TestB, test_ctx, ASYNC);
        ControlQ2.insert(d_ControlB, control_ctx, ASYNC);

        TestQ2.refill_batch(d_TestB2, test_ctx, ASYNC);
        ControlQ2.refill_batch(d_ControlB2, control_ctx, ASYNC);
        
        //Copy to host to test equality
        copy(h_TestB, d_TestB2, test_ctx, ASYNC);
        copy(h_ControlB, d_ControlB2, control_ctx, ASYNC);
        //Wait for operations on streams to complete.
        test_ctx.wait(); control_ctx.wait();


        //Test equality.
        bool pass = h_TestB == h_ControlB;
        if (!pass) {
            //std::cout << h_TestB;
        }
        return pass;
    };

    auto CUBIC_TEST = [&](){
        IsomerQueue TestQ(N,0);
        IsomerQueue TestQ2(N,0);
        IsomerQueue ControlQ(N,0);
        IsomerQueue ControlQ2(N,0);
        IsomerBatch<GPU> d_TestB(N,sample_size, 0);
        IsomerBatch<GPU> d_TestB2(N,sample_size, 0);
        IsomerBatch<GPU> d_ControlB(N,sample_size, 0);
        IsomerBatch<GPU> d_ControlB2(N,sample_size, 0);

        IsomerBatch<CPU> h_TestB(N,sample_size, 0);
        IsomerBatch<CPU> h_ControlB(N,sample_size, 0);

        LaunchCtx test_ctx(0);
        LaunchCtx control_ctx(0);

        G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
        G.N = Nf;
        
        for (int i = 0; i < sample_size; ++i){
            for (size_t j = 0; j < Nf; j++){
                G.neighbours[j].clear();
                for (size_t k = 0; k < 6; k++) {
                    auto u = dual_neighbours[id_subset[i]*Nf*6 + j*6 +k];
                    if(u != UINT16_MAX) G.neighbours[j].push_back(u);
                }
            }
            //Insert isomers in both the test and control queue.
            TestQ.insert(G,i, test_ctx, ASYNC);
            ControlQ.insert(G,i, control_ctx, ASYNC);
        }
        //Refill batches 
        TestQ.refill_batch(d_TestB, test_ctx, ASYNC);
        ControlQ.refill_batch(d_ControlB, control_ctx, ASYNC);

        isomerspace_dual::dualize(d_TestB, test_ctx, ASYNC);
        isomerspace_tutte::tutte_layout(d_TestB, 10000000, test_ctx, ASYNC);
        isomerspace_X0::zero_order_geometry(d_TestB, 4.0, test_ctx, ASYNC);
        isomerspace_dual::dualize(d_ControlB, control_ctx, ASYNC);
        isomerspace_tutte::tutte_layout(d_ControlB, 10000000, control_ctx, ASYNC);
        isomerspace_X0::zero_order_geometry(d_ControlB, 4.0, control_ctx, ASYNC);

        TestQ2.insert(d_TestB, test_ctx, ASYNC);
        ControlQ2.insert(d_ControlB, control_ctx, ASYNC);

        TestQ2.refill_batch(d_TestB2, test_ctx, ASYNC);
        ControlQ2.refill_batch(d_ControlB2, control_ctx, ASYNC);
        
        //Copy to host to test equality
        copy(h_TestB, d_TestB2, test_ctx, ASYNC);
        copy(h_ControlB, d_ControlB2, control_ctx, ASYNC);
        //Wait for operations on streams to complete.
        test_ctx.wait(); control_ctx.wait();


        //Test equality.
        bool pass = h_TestB == h_ControlB;
        if (!pass) {
            //std::cout << h_TestB;
        }
        return pass;
    };

    auto OPTIMIZE_TEST = [&](){
        IsomerQueue TestQ(N,0);
        IsomerQueue TestQ2(N,0);
        IsomerQueue ControlQ(N,0);
        IsomerQueue ControlQ2(N,0);
        IsomerBatch<GPU> d_TestB(N,sample_size*2, 0);
        IsomerBatch<GPU> d_TestB2(N,sample_size, 0);
        IsomerBatch<GPU> d_ControlB(N,sample_size*2, 0);
        IsomerBatch<GPU> d_ControlB2(N,sample_size, 0);

        IsomerBatch<CPU> h_TestB(N,sample_size, 0);
        IsomerBatch<CPU> h_ControlB(N,sample_size, 0);

        LaunchCtx test_ctx(0);
        LaunchCtx control_ctx(0);

        G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
        G.N = Nf;
        
        for (int i = 0; i < sample_size; ++i){
            for (size_t j = 0; j < Nf; j++){
                G.neighbours[j].clear();
                for (size_t k = 0; k < 6; k++) {
                    auto u = dual_neighbours[id_subset[i]*Nf*6 + j*6 +k];
                    if(u != UINT16_MAX) G.neighbours[j].push_back(u);
                }
            }
            //Insert isomers in both the test and control queue.
            TestQ.insert(G,i*2, test_ctx, ASYNC);
            TestQ.insert(G,i*2+1, test_ctx, ASYNC);
            ControlQ.insert(G,i*2, control_ctx, ASYNC);
            ControlQ.insert(G,i*2+1, control_ctx, ASYNC);
        }
        //Refill batches 
        TestQ.refill_batch(d_TestB, test_ctx, ASYNC);
        ControlQ.refill_batch(d_ControlB, control_ctx, ASYNC);

        isomerspace_dual::dualize(d_TestB, test_ctx, ASYNC);
        isomerspace_tutte::tutte_layout(d_TestB, 10000000, test_ctx, ASYNC);
        isomerspace_X0::zero_order_geometry(d_TestB, 4.0, test_ctx, ASYNC);
        reset_convergence_statuses(d_TestB, test_ctx, ASYNC);
        isomerspace_dual::dualize(d_ControlB, control_ctx, ASYNC);
        isomerspace_tutte::tutte_layout(d_ControlB, 10000000, control_ctx, ASYNC);
        isomerspace_X0::zero_order_geometry(d_ControlB, 4.0, control_ctx, ASYNC);
        reset_convergence_statuses(d_ControlB, control_ctx, ASYNC);
        TestQ2.insert(d_TestB, test_ctx, ASYNC);
        ControlQ2.insert(d_ControlB, control_ctx, ASYNC);

        
        TestQ2.refill_batch(d_TestB2, test_ctx, ASYNC);
        ControlQ2.refill_batch(d_ControlB2, control_ctx, ASYNC);
        isomerspace_forcefield::optimize<PEDERSEN>(d_TestB2, N, N*2,test_ctx, ASYNC);
        isomerspace_forcefield::optimize<PEDERSEN>(d_ControlB2, N, N*2,control_ctx, ASYNC);
        TestQ2.refill_batch(d_TestB2, test_ctx, ASYNC);
        ControlQ2.refill_batch(d_ControlB2, control_ctx, ASYNC);
        TestQ2.insert(d_TestB, test_ctx, ASYNC);
        ControlQ2.insert(d_ControlB, control_ctx, ASYNC);
        isomerspace_forcefield::optimize<PEDERSEN>(d_TestB2, N, N*2,test_ctx, ASYNC);
        isomerspace_forcefield::optimize<PEDERSEN>(d_ControlB2, N, N*2,control_ctx, ASYNC);
        TestQ2.refill_batch(d_TestB2, test_ctx, ASYNC);
        ControlQ2.refill_batch(d_ControlB2, control_ctx, ASYNC);

        //Copy to host to test equality
        copy(h_TestB, d_TestB2, test_ctx, ASYNC);
        copy(h_ControlB, d_ControlB2, control_ctx, ASYNC);
        //Wait for operations on streams to complete.
        test_ctx.wait(); control_ctx.wait();


        //Test equality.
        bool pass = h_TestB == h_ControlB;
        if (!pass) {
        }
        return pass;
    };

    auto ORDERED_TEST = [&](){
        IsomerQueue TestQ(N,0);
        IsomerQueue TestQ2(N,0);
        IsomerBatch<GPU> d_TestB(N,sample_size*2, 0);
        IsomerBatch<GPU> d_TestB2(N,sample_size, 0);

        LaunchCtx test_ctx(0);

        G.neighbours = neighbours_t(Nf, std::vector<node_t>(6));
        G.N = Nf;
        
        for (int i = 0; i < sample_size; ++i){
            for (size_t j = 0; j < Nf; j++){
                G.neighbours[j].clear();
                for (size_t k = 0; k < 6; k++) {
                    auto u = dual_neighbours[id_subset[i]*Nf*6 + j*6 +k];
                    if(u != UINT16_MAX) G.neighbours[j].push_back(u);
                }
            }
            //Insert isomers in both the test and control queue.
            TestQ.insert(G,i*2, test_ctx, ASYNC);
            TestQ.insert(G,i*2+1, test_ctx, ASYNC);
        }
        //Refill batches 
        std::cout << "Counters: "<< TestQ.get_front()  << " , "<< TestQ.get_back() << std::endl;
        TestQ.refill_batch(d_TestB, test_ctx, SYNC);
        std::cout << "Counters: "<< TestQ.get_front()  << " , "<< TestQ.get_back() << std::endl;

        isomerspace_dual::dualize(d_TestB, test_ctx, SYNC);
        isomerspace_tutte::tutte_layout(d_TestB, 10000000, test_ctx, ASYNC);
        isomerspace_X0::zero_order_geometry(d_TestB, 4.0, test_ctx, ASYNC);
        reset_convergence_statuses(d_TestB, test_ctx, SYNC);

        //std::cout << d_TestB << std::endl;
        //std::cout << TestQ2.device_batch << std::endl;
        std::cout << "===============================================================" << endl;
        TestQ2.insert(d_TestB, test_ctx, SYNC);
        //std::cout << TestQ2.device_batch << std::endl;
        std::cout << "===============================================================" << endl;
        TestQ2.refill_batch(d_TestB2, test_ctx, SYNC);
        //std::cout << d_TestB2 << std::endl;
        std::cout << "===============================================================" << endl;

        isomerspace_forcefield::optimize<PEDERSEN>(d_TestB2, N*3, N*4,test_ctx, SYNC);
        //std::cout<< d_TestB2 << endl;
        std::cout << "===============================================================" << endl;

        TestQ2.refill_batch(d_TestB2, test_ctx, SYNC);
        //std::cout<< d_TestB2 << endl;
        std::cout << "===============================================================" << endl;



        isomerspace_forcefield::optimize<PEDERSEN>(d_TestB2, N*1, N*4,test_ctx, SYNC);
        TestQ2.refill_batch(d_TestB2, test_ctx, SYNC);
        isomerspace_forcefield::optimize<PEDERSEN>(d_TestB2, N*4, N*4,test_ctx, SYNC);
        TestQ2.refill_batch(d_TestB2, test_ctx, SYNC);

        //std::cout << d_TestB2 << endl;

        std::vector<size_t> id_list(n_samples);
        std::iota(id_list.begin(), id_list.end(), n_samples);

        //Wait for operations on streams to complete.
        test_ctx.wait();

        //Test equality.
        bool pass = true;

//        for(int i= 0; i < h_TestB.isomer_capacity; ++i){
//            pass = pass && id_list[i] == h_TestB.IDs[i];
//            if (!pass) {std::cout << "Failed at: " << id_list[i] << " , " << h_TestB.IDs[i] << std::endl; break;}
//        }
        //std::cout << h_TestB;
        return pass;
    };

    auto test_iterations = 5;
    std::array<bool,5> test_passed({true, true, true, true, true});
    for(int i = 0; i < test_iterations; ++i){test_passed[0] = test_passed[0] && SEQUENTIAL_INSERT_AND_FILL_TEST();} if(test_passed[0]) std::cout << "Test 0 Passed!" << std::endl; else std::cout << "Test 0 Failed!" << std::endl;
    for(int i = 0; i < test_iterations; ++i){test_passed[1] = test_passed[1] && BATCH_INSERT_AND_FILL_TEST();} if(test_passed[1]) std::cout << "Test 1 Passed!" << std::endl; else std::cout << "Test 1 Failed!" << std::endl;
    for(int i = 0; i < test_iterations; ++i){test_passed[2] = test_passed[2] && CUBIC_TEST();} if(test_passed[2]) std::cout << "Test 2 Passed!" << std::endl; else std::cout << "Test 2 Failed!" << std::endl;
    for(int i = 0; i < test_iterations; ++i){test_passed[3] = test_passed[3] && OPTIMIZE_TEST();} if(test_passed[3]) std::cout << "Test 3 Passed!" << std::endl; else std::cout << "Test 3 Failed!" << std::endl;
    for(int i = 0; i < 1; ++i){test_passed[4] = test_passed[4] && ORDERED_TEST();} if(test_passed[4]) std::cout << "Test 4 Passed!" << std::endl; else std::cout << "Test 4 Failed!" << std::endl;
    

}
