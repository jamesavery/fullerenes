#include <sys/stat.h>
#include <limits.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <thread>
#include <future>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/isomerspace_tutte.hh"

#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/kernels.hh"

#include "numeric"

using namespace std;
using namespace std::chrono;

#define TEST_SAMPLES 2048


int main(int ac, char **argv)
{
    int N = 60;
    IsomerspaceTutte tutte_kernel   = IsomerspaceTutte(N);
    
    BuckyGen::buckygen_queue Q = BuckyGen::start(N,false,false);
    
    IsomerBatch d_in(N,TEST_SAMPLES,DEVICE_BUFFER);
    IsomerBatch d_opt(N,TEST_SAMPLES,DEVICE_BUFFER);
    IsomerBatch h_out(N,10,HOST_BUFFER);
    IsomerBatch h_out_small(N,10,HOST_BUFFER);
    FullereneDual dualG;
    cuda_io::IsomerQueue batch_queue(N);
    cuda_io::IsomerQueue batch_queue_2(N);

    auto filler = [&](int samples){
    for (size_t i = 0; i < samples; i++)
    {
        static int I = 0;
        static bool more_to_generate = true;
        if (!more_to_generate){break;}
        more_to_generate &= BuckyGen::next_fullerene(Q,dualG);
        if (!more_to_generate){break;}
        dualG.update();
        FullereneGraph G =  dualG.dual_graph();
        batch_queue.insert(G, I, LaunchCtx(), LaunchPolicy::SYNC, false);
        I++;
    }};
    std::queue<std::tuple<Polyhedron,size_t,IsomerStatus>> out_queue;
    filler(TEST_SAMPLES);
    batch_queue.refill_batch(d_in);
    cuda_io::copy(h_out_small,d_in, LaunchCtx(), LaunchPolicy::SYNC, {0,10}, {10,20});
    cuda_io::copy(h_out,d_in, LaunchCtx(), LaunchPolicy::SYNC, {0,10}, {10,20});

    assert(h_out == h_out_small);
    std::cout << h_out_small;

}