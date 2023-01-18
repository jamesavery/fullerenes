
#include <chrono>
#include <fstream>
#include <stdio.h>

#include "fullerenes/gpu/isomer_queue.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
using namespace cuda_benchmark;
using namespace chrono;
using namespace chrono_literals;

int main(int argc, char** argv){
    size_t N_start                = argc > 1 ? strtol(argv[1],0,0) : (size_t)1;   // Argument 1: Number of vertices N
    size_t N_end                = argc > 2 ? strtol(argv[2],0,0) : (size_t)1024;   // Argument 1: Number of vertices N
    size_t N_iter                = argc > 3 ? strtol(argv[3],0,0) : 1000;     // Argument 1: Number of vertices N

    auto n_samples = 1;
    for (int i = N_start; i < N_end + 1; i++){
        
        bool test = test_global_scan(i,N_iter);
        std::string state = test? "Passed!" : "Failed!";
        std::cout << state  << std::endl;
    }

}