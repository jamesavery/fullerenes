
#include <chrono>
#include <fstream>
#include <stdio.h>
#include "fullerenes/isomer_queue.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/device_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
using namespace cuda_benchmark;
using namespace chrono;
using namespace chrono_literals;

int main(int argc, char** argv){
    size_t N_start                  = argc > 1 ? strtol(argv[1],0,0) : (size_t)1;       // Argument 1: Number of vertices N
    size_t N_end                    = argc > 2 ? strtol(argv[2],0,0) : (size_t)1024;    // Argument 2: Number of vertices N
    size_t N_iter                   = argc > 3 ? strtol(argv[3],0,0) : 100;             // Argument 3: Number of vertices N

    auto n_samples = 10;

    std::vector<std::chrono::nanoseconds> 
            T_seq(n_samples),
            T_blelloch(n_samples),
            T_warp(n_samples);

    ofstream out_file("scan_benchmark_" + to_string(N_start) + "_" + to_string(N_end) + ".txt");
    ofstream std_file("scan_benchmark_std_" + to_string(N_start) + "_" + to_string(N_end) + ".txt");

    using namespace device_io;
    for (int i = N_start; i < N_end + 1; i++){
        
        unsigned int actual_blocksize = max(64, ( ( (i -1) >> 5 ) + 1) * 32);
        for (int j = 0; j < n_samples; j++){
            T_seq[j] = chrono::nanoseconds( benchmark_scan(i,N_iter,0));
            T_blelloch[j] = chrono::nanoseconds( benchmark_scan(i,N_iter,1));
            T_warp[j] = chrono::nanoseconds( benchmark_scan(i,N_iter,2));
        };
        out_file <<  i << "," << 68* (1024/actual_blocksize) << "," << mean(T_seq)/1us <<  "," << mean(T_blelloch)/1us << "," << mean(T_warp)/1us << "\n";
        std_file <<  i << "," << 68* (1024/actual_blocksize) << "," << sdev(T_seq)/1us <<  "," << sdev(T_blelloch)/1us << "," << sdev(T_warp)/1us << "\n";
    }

}
