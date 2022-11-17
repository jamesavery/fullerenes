
#include <chrono>
#include <fstream>
#include <stdio.h>

#include "fullerenes/gpu/batch_queue.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/benchmark_functions.hh"
using namespace cuda_benchmark;
using namespace chrono;
using namespace chrono_literals;

int main(int argc, char** argv){
    size_t N_start                  = argc > 1 ? strtol(argv[1],0,0) : (size_t)1;   // Argument 1: Number of vertices N
    size_t N_end                    = argc > 2 ? strtol(argv[2],0,0) : (size_t)1024;   // Argument 1: Number of vertices N
    size_t N_iter                   = argc > 3 ? strtol(argv[3],0,0) : 100;     // Argument 1: Number of vertices N

    auto n_samples = 5;

    auto mean = [&](std::vector<std::chrono::nanoseconds> &input, const int samples){
            auto result = std::chrono::nanoseconds(0);
            for (int i = 0; i < samples; i++){
                result += input[i];
            }
            return result / samples;
        };

    auto standard_deviation = [&](std::vector<std::chrono::nanoseconds> &input, const int samples){
        auto mean_ = mean(input, samples);
        auto result = std::chrono::nanoseconds(0);
        for (int i = 0; i < samples; i++){
            result += (input[i] - mean_) *  (input[i] - mean_).count();
        }
        return std::chrono::nanoseconds((int)std::sqrt( (result / samples).count()));
    };

    std::vector<std::chrono::nanoseconds> 
            T_seq(n_samples),
            T_interleave(n_samples),
            T_noconfl(n_samples),
            T_warp(n_samples),
            T_warp_manual(n_samples),
            T_warp_atomic(n_samples),
            T_warp_other_out(n_samples),
            T_warp_other_full(n_samples);

    ofstream out_file("reduction_benchmark_" + to_string(N_start) + "_" + to_string(N_end) + ".txt");
    ofstream std_file("reduction_benchmark_std_" + to_string(N_start) + "_" + to_string(N_end) + ".txt");


    for (int i = N_start; i < N_end + 1; i++){
        unsigned int actual_blocksize = max(64, ( ( (i -1) >> 5 ) + 1) * 32);
        size_t N_blocks = n_blocks(i);
        for (int j = 0; j < n_samples; j++){
            T_seq[j] = benchmark_reduction(i,N_iter,0);
            T_interleave[j] = benchmark_reduction(i,N_iter,1);
            T_noconfl[j] = benchmark_reduction(i,N_iter,2);
            T_warp[j] = benchmark_reduction(i,N_iter,3);
            T_warp_manual[j] = benchmark_reduction(i,N_iter,4);
            T_warp_atomic[j] = benchmark_reduction(i,N_iter,5);
            T_warp_other_out[j] = benchmark_reduction(i,N_iter,6);
            T_warp_other_full[j] = benchmark_reduction(i,N_iter,7);
        };
        out_file <<  i << "," << N_blocks << "," << mean(T_seq,n_samples)/1ns <<  "," << mean(T_interleave,n_samples)/1ns << "," << mean(T_noconfl,n_samples)/1ns << "," << mean(T_warp,n_samples)/1ns << "," << mean(T_warp_manual,n_samples)/1ns << "," << mean(T_warp_atomic,n_samples)/1ns << "," << mean(T_warp_other_out,n_samples)/1ns << "," << mean(T_warp_other_full,n_samples)/1ns << "\n";
        std_file <<  i << "," << N_blocks << "," << standard_deviation(T_seq,n_samples)/1ns <<  "," << standard_deviation(T_interleave,n_samples)/1ns << "," << standard_deviation(T_noconfl,n_samples)/1ns << "," << standard_deviation(T_warp,n_samples)/1ns << "," << standard_deviation(T_warp_manual,n_samples)/1ns << "," << standard_deviation(T_warp_atomic, n_samples)/1ns << "," << standard_deviation(T_warp_other_out,n_samples)/1ns << "," << standard_deviation(T_warp_other_full,n_samples)/1ns <<"\n";
    }

}