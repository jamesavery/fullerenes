#include <fullerenes/graph.hh>
#include <fullerenes/sycl-headers/all-kernels.hh>
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <string>
#include <algorithm>
using namespace sycl;

template<typename T>
T stddev(const std::vector<T>& data)
{
  if(data.size() <= 1) return 0;
  
    // Calculate the mean
    T mn = mean(data);

    // Calculate the sum of squared differences from the mean
    T sum_of_squares = 0.0;
    for (const T& value : data)
    {
        T diff = value - mn;
        sum_of_squares += diff * diff;
    }

    // Calculate the variance and return the square root
    T variance = sum_of_squares / (data.size() - 1);
    return std::sqrt(variance);
}

int main(int argc, char** argv) {
    typedef float real_t;
    typedef uint16_t node_t;

    size_t N  = argc>1 ? strtol(argv[1],0,0) : 80;
    size_t Nf = N/2 + 2;
    size_t NumNodes = argc>2 ? strtol(argv[2],0,0) : 2000000;
    std::string device_type = argc>3 ? argv[3] : "gpu";
    size_t Nruns = argc>4 ? strtol(argv[4],0,0) : 10;
    size_t Nwarmup = argc>5 ? strtol(argv[5],0,0) : 1;
    auto dual_version = argc>6 ? strtol(argv[6],0,0) : 0;

    size_t BatchSize = std::ceil((real_t)NumNodes/(real_t)N);


    //auto selector =  device_type == "cpu" ? sycl::cpu_selector_v : sycl::gpu_selector_v;

    SyclQueue Q(device_type);
    
    FullereneBatch<real_t,node_t> batch(N, BatchSize);
    Graph G(N);
    auto fill = [&](FullereneBatch<real_t,node_t>& batch)
    {
    BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, 0, 0);

    auto acc_dual = batch.d_.A_dual_;
    auto acc_degs = batch.d_.deg_;
    auto acc_status = batch.m_.flags_;
    for (size_t ii = 0; ii < BatchSize; ii++)
    {
        auto more = BuckyGen::next_fullerene(BuckyQ, G);
        if(!more) break;

        for (size_t j = 0; j < Nf; j++)
        {
            for(size_t k = 0; k < G.neighbours[j].size(); k++)
            {
                acc_dual[ii*Nf + j][k] = G.neighbours[j][k];
            } 
            if(G.neighbours[j].size() == 5){
                acc_dual[ii*Nf + j][5] = std::numeric_limits<node_t>::max();
                acc_degs[ii*Nf + j] = 5;
            } else {
                acc_degs[ii*Nf + j] = 6;
            }   

        }
        acc_status[ii] = IsomerStatus::NOT_CONVERGED;
    }

    BuckyGen::stop(BuckyQ);
    };
    
    vector<double> times_generate(Nruns); //Times in nanoseconds.
    vector<double> times_memcpy(Nruns); //Times in nanoseconds.
    vector<double> times_dual(Nruns); //Times in nanoseconds.

    DualizeFunctor<real_t, node_t> dualize;

    for(int i = 0; i < Nruns + Nwarmup; i++){
        auto start = std::chrono::steady_clock::now();
        fill(batch);
        auto T0 = std::chrono::steady_clock::now(); times_generate[i] = std::chrono::duration<double, std::nano>(T0 - start).count();
        //nop_kernel(Q, batch, LaunchPolicy::SYNC);
        auto T1 = std::chrono::steady_clock::now(); times_memcpy[i] = std::chrono::duration<double, std::nano>(T1 - T0).count();
        switch(dual_version){
            case 0:
                dualize(Q, batch, LaunchPolicy::SYNC);
                break;
            default:
                dualize(Q, batch, LaunchPolicy::SYNC);
                break;
        }
        if (i >= Nwarmup) {auto T2 = std::chrono::steady_clock::now(); times_dual[i- Nwarmup] = std::chrono::duration<double, std::nano>(T2 - T1).count();}
    }

    std::cout << "N, Nf, BatchSize, device_type" << std::endl;
    std::cout << "N: " << N << ", Nf: " << Nf << ", BatchSize: " << BatchSize << ", device_type: " << device_type << "\n";
    std::cout << "Memory copy: " << mean(times_memcpy)/BatchSize << ", " << stddev(times_memcpy)/BatchSize << " ns \n";
    std::cout << "Dual: " << mean(times_dual)/BatchSize << ", " << stddev(times_dual)/BatchSize << " ns \n";
    return 0;
}
