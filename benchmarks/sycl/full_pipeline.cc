#include <fullerenes/graph.hh>
#include <fullerenes/sycl-kernels.hh>
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

    size_t N  = argc>1 ? strtol(argv[1],0,0) : 20;
    size_t Nf = N/2 + 2;
    size_t BatchSize = argc>2 ? strtol(argv[2],0,0) : 1;
    std::string device_type = argc>3 ? argv[3] : "gpu";
    size_t Nruns = argc>4 ? strtol(argv[4],0,0) : 10;

    auto selector =  device_type == "cpu" ? sycl::cpu_selector_v : sycl::gpu_selector_v;

    sycl::queue Q = sycl::queue(selector, sycl::property::queue::in_order{});
    
    IsomerBatch<real_t,node_t> batch(N, BatchSize);
    Graph G(N);
    auto fill = [&](IsomerBatch<real_t,node_t>& batch)
    {
    BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, 0, 0);

    sycl::host_accessor acc_dual(batch.dual_neighbours, sycl::write_only);
    sycl::host_accessor acc_degs(batch.face_degrees, sycl::write_only);
    sycl::host_accessor acc_status (batch.statuses, sycl::write_only);
    for (size_t ii = 0; ii < BatchSize; ii++)
    {
        auto more = BuckyGen::next_fullerene(BuckyQ, G);
        if(!more) break;

        for (size_t j = 0; j < Nf; j++)
        {
            for(size_t k = 0; k < G.neighbours[j].size(); k++)
            {
                acc_dual[ii*Nf*6 + j*6 + k] = G.neighbours[j][k];
            } 
            if(G.neighbours[j].size() == 5){
                acc_dual[ii*Nf*6 + j*6 + 5] = std::numeric_limits<node_t>::max();
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
    vector<double> times_tutte(Nruns); //Times in nanoseconds.
    vector<double> times_project(Nruns); //Times in nanoseconds.
    vector<double> times_opt(Nruns); //Times in nanoseconds.

    for(int i = 0; i < Nruns; i++){
        auto start = std::chrono::steady_clock::now();
        fill(batch);
        auto T0 = std::chrono::steady_clock::now(); times_generate[i] = std::chrono::duration<double, std::nano>(T0 - start).count();
        nop_kernel(Q, batch, LaunchPolicy::SYNC);
        auto T1 = std::chrono::steady_clock::now(); times_memcpy[i] = std::chrono::duration<double, std::nano>(T1 - T0).count();
        dualise(Q, batch, LaunchPolicy::SYNC);
        auto T2 = std::chrono::steady_clock::now(); times_dual[i] = std::chrono::duration<double, std::nano>(T2 - T1).count();
        tutte_layout(Q, batch, LaunchPolicy::SYNC);
        auto T3 = std::chrono::steady_clock::now(); times_tutte[i] = std::chrono::duration<double, std::nano>(T3 - T2).count();
        spherical_projection(Q, batch, LaunchPolicy::SYNC);
        auto T4 = std::chrono::steady_clock::now(); times_project[i] = std::chrono::duration<double, std::nano>(T4 - T3).count();
        forcefield_optimise(Q, batch, 5*N, 5*N, LaunchPolicy::SYNC);
        auto T5 = std::chrono::steady_clock::now(); times_opt[i] = std::chrono::duration<double, std::nano>(T5 - T4).count();
        


    }

    std::cout << "N, Nf, BatchSize, device_type, generate, memcpy, dual, tutte, project, opt" << std::endl;
    std::cout << "N: " << N << ", Nf: " << Nf << ", BatchSize: " << BatchSize << ", device_type: " << device_type << "\n";
    std::cout << "Generate: " << mean(times_generate)/BatchSize << ", " << stddev(times_generate)/BatchSize << " ns \n";
    std::cout << "Memcpy: " << mean(times_memcpy)/BatchSize << ", " << stddev(times_memcpy)/BatchSize << " ns \n";
    std::cout << "Dual: " << mean(times_dual)/BatchSize << ", " << stddev(times_dual)/BatchSize << " ns \n";
    std::cout << "Tutte: " << mean(times_tutte)/BatchSize << ", " << stddev(times_tutte)/BatchSize << " ns \n";
    std::cout << "Project: " << mean(times_project)/BatchSize << ", " << stddev(times_project)/BatchSize << " ns \n";
    std::cout << "Opt: " << mean(times_opt)/BatchSize << ", " << stddev(times_opt)/BatchSize << " ns \n";

    return 0;
}
