#include <fullerenes/graph.hh>
#include "fullerenes/polyhedron.hh"
#include <fullerenes/sycl-headers/all-kernels.hh>
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <string>
#include <algorithm>

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

    size_t BatchSize = 1280; //std::ceil((real_t)NumNodes/(real_t)N);

    //auto selector =  device_type == "cpu" ? sycl::cpu_selector_v : sycl::gpu_selector_v;
    auto device = Device(device_type);

    //sycl::queue Q = sycl::queue(selector, sycl::property::queue::in_order{});
    SyclQueue Q(device);
    
    FullereneBatch<real_t,node_t> batch(N, BatchSize);
    Graph G(Nf, true);
    auto fill_and_dualize = [&](FullereneBatch<real_t,node_t>& batch, double& filltime, double& dualtime)
    {
    BuckyGen::buckygen_queue BuckyQ = BuckyGen::start(N, 0, 0);

    double ftime = 0; double dtime = 0;
    for (size_t ii = 0; ii < BatchSize; ii++)
    {   
        auto start = std::chrono::steady_clock::now();
        auto more = BuckyGen::next_fullerene(BuckyQ, G);
        //batch.push_back(G, ii);
        FullereneDual FD(G);
        auto T0 = std::chrono::steady_clock::now(); ftime += std::chrono::duration<double, std::nano>(T0 - start).count();
        FD.update();
        PlanarGraph pG = FD.dual_graph();
        batch.push_back(pG, ii);
        auto T1 = std::chrono::steady_clock::now(); dtime += std::chrono::duration<double, std::nano>(T1 - T0).count();
        if(!more) break;

        
        batch.m_.flags_[ii] |= StatusFlag::FULLERENEGRAPH_PREPARED;
    }
    BuckyGen::stop(BuckyQ);
    filltime = ftime;
    dualtime = dtime;
    };
    
    vector<double> times_generate(Nruns); //Times in nanoseconds.
    vector<double> times_memcpy(Nruns); //Times in nanoseconds.
    vector<double> times_dual(Nruns); //Times in nanoseconds.
    vector<double> times_tutte(Nruns); //Times in nanoseconds.
    vector<double> times_project(Nruns); //Times in nanoseconds.
    vector<double> times_opt(Nruns); //Times in nanoseconds.

    TutteFunctor<real_t,node_t> tutte_layout;
    SphericalProjectionFunctor<real_t,node_t> spherical_projection;
    ForcefieldOptimizeFunctor<PEDERSEN,real_t,node_t> forcefield_optimize;


    for(int i = 0; i < Nruns; i++){
        auto start = std::chrono::steady_clock::now();
        FullereneBatch<real_t,node_t> batch(N, BatchSize);
        fill_and_dualize(batch, times_generate[i], times_dual[i]);
        auto T0 = std::chrono::steady_clock::now();
        //nop_kernel(Q, batch, LaunchPolicy::SYNC);
        auto T1 = std::chrono::steady_clock::now(); times_memcpy[i] = std::chrono::duration<double, std::nano>(T1 - T0).count();
        tutte_layout(Q, batch, LaunchPolicy::SYNC);
        auto T2 = std::chrono::steady_clock::now(); times_tutte[i] = std::chrono::duration<double, std::nano>(T2 - T1).count();
        spherical_projection(Q, batch, LaunchPolicy::SYNC);
        auto T3 = std::chrono::steady_clock::now(); times_project[i] = std::chrono::duration<double, std::nano>(T3 - T2).count();
        forcefield_optimize(Q, batch, LaunchPolicy::SYNC, 5*N, 5*N);
        auto T4 = std::chrono::steady_clock::now(); times_opt[i] = std::chrono::duration<double, std::nano>(T4 - T3).count();
        


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
