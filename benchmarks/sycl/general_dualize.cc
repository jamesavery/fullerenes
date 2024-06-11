#include <fullerenes/graph.hh>
#include <fullerenes/planargraph.hh>
#include <fullerenes/spiral.hh>
#include <fullerenes/sycl-kernels.hh>
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/argparser.hh>
#include <string>
#include <algorithm>
#include <fstream>
#include <fullerenes/triangulation.hh>
using namespace sycl;
    // Function to parse command line arguments

int main(int argc, char** argv) {
    CmdArgs args;
    parseArguments(argc, argv, args);
    typedef float real_t;
    typedef uint16_t node_t;

    size_t N = args.natoms;
    size_t BatchSize = args.nisomers;
    std::string device_type = args.device_type;
    bool use_double_precision = args.use_double_precision;
    size_t Nlanczos = args.nlanczos;
    std::string output_file = args.output_file;
    size_t print_vector = args.print_vector;
    size_t device_id = args.device_id;
    bool print_out = false;


    size_t Nf = N/2 + 2;


    auto selector =  device_type == "cpu" ? sycl::cpu_selector_v : sycl::gpu_selector_v;

    std::vector<sycl::device> devices;
    for (auto &device : sycl::device::get_devices())
    {
        if (device_type == "cpu" && device.is_cpu())
        {
            devices.push_back(device);
        }
        else if (device_type == "gpu" && device.is_gpu())
        {
            devices.push_back(device);
        }
    }
    sycl::queue Q = sycl::queue(devices[device_id]);

    Triangulation T(spiral_nomenclature("C49056=[4374,7762,8678,12888,17567,18028,18198,18827,20169,22581,22707,23966]-fullerene"));

    auto begin = std::chrono::high_resolution_clock::now();
    T.dual_graph();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Dual graph time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms" << std::endl;
    

    sycl::buffer<node_t> G_in = sycl::buffer<node_t>(T.N*6);
    sycl::buffer<node_t> Deg_in = sycl::buffer<node_t>(T.N);
    sycl::buffer<node_t> G_out = sycl::buffer<node_t>((T.N-2)*2);
    sycl::buffer<node_t> Deg_out = sycl::buffer<node_t>((T.N-2)*2);
    {
        host_accessor G_in_acc(G_in, write_only);
        host_accessor Deg_in_acc(Deg_in, write_only);

        for (int i = 0; i < T.N; i++) {
            Deg_in_acc[i] = T.neighbours[i].size();
            for (int j = 0; j < T.neighbours[i].size(); j++) {
                G_in_acc[i*6+j] = T.neighbours[i][j];
            }            
        }
    }
    


    dualize_general<6,3>(Q, G_in, Deg_in, G_out, Deg_out, T.N, (T.N-2)*2, LaunchPolicy::SYNC);
/* 
    neighbours_t neighbours((T.N-2)*2);
    {
        host_accessor G_out_acc(G_out, read_only);
        host_accessor Deg_out_acc(Deg_out, read_only);

        for (int i = 0; i < (T.N-2)*2; i++) {
            neighbours[i].resize(Deg_out_acc[i]);
            for (int j = 0; j < Deg_out_acc[i]; j++) {
                neighbours[i][j] = G_out_acc[i*6+j];
            }
        }
    }   

    PlanarGraph dual(Graph(neighbours,true)); */

}