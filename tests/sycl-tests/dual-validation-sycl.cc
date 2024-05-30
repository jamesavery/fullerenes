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
#include <fullerenes/isomerdb.hh>
using namespace sycl;
    // Function to parse command line arguments

int main(int argc, char** argv) {
    CmdArgs args;
    parseArguments(argc, argv, args);
    typedef float real_t;
    typedef uint16_t node_t;

    size_t N = args.natoms;
    size_t BatchSize = args.nisomers == 0 ? IsomerDB::number_isomers(N) : args.nisomers;
    std::string device_type = args.device_type;
    size_t device_id = args.device_id;

    IsomerBatch<real_t,node_t> batch(N, BatchSize);
    fill(batch);

    auto selector =  device_type == "cpu" ? sycl::cpu_selector_v : sycl::gpu_selector_v;
    size_t Nf = N/2 + 2;
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

    sycl::buffer<node_t,1> G_in(sycl::range<1>(Nf*6));
    sycl::buffer<node_t,1> Deg_in(sycl::range<1>(Nf*1));
    sycl::buffer<node_t,1> G_out(sycl::range<1>(N*3));
    sycl::buffer<node_t,1> Deg_out(sycl::range<1>(N*1));
    neighbours_t input_neighbours(Nf);
    neighbours_t output_neighbours(N);

    for (size_t i = 0; i < BatchSize; i++) {
        {
            host_accessor G_in_acc(G_in, read_write);
            host_accessor Deg_in_acc(Deg_in, read_write);
            host_accessor G_out_acc(G_out, read_write);
            host_accessor Deg_out_acc(Deg_out, read_write);
            host_accessor batch_Gdual_acc(batch.dual_neighbours, read_only);
            host_accessor batch_GdualDeg_acc(batch.face_degrees , read_only);
            for (size_t j = 0; j < Nf*6; j++) {
                G_in_acc[j] = batch_Gdual_acc[i*Nf*6 + j];
            }
            for (size_t j = 0; j < Nf; j++) {
                Deg_in_acc[j] = batch_GdualDeg_acc[i*Nf + j];
                if (Deg_in_acc[j] == 6){
                    input_neighbours[j] = {G_in_acc[j*6], G_in_acc[j*6 + 1], G_in_acc[j*6 + 2], G_in_acc[j*6 + 3], G_in_acc[j*6 + 4], G_in_acc[j*6 + 5]};
                } else {
                    input_neighbours[j] = {G_in_acc[j*6], G_in_acc[j*6 + 1], G_in_acc[j*6 + 2], G_in_acc[j*6 + 3], G_in_acc[j*6 + 4]};
                }
            }
        }
        dualize_general<6,3>(Q, G_in, Deg_in, G_out, Deg_out, Nf, N, LaunchPolicy::SYNC);
        {
            host_accessor G_out_acc(G_out, read_only);
            for (size_t j = 0; j < N; j++) {
                output_neighbours[j] = {G_out_acc[j*3], G_out_acc[j*3 + 1], G_out_acc[j*3 + 2]};
            }
        }        
        PlanarGraph pg(input_neighbours);
        PlanarGraph pg_dual(output_neighbours);
        spiral_nomenclature Sin(pg, spiral_nomenclature::CAGE, spiral_nomenclature::TRIANGULATION);
        spiral_nomenclature Sout(pg_dual, spiral_nomenclature::CAGE, spiral_nomenclature::CUBIC);
        std::string true_false = (Sin.to_string() == Sout.to_string()) ? "True" : "False";
        std::cout << "Isomer " << i << " Is isomorphic: " << true_false << std::endl;
        if (Sin.to_string() != Sout.to_string()) {
            std::cout << "Input: " << Sin.to_string() << std::endl;
            std::cout << "Output: " << Sout.to_string() << std::endl;
            std::cerr << "Isomer " << i << " failed" << std::endl;
        }
    }


}