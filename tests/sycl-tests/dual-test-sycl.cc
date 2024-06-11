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
    
    /* if (use_double_precision) {
        std::cout << "Using double precision" << std::endl;
        using real_t = double;
    } else {
        std::cout << "Using single precision" << std::endl;
        using real_t = float;
    } */


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
    
    IsomerBatch<real_t,node_t> batch(N, BatchSize);
    fill(batch);
/*     dualize(Q, batch, LaunchPolicy::SYNC);
    tutte_layout(Q, batch, LaunchPolicy::SYNC);
    {
        host_accessor positions_acc(batch.xys, read_only);
        std::vector<std::array<real_t, 2>> positions(N);
        for (int i = 0; i < N; i++) {
            positions[i] = positions_acc[i];
        }
        ofstream myfile("tutte_embedding_N=" + std::to_string(N) + "_dims_" + std::to_string(N) + "_X_2" ".xyz");
        myfile.write((char*)&positions , N * sizeof(std::array<real_t, 2>));
    }
 */
    

    sycl::buffer<uint16_t, 1> cubic_degrees(N); //Will be filled with 3s if working correctly
    {
        auto dual_neighbours_acc = batch.dual_neighbours.get_access<sycl::access::mode::read>();
        auto face_degrees_acc = batch.face_degrees.get_access<sycl::access::mode::read>();
        std::cout << "Dual neighbours [From Buckygen]" << std::endl;
        for (int i = 0; i < Nf; i++) {
            for (int j = 0; j < 6; j++) {
                std::cout << dual_neighbours_acc[i*6 + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    dualize_general<6, 3>(Q, batch.dual_neighbours, batch.face_degrees, batch.cubic_neighbours, cubic_degrees, Nf, N, LaunchPolicy::SYNC);  
    {
        auto cubic_neighbours_acc = batch.cubic_neighbours.get_access<sycl::access::mode::read>();
        auto cubic_degrees_acc = cubic_degrees.get_access<sycl::access::mode::read>();
        std::cout << "Cubic neighbours [G*]" << std::endl;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << cubic_neighbours_acc[i*3 + j] << " ";
            }
            std::cout << std::endl;
        }
    }
    dualize_general<3, 6>(Q, batch.cubic_neighbours, cubic_degrees, batch.dual_neighbours, batch.face_degrees, N, Nf, LaunchPolicy::SYNC);
    {
        auto dual_neighbours_acc = batch.dual_neighbours.get_access<sycl::access::mode::read>();
        auto face_degrees_acc = batch.face_degrees.get_access<sycl::access::mode::read>();
        std::cout << "Dual neighbours [G**]" << std::endl;
        for (int i = 0; i < Nf; i++) {
            for (int j = 0; j < 6; j++) {
                std::cout << dual_neighbours_acc[i*6 + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    neighbours_t c_neighbours(N);
    neighbours_t d_neighbours(Nf);
    {
        auto cubic_neighbours_acc = batch.cubic_neighbours.get_access<sycl::access::mode::read>();
        auto dual_neighbours_acc = batch.dual_neighbours.get_access<sycl::access::mode::read>();
        auto face_degrees_acc = batch.face_degrees.get_access<sycl::access::mode::read>();
        for (int i = 0; i < N; i++) {
            c_neighbours[i] = {cubic_neighbours_acc[i*3], cubic_neighbours_acc[i*3 + 1], cubic_neighbours_acc[i*3 + 2]};
        }
        for (int i = 0; i < Nf; i++) {
            if (face_degrees_acc[i] == 6) {
                d_neighbours[i] = {dual_neighbours_acc[i*6], dual_neighbours_acc[i*6 + 1], dual_neighbours_acc[i*6 + 2], dual_neighbours_acc[i*6 + 3], dual_neighbours_acc[i*6 + 4], dual_neighbours_acc[i*6 + 5]};
            } else if (face_degrees_acc[i] == 5) {
                d_neighbours[i] = {dual_neighbours_acc[i*6], dual_neighbours_acc[i*6 + 1], dual_neighbours_acc[i*6 + 2], dual_neighbours_acc[i*6 + 3], dual_neighbours_acc[i*6 + 4]};
            }
        }
    }
    PlanarGraph G(c_neighbours);
    PlanarGraph Gstar(d_neighbours);

    spiral_nomenclature spiralG(G,spiral_nomenclature::naming_scheme_t::CAGE, spiral_nomenclature::construction_scheme_t::CUBIC);
    spiral_nomenclature spiralGstar(Gstar,spiral_nomenclature::naming_scheme_t::CAGE, spiral_nomenclature::construction_scheme_t::TRIANGULATION);

    std::cout << "Spiral G: " << spiralG << std::endl;
    std::cout << "Spiral G*: " << spiralGstar << std::endl;

    dualize_general<6, 3>(Q, batch.dual_neighbours, batch.face_degrees, batch.cubic_neighbours, cubic_degrees, Nf, N, LaunchPolicy::SYNC);  
    {
        auto cubic_neighbours_acc = batch.cubic_neighbours.get_access<sycl::access::mode::read>();
        std::cout << "Cubic neighbours [G***]" << std::endl;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << cubic_neighbours_acc[i*3 + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    dualize_general<3, 6>(Q, batch.cubic_neighbours, cubic_degrees, batch.dual_neighbours, batch.face_degrees, N, Nf, LaunchPolicy::SYNC);
    {
        auto dual_neighbours_acc = batch.dual_neighbours.get_access<sycl::access::mode::read>();
        auto face_degrees_acc = batch.face_degrees.get_access<sycl::access::mode::read>();
        std::cout << "Dual neighbours [G***]" << std::endl;
        for (int i = 0; i < Nf; i++) {
            for (int j = 0; j < 6; j++) {
                std::cout << dual_neighbours_acc[i*6 + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }


}