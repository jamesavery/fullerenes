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
    size_t Nisomers = IsomerDB::number_isomers(N);
    size_t BatchSize = (args.nisomers == 0 || args.nisomers > Nisomers)? Nisomers : args.nisomers;
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
    sycl::buffer<node_t,1> G_back(sycl::range<1>(Nf*6));
    sycl::buffer<node_t,1> Deg_back(sycl::range<1>(Nf*1));

    neighbours_t input_neighbours(Nf);
    neighbours_t output_neighbours(N);
    neighbours_t back_neighbours(Nf);      

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

        dualize_general<3,6>(Q, G_out, Deg_out, G_back, Deg_back, N, Nf, LaunchPolicy::SYNC);
        {
            host_accessor G_back_acc(G_back, read_only);
            host_accessor Deg_back_acc(Deg_back, read_write);            
            for (size_t j = 0; j < Nf; j++){
                if (Deg_back_acc[j] == 6){
                    back_neighbours[j] = {G_back_acc[j*6], G_back_acc[j*6 + 1], G_back_acc[j*6 + 2], G_back_acc[j*6 + 3], G_back_acc[j*6 + 4], G_back_acc[j*6 + 5]};
                } else {
                    back_neighbours[j] = {G_back_acc[j*6], G_back_acc[j*6 + 1], G_back_acc[j*6 + 2], G_back_acc[j*6 + 3], G_back_acc[j*6 + 4]};
                }
            }
        }
        Triangulation pg(input_neighbours, true);
        PlanarGraph   pg_dref(pg.dual_graph());
        PlanarGraph   pg_dual(Graph(output_neighbours, true));
        Triangulation pg_back(back_neighbours, true);        

        spiral_nomenclature Sin(pg,       spiral_nomenclature::FULLERENE, spiral_nomenclature::TRIANGULATION, true);
        spiral_nomenclature Sout(pg_dual, spiral_nomenclature::FULLERENE, spiral_nomenclature::CUBIC, true);
        spiral_nomenclature Sref(pg_dref, spiral_nomenclature::FULLERENE, spiral_nomenclature::CUBIC, true);        
        spiral_nomenclature Sback(pg_back,spiral_nomenclature::FULLERENE, spiral_nomenclature::TRIANGULATION, true);
        
        bool validates_all_checks = (Sout.spiral == Sref.spiral) 
                                 && (Sin.spiral  == Sback.spiral);

        //std::cout << "Isomer " << i << " Is isomorphic: " << (validates_all_checks? "True":"False") << std::endl;
        if (!validates_all_checks) {
                     std::cout << "Input:  " << Sin.spiral  << std::endl;
            std::cout << "Output: " << Sout.spiral << std::endl;
            std::cout << "Ref:    " << Sref.spiral << std::endl;            
            std::cerr << "Isomer " << i << " failed" << std::endl;
        }
        if(i % 1000 == 999) {
            std::cout << (i+1) << " isomers passed isomorphism check.\n";
        }    
    }
    std::cout << "All " << BatchSize << " isomers passed isomorphism check.\n";

}