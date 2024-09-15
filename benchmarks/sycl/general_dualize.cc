#include <fullerenes/graph.hh>
#include <fullerenes/planargraph.hh>
#include <fullerenes/spiral.hh>
#include <fullerenes/sycl-headers/all-kernels.hh>
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/argparser.hh>
#include <string>
#include <algorithm>
#include <fstream>
#include <fullerenes/triangulation.hh>
#include <unistd.h>
#include <sys/wait.h>
using namespace sycl;
    // Function to parse command line arguments

int main(int argc, char** argv) {
    CmdArgs args;
    parseArguments(argc, argv, args);
    typedef float real_t;
    typedef uint16_t node_t;
    size_t Nruns = args.nruns;
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


    
/* 
    std::vector<Triangulation> Ts;
    Ts.push_back(Triangulation(spiral_nomenclature("C21888-[3300,3494,4574,4792,8540,8740,9990,10096,10295,10395,10491,10942]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C24192-[2052,3078,3087,4251,9613,10261,10900,10912,11644,11948,12044,12056]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C25460-[1,183,192,7835,7902,8267,9450,9737,9746,9815,9882,11132]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C25728-[4370,6090,6364,8015,9745,10208,10671,11000,11457,11564,12338,12746]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C28880-[80,91,102,9884,9966,10048,11108,11190,11272,11362,11442,11522]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C30020-[1,229,269,9078,10044,10120,10697,10769,10845,10931,11795,12930]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C45024-[1976,4611,5316,8071,15966,17118,18276,18679,19817,20761,20913,22162]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C47120-[150,165,180,15909,16013,16117,16543,16647,16751,19231,19326,19421]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C47616-[9633,9987,13240,13634,18080,19656,20118,20755,21212,21385,22716,23136]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C47656-[1220,7081,7370,11493,13992,14607,15827,18376,19073,20339,22930,23378]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C47656-[3,679,696,8205,8515,15716,17645,18690,21368,22155,22605,23795]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C48020-[1,7722,7778,7834,7890,7946,17543,17599,17655,17711,17767,24012]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C48840-[8965,9285,9610,12410,18096,18714,21201,21799,23834,23956,24074,24418]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C48860-[1,8007,8064,8121,8178,8235,16868,16925,16982,17039,17096,24432]-fullerene")));
    Ts.push_back(Triangulation(spiral_nomenclature("C49056=[4374,7762,8678,12888,17567,18028,18198,18827,20169,22581,22707,23966]-fullerene")));
    auto maxThreads = Q.get_device().get_info<info::device::max_compute_units>() * Q.get_device().get_info<info::device::max_work_group_size>();
    for (auto &T : Ts) {
        std::cout << "                                Nin/Nout: \t" << T.N << "/" << (T.N-2)*2 <<  "                                   " << std::endl;
        std::cout << "                                Occupancy:\t" << (float)T.N*100 / maxThreads << "% / " << (float)((T.N-2)*2)*100. / maxThreads << "%" << std::endl;
        sycl::buffer<node_t> G_in = sycl::buffer<node_t>(T.N*6);
        sycl::buffer<node_t> Deg_in = sycl::buffer<node_t>(T.N);
        sycl::buffer<node_t> G_out = sycl::buffer<node_t>(((T.N-2)*2)*3);
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
        auto Timpl0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i  < Nruns; i++) {
            dualize_general<6,3>(Q, G_in, Deg_in, G_out, Deg_out, T.N, (T.N-2)*2, LaunchPolicy::ASYNC);
        }
        Q.wait();
        auto ElapsedMicro = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - Timpl0).count() / Nruns;   



        auto T0 = std::chrono::high_resolution_clock::now();
        PlanarGraph Gref = T.dual_graph();
        auto ElapsedRefMicro = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - T0).count();

        neighbours_t neighbours((T.N-2)*2);
        {
            host_accessor G_out_acc(G_out, read_only);
            host_accessor Deg_out_acc(Deg_out, read_only);

            for (int i = 0; i < (T.N-2)*2; i++) {
                for (int j = 0; j < Deg_out_acc[i]; j++) {
                    neighbours[i].push_back(G_out_acc[i*3+j]);
                }
            }
        }   


        PlanarGraph G_star(Graph(neighbours,true));
        spiral_nomenclature refSpiral(Gref, spiral_nomenclature::FULLERENE, spiral_nomenclature::CUBIC);
        spiral_nomenclature outputSpiral(G_star, spiral_nomenclature::FULLERENE, spiral_nomenclature::CUBIC);
        auto is_isomorphic = (refSpiral.spiral.spiral_code == outputSpiral.spiral.spiral_code) ? "Yes" : "No";
        std::cout << "Are G* and reference G* isomorphic? : " << is_isomorphic << std::endl;
        std::cout << "Triangulation Specialized Reference Routine: " << ElapsedRefMicro << " µs, General Lockstep Implementation: " << ElapsedMicro << " µs" << std::endl;
        std::cout << "Speedup: " << (double)ElapsedRefMicro/ElapsedMicro << std::endl;
        

        dualize_general<3,6>(Q, G_out, Deg_out, G_in, Deg_in, (T.N-2)*2, T.N, LaunchPolicy::SYNC);
        auto Tref1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i  < Nruns; i++) {
            dualize_general<3,6>(Q, G_out, Deg_out, G_in, Deg_in, (T.N-2)*2, T.N, LaunchPolicy::ASYNC);
        }
        Q.wait();
        auto ElapsedMicroSeconds1 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - Tref1).count() / Nruns;

        auto T2 = std::chrono::high_resolution_clock::now();
        PlanarGraph Gref2 = Gref.dual_graph();
        auto ElapsedRefMicro2 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - T2).count();

        neighbours_t neighbours2(T.N);
        {
            host_accessor G_in_acc(G_in, read_only);
            host_accessor Deg_in_acc(Deg_in, read_only);

            for (int i = 0; i < T.N; i++) {
                for (int j = 0; j < Deg_in_acc[i]; j++) {
                    neighbours2[i].push_back(G_in_acc[i*6+j]);
                }
            }
        }

        PlanarGraph G_star2(Graph(neighbours2,true));
        spiral_nomenclature refSpiral2(Gref2, spiral_nomenclature::FULLERENE, spiral_nomenclature::TRIANGULATION);
        spiral_nomenclature outputSpiral2(G_star2, spiral_nomenclature::FULLERENE, spiral_nomenclature::TRIANGULATION);
        auto is_isomorphic2 = (refSpiral2.spiral.spiral_code == outputSpiral2.spiral.spiral_code) ? "Yes" : "No";
        std::cout << "Are G** and reference G** isomorphic? : " << is_isomorphic2 << std::endl;
        std::cout << "General Reference Routine: " << ElapsedRefMicro2 << " µs, General Lockstep Implementation: " << ElapsedMicroSeconds1 << " µs" << std::endl;
        std::cout << "Speedup: " << (double)ElapsedRefMicro2/ElapsedMicroSeconds1 << std::endl;
        std::cout << "=====================================================================================================" << std::endl;

    } */
    //QueueWrapper Q(device_id, device_type == "GPU" ? DeviceType::GPU : DeviceType::CPU);
    /* 
    
    BuckyGen::next_fullerene(BQ, Gr);
    std::cout << "Gr has " << Gr.neighbours.size() << " vertices" << std::endl;
 */
    //BuckyGen::buckygen_queue BQ = BuckyGen::start(20, false, false);
    pid_t pid = fork();
    if(pid == 0) {
        std::cout << "Child Process" << std::endl;
        return 0;
    } else {
        std::cout << "Parent Process" << std::endl;
        sleep(5);
        std::cout << "Done Waiting" << std::endl;
    }

    
    //DeviceWrapper D = Q.get_device();
    //Graph Gr(neighbours_t(Nf), true);
    //FullereneIsomer<float,uint16_t> isomer(Gr, false, false);
    //std::cout << "Device : " << D.get_name() <<  " Compute Units: " << D.get_property(DeviceProperty::MAX_COMPUTE_UNITS) << " Max Work Group Size: " << D.get_property(DeviceProperty::MAX_WORK_GROUP_SIZE) << std::endl;

    //std::cout <<  "Device : " << Q.get_device().get_name() << std::endl;
 //   dualize(QQ,isomer,LaunchPolicy::SYNC);
    

    waitpid(pid, NULL, 0);
}