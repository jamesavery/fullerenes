#include <fullerenes/graph.hh>
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
    
    dualize(Q, batch, LaunchPolicy::SYNC);
    tutte_layout(Q, batch, LaunchPolicy::SYNC);
    spherical_projection(Q, batch, LaunchPolicy::SYNC);

    IsomerBatch<double,node_t> batch_double(N, BatchSize);
    {
        auto acc_double = sycl::host_accessor(batch_double.X);
        auto acc = sycl::host_accessor(batch.X);
        auto acc_double_cubic_neighbours = sycl::host_accessor(batch_double.cubic_neighbours);
        auto acc_cubic_neighbours = sycl::host_accessor(batch.cubic_neighbours);

        for (size_t i = 0; i < N*BatchSize; i++)
        {
            for (size_t j = 0; j < 3; j++)
            {
                acc_double_cubic_neighbours[i*3 + j] = acc_cubic_neighbours[i*3 + j];
                acc_double[i][j] = acc[i][j];
            }
        }
    }

    forcefield_optimize(Q, batch_double, 5*N, 5*N, LaunchPolicy::SYNC);

    sycl::buffer<double, 1> hessians(range<1>(N*90*BatchSize));
    sycl::buffer<node_t, 1> cols(range<1>(N*90*BatchSize));
    sycl::buffer<double, 1> eigenvalues(range<1>(N*3*BatchSize));
    sycl::buffer<double, 1> eigenvalue_ends(range<1>(BatchSize*2));
    sycl::buffer<double, 1> eigenvectors(range<1>(2*N*3*N*3*BatchSize));
    compute_hessians(Q, batch_double, hessians, cols, LaunchPolicy::SYNC);
    //eigensolve<EigensolveMode::FULL_SPECTRUM>(Q, batch, hessians, cols, eigenvalues);
    eigensolve<EigensolveMode::FULL_SPECTRUM_VECTORS>(Q, batch_double, hessians, cols, eigenvalues, LaunchPolicy::SYNC, Nlanczos, eigenvectors);
    //eigensolve<EigensolveMode::ENDS>(Q, batch, hessians, cols, eigenvalue_ends, LaunchPolicy::SYNC, Nlanczos);

    //Print First Eigenvector
    {
        auto acc_eigenvectors = sycl::host_accessor(eigenvectors);
        auto acc_eigenvalues = sycl::host_accessor(eigenvalues);
        auto acc_X = sycl::host_accessor(batch.X);
        {
            std::cout << "Eigenvector " << print_vector << " : " << std::endl;
            for (size_t i = 0; i < N; i++)
            {
                std::cout << acc_eigenvectors[print_vector*N*3 + i*3] << " " << acc_eigenvectors[print_vector*N*3 + i*3 + 1] << " " << acc_eigenvectors[print_vector*N*3 + i*3 + 2] << std::endl;
            }
            std::cout << "Eigenvalue "<< print_vector << " : " << acc_eigenvalues[print_vector] << std::endl;
        }
    }

    
    {   
        std::vector <double> matrices(N*3*N*3*BatchSize);
        std::vector <double> vect_eigenvalues(N*3*BatchSize);
        std::vector <std::array<double,3>> vect_X(N*BatchSize);
        auto acc_eigenvalues    = sycl::host_accessor(eigenvalues);
        auto acc_hessians       = sycl::host_accessor(hessians);
        auto acc_cols           = sycl::host_accessor(cols);
        auto acc_X              = sycl::host_accessor(batch_double.X);

        for (size_t ii = 0; ii < BatchSize; ii++)
        {
            //Create the matrices (Densely stored) from the hessians and cols
            for (int i = 0; i < N*3; i++){
                for (int j = 0; j < 30; j++){
                    matrices[ii*N*3*N*3 + i*N*3 + acc_cols[ii*90*N + i*30 + j]] = acc_hessians[ii*90*N + i*30 + j];
                }
            }

            //Store the eigenvalues
            for (int i = 0; i < N*3; i++){
                vect_eigenvalues[ii*N*3 + i] = acc_eigenvalues[ii*N*3 + i];
            }
            for (int i = 0; i < N; i++){
                vect_X[ii*N + i] = acc_X[ii*N + i];
            }
        }

        std::ofstream out_matrix("matrices.float64", ios::out | ios::binary);
        std::ofstream out_eigenvalues("eigenvalues.float64", ios::out | ios::binary);
        std::ofstream out_X("X.float64", ios::out | ios::binary);
        out_matrix.write((char*)&matrices[0], matrices.size()*sizeof(double));
        out_eigenvalues.write((char*)&vect_eigenvalues[0], vect_eigenvalues.size()*sizeof(double));
        out_X.write((char*)&vect_X[0], vect_X.size()*sizeof(std::array<double,3>));
    }
    
    return 0;
}
