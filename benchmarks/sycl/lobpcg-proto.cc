#include <fullerenes/argparser.hh>
#include <fullerenes/sycl-headers/fill.hh>
#include <fullerenes/kernel-headers/all-kernels.hh>

template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG_V1(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest);

template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG_V2(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest);

template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG_V3(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest, Span<T> eigvects, Span<T> eigvals);

template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG_V4(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest, Span<T> eigvects, Span<T> eigvals);
int main(int argc, char** argv){

    CmdArgs args;
    parseArguments(argc, argv, args);
    size_t N = args.natoms;
    size_t BatchSize = args.nisomers;
    std::string device_type = args.device_type;
    size_t maxiter = args.nlanczos;
    std::cout << "N = " << N << " BatchSize = " << BatchSize << " Device Type = " << device_type << " Max Iterations = " << maxiter << std::endl;
    auto queue = SyclQueue(device_type);

    FullereneBatch<float, uint16_t> batch(N, BatchSize);
    DualizeFunctor<float, uint16_t> dualize;
    TutteFunctor<float, uint16_t> tutte_layout;
    SphericalProjectionFunctor<float, uint16_t> spherical_projection;
    ForcefieldOptimizeFunctor<PEDERSEN, float, uint16_t> forcefield_optimize;
    ForcefieldOptimizeFunctor<PEDERSEN, double, uint16_t> forcefield_optimize_double;

    HessianFunctor<PEDERSEN, float, uint16_t> compute_hessians;
    HessianFunctor<PEDERSEN, double, uint16_t> compute_hessians_double;

    EigenFunctor<EigensolveMode::FULL_SPECTRUM, float, uint16_t> eigensolve;
    EigenFunctor<EigensolveMode::FULL_SPECTRUM, double, uint16_t> eigensolve_double;

    fill(batch);
    dualize(queue, batch, LaunchPolicy::SYNC);
    tutte_layout(queue, batch, LaunchPolicy::SYNC);
    //spherical_projection(queue, batch, LaunchPolicy::SYNC);
    for(int i = 0; i < batch.size(); i++){
        spherical_projection(queue, batch[i], LaunchPolicy::SYNC);
    }
    //The batched version of spherical projection is not deterministic (Floating point associativity of atomic operations)
    FullereneBatch<double, uint16_t> batch_double(N, BatchSize);
    /* {
        auto batch_acc_X = batch.d_.X_cubic_;
        auto batch_double_acc_X = batch_double.d_.X_cubic_;
        auto batch_acc_cubic_neighbours = batch.d_.A_cubic_;
        auto batch_double_acc_cubic_neighbours = batch_double.d_.A_cubic_;

        for(int i = 0; i < N*BatchSize; i++){
            for(int j = 0; j < 3; j++){
                batch_double_acc_X[i][j] = batch_acc_X[i][j];
                batch_double_acc_cubic_neighbours[i][j] = batch_acc_cubic_neighbours[i][j];
            }
        }
    } */
    for (int i = 0; i < batch.size(); i++){
        batch_double.push_back(batch[i]);
    }
    forcefield_optimize_double(queue, batch_double, LaunchPolicy::SYNC, 5*N, 5*N);
    forcefield_optimize(queue, batch, LaunchPolicy::SYNC, 5*N, 5*N);
    

    SyclVector<float> hessians((N*90*BatchSize));
    SyclVector<double> hessians_double((N*90*BatchSize));
    SyclVector<uint16_t> cols((N*90*BatchSize));
    SyclVector<float> eigenvalues((21*BatchSize));
    SyclVector<double> eigenvalues_double((21*BatchSize));
    SyclVector<float> eigenvects((N*3*21*BatchSize));
    SyclVector<double> eigenvects_double;

    compute_hessians_double(queue, batch_double, LaunchPolicy::SYNC, hessians_double, cols);
    compute_hessians(queue, batch, LaunchPolicy::SYNC, hessians, cols);

    //eigensolve(queue, batch, LaunchPolicy::SYNC, hessians, cols, N*3 - 6, eigenvalues, eigenvects);
    //eigensolve_double(queue, batch_double, LaunchPolicy::SYNC, hessians_double, cols, N*3 - 6, eigenvalues_double, eigenvects_double);



    int m = 3;
    int n = 3;
    //std::vector<float> A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    //std::vector<int> cols = {0, 1, 2, 1, 2, 0, 2, 0, 1};


    /* 
    std::cout << "Starting LOBPCG-V1" << std::endl;
    auto T0 = std::chrono::high_resolution_clock::now();
    LOBPCG_V1<float, uint16_t, 21, 30>(queue, hessians, cols, BatchSize, (int)N*3, maxiter,true);
    auto T1 = std::chrono::high_resolution_clock::now();

 */
    //std::cout << "Starting LOBPCG-V2" << std::endl;
    //auto T2 = std::chrono::high_resolution_clock::now();
    //LOBPCG_V2<float, uint16_t, 21, 30>(queue, hessians, cols, BatchSize, (int)N*3, maxiter, true);
    //auto T3 = std::chrono::high_resolution_clock::now();
    std::cout << "Starting LOBPCG-V3" << std::endl;
     auto T4 = std::chrono::high_resolution_clock::now();
    //LOBPCG_V3<float, uint16_t, 21, 30>(queue, hessians, cols, BatchSize, (int)N*3, maxiter, true, eigenvects, eigenvalues);
    auto T5 = std::chrono::high_resolution_clock::now();

    std::cout << "Starting LOBPCG-V4" << std::endl;
    //auto T6 = std::chrono::high_resolution_clock::now();
    LOBPCG_V4<float, uint16_t, 3, 30>( queue, 
                                        hessians, 
                                        cols, 
                                        BatchSize, 
                                        (int)N*3, 
                                        maxiter, 
                                        true, 
                                        eigenvects, 
                                        eigenvalues);
    //auto T7 = std::chrono::high_resolution_clock::now();

    std::cout << "Eigenvalues: " << eigenvalues << std::endl;

    //std::cout << "LOBPCG-V1: " << float(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count())/float(BatchSize) << " ms/isomer" << std::endl;
    //std::cout << "LOBPCG-V2: " << float(std::chrono::duration_cast<std::chrono::milliseconds>(T3 - T2).count())/float(BatchSize) << " ms/isomer" << std::endl;
    //std::cout << "LOBPCG-V3: " << float(std::chrono::duration_cast<std::chrono::milliseconds>(T5 - T4).count())/float(BatchSize) << " ms/isomer" << std::endl;
    //std::cout << "LOBPCG: " << std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count()/10 << " ms" << std::endl;
    //LOBPCG<double, uint16_t, 21, 30>(queue, hessians_double, cols, BatchSize, N*3, maxiter);

    std::vector <float> matrices(N*3*N*3*BatchSize);
    std::vector <double> matrices_double(N*3*N*3*BatchSize);
        std::vector <float> vect_eigenvalues(21*BatchSize);
        std::vector <double> vect_eigenvalues_double(21*BatchSize);
        std::vector <std::array<float,3>> vect_X(N*BatchSize);
        std::vector <std::array<double,3>> vect_X_double(N*BatchSize);
        auto acc_eigenvalues    = eigenvalues;
        auto acc_hessians       = hessians;
        auto acc_cols           = cols;
        auto acc_X              = batch.d_.X_cubic_;

        auto acc_eigenvalues_double    = eigenvalues_double;
        auto acc_hessians_double       = hessians_double;
        auto acc_X_double              = batch_double.d_.X_cubic_;

        for (size_t ii = 0; ii < BatchSize; ii++)
        {
            //Create the matrices (Densely stored) from the hessians and cols
            for (int i = 0; i < N*3; i++){
                for (int j = 0; j < 30; j++){
                    matrices[ii*N*3*N*3 + i*N*3 + acc_cols[ii*90*N + i*30 + j]] = acc_hessians[ii*90*N + i*30 + j];
                    matrices_double[ii*N*3*N*3 + i*N*3 + acc_cols[ii*90*N + i*30 + j]] = acc_hessians_double[ii*90*N + i*30 + j];
                }
            }

            //Store the eigenvalues
            for (int i = 0; i < 21; i++){
                vect_eigenvalues[ii*21 + i] = acc_eigenvalues[ii*21 + i];
                vect_eigenvalues_double[ii*21 + i] = acc_eigenvalues_double[ii*21 + i];
            }
            for (int i = 0; i < N; i++){
                vect_X[ii*N + i] = acc_X[ii*N + i];
                vect_X_double[ii*N + i] = acc_X_double[ii*N + i];
            }
        }

        std::ofstream out_matrix("matrices.float32", ios::out | ios::binary);
        std::ofstream out_eigenvalues("eigenvalues.float32", ios::out | ios::binary);
        std::ofstream out_X("X.float32", ios::out | ios::binary);
        out_matrix.write((char*)&matrices[0], matrices.size()*sizeof(float));
        out_eigenvalues.write((char*)&vect_eigenvalues[0], vect_eigenvalues.size()*sizeof(float));
        out_X.write((char*)&vect_X[0], vect_X.size()*sizeof(std::array<float,3>));

        std::ofstream out_matrix_double("matrices.float64", ios::out | ios::binary);
        std::ofstream out_eigenvalues_double("eigenvalues.float64", ios::out | ios::binary);
        std::ofstream out_X_double("X.float64", ios::out | ios::binary);
        out_matrix_double.write((char*)&matrices_double[0], matrices_double.size()*sizeof(double));
        out_eigenvalues_double.write((char*)&vect_eigenvalues_double[0], vect_eigenvalues_double.size()*sizeof(double));
        out_X_double.write((char*)&vect_X_double[0], vect_X_double.size()*sizeof(std::array<double,3>));
    return 0;
}

