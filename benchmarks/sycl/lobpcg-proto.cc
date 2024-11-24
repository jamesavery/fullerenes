#include <fullerenes/argparser.hh>
#include <fullerenes/sycl-headers/fill.hh>
#include <fullerenes/kernel-headers/all-kernels.hh>

template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG(SyclQueue &ctx, Span<T> A, Span<K> cols, int m, size_t maxiters);


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
    fill(batch);
    dualize(queue, batch, LaunchPolicy::SYNC);
    tutte_layout(queue, batch, LaunchPolicy::SYNC);
    spherical_projection(queue, batch, LaunchPolicy::SYNC);
    FullereneBatch<double, uint16_t> batch_double(N, BatchSize);
    {
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
    }
    forcefield_optimize_double(queue, batch_double, LaunchPolicy::SYNC, 5*N, 5*N);
    forcefield_optimize(queue, batch, LaunchPolicy::SYNC, 5*N, 5*N);
    

    SyclVector<float> hessians((N*90*BatchSize));
    SyclVector<double> hessians_double((N*90*BatchSize));
    SyclVector<uint16_t> cols((N*90*BatchSize));

    compute_hessians_double(queue, batch_double, LaunchPolicy::SYNC, hessians_double, cols);
    compute_hessians(queue, batch, LaunchPolicy::SYNC, hessians, cols);



    int m = 3;
    int n = 3;
    //std::vector<float> A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    //std::vector<int> cols = {0, 1, 2, 1, 2, 0, 2, 0, 1};

    
    LOBPCG<float, uint16_t, 3, 30>(queue, hessians, cols, (int)N*3, maxiter);
    //LOBPCG<double, uint16_t, 3, 30>(queue, hessians_double, cols, N*3, maxiter);
    return 0;
}

