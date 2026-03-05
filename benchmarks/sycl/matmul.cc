
//Test the performance of various ways to multiply a sparse Matrix with a Dense Matrix

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include <numeric>
#include <sycl/sycl.hpp>


#define BDIM 256
#define PRINTF(fmt, ...) sycl::ext::oneapi::experimental::printf(fmt, __VA_ARGS__)
#define PRINTS(str) sycl::ext::oneapi::experimental::printf(str)

/* void load_A(int tid, int i, int j, int M, int K, int TILE_M, int TILE_K, const sycl::half* Aacc, sycl::local_accessor<sycl::half, 1> localA){
    int offset = i*K + j;
    auto limit = sycl::min(TILE_M*TILE_K, sycl::max(M*K - offset, 0));
    auto divisor = sycl::min(K, TILE_K);
    for(int i = tid; i < limit; i+= BDIM){
        auto row = i/divisor;
        auto col = i%divisor;
        localA[row*TILE_K + col] = Aacc[offset + row*K + col];
    }
}

void load_B(int tid, int i, int j, int N, int K, int TILE_N, int TILE_K, const sycl::half* Bacc, sycl::local_accessor<sycl::half, 1> localB){
    int offset = i*N + j;
    auto limit = sycl::min(TILE_N*TILE_K, sycl::max(N*K - offset, 0));
    auto divisor = sycl::min(N, TILE_N);
    for(int i = tid; i < limit; i+= BDIM){
        auto row = i/divisor;
        auto col = i%divisor;
        localB[row*TILE_N + col] = Bacc[offset + row*N + col];
    }
} */

int main(int argc, char** argv){
    using namespace sycl;
    using use = ext::intel::experimental::matrix::use;
    using layout = ext::oneapi::experimental::matrix::layout;
    using bfloat16 = ext::oneapi::bfloat16;
    using float16 = sycl::float16;
    using namespace sycl::ext::oneapi::experimental::matrix;


    
    queue Q = queue(gpu_selector_v, property::queue::in_order{});
    int M = argc > 1 ? std::stoi(argv[1]) : 16;
    int N = argc > 2 ? std::stoi(argv[2]) : 16;
    int K = argc > 3 ? std::stoi(argv[3]) : 16;
    int TILE_M = argc > 4 ? std::stoi(argv[4]) : 32;
    int TILE_N = argc > 5 ? std::stoi(argv[5]) : 32;
    int TILE_K = argc > 6 ? std::stoi(argv[6]) : 32;
    int Nmatrices = argc > 7 ? std::stoi(argv[7]) : 1;
    float sparsity = argc > 8 ? std::stof(argv[8]) : 0.5;
    int Nruns = argc > 9 ? std::stoi(argv[9]) : 1;
    int SZ = Q.get_device().get_info<info::device::sub_group_sizes>()[0];


    int numcols = int((1.-sparsity)*N);
    std::vector<half> A(M*K*Nmatrices);
    std::vector<int> Cols(M*numcols*Nmatrices);
    std::vector<half> B(N*K*Nmatrices);
    std::vector<float> C(M*N*Nmatrices);

    std::default_random_engine generator(0);
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    std::vector<int> linspace(K);
    std::iota(linspace.begin(), linspace.end(), 0);
    //
    for (int i = 0; i < M; i++){
        //Same sparsity in each row, but different pattern. This is intentional to mimic real use case.
        std::shuffle(linspace.begin(), linspace.end(), generator);
        std::sort(linspace.begin(), linspace.begin() + numcols);
        std::copy(linspace.begin(), linspace.begin() + numcols, Cols.begin() + i*numcols);
    }
    

    //Fill randum values in A
    std::transform(A.begin(), A.end(), A.begin(), [&](float x){return distribution(generator);});

    //Fill random values in B
    std::transform(B.begin(), B.end(), B.begin(), [&](float x){return distribution(generator);});

    sycl::buffer<half, 1> A_buf(A.data(), range<1>(M*K*Nmatrices));
    sycl::buffer<int, 1> Cols_buf(Cols.data(), range<1>(M*numcols*Nmatrices));
    sycl::buffer<half, 1> B_buf(B.data(), range<1>(N*K*Nmatrices));
    sycl::buffer<float, 1> C_buf(C.data(), range<1>(M*N*Nmatrices));

    /* Q.submit([&](sycl::handler& h){
        auto Aacc = A_buf.get_access<sycl::access::mode::read>(h);
        auto Bacc = B_buf.get_access<sycl::access::mode::read>(h);
        auto Cacc = C_buf.get_access<sycl::access::mode::write>(h);
        auto Colsacc = Cols_buf.get_access<sycl::access::mode::read>(h);
        h.parallel_for<class matmul>(range<1>(N*N*Nmatrices), [=](sycl::id<1> item){
            int i = item[0]/N;
            int j = item[0]%N;
            float sum = 0;
            for (int k = 0; k < numcols; k++){
                sum += Aacc[i*numcols + k]*Bacc[Colsacc[i*numcols + k]*N + j];
            }
            Cacc[item] = sum;
        });
    });

    Q.wait(); */

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    Q.submit([&](sycl::handler& h){
        auto Aacc = A_buf.get_access<sycl::access::mode::read>(h);
        auto Bacc = B_buf.get_access<sycl::access::mode::read>(h);
        auto Cacc = C_buf.get_access<sycl::access::mode::write>(h);
        auto Colsacc = Cols_buf.get_access<sycl::access::mode::read>(h);
        auto localA = local_accessor<half, 1>(TILE_M*TILE_K, h);
        auto localB = local_accessor<half, 1>(TILE_N*TILE_K, h);
        auto localC = local_accessor<float, 1>(TILE_M*TILE_N, h);

        h.parallel_for<class matmul2>(nd_range<1>(BDIM*Nmatrices, BDIM), [=](sycl::nd_item<1> item){
            auto cta = item.get_group();
            int bid = item.get_group_linear_id();
            int tid = cta.get_local_linear_id();
            int Aoffset = bid*numcols*M;
            int Boffset = bid*K*N;
            auto col_ptr = Colsacc.get_pointer() + Aoffset;
            //auto A_ptr = Aacc.get_pointer() + Aoffset;
            auto B_ptr = Bacc.get_pointer() + Boffset;
            auto A_ptr = Aacc.get_pointer() + bid*M*K;
            auto C_ptr = Cacc.get_pointer() + bid*M*N;
            auto Nwarps = cta.get_local_linear_range()/SZ;
            auto subgroup_id = tid/SZ;

            sycl::sub_group sg = item.get_sub_group();
            joint_matrix<sycl::sub_group, half, use::a, 16, 16, layout::row_major> sub_a;
            joint_matrix<sycl::sub_group, half, use::b, 16, 16, layout::row_major> sub_b;
            joint_matrix<sycl::sub_group, float, use::accumulator, 16, 16> sub_c;

            /* auto load_Asparse = [&](int i, int j){
                memset(localA.get_pointer(), 0, TILE_M*TILE_K*sizeof(half));
                auto offset = i*numcols;
                auto limit = sycl::min(numcols*TILE_M, sycl::max(M*numcols - offset, 0));
                for(int i = tid; i < limit; i+= BDIM){
                    auto row = i/numcols;
                    auto col = i%numcols;
                    auto idx = offset + row*numcols + col;
                    auto c = col_ptr[idx] - j;
                    if (c >= 0 && c < TILE_K)
                        localA[row*TILE_K + c] = A_ptr[idx];
                }
            }; */

            auto load_A = [&](int i, int j){
                int offset = i*K + j;
                auto limit = sycl::min(TILE_M*TILE_K, sycl::max(M*K - offset, 0));
                auto divisor = sycl::min(K, TILE_K);
                for(int i = tid; i < limit; i+= BDIM){
                    auto row = i/divisor;
                    auto col = i%divisor;
                    localA[row*TILE_K + col] = A_ptr[offset + row*K + col];
                }
            };

            auto load_B = [&](int i, int j){
                int offset = i*N + j;
                auto limit = sycl::min(TILE_N*TILE_K, sycl::max(N*K - offset, 0));
                auto divisor = sycl::min(N, TILE_N);
                for(int i = tid; i < limit; i+= BDIM){
                    auto row = i/divisor;
                    auto col = i%divisor;
                    localB[row*TILE_N + col] = B_ptr[offset + row*N + col];
                }
            };
            auto store_C = [&](int i, int j){
                int offset = i*N + j;
                auto limit = sycl::min(TILE_N*TILE_M, sycl::max(N*M - offset, 0));
                auto divisor = sycl::min(N, TILE_N);
                for(int i = tid; i < limit; i+= BDIM){
                    auto row = i/divisor;
                    auto col = i%divisor;
                    C_ptr[offset + row*N + col] = localC[row*TILE_N + col];
                }
            };
            int NsubM = TILE_M/16;
            int NsubN = TILE_N/16;
            int NsubK = TILE_K/16;
    for(int iii = 0; iii < Nruns; iii++){
            //B is more memory intensive (it is not sparse) , so we wish to reuse it as much as possible
            for (int i = 0; i < M; i+=TILE_M){
                for (int j = 0; j < N; j+=TILE_N){
                    //We need to reset the local C block
                    sycl::group_barrier(cta);
                    if(tid==0) memset(localC.get_pointer(), 0, TILE_M*TILE_N*sizeof(float));

                    for(int k = 0; k  < K; k+=TILE_K){
                        load_A(i, k);
                        load_B(k, j);
                        sycl::group_barrier(cta);
                        for (int kk = subgroup_id; kk < NsubM*NsubN; kk+=Nwarps){
                            auto cr = kk/NsubN; //Row in the local C block
                            auto cc = kk%NsubN; //Column in the local C block
                            //Iterate over the rows of A and the columns of B
                            joint_matrix_load(sg, sub_c, localC.get_pointer() + cr*TILE_N*16 + cc*16, TILE_N, layout::row_major);
                            for (int ii = 0; ii < NsubK; ii++)
                            {   
                                joint_matrix_load(sg, sub_a, localA.get_pointer() + cr*TILE_K*16 + ii*16, TILE_K);
                                joint_matrix_load(sg, sub_b, localB.get_pointer() + cc*16 + ii*TILE_N*16, TILE_N);
                                sub_c = joint_matrix_mad(sg, sub_a, sub_b, sub_c);
                            }
                            joint_matrix_store(sg, sub_c, localC.get_pointer() + cr*TILE_N*16 + cc*16, TILE_N, layout::row_major);
                        }
                        sycl::group_barrier(cta);
                    }
                    store_C(i, j);
                    //cta.async_work_group_copy(C_ptr + i*N + j, localC.get_pointer(), TILE_M*TILE_N);
                }
            }
            //Cacc[item] = sum;
    }
            //Cacc[Boffset + item.get_local_linear_id()] = 2;
            
        });
    });
    Q.wait();

    auto elapsed_nano = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start).count();
    {
        auto Cacc = C_buf.get_access<sycl::access::mode::read>();
        auto Aacc = A_buf.get_access<sycl::access::mode::read>();
        auto Bacc = B_buf.get_access<sycl::access::mode::read>();

        /* std::cout << "C: =======================" << std::endl;
        for (int i = 0; i < M; i++){
            for(int j = 0; j < N; j++){
                printf("%f ", Cacc[i*N + j]);
            }
            std::cout << std::endl;
        } */
        /* std::cout << "\nA: =======================" << std::endl;
        for (int i = 0; i < M; i++){
            int col = 0;
            for(int j = 0; j < K; j++){
                if (Cols[i*numcols + col] == j){
                    printf("%f ", Aacc[i*numcols + col]);
                    col++;
                }
                else{
                    printf("%f ", 0.0);
                }
            }
            std::cout << std::endl;
        } */
        

        /* std::cout << "\nB: =======================" << std::endl;
        for (int i = 0; i < K; i++){
            for(int j = 0; j < N; j++){
                printf("%f ", Bacc[i*N + j]);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl; */
    }

    //Validate
    {
        auto Cacc = C_buf.get_access<sycl::access::mode::read>();
        auto Aacc = A_buf.get_access<sycl::access::mode::read>();
        auto Bacc = B_buf.get_access<sycl::access::mode::read>();
        float abs_err = 0;
        float rel_err = 0;
        std::cout << "C Err Matrix: =======================" << std::endl;
        for (int i = 0; i < M; i++){
            for(int j = 0; j < N; j++){
                float sum = 0;
                for (int k = 0; k < K; k++){
                    //sum += Aacc[i*numcols + k]*Bacc[Cols[i*numcols + k]*N + j];
                    sum += Aacc[i*K + k]*Bacc[k*N + j];
                }
                printf("%f ", std::abs(Cacc[i*N + j] - sum)/std::abs(sum));
                abs_err += std::abs(Cacc[i*N + j] - sum);
                rel_err += std::abs(Cacc[i*N + j] - sum)/std::abs(sum);
            }
            std::cout << std::endl;
        }
        std::cout << "Average Relative Error: " << (rel_err/(N*M)) * 100. << " %" << std::endl;
        std::cout << "Average Absolute Error: " << abs_err/(N*M) << std::endl;
    }

    std::cout << "Time: " << elapsed_nano << " ms" << std::endl;
    auto denseGFLOPS = Nruns*2.*M*N*K*Nmatrices/elapsed_nano/1e6;
    auto sparseGFLOPS = Nruns*2.*numcols*N*K*Nmatrices/elapsed_nano/1e6;
    std::cout << "Dense Performance  : " << denseGFLOPS << " GFLOPS\n";
    std::cout << "Sparse Performance : " << sparseGFLOPS << " GFLOPS\n";
    std::cout << "Bandwidth: " << Nruns*1.*(TILE_K*TILE_N + TILE_K*TILE_M + TILE_N*TILE_M)*Nmatrices* (M/TILE_M)*(N/TILE_N)*(K/TILE_K)/elapsed_nano/1e6 << " GB/s" << std::endl;
    

}