#include <array>
#include <algorithm>
#include <fstream>
#include "sycl/sycl.hpp"
#include <oneapi/dpl/random>
#include <oneapi/dpl/iterator>
#include <fullerenes/graph.hh>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include "queue-impl.cc"
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/argparser.hh>
#include <fullerenes/sycl-headers/fill.hh>
#include <fullerenes/sycl-headers/sycl-mdspan.hh>
#include "lanczos.cc"
#define CEIL_DIV(x, y) (((x) + (y) - 1) / (y))

#if defined(SOLVER_BACKEND) && (SOLVER_BACKEND == CUSOLVER_BACKEND)
    #define USE_CUSOLVER
    #include <cusolverDn.h>
    #include <cublas_v2.h>
    #include <cusparse.h>
    #include <cuda_runtime.h>
    #include <cuda_runtime_api.h>
    #define handle_t cusolverDnHandle_t
    #define status_t cusolverStatus_t
    #define stream_t cudaStream_t
    #define params_t cusolverDnParams_t
    #define props_t cudaDeviceProp
    #define synchronize cudaDeviceSynchronize
    #define sync_stream cudaStreamSynchronize
    #define create_handle cusolverDnCreate
    #define destroy_handle cusolverDnDestroy
    #define set_stream cusolverDnSetStream
    #define create_stream cudaStreamCreate
    #define destroy_stream cudaStreamDestroy
    #define get_device_properties cudaGetDeviceProperties
    #define create_params cusolverDnCreateParams
    #define EIG_SOLVE_MODE CUSOLVER_EIG_MODE_VECTOR
    #define FILL_MODE_LOWER CUBLAS_FILL_MODE_LOWER
    #define SsyevBatched cusolverDnXsyevBatched
    #define SsyevBatched_bufferSize cusolverDnXsyevBatched_bufferSize
    #define SOLVER_COMPUTE_TYPE CUDA_R_32F
#elif defined(SOLVER_BACKEND) && (SOLVER_BACKEND == ROCSOLVER_BACKEND)
    #define USE_ROCSOLVER
    #include <rocsolver.h>
    #define handle_t rocsolver_handle
    #define status_t rocsolver_status
    #define stream_t hipStream_t
    #define synchronize hipDeviceSynchronize
    #define create_handle rocsolver_create_handle
    #define destroy_handle rocsolver_destroy_handle
    #define set_stream rocsolver_set_stream
    #define create_stream hipStreamCreate
    #define destroy_stream hipStreamDestroy
#endif

#define DECLARE_SUBSPACE\
    auto subspace = S_acc.subspan(bid * BlockVectors*9 * m, BlockVectors * m * 9);\
    auto blockX = S_acc.subspan(bid * BlockVectors*9 * m, BlockVectors*m);\
    auto blockR = S_acc.subspan(bid * BlockVectors*9 * m + BlockVectors*m, BlockVectors*m);\
    auto blockP = S_acc.subspan(bid * BlockVectors*9 * m + 2*BlockVectors*m, BlockVectors*m);\
    auto blockXtemp = S_acc.subspan(bid * BlockVectors*9 * m + 3*BlockVectors*m, BlockVectors*m);\
    auto blockRtemp = S_acc.subspan(bid * BlockVectors*9 * m + 4*BlockVectors*m, BlockVectors*m);\
    auto blockPtemp = S_acc.subspan(bid * BlockVectors*9 * m + 5*BlockVectors*m, BlockVectors*m);\
    auto blockAX = S_acc.subspan(bid * BlockVectors*9 * m + 6*BlockVectors*m, BlockVectors*m);\
    auto blockAR = S_acc.subspan(bid * BlockVectors*9 * m + 7*BlockVectors*m, BlockVectors*m);\
    auto blockAP = S_acc.subspan(bid * BlockVectors*9 * m + 8*BlockVectors*m, BlockVectors*m);

template <typename T, typename K, int BlockVectors, int NZ> class LOBPCGg {};

using namespace sycl;
//LOBPCG Algorithm implemented in SYCL
//Version 1.0
//Naive Approach: All matrices, and block vectors are stored in Global Memory


//Finds the BlockVectors- smallest or largest eigenvalues of the matrix A
template <typename T, int N, int BlockVectors>
void compute_k_eigenpairs(sycl::group<1>& cta, const Span<T>& A, const Span<T>& eigvects, const Span<T>& lambdas, T* working_space, const local_accessor<std::byte,1>& sort_scratch, const bool largest){
    auto tid = cta.get_local_linear_id();
    auto bdim = cta.get_local_linear_range();

    T* LanczosVectors = working_space;
    //The memory used for X, alphas and betas can be reused after lanczos is done, LanczosVectors must be preserved for eigenvector computation
    T* X = LanczosVectors + N*N;
    T* alphas = X + N;
    T* betas = alphas + N + 1;

    T* D = alphas;
    T* L = betas;
    T* U = L + N + 1;
    T* V = U + N*2 + 1;
    T* Q = V + N*2 + 1;
    int* k_indices = reinterpret_cast<int*>(betas);
    lanczos<T, N>(cta, A, Span<T>(X,N), Span<T>(alphas,N+1), Span<T>(betas,N+1), Span<T>(LanczosVectors,N*N));
    diagonalize<T, N>(cta, U, L, D, V, Q);
    if(tid < N) k_indices[tid] = tid;
    sycl::group_barrier(cta);

    
    const auto bytes = sycl::ext::oneapi::experimental::default_sorters::joint_sorter<>::memory_required<T>(sycl::memory_scope::work_group, N);

    sycl::ext::oneapi::experimental::joint_sort(ext::oneapi::experimental::group_with_scratchpad(cta, sycl::span{(std::byte*)sort_scratch.get_pointer(), bytes}), 
                                                k_indices, 
                                                k_indices + N,
                                                [&](auto x, auto y){ return D[x] < D[y]; });



    

    if(tid < N){
        lambdas[tid] = largest ? D[k_indices[N - 1 - tid]] : D[k_indices[tid]];
    }
    //if (tid < BlockVectors) sycl::ext::oneapi::experimental::printf("Lambdas[%d] = %f\n", tid, lambdas[tid]);

    if(tid < N){
        for(int i = 0; i < BlockVectors; i++){
            auto index = largest ? k_indices[N -1 - i] : k_indices[i];
            eigvects[i*N + tid] = 0;
            for (int j = 0; j < N; j++){
                eigvects[i*N + tid] += LanczosVectors[j*N + tid] * Q[index*N + j];
            }
        }
    }

    //Normalize eigenvectors
    for(int i = 0; i < BlockVectors; i++){
        T val = tid < N ? eigvects[i*N + tid] : 0;
        T rnorm = rsqrt(reduce_over_group(cta, val * val, sycl::plus<T>{}));

        //Invert the sign of the eigenvector if the first element is negative
        T first = tid < N ? eigvects[i*N] : 1;
        group_barrier(cta);
        if(tid<N) eigvects[i*N + tid] = sign(first) * val * rnorm;
    }

    sycl::group_barrier(cta);
} 

//Device Side Matmul
//A: m x k. Assume A is stored in thread local memory for now. I.E. Each thread has k*RowsPerThread elements of A
//blockX: m x n
//m is the number of rows in A and blockX
//k is the number of non-zero elements in each row of A, A is Square, Symmetric, and Positive Definite
//SN is the number of columns in the subspace i.e. the number of block vectors. (Subspace Dimension)
template <typename T, typename K, int RowsPerThread, int SN>
void matBlockVector(sycl::group<1>& cta ,const private_ptr<T> A, const private_ptr<K> cols, const Span<T> blockX, Span<T> AX, int n, int m){    
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
    group_barrier(cta);
    if(tid < m){
        for(int i = 0; i < RowsPerThread; i++){
            for(int ii = 0; ii < SN; ii++){
                T sum = 0;
                for(int j = 0; j < n; j++){
                    sum += A[i*n + j] * blockX[cols[i*n + j] + ii*m];
                }
                AX[ii*m + (tid + i*bdim)] = sum;
            }
        }

    }
}

template <typename T, int SN>
void compute_gram_matrix(sycl::group<1>& cta, T* St, T* S, T* StS, int m){
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
    for(int i = 0; i < SN; i++){
        for(int j = i; j < SN; j++){
            //StS[i*SN + j] = joint_reduce(cta, St + i*m, S + j*m)
        }
    }
    for (int ij = tid; ij < SN*SN; ij += bdim){
        int i = ij / SN;
        int j = ij % SN;
        StS[j*SN + i] = StS[i*SN + j];
    }
}

//Projection of A onto the krylovspace spanned by S and then 
template <typename T, typename K, int RowsPerThread, int SN>
void STAS(sycl::group<1>& cta, const private_ptr<T> A, const private_ptr<K> cols, Span<T> S, int m, int n, Span<T> StAS, local_ptr<T> Scache){
//S^T * A * S
//To produce the i,j-th element of StAS, we need to compute the dot product of the i-th row of S^T and the j-th column of A*S
//The j-th column of A*S is the matrix vector product of A and the j-th column of S
//For this reason it makes sense to store the j-th column of S in shared memory and tcd .he i-th row of S^T in shared memory (the i-th column of S)

//global_ptr<T> S_global(S);
//local_ptr<T> S_j(Scache, m);
    auto tid = cta.get_local_id(0);
    for(int j = 0; j < SN; j++){
        //Load the j-th column of S into shared memory
        for (int ix = tid; ix < m; ix += cta.get_local_range(0)) Scache[ix] = S[j*m + ix];
        sycl::group_barrier(cta);
        for(int i = j; i < SN; i++){
            //Load the i-th row of S^T into registers
            T S_itid = S[i*m + tid];
            //Compute the dot product of the i-th row of S^T and the j-th column of A*S
            T AS_tidj = 0;
            for(int k = 0; k < n; k++){
                AS_tidj += A[k] * Scache[cols[k]];
            }
            //Compute the i,j-th element of StAS, the dot product of the i-th row of S^T and the j-th column of A*S
            StAS[i*SN + j] = reduce_over_group(cta, S_itid * AS_tidj, sycl::plus<T>{});
            if(tid==0) StAS[j*SN + i] = StAS[i*SN + j];
        }
        sycl::group_barrier(cta);
    }
}

//M != N = K in place matrix matrix multiplication
template <typename T, int SN>
void inPlaceMatMatMul(sycl::group<1>& cta, Span<T> A, Span<T> B, int m){
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
    T Alocal[SN];
    for(int j = 0; j < SN; j++){
        T sum = 0;
        for(int k = 0; k < SN; k++){
            sum += A[k*m + tid] * B[j*SN + k];
        }
        Alocal[j] = sum;
    }
    sycl::group_barrier(cta);
    for(int j = 0; j < SN; j++){
        A[j*m + tid] = Alocal[j];
    }   
}

template <typename T, int SN>
void orthonormalize(sycl::group<1>& cta, Span<T> S, int m){
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
    sycl::group_barrier(cta);
    for(int i = 0; i < SN; i++){
        T* S_i = S.data() + i*m;
        T norm = sycl::sqrt(reduce_over_group(cta, S_i[tid]*S_i[tid], sycl::plus<T>{}));
        S_i[tid] /= norm;
        sycl::group_barrier(cta);
        for(int j = i+1; j < SN; j++){
            T* S_j = S.data() + j*m;
            T projection = reduce_over_group(cta, S_i[tid]*S_j[tid], sycl::plus<T>{}) * S_i[tid];
            sycl::group_barrier(cta);
            S_j[tid] -= projection;
        }
        sycl::group_barrier(cta);
    }
}

template<typename T, int SN>
void applyConstraints(sycl::group<1>& cta, Span<T> S, const Span<T> C, int m){
    //S = S - C @ (C^T @ S)
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
    for(int i = 0; i < SN; i++){
        T* S_i = S.data() + i*m;
        T projection = 0;
        for(int j = 0; j < SN; j++){
            T* C_j = C.data() + i*m;
            projection += sycl::reduce_over_group(cta, S_i[tid]*C_j[tid], sycl::plus<T>{}) * C_j[tid];
        }
        sycl::group_barrier(cta);
        S_i[tid] -= projection;
    }
}

template <typename T, int SN, int BlockVectors>
void update_vectors(sycl::group<1>& cta, Span<T> X, Span<T> R, Span<T> P, Span<T> AX, Span<T> AR, Span<T> AP, Span<T> eigenvects, int m, bool restart){
    auto tid = cta.get_local_id(0);
    group_barrier(cta);
    T Ps[BlockVectors], APs[BlockVectors], Xs[BlockVectors], AXs[BlockVectors];
    for(int i = 0; i < BlockVectors; i++){
        T* X_i = X.data() + i*m; T Xl = X_i[tid];
        T* R_i = R.data() + i*m; T Rl = R_i[tid];
        T* P_i = P.data() + i*m; T Pl = restart ? 0 : P_i[tid];
        T* AX_i = AX.data() + i*m; T AXl = AX_i[tid];
        T* AR_i = AR.data() + i*m; T ARl = AR_i[tid];
        T* AP_i = AP.data() + i*m; T APl = restart ? 0 : AP_i[tid];
        T* eig_X = eigenvects.data() + i*(SN);
        T* eig_R = eigenvects.data() + i*(SN) + BlockVectors;
        T* eig_P = eigenvects.data() + i*(SN) + 2*BlockVectors;
        T tempP = 0;
        T tempAP = 0;
        T tempX = 0;
        T tempAX = 0;
        if (restart){
            for(int k = 0; k < BlockVectors; k++){
                tempP += eig_R[k] * R[k*m + tid];
                tempAP += eig_R[k] * AR[k*m + tid];
            }
        } else {
            for(int k = 0; k < BlockVectors; k++){
                tempP += eig_R[k] * R[k*m + tid] + eig_P[k] * P[k*m + tid];
                tempAP += eig_R[k] * AR[k*m + tid] + eig_P[k] * AP[k*m + tid];
        }}
        for(int k = 0; k < BlockVectors; k++){
            tempX += eig_X[k] * X[k*m + tid];
            tempAX += eig_X[k] * AX[k*m + tid];
        }

        Ps[i] = tempP;
        APs[i] = tempAP;
        Xs[i] = tempX + tempP;
        AXs[i] = tempAX + tempAP;
        //Thread 0 print temporary values
    }
    group_barrier(cta);
    for(int i = 0; i < BlockVectors; i++){
        P[i*m + tid] = Ps[i];
        AP[i*m + tid] = APs[i];
        X[i*m + tid] = Xs[i];
        AX[i*m + tid] = AXs[i];
    }
    group_barrier(cta);
}

template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG_V1(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest){
    static SyclVector<T> S(batch_size * BlockVectors*9 * m);
    static SyclVector<T> X0(batch_size * BlockVectors * m);
    static SyclVector<T> LastEigVects(batch_size * BlockVectors * BlockVectors*3);
    static SyclVector<T> LastGram(batch_size * 3 * BlockVectors * 3 * BlockVectors);
    static SyclVector<T> vals(batch_size * BlockVectors*3);
    static SyclVector<int> indices(batch_size * BlockVectors);

    auto S_acc = S.to_span();
    auto X0_acc = X0.to_span();
    auto LastEigVects_acc = LastEigVects.to_span();
    auto LastGram_acc = LastGram.to_span();
    auto keys_acc = indices.to_span();
    auto vals_acc = vals.to_span();
    auto indices_acc = indices.to_span();

    const T tol = std::numeric_limits<T>::epsilon() * sycl::sqrt(float(m));
    //
    ctx -> submit([&](sycl::handler& h){
        using TupleType = typename std::iterator_traits<oneapi::dpl::zip_iterator<T*, int*>>::value_type;
        auto StAS_smem      = sycl::local_accessor<T, 1>(BlockVectors*3 * BlockVectors*3, h);
        auto eigvects_smem  = sycl::local_accessor<T, 1>(BlockVectors * BlockVectors*3,h);
        auto lambdas_smem   = sycl::local_accessor<T, 1>(BlockVectors, h);
        auto working_space  = sycl::local_accessor<T, 1>(
                        BlockVectors*3 * BlockVectors*3 * 2 + //Lanczos Vectors and transformation matrix Q
                        BlockVectors*3*8                      /*Everything else */, h);
        auto Scache         = sycl::local_accessor<T, 1>(m, h);
        
        const auto bytes = sycl::ext::oneapi::experimental::default_sorters::joint_sorter<>::memory_required<TupleType>(sycl::memory_scope::work_group, BlockVectors*3);
        auto sort_scratchpad = local_accessor<std::byte, 1>(bytes, h);

        h.parallel_for<LOBPCGg<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            sycl::group<1> cta = item.get_group();
            auto bid = item.get_group_linear_id();
            constexpr auto SN = BlockVectors*3;
            std::bitset<BlockVectors> converged;
            oneapi::dpl::uniform_real_distribution<T> distr(0.0, 1.0);            
            oneapi::dpl::minstd_rand engine(42, tid);
            auto A_acc = A.subspan(bid * m * NZ, m * NZ);
            auto cols_acc = cols.subspan(bid * m * NZ, m * NZ);
            auto StAS = Span<T>(StAS_smem.get_pointer(), BlockVectors*3 * BlockVectors*3);
            auto eigvects = Span<T>(eigvects_smem.get_pointer(), BlockVectors * BlockVectors*3);
            auto lambdas = Span<T>(lambdas_smem.get_pointer(), BlockVectors*3);


            //Load the i-th row of A into registers
            T A_tid[NZ];
            for(int i = 0; i < NZ; i++){
                A_tid[i] = A_acc[tid*NZ + i];
            }
            //Load the i-th row of cols into registers
            K cols_tid[NZ]; 
            for(int i = 0; i < NZ; i++){
                cols_tid[i] = cols_acc[tid*NZ + i];
            }


            //X^T @ X  = (X^T @ X)^T = X @ X^T
            //iff A^T = A then 
            //X^T @ A @ X = (X^T @ A @ X)^T
            DECLARE_SUBSPACE

            for (int i = tid + m*BlockVectors; i < 9*m*BlockVectors; i+=cta.get_local_range(0)){
                subspace[i] = T(0);
            }
            for (int i = tid; i < m*BlockVectors; i+=cta.get_local_range(0)){
                blockX[i] = distr(engine);
                //X0_acc[i] = blockX[i];
                X0_acc[bid*m*BlockVectors + i] = blockX[i];
            }
             
            group_barrier(cta);
            //Normalize S vectors
            orthonormalize<T, BlockVectors>(cta, blockX, m); //Modified Gram-Schmidt

            //A * X
            //matBlockVector<T, K, 1, N>(cta, A_tid, cols_tid, blockX, blockAX, NZ, m); 
            //X^T * A * X
            
            STAS<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockX, m, NZ, StAS, Scache.get_pointer());
            //Fill in with simple 3x3 matrix
            compute_k_eigenpairs<T, BlockVectors, BlockVectors>(cta, StAS, eigvects, lambdas, working_space.get_pointer(), sort_scratchpad, largest);
            matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockX, blockAX, NZ, m);
            inPlaceMatMatMul<T, BlockVectors>(cta, blockX, eigvects, m);
            inPlaceMatMatMul<T, BlockVectors>(cta, blockAX, eigvects, m);

            //Compute the residual R = A*X - X*Lambda
            int iter = 0;
            bool restart = true;
            while(!converged.all() && iter < maxiters){
                //R = A*X - X*Lambda
                for(int i = 0; i < BlockVectors; i++) blockR[i*m + tid] = blockAX[i*m + tid] - lambdas[i] * blockX[i*m + tid];
                //Convergence Check
                for(int i = 0; i < BlockVectors; i++) {if(converged[i]) continue; converged[i] = sycl::sqrt(std::abs(reduce_over_group(cta, blockR[i*m + tid]*blockR[i*m + tid], sycl::plus<T>{}))) < tol;}
                sycl::group_barrier(cta);


                //R = R - X * (X^T * R)
                applyConstraints<T, BlockVectors>(cta, blockR, blockX, m);
                sycl::group_barrier(cta);
                orthonormalize<T, BlockVectors>(cta, blockR, m);
                sycl::group_barrier(cta);
                matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockR, blockAR, NZ, m);
                /* 
                 */

                if (!restart) {
                     orthonormalize<T, BlockVectors>(cta, blockP, m);
                    for(int i = 0; i < BlockVectors; i++){ 
                        blockXtemp[i*m + tid] = blockX[i*m + tid];
                        blockRtemp[i*m + tid] = blockR[i*m + tid];
                        blockPtemp[i*m + tid] = blockP[i*m + tid];
                    }
                    matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockP, blockAP, NZ, m);
                    orthonormalize<T, BlockVectors*3>(cta, blockX, m);
                    STAS<T, K, 1, BlockVectors*3>(cta, A_tid, cols_tid, subspace, m, NZ, StAS, Scache.get_pointer());
                    compute_k_eigenpairs<T, BlockVectors*3, BlockVectors>(cta, StAS, eigvects, lambdas, working_space.get_pointer(), sort_scratchpad, largest);
                    
                    for (int i = 0; i < BlockVectors; i++){ 
                        blockX[i*m + tid] = blockXtemp[i*m + tid];
                        blockR[i*m + tid] = blockRtemp[i*m + tid];
                        blockP[i*m + tid] = blockPtemp[i*m + tid];
                    } 
                    group_barrier(cta);
                    update_vectors<T, BlockVectors*3, BlockVectors>(cta, blockX, blockR, blockP, blockAX, blockAR, blockAP, eigvects, m, false);
                    group_barrier(cta);

                } else {
                    for (int i = 0; i < BlockVectors; i++){
                        blockXtemp[i*m + tid] = blockX[i*m + tid];
                        blockRtemp[i*m + tid] = blockR[i*m + tid];
                    }
                    orthonormalize<T, BlockVectors*2>(cta, blockX, m);
                    STAS<T, K, 1, BlockVectors*2>(cta, A_tid, cols_tid, subspace, m, NZ, StAS, Scache.get_pointer());
                    compute_k_eigenpairs<T, BlockVectors*2, BlockVectors>(cta, StAS, eigvects, lambdas, working_space.get_pointer(), sort_scratchpad, largest);
                    for (int i = 0; i < BlockVectors; i++){
                        blockX[i*m + tid] = blockXtemp[i*m + tid];
                        blockR[i*m + tid] = blockRtemp[i*m + tid];
                    }
                    group_barrier(cta);
                    update_vectors<T, BlockVectors*2, BlockVectors>(cta, blockX, blockR, blockP, blockAX, blockAR, blockAP, eigvects, m, true);
                    restart = false;
                }
                sycl::group_barrier(cta);
                iter++;
            }
            for (int i = tid; i < BlockVectors*BlockVectors*3; i+=cta.get_local_range(0)){
                LastEigVects_acc[i] = eigvects[i];
            }
            for (int i = tid; i < 3*BlockVectors*3*BlockVectors; i+=cta.get_local_range(0)){
                LastGram_acc[i] = StAS[i];
            }
            if(tid < BlockVectors*3)   vals_acc[tid] = lambdas[tid];
        });
    });
    std::string dtype = std::is_same<T, float>::value ? ".float32" : ".float64";
    ctx.wait();
    /* int SD = std::min(int(maxiters+1), 3);
    std::ofstream file("Subspace_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        file.write(reinterpret_cast<const char*>(S_acc.data()), SD * BlockVectors * m * sizeof(T));
    }

    std::ofstream file2("AXARAP_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        file2.write(reinterpret_cast<const char*>(S_acc.data() + 6 * BlockVectors * m), 3 * BlockVectors * m * sizeof(T));
    }

    std::ofstream file3("Initial_X_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + dtype, ios::binary);
    {
        file3.write(reinterpret_cast<const char*>(X0_acc.data()), BlockVectors * m * sizeof(T));
    }

    std::ofstream file4("LastEigVects_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        file4.write(reinterpret_cast<const char*>(LastEigVects_acc.data()), BlockVectors * BlockVectors * SD * sizeof(T));
    }

    std::ofstream file5("LastGram_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        file5.write(reinterpret_cast<const char*>(LastGram_acc.data()), SD * BlockVectors * SD * BlockVectors * sizeof(T));
    }

    std::ofstream file6("Lambdas_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        file6.write(reinterpret_cast<const char*>(vals_acc.data()), BlockVectors * sizeof(T));
    }
    */
    std::cout << vals_acc << std::endl; 
}

template void LOBPCG_V1<float, uint16_t, 3, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest);
template void LOBPCG_V1<float, uint16_t, 21, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest);

//template <typename T, typename K, int BlockVectors, int NZ>
//void LOBPCG_V3(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest){}
//
//template void LOBPCG_V3<float, uint16_t, 21, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest);
//template void LOBPCG_V3<float, uint16_t, 3, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest);

//template void LOBPCG<float, uint16_t, 9, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters);
//template void LOBPCG<double, uint16_t, 21, 30>(SyclQueue &ctx, Span<double> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters);
//template void LOBPCG<double, uint16_t, 3, 30>(SyclQueue &ctx, Span<double> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters);
//template void LOBPCG<double, uint16_t, 9, 30>(SyclQueue &ctx, Span<double> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters);

//template void LOBPCG_V1<float, uint16_t, 21, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters);