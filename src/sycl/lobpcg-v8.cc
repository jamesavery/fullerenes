#include <array>
#include <algorithm>
#include <fstream>
#include <sycl/sycl.hpp>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <oneapi/dpl/random>
#include <oneapi/dpl/iterator>
#include <oneapi/mkl/lapack.hpp>
#include <oneapi/mkl/spblas.hpp>
#include <fullerenes/graph.hh>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include "queue-impl.cc"
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/argparser.hh>
#include <fullerenes/sycl-headers/fill.hh>
#include <cusolverDn.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include "primitives.cc"
#include <iomanip>
#include <fullerenes/mempool.hh>
#include "../linalg/linalg-impl.hh"
#include "linalg/batched/QR/chol_qr.cc"


#define CEIL_DIV(x, y) (((x) + (y) - 1) / (y))
#define DECLARE_SUBSPACE\
    auto subspace = S_acc.subspan(bid * BlockVectors*3 * m, BlockVectors * m * 3);\
    auto blockX = S_acc.subspan(bid * BlockVectors*3 * m + 0*BlockVectors*m, BlockVectors*m);\
    auto blockP = S_acc.subspan(bid * BlockVectors*3 * m + 1*BlockVectors*m, BlockVectors*m);\
    auto blockR = S_acc.subspan(bid * BlockVectors*3 * m + 2*BlockVectors*m, BlockVectors*m);\
    auto blockXtemp = S_temp_acc.subspan(bid * BlockVectors*3 * m + 0*BlockVectors*m, BlockVectors*m);\
    auto blockPtemp = S_temp_acc.subspan(bid * BlockVectors*3 * m + 1*BlockVectors*m, BlockVectors*m);\
    auto blockRtemp = S_temp_acc.subspan(bid * BlockVectors*3 * m + 2*BlockVectors*m, BlockVectors*m);\
    auto blockAX = AS_acc.subspan(bid * BlockVectors*3 * m + 0*BlockVectors*m, BlockVectors*m);\
    auto blockAP = AS_acc.subspan(bid * BlockVectors*3 * m + 1*BlockVectors*m, BlockVectors*m);\
    auto blockAR = AS_acc.subspan(bid * BlockVectors*3 * m + 2*BlockVectors*m, BlockVectors*m);

#define CHECK_CUSPARSE(func)                                                   \
{                                                                              \
    cusparseStatus_t cusparse_status = (func);                                          \
    if (cusparse_status != CUSPARSE_STATUS_SUCCESS) {                                   \
        printf("CUSPARSE API failed at line %d with error: %s (%d)\n",         \
               __LINE__, cusparseGetErrorString(cusparse_status), cusparse_status);              \
        assert(false);                                                   \
    }                                                                          \
}



template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V8_1 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V8_2 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V8_3 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V8_4 {};

using namespace sycl;

template <typename T, int SN>
void inPlaceMatMatMul(sycl::group<1>& cta, Span<T> A, Span<T> B, int m, bool largest){
    auto tid = cta.get_local_id(0);
    T Alocal[SN];
    for(int j = 0; j < SN; j++){
        auto offset = SN * (largest ? (SN -j - 1) : j);
        T sum = 0;
        for(int k = 0; k < SN; k++){
            sum += A[k*m + tid] * B[offset + k];
        }
        Alocal[j] = sum;
    }
    sycl::group_barrier(cta);
    for(int j = 0; j < SN; j++){
        A[j*m + tid] = Alocal[j];
    }   
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
void SVQB(SyclQueue &ctx,   Span<T> S /* Subspace */,
                            Span<T> StS /* Subspace Gram Matrix*/,
                            Span<T> D /* D = diag(StS)^-1/2 */,
                            Span<T> U /* Output Orthonormal Basis */,
                            size_t m,
                            size_t Sstride,
                            size_t batch_size){
    //Compute StS = S^T * S
    static cublasHandle_t handle;
    static cusolverDnHandle_t cusolver_handle;
    static cusolverDnParams_t params;
    static SyclVector<T> d_solver_scratchpad;
    static SyclVector<T> h_solver_scratchpad;
    static SyclVector<T> lambdas;
    static SyclVector<int> d_info(batch_size);
    static size_t d_scratchpad_size = 0;
    static size_t h_scratchpad_size = 0;
    constexpr T alpha = 1.0;
    constexpr T beta = 0.0;
    static size_t stored_SN = 0;
    if (stored_SN != SN){
        cublasCreate(&handle);
        cusolverDnCreate(&cusolver_handle);
        cusolverDnCreateParams(&params);
        
        lambdas = SyclVector<T>(batch_size * SN);
        cusolverDnXsyevBatched_bufferSize(cusolver_handle,
                                            params,
                                            CUSOLVER_EIG_MODE_VECTOR,
                                            CUBLAS_FILL_MODE_LOWER,
                                            SN,
                                            CUDA_R_32F,
                                            StS.data(),
                                            SN,
                                            CUDA_R_32F,
                                            lambdas.data(),
                                            CUDA_R_32F,
                                            &d_scratchpad_size,
                                            &h_scratchpad_size,
                                            batch_size);
        d_solver_scratchpad = SyclVector<T>(d_scratchpad_size);
        h_solver_scratchpad = SyclVector<T>(h_scratchpad_size);
        stored_SN = SN;
    }

    //Compute StS = S^T * S
    cublasSgemmStridedBatched(handle, CUBLAS_OP_T, CUBLAS_OP_N,
                                SN, /* m */
                                SN, /* n */
                                m, /* k */
                                &alpha, /* alpha */
                                S.data(), /* A */
                                m, /* lda */
                                Sstride, /* strideA */
                                S.data(), /* B */
                                m, /* ldb */
                                Sstride, /* strideB */
                                &beta, /* beta */
                                StS.data(), /* C */
                                SN, /* ldc */
                                SN*SN, /* strideC */
                                batch_size);
    cudaDeviceSynchronize();

    //Compute D = diag(StS)^-1/2
    ctx -> submit([&](sycl::handler& h){
        h.parallel_for(nd_range<1>(sycl::range{size_t(batch_size*SN)}, sycl::range{size_t(SN)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            auto bid = item.get_group_linear_id();
            auto D_acc = D.subspan(bid * SN, SN);
            auto StS_acc = StS.subspan(bid * SN * SN, SN * SN);
            D_acc[tid] = sycl::rsqrt(StS_acc[tid * SN + tid]);
        });
    });

    //Compute StS = D * StS * D
    ctx -> submit([&](sycl::handler& h){
        auto D_local = sycl::local_accessor<T, 1>(SN, h);
        h.parallel_for(nd_range<1>(sycl::range{size_t(batch_size*SN)}, sycl::range{size_t(SN)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            auto bid = item.get_group_linear_id();
            auto cta = item.get_group();
            auto D_acc = D.subspan(bid * SN, SN);
            D_local[tid] = D_acc[tid];
            sycl::group_barrier(cta);
            auto StS_acc = StS.subspan(bid * SN * SN, SN * SN);
            auto D_tid = D_local[tid];
            for(int i = 0; i < SN; i++){
                auto D_i = D_local[i];
                StS_acc[i * SN + tid] *= D_tid * D_i;
            }
        });
    });

    ctx.wait();

    //Compute the eigenvalues of D StS D
    auto solve_status = cusolverDnXsyevBatched(cusolver_handle,
                            params,
                            CUSOLVER_EIG_MODE_VECTOR,
                            CUBLAS_FILL_MODE_LOWER,
                            SN,
                            CUDA_R_32F,
                            StS.data(),
                            SN,
                            CUDA_R_32F,
                            lambdas.data(),
                            CUDA_R_32F,
                            d_solver_scratchpad.data(),
                            d_scratchpad_size,
                            h_solver_scratchpad.data(),
                            h_scratchpad_size,
                            d_info.data(),
                            batch_size);

    cudaDeviceSynchronize();
    //Compute Q = S * D * EigenVectors * Lambda^-1/2

    //First Compute D * EigenVectors * Lambda^-1/2
    ctx -> submit([&](sycl::handler& h){
        auto D_local = sycl::local_accessor<T, 1>(SN, h);
        auto lambdas_acc = lambdas.to_span();
        h.parallel_for(nd_range<1>(sycl::range{size_t(batch_size*SN)}, sycl::range{size_t(SN)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            auto bid = item.get_group_linear_id();
            auto cta = item.get_group();
            auto D_acc = D.subspan(bid * SN, SN);
            D_local[tid] = D_acc[tid];
            sycl::group_barrier(cta);
            auto StS_acc = StS.subspan(bid * SN * SN, SN * SN);
            //auto lambdas_i = lambdas_acc[bid * SN + tid];
            auto D_tid = D_local[tid];
            auto tau = std::numeric_limits<T>::epsilon() * std::abs(lambdas_acc[bid * SN + SN - 1]);
            for(int i = 0; i < SN; i++){
                //auto D_i = D_local[i];
                auto lambda_i = lambdas_acc[bid * SN + i] < tau ? tau : lambdas_acc[bid * SN + i];
                StS_acc[i * SN + tid] *= D_tid * sycl::rsqrt(std::abs(lambda_i));
            }
        });
    });
    ctx.wait();
    
    //Compute Q = S * D * EigenVectors * Lambda^-1/2
    cublasSgemmStridedBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                m, /* m */
                                SN, /* n */
                                SN, /* k */
                                &alpha, /* alpha */
                                S.data(), /* A */ 
                                m, /* lda */
                                m*SN, /* StrideA */
                                StS.data(), /* B */
                                SN, /* ldb */
                                SN*SN, /* StrideB */
                                &beta, /* beta */
                                U.data(), /* C */
                                m, /* ldc */
                                Sstride, /* StrideC */
                                batch_size);
    
    cudaDeviceSynchronize();
}

template <typename T, int ExtDim /* Number of vectors in External Basis */, int CandDim /* Number of vectors in Candidate Basis */>
void ExternalOrthogonalization(SyclQueue& ctx,  Span<T> S  /* Candidate Basis [k,CandDim] */, 
                                                Span<T> E  /* External Basis [k,ExtDim]*/,
                                                Span<T> U  /* Workspace [ExtDim * CandDim] */,
                                                size_t k,
                                                size_t ExtStride,
                                                size_t CandStride,
                                                size_t batch_size){
    static cublasHandle_t handle;
    static bool initialized = false;
    if (!initialized){
        cublasCreate(&handle);
    }

    constexpr T alpha0 = 1.0;
    constexpr T beta0 = 0.0;
    //Compute U = E^T * S
    cublasSgemmStridedBatched(  handle, CUBLAS_OP_T, CUBLAS_OP_N,
                                ExtDim, /* m */
                                CandDim, /* n */
                                k, /* k */
                                &alpha0, /* alpha */
                                E.data(), /* A */
                                k, /* lda */
                                ExtStride, /* strideA */
                                S.data(), /* B */
                                k, /* ldb */
                                CandStride, /* strideB */
                                &beta0, /* beta */
                                U.data(), /* C */
                                ExtDim, /* ldc */
                                ExtDim*CandDim, /* strideC */
                                batch_size);
    cudaDeviceSynchronize();

    //Compute S = S - E * (E^T * S)
    constexpr T alpha1 = -1.0;
    constexpr T beta1 = 1.0;

    cublasSgemmStridedBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                k, /* m */
                                CandDim, /* n */
                                ExtDim, /* k */
                                &alpha1, /* alpha */
                                E.data(), /* A */
                                k, /* lda */
                                ExtStride, /* strideA */
                                U.data(), /* B */
                                ExtDim, /* ldb */
                                ExtDim*CandDim, /* strideB */
                                &beta1, /* beta */
                                S.data(), /* C */
                                k, /* ldc */
                                CandStride, /* strideC */
                                batch_size);
    cudaDeviceSynchronize();
}

template <typename T, int ExtDim, int CandDim>
void Ortho(SyclQueue& ctx, Span<T> S /* Candidate Basis [k,CandDim] */, 
                            Span<T> E /* External Basis [k,ExtDim] */,
                            Span<T> U /* Workspace [ExtDim * CandDim] */,
                            size_t k,
                            size_t ExtStride,
                            size_t CandStride,
                            size_t batch_size){
    for(int i = 0; i < 2; i++){
        ExternalOrthogonalization<T, ExtDim, CandDim>(ctx, S, E, U, k, ExtStride, CandStride, batch_size);
        //Chol2QR<T, CandDim>(ctx, S_ptrs, StS_ptrs, k, batch_size);
        chol2_qr_batched(ctx, false, batch_size, CandStride, k, CandDim, S, U.template as_span<std::byte>());
    }
}
                            

template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG_V8(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest, Span<T> out_eigvects, Span<T> out_eigvals){
    SyclVector<T> S(batch_size * BlockVectors*3 * m);
    SyclVector<T> S_temp(batch_size * BlockVectors*3 * m);
    SyclVector<T> AS(batch_size * BlockVectors*3 * m);
    SyclVector<T> X0(batch_size * BlockVectors * m);
    SyclVector<T> LastEigVects(batch_size * BlockVectors * BlockVectors*3);
    SyclVector<T> LastGram(batch_size * 3 * BlockVectors * 3 * BlockVectors);
    SyclVector<T> vals(batch_size * BlockVectors);
    SyclVector<T> StAS(batch_size * BlockVectors*3 * BlockVectors*3, 0);
    auto S_acc = S.to_span(); 
    auto S_temp_acc = S_temp.to_span();
    auto AS_acc = AS.to_span();
    auto X0_acc = X0.to_span();
    auto StAS_acc = StAS.to_span();
    constexpr T alpha = 1.0;
    constexpr T beta = 0.0;
    auto SubspaceStride = m*BlockVectors*3;

    SyclVector<int> cols_32I(cols.size());
    primitives::transform(ctx, cols, cols_32I, [](K i){return static_cast<int>(i);});
    
    size_t syevd_scratchpad_size = 0;
    size_t syevd_scratchpad_host_size = 0;
    SyclVector<int> syevd_info(batch_size);
    SyclVector<T> lambdas(batch_size * BlockVectors*3);
    SyclVector<int> csr_row_offsets(m + 1);
    
    //csr_row_offsets = [0, NZ, 2*NZ, 3*NZ, ...]
    primitives::iota(ctx, csr_row_offsets, 0);
    primitives::transform(ctx, csr_row_offsets, csr_row_offsets, [&](int i){return i * NZ;}); 

    auto lambda_acc = lambdas.to_span();
    
    const T tol = sycl::sqrt(std::numeric_limits<T>::epsilon()) * 1e-2 * T(m);
    
    SyclVector<T> Diag(3*BlockVectors * batch_size, 0);
    SyclVector<T> R_inv(batch_size * 3 * BlockVectors * 3 * BlockVectors, 0);
    SyclVector<T> U(batch_size * BlockVectors * 3 * m, 0);
    SyclVector<T> C_p(batch_size * BlockVectors * BlockVectors * 3, 0); //At most it will contain BlockVectors eigenvectors of size 3*BlockVectors
    SyclVector<T> QRworkspace(batch_size * BlockVectors * 3 * m, 0);
    //Setup
    ctx -> submit([&](sycl::handler& h ){
        auto Scache = sycl::local_accessor<T, 1>(m, h);

        h.parallel_for<LOBPCG_V8_1<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            sycl::group<1> cta = item.get_group();
            auto bid = item.get_group_linear_id();
            oneapi::dpl::uniform_real_distribution<T> distr(0.0, 1.0);            
            oneapi::dpl::minstd_rand engine(42, tid);
            DECLARE_SUBSPACE;

            sycl::group_barrier(cta);
            for (int i = tid; i < m*BlockVectors; i+=cta.get_local_range(0)){
                blockX[i] = distr(engine);
                //X0_acc[i] = blockX[i];
                X0_acc[bid*m*BlockVectors + i] = blockX[i];
            }
            //Normalize S vectors
            //orthonormalize<T, BlockVectors>(cta, blockX, m); //Modified Gram-Schmidt
            
        });
    });
    ctx.wait();
    //CholQR<T, BlockVectors>(ctx, S_ptrs, STS_ptrs, m, batch_size);
    chol2_qr_batched(ctx, false, batch_size, m*3*BlockVectors, m, BlockVectors, S.to_span(), QRworkspace.to_span().template as_span<std::byte>());
    linalg::SparseMatHandle<T, Format::CSR, BatchType::Batched> handleA(A.data(), csr_row_offsets.data(), cols_32I.data(), m*NZ, m, m, m*NZ, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleS(S.data(), m, BlockVectors*3, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleX(S.data(), m, BlockVectors, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleR(S.data() + 2*BlockVectors*m, m, BlockVectors, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleP(S.data() + BlockVectors*m, m, BlockVectors, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleAS(AS.data(), m, BlockVectors*3, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleAX(AS.data(), m, BlockVectors, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleAR(AS.data() + 2*BlockVectors*m, m, BlockVectors, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleAP(AS.data() + BlockVectors*m, m, BlockVectors, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleXtAX(StAS.data(), BlockVectors, BlockVectors, BlockVectors, BlockVectors*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleStAS(StAS.data(), BlockVectors*3, BlockVectors*3, BlockVectors*3, BlockVectors*3*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleS_restart(S.data(), m, BlockVectors*2, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleAS_restart(AS.data(), m, BlockVectors*2, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleStAS_restart(StAS.data(), BlockVectors*2, BlockVectors*2, BlockVectors*2, BlockVectors*2*2*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleU(U.data(), m, BlockVectors*3, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleStemp(S_temp.data(), m, BlockVectors*3, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleC_p(C_p.data(), BlockVectors*3, BlockVectors, BlockVectors*3, BlockVectors*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleC_p_restart(C_p.data(), BlockVectors*2, BlockVectors, BlockVectors*2, BlockVectors*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleP_new(U.data() + m*BlockVectors, m, BlockVectors, m, m*3*BlockVectors, batch_size);
    linalg::DenseMatHandle<T, BatchType::Batched> handleAP_new(S_temp.data() + m*BlockVectors, m, BlockVectors, m, m*3*BlockVectors, batch_size);
    SyclVector<std::byte> spmm_workspace(linalg::spmm_buffer_size<Backend::CUDA>(ctx, handleA, handleS, handleAS, alpha, beta, Transpose::NoTrans, Transpose::NoTrans));
    SyclVector<std::byte> syev_workspace(linalg::syev_buffer_size<Backend::CUDA>(ctx, handleStAS, lambdas.to_span(), Uplo::Lower));


    //Compute AX
    linalg::spmm<Backend::CUDA>(ctx, handleA, handleX, handleAX, alpha, beta, Transpose::NoTrans, Transpose::NoTrans, spmm_workspace.to_span());
    
    //Compute X^T AX
    linalg::gemm<Backend::CUDA>(ctx, handleX, handleAX, handleXtAX, alpha, beta, Transpose::Trans, Transpose::NoTrans);

    //Solve the eigenvalue problem
    linalg::syev<Backend::CUDA>(ctx, handleXtAX, lambdas.to_span(), Uplo::Lower, syev_workspace.to_span());
    ctx.wait();


    ctx -> submit([&](sycl::handler& h){
        auto Scache = sycl::local_accessor<T, 1>(m, h);

        h.parallel_for<LOBPCG_V8_2<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            sycl::group<1> cta = item.get_group();
            auto bid = item.get_group_linear_id();

            auto A_acc = A.subspan(bid * m * NZ, m * NZ);
            auto cols_acc = cols.subspan(bid * m * NZ, m * NZ);

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

            DECLARE_SUBSPACE;

            matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockX, blockAX, NZ, m);
            //cuSolver stores the eigenvectors in-place in the input matrix
            auto eigvects = StAS_acc.subspan(bid * BlockVectors * BlockVectors, BlockVectors * BlockVectors);
            inPlaceMatMatMul<T, BlockVectors>(cta, blockX, eigvects, m , largest);
            inPlaceMatMatMul<T, BlockVectors>(cta, blockAX, eigvects, m, largest);
            //Compute the residual R = A*X - X*Lambda
        });
    });
    ctx->wait();



    SyclVector<std::bitset<BlockVectors>> converged(batch_size);
    SyclVector<T> residuals(batch_size * BlockVectors);
    auto converged_acc = converged.to_span();
    auto residuals_acc = residuals.to_span();

    auto iteration = [&](SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t iter, bool restart){
        
        auto start = std::chrono::high_resolution_clock::now();
        //Setup
        ctx -> submit([&](sycl::handler& h){
            auto Scache = sycl::local_accessor<T, 1>(m, h);

            h.parallel_for<LOBPCG_V8_3<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
                auto tid = item.get_local_linear_id();
                auto cta = item.get_group();
                auto bid = item.get_group_linear_id();
                auto& converged = converged_acc[bid];
                auto residuals = residuals_acc.subspan(bid * BlockVectors, BlockVectors);
                auto num_eigvals = iter < 2 ? (iter+1) * BlockVectors : 3*BlockVectors;
                auto lambdas = lambda_acc.subspan(bid * num_eigvals, num_eigvals);
                DECLARE_SUBSPACE;

                //R = A*X - X*Lambda
                for(int i = 0; i < BlockVectors; i++) blockR[i*m + tid] = blockAX[i*m + tid] - lambdas[largest ? (num_eigvals - 1 - i) : i] * blockX[i*m + tid];
                //Convergence Check
                for(int i = 0; i < BlockVectors; i++) {
                    //if(converged[i]) continue; 
                    residuals[i] = sycl::sqrt(reduce_over_group(cta, blockR[i*m + tid]*blockR[i*m + tid], sycl::plus<T>{}));
                    residuals[i] /= sycl::sqrt(reduce_over_group(cta, blockX[i*m + tid]*blockX[i*m + tid], sycl::plus<T>{})) * lambdas[largest ? (num_eigvals - 1 - i) : i];
                    converged[i] = residuals[i] < tol;
                }
            });
        });
        ctx -> wait();
        
        start = std::chrono::high_resolution_clock::now();
        if (restart){
            Ortho<T, BlockVectors, BlockVectors>(ctx, S.subspan(2*m*BlockVectors), S, U, m, SubspaceStride, SubspaceStride, batch_size);
        } else {
            Ortho<T, 2*BlockVectors, BlockVectors>(ctx, S.subspan(2*m*BlockVectors), S, U, m, SubspaceStride, SubspaceStride, batch_size);
        }
        //Strided Copy Kernel:
        if (restart){
            ctx -> submit([&](sycl::handler& h){
                auto src = S.to_span();
                auto dst = S.to_span();
                h.parallel_for(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
                    auto tid = item.get_local_linear_id();
                    auto bid = item.get_group_linear_id();
                    auto cta = item.get_group();
                    auto block_src = src.subspan((bid * 3 + 2) * m * BlockVectors, m * BlockVectors);
                    auto block_dst = dst.subspan((bid * 3 + 1) * m * BlockVectors, m * BlockVectors);
                    for(int i = tid; i < m*BlockVectors; i+=cta.get_local_range(0)){
                        block_dst[i] = block_src[i];
                    }
                });
            });
            ctx.wait();
        }
        //Compute AR
        linalg::spmm<Backend::CUDA>(ctx, handleA, restart ? handleP : handleR, restart ? handleAP : handleAR, alpha, beta, Transpose::NoTrans, Transpose::NoTrans, spmm_workspace.to_span());
        //Compute S^T A S
        linalg::gemm<Backend::CUDA>(ctx, restart ? handleS_restart : handleS, restart ? handleAS_restart : handleAS, restart ? handleStAS_restart : handleStAS, alpha, beta, Transpose::Trans, Transpose::NoTrans);

        //Solve the eigenvalue problem
        linalg::syev<Backend::CUDA>(ctx, restart ? handleStAS_restart : handleStAS, lambdas.to_span(), Uplo::Lower, syev_workspace.to_span());
        ctx.wait();
        //cudaMemcpy(R_inv.data(), StAS.data(), StAS.size() * sizeof(T), cudaMemcpyDeviceToDevice);
        R_inv = StAS;
        //cudaDeviceSynchronize();

        //If largest = true, then the order of the eigenvectors is reversed
        if (largest){
            ctx -> submit([&](sycl::handler& h){
                auto cols = restart ? BlockVectors*2 : BlockVectors*3;
                auto R_inv_acc = R_inv.to_span();
                h.parallel_for(nd_range<1>(sycl::range{size_t(batch_size * 256)}, sycl::range{size_t(256)}), [=](sycl::nd_item<1> item){
                    auto tid = item.get_local_linear_id();
                    auto bid = item.get_group_linear_id();
                    auto cta = item.get_group();
                    auto dst_acc = StAS_acc.subspan(bid * cols * cols, cols * cols);
                    auto src_acc = R_inv_acc.subspan(bid * cols * cols, cols * cols);
                    for(int i = tid; i < cols*cols; i+=cta.get_local_range(0)){
                        auto src_ix = i % cols + (cols - 1 - i / cols)*cols ; 
                        dst_acc[i] = src_acc[src_ix];
                    }
                });
            });
            ctx -> wait();
        }

        //X(i+1) =  [X(i), R(i), P(i)] * [Zx, Zr, Zp]^T

        //Compute C_p = C_x - [I 
        //                     0 
        //                     0]

        ctx -> submit([&](sycl::handler& h){
            auto C_p_acc = C_p.to_span();
            auto C_x_acc = StAS.to_span();
            auto Cdim = BlockVectors * (restart ? 2 : 3);
            auto Cstride = BlockVectors * BlockVectors * 3;
            auto Nactive = BlockVectors; //When we start soft-locking the eigenvectors, we will need to update this
            h.parallel_for(nd_range<1>(sycl::range{size_t(batch_size*Nactive)}, sycl::range{size_t(Nactive)}), [=](sycl::nd_item<1> item){
                auto cta = item.get_group();
                auto tid = item.get_local_linear_id();
                auto bid = item.get_group_linear_id();

                auto C_x = C_x_acc.subspan(bid * Cdim * Cdim, Cdim * Cdim);
                auto C_p = C_p_acc.subspan(bid * Cstride, Cdim * Nactive);

                for(int i = tid; i < Nactive * Cdim; i+=cta.get_local_range(0)) {C_p[i] = C_x[i];}
                
                sycl::group_barrier(cta);
                C_p[tid * Cdim + tid] -= 1;
            });
        });
        ctx -> wait();
        //Compute X = [X, R, P] * [Zx, Zr, Zp]^T
        linalg::gemm<Backend::CUDA>(ctx, restart ? handleS_restart : handleS, restart ? handleStAS_restart : handleStAS, handleU, alpha, beta, Transpose::NoTrans, Transpose::NoTrans); 
        //Make an implicit update of AX:
        linalg::gemm<Backend::CUDA>(ctx, restart ? handleAS_restart : handleAS, restart ? handleStAS_restart : handleStAS, handleStemp, alpha, beta, Transpose::NoTrans, Transpose::NoTrans);
        ctx.wait();
        Ortho<T, BlockVectors, BlockVectors>(ctx, C_p, StAS, QRworkspace, restart ? BlockVectors*2 : BlockVectors*3, BlockVectors*BlockVectors * (restart? (2*2) : (3*3)), BlockVectors*BlockVectors *  3, batch_size);
        //First compute P = [X, R, P] * C_p

        linalg::gemm<Backend::CUDA>(ctx, restart ? handleS_restart : handleS, restart ? handleC_p_restart : handleC_p, handleP_new, alpha, beta, Transpose::NoTrans, Transpose::NoTrans);
        //Make an implicit update of AP:
        linalg::gemm<Backend::CUDA>(ctx, restart ? handleAS_restart : handleAS, restart ? handleC_p_restart : handleC_p, handleAP_new, alpha, beta, Transpose::NoTrans, Transpose::NoTrans);
        ctx.wait();
        //AS = S_temp;
        //AS.swap(S_temp);
        AS = S_temp;
        S = U;
    };

    auto T0 = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < maxiters; i++){
        iteration(ctx, A, cols, batch_size, m, i, i < 1 ? true : false);
    }
    ctx.wait();   
    for(int i = 0; i < batch_size; i++){
        primitives::copy(ctx, lambdas.subspan(i * (BlockVectors *3) + BlockVectors * 2, BlockVectors), out_eigvals.subspan(i * BlockVectors, BlockVectors));
        //primitives::copy(ctx, S_acc.subspan(i * BlockVectors * 3 * m, BlockVectors * m), eigvects.subspan(i * BlockVectors * m, BlockVectors * m));
    }
}

template void LOBPCG_V8<float, uint16_t, 3, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest, Span<float> eigvects, Span<float> eigvals);
template void LOBPCG_V8<float, uint16_t, 21, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest, Span<float> eigvects, Span<float> eigvals);