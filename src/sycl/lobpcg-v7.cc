#include <array>
#include <algorithm>
#include <fstream>
#include <sycl/sycl.hpp>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <oneapi/dpl/random>
#include <oneapi/dpl/iterator>
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



template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V7_1 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V7_2 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V7_3 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V7_4 {};

using namespace sycl;

template <typename T>
void print_matrix(Span<T> matrix, int rows, int cols, bool transpose = false) {
    std::cout << std::fixed << std::setprecision(6);
    if (!transpose) {
        std::cout << "[\n";
        for (int i = 0; i < rows; i++) {
            std::cout << "  [";
            for (int j = 0; j < cols; j++) {
                std::cout << std::setw(10) << matrix[i * cols + j];
                if (j < cols - 1) std::cout << ", ";
            }
            std::cout << "]";
            if (i < rows - 1) std::cout << ",";
            std::cout << "\n";
        }
        std::cout << "]\n";
    } else {
        std::cout << "[\n";
        for (int i = 0; i < cols; i++) { // Outer loop is cols for correct transposition
            std::cout << "  [";
            for (int j = 0; j < rows; j++) {
                std::cout << std::setw(10) << matrix[j * cols + i]; // Transposed indexing
                if (j < rows - 1) std::cout << ", ";
            }
            std::cout << "]";
            if (i < cols - 1) std::cout << ",";
            std::cout << "\n";
        }
        std::cout << "]\n";
    }   
}

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

template <typename T, int SN>
void ComputeGramMatrix(SyclQueue& ctx, Span<T> A,
                                        Span<T> Gram,
                                        size_t m,
                                        size_t Astride,
                                        size_t batch_size){
    static cublasHandle_t handle;
    static size_t stored_SN = 0;
    if (stored_SN != SN){
        cublasCreate(&handle);
        stored_SN = SN;
    }
    constexpr T alpha = 1.0;
    constexpr T beta = 0.0;
    cublasSgemmStridedBatched(handle, CUBLAS_OP_T, CUBLAS_OP_N,
                                SN, /* m */
                                SN, /* n */
                                m, /* k */
                                &alpha, /* alpha */
                                A.data(), /* A */
                                m, /* lda */
                                Astride, /* strideA */
                                A.data(), /* B */
                                m, /* ldb */
                                Astride, /* strideB */
                                &beta, /* beta */
                                Gram.data(), /* C */
                                SN, /* ldc */
                                SN*SN, /* strideC */
                                batch_size);
    cudaDeviceSynchronize();
}


template <typename T, int SN>
void CholQR(SyclQueue &ctx, Span<T*> S_ptrs /* Subspace */,
                            Span<T*> STS_ptrs /* Subspace Gram Matrix*/,
                            size_t m,
                            size_t batch_size){
    static cublasHandle_t blas_handle;
    static cusolverDnHandle_t solver_handle;
    static SyclVector<int> d_info(batch_size);    
    static bool initialized = false;
    if (!initialized){
        cublasCreate(&blas_handle);
        cusolverDnCreate(&solver_handle);
        initialized = true;
    }

    auto S_stride = batch_size > 1 ? std::distance(S_ptrs[0], S_ptrs[1]) : m*SN;
    auto STS_stride = batch_size > 1 ? std::distance(STS_ptrs[0], STS_ptrs[1]) : SN*SN;

    constexpr T alpha = 1.0;
    constexpr T beta = 0.0;

    //Compute StS = S^T * S
    cublasSgemmStridedBatched(blas_handle, CUBLAS_OP_T, CUBLAS_OP_N,
                                SN, /* m */
                                SN, /* n */
                                m, /* k */
                                &alpha, /* alpha */
                                S_ptrs.front(), /* A */
                                m, /* lda */
                                S_stride, /* strideA */
                                S_ptrs.front(), /* B */
                                m, /* ldb */
                                S_stride, /* strideB */
                                &beta, /* beta */
                                STS_ptrs.front(), /* C */
                                SN, /* ldc */
                                STS_stride, /* strideC */
                                batch_size);


    cudaDeviceSynchronize();

    //Compute the Cholesky Factorization of StS
    cusolverDnSpotrfBatched(solver_handle, CUBLAS_FILL_MODE_LOWER, SN, STS_ptrs.data(), SN, d_info.data(), batch_size);
    //cublasSgetrfBatched(blas_handle, SN, Sts_ptrs.data(), SN, NULL, d_info.data(), batch_size);
    cudaDeviceSynchronize();
    //Compute Q = S * StS^-1
    cublasStrsmBatched( blas_handle,
                        CUBLAS_SIDE_RIGHT, 
                        CUBLAS_FILL_MODE_LOWER, 
                        CUBLAS_OP_T, 
                        CUBLAS_DIAG_NON_UNIT, 
                        m, 
                        SN, 
                        &alpha, 
                        STS_ptrs.data(), 
                        SN, 
                        S_ptrs.data(), 
                        m, batch_size);
    cudaDeviceSynchronize();
}

template <typename T, int SN>
void Chol2QR(SyclQueue &ctx, Span<float*> S_ptrs /* Subspace */,
                            Span<float*> STS_ptrs /* Subspace Gram Matrix*/,
                            size_t m,
                            size_t batch_size){
    CholQR<T,SN>(ctx, S_ptrs, STS_ptrs, m, batch_size);
    CholQR<T,SN>(ctx, S_ptrs, STS_ptrs, m, batch_size);
}

template <typename T, int SN>
void ShiftCholQR3(SyclQueue &ctx, Span<T*> S_ptrs /* Subspace */,
                            Span<T*> STS_ptrs /* Subspace Gram Matrix*/,
                            size_t m,
                            size_t batch_size){
    static cublasHandle_t blas_handle;
    static cusolverDnHandle_t solver_handle;
    static SyclVector<int> d_info(batch_size);    
    static bool initialized = false;
    if (!initialized){
        cublasCreate(&blas_handle);
        cusolverDnCreate(&solver_handle);
        initialized = true;
    }

    auto S_stride = batch_size > 1 ? std::distance(S_ptrs[0], S_ptrs[1]) : m*SN;
    auto STS_stride = batch_size > 1 ? std::distance(STS_ptrs[0], STS_ptrs[1]) : SN*SN;

    constexpr T alpha = 1.0;
    constexpr T beta = 0.0;

    //Compute StS = S^T * S
    cublasSgemmStridedBatched(blas_handle, CUBLAS_OP_T, CUBLAS_OP_N,
                                SN, /* m */
                                SN, /* n */
                                m, /* k */
                                &alpha, /* alpha */
                                S_ptrs.front(), /* A */
                                m, /* lda */
                                S_stride, /* strideA */
                                S_ptrs.front(), /* B */
                                m, /* ldb */
                                S_stride, /* strideB */
                                &beta, /* beta */
                                STS_ptrs.front(), /* C */
                                SN, /* ldc */
                                STS_stride, /* strideC */
                                batch_size);

    ctx -> submit([&](sycl::handler& h){
        h.parallel_for(nd_range<1>(sycl::range{size_t(batch_size*SN)}, sycl::range{size_t(SN)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            auto bid = item.get_group_linear_id();
            auto cta = item.get_group();
            auto STS_acc = STS_ptrs[bid];
            
            auto g_norm = sycl::reduce_over_group(cta, sycl::sqrt(STS_acc[tid*SN + tid]), sycl::maximum<T>());
            auto eps = std::numeric_limits<T>::epsilon();
            auto shift = T(11.0) * (m*SN*eps +(SN + 1)*SN*eps) * g_norm;
            STS_acc[tid*SN + tid] += shift;
        });
    });
    ctx.wait();

    //Compute the Cholesky Factorization of StS
    cusolverDnSpotrfBatched(solver_handle, CUBLAS_FILL_MODE_LOWER, SN, STS_ptrs.data(), SN, d_info.data(), batch_size);
    //cublasSgetrfBatched(blas_handle, SN, Sts_ptrs.data(), SN, NULL, d_info.data(), batch_size);
    cudaDeviceSynchronize();
    //Compute Q = S * StS^-1
    cublasStrsmBatched( blas_handle,
                        CUBLAS_SIDE_RIGHT, 
                        CUBLAS_FILL_MODE_LOWER, 
                        CUBLAS_OP_T, 
                        CUBLAS_DIAG_NON_UNIT, 
                        m, 
                        SN, 
                        &alpha, 
                        STS_ptrs.data(), 
                        SN, 
                        S_ptrs.data(), 
                        m, batch_size);
    cudaDeviceSynchronize();

    Chol2QR<T,SN>(ctx, S_ptrs, STS_ptrs, m, batch_size);
}

template <typename T, int ExtDim, int CandDim>
void Ortho(SyclQueue& ctx, Span<T*> S_ptrs /* Candidate Basis [k,CandDim] */, 
                            Span<T*> E_ptrs /* External Basis [k,ExtDim] */,
                            Span<T*> StS_ptrs /* Subspace Gram Matrix [k, CandDim] */,
                            Span<T> U /* Workspace [ExtDim * CandDim] */,
                            size_t k,
                            size_t ExtStride,
                            size_t CandStride,
                            size_t batch_size){
    for(int i = 0; i < 2; i++){
        ExternalOrthogonalization<T, ExtDim, CandDim>(ctx, Span<T>(S_ptrs.front(),CandStride*batch_size), Span<T>(E_ptrs.front(),ExtStride*batch_size), U, k, ExtStride, CandStride, batch_size);
        CholQR<T, CandDim>(ctx, S_ptrs, StS_ptrs, k, batch_size);
    }
}
                            

template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG_V7(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest, Span<T> out_eigvects, Span<T> out_eigvals){
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

    cusolverDnHandle_t handle{};
    cudaStream_t stream{};
    cudaDeviceProp props{};
    cusolverDnParams_t params{};
    cusolverDnCreate(&handle);
    cudaStreamCreate(&stream);
    cusolverDnSetStream(handle, stream);
    cusolverDnCreateParams(&params);
    /* LOBPCG_Helper<T, K> helper(stream, BlockVectors,
                                m, NZ, batch_size, A.data(), cols.data(), csr_row_offsets.data(),
                                S.data(), StAS.data(), S.data() + BlockVectors*9*m, lambdas.data()); */
     
    #if (SOLVER_BACKEND == CUSOLVER_BACKEND)
        T alpha = 1.0;
        T beta = 0.0;
        cusolverStatus_t status = cusolverDnXsyevBatched_bufferSize(handle, 
                                params, 
                                CUSOLVER_EIG_MODE_VECTOR, 
                                CUBLAS_FILL_MODE_LOWER,
                                m,
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                A.data(),
                                m,
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                lambdas.data(),
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                &syevd_scratchpad_size,
                                &syevd_scratchpad_host_size,
                                batch_size);

        cusparseHandle_t SpMM_AS_handle;
        CHECK_CUSPARSE(cusparseCreate(&SpMM_AS_handle))
        CHECK_CUSPARSE(cusparseSetStream(SpMM_AS_handle, stream))
        cusparseSpMatDescr_t descrA;    //Left matrix
        cusparseDnMatDescr_t descrX;    //Right matrix
        cusparseDnMatDescr_t descrR;    //Block Residual matrix
        cusparseDnMatDescr_t descrP;    //Block Previous matrix
        cusparseDnMatDescr_t descrS;    //The full subspace matrix [X, R, P]
        cusparseDnMatDescr_t descrAX;  //The result of the matrix matrix multiplication
        cusparseDnMatDescr_t descrAR; //The result of the matrix matrix multiplication
        cusparseDnMatDescr_t descrAS;  //The result of the matrix matrix multiplication
        cusparseDnMatDescr_t descrAP;  //The result of the matrix matrix multiplication
        CHECK_CUSPARSE(cusparseCreateCsr(  &descrA,    m,  m, m*NZ, csr_row_offsets.data(), cols_32I.data(), A.data(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrX,    m,  BlockVectors,   m,  S.data(),   CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrR,    m,  BlockVectors,   m,  S.data() +  2*BlockVectors*m,   CUDA_R_32F, CUSPARSE_ORDER_COL)) 
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrP,    m,  BlockVectors,   m,  S.data() +  1*BlockVectors*m,   CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrS,    m,  BlockVectors*3, m,  S.data(),   CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrAX,   m,  BlockVectors,   m,  AS.data(),  CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrAR,   m,  BlockVectors,   m,  AS.data() + 2*BlockVectors*m,  CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrAP,   m,  BlockVectors,   m,  AS.data() + 1*BlockVectors*m,  CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrAS,   m,  BlockVectors*3, m,  AS.data(),   CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCsrSetStridedBatch(     descrA, batch_size, m*NZ, m*NZ))
        auto SubspaceStride = m*BlockVectors*3; //[X_0, R_0, P_0, X_1, R_1, P_1, ...]

        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrX, batch_size,     SubspaceStride))
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrR, batch_size,     SubspaceStride))
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrS, batch_size,     SubspaceStride))
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrP, batch_size,     SubspaceStride))
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrAX, batch_size,    SubspaceStride)) //AX = A * X
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrAR, batch_size,    SubspaceStride)) //AR = A * R
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrAP, batch_size,    SubspaceStride)) //AP = A * P
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrAS, batch_size,    SubspaceStride)) //AS = A * [X, R, P]

        size_t SpMM_buffer_size = 0;
        CHECK_CUSPARSE(cusparseSpMM_bufferSize(SpMM_AS_handle,
                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha,
                                descrA,
                                descrR,
                                &beta,
                                descrAR,
                                CUDA_R_32F,
                                CUSPARSE_SPMM_CSR_ALG1,
                                &SpMM_buffer_size))
        
        
        SyclVector<T> SpMM_buffer(SpMM_buffer_size);

        CHECK_CUSPARSE(cusparseSpMM_preprocess(SpMM_AS_handle,
                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha,
                                descrA,
                                descrR,
                                &beta,
                                descrAR,
                                CUDA_R_32F,
                                CUSPARSE_SPMM_CSR_ALG1,
                                SpMM_buffer.data()))

        cublasHandle_t GEMM_STAS_handle;
        cublasCreate(&GEMM_STAS_handle);
        cublasSetStream(GEMM_STAS_handle, stream);
    

    #elif defined(USE_ROCSOLVER)
        status_t status = rocsolver_ssyevd_bufferSize(handle, params, EIG_SOLVE_MODE, FILL_MODE_LOWER, m, A.data(), m, lambdas.data(), &syevd_scratchpad_size, &syevd_scratchpad_host_size, batch_size);
    #endif

    SyclVector<T> syevd_scratchpad (syevd_scratchpad_size);
    SyclVector<T> syevd_scratchpad_host (syevd_scratchpad_host_size);

    #if (SOLVER_BACKEND == CUSOLVER_BACKEND)
        auto compute_eigenpairs = [&](T* A, T* lambdas, int m, int batch_size){
            cusolverStatus_t status = cusolverDnXsyevBatched(handle, 
                                params, 
                                CUSOLVER_EIG_MODE_VECTOR, 
                                CUBLAS_FILL_MODE_LOWER,
                                m,
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                A,
                                m,
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                lambdas,
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                syevd_scratchpad.data(),
                                syevd_scratchpad_size,
                                syevd_scratchpad_host.data(),
                                syevd_scratchpad_host_size,
                                syevd_info.data(),
                                batch_size);
            if(status != CUSOLVER_STATUS_SUCCESS){
                std::cerr << "Error in cusolverDnXsyevBatched: " << status << std::endl;
            }
        };
    #elif defined(USE_ROCSOLVER)
        auto compute_eigenpairs = [&](T* A, T* lambdas, T* eigenvectors, int m, int batch_size){
            status_t status = rocsolver_ssyevd(handle, params, EIG_SOLVE_MODE, FILL_MODE_LOWER, m, A, m, lambdas, eigenvectors, m, syevd_scratchpad, syevd_scratchpad_size, syevd_scratchpad_host_size, syevd_info.data(), batch_size);
        };
    #endif
    const T tol = sycl::sqrt(std::numeric_limits<T>::epsilon()) * 1e-2 * T(m);
    
    SyclVector<T> Diag(3*BlockVectors * batch_size, 0);
    SyclVector<T> R_inv(batch_size * 3 * BlockVectors * 3 * BlockVectors, 0);
    SyclVector<T> U(batch_size * BlockVectors * 3 * m, 0);
    SyclVector<T> C_p(batch_size * BlockVectors * BlockVectors * 3, 0); //At most it will contain BlockVectors eigenvectors of size 3*BlockVectors
    SyclVector<T> QRworkspace(batch_size * BlockVectors * 3 * m, 0);
    SyclVector<T*> S_ptrs(batch_size);
    SyclVector<T*> R_ptrs(batch_size);
    SyclVector<T*> STS_ptrs(batch_size);
    SyclVector<T*> C_p_ptrs(batch_size);
    SyclVector<T*> C_x_restart_ptrs(batch_size);
    SyclVector<T*> R_inv_ptrs(batch_size);

    for(int i = 0; i < batch_size; i++){
        S_ptrs[i] = S.data() + i * SubspaceStride;
        R_ptrs[i] = S.data() + i * (BlockVectors*3) * m + 2 * BlockVectors*m;
        STS_ptrs[i] = StAS.data() + i * BlockVectors * 3 * BlockVectors * 3;
        C_p_ptrs[i] = C_p.data() + i * BlockVectors * BlockVectors * 3;
        C_x_restart_ptrs[i] = StAS.data() + i * BlockVectors * 2 * BlockVectors * 2;
        R_inv_ptrs[i] = R_inv.data() + i * 3 * BlockVectors * 3 * BlockVectors;
    }
    //Setup
    ctx -> submit([&](sycl::handler& h ){
        auto Scache = sycl::local_accessor<T, 1>(m, h);

        h.parallel_for<LOBPCG_V7_1<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            sycl::group<1> cta = item.get_group();
            auto bid = item.get_group_linear_id();
            std::bitset<BlockVectors> converged;
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
    CholQR<T, BlockVectors>(ctx, S_ptrs, STS_ptrs, m, batch_size);

    //Compute AX
    auto SpMM_status = cusparseSpMM(SpMM_AS_handle, 
                                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    &alpha,
                                    descrA,
                                    descrX,
                                    &beta,
                                    descrAX,
                                    CUDA_R_32F,
                                    CUSPARSE_SPMM_CSR_ALG1,
                                    SpMM_buffer.data());

    cudaStreamSynchronize(stream);
    
    //Compute X^T AX

    auto gemm_status = cublasSgemmStridedBatched(GEMM_STAS_handle,
                                                CUBLAS_OP_T,
                                                CUBLAS_OP_N,
                                                BlockVectors,
                                                BlockVectors,
                                                m,
                                                &alpha,
                                                S.data(), 
                                                m, m*3*BlockVectors,
                                                AS.data(),
                                                m, m*3*BlockVectors,
                                                &beta,
                                                StAS.data(),
                                                BlockVectors, 
                                                BlockVectors * BlockVectors,
                                                batch_size); 
    cudaStreamSynchronize(stream);
    //ComputeGramMatrix<T, BlockVectors>(ctx, S, StAS, m, BlockVectors*3*m, batch_size);
    compute_eigenpairs(StAS.data(), lambdas.data(), BlockVectors, batch_size);
    cudaStreamSynchronize(stream);

    ctx -> submit([&](sycl::handler& h){
        auto Scache = sycl::local_accessor<T, 1>(m, h);

        h.parallel_for<LOBPCG_V7_2<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
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


    auto Tresiduals = std::chrono::duration<double>(0);
    auto Tortho = std::chrono::duration<double>(0);
    auto Tgemm = std::chrono::duration<double>(0);
    auto Tspmm = std::chrono::duration<double>(0);
    auto Teigen = std::chrono::duration<double>(0);
    auto Tupdate = std::chrono::duration<double>(0);
    auto Tmemcpy = std::chrono::duration<double>(0);

    auto iteration = [&](SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t iter, bool restart){
        
        auto start = std::chrono::high_resolution_clock::now();
        //Setup
        ctx -> submit([&](sycl::handler& h){
            auto Scache = sycl::local_accessor<T, 1>(m, h);

            h.parallel_for<LOBPCG_V7_3<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
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

        /* std::cout << "Residuals: " << residuals << std::endl;
        std::cout << "Converged: ";
        for(int i = 0; i < BlockVectors; i++) std::cout << converged_acc[0][i] << " ";
        std::cout << std::endl; */

        auto end = std::chrono::high_resolution_clock::now();
        Tresiduals += end - start;

        
        start = std::chrono::high_resolution_clock::now();
        if (restart){
            Ortho<T, BlockVectors, BlockVectors>(ctx, R_ptrs, S_ptrs, STS_ptrs, U, m, SubspaceStride, SubspaceStride, batch_size);
        } else {
            Ortho<T, 2*BlockVectors, BlockVectors>(ctx, R_ptrs, S_ptrs, STS_ptrs, U, m, SubspaceStride, SubspaceStride, batch_size);
        }
        end = std::chrono::high_resolution_clock::now();
        Tortho += end - start;
        /* ComputeGramMatrix<T, BlockVectors*3>(ctx, S_acc, StAS_acc, m, SubspaceStride, batch_size);
        for (int i = 0; i < batch_size; i++ ){
            std::cout << "StS matrix " << i << std::endl;
            print_matrix(StAS.subspan(i * BlockVectors * 3 * BlockVectors * 3, BlockVectors * 3 * BlockVectors * 3), 3*BlockVectors, 3*BlockVectors, false);
            
        } */

        start = std::chrono::high_resolution_clock::now();
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
        end = std::chrono::high_resolution_clock::now();
        Tmemcpy += end - start;
        //ComputeGramMatrix<T,BlockVectors*3>(ctx, S, R_inv, m, BlockVectors*m*3, batch_size);
        //print_matrix(R_inv.to_span(), BlockVectors*3 , BlockVectors*3 , false);

        start = std::chrono::high_resolution_clock::now();
        //Compute AR
        auto SpMM_status = cusparseSpMM(SpMM_AS_handle, 
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        &alpha,
                                        descrA,
                                        restart ? descrP : descrR,
                                        &beta,
                                        restart ? descrAP: descrAR,
                                        CUDA_R_32F,
                                        CUSPARSE_SPMM_CSR_ALG1,
                                        SpMM_buffer.data());

        cudaStreamSynchronize(stream);
        end = std::chrono::high_resolution_clock::now();
        Tspmm += end - start;        

        start = std::chrono::high_resolution_clock::now();
        //Compute X^T AX
        auto gemm_status = cublasSgemmStridedBatched(GEMM_STAS_handle,
                                                    CUBLAS_OP_T,
                                                    CUBLAS_OP_N,
                                                    BlockVectors * (restart ? 2 : 3),
                                                    BlockVectors * (restart ? 2 : 3),
                                                    m,
                                                    &alpha,
                                                    S.data(), 
                                                    m, m*3*BlockVectors,
                                                    AS.data(),
                                                    m, m*3*BlockVectors,
                                                    &beta,
                                                    StAS.data(),
                                                    BlockVectors * (restart ? 2: 3), 
                                                    BlockVectors * (restart ? 2*2: 3*3) * BlockVectors,
                                                    batch_size); 

        //Solve the eigenvalue problem
        cudaStreamSynchronize(stream);
        end = std::chrono::high_resolution_clock::now();
        Tgemm += end - start;

        start = std::chrono::high_resolution_clock::now();
        compute_eigenpairs(StAS.data(), lambdas.data(), restart ? BlockVectors*2 : BlockVectors*3, batch_size);
        cudaStreamSynchronize(stream);
        end = std::chrono::high_resolution_clock::now();
        Teigen += end - start;

        start = std::chrono::high_resolution_clock::now();
        cudaMemcpy(R_inv.data(), StAS.data(), StAS.size() * sizeof(T), cudaMemcpyDeviceToDevice);
        cudaDeviceSynchronize();

        //If largest = true, then the order of the eigenvectors is reversed
        if (largest){
            ctx -> submit([&](sycl::handler& h){
                auto cols = restart ? BlockVectors*2 : BlockVectors*3;
                auto R_inv_acc = R_inv.to_span();
                h.parallel_for(nd_range<1>(sycl::range{size_t(batch_size*256)}, sycl::range{size_t(256)}), [=](sycl::nd_item<1> item){
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
        end = std::chrono::high_resolution_clock::now();
        Tmemcpy += end - start;
        start = std::chrono::high_resolution_clock::now();
        //Compute X = [X, R, P] * [Zx, Zr, Zp]^T
        
        cublasSgemmStridedBatched(GEMM_STAS_handle,
                                CUBLAS_OP_N,
                                CUBLAS_OP_N,
                                m,
                                BlockVectors,
                                BlockVectors * (restart ? 2 : 3),
                                &alpha,
                                S.data(),
                                m, m*3*BlockVectors,
                                StAS.data(),
                                BlockVectors * (restart ? 2 : 3),
                                BlockVectors * (restart ? 2*2 : 3*3) * BlockVectors,
                                &beta,
                                U.data(),
                                m,
                                m*3*BlockVectors,
                                batch_size);

        //Make an implicit update of AX:
        cublasSgemmStridedBatched(  GEMM_STAS_handle,
                                    CUBLAS_OP_N,
                                    CUBLAS_OP_N,
                                    m,
                                    BlockVectors,
                                    BlockVectors * (restart ? 2 : 3),
                                    &alpha,
                                    AS.data(),
                                    m, m*3*BlockVectors,
                                    StAS.data(),
                                    BlockVectors * (restart ? 2 : 3),
                                    BlockVectors * (restart ? 2*2 : 3*3) * BlockVectors,
                                    &beta,
                                    S_temp.data(),
                                    m,
                                    m*3*BlockVectors,
                                    batch_size);

        cudaStreamSynchronize(stream);
        end = std::chrono::high_resolution_clock::now();
        Tgemm += end - start;

        start = std::chrono::high_resolution_clock::now();
        ExternalOrthogonalization<T, BlockVectors, BlockVectors>(ctx, C_p, StAS, R_inv, restart ? BlockVectors*2 : BlockVectors*3, BlockVectors*BlockVectors * (restart? (2*2) : (3*3)), BlockVectors*BlockVectors *  3, batch_size);
        CholQR<T, BlockVectors>(ctx, C_p_ptrs, STS_ptrs, BlockVectors * (restart ? 2 : 3), batch_size); //Since we're done using StAS (eigenvectors of previous iteration), we can use it as a scratchpad
        //ExternalOrthogonalization<T, BlockVectors, BlockVectors>(ctx, C_p, StAS, R_inv, restart ? BlockVectors*2 : BlockVectors*3, BlockVectors*BlockVectors * (restart? (2*2) : (3*3)), BlockVectors*BlockVectors *  3, batch_size);
        //CholQR<T, BlockVectors>(ctx, C_p_ptrs, STS_ptrs, BlockVectors * (restart ? 2 : 3), batch_size); //Since we're done using StAS (eigenvectors of previous iteration), we can use it as a scratchpad
        //Ortho<T, BlockVectors, BlockVectors>(ctx, C_p_ptrs, restart ? C_x_restart_ptrs : STS_ptrs, R_inv_ptrs, QRworkspace, BlockVectors * (restart ? 2 : 3), BlockVectors*BlockVectors * (restart? (2*2) : (3*3)), BlockVectors*BlockVectors *  3, batch_size);


        end = std::chrono::high_resolution_clock::now();
        Tortho += end - start;

        start = std::chrono::high_resolution_clock::now();
        end = std::chrono::high_resolution_clock::now();
        Tmemcpy += end - start;
        //First compute P = [X, R, P] * C_p^T
        start = std::chrono::high_resolution_clock::now();
        cublasSgemmStridedBatched(GEMM_STAS_handle,
                                CUBLAS_OP_N,
                                CUBLAS_OP_N,
                                m, /* m */
                                BlockVectors, /* n */
                                BlockVectors * (restart ? 2 : 3), /* k */
                                &alpha, /* alpha */
                                S.data(),    /* A */
                                m, m*3*BlockVectors, /* lda, strideA */
                                C_p.data(), /* B */
                                BlockVectors * (restart ? 2 : 3), /* ldb */
                                BlockVectors * 3 * BlockVectors, /* strideB */
                                &beta, /* beta */
                                U.data() + BlockVectors*1*m, /* C */
                                m, /* ldc */
                                m*3*BlockVectors, /* strideC */
                                batch_size);
        
        //Make an implicit update of AP:
        cublasSgemmStridedBatched(  GEMM_STAS_handle,
                                    CUBLAS_OP_N,
                                    CUBLAS_OP_N,
                                    m,
                                    BlockVectors,
                                    BlockVectors * (restart ? 2 : 3),
                                    &alpha,
                                    AS.data(),
                                    m, m*3*BlockVectors,
                                    C_p.data(),
                                    BlockVectors * (restart ? 2 : 3),
                                    BlockVectors * 3 * BlockVectors,
                                    &beta,
                                    S_temp.data() + m*BlockVectors*1,
                                    m,
                                    m*3*BlockVectors,
                                    batch_size);

        cudaStreamSynchronize(stream);
        end = std::chrono::high_resolution_clock::now();
        Tgemm += end - start;

        start = std::chrono::high_resolution_clock::now();
        cudaMemcpy(AS.data(), S_temp.data(), S_temp.size() * sizeof(T), cudaMemcpyDeviceToDevice);
        cudaMemcpy(S.data(), U.data(), U.size() * sizeof(T), cudaMemcpyDeviceToDevice);
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        Tmemcpy += end - start;
        //std::cout << "Post Update Gram Matrix" << std::endl;
        //ComputeGramMatrix<T, BlockVectors*3>(ctx, S, StAS, m, BlockVectors*3*m, batch_size);
        //print_matrix(StAS.to_span(), BlockVectors * 3, BlockVectors * 3, false);
    };
    auto T0 = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < maxiters; i++){
        iteration(ctx, A, cols, batch_size, m, i, i < 1 ? true : false);
    }
    ctx.wait();   
    auto T1 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(T1 - T0).count();
    std::cout << "LOBPCG V4 Total Iteration Time: " << duration / batch_size << " µs / fullerene" << std::endl;
    std::cout << "Residuals Time: " << std::chrono::duration_cast<std::chrono::microseconds>(Tresiduals).count() / batch_size << " µs / fullerene" << std::endl;
    std::cout << "Orthogonalization Time: " << std::chrono::duration_cast<std::chrono::microseconds>(Tortho).count() / batch_size << " µs / fullerene" << std::endl;
    std::cout << "GEMM Time: " << std::chrono::duration_cast<std::chrono::microseconds>(Tgemm).count() / batch_size << " µs / fullerene" << std::endl;
    std::cout << "SpMM Time: " << std::chrono::duration_cast<std::chrono::microseconds>(Tspmm).count() / batch_size << " µs / fullerene" << std::endl;
    std::cout << "Eigenvalue Time: " << std::chrono::duration_cast<std::chrono::microseconds>(Teigen).count() / batch_size << " µs / fullerene" << std::endl;
    std::cout << "Update Time: " << std::chrono::duration_cast<std::chrono::microseconds>(Tupdate).count() / batch_size << " µs / fullerene" << std::endl;
    std::cout << "Memcpy Time: " << std::chrono::duration_cast<std::chrono::microseconds>(Tmemcpy).count() / batch_size << " µs / fullerene" << std::endl;


    //std::cout << "Eigenvalues: " << lambdas << std::endl;
    //std::cout << "Residuals: " << residuals << std::endl;

    for(int i = 0; i < batch_size; i++){
        primitives::copy(ctx, lambdas.subspan(i * (BlockVectors *3) + BlockVectors * 2, BlockVectors), out_eigvals.subspan(i * BlockVectors, BlockVectors));
        //primitives::copy(ctx, S_acc.subspan(i * BlockVectors * 3 * m, BlockVectors * m), eigvects.subspan(i * BlockVectors * m, BlockVectors * m));
    }
}

template void LOBPCG_V7<float, uint16_t, 3, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest, Span<float> eigvects, Span<float> eigvals);
template void LOBPCG_V7<float, uint16_t, 21, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest, Span<float> eigvects, Span<float> eigvals);