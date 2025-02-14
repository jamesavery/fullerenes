#include <fullerenes/linalg/qr/batched_qr.hh>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/mempool.hh>
#include "../../../primitives.cc"
#include "../../../../linalg/linalg-impl.hh"

using namespace linalg;

//Non-owning memory pool

template <typename T>
void chol_qr_batched(   SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace) 
{
    BumpAllocator pool(workspace.data(), workspace.size());
    //static SyclVector<int> d_info(batch_size);    
    auto ATA =      pool.allocate<T>(ctx, n*n*batch_size);
    auto ATA_stride = n*n;

    constexpr T alpha = 1.0;
    constexpr T beta = 0.0;
    auto descrA = DenseMatHandle<T, BatchType::Batched>(A.data(), m, n, m, Astride, batch_size);
    auto descrC = DenseMatHandle<T, BatchType::Batched>(ATA.data(), n, n, n, ATA_stride, batch_size);
    
    auto potrf_workspace = pool.allocate<std::byte>(ctx, potrf_buffer_size<Backend::CUDA>(ctx, descrC, Uplo::Lower));
    //Compute StS = S^T * S
    gemm<Backend::CUDA>(ctx, descrA, descrA, descrC, alpha, beta, Transpose::Trans, Transpose::NoTrans);
    //Compute the Cholesky Factorization of StS
    potrf<Backend::CUDA>(ctx, descrC, Uplo::Lower, potrf_workspace);
    //Compute Q = S * StS^-1 (S is overwritten with Q)
    trsm<Backend::CUDA>(ctx, descrC, descrA, Side::Right, Uplo::Lower, Transpose::Trans, Diag::NonUnit, alpha);
    //Compute the QR factorization of Q
    gemm<Backend::CUDA>(ctx, descrA, descrA, descrC, alpha, beta, Transpose::Trans, Transpose::NoTrans);
    
    print_matrix(ATA, n, n);
    ctx.wait();
}

template <typename T>
size_t chol_qr_batched_buffer_size(SyclQueue& ctx,
                                bool transpose_op,
                                int batch_size,
                                int Astride,
                                int m,
                                int n,
                                Span<T> A) {
    return  BumpAllocator::allocation_size<T>(ctx, batch_size) + 
            2*BumpAllocator::allocation_size<T*>(ctx, batch_size) + 
            BumpAllocator::allocation_size<T>(ctx, n*n*batch_size);
}

template <typename T>
void chol2_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace){
    chol_qr_batched(ctx, transpose_op, batch_size, Astride, m, n, A, workspace);
    chol_qr_batched(ctx, transpose_op, batch_size, Astride, m, n, A, workspace);
}
template <typename T>
size_t chol2_qr_batched_buffer_size(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A){
    return chol_qr_batched_buffer_size(ctx, transpose_op, batch_size, Astride, m, n, A);
}

template <typename T>
void shift_chol3_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace){
    constexpr T alpha = 1.0;
    constexpr T beta = 0.0;

    //Compute StS = S^T * S
    //cublasSgemmStridedBatched(blas_handle, CUBLAS_OP_T, CUBLAS_OP_N,
    //                            SN, /* m */
    //                            SN, /* n */
    //                            m, /* k */
    //                            &alpha, /* alpha */
    //                            S_ptrs.front(), /* A */
    //                            m, /* lda */
    //                            S_stride, /* strideA */
    //                            S_ptrs.front(), /* B */
    //                            m, /* ldb */
    //                            S_stride, /* strideB */
    //                            &beta, /* beta */
    //                            STS_ptrs.front(), /* C */
    //                            SN, /* ldc */
    //                            STS_stride, /* strideC */
    //                            batch_size);

    /* ctx -> submit([&](sycl::handler& h){
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



    chol_qr_batched(ctx, transpose_op, batch_size, Astride, m, n, A, workspace);
    chol_qr_batched(ctx, transpose_op, batch_size, Astride, m, n, A, workspace); */
} 

//Instantiate the template
template void chol_qr_batched<float>(SyclQueue&,bool,int,int,int,int,Span<float>,Span<std::byte> workspace);

//Instantiate the template
template size_t chol_qr_batched_buffer_size<float>(SyclQueue&,bool,int,int,int,int,Span<float>);