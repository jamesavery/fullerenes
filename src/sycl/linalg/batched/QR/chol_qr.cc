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
    auto ATA =          pool.allocate<T>(ctx, n*n*batch_size);
    auto matAmem =      pool.allocate<T*>(ctx, batch_size);
    auto matATAmem =    pool.allocate<T*>(ctx, batch_size);
    auto ATA_stride = n*n;

    constexpr T alpha = 1.0;
    constexpr T beta = 0.0;
    auto descrA = DenseMatView<T, BatchType::Batched>(A.data(), m, n, m, Astride, batch_size, matAmem);
    auto descrC = DenseMatView<T, BatchType::Batched>(ATA.data(), n, n, n, ATA_stride, batch_size, matATAmem);
    
    auto potrf_workspace = pool.allocate<std::byte>(ctx, potrf_buffer_size<Backend::CUDA>(ctx, descrC, Uplo::Lower));
    //Compute StS = S^T * S
    gemm<Backend::CUDA>(ctx, descrA, descrA, descrC, alpha, beta, Transpose::Trans, Transpose::NoTrans);
    //Compute the Cholesky Factorization of StS
    potrf<Backend::CUDA>(ctx, descrC, Uplo::Lower, potrf_workspace);
    //Compute Q = S * StS^-1 (S is overwritten with Q)
    trsm<Backend::CUDA>(ctx, descrC, descrA, Side::Right, Uplo::Lower, Transpose::Trans, Diag::NonUnit, alpha);
    //Compute the QR factorization of Q
    gemm<Backend::CUDA>(ctx, descrA, descrA, descrC, alpha, beta, Transpose::Trans, Transpose::NoTrans);
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
    return  BumpAllocator::allocation_size<std::byte>(ctx, potrf_buffer_size<Backend::CUDA>(ctx, DenseMatView<T, BatchType::Batched>(A.data(), m, n, m, Astride, batch_size, Span<T*>{}), Uplo::Lower)) + 
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
    
    BumpAllocator pool(workspace);
    auto ATA = pool.allocate<T>(ctx, n*n*batch_size);
    auto matAmem = pool.allocate<T*>(ctx, batch_size);
    auto matATAmem = pool.allocate<T*>(ctx, batch_size);
    
    auto ATA_stride = n*n;
    auto descrA = DenseMatView<T, BatchType::Batched>(A.data(), m, n, m, Astride, batch_size, matAmem);
    auto descrC = DenseMatView<T, BatchType::Batched>(ATA.data(), n, n, n, ATA_stride, batch_size, matATAmem);
    auto potrf_workspace = pool.allocate<std::byte>(ctx, potrf_buffer_size<Backend::CUDA>(ctx, descrC, Uplo::Lower));

    gemm<Backend::CUDA>(ctx, descrA, descrA, descrC, alpha, beta, Transpose::Trans, Transpose::NoTrans);

    auto ATA_ptrs = get_data(descrC);
    ctx -> submit([&](sycl::handler& h){
        h.parallel_for(sycl::nd_range<1>(sycl::range{size_t(batch_size*n)}, sycl::range{size_t(n)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            auto bid = item.get_group_linear_id();
            auto cta = item.get_group();
            auto ATA_acc = ATA_ptrs + bid * ATA_stride;
            
            auto g_norm = sycl::reduce_over_group(cta, sycl::sqrt(ATA_acc[tid*n + tid]), sycl::maximum<T>());
            auto eps = std::numeric_limits<T>::epsilon();
            auto shift = T(11.0) * (m*n*eps + (n + 1)*n*eps) * g_norm;
            ATA_acc[tid*n + tid] += shift;
        });
    });
    ctx.wait();

    //Compute the Cholesky Factorization of StS
    potrf<Backend::CUDA>(ctx, descrC, Uplo::Lower, potrf_workspace);
    trsm<Backend::CUDA>(ctx, descrC, descrA, Side::Right, Uplo::Lower, Transpose::Trans, Diag::NonUnit, alpha);

    chol_qr_batched(ctx, transpose_op, batch_size, Astride, m, n, A, workspace);
    chol_qr_batched(ctx, transpose_op, batch_size, Astride, m, n, A, workspace);
} 


template <typename T>
size_t shift_chol3_qr_batched_buffer_size(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A){
    return chol_qr_batched_buffer_size(ctx, transpose_op, batch_size, Astride, m, n, A) + 
            BumpAllocator::allocation_size<T>(ctx, n*n*batch_size);
}

//Instantiate the template
template void chol_qr_batched<float>(SyclQueue&,bool,int,int,int,int,Span<float>,Span<std::byte> workspace);

//Instantiate the template
template size_t chol_qr_batched_buffer_size<float>(SyclQueue&,bool,int,int,int,int,Span<float>);

//Instantiate the template
template void chol2_qr_batched<float>(SyclQueue&,bool,int,int,int,int,Span<float>,Span<std::byte> workspace);

//Instantiate the template
template size_t chol2_qr_batched_buffer_size<float>(SyclQueue&,bool,int,int,int,int,Span<float>);

//Instantiate the template
template void shift_chol3_qr_batched<float>(SyclQueue&,bool,int,int,int,int,Span<float>,Span<std::byte> workspace);

//Instantiate the template
template size_t shift_chol3_qr_batched_buffer_size<float>(SyclQueue&,bool,int,int,int,int,Span<float>);