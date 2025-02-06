#include <fullerenes/linalg/qr/batched_qr.hh>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include "../../../primitives.cc"
#include "../../../../linalg/linalg-impl.hh"

using namespace linalg;

//Non-owning memory pool
struct Mempool {
    void* data;
    size_t byte_size;

    Mempool(void* data, size_t byte_size): data(data), byte_size(byte_size){}

    template<typename T>
    static size_t allocation_size(size_t size) {
        std::uintptr_t alignment_padding = alignof(T) - 1;
        return (size * sizeof(T)) + alignment_padding;
    }

    template<typename T>
    T* allocate(size_t size){
        size_t alloc_size = allocation_size<T>(size);
        if (alloc_size > byte_size){
            throw std::runtime_error("Out of memory");
        }
        std::uintptr_t addr = reinterpret_cast<std::uintptr_t>(data);
        std::uintptr_t aligned = (addr + alignof(T) - 1) & ~(alignof(T) - 1);
        T* ptr = reinterpret_cast<T*>(aligned);

        data = reinterpret_cast<void*>(ptr + size);
        byte_size -= (reinterpret_cast<char*>(data) - reinterpret_cast<char*>(ptr));

        return ptr;
    }
};


template <typename T>
void chol_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace) 
{
    Mempool pool(workspace.data(), workspace.size());
    //static SyclVector<int> d_info(batch_size);    
    auto d_info =   Span(pool.allocate<int>(batch_size), batch_size);
    auto A_ptrs =   Span(pool.allocate<T*>(batch_size), batch_size);
    auto ATA_ptrs = Span(pool.allocate<T*>(batch_size), batch_size);
    auto ATA =      Span(pool.allocate<T>(n*n*batch_size), n*n*batch_size);

    auto ATA_stride = n*n;
    
    ctx -> parallel_for(batch_size, [=](sycl::id<1> id){
        A_ptrs[id] = A.data() + id*Astride;
        ATA_ptrs[id] = ATA.data() + id*ATA_stride;
    });
    ctx.wait();

    static bool initialized = false;
    static LinalgHandle<Backend::CUDA> handle;
    handle.setStream(ctx);

    std::cout << "A size: " << A.size() << std::endl;
    std::cout << "Workspace size: " << workspace.size() << std::endl;


/*     cublasHandle_t handle; 
    auto blas_status = cublasCreate(&handle);
    if (blas_status != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "cuBLAS initialization failed (in chol_qr.cc) with status: " << blas_status << std::endl;
        throw std::runtime_error("cuBLAS initialization failed");
    }
 */    //cublasSetStream(handle, sycl::get_native<sycl::backend::ext_oneapi_cuda>(*ctx));

    /* cusolverDnHandle_t solver_handle; 
    auto solver_status = cusolverDnCreate(&solver_handle);
    if (solver_status != CUSOLVER_STATUS_SUCCESS) {
        std::cerr << "cusolver initialization failed (in chol_qr.cc) with status: " << solver_status << std::endl;
        throw std::runtime_error("cusolver initialization failed");
    } */
    //cusolverDnSetStream(solver_handle, sycl::get_native<sycl::backend::ext_oneapi_cuda>(*ctx));

    constexpr T alpha = 1.0;
    constexpr T beta = 0.0;
    auto descrA = MatHandle<T, BatchType::Batched>(A_ptrs.front(), m, n, m, Astride, batch_size);
    auto descrC = MatHandle<T, BatchType::Batched>(ATA_ptrs.front(), n, n, n, ATA_stride, batch_size);

    //Compute StS = S^T * S
    gemm<Backend::CUDA>(ctx, descrA, descrA, descrC, alpha, beta, Transpose::Trans, Transpose::NoTrans);
    

    //cublasSgemmStridedBatched(blas_handle, CUBLAS_OP_T, CUBLAS_OP_N,
    //                            n, /* m */
    //                            n, /* n */
    //                            m, /* k */
    //                            &alpha, /* alpha */
    //                            A_ptrs.front(), /* A */
    //                            m, /* lda */
    //                            Astride, /* strideA */
    //                            A_ptrs.front(), /* B */
    //                            m, /* ldb */
    //                            Astride, /* strideB */
    //                            &beta, /* beta */
    //                            ATA_ptrs.front(), /* C */
    //                            n, /* ldc */
    //                            ATA_stride, /* strideC */
    //                            batch_size);

    //Compute the Cholesky Factorization of StS
    cusolverDnSpotrfBatched(handle, CUBLAS_FILL_MODE_LOWER, n, ATA_ptrs.data(), n, d_info.data(), batch_size);
    //cublasSgetrfBatched(blas_handle, SN, ATA_ptrs.data(), SN, NULL, d_info.data(), batch_size);

    //Compute Q = S * StS^-1
    cublasStrsmBatched( handle,
                        CUBLAS_SIDE_RIGHT, 
                        CUBLAS_FILL_MODE_LOWER, 
                        CUBLAS_OP_T, 
                        CUBLAS_DIAG_NON_UNIT, 
                        m, 
                        n, 
                        &alpha, 
                        ATA_ptrs.data(), 
                        n, 
                        A_ptrs.data(), 
                        m, batch_size);

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
    return Mempool::allocation_size<T>(batch_size) + 2*Mempool::allocation_size<T*>(batch_size) + Mempool::allocation_size<T>(n*n*batch_size);
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