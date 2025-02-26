#include "../../../linalg-impl.hh"
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <sycl/sycl.hpp>

//Signature adapter
namespace linalg {

    template <Backend B, typename T, BatchType BT>
    SyclEvent gemm(SyclQueue& ctx,
                   DenseMatHandle<T,BT>& descrA,
                   DenseMatHandle<T,BT>& descrB,
                   DenseMatHandle<T,BT>& descrC,
                   T alpha,
                   T beta,
                   Transpose transA,
                   Transpose transB,
                   ComputePrecision precision) {
        // Call cuBLAS
        static LinalgHandle<B> handle;
        handle.setStream(ctx);
        
        auto [m, k] = get_effective_dims(descrA, transA);
        auto [kB, n] = get_effective_dims(descrB, transB);

        if constexpr (BT == BatchType::Single) {
            //Can't really use the call_backend function here, because cublasGemmEx is an overloaded function
            cublasGemmEx(handle,
                enum_convert<BackendLibrary::CUBLAS>(transA), enum_convert<BackendLibrary::CUBLAS>(transB),
                m, n, k,
                &alpha,
                descrA.data_, BackendScalar<T,B>::type, descrA.ld_,
                descrB.data_, BackendScalar<T,B>::type, descrB.ld_,
                &beta,
                descrC.data_, BackendScalar<T,B>::type, descrC.ld_,
                enum_convert<BackendLibrary::CUBLAS, T>(precision),
                CUBLAS_GEMM_DFALT);
        } else {
            //Can't really use the call_backend function here, because cublasGemmStridedBatchedEx is an overloaded function
            cublasGemmStridedBatchedEx(handle,
                enum_convert<BackendLibrary::CUBLAS>(transA), enum_convert<BackendLibrary::CUBLAS>(transB),
                m, n, k,
                &alpha,
                descrA.data_, BackendScalar<T,B>::type, descrA.ld_, descrA.stride_,
                descrB.data_, BackendScalar<T,B>::type, descrB.ld_, descrB.stride_,
                &beta,
                descrC.data_, BackendScalar<T,B>::type, descrC.ld_, descrC.stride_,
                descrA.batch_size_,
                enum_convert<BackendLibrary::CUBLAS, T>(precision),
                CUBLAS_GEMM_DFALT);
        }
        return ctx.get_event();
    }

    template <Backend B, typename T, BatchType BT>
    SyclEvent trsm(SyclQueue& ctx,
                   DenseMatHandle<T,BT>& descrA,
                   DenseMatHandle<T,BT>& descrB,
                   Side side,
                   Uplo uplo,
                   Transpose transA,
                   Diag diag,
                   T alpha) {
        static LinalgHandle<B> handle;
        handle.setStream(ctx);
        auto [kB, n] = get_effective_dims(descrB, Transpose::NoTrans);

        if constexpr (BT == BatchType::Single) {
            call_backend<T, BackendLibrary::CUBLAS, B>(cublasStrsm, cublasDtrsm, cublasCtrsm, cublasZtrsm, 
                handle, side, uplo, transA, diag, kB, n, &alpha, get_data(descrA), descrA.ld_, get_data(descrB), descrB.ld_); 
        } else {
            call_backend<T, BackendLibrary::CUBLAS, B>(cublasStrsmBatched, cublasDtrsmBatched, cublasCtrsmBatched, cublasZtrsmBatched, 
                handle, side, uplo, transA, diag, kB, n, &alpha, get_data(descrA), descrA.ld_, get_data(descrB), descrB.ld_, descrA.batch_size_);
        }
        return ctx.get_event();
    }

    template <Backend B, typename T, BatchType BT>
    size_t potrf_buffer_size(SyclQueue& ctx,
                                DenseMatHandle<T,BT>& A,
                                Uplo uplo) {
        static LinalgHandle<B> handle;
        handle.setStream(ctx);
        int size = 0;
        if constexpr (BT == BatchType::Single) {
            call_backend<T, BackendLibrary::CUSOLVER, B>(cusolverDnSpotrf_bufferSize, cusolverDnDpotrf_bufferSize, cusolverDnCpotrf_bufferSize, cusolverDnZpotrf_bufferSize,
                handle, uplo, A.rows_, get_data(A), A.ld_, &size);
            size = BumpAllocator::allocation_size<std::byte>(ctx, size);
        } else {
            size =  BumpAllocator::allocation_size<int>(ctx, A.batch_size_);
        }
        return size;
    }

    template <Backend B, typename T, BatchType BT>
    SyclEvent potrf(SyclQueue& ctx,
                    DenseMatHandle<T,BT>& descrA,
                    Uplo uplo,
                    Span<std::byte> workspace) {        
        static LinalgHandle<B> handle;
        handle.setStream(ctx);
        BumpAllocator pool(workspace);
        auto Lwork = potrf_buffer_size<B>(ctx, descrA, uplo);
        
        if constexpr (BT == BatchType::Single) {
            int info;
            auto potrf_span = pool.allocate<std::byte>(ctx, Lwork);
            call_backend<T, BackendLibrary::CUSOLVER, B>(cusolverDnSpotrf, cusolverDnDpotrf, cusolverDnCpotrf, cusolverDnZpotrf,
                handle, uplo, descrA.rows_, get_data(descrA), descrA.ld_, reinterpret_cast<T*>(potrf_span.data()), Lwork, &info);
        } else {
            auto info = pool.allocate<int>(ctx, descrA.batch_size_);
            call_backend<T, BackendLibrary::CUSOLVER, B>(cusolverDnSpotrfBatched, cusolverDnDpotrfBatched, cusolverDnCpotrfBatched, cusolverDnZpotrfBatched,
                handle, uplo, descrA.rows_, get_data(descrA), descrA.ld_, info.data(), descrA.batch_size_);
        }
        return ctx.get_event();
        
    }
        }
        return ctx.get_event();
    }

    


                   
    #define GEMM_INSTANTIATE(fp, BT) \
    template SyclEvent gemm<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatHandle<fp, BT>&, \
        DenseMatHandle<fp, BT>&, \
        DenseMatHandle<fp, BT>&, \
        fp, fp, Transpose, Transpose, ComputePrecision);

    #define TRSM_INSTANTIATE(fp, BT) \
    template SyclEvent trsm<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatHandle<fp, BT>&, \
        DenseMatHandle<fp, BT>&, \
        Side, Uplo, Transpose, Diag, fp);

    #define POTRF_INSTANTIATE(fp, BT) \
    template SyclEvent potrf<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatHandle<fp, BT>&, \
        Uplo, \
        Span<std::byte>);
    
    #define POTRF_BUFFER_SIZE_INSTANTIATE(fp, BT) \
    template size_t potrf_buffer_size<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatHandle<fp, BT>&, \
        Uplo);
    
    #define BLAS_INSTANTIATE(fp, BT) \
        GEMM_INSTANTIATE(fp, BT) \
        TRSM_INSTANTIATE(fp, BT) \
        POTRF_INSTANTIATE(fp, BT) \
        POTRF_BUFFER_SIZE_INSTANTIATE(fp, BT)
    
    // Macro that covers all layout and batch type combinations for a given floating-point type.
    #define GEMM_INSTANTIATE_FOR_FP(fp)                \
        BLAS_INSTANTIATE(fp, BatchType::Batched)        \
        BLAS_INSTANTIATE(fp, BatchType::Single)        

    // Instantiate for the floating-point types of interest.
    GEMM_INSTANTIATE_FOR_FP(float)
    GEMM_INSTANTIATE_FOR_FP(double)
    GEMM_INSTANTIATE_FOR_FP(std::complex<float>)
    GEMM_INSTANTIATE_FOR_FP(std::complex<double>)

    #undef GEMM_INSTANTIATE
    #undef GEMM_INSTANTIATE_FOR_FP
}


