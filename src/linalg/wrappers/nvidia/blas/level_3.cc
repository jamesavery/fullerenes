#include "../../../linalg-impl.hh"
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <sycl/sycl.hpp>

//Signature adapter
namespace linalg {

    template <Backend B, typename T, BatchType BT>
    SyclEvent gemm(SyclQueue& ctx,
                   DenseMatView<T,BT> descrA,
                   DenseMatView<T,BT> descrB,
                   DenseMatView<T,BT> descrC,
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
                   DenseMatView<T,BT> descrA,
                   DenseMatView<T,BT> descrB,
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
                handle, side, uplo, transA, diag, kB, n, &alpha, get_ptr_arr(ctx, descrA), descrA.ld_, get_ptr_arr(ctx, descrB), descrB.ld_, descrA.batch_size_);
        }
        return ctx.get_event();
    }

    template <Backend B, typename T, BatchType BT>
    size_t potrf_buffer_size(SyclQueue& ctx,
                                DenseMatView<T,BT> A,
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
                    DenseMatView<T,BT> descrA,
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
                handle, uplo, descrA.rows_, get_ptr_arr(ctx, descrA), descrA.ld_, info.data(), descrA.batch_size_);
        }
        return ctx.get_event();
        
    }

    template <Backend B, typename T, Format F, BatchType BT>
    SyclEvent spmm(SyclQueue& ctx,
                   SparseMatHandle<T, F, BT>& descrA,
                   DenseMatView<T,BT> descrB,
                   DenseMatView<T,BT> descrC,
                   T alpha,
                   T beta,
                   Transpose transA,
                   Transpose transB,
                   Span<std::byte> workspace) {
        // Call cuBLAS
        static LinalgHandle<B> handle;
        handle.setStream(ctx);

        BumpAllocator pool(workspace);
        auto buffer_size = spmm_buffer_size<B>(ctx, descrA, descrB, descrC, alpha, beta, transA, transB);
        auto buffer = pool.allocate<std::byte>(ctx, buffer_size);

        cusparseSpMM(
            handle,
            enum_convert<BackendLibrary::CUSPARSE>(transA),//transA == Transpose::NoTrans ? CUSPARSE_OPERATION_NON_TRANSPOSE : CUSPARSE_OPERATION_TRANSPOSE,
            enum_convert<BackendLibrary::CUSPARSE>(transB),//transB == Transpose::NoTrans ? CUSPARSE_OPERATION_NON_TRANSPOSE : CUSPARSE_OPERATION_TRANSPOSE,
            &alpha,
            *descrA,
            *descrB,
            &beta,
            *descrC,
            BackendScalar<T,B>::type,
            CUSPARSE_SPMM_ALG_DEFAULT,
            buffer.data()
        );
        return ctx.get_event();
    }

    template <Backend B, typename T, Format F, BatchType BT>
    size_t spmm_buffer_size(SyclQueue& ctx,
                               SparseMatHandle<T, F, BT>& descrA,
                               DenseMatView<T,BT> descrB,
                               DenseMatView<T,BT> descrC,
                               T alpha,
                               T beta,
                               Transpose transA,
                               Transpose transB) {
        // Call cuBLAS
        static LinalgHandle<B> handle;
        handle.setStream(ctx);

        size_t size = 0;

        cusparseSpMM_bufferSize(
            handle,
            enum_convert<BackendLibrary::CUSPARSE>(transA),
            enum_convert<BackendLibrary::CUSPARSE>(transB),
            &alpha,
            *descrA,
            *descrB,
            &beta,
            *descrC,
            BackendScalar<T,B>::type,
            CUSPARSE_SPMM_ALG_DEFAULT,
            &size
        );
        return size;
    }

    template <Backend B, typename T, BatchType BT>
    SyclEvent syev(SyclQueue& ctx,
                   DenseMatView<T,BT> descrA,
                   Span<T> eigenvalues,
                   Uplo uplo,
                   Span<std::byte> workspace) {
        static LinalgHandle<B> handle;
        handle.setStream(ctx);
        BumpAllocator pool(workspace);

        size_t l_work_device;
        size_t l_work_host;
        
        cusolverDnParams_t params;
        cusolverDnCreateParams(&params);
        if constexpr (BT == BatchType::Single){
            cusolverDnXsyevd_bufferSize(handle, params, CUSOLVER_EIG_MODE_VECTOR, enum_convert<BackendLibrary::CUSOLVER>(uplo), descrA.rows_, BackendScalar<T,B>::type, descrA.data_, descrA.ld_, BackendScalar<T,B>::type, eigenvalues.data(), BackendScalar<T,B>::type, &l_work_device, &l_work_host);
        } else {
            cusolverDnXsyevBatched_bufferSize(handle, params, CUSOLVER_EIG_MODE_VECTOR, enum_convert<BackendLibrary::CUSOLVER>(uplo), descrA.rows_, BackendScalar<T,B>::type, descrA.data_, descrA.ld_, BackendScalar<T,B>::type, eigenvalues.data(), BackendScalar<T,B>::type, &l_work_device, &l_work_host, descrA.batch_size_);
        }
        
        auto host_workspace = pool.allocate<std::byte>(ctx, l_work_host);
        auto device_workspace = pool.allocate<std::byte>(ctx, l_work_device);
                                                       //auto syev_span = pool.allocate<std::byte>(ctx, Lwork);
        if constexpr (BT == BatchType::Single) {

            int info;
            cusolverDnXsyevd(
                handle,
                params,
                CUSOLVER_EIG_MODE_VECTOR,
                enum_convert<BackendLibrary::CUSOLVER>(uplo),
                descrA.rows_,
                BackendScalar<T,B>::type,
                get_data(descrA),
                descrA.ld_,
                BackendScalar<T,B>::type,
                eigenvalues.data(),
                BackendScalar<T,B>::type,
                device_workspace.data(),
                l_work_device,
                host_workspace.data(),
                l_work_host,
                &info);
        } else {
            auto info = pool.allocate<int>(ctx, descrA.batch_size_);
            cusolverDnXsyevBatched(
                handle,
                params,
                CUSOLVER_EIG_MODE_VECTOR,
                enum_convert<BackendLibrary::CUSOLVER>(uplo),
                descrA.rows_,
                BackendScalar<T,B>::type,
                descrA.data_,
                descrA.ld_,
                BackendScalar<T,B>::type,
                eigenvalues.data(),
                BackendScalar<T,B>::type,
                device_workspace.data(),
                l_work_device,
                host_workspace.data(),
                l_work_host,
                info.data(),
                descrA.batch_size_);
        }
        return ctx.get_event();
    }

    template <Backend B, typename T, BatchType BT>
    size_t syev_buffer_size(SyclQueue& ctx,
                            DenseMatView<T,BT> descrA,                            
                            Span<T> eigenvalues,
                            Uplo uplo) {
        static LinalgHandle<B> handle;
        handle.setStream(ctx);
        size_t l_work_device;
        size_t l_work_host;
        size_t total = 0;
        cusolverDnParams_t params;
        cusolverDnCreateParams(&params);
        if constexpr (BT == BatchType::Single){
            cusolverDnXsyevd_bufferSize(handle, params, CUSOLVER_EIG_MODE_VECTOR, enum_convert<BackendLibrary::CUSOLVER>(uplo), descrA.rows_, BackendScalar<T,B>::type, descrA.data_, descrA.ld_, BackendScalar<T,B>::type, eigenvalues.data(), BackendScalar<T,B>::type, &l_work_device, &l_work_host);
        } else {
             cusolverDnXsyevBatched_bufferSize(handle, params, CUSOLVER_EIG_MODE_VECTOR, enum_convert<BackendLibrary::CUSOLVER>(uplo), descrA.rows_, BackendScalar<T,B>::type, descrA.data_, descrA.ld_, BackendScalar<T,B>::type, eigenvalues.data(), BackendScalar<T,B>::type, &l_work_device, &l_work_host, descrA.batch_size_);
             total = BumpAllocator::allocation_size<int>(ctx, descrA.batch_size_);
        }

        return BumpAllocator::allocation_size<std::byte>(ctx, l_work_host) + BumpAllocator::allocation_size<std::byte>(ctx, l_work_device) + total;
    }

    template <Backend B, typename T, BatchType BT>
    SyclEvent ortho(SyclQueue& ctx,
                    DenseMatView<T,BT> A,
                    Transpose transA,
                    Span<std::byte> workspace,
                    OrthoAlgorithm algo) {
        static LinalgHandle<B> handle;
        handle.setStream(ctx);
        BumpAllocator pool(workspace);
        //auto buffer_size = ortho_buffer;
        size_t size = 0;


    }

                   
    #define GEMM_INSTANTIATE(fp, BT) \
    template SyclEvent gemm<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatView<fp, BT>, \
        DenseMatView<fp, BT>, \
        DenseMatView<fp, BT>, \
        fp, fp, Transpose, Transpose, ComputePrecision);

    #define TRSM_INSTANTIATE(fp, BT) \
    template SyclEvent trsm<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatView<fp, BT>, \
        DenseMatView<fp, BT>, \
        Side, Uplo, Transpose, Diag, fp);

    #define POTRF_INSTANTIATE(fp, BT) \
    template SyclEvent potrf<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatView<fp, BT>, \
        Uplo, \
        Span<std::byte>);
    
    #define POTRF_BUFFER_SIZE_INSTANTIATE(fp, BT) \
    template size_t potrf_buffer_size<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatView<fp, BT>, \
        Uplo);

    #define SPMM_INSTANTIATE(fp, F, BT) \
    template SyclEvent spmm<Backend::CUDA, fp, F, BT>( \
        SyclQueue&, \
        SparseMatHandle<fp, F, BT>&, \
        DenseMatView<fp, BT>, \
        DenseMatView<fp, BT>, \
        fp, fp, Transpose, Transpose, Span<std::byte>);
    
    #define SPMM_BUFFER_SIZE_INSTANTIATE(fp, F, BT) \
    template size_t spmm_buffer_size<Backend::CUDA, fp, F, BT>( \
        SyclQueue&, \
        SparseMatHandle<fp, F, BT>&, \
        DenseMatView<fp, BT>, \
        DenseMatView<fp, BT>, \
        fp, fp, Transpose, Transpose);

    #define SYEV_INSTANTIATE(fp, BT) \
    template SyclEvent syev<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatView<fp, BT>, \
        Span<fp>, \
        Uplo, \
        Span<std::byte>);

    #define SYEV_BUFFER_SIZE_INSTANTIATE(fp, BT) \
    template size_t syev_buffer_size<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatView<fp, BT>, \
        Span<fp>, \
        Uplo);
    
    #define BLAS_INSTANTIATE(fp, BT) \
        GEMM_INSTANTIATE(fp, BT) \
        TRSM_INSTANTIATE(fp, BT) \
        POTRF_INSTANTIATE(fp, BT) \
        POTRF_BUFFER_SIZE_INSTANTIATE(fp, BT) \
        SPMM_INSTANTIATE(fp, Format::CSR, BT) \
        SPMM_BUFFER_SIZE_INSTANTIATE(fp, Format::CSR, BT) \
        SYEV_INSTANTIATE(fp, BT) \
        SYEV_BUFFER_SIZE_INSTANTIATE(fp, BT) 
        
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


