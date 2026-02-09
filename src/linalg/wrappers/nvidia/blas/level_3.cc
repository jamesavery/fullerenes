#include "../../../linalg-impl.hh"
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <sycl/sycl.hpp>
#include <complex>

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
        //If transA == NoTrans we find the orthonormal basis of the column space
        //Else we find the orthonormal basis of the row space

        static LinalgHandle<B> handle;
        auto batch_size = get_batch_size(A);
        handle.setStream(ctx);
        BumpAllocator pool(workspace);
        //auto buffer_size = ortho_buffer;
        auto [m, k] = get_effective_dims(A, transA);
        Transpose inv_trans = transA == Transpose::Trans ? Transpose::NoTrans : Transpose::Trans;
        assert(k <= m);
        //If k > m && transA == NoTrans the columns of A are linearly dependent
        //Else if k > m && transA == Trans the rows of A are linearly dependent
        auto ATA =          pool.allocate<T>(ctx, k*k*batch_size);
        auto matAmem =      pool.allocate<T*>(ctx, batch_size);
        auto matATAmem =    pool.allocate<T*>(ctx, batch_size);
        auto ATA_stride = k*k;
        auto C = create_view<T, BT>(ATA.data(), k, k, k, ATA_stride, batch_size, matATAmem);
        auto potrf_workspace = pool.allocate<std::byte>(ctx, potrf_buffer_size<Backend::CUDA>(ctx, C, Uplo::Lower));
        
        auto chol_alg = [&](){
            constexpr T alpha = 1.0;
            constexpr T beta = 0.0;
            
            //Compute StS = S^T * S or StS = S * S^T (depending on transA)
            gemm<Backend::CUDA>(ctx, A, A, C, T(1.0), T(0.0), inv_trans, transA);
            //Compute the Cholesky Factorization of StS
            potrf<Backend::CUDA>(ctx, C, Uplo::Lower, potrf_workspace);
            //Compute Q = S * StS^-1 (S is overwritten with Q)
            trsm<Backend::CUDA>(ctx, C, A, Side::Right, Uplo::Lower, inv_trans, Diag::NonUnit, alpha);
            //Compute the QR factorization of Q
            gemm<Backend::CUDA>(ctx, A, A, C, T(1.0), T(0.0), inv_trans, transA);
        };

        if (algo == OrthoAlgorithm::Cholesky){
            chol_alg();
        } else if (algo == OrthoAlgorithm::Chol2){
            chol_alg();
            chol_alg();
        } else if (algo == OrthoAlgorithm::ShiftChol3) {
            gemm<Backend::CUDA>(ctx, A, A, C, T(1.0), T(0.0), inv_trans, transA);

            auto ATA_ptr = get_data(C);
            ctx -> submit([&](sycl::handler& h){
                h.parallel_for(sycl::nd_range<1>(sycl::range{size_t(batch_size * k)}, sycl::range{size_t(k)}), [=](sycl::nd_item<1> item){
                    auto tid = item.get_local_linear_id();
                    auto bid = item.get_group_linear_id();
                    auto cta = item.get_group();
                    auto ATA_acc = ATA_ptr + bid * ATA_stride;
                    T g_norm = 0.0;
                    if constexpr (sycl::detail::is_complex<T>::value){
                        g_norm = sycl::reduce_over_group(cta, std::sqrt(ATA_acc[tid * k + tid].real()), sycl::maximum<typename T::value_type>());
                    } else {
                        g_norm = sycl::reduce_over_group(cta, std::sqrt(ATA_acc[tid * k + tid]), sycl::maximum<T>());
                    }
                    auto eps = std::numeric_limits<T>::epsilon();
                    auto shift = T(11.0) * T( T(m * k) * T(eps) + T(k + 1) * T(k) * T(eps)) * g_norm;
                    ATA_acc[tid * k + tid] += shift;
                });
            });

            //Compute the Cholesky Factorization of StS
            potrf<Backend::CUDA>(ctx, C, Uplo::Lower, potrf_workspace);
            trsm<Backend::CUDA>(ctx, C, A, Side::Right, Uplo::Lower, inv_trans, Diag::NonUnit, T(1.0));
            chol_alg();
            chol_alg();
        }
        
        return ctx.get_event();
    }

    template <Backend B, typename T, BatchType BT>
    size_t ortho_buffer_size(SyclQueue& ctx,
                             DenseMatView<T,BT> A,
                             Transpose transA,
                             OrthoAlgorithm algo) {
        size_t size = 0;
        auto [m, k] = get_effective_dims(A, transA);
        return  BumpAllocator::allocation_size<std::byte>(ctx, potrf_buffer_size<Backend::CUDA>(ctx, DenseMatView<T, BatchType::Batched>(A.data_, m, k, m, m*k, get_batch_size(A), Span<T*>{}), Uplo::Lower)) + 
        2*BumpAllocator::allocation_size<T*>(ctx, get_batch_size(A)) + 
        BumpAllocator::allocation_size<T>(ctx, k*k*get_batch_size(A));
    }

    template <Backend B, typename T, BatchType BT>
    SyclEvent ortho(SyclQueue& ctx,
                    DenseMatView<T,BT> A,
                    DenseMatView<T,BT> M,
                    Transpose transA,
                    Transpose transM,
                    Span<std::byte> workspace,
                    OrthoAlgorithm algo,
                    size_t iterations) {
        BumpAllocator pool(workspace);
        //When orthogonalizing against an external basis M,
        //M must be an orthonormal basis
        //Both A and M must be either tall-and-skinny or short-and-fat
        //Furthermore the number of vectors in A and M must sum to at most the dimension of these vectors 
        auto nM = transM == Transpose::NoTrans ? M.cols_ : M.rows_;
        auto nA = transA == Transpose::NoTrans ? A.cols_ : A.rows_;
        auto k = transA == Transpose::NoTrans ? A.rows_ : A.cols_;
        assert(nA + nM <= k);
        assert(k == (transM == Transpose::NoTrans ? M.rows_ : M.cols_));

        auto inv_transA = transA == Transpose::Trans ? Transpose::NoTrans : Transpose::Trans;
        auto inv_transM = transM == Transpose::Trans ? Transpose::NoTrans : Transpose::Trans;
        auto MAmem =          pool.allocate<T>(ctx, nM*nA * get_batch_size(A));
        auto orthoworkspace = pool.allocate<std::byte>(ctx, ortho_buffer_size<Backend::CUDA>(ctx, A, transA, algo));
        auto descrMA = create_view<T, BT>(MAmem.data(), nM, nA, nM, nM*nA, get_batch_size(A), Span<T*>{});
        auto isAtrans = transA == Transpose::Trans;
        auto is_first_transposed = static_cast<Transpose>(((transA == Transpose::Trans) || (transM == Transpose::Trans)));
        auto is_second_transposed = static_cast<Transpose>(((transA == Transpose::Trans) && (transM == Transpose::NoTrans)));
        for (size_t i = 0; i < iterations; i++){
            gemm<Backend::CUDA>(ctx, M, A, descrMA, T(1.0), T(0.0), inv_transM, transA);
            gemm<Backend::CUDA>(ctx, isAtrans ? descrMA : M, isAtrans ? M : descrMA, A, T(-1.0), T(1.0), is_first_transposed, is_second_transposed);

            ortho<Backend::CUDA>(ctx, A, transA, orthoworkspace, algo);
        }
        return ctx.get_event();
    }

    template <Backend B, typename T, BatchType BT>
    size_t ortho_buffer_size(SyclQueue& ctx,
                             DenseMatView<T,BT> A,
                             DenseMatView<T,BT> M,
                             Transpose transA,
                             Transpose transM,
                             OrthoAlgorithm algo,
                             size_t iterations) {
        auto nM = transM == Transpose::NoTrans ? M.cols_ : M.rows_;
        auto nA = transA == Transpose::NoTrans ? A.cols_ : A.rows_;
        return  BumpAllocator::allocation_size<std::byte>(ctx, ortho_buffer_size<Backend::CUDA>(ctx, A, transA, algo)) +
                BumpAllocator::allocation_size<T>(ctx, nM*nA * get_batch_size(A));
    }   


    #define ORTHO_INSTANTIATE(fp, BT) \
    template SyclEvent ortho<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatView<fp, BT>, \
        Transpose, \
        Span<std::byte>, \
        OrthoAlgorithm); \
    template SyclEvent ortho<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatView<fp, BT>, \
        DenseMatView<fp, BT>, \
        Transpose, \
        Transpose, \
        Span<std::byte>, \
        OrthoAlgorithm, \
        size_t); \
    template size_t ortho_buffer_size<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatView<fp, BT>, \
        Transpose, \
        OrthoAlgorithm); \
    template size_t ortho_buffer_size<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        DenseMatView<fp, BT>, \
        DenseMatView<fp, BT>, \
        Transpose, \
        Transpose, \
        OrthoAlgorithm, \
        size_t);
                   
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
        ORTHO_INSTANTIATE(fp, BT) \
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


