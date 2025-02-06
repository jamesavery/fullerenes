#include "../../../linalg-impl.hh"
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <CL/sycl.hpp>

//Signature adapter


namespace linalg {
    template <Backend B, typename T, BatchType BT>
    SyclEvent gemm(SyclQueue& ctx,
                   MatHandle<T,BT>& descrA,
                   MatHandle<T,BT>& descrB,
                   MatHandle<T,BT>& descrC,
                   T alpha,
                   T beta,
                   Transpose transA,
                   Transpose transB) {
        // Call cuBLAS
        static LinalgHandle<B> handle;
        handle.setStream(ctx);
        //cublasSetStream(handle, sycl::get_native<sycl::backend::ext_oneapi_cuda>(*ctx));

        const auto opA = backendTransposeOp(transA);
        const auto opB = backendTransposeOp(transB);
        const auto m = transA == Transpose::NoTrans ? descrA.rows_ : descrA.cols_;
        const auto n = transB == Transpose::NoTrans ? descrB.cols_ : descrB.rows_;
        const auto k = transA == Transpose::NoTrans ? descrA.cols_ : descrA.rows_;        

        if constexpr (BT == BatchType::Single) {
            auto status = cublasGemmEx(handle,
                         opA, opB,
                            m, n, k,
                         &alpha,
                         descrA.data_, BackendScalar<T,B>::type, descrA.ld_,
                         descrB.data_, BackendScalar<T,B>::type, descrB.ld_,
                         &beta,
                         descrC.data_, BackendScalar<T,B>::type, descrC.ld_,
                         BackendScalar<T,B>::type,
                         CUBLAS_GEMM_DFALT);

            if (status != CUBLAS_STATUS_SUCCESS) {
                std::cerr << "cuBLAS GEMM failed with status: " << status << std::endl;
                throw std::runtime_error("cuBLAS GEMM failed");
            }
        } else {
            auto status = cublasGemmStridedBatchedEx(handle,
                                       opA, opB,
                                        m, n, k,
                                       &alpha,
                                       descrA.data_, BackendScalar<T,B>::type, descrA.ld_, descrA.stride_,
                                       descrB.data_, BackendScalar<T,B>::type, descrB.ld_, descrB.stride_,
                                       &beta,
                                       descrC.data_, BackendScalar<T,B>::type, descrC.ld_, descrC.stride_,
                                       descrA.batch_size_,
                                       BackendScalar<T,B>::type,
                                       CUBLAS_GEMM_DFALT);
            if (status != CUBLAS_STATUS_SUCCESS) {
                std::cerr << "cuBLAS GEMM failed with status: " << status << std::endl;
                throw std::runtime_error("cuBLAS GEMM failed");
            }
        }

        return ctx.get_event();
    }

    template <Backend B, typename T, BatchType BT>
    SyclEvent trsm(SyclQueue& ctx,
                   MatHandle<T,BT>& descrA,
                   MatHandle<T,BT>& descrB,
                   Side side,
                   Uplo uplo,
                   Transpose transA,
                   Diag diag,
                   T alpha) {}


                   
    #define GEMM_INSTANTIATE(fp, BT) \
    template SyclEvent gemm<Backend::CUDA, fp, BT>( \
        SyclQueue&, \
        MatHandle<fp, BT>&, \
        MatHandle<fp, BT>&, \
        MatHandle<fp, BT>&, \
        fp, fp, Transpose, Transpose);

    // Macro that covers all layout and batch type combinations for a given floating-point type.
    #define GEMM_INSTANTIATE_FOR_FP(fp)                \
        GEMM_INSTANTIATE(fp, BatchType::Single)  \
        GEMM_INSTANTIATE(fp, BatchType::Batched)

    // Instantiate for the floating-point types of interest.
    GEMM_INSTANTIATE_FOR_FP(float)
    GEMM_INSTANTIATE_FOR_FP(double)
    GEMM_INSTANTIATE_FOR_FP(std::complex<float>)
    GEMM_INSTANTIATE_FOR_FP(std::complex<double>)

    #undef GEMM_INSTANTIATE
    #undef GEMM_INSTANTIATE_FOR_FP
}


