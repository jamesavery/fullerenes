#include <fullerenes/linalg.hh>
#include "../sycl/queue-impl.cc"
#include <CL/sycl.hpp>

#define USE_CUDA
#ifdef USE_CUDA
    #include <cuda_runtime.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
    #include <cusparse.h>
    #include <cusolverDn.h>
#endif

#ifdef USE_ROCM
    #include <hip/hip_runtime.h>
    #include <hip/hip_runtime_api.h>
    #include <rocblas.h>
#endif

#ifdef USE_MAGMA
    #include <magma_v2.h>
#endif

#ifdef USE_MKL
    #include <mkl.h>
#endif



namespace linalg{

    template<typename T>
    struct blas_precision;

    template<>
    struct blas_precision<float> {
        static constexpr char value = 'S';
    };

    template<>
    struct blas_precision<double> {
        static constexpr char value = 'D';
    };

    template<>
    struct blas_precision<std::complex<float>> {
        static constexpr char value = 'C';
    };

    template<>
    struct blas_precision<std::complex<double>> {
        static constexpr char value = 'Z';
    };

    // Function existence checker
    template<typename T>
    struct function_exists {
        template<typename U>
        static auto test(U*) -> decltype(std::declval<U>()(), std::true_type{});
        template<typename U>
        static auto test(...) -> std::false_type;
        static constexpr bool value = decltype(test<T>(nullptr))::value;
    };

    // Function name generator
    template<typename T, Backend B, bool Batched>
    struct blas_function {
        template<const char* Name>
        static constexpr auto get() {
            if constexpr (B == Backend::CUDA) {
                using namespace std::string_literals;
                constexpr std::string_view prefix = "cublas";
                constexpr std::string_view batched_suffix = Batched ? "Batched" : "";
                return prefix + blas_precision<T>::value + Name + batched_suffix;
            }
            // Add similar cases for other backends
            else {
                static_assert(std::false_type{}, "Unsupported backend");
            }
        }

        template<const char* Name>
        static constexpr auto get_fallback() {
            if constexpr (Batched) {
                return get<Name, false>(); // Fallback to non-batched version
            } else {
                return nullptr; // No fallback available
            }
        }
    };


    // Macro for function retrieval with fallback
    #define GET_BLAS_FUNC(backend, name, T, batched) \
    []() { \
        constexpr auto func = blas_function<T, backend, batched>::get<#name>(); \
        if constexpr (function_exists<decltype(func)>::value) { \
            return func; \
        } else { \
            constexpr auto fallback = blas_function<T, backend, batched>::get_fallback<#name>(); \
            static_assert(function_exists<decltype(fallback)>::value, \
                         "Neither primary nor fallback function exists"); \
            return fallback; \
        } \
    }()


    template <typename T, Backend B>
    struct BackendScalar;

    template <Transpose T, Backend B>
    struct BackendTranspose;

    template <Backend B>
    struct LinalgHandle;



    #ifdef USE_CUDA
        template <>
        struct BackendScalar<float, Backend::CUDA> {
            static constexpr cudaDataType_t type = CUDA_R_32F;
        };

        template <>
        struct BackendScalar<double, Backend::CUDA> {
            static constexpr cudaDataType_t type = CUDA_R_64F;
        };

        template <>
        struct BackendScalar<std::complex<float>, Backend::CUDA> {
            static constexpr cudaDataType_t type = CUDA_C_32F;
        };

        template <>
        struct BackendScalar<std::complex<double>, Backend::CUDA> {
            static constexpr cudaDataType_t type = CUDA_C_64F;
        };
        /* template <>
        struct BackendTranspose<Transpose::NoTrans, Layout::ColMajor, Backend::CUDA> {
            static constexpr cublasOperation_t type = CUBLAS_OP_N;
        };

        template <>
        struct BackendTranspose<Transpose::NoTrans, Layout::RowMajor, Backend::CUDA> {
            static constexpr cublasOperation_t type = CUBLAS_OP_T;
        };  

        template <>
        struct BackendTranspose<Transpose::Trans, Layout::ColMajor, Backend::CUDA> {
            static constexpr cublasOperation_t type = CUBLAS_OP_T;
        };

        template <>
        struct BackendTranspose<Transpose::Trans, Layout::RowMajor, Backend::CUDA> {
            static constexpr cublasOperation_t type = CUBLAS_OP_N;
        }; */

        inline cublasOperation_t backendTransposeOp(Transpose t) {
            return (t == Transpose::NoTrans) ? CUBLAS_OP_N : CUBLAS_OP_T;
        }

        template <>
        struct LinalgHandle<Backend::CUDA> {
            cublasHandle_t blas_handle_;
            cusparseHandle_t sparse_handle_;
            cusolverDnHandle_t solver_handle_;

            LinalgHandle() {
                auto blas_status = cublasCreate(&blas_handle_);
                if (blas_status != CUBLAS_STATUS_SUCCESS) {
                    std::cerr << "CUBLAS initialization failed with status: " << blas_status << std::endl;
                    throw std::runtime_error("CUBLAS initialization failed");
                }

                auto sparse_status = cusparseCreate(&sparse_handle_);
                if (sparse_status != CUSPARSE_STATUS_SUCCESS) {
                    std::cerr << "CUSPARSE initialization failed with status: " << sparse_status << std::endl;
                    throw std::runtime_error("CUSPARSE initialization failed");
                }

                auto solver_status = cusolverDnCreate(&solver_handle_);
                if (solver_status != CUSOLVER_STATUS_SUCCESS) {
                    std::cerr << "CUSOLVER initialization failed with status: " << solver_status << std::endl;
                    throw std::runtime_error("CUSOLVER initialization failed");
                }
            }

            ~LinalgHandle() {
                cublasDestroy(blas_handle_);
                cusparseDestroy(sparse_handle_);
                cusolverDnDestroy(solver_handle_);
            }

            constexpr inline operator cublasHandle_t() const {
                return blas_handle_;
            }

            constexpr inline operator cusparseHandle_t() const {
                return sparse_handle_;
            }

            constexpr inline operator cusolverDnHandle_t() const {
                return solver_handle_;
            }

            void setStream(const SyclQueue& ctx) {
                cudaStream_t stream = sycl::get_native<sycl::backend::ext_oneapi_cuda>(*ctx);
                cublasSetStream(blas_handle_, stream);
                cusparseSetStream(sparse_handle_, stream);
                cusolverDnSetStream(solver_handle_, stream);
            }
        };
    #endif

    
    

}