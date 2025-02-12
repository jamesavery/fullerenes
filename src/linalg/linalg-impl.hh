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


    // Helper type traits to identify our enum types
    template<typename T>
    struct is_linalg_enum : std::false_type {};
    
    template<> struct is_linalg_enum<Side> : std::true_type {};
    template<> struct is_linalg_enum<Uplo> : std::true_type {};
    template<> struct is_linalg_enum<Transpose> : std::true_type {};
    template<> struct is_linalg_enum<Diag> : std::true_type {};

    template<class T>
    struct is_complex_or_floating_point : std::is_floating_point<T> { };

    template<class T>
    struct is_complex_or_floating_point<std::complex<T>> : std::is_floating_point<T> { };

    template<typename T>
    struct base_type {
        using type = T;
    };

    template<typename T>
    struct base_type<std::complex<T>> {
        using type = T;
    };

    // Individual enum conversions for CUDA backend
    template<Backend B>
    constexpr auto enum_convert(Side side) {
        if constexpr (B == Backend::CUDA) {
            return static_cast<cublasSideMode_t>(
                side == Side::Left ? CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT
            );
        }
    }

    template<Backend B>
    constexpr auto enum_convert(Uplo uplo) {
        if constexpr (B == Backend::CUDA) {
            return static_cast<cublasFillMode_t>(
                uplo == Uplo::Upper ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER
            );
        }
    }

    template<Backend B>
    constexpr auto enum_convert(Transpose trans) {
        if constexpr (B == Backend::CUDA) {
            return static_cast<cublasOperation_t>(
                trans == Transpose::NoTrans ? CUBLAS_OP_N : CUBLAS_OP_T
            );
        }
    }

    template<Backend B>
    constexpr auto enum_convert(Diag diag) {
        if constexpr (B == Backend::CUDA) {
            return static_cast<cublasDiagType_t>(
                diag == Diag::Unit ? CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT
            );
        }
    }

    template<Backend B, typename T>
    constexpr auto enum_convert(ComputePrecision precision) {
        if constexpr (B == Backend::CUDA) {
            using BaseType = typename base_type<T>::type;
            
            // Handle default precision first
            if (precision == ComputePrecision::Default) {
                if constexpr (std::is_same_v<BaseType, float>) {
                    return CUBLAS_COMPUTE_32F;
                } else if constexpr (std::is_same_v<BaseType, double>) {
                    return CUBLAS_COMPUTE_64F;
                }
            }

            // Handle specific precision requests
            if constexpr (std::is_same_v<BaseType, float>) {
                switch (precision) {
                    case ComputePrecision::F32:  return CUBLAS_COMPUTE_32F;
                    case ComputePrecision::F16:  return CUBLAS_COMPUTE_32F_FAST_16F;
                    case ComputePrecision::BF16: return CUBLAS_COMPUTE_32F_FAST_16BF;
                    case ComputePrecision::TF32: return CUBLAS_COMPUTE_32F_FAST_TF32;
                    default:
                        throw std::runtime_error("Unsupported precision for single precision type");
                }
            } 
            else if constexpr (std::is_same_v<BaseType, double>) {
                if (precision == ComputePrecision::F64) {
                    return CUBLAS_COMPUTE_64F;
                }
                throw std::runtime_error("Only F64 precision supported for double precision type");
            }
        }
        throw std::runtime_error("Unsupported backend or type combination");
    }

    template<Backend B, typename T>
    constexpr auto ptr_convert(T** ptr) {
        if constexpr (B == Backend::CUDA) {
            if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
                return ptr; // No conversion needed
            } else if constexpr (std::is_same_v<T, std::complex<float>>) {
                return reinterpret_cast<cuComplex**>(ptr);
            } else if constexpr (std::is_same_v<T, std::complex<double>>) {
                return reinterpret_cast<cuDoubleComplex**>(ptr);
            }
        }
    }

    


    template<Backend B, typename T>
    constexpr auto ptr_convert(T* ptr) {
        if constexpr (B == Backend::CUDA) {
            if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
                return ptr; // No conversion needed
            } else if constexpr (std::is_same_v<T, std::complex<float>>) {
                return reinterpret_cast<cuComplex*>(ptr);
            } else if constexpr (std::is_same_v<T, std::complex<double>>) {
                return reinterpret_cast<cuDoubleComplex*>(ptr);
            }
        }
    }

    // Const pointer version
    template<Backend B, typename T>
    constexpr auto ptr_convert(const T* ptr) {
        if constexpr (B == Backend::CUDA) {
            if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
                return ptr; // No conversion needed
            } else if constexpr (std::is_same_v<T, std::complex<float>>) {
                return reinterpret_cast<const cuComplex*>(ptr);
            } else if constexpr (std::is_same_v<T, std::complex<double>>) {
                return reinterpret_cast<const cuDoubleComplex*>(ptr);
            }
        }
    }

    /* template<Backend B, typename T>
    constexpr auto ptr_convert(T* ptr) {
        if constexpr (B == Backend::CUDA) {
            using BaseT = std::remove_const_t<T>;
            using CudaT = std::conditional_t<std::is_same_v<BaseT, float> || std::is_same_v<BaseT, double>,
                                           BaseT,
                                           std::conditional_t<std::is_same_v<BaseT, std::complex<float>>,
                                                            cuComplex,
                                                            cuDoubleComplex>>;
            using ReturnT = std::conditional_t<std::is_const_v<T>,
                                             const CudaT*,
                                             CudaT*>;
            
            if constexpr (std::is_same_v<BaseT, float> || std::is_same_v<BaseT, double>) {
                return ptr; // No conversion needed
            } else {
                return reinterpret_cast<ReturnT>(ptr);
            }
        }
    } */

    // Variadic template for converting multiple pointers
    namespace detail {
        template <Backend B, typename T>
        constexpr auto convert_arg(T&& arg) {
            if constexpr (is_linalg_enum<std::remove_reference_t<T>>::value) {
                return enum_convert<B>(std::forward<T>(arg));
            } else if constexpr (std::is_integral_v<std::remove_reference_t<T>>) {
                return std::forward<T>(arg);
            } else if constexpr (std::is_integral_v<std::remove_pointer_t<std::remove_reference_t<T>>>) {
                return std::forward<T>(arg);
            } else if constexpr (std::is_pointer_v<std::remove_reference_t<T>>) {
                return ptr_convert<B>(std::forward<T>(arg));
            } else {
                return std::forward<T>(arg);
            }
        }

        template <typename F, typename Tuple, std::size_t... I>
        auto apply_tuple_impl(F&& f, Tuple&& t, std::index_sequence<I...>) {
            return f(std::get<I>(std::forward<Tuple>(t))...);
        }

        template <typename F, typename Tuple>
        auto apply_tuple(F&& f, Tuple&& t) {
            return apply_tuple_impl(
                std::forward<F>(f),
                std::forward<Tuple>(t),
                std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>{}
            );
        }
    }

    // Combined enum and pointer conversion
    template <Backend B, typename... Args>
    auto backend_convert(Args&&... args) {
        return std::make_tuple(detail::convert_arg<B>(std::forward<Args>(args))...);
    }
    
    template <typename T, Backend B, typename Fun1, typename Fun2, typename Fun3, typename Fun4, typename... Args>
    auto call_backend(const Fun1& fun1, const Fun2& fun2, const Fun3& fun3, const Fun4& fun4, const LinalgHandle<B>& handle, Args&&... args) {
        if constexpr (std::is_same_v<T,float>) {
            return std::apply(fun1, std::tuple_cat(std::forward_as_tuple(handle), backend_convert<B>(std::forward<Args>(args)...)));
        } else if constexpr (std::is_same_v<T,double>) {
            return std::apply(fun2, std::tuple_cat(std::forward_as_tuple(handle), backend_convert<B>(std::forward<Args>(args)...)));
        } else if constexpr (std::is_same_v<T,std::complex<float>>) {
            return std::apply(fun3, std::tuple_cat(std::forward_as_tuple(handle), backend_convert<B>(std::forward<Args>(args)...)));
        } else if constexpr (std::is_same_v<T,std::complex<double>>) {
            return std::apply(fun4, std::tuple_cat(std::forward_as_tuple(handle), backend_convert<B>(std::forward<Args>(args)...)));
        }
    }

    // Variadic template for converting multiple enums
/*     template <Backend B, typename... Args>
    auto enum_convert(Args&&... args) {
        static_assert((is_linalg_enum<std::remove_reference_t<Args>>::value && ...), 
                     "All arguments must be linalg enum types");
        return std::make_tuple(enum_convert<B>(std::forward<Args>(args))...);
    }
 */


        


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

        template <typename T>
        struct BlasComputeType<T, ComputePrecision::Default, Backend::CUDA> {
            static constexpr cublasComputeType_t type = (std::is_same_v<T, float> || std::is_same_v<T, std::complex<float>>) ? CUBLAS_COMPUTE_32F : CUBLAS_COMPUTE_64F;
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