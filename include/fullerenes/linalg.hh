#pragma once
#include <complex>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>
#include <memory>
#include <array>  // Add this include at the top with other includes

namespace linalg{
    enum class Backend {
        AUTO,
        CUDA,
        ROCM,
        MKL,
        MAGMA,
        SYCL,
        NETLIB
        // Add more as needed
    };

    enum class BackendLibrary {
        CUBLAS,     //Belongs to CUDA backend
        CUSPARSE,   //Belongs to CUDA backend
        CUSOLVER,   //Belongs to CUDA backend
        ROCBLAS,    //Belongs to ROCM backend
        ROCSPARSE,  //Belongs to ROCM backend
        ROCSOLVER,  //Belongs to ROCM backend
        MAGMA,      //Belongs to MAGMA backend
        MKL,        //Belongs to MKL backend
        NETLIB,     //Belongs to NETLIB backend
        LAPACK      //Belongs to NETLIB backend
    };

    enum class BatchType {
        Single,
        Batched
    };
    
    enum class Transpose {
        NoTrans,
        Trans
    };

    enum class Uplo {
        Upper,
        Lower
    };

    enum class Diag {
        NonUnit,
        Unit
    };

    enum class Side {
        Left,
        Right
    };
    
    enum class OrthoAlgorithm {
        Chol2,      //Default
        Cholesky,   //Rarely sufficient
        ShiftChol3, //More stable than Chol2
        Householder, //Most numerically stable
        ModifiedGramSchmidt
    };
    
    //Some of the types are not supported by all backends, compilation errors will make this apparent
    enum class ComputePrecision {
        Default, //Use same precision as input
        F32,
        F64,
        F16,
        BF16,
        TF32
    };


    enum class Format {
        COO,        //Coordinate
        CSR,        //Compressed Sparse Row
        CSC,        //Compressed Sparse Column
        SELL,       //Sliced ELLPACK
        BSR,        //Blocked Sparse Row
        BLOCKED_ELL //Blocked ELLPACK
    };

    enum class Layout {
        RowMajor,
        ColMajor
    };

    namespace detail {
        // Tag types for overload resolution
        template<Backend B> struct backend_tag {};
        struct fallback_tag : backend_tag<Backend::AUTO> {};
    }

    struct BackendSelector {
        // Count available backends at compile time
        static constexpr size_t num_backends = 1   // NETLIB is always available
            #ifdef USE_CUDA
                + 1
            #endif
            #ifdef USE_ROCM
                + 1
            #endif
            #ifdef USE_MKL
                + 1
            #endif
            #ifdef USE_MAGMA
                + 1
            #endif
            #ifdef USE_SYCL
                + 1
            #endif
        ;

        static constexpr std::array<Backend, num_backends> available_backends = {{
                #ifdef USE_CUDA
                    Backend::CUDA,
                #endif
                #ifdef USE_ROCM
                    Backend::ROCM,
                #endif
                #ifdef USE_MKL
                    Backend::MKL,
                #endif
                #ifdef USE_MAGMA
                    Backend::MAGMA,
                #endif
                #ifdef USE_SYCL
                    Backend::SYCL,
                #endif
                Backend::NETLIB  // Always available as fallback
        }};

        static constexpr Backend get(SyclQueue& ctx) {
            auto device = ctx.device();
            switch (device.type) {
                case DeviceType::CPU:
                    #ifdef USE_MKL
                        return Backend::MKL;
                    #else
                        return Backend::NETLIB;
                    #endif

                case DeviceType::GPU:
                    if (device.get_vendor() == "NVIDIA Corporation") {
                        #ifdef USE_CUDA
                            return Backend::CUDA;
                        #elif defined(USE_MAGMA)
                            return Backend::MAGMA;
                        #else
                            return Backend::SYCL;  // SYCL fallback for NVIDIA
                        #endif
                    }
                    else if (device.get_vendor() == "Advanced Micro Devices, Inc.") {
                        #ifdef USE_ROCM
                            return Backend::ROCM;
                        #else
                            return Backend::SYCL;  // SYCL fallback for AMD
                        #endif
                    }
                    return Backend::SYCL;  // SYCL fallback for other GPUs

                case DeviceType::ACCELERATOR:
                    return Backend::SYCL;  // SYCL is primary backend for accelerators

                case DeviceType::HOST:
                    return Backend::NETLIB;

                default:
                    return Backend::AUTO;
            }
        }

        template <typename... Args>
        static constexpr auto select(SyclQueue& ctx, Args&&... args) {
            constexpr size_t num_args = sizeof...(Args);
            static_assert(num_args > 0, "At least one argument must be provided");
            
            // Convert args to array for indexed access
            std::array<std::common_type_t<Args...>, num_args> implementations{std::forward<Args>(args)...};
            
            // Get the backend for this device
            Backend backend = get(ctx);
            
            // Find the implementation index based on backend priority
            size_t impl_index = 0;
            for(size_t i = 0; i < available_backends.size() && i < num_args; ++i) {
                if(available_backends[i] == backend) {
                    impl_index = i;
                    break;
                }
            }
            
            // If no matching backend found, use last argument as fallback
            if(impl_index >= num_args) {
                impl_index = num_args - 1;
            }
            
            return implementations[impl_index];
        }
    };

    //Forward declarations
    template <typename T, BatchType BT>
    struct DenseMatHandle;

    template <typename T, Format F, BatchType BT>
    struct SparseMatHandle;

    template <typename T>
    struct BackendDenseMatrixHandle;

    template <typename T, Format F>
    struct BackendSparseMatrixHandle;

    template <typename T>
    struct BackendDenseVectorHandle;

    template <typename T>
    struct BackendSparseVectorHandle;
    
    template <typename T, BatchType BT>
    struct DenseMatView;

    template <typename T>
    struct DenseMatHandle<T, BatchType::Single> {
        DenseMatHandle(T* data, int rows, int cols, int ld);
        ~DenseMatHandle();
        void init(SyclQueue& ctx);
        void init_backend();
        // Accessors...
        T* data_;
        int rows_, cols_, ld_;
        Layout layout_ = Layout::ColMajor; //Most backends don't support row-major dense matrices

        BackendDenseMatrixHandle<T>* operator->();
        BackendDenseMatrixHandle<T>& operator*();
        
        DenseMatView<T, BatchType::Single> operator()() const {
            return DenseMatView<T, BatchType::Single>(*this);
        }

        private:
            std::unique_ptr<BackendDenseMatrixHandle<T>> backend_handle_;
    };

    template <typename T>
    struct DenseMatHandle<T, BatchType::Batched> {
        DenseMatHandle(T* data, int rows, int cols, int ld, int stride, int batch_size);
        ~DenseMatHandle();
        void init(SyclQueue& ctx);
        void init_backend();

        // Accessors...
        T* data_;
        SyclVector<T*> data_ptrs_;
        int rows_, cols_, ld_, stride_, batch_size_;
        Layout layout_ = Layout::ColMajor; //Most backends don't support row-major dense matrices

        BackendDenseMatrixHandle<T>* operator->();
        BackendDenseMatrixHandle<T>& operator*();

        DenseMatView<T, BatchType::Batched> operator()() const {
            return DenseMatView<T, BatchType::Batched>(*this);
        }

        private:
            std::unique_ptr<BackendDenseMatrixHandle<T>> backend_handle_;
    };

    

    template <typename T>
    struct DenseMatView<T, BatchType::Batched> {
        DenseMatView(T* data, int rows, int cols, int ld, int stride, int batch_size, Span<T*> data_ptrs);
        DenseMatView(const DenseMatView<T, BatchType::Batched>& view);
        DenseMatView(DenseMatView<T, BatchType::Batched>&& view);
        DenseMatView<T, BatchType::Batched>& operator=(const DenseMatView<T, BatchType::Batched>& view);
        DenseMatView<T, BatchType::Batched>& operator=(DenseMatView<T, BatchType::Batched>&& view);
        // Allow lvalue reference construction but explicitly delete rvalue reference constructor
        DenseMatView(const DenseMatHandle<T, BatchType::Batched>& handle);
        // Allow for the view to be a reinterpreted view of the matrices
        DenseMatView(const DenseMatHandle<T, BatchType::Batched>& handle, int rows, int cols, int ld, int stride, int batch_size);
        DenseMatView(DenseMatHandle<T, BatchType::Batched>&& handle) = delete;

        // Allow lvalue reference assignment but explicitly delete rvalue reference assignment
        DenseMatView<T, BatchType::Batched>& operator=(const DenseMatHandle<T, BatchType::Batched>& handle);
        DenseMatView<T, BatchType::Batched>& operator=(DenseMatHandle<T, BatchType::Batched>&& handle) = delete;
        
        ~DenseMatView();
        void init(SyclQueue& ctx);
        void init_backend();
        // Accessors...
        T* data_; 
        Span<T*> data_ptrs_;
        int rows_, cols_, ld_, stride_, batch_size_;
        Layout layout_ = Layout::ColMajor; //Most backends don't support row-major dense matrices

        BackendDenseMatrixHandle<T>* operator->();
        BackendDenseMatrixHandle<T>& operator*();

        private:
            std::shared_ptr<BackendDenseMatrixHandle<T>> backend_handle_;
    };
    
    

    template <typename T>
    struct DenseMatView<T, BatchType::Single> {
        DenseMatView(T* data, int rows, int cols, int ld);
        DenseMatView(const DenseMatView<T, BatchType::Single>& view);
        DenseMatView(DenseMatView<T, BatchType::Single>&& view);
        DenseMatView<T, BatchType::Single>& operator=(const DenseMatView<T, BatchType::Single>& view);
        DenseMatView<T, BatchType::Single>& operator=(DenseMatView<T, BatchType::Single>&& view);
        
        //These handles are already non-owning, but for consistency, delete rvalue reference constructor
        DenseMatView(const DenseMatHandle<T, BatchType::Single>& handle);
        DenseMatView(const DenseMatHandle<T, BatchType::Single>& handle, int rows, int cols, int ld);
        DenseMatView(DenseMatHandle<T, BatchType::Single>&& handle) = delete;

        //These handles are already non-owning, but for consistency, delete rvalue reference assignment
        DenseMatView<T, BatchType::Single>& operator=(const DenseMatHandle<T, BatchType::Single>& handle);
        DenseMatView<T, BatchType::Single>& operator=(DenseMatHandle<T, BatchType::Single>&& handle) = delete;

        ~DenseMatView();
        void init(SyclQueue& ctx);
        void init_backend();

        // Accessors...
        T* data_;
        int rows_, cols_, ld_;
        Layout layout_ = Layout::ColMajor;

        BackendDenseMatrixHandle<T>* operator->();
        BackendDenseMatrixHandle<T>& operator*();

        private:
            std::shared_ptr<BackendDenseMatrixHandle<T>> backend_handle_;
    };

    // Deduction guides for DenseMatView

    template <typename T>
    DenseMatView(const DenseMatHandle<T, BatchType::Single>&) -> DenseMatView<T, BatchType::Single>;

    template <typename T>
    DenseMatView(const DenseMatHandle<T, BatchType::Batched>&) -> DenseMatView<T, BatchType::Batched>;

    template <typename T>
    struct SparseMatHandle<T, Format::CSR, BatchType::Single> {
        /**
         * @brief Constructs a single CSR sparse matrix handle
         * @param data Array of non-zero values [nnz]
         * @param row_offsets Array of row offsets [rows + 1]
         * @param col_indices Array of column indices [nnz]
         * @param nnz Number of non-zero elements
         * @param rows Number of rows
         * @param cols Number of columns
         * @param layout Layout of the matrix
         */
        SparseMatHandle(T* data, int* row_offsets, int* col_indices, int nnz, int rows, int cols, Layout layout = Layout::RowMajor);
        ~SparseMatHandle();
        void init(SyclQueue& ctx);
        void init_backend();

        // Raw pointers to externally owned memory
        T* data_;              // [nnz] non-zero values
        int* row_offsets_;     // [rows + 1] row offsets
        int* col_indices_;     // [nnz] column indices
        int nnz_, rows_, cols_;
        Layout layout_ = Layout::RowMajor;

        BackendSparseMatrixHandle<T, Format::CSR>* operator->();
        BackendSparseMatrixHandle<T, Format::CSR>& operator*();

        private:
            std::unique_ptr<BackendSparseMatrixHandle<T, Format::CSR>> backend_handle_;

        
    };

    template <typename T>
    struct SparseMatHandle<T, Format::CSR, BatchType::Batched> {
        /**
         * @brief Constructs a batched CSR sparse matrix handle
         * @param data Array of non-zero values [batch_size * stride]
         * @param row_offsets Array of row offsets [batch_size * (rows + 1)]
         * @param col_indices Array of column indices [batch_size * nnz]
         * @param nnz Number of non-zero elements per matrix
         * @param rows Number of rows per matrix
         * @param cols Number of columns per matrix
         * @param stride Stride between consecutive matrices in data array
         * @param batch_size Number of matrices in batch
         */
        SparseMatHandle(T* data, int* row_offsets, int* col_indices, 
            int nnz, int rows, int cols, int stride, int batch_size);

        ~SparseMatHandle();
        void init(SyclQueue& ctx);
        void init_backend();

        // Raw pointers to externally owned memory
        T* data_;              // [batch_size * stride] non-zero values
        int* row_offsets_;     // [batch_size * (rows + 1)] row offsets
        int* col_indices_;     // [batch_size * nnz] column indices
        
        int nnz_, rows_, cols_, stride_, batch_size_;
        Layout layout_ = Layout::RowMajor;

        BackendSparseMatrixHandle<T, Format::CSR>* operator->();
        BackendSparseMatrixHandle<T, Format::CSR>& operator*();

        private:
            std::unique_ptr<BackendSparseMatrixHandle<T, Format::CSR>> backend_handle_;
    };

    //Uniform accessor for data
    template <template <typename, BatchType> class Handle, typename T, BatchType BT>
    auto get_data(Handle<T,BT>& handle) {
        return handle.data_;
    }
    template <template <typename, BatchType> class Handle, typename T, BatchType BT, std::enable_if_t<BT == BatchType::Batched, int> = 0>
    auto get_ptr_arr(SyclQueue& ctx, Handle<T,BT>& handle) {
        handle.init(ctx);
        return handle.data_ptrs_.data();

    }

    //Uniform accessor for data
    template <template <typename, Format, BatchType> class Handle, typename T, Format F, BatchType BT>
    auto get_data(Handle<T,F,BT>& handle) {
        return handle.data_;
    }


    template <typename T, BatchType BT>
    struct DenseVecHandle;

    template <typename T, BatchType BT>
    struct SparseVecHandle;

    template <typename T>
    struct DenseVecHandle<T, BatchType::Single> {
        DenseVecHandle(T* data, int size, int ldc) : data_(data), size_(size), ldc_(ldc) {}
        // Accessors...
        T* data_;
        SyclVector<T*> data_ptrs_;
        int size_, ldc_;

        private:
            std::unique_ptr<BackendDenseVectorHandle<T>> backend_handle_;
    };

    template <typename T>
    struct DenseVecHandle<T, BatchType::Batched> {
        DenseVecHandle(T* data, int size, int ldc, int stride, int batch_size) 
            : data_(data), size_(size), ldc_(ldc), stride_(stride), batch_size_(batch_size) {}
        // Accessors...
        T* data_;
        int size_, ldc_, stride_, batch_size_;

        private:
            std::unique_ptr<BackendDenseVectorHandle<T>> backend_handle_;
    };

    template <typename T>
    struct SparseVecHandle<T, BatchType::Single> {
        SparseVecHandle(T* data, int* indices, int size) : data_(data), indices_(indices), size_(size) {}
        // Accessors...
        T* data_;
        int* indices_;
        int size_;

        private:
            std::unique_ptr<BackendDenseVectorHandle<T>> backend_handle_;
    };

    template <typename T>
    struct SparseVecHandle<T, BatchType::Batched> {
        SparseVecHandle(T* data, int* indices, int size, int stride, int batch_size) 
            : data_(data), indices_(indices), size_(size), stride_(stride), batch_size_(batch_size) {}
        // Accessors...
        T* data_;
        int* indices_;
        int size_, stride_, batch_size_;

        private:
            std::unique_ptr<BackendSparseVectorHandle<T>> backend_handle_;
    };


    namespace detail {
        // Type trait to check for complex or floating point types
        template<typename T>
        struct is_complex_or_floating_point : 
            std::bool_constant<std::is_floating_point_v<T> || 
                             std::is_same_v<T, std::complex<float>> || 
                             std::is_same_v<T, std::complex<double>>> {};

        template<Backend B>
        [[noreturn]] void throw_unsupported() {
            throw std::runtime_error("Operation not supported for selected backend: " + std::to_string(static_cast<int>(B)));
        }

        template<typename T>
        using enable_if_scalar_t = typename std::enable_if<
            is_complex_or_floating_point<T>::value
        >::type;

        // Base template for CRTP-style function declarations
        /* template<typename Name, typename RetType>
        struct LinalgFuncBase {
            // Default implementation that throws
            template <Backend B, typename T, BatchType BT, typename... Args>
            static RetType impl(detail::fallback_tag, Args&&... args) {
                detail::throw_unsupported<B>();
            }

            // Backend-specific version declaration
            template <Backend B, typename T, BatchType BT, typename... Args>
            static RetType impl(detail::backend_tag<B>, Args&&... args);

            // Wrapper that forwards to implementation with appropriate tag
            template <Backend B, typename T, BatchType BT, typename... Args>
            static RetType call(Args&&... args) {
                return impl<B, T, BT>(detail::backend_tag<B>{}, std::forward<Args>(args)...);
            }

            // Auto-dispatching version
            template <typename T, BatchType BT, typename... Args>
            static RetType dispatch(SyclQueue& ctx, Args&&... args) {
                Backend backend = BackendSelector::get(ctx);
                
                switch (backend) {
                    case Backend::CUDA:
                        #ifdef USE_CUDA
                            return impl<Backend::CUDA, T, BT>(detail::backend_tag<Backend::CUDA>{}, std::forward<Args>(args)...);
                        #endif
                    case Backend::ROCM:
                        #ifdef USE_ROCM
                            return impl<Backend::ROCM, T, BT>(detail::backend_tag<Backend::ROCM>{}, std::forward<Args>(args)...);
                        #endif
                    case Backend::MKL:
                        #ifdef USE_MKL
                            return impl<Backend::MKL, T, BT>(detail::backend_tag<Backend::MKL>{}, std::forward<Args>(args)...);
                        #endif
                    case Backend::MAGMA:
                        #ifdef USE_MAGMA
                            return impl<Backend::MAGMA, T, BT>(detail::backend_tag<Backend::MAGMA>{}, std::forward<Args>(args)...);
                        #endif
                    case Backend::SYCL:
                        #ifdef USE_SYCL
                            return impl<Backend::SYCL, T, BT>(detail::backend_tag<Backend::SYCL>{}, std::forward<Args>(args)...);
                        #endif
                    case Backend::NETLIB:
                        return impl<Backend::NETLIB, T, BT>(detail::backend_tag<Backend::NETLIB>{}, std::forward<Args>(args)...);
                    default:
                        return impl<Backend::AUTO, T, BT>(detail::fallback_tag{}, std::forward<Args>(args)...);
                }
            }
        }; */
    }
    /* struct Gemm : detail::LinalgFuncBase<Gemm, SyclEvent> {};
    struct Trsm : detail::LinalgFuncBase<Trsm, SyclEvent> {};
    struct Potrf : detail::LinalgFuncBase<Potrf, SyclEvent> {};
    struct PotrfBufferSize : detail::LinalgFuncBase<PotrfBufferSize, size_t> {}; */


    template <Backend B, typename T, BatchType BT>
    SyclEvent gemm(SyclQueue& ctx,
        DenseMatView<T,BT> descrA,
        DenseMatView<T,BT> descrB,
        DenseMatView<T,BT> descrC,
        T alpha,
        T beta,
        Transpose transA,        
        Transpose transB,
        ComputePrecision precision = ComputePrecision::Default);

    
    template <Backend B, typename T, Format F, BatchType BT>
    SyclEvent spmm(SyclQueue& ctx,
        SparseMatHandle<T, F, BT>& descrA,
        DenseMatView<T,BT> descrB,
        DenseMatView<T,BT> descrC,
        T alpha,
        T beta,
        Transpose transA,
        Transpose transB,
        Span<std::byte> workspace);

    template <Backend B, typename T, Format F, BatchType BT>
    size_t spmm_buffer_size(SyclQueue& ctx,
        SparseMatHandle<T, F, BT>& descrA,
        DenseMatView<T,BT> descrB,
        DenseMatView<T,BT> descrC,
        T alpha,
        T beta,
        Transpose transA,        
        Transpose transB);

    template <Backend B, typename T, BatchType BT>
    SyclEvent trsm(SyclQueue& ctx,
        DenseMatView<T,BT> descrA,
        DenseMatView<T,BT> descrB,
        Side side,
        Uplo uplo,
        Transpose transA,
        Diag diag,
        T alpha);

    template <Backend B, typename T, BatchType BT>
    size_t potrf_buffer_size(SyclQueue& ctx,
                        DenseMatView<T,BT> A,
                        Uplo uplo);

    template <Backend B, typename T, BatchType BT>
    SyclEvent potrf(SyclQueue& ctx,
            DenseMatView<T,BT> descrA,
            Uplo uplo,
            Span<std::byte> workspace); 

    template <Backend B, typename T, BatchType BT>
    SyclEvent syev(SyclQueue& ctx,
            DenseMatView<T,BT> descrA, //A is overwritten with eigenvectors
            Span<T> eigenvalues,
            Uplo uplo,
            Span<std::byte> workspace);

    template <Backend B, typename T, BatchType BT>
    size_t syev_buffer_size(SyclQueue& ctx,
            DenseMatView<T,BT> A,
            Span<T> eigenvalues,
            Uplo uplo);

    //Produces an orthonormal matrix in-place that spans the column space of A (if transA == NoTrans) or the row space of A (if transA == Trans)
    template <Backend B, typename T, BatchType BT>
    SyclEvent ortho(SyclQueue& ctx,
            DenseMatView<T,BT> A, //A is overwritten with orthogonal vectors
            Transpose transA,
            Span<std::byte> workspace,
            OrthoAlgorithm algo = OrthoAlgorithm::Chol2);

    //Produces an orthonormal matrix in-place that spans the column space of A (if transA == NoTrans) or the row space of A (if transA == Trans)
    //The matrix A is orthogonalized with respect to the columns of M (if transM == NoTrans) or the rows of M (if transM == Trans)
    template <Backend B, typename T, BatchType BT>
    SyclEvent ortho(SyclQueue& ctx,
            DenseMatView<T,BT> A, //A is overwritten with orthogonal vectors
            DenseMatView<T,BT> M, //External metric
            Transpose transA,
            Transpose transM,
            Span<std::byte> workspace,
            OrthoAlgorithm algo = OrthoAlgorithm::Chol2);
    
    template <Backend B, typename T, BatchType BT>
    size_t ortho_buffer_size(SyclQueue& ctx,
            DenseMatView<T,BT> A,
            Transpose transA,
            OrthoAlgorithm algo = OrthoAlgorithm::Chol2);

    template <Backend B, typename T, BatchType BT>
    size_t ortho_buffer_size(SyclQueue& ctx,
            DenseMatView<T,BT> A,
            DenseMatView<T,BT> M,
            Transpose transA,
            Transpose transM,
            OrthoAlgorithm algo = OrthoAlgorithm::Chol2);
}
