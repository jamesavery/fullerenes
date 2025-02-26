#pragma once
#include <complex>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

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
    struct DenseMatHandle<T, BatchType::Single> {
        DenseMatHandle(T* data, int rows, int cols, int ld) 
            : data_(data), rows_(rows), cols_(cols), ld_(ld) {}
        // Accessors...
        T* data_;
        int rows_, cols_, ld_;
        Layout layout_ = Layout::ColMajor; //Most backends don't support row-major dense matrices
    };

    template <typename T>
    struct DenseMatHandle<T, BatchType::Batched> {
        DenseMatHandle(T* data, int rows, int cols, int ld, int stride, int batch_size)
            : data_(data), rows_(rows), cols_(cols), ld_(ld), stride_(stride), batch_size_(batch_size), data_ptrs_(batch_size) {
                for (int i = 0; i < batch_size; i++) {
                    data_ptrs_[i] = data + i * stride;
                }
            }
        // Accessors...
        T* data_;
        SyclVector<T*> data_ptrs_;
        int rows_, cols_, ld_, stride_, batch_size_;
        Layout layout_ = Layout::ColMajor; //Most backends don't support row-major dense matrices

        private:
            void* backend_handle_;
    };

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
        SparseMatHandle(T* data, int* row_offsets, int* col_indices, int nnz, int rows, int cols, Layout layout = Layout::RowMajor)
            : data_(data), row_offsets_(row_offsets), col_indices_(col_indices), nnz_(nnz), rows_(rows), cols_(cols), layout_(layout) {
            assert(data && row_offsets && col_indices && "Null pointer provided");
            assert(nnz > 0 && rows > 0 && cols > 0 && "Invalid dimensions");
        }

        // Raw pointers to externally owned memory
        T* data_;              // [nnz] non-zero values
        int* row_offsets_;     // [rows + 1] row offsets
        int* col_indices_;     // [nnz] column indices
        int nnz_, rows_, cols_;
        Layout layout_ = Layout::RowMajor;

        private:
            void* backend_handle_;

        
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
            int nnz, int rows, int cols, int stride, int batch_size)
            :    data_(data), row_offsets_(row_offsets), col_indices_(col_indices),
                nnz_(nnz), rows_(rows), cols_(cols), stride_(stride), batch_size_(batch_size) {
            assert(data && row_offsets && col_indices && "Null pointer provided");
            assert(nnz > 0 && rows > 0 && cols > 0 && "Invalid dimensions");
            assert(stride >= nnz && "Stride must be >= nnz");
            assert(batch_size > 0 && "Batch size must be positive");
            
            for (int i = 0; i < batch_size; i++) {
                data_ptrs_[i] = data + i * stride;
                row_offsets_ptrs_[i] = row_offsets + i * (rows + 1);
                col_indices_ptrs_[i] = col_indices + i * nnz;
            }
        }

        // Raw pointers to externally owned memory
        T* data_;              // [batch_size * stride] non-zero values
        int* row_offsets_;     // [batch_size * (rows + 1)] row offsets
        int* col_indices_;     // [batch_size * nnz] column indices
        
        // Per-matrix pointers for batch processing
        SyclVector<T*> data_ptrs_;
        SyclVector<int*> row_offsets_ptrs_, col_indices_ptrs_;
        
        int nnz_, rows_, cols_, stride_, batch_size_;
        Layout layout_ = Layout::RowMajor;

        private:
            void* backend_handle_;
    };

    //Uniform accessor for data
    template <template <typename, BatchType> class Handle, typename T, BatchType BT>
    auto get_data(Handle<T,BT>& handle) {
        if constexpr (BT == BatchType::Single) {
            return handle.data_;
        } else {
            return handle.data_ptrs_.data();
        }
    }

    //Uniform accessor for data
    template <template <typename, Format, BatchType> class Handle, typename T, Format F, BatchType BT>
    auto get_data(Handle<T,F,BT>& handle) {
        if constexpr (BT == BatchType::Single) {
            return handle.data_;
        } else {
            return handle.data_ptrs_.data();
        }
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
            void* backend_handle_;
    };

    template <typename T>
    struct DenseVecHandle<T, BatchType::Batched> {
        DenseVecHandle(T* data, int size, int ldc, int stride, int batch_size) 
            : data_(data), size_(size), ldc_(ldc), stride_(stride), batch_size_(batch_size) {}
        // Accessors...
        T* data_;
        int size_, ldc_, stride_, batch_size_;

        private:
            void* backend_handle_;
    };

    template <typename T>
    struct SparseVecHandle<T, BatchType::Single> {
        SparseVecHandle(T* data, int* indices, int size) : data_(data), indices_(indices), size_(size) {}
        // Accessors...
        T* data_;
        int* indices_;
        int size_;

        private:
            void* backend_handle_;
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
            void* backend_handle_;
    };

    template <Backend B, typename T, BatchType BT>
    SyclEvent gemm(SyclQueue& ctx,
                    DenseMatHandle<T,BT>& descrA,
                    DenseMatHandle<T,BT>& descrB,
                    DenseMatHandle<T,BT>& descrC,
                    T alpha,
                    T beta,
                    Transpose transA,
                    Transpose transB,
                    ComputePrecision precision = ComputePrecision::Default);
    
    template <Backend B, typename T, BatchType BT>
    SyclEvent trsm(SyclQueue& ctx,
                    DenseMatHandle<T,BT>& descrA,
                    DenseMatHandle<T,BT>& descrB,
                    Side side,
                    Uplo uplo,
                    Transpose transA,
                    Diag diag,
                    T alpha);
    
    template <Backend B, typename T, BatchType BT>
    size_t potrf_buffer_size(SyclQueue& ctx,
                                DenseMatHandle<T,BT>& A,
                                Uplo uplo);
 
    template <Backend B, typename T, BatchType BT>
    SyclEvent potrf(SyclQueue& ctx,
                    DenseMatHandle<T,BT>& descrA,
                    Uplo uplo,
                    Span<std::byte> workspace);

}
