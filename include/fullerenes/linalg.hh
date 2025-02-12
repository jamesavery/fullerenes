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
        SYCL
        // Add more as needed
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



    template <typename T, BatchType BT>
    struct MatHandle;


    template <typename T, Format F, BatchType BT>
    struct SparseMatHandle;


    template <typename T>
    struct MatHandle<T, BatchType::Single> {
        MatHandle(T* data, int rows, int cols, int ld) 
            : data_(data), rows_(rows), cols_(cols), ld_(ld) {}
        // Accessors...
        T* data_;
        int rows_, cols_, ld_;
    };

    template <typename T>
    struct MatHandle<T, BatchType::Batched> {
        MatHandle(T* data, int rows, int cols, int ld, int stride, int batch_size)
            : data_(data), rows_(rows), cols_(cols), ld_(ld), stride_(stride), batch_size_(batch_size), data_ptrs_(batch_size) {
                for (int i = 0; i < batch_size; i++) {
                    data_ptrs_[i] = data + i * stride;
                }
            }
        // Accessors...
        T* data_;
        SyclVector<T*> data_ptrs_;
        int rows_, cols_, ld_, stride_, batch_size_;
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
         */
        SparseMatHandle(T* data, int* row_offsets, int* col_indices, int nnz, int rows, int cols)
            : data_(data), row_offsets_(row_offsets), col_indices_(col_indices), nnz_(nnz), rows_(rows), cols_(cols) {
            assert(data && row_offsets && col_indices && "Null pointer provided");
            assert(nnz > 0 && rows > 0 && cols > 0 && "Invalid dimensions");
        }

        // Raw pointers to externally owned memory
        T* data_;              // [nnz] non-zero values
        int* row_offsets_;     // [rows + 1] row offsets
        int* col_indices_;     // [nnz] column indices
        int nnz_, rows_, cols_;
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
            : data_(data), row_offsets_(row_offsets), col_indices_(col_indices),
              nnz_(nnz), rows_(rows), cols_(cols), stride_(stride), batch_size_(batch_size),
              data_ptrs_(batch_size), row_offsets_ptrs_(batch_size), col_indices_ptrs_(batch_size) {
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
    struct VecHandle;

    template <typename T>
    struct VecHandle<T, BatchType::Single> {
        VecHandle(T* data, int size, int ldc) : data_(data), size_(size), ldc_(ldc) {}
        // Accessors...
        T* data_;
        SyclVector<T*> data_ptrs_;
        int size_, ldc_;
    };

    template <typename T>
    struct VecHandle<T, BatchType::Batched> {
        VecHandle(T* data, int size, int ldc, int stride, int batch_size) 
            : data_(data), size_(size), ldc_(ldc), stride_(stride), batch_size_(batch_size) {}
        // Accessors...
        T* data_;
        int size_, ldc_, stride_, batch_size_;
    };

    template <Backend B, typename T, BatchType BT>
    SyclEvent gemm(SyclQueue& ctx,
                    MatHandle<T,BT>& descrA,
                    MatHandle<T,BT>& descrB,
                    MatHandle<T,BT>& descrC,
                    T alpha,
                    T beta,
                    Transpose transA,
                    Transpose transB,
                    ComputePrecision precision = ComputePrecision::Default);
    
    template <Backend B, typename T, BatchType BT>
    SyclEvent trsm(SyclQueue& ctx,
                    MatHandle<T,BT>& descrA,
                    MatHandle<T,BT>& descrB,
                    Side side,
                    Uplo uplo,
                    Transpose transA,
                    Diag diag,
                    T alpha);
    
    template <Backend B, typename T, BatchType BT>
    size_t potrf_buffer_size(SyclQueue& ctx,
                                MatHandle<T,BT>& A,
                                Uplo uplo);

    template <Backend B, typename T, BatchType BT>
    SyclEvent potrf(SyclQueue& ctx,
                    MatHandle<T,BT>& descrA,
                    Uplo uplo,
                    Span<std::byte> workspace);

}
