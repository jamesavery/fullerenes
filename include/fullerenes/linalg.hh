#pragma once
#include <complex>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>

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



    template <typename T, BatchType BT>
    struct MatHandle;

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
            : data_(data), rows_(rows), cols_(cols), ld_(ld), stride_(stride), batch_size_(batch_size) {}
        // Accessors...
        T* data_;
        int rows_, cols_, ld_, stride_, batch_size_;
    };


    template <typename T, BatchType BT>
    struct VecHandle;

    template <typename T>
    struct VecHandle<T, BatchType::Single> {
        VecHandle(T* data, int size, int ldc) : data_(data), size_(size), ldc_(ldc) {}
        // Accessors...
        T* data_;
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
                    Transpose transB);
    
    template <Backend B, typename T, BatchType BT>
    SyclEvent trsm(SyclQueue& ctx,
                    Side side,
                    Uplo uplo,
                    Diag diag,
                    MatHandle<T,BT>& descrA,
                    MatHandle<T,BT>& descrB,
                    T alpha);

}
