#include <sycl/sycl.hpp>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-device-queue.hh>
#include <iostream>
#include <iomanip>

template <typename T>
void print_matrix(Span<T> matrix, int rows, int cols, bool transpose = false) {
    std::cout << std::fixed << std::setprecision(6);
    if (!transpose) {
        std::cout << "[\n";
        for (int i = 0; i < rows; i++) {
            std::cout << "  [";
            for (int j = 0; j < cols; j++) {
                std::cout << std::setw(10) << matrix[i * cols + j];
                if (j < cols - 1) std::cout << ", ";
            }
            std::cout << "]";
            if (i < rows - 1) std::cout << ",";
            std::cout << "\n";
        }
        std::cout << "]\n";
    } else {
        std::cout << "[\n";
        for (int i = 0; i < cols; i++) { // Outer loop is cols for correct transposition
            std::cout << "  [";
            for (int j = 0; j < rows; j++) {
                std::cout << std::setw(10) << matrix[j * cols + i]; // Transposed indexing
                if (j < rows - 1) std::cout << ", ";
            }
            std::cout << "]";
            if (i < cols - 1) std::cout << ",";
            std::cout << "\n";
        }
        std::cout << "]\n";
    }   
}



template <typename T>
void chol_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace);

template <typename T>
size_t chol_qr_batched_buffer_size(SyclQueue& ctx,
                                bool transpose_op,
                                int batch_size,
                                int Astride,
                                int m,
                                int n,
                                Span<T> A);

template <typename T>
void chol2_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace);
template <typename T>
size_t chol2_qr_batched_buffer_size(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A);

template <typename T>
void shift_chol3_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace);

template <typename T>
size_t shift_chol3_qr_batched_buffer_size(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A);

template <typename T>
void svqb_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace);

template <typename T>
size_t svqb_qr_batched_buffer_size(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A);

template <typename T>
void mgs_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace);

template <typename T>
size_t mgs_qr_batched_buffer_size(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A);

template <typename T>
void blocked_mgs_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace);

template <typename T>
size_t blocked_mgs_qr_batched_buffer_size(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A);

template <typename T>
void householder_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace);

template <typename T>
size_t householder_qr_batched_buffer_size(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A);

template <typename T>
void givens_qr_batched(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A,
                        Span<std::byte> workspace);

template <typename T>
size_t givens_qr_batched_buffer_size(SyclQueue& ctx,
                        bool transpose_op,
                        int batch_size,
                        int Astride,
                        int m,
                        int n,
                        Span<T> A);
