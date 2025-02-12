#include <array>
#include <algorithm>
#include <fstream>
#include <sycl/sycl.hpp>
#include <fullerenes/sycl-headers/sycl-span.hh>
#include <oneapi/dpl/random>
#include <oneapi/dpl/iterator>
#include <oneapi/mkl/lapack.hpp>
#include <oneapi/mkl/spblas.hpp>
#include <fullerenes/graph.hh>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include "queue-impl.cc"
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/argparser.hh>
#include <fullerenes/sycl-headers/fill.hh>
#include <cusolverDn.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include "primitives.cc"
#include <iomanip>
#define CEIL_DIV(x, y) (((x) + (y) - 1) / (y))
#define DECLARE_SUBSPACE\
    auto subspace = S_acc.subspan(bid * BlockVectors*3 * m, BlockVectors * m * 3);\
    auto blockX = S_acc.subspan(bid * BlockVectors*3 * m + 0*BlockVectors*m, BlockVectors*m);\
    auto blockR = S_acc.subspan(bid * BlockVectors*3 * m + 1*BlockVectors*m, BlockVectors*m);\
    auto blockP = S_acc.subspan(bid * BlockVectors*3 * m + 2*BlockVectors*m, BlockVectors*m);\
    auto blockXtemp = S_temp_acc.subspan(bid * BlockVectors*3 * m + 0*BlockVectors*m, BlockVectors*m);\
    auto blockRtemp = S_temp_acc.subspan(bid * BlockVectors*3 * m + 1*BlockVectors*m, BlockVectors*m);\
    auto blockPtemp = S_temp_acc.subspan(bid * BlockVectors*3 * m + 2*BlockVectors*m, BlockVectors*m);\
    auto blockAX = AS_acc.subspan(bid * BlockVectors*3 * m + 0*BlockVectors*m, BlockVectors*m);\
    auto blockAR = AS_acc.subspan(bid * BlockVectors*3 * m + 1*BlockVectors*m, BlockVectors*m);\
    auto blockAP = AS_acc.subspan(bid * BlockVectors*3 * m + 2*BlockVectors*m, BlockVectors*m);

#define CHECK_CUSPARSE(func)                                                   \
{                                                                              \
    cusparseStatus_t cusparse_status = (func);                                          \
    if (cusparse_status != CUSPARSE_STATUS_SUCCESS) {                                   \
        printf("CUSPARSE API failed at line %d with error: %s (%d)\n",         \
               __LINE__, cusparseGetErrorString(cusparse_status), cusparse_status);              \
        assert(false);                                                   \
    }                                                                          \
}



template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V3_1 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V3_2 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V3_3 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V3_4 {};

using namespace sycl;

//Device Side Matmul
//A: m x k. Assume A is stored in thread local memory for now. I.E. Each thread has k*RowsPerThread elements of A
//blockX: m x n
//m is the number of rows in A and blockX
//k is the number of non-zero elements in each row of A, A is Square, Symmetric, and Positive Definite
//SN is the number of columns in the subspace i.e. the number of block vectors. (Subspace Dimension)
template <typename T, typename K, int RowsPerThread, int SN>
void matBlockVector(sycl::group<1>& cta ,const private_ptr<T> A, const private_ptr<K> cols, const Span<T> blockX, Span<T> AX, int n, int m){    
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
    group_barrier(cta);
    if(tid < m){
        for(int i = 0; i < RowsPerThread; i++){
            for(int ii = 0; ii < SN; ii++){
                T sum = 0;
                for(int j = 0; j < n; j++){
                    sum += A[i*n + j] * blockX[cols[i*n + j] + ii*m];
                }
                AX[ii*m + (tid + i*bdim)] = sum;
            }
        }

    }
}

//Projection of A onto the krylovspace spanned by S and then 
template <typename T, typename K, int RowsPerThread, int SN>
void STAS(sycl::group<1>& cta, const private_ptr<T> A, const private_ptr<K> cols, Span<T> S, int m, int n, Span<T> StAS, local_ptr<T> Scache){
//S^T * A * S
//To produce the i,j-th element of StAS, we need to compute the dot product of the i-th row of S^T and the j-th column of A*S
//The j-th column of A*S is the matrix vector product of A and the j-th column of S
//For this reason it makes sense to store the j-th column of S in shared memory and tcd .he i-th row of S^T in shared memory (the i-th column of S)

//global_ptr<T> S_global(S);
//local_ptr<T> S_j(Scache, m);
    auto tid = cta.get_local_id(0);
    for(int j = 0; j < SN; j++){
        //Load the j-th column of S into shared memory
        for (int ix = tid; ix < m; ix += cta.get_local_range(0)) Scache[ix] = S[j*m + ix];
        sycl::group_barrier(cta);
        for(int i = j; i < SN; i++){
            //Load the i-th row of S^T into registers
            T S_itid = S[i*m + tid];
            //Compute the dot product of the i-th row of S^T and the j-th column of A*S
            T AS_tidj = 0;
            for(int k = 0; k < n; k++){
                AS_tidj += A[k] * Scache[cols[k]];
            }
            //Compute the i,j-th element of StAS, the dot product of the i-th row of S^T and the j-th column of A*S
            StAS[i*SN + j] = reduce_over_group(cta, S_itid * AS_tidj, sycl::plus<T>{});
            if(tid==0) StAS[j*SN + i] = StAS[i*SN + j];
        }
        sycl::group_barrier(cta);
    }
}

template <typename T, int SN>
void inPlaceMatMatMul(sycl::group<1>& cta, Span<T> A, Span<T> B, int m, bool largest){
    auto tid = cta.get_local_id(0);
    T Alocal[SN];
    for(int j = 0; j < SN; j++){
        auto offset = SN * (largest ? (SN -j - 1) : j);
        T sum = 0;
        for(int k = 0; k < SN; k++){
            sum += A[k*m + tid] * B[offset + k];
        }
        Alocal[j] = sum;
    }
    sycl::group_barrier(cta);
    for(int j = 0; j < SN; j++){
        A[j*m + tid] = Alocal[j];
    }   
}

template <typename T, int SN, int BlockVectors>
void update_vectors(sycl::group<1>& cta, Span<T> X, Span<T> R, Span<T> P, Span<T> AX, Span<T> AR, Span<T> AP, Span<T> eigenvects, int m, bool restart, bool largest){
    auto tid = cta.get_local_id(0);
    group_barrier(cta);
    T Ps[BlockVectors], APs[BlockVectors], Xs[BlockVectors], AXs[BlockVectors];
    for(int i = 0; i < BlockVectors; i++){
        T* X_i = X.data() + i*m; T Xl = X_i[tid];
        T* R_i = R.data() + i*m; T Rl = R_i[tid];
        T* P_i = P.data() + i*m; T Pl = restart ? 0 : P_i[tid];
        T* AX_i = AX.data() + i*m; T AXl = AX_i[tid];
        T* AR_i = AR.data() + i*m; T ARl = AR_i[tid];
        T* AP_i = AP.data() + i*m; T APl = restart ? 0 : AP_i[tid];
        T* eig_X = eigenvects.data() +   (largest ? ((SN-1)*SN - i*SN ) : i*SN);
        T* eig_R = eigenvects.data() +   (largest ? ((SN-1)*SN - i*SN ) : i*SN) + BlockVectors;
        T* eig_P = eigenvects.data() +   (largest ? ((SN-1)*SN - i*SN ) : i*SN) + 2*BlockVectors;
        T tempP = 0;
        T tempAP = 0;
        T tempX = 0;
        T tempAX = 0;
        if (restart){
            for(int k = 0; k < BlockVectors; k++){
                tempP += eig_R[k] * R[k*m + tid];
                tempAP += eig_R[k] * AR[k*m + tid];
            }
        } else {
            for(int k = 0; k < BlockVectors; k++){
                tempP += eig_R[k] * R[k*m + tid] + eig_P[k] * P[k*m + tid];
                tempAP += eig_R[k] * AR[k*m + tid] + eig_P[k] * AP[k*m + tid];
        }}
        for(int k = 0; k < BlockVectors; k++){
            tempX += eig_X[k] * X[k*m + tid];
            tempAX += eig_X[k] * AX[k*m + tid];
        }

        Ps[i] = tempP;
        APs[i] = tempAP;
        Xs[i] = tempX + tempP;
        AXs[i] = tempAX + tempAP;
        //Thread 0 print temporary values
    }
    group_barrier(cta);
    for(int i = 0; i < BlockVectors; i++){
        P[i*m + tid] = Ps[i];
        AP[i*m + tid] = APs[i];
        X[i*m + tid] = Xs[i];
        AX[i*m + tid] = AXs[i];
    }
    group_barrier(cta);

}

template <typename T, int SN>
void orthonormalize(sycl::group<1>& cta, Span<T> S, int m){
    auto tid = cta.get_local_id(0);
    sycl::group_barrier(cta);
    for(int i = 0; i < SN; i++){
        T* S_i = S.data() + i*m;
        T norm = sycl::sqrt(reduce_over_group(cta, S_i[tid]*S_i[tid], sycl::plus<T>{}));
        S_i[tid] /= norm;
        sycl::group_barrier(cta);
        for(int j = i+1; j < SN; j++){
            T* S_j = S.data() + j*m;
            T projection = reduce_over_group(cta, S_i[tid]*S_j[tid], sycl::plus<T>{}) * S_i[tid];
            sycl::group_barrier(cta);
            S_j[tid] -= projection;
        }
        sycl::group_barrier(cta);
    }
}

template<typename T, int SN>
void applyConstraints(sycl::group<1>& cta, Span<T> S, const Span<T> C, int m){
    //S = S - C @ (C^T @ S)
    auto tid = cta.get_local_id(0);
    for(int i = 0; i < SN; i++){
        T* S_i = S.data() + i*m;
        T projection = 0;
        for(int j = 0; j < SN; j++){
            T* C_j = C.data() + i*m;
            projection += sycl::reduce_over_group(cta, S_i[tid]*C_j[tid], sycl::plus<T>{}) * C_j[tid];
        }
        sycl::group_barrier(cta);
        S_i[tid] -= projection;
    }
}

template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG_V3(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest, Span<T> out_eigvects, Span<T> out_eigvals){
    SyclVector<T> S(batch_size * BlockVectors*3 * m);
    SyclVector<T> S_temp(batch_size * BlockVectors*3 * m);
    SyclVector<T> AS(batch_size * BlockVectors*3 * m);
    SyclVector<T> X0(batch_size * BlockVectors * m);
    SyclVector<T> LastEigVects(batch_size * BlockVectors * BlockVectors*3);
    SyclVector<T> LastGram(batch_size * 3 * BlockVectors * 3 * BlockVectors);
    SyclVector<T> vals(batch_size * BlockVectors);
    SyclVector<T> StAS(batch_size * BlockVectors*3 * BlockVectors*3, 0);
    auto S_acc = S.to_span(); 
    auto S_temp_acc = S_temp.to_span();
    auto AS_acc = AS.to_span();
    auto X0_acc = X0.to_span();
    auto StAS_acc = StAS.to_span();

    SyclVector<int> cols_32I(cols.size());
    primitives::transform(ctx, cols, cols_32I, [](K i){return static_cast<int>(i);});
    
    size_t syevd_scratchpad_size = 0;
    size_t syevd_scratchpad_host_size = 0;
    SyclVector<int> syevd_info(batch_size);
    SyclVector<T> lambdas(batch_size * BlockVectors*3);
    SyclVector<int> csr_row_offsets(m + 1);
    
    //csr_row_offsets = [0, NZ, 2*NZ, 3*NZ, ...]
    primitives::iota(ctx, csr_row_offsets, 0);
    primitives::transform(ctx, csr_row_offsets, csr_row_offsets, [&](int i){return i * NZ;}); 

    auto lambda_acc = lambdas.to_span();

    cusolverDnHandle_t handle{};
    cudaStream_t stream{};
    cudaDeviceProp props{};
    cusolverDnParams_t params{};
    cusolverDnCreate(&handle);
    cudaStreamCreate(&stream);
    cusolverDnSetStream(handle, stream);
    cusolverDnCreateParams(&params);
    /* LOBPCG_Helper<T, K> helper(stream, BlockVectors,
                                m, NZ, batch_size, A.data(), cols.data(), csr_row_offsets.data(),
                                S.data(), StAS.data(), S.data() + BlockVectors*9*m, lambdas.data()); */
     
    #if (SOLVER_BACKEND == CUSOLVER_BACKEND)
        T alpha = 1.0;
        T beta = 0.0;
        cusolverStatus_t status = cusolverDnXsyevBatched_bufferSize(handle, 
                                params, 
                                CUSOLVER_EIG_MODE_VECTOR, 
                                CUBLAS_FILL_MODE_LOWER,
                                m,
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                A.data(),
                                m,
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                lambdas.data(),
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                &syevd_scratchpad_size,
                                &syevd_scratchpad_host_size,
                                batch_size);

        cusparseHandle_t SpMM_AS_handle;
        CHECK_CUSPARSE(cusparseCreate(&SpMM_AS_handle))
        CHECK_CUSPARSE(cusparseSetStream(SpMM_AS_handle, stream))
        cusparseSpMatDescr_t descrA;    //Left matrix
        cusparseDnMatDescr_t descrX;    //Right matrix
        cusparseDnMatDescr_t descrXR;   //Concatination matrix [X, R]
        cusparseDnMatDescr_t descrS;    //The full subspace matrix [X, R, P]
        cusparseDnMatDescr_t descrAX;  //The result of the matrix matrix multiplication
        cusparseDnMatDescr_t descrAXR; //The result of the matrix matrix multiplication
        cusparseDnMatDescr_t descrAS;  //The result of the matrix matrix multiplication
        CHECK_CUSPARSE(cusparseCreateCsr(  &descrA,    m,  m, m*NZ, csr_row_offsets.data(), cols_32I.data(), A.data(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrX,    m,  BlockVectors,   m,  S.data(),   CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrXR,   m,  BlockVectors*2, m,  S.data(),   CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrS,    m,  BlockVectors*3, m,  S.data(),   CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrAX,   m,  BlockVectors,   m,  AS.data(),  CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrAXR,  m,  BlockVectors*2, m,  AS.data(),   CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCreateDnMat(&descrAS,   m,  BlockVectors*3, m,  AS.data(),   CUDA_R_32F, CUSPARSE_ORDER_COL))
        CHECK_CUSPARSE(cusparseCsrSetStridedBatch(     descrA, batch_size, m*NZ, m*NZ))
        auto SubspaceStride = m*BlockVectors*3; //[X_0, R_0, P_0, X_1, R_1, P_1, ...]

        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrX, batch_size,     SubspaceStride))
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrXR, batch_size,    SubspaceStride))
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrS, batch_size,     SubspaceStride))
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrAX, batch_size,    SubspaceStride)) //AX = A * X
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrAXR, batch_size,   SubspaceStride)) //AXR = A * [X, R]
        CHECK_CUSPARSE(cusparseDnMatSetStridedBatch(   descrAS, batch_size,    SubspaceStride)) //AS = A * [X, R, P]

        size_t SpMM_buffer_size = 0;
        CHECK_CUSPARSE(cusparseSpMM_bufferSize(SpMM_AS_handle,
                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha,
                                descrA,
                                descrS,
                                &beta,
                                descrAS,
                                CUDA_R_32F,
                                CUSPARSE_SPMM_ALG_DEFAULT,
                                &SpMM_buffer_size))
        
        SyclVector<T> SpMM_buffer(SpMM_buffer_size);

        cublasHandle_t GEMM_STAS_handle;
        cublasCreate(&GEMM_STAS_handle);
        cublasSetStream(GEMM_STAS_handle, stream);
    

    #elif defined(USE_ROCSOLVER)
        status_t status = rocsolver_ssyevd_bufferSize(handle, params, EIG_SOLVE_MODE, FILL_MODE_LOWER, m, A.data(), m, lambdas.data(), &syevd_scratchpad_size, &syevd_scratchpad_host_size, batch_size);
    #endif

    SyclVector<T> syevd_scratchpad (syevd_scratchpad_size);
    SyclVector<T> syevd_scratchpad_host (syevd_scratchpad_host_size);

    #if (SOLVER_BACKEND == CUSOLVER_BACKEND)
        auto compute_eigenpairs = [&](T* A, T* lambdas, int m, int batch_size){
            cusolverStatus_t status = cusolverDnXsyevBatched(handle, 
                                params, 
                                CUSOLVER_EIG_MODE_VECTOR, 
                                CUBLAS_FILL_MODE_LOWER,
                                m,
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                A,
                                m,
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                lambdas,
                                std::is_same_v<T, float> ? CUDA_R_32F : CUDA_R_64F,
                                syevd_scratchpad.data(),
                                syevd_scratchpad_size,
                                syevd_scratchpad_host.data(),
                                syevd_scratchpad_host_size,
                                syevd_info.data(),
                                batch_size);
            if(status != CUSOLVER_STATUS_SUCCESS){
                std::cerr << "Error in cusolverDnXsyevBatched: " << status << std::endl;
            }
        };
    #elif defined(USE_ROCSOLVER)
        auto compute_eigenpairs = [&](T* A, T* lambdas, T* eigenvectors, int m, int batch_size){
            status_t status = rocsolver_ssyevd(handle, params, EIG_SOLVE_MODE, FILL_MODE_LOWER, m, A, m, lambdas, eigenvectors, m, syevd_scratchpad, syevd_scratchpad_size, syevd_scratchpad_host_size, syevd_info.data(), batch_size);
        };
    #endif
    const T tol = sycl::sqrt(std::numeric_limits<T>::epsilon()) * T(m);
    auto print_matrix = [&](Span<T> matrix, int rows, int cols, bool transpose = false) {
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
};
    //Setup
    ctx -> submit([&](sycl::handler& h ){
        auto Scache = sycl::local_accessor<T, 1>(m, h);

        h.parallel_for<LOBPCG_V3_1<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            sycl::group<1> cta = item.get_group();
            auto bid = item.get_group_linear_id();
            std::bitset<BlockVectors> converged;
            oneapi::dpl::uniform_real_distribution<T> distr(0.0, 1.0);            
            oneapi::dpl::minstd_rand engine(42, tid);
            auto A_acc = A.subspan(bid * m * NZ, m * NZ);
            auto cols_acc = cols.subspan(bid * m * NZ, m * NZ);

            //Load the i-th row of A into registers
            T A_tid[NZ];
            for(int i = 0; i < NZ; i++){
                A_tid[i] = A_acc[tid*NZ + i];
            }
            //Load the i-th row of cols into registers
            K cols_tid[NZ]; 
            for(int i = 0; i < NZ; i++){
                cols_tid[i] = cols_acc[tid*NZ + i];
            }
            DECLARE_SUBSPACE;

            for (int i = tid; i < 3*m*BlockVectors; i+=cta.get_local_range(0)){
                S_acc[i] = T(0);
                S_temp_acc[i] = T(0);
                AS_acc[i] = T(0);
            }
            sycl::group_barrier(cta);
            for (int i = tid; i < m*BlockVectors; i+=cta.get_local_range(0)){
                blockX[i] = distr(engine);
                //X0_acc[i] = blockX[i];
                X0_acc[bid*m*BlockVectors + i] = blockX[i];
            }
            group_barrier(cta);
            //Normalize S vectors
            orthonormalize<T, BlockVectors>(cta, blockX, m); //Modified Gram-Schmidt
            
        });
    });
    ctx.wait();

    //Compute AX
    auto SpMM_status = cusparseSpMM(SpMM_AS_handle, 
                                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    &alpha,
                                    descrA,
                                    descrX,
                                    &beta,
                                    descrAX,
                                    CUDA_R_32F,
                                    CUSPARSE_SPMM_ALG_DEFAULT,
                                    SpMM_buffer.data());

    cudaStreamSynchronize(stream);
    
    //Compute X^T AX
    auto gemm_status = cublasSgemmStridedBatched(GEMM_STAS_handle,
                                                CUBLAS_OP_T,
                                                CUBLAS_OP_N,
                                                BlockVectors, BlockVectors, m,
                                                &alpha,
                                                S.data(), 
                                                m, m*3*BlockVectors,
                                                AS.data(),
                                                m, m*3*BlockVectors,
                                                &beta,
                                                StAS.data(),
                                                BlockVectors, BlockVectors*BlockVectors,
                                                batch_size);

    compute_eigenpairs(StAS.data(), lambdas.data(), BlockVectors, batch_size);
    cudaStreamSynchronize(stream);

    ctx -> submit([&](sycl::handler& h){
        auto Scache = sycl::local_accessor<T, 1>(m, h);

        h.parallel_for<LOBPCG_V3_2<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            sycl::group<1> cta = item.get_group();
            auto bid = item.get_group_linear_id();

            auto A_acc = A.subspan(bid * m * NZ, m * NZ);
            auto cols_acc = cols.subspan(bid * m * NZ, m * NZ);

            //Load the i-th row of A into registers
            T A_tid[NZ];
            for(int i = 0; i < NZ; i++){
                A_tid[i] = A_acc[tid*NZ + i];
            }
            //Load the i-th row of cols into registers
            K cols_tid[NZ]; 
            for(int i = 0; i < NZ; i++){
                cols_tid[i] = cols_acc[tid*NZ + i];
            }


            //X^T @ X  = (X^T @ X)^T = X @ X^T
            //iff A^T = A then 
            //X^T @ A @ X = (X^T @ A @ X)^T

            DECLARE_SUBSPACE;

            matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockX, blockAX, NZ, m);
            //cuSolver stores the eigenvectors in-place in the input matrix
            auto eigvects = StAS_acc.subspan(bid * BlockVectors * BlockVectors, BlockVectors * BlockVectors);
            inPlaceMatMatMul<T, BlockVectors>(cta, blockX, eigvects, m , largest);
            inPlaceMatMatMul<T, BlockVectors>(cta, blockAX, eigvects, m, largest);
            //Compute the residual R = A*X - X*Lambda
        });
    });
    ctx->wait();

    SyclVector<std::bitset<BlockVectors>> converged(batch_size);
    SyclVector<T> residuals(batch_size * BlockVectors);
    auto converged_acc = converged.to_span();
    auto residuals_acc = residuals.to_span();

    auto iteration = [&](SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t iter, bool restart){
        
        //Setup
        ctx -> submit([&](sycl::handler& h){
            auto Scache = sycl::local_accessor<T, 1>(m, h);

            h.parallel_for<LOBPCG_V3_3<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
                auto tid = item.get_local_linear_id();
                auto cta = item.get_group();
                auto bid = item.get_group_linear_id();
                auto& converged = converged_acc[bid];
                auto residuals = residuals_acc.subspan(bid * BlockVectors, BlockVectors);
                auto num_eigvals = iter < 2 ? (iter+1) * BlockVectors : 3*BlockVectors;
                auto lambdas = lambda_acc.subspan(bid * num_eigvals, num_eigvals);


                auto A_acc = A.subspan(bid * m * NZ, m * NZ);
                auto cols_acc = cols.subspan(bid * m * NZ, m * NZ);

                //Load the i-th row of A into registers
                T A_tid[NZ];
                for(int i = 0; i < NZ; i++){
                    A_tid[i] = A_acc[tid*NZ + i];
                }
                //Load the i-th row of cols into registers
                K cols_tid[NZ]; 
                for(int i = 0; i < NZ; i++){
                    cols_tid[i] = cols_acc[tid*NZ + i];
                }
            
                DECLARE_SUBSPACE;

                //R = A*X - X*Lambda
                for(int i = 0; i < BlockVectors; i++) blockR[i*m + tid] = blockAX[i*m + tid] - lambdas[largest ? (num_eigvals - 1 - i) : i] * blockX[i*m + tid];
                //Convergence Check
                for(int i = 0; i < BlockVectors; i++) {
                    //if(converged[i]) continue; 
                    residuals[i] = sycl::sqrt(reduce_over_group(cta, blockR[i*m + tid]*blockR[i*m + tid], sycl::plus<T>{}));
                    residuals[i] /= sycl::sqrt(reduce_over_group(cta, blockX[i*m + tid]*blockX[i*m + tid], sycl::plus<T>{})) * lambdas[largest ? (num_eigvals - 1 - i) : i];
                    converged[i] = residuals[i] < tol;
                }
                sycl::group_barrier(cta);
                sycl::group_barrier(cta);
                auto StAS_offset = restart ? bid * BlockVectors*2 * BlockVectors*2 : bid * BlockVectors*3 * BlockVectors*3;
                auto StAS_size = restart ? BlockVectors*2 * BlockVectors*2 : BlockVectors*3 * BlockVectors*3;
                auto StAS = StAS_acc.subspan(StAS_offset, StAS_size);

                auto AS_block = AS_acc.subspan(bid * m * BlockVectors*3, m * BlockVectors*3);
                auto S_block = S_acc.subspan(bid * m * BlockVectors*3, m * BlockVectors*3);

                if (!restart) {
                    orthonormalize<T, BlockVectors*3>(cta, blockX, m);
                } else {
                    orthonormalize<T, BlockVectors*2>(cta, blockX, m);
                }
            });
        });
        ctx -> wait();

        //Compute AX
        auto SpMM_status = cusparseSpMM(SpMM_AS_handle, 
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        &alpha,
                                        descrA,
                                        restart ? descrXR : descrS,
                                        &beta,
                                        restart ? descrAXR : descrAS,
                                        CUDA_R_32F,
                                        CUSPARSE_SPMM_ALG_DEFAULT,
                                        SpMM_buffer.data());

        cudaStreamSynchronize(stream);

        //Compute X^T AX
        auto gemm_status = cublasSgemmStridedBatched(GEMM_STAS_handle,
                                                    CUBLAS_OP_T,
                                                    CUBLAS_OP_N,
                                                    BlockVectors * (restart ? 2 : 3),
                                                    BlockVectors * (restart ? 2 : 3),
                                                    m,
                                                    &alpha,
                                                    S.data(), 
                                                    m, m*3*BlockVectors,
                                                    AS.data(),
                                                    m, m*3*BlockVectors,
                                                    &beta,
                                                    StAS.data(),
                                                    BlockVectors * (restart ? 2: 3), 
                                                    BlockVectors * (restart ? 2*2: 3*3) * BlockVectors,
                                                    batch_size); 

        //Solve the eigenvalue problem
        cudaStreamSynchronize(stream);
        compute_eigenpairs(StAS.data(), lambdas.data(), restart ? BlockVectors*2 : BlockVectors*3, batch_size);
        cudaStreamSynchronize(stream);

        ctx -> submit([&](sycl::handler& h){
            auto Scache = sycl::local_accessor<T, 1>(m, h);
            h.parallel_for<LOBPCG_V3_4<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
                sycl::group<1> cta = item.get_group();
                auto bid = item.get_group_linear_id();
            
                DECLARE_SUBSPACE;
                auto eigvects_offset = restart ? bid * BlockVectors * 2 * BlockVectors*2 : bid * BlockVectors * 3 * BlockVectors*3;
                auto eigvects_size = restart ? BlockVectors * 2 * BlockVectors*2 : BlockVectors * 3 * BlockVectors*3;
                auto eigvects = StAS_acc.subspan(eigvects_offset, eigvects_size);

                if (!restart) {
                    update_vectors<T, BlockVectors*3, BlockVectors>(cta, blockX, blockR, blockP, blockAX, blockAR, blockAP, eigvects, m, false, largest);
                } else {
                    update_vectors<T, BlockVectors*2, BlockVectors>(cta, blockX, blockR, blockP, blockAX, blockAR, blockAP, eigvects, m, true, largest);
                }
            });
        });
    };
    auto T0 = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < maxiters; i++){
        iteration(ctx, A, cols, batch_size, m, i, i < 1 ? true : false);
    }
    ctx.wait();   
    auto T1 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(T1 - T0).count();
    std::cout << "LOBPCG V3 Kernel Time: " << duration / batch_size << " µs / fullerene" << std::endl;

    //std::cout << "Eigenvalues: " << lambdas << std::endl;
    //std::cout << "Residuals: " << residuals << std::endl;

    for(int i = 0; i < batch_size; i++){
        primitives::copy(ctx, lambdas.subspan(i * (BlockVectors *3) + BlockVectors * 2, BlockVectors), out_eigvals.subspan(i * BlockVectors, BlockVectors));
        //primitives::copy(ctx, S_acc.subspan(i * BlockVectors * 3 * m, BlockVectors * m), eigvects.subspan(i * BlockVectors * m, BlockVectors * m));
    }
}

template void LOBPCG_V3<float, uint16_t, 3, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest, Span<float> eigvects, Span<float> eigvals);
template void LOBPCG_V3<float, uint16_t, 21, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest, Span<float> eigvects, Span<float> eigvals);