#include <array>
#include <algorithm>
#include <fstream>
#include "sycl/sycl.hpp"
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
#include <fullerenes/sycl-headers/sycl-mdspan.hh>
#define CEIL_DIV(x, y) (((x) + (y) - 1) / (y))

#if defined(SOLVER_BACKEND) && (SOLVER_BACKEND == CUSOLVER_BACKEND)
    #define USE_CUSOLVER
    #include <cusolverDn.h>
    #include <cublas_v2.h>
    #include <cusparse.h>
    #include <cuda_runtime.h>
    #include <cuda_runtime_api.h>
    #define handle_t cusolverDnHandle_t
    #define status_t cusolverStatus_t
    #define stream_t cudaStream_t
    #define params_t cusolverDnParams_t
    #define props_t cudaDeviceProp
    #define synchronize cudaDeviceSynchronize
    #define sync_stream cudaStreamSynchronize
    #define create_handle cusolverDnCreate
    #define destroy_handle cusolverDnDestroy
    #define set_stream cusolverDnSetStream
    #define create_stream cudaStreamCreate
    #define destroy_stream cudaStreamDestroy
    #define get_device_properties cudaGetDeviceProperties
    #define create_params cusolverDnCreateParams
    #define EIG_SOLVE_MODE CUSOLVER_EIG_MODE_VECTOR
    #define FILL_MODE_LOWER CUBLAS_FILL_MODE_LOWER
    #define SsyevBatched cusolverDnXsyevBatched
    #define SsyevBatched_bufferSize cusolverDnXsyevBatched_bufferSize
    #define SOLVER_COMPUTE_TYPE CUDA_R_32F
#elif defined(SOLVER_BACKEND) && (SOLVER_BACKEND == ROCSOLVER_BACKEND)
    #define USE_ROCSOLVER
    #include <rocsolver.h>
    #define handle_t rocsolver_handle
    #define status_t rocsolver_status
    #define stream_t hipStream_t
    #define synchronize hipDeviceSynchronize
    #define create_handle rocsolver_create_handle
    #define destroy_handle rocsolver_destroy_handle
    #define set_stream rocsolver_set_stream
    #define create_stream hipStreamCreate
    #define destroy_stream hipStreamDestroy
#endif

#define DECLARE_SUBSPACE\
    auto subspace = S_acc.subspan(bid * BlockVectors*9 * m, BlockVectors * m * 9);\
    auto blockX = S_acc.subspan(bid * BlockVectors*9 * m, BlockVectors*m);\
    auto blockR = S_acc.subspan(bid * BlockVectors*9 * m + BlockVectors*m, BlockVectors*m);\
    auto blockP = S_acc.subspan(bid * BlockVectors*9 * m + 2*BlockVectors*m, BlockVectors*m);\
    auto blockXtemp = S_acc.subspan(bid * BlockVectors*9 * m + 3*BlockVectors*m, BlockVectors*m);\
    auto blockRtemp = S_acc.subspan(bid * BlockVectors*9 * m + 4*BlockVectors*m, BlockVectors*m);\
    auto blockPtemp = S_acc.subspan(bid * BlockVectors*9 * m + 5*BlockVectors*m, BlockVectors*m);\
    auto blockAX = S_acc.subspan(bid * BlockVectors*9 * m + 6*BlockVectors*m, BlockVectors*m);\
    auto blockAR = S_acc.subspan(bid * BlockVectors*9 * m + 7*BlockVectors*m, BlockVectors*m);\
    auto blockAP = S_acc.subspan(bid * BlockVectors*9 * m + 8*BlockVectors*m, BlockVectors*m);

template <typename T, typename K, int BlockVectors, int NZ> class LOBPCGg {};

template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V1_1 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V1_2 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V1_3 {};
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCG_V1_4 {};



using namespace sycl;
//LOBPCG Algorithm implemented in SYCL
//Version 1.0
//Naive Approach: All matrices, and block vectors are stored in Global Memory

template <typename T>
void apply_all_reflections(const sycl::group<1> &cta, T* V, const int n, const int m, T* Q)
{   
    static_assert(std::is_floating_point<T>::value, "T must be floating point");
    auto tid = cta.get_local_linear_id();
    auto bdim = cta.get_local_linear_range();
    for(int k=0;k<n;k++){
        const T &v0 = V[2*k], &v1 = V[2*k+1];      
        // Udrullet:
        //       apply_reflection(Q({k,k+2},{0,m}), v);
        for(int l=tid;l<m; l+=bdim){
            T &q0 = Q[k*m+l], &q1 = Q[(k+1)*m+l];
            T vTA = q0*v0 + q1*v1;
            q0 -= 2*v0*vTA;
            q1 -= 2*v1*vTA;
        }      
    }  
}
template <typename T>
void T_QTQ(sycl::group<1>& cta, const int n, T* D, T* L, T* U, T* Vout, T shift=0){
    int tix = cta.get_local_linear_id();
    int bdim = cta.get_local_linear_range();
    FLOAT_TYPEDEFS(T);
    //  QTQ_calls ++;
    // Unrolled
    //  real_t numerical_zero = T.max_norm()*10*std::numeric_limits<real_t>::epsilon();
    // specialized max_norm = max(sum(abs(A),axis=1)) for tridiagonal matrix. 
    real_t local_max = real_t(0.);
    for (int i = tix; i < n; i += bdim){
        local_max = std::max(local_max, std::abs(D[i]) + 2*std::abs(L[i]));
    }
    real_t max_norm = reduce_over_group(cta, local_max, sycl::maximum<real_t>());
    real_t numerical_zero = 10*std::numeric_limits<real_t>::epsilon();
    real_t d_n, l_n, l_nm1;
    d_n = D[n]; l_n = L[n]; l_nm1 = L[n-1];
    sycl::group_barrier(cta);
    //real_t a[2], v[2], D[n+1], L[n+1], U[2*(n+1)];
    real_t a[2], v[2];//, D[n+1], L[n+1], U[2*(n+1)];
    for(int k = tix; k < n + 1; k += bdim){
        D[k] -= shift;
        U[n+1 + k] = real_t(0.);
        if(k < n-1){
            U[k] = L[k];
            Vout[2*k] = real_t(0.); Vout[2*k+1] = real_t(0.);
        } else {
            L[k] = real_t(0.);
            U[k] = real_t(0.);
        }
    }

    sycl::group_barrier(cta);
    if(tix == 0)
        for(int k=0;k<n-1;k++){
            if (std::abs(L[k]) > numerical_zero){
            a[0] = D[k]; a[1] = L[k];       // a = T[k:k+2,k] is the vector of nonzeros in kth subdiagonal column.
            
            real_t anorm = sycl::sqrt(a[0]*a[0] + a[1]*a[1]); 

            // // Udrullet
            // //    reflection_vector(a,anorm,v);
            v[0] = D[k]; v[1] = L[k];
            real_t alpha = -sycl::copysign(anorm,a[0]); // Koster ingenting
            v[0] -= alpha;

            real_t vnorm = sycl::sqrt(v[0]*v[0]+v[1]*v[1]);
            real_t norm_inv = real_t(1.)/vnorm;               //Normalize
            v[0] *= norm_inv;  v[1] *= norm_inv;

            Vout[2*k] = v[0]; Vout[2*k+1] = v[1];
            
            // // Udrullet 
            // //    apply_reflection(T({k,k+2},{k,k+3}),v);
            // //      if(k+1<n){			// k=n-1 case handled by padding with zeros
            coord3d vTA = { D[ k ]*v[0] + L[ k ]*v[1],  // T(k,k  )*v[0] + T(k+1,k  )*v[1]
                            U[ k ]*v[0] + D[k+1]*v[1],  // T(k,k+1)*v[0] + T(k+1,k+1)*v[1]
                            U[(n+1)+k]*v[0] + U[k+1]*v[1]}; // T(k,k+2)*v[0] + T(k+1,k+2)*v[1]

        
            D[k]     -= real_t(2.)*v[0]*vTA[0];
            L[k]     -= real_t(2.)*v[1]*vTA[0];
            U[k]     -= real_t(2.)*v[0]*vTA[1];
            D[k+1]     -= real_t(2.)*v[1]*vTA[1];
            U[(n+1)+k] -= real_t(2.)*v[0]*vTA[2];
            U[k+1]     -= real_t(2.)*v[1]*vTA[2];
            }
        }

    if(tix == 0)
    { // Transform from the right = transform columns of the transpose.
        int k = 0;
        const real_t *v = &Vout[0];
        real_t   vTA[2] = {D[ k ]*v[0] + U[  k  ]*v[1],  // T(k,k  )*v[0] + T(k,  k+1)*v[1]
                    0        + D[ k+1 ]*v[1]}; // T(k+1,k)*v[0] + T(k+1,k+1)*v[1]. Lower subdiagonal is zero at this stage.
        
        D[k]       -= real_t(2.)*v[0]*vTA[0]; // T(k,k)     -= 2*v[0]*vTA[0]
        U[k]       -= real_t(2.)*v[1]*vTA[0]; // T(k,k+1)   -= 2*v[1]*vTA[0]
        L[k]       -= real_t(2.)*v[0]*vTA[1]; // T(k+1,k)   -= 2*v[0]*vTA[1]
        D[k+1]     -= real_t(2.)*v[1]*vTA[1]; // T(k+1,k+1) -= 2*v[1]*vTA[1]        
    }
    sycl::group_barrier(cta);


    if(tix == 0){
        for(int k=1;k<n-1;k++){
            const real_t *v = &Vout[2*k];
            coord3d vTA = {U[k-1]*v[0] + U[(n+1)+k-1]*v[1], // T(k-1,k)*v[0] + T(k-1,k+1)*v[1]  
                            D[ k ]*v[0] + U[  k  ]*v[1],     // T(k,k  )*v[0] + T(k,  k+1)*v[1]
                            L[ k ]*v[0] + D[ k+1 ]*v[1]};    // T(k+1,k)*v[0] + T(k+1,k+1)*v[1]. Lower subdiagonal is zero at this stage

            U[k-1]     -= real_t(2.)*v[0]*vTA[0];     // T(k-1,k)   -= 2*v[0]*vTA[0]
            U[(n+1)+(k-1)] -= real_t(2.)*v[1]*vTA[0]; // T(k-1,k+1) -= 2*v[1]*vTA[0]
            U[k]       -= real_t(2.)*v[1]*vTA[1];     // T(k,  k+1)   -= 2*v[1]*vTA[1]
            D[k]       -= real_t(2.)*v[0]*vTA[1];     // T(k,  k)     -= 2*v[0]*vTA[1]
            L[k]       -= real_t(2.)*v[0]*vTA[2];     // T(k+1,k)   -= 2*v[0]*vTA[2]
            D[k+1]     -= real_t(2.)*v[1]*vTA[2];     // T(k+1,k+1) -= 2*v[1]*vTA[2]        
        }
    }


    sycl::group_barrier(cta);
    for (int k = tix; k<n; k+=bdim){  // Copy working diagonals to output
        D[k] += shift;
        if(k < n-1){
            L[k] = U[k];
        }
    }
    sycl::group_barrier(cta);
    if (tix==0){
        D[n] = d_n;
        L[n-1] = l_nm1;
        L[n] = l_n;
    }
    sycl::group_barrier(cta);

}
template <typename T>
std::array<T,2> eigvalsh2x2(const std::array<T,4> &A){
    auto [a,b,c,d] = A;
    T D = sycl::sqrt(4*b*c+(a-d)*(a-d));
    return {(a+d-D)/2, (a+d+D)/2};
}
//Assumes A is symmetric
template <typename T, int N> 
void lanczos(sycl::group<1>& cta, Span<T> A, Span<T> X, Span<T> alphas, Span<T> betas, Span<T> Vspan){
    auto tid = cta.get_local_linear_id();
    T* V = Vspan.data() + tid;
    oneapi::dpl::uniform_real_distribution<T> distr(0.0, 1.0);            
    oneapi::dpl::minstd_rand engine(42, tid);

    //Assumes A is column major, gives the best access pattern
    auto mat_vect = [&](const T x){
        T result = T(0);
        if(tid < N) X[tid] = x;
        sycl::group_barrier(cta);
        if (tid < N){
        #pragma unroll
        for (int j = 0; j < N; j++){
            result += A[j*N + tid] * X[j];
        }
        }
        sycl::group_barrier(cta);
        return result;
    };

    //Modified Gram-Schmidt
    auto MGS = [&](int index){
        sycl::group_barrier(cta);
        T result = (tid < N) ?  V[index*N] : T(0);

        #pragma unroll
        for (int j = 0; j < index; j++){
            auto other = (tid < N) ? V[j*N] : T(0);
            auto proj = reduce_over_group(cta, result * other, sycl::plus<T>{}) * other;
            result -= (tid < N) ? proj : T(0);
        }
        result /= sycl::sqrt(reduce_over_group(cta, result * result, sycl::plus<T>{}));
        return result;
    };

    auto v0 = (tid < N) ? distr(engine) : T(0);
    auto norm = sycl::sqrt(reduce_over_group(cta, v0 * v0, sycl::plus<T>{}));
    if (tid < N)  V[0*N] = v0 / norm;
    v0 = MGS(0);
    if (tid < N) V[0*N] = v0;
    for (int i = 0; i < N; i++){
        if (i > 0){
            T v_i = MGS(i);
            if(tid < N) V[i*N] = v_i;
        }
        T vect_input = (tid < N) ? V[i*N] : T(0);
        T v = mat_vect(vect_input);
        T alpha = reduce_over_group(cta, v * vect_input,sycl::plus<T>{});
        if (tid == i) alphas[i] = alpha;
        if (i == 0){
            v -= alpha * vect_input;
        } else {
            T v_minus_one = (tid < N) ? V[(i-1)*N] : T(0);
            v -= betas[i- 1] * v_minus_one + alpha * vect_input;
        }
        v = (tid < N) ? v : T(0);
        T beta = sycl::sqrt(reduce_over_group(cta, v * v, sycl::plus<T>{}));
        if (tid == i) betas[i] = beta;
        if ((i < N-1) && (tid  < N)) V[(i+1)*N] = v / beta;
    } 
}
/* 
template <typename real_t>
void reflect_region(group<1>& cta,
    local_accessor<real_t,1> &A,int N,           // Matrix top transform
    int i0, int j0, int m, int n,                // Region of A to transform
    const real_t &v_i,                           // Thread's element of the reflection vector v
    int cols,                                    // Left- or right-reflection (COLS or ROWS)
    local_accessor<real_t,1>& vHA)               // Length-n working memory for v^H A[region] 
{
    // TODO: Implement segmented sum to get full m*n parallelism.
    int stride[2] = {N * (!cols) + 1 * cols, N * cols + 1 * (!cols)};

    int i_tid = cta.get_local_id(0);
    for (int j = 0; j < n; j++) 
    {
        size_t IJ = (i_tid + i0) * stride[0] + (j + j0) * stride[1];
        real_t vAij = v_i*A[IJ];
        real_t sum = reduce_over_group(cta, vAij, plus<real_t>());
        vHA[j] = sum;
    }

    for (size_t j = 0; j < n; j++) // A += -2*outer(v,vTA)
    {
        size_t IJ = (i_tid + i0) * stride[0] + (j + j0) * stride[1];
        A[IJ] -= 2 * v_i * vHA[j];
    }
}

template <typename real_t>
real_t max_norm(const group<1>& cta, 
               const local_accessor<real_t,1>& A, const int m, const int n)
{//TODO: Implement regular segmented reduce -> m*n parallel
  real_t mx = 0;
  int j = cta.get_local_id(0);
  for(int i=0;i<m;i++){ 	
    real_t row_norm = reduce_over_group(cta, abs(A[i*n+j]), plus<real_t>());
    mx = max(mx,row_norm);
  }
  return mx;  
}
template <typename real_t>
real_t reflection_vector(const group<1>& cta,
                         const real_t& a_i,const real_t& anorm)
{
    int i_tid = cta.get_local_id(0);
    real_t alpha = -sycl::copysign(anorm,a_i);
    real_t v_i = a_i + (i_tid==0)*alpha; // TODO: Check fortegn 
    real_t vnorm = sqrt(reduce_over_group(cta, v_i*v_i, plus<real_t>()));
    return v_i / vnorm;
}
template <typename real_t>
void QHQ(const group<1>& cta,
        local_accessor<real_t,1>& A, int n, //in/out
        local_accessor<real_t,1>& Q,        //in/out
        local_accessor<real_t,1>& vHA_wm,   //workmem
         bool compute_eigenvectors=true)
{
  real_t numerical_zero = max_norm(A,n,n)*10*std::numeric_limits<real_t>::epsilon();
  
  int j_tid = cta.get_local_id(0);
  for(int k=0;k<n-1;k++){ // Sequential loop
    int l = n-(k+1);		                                                //length of kth postdiagonal row a
    real_t a_i   = j_tid<l? A[k*n+(k+1)+j_tid] : 0;                         //a = A[k,(k+1):n], kth postdiagonal row.
    real_t anorm = sqrt(reduce_over_group(cta, a_i*a_i, plus<real_t>()));   //Norm of a

    if(anorm < numerical_zero) continue;                                    //Already eliminated, don't divide by 0

    real_t v_i =  reflection_vector(cta, a_i, anorm);                       //Vector definining elimination operations
    
    reflect_region(cta,A,n,k+1,k,l,l+1, v_i,2, ROWS, vHA_wm);
    reflect_region(cta,A,n,k+1,k,l,l+1, v_i,2, COLS, vHA_wm);

    if(compute_eigenvectors) reflect_region(cta,Q,n, k+1,0, l,n, v_i, sigma, ROWS, vHA_wm); 
  }
} */

template <typename T, int N>
void diagonalize(sycl::group<1>& cta, T* U, T* L, T* D, T* V, T* Q){
    auto tid = cta.get_local_linear_id();
    auto bid = cta.get_group_linear_id();
    auto bdim = cta.get_local_linear_range();

    for (int i = tid; i < N*N; i += bdim){
        Q[i] = 0;
    }
    sycl::group_barrier(cta);
    for (int i = tid; i < N; i += bdim){
        U[i]        = L[i]; //Matrix is symmetric
        Q[i*(N+1)] = 1;     //Initialize to identity matrix
    }
    sycl::group_barrier(cta);
    for(int k=N-1;k>=0;k--){
        T d = D[k];
        T shift = d;

        int i = 0;
        T GR = (k>0?std::abs(L[k-1]):0)+std::abs(L[k]);
        int not_done = 1;
        while (not_done > 0){
            i++;
            T_QTQ(cta, k+1, D, L, U, V, shift);
            apply_all_reflections(cta, V,k,N,Q);
            GR = (k>0?std::abs(L[k-1]):0)+(k+1<N?std::abs(L[k]):0);

            if(k>0){
                std::array<T,4> args = {D[k-1], L[k-1], L[k-1], D[k]};
                auto [l0, l1] = eigvalsh2x2(args);
                shift = std::abs(l0-d) < std::abs(l1-d)? l0:l1;
            } else {shift = D[k];}
        
            if(GR <= std::numeric_limits<T>::epsilon()*T(10.)) not_done--; // Do one (or optionally more) steps after reaching tolerance, to get all off-diagonals below.
                                            // GPU NB: Se GPU NB ovenfor.
            if(i>10){
                //printf("%dth run: Cannot converge eigenvalue %d to tolerance " G " using machine precision %g (d=%g, shift=%g, G=%g)\n" "D[k] = %g, L[k-1] = %g, L[k] = %g\n", nth_time,k,tolerance, std::numeric_limits<T>::epsilon(),d,shift,GR, D[k], (k>0)?L[k-1]:0, (k+1<n)?L[k]:0);
                auto max_error = sycl::max(std::numeric_limits<T>::epsilon()*T(10.),GR);
                break;
            }
        }
    }
    sycl::group_barrier(cta);
    
    //Eigenvalues now reside in D
}
//Finds the BlockVectors- smallest or largest eigenvalues of the matrix A
template <typename T, int N, int BlockVectors>
void compute_k_eigenpairs(sycl::group<1>& cta, const Span<T>& A, const Span<T>& eigvects, const Span<T>& lambdas, T* working_space, const local_accessor<std::byte,1>& sort_scratch, const bool largest){
    auto tid = cta.get_local_linear_id();
    auto bdim = cta.get_local_linear_range();

    T* LanczosVectors = working_space;
    //The memory used for X, alphas and betas can be reused after lanczos is done, LanczosVectors must be preserved for eigenvector computation
    T* X = LanczosVectors + N*N;
    T* alphas = X + N;
    T* betas = alphas + N + 1;

    T* D = alphas;
    T* L = betas;
    T* U = L + N + 1;
    T* V = U + N*2 + 1;
    T* Q = V + N*2 + 1;
    int* k_indices = reinterpret_cast<int*>(betas);
    lanczos<T, N>(cta, A, Span<T>(X,N), Span<T>(alphas,N+1), Span<T>(betas,N+1), Span<T>(LanczosVectors,N*N));
    diagonalize<T, N>(cta, U, L, D, V, Q);
    if(tid < N) k_indices[tid] = tid;
    sycl::group_barrier(cta);

    
    const auto bytes = sycl::ext::oneapi::experimental::default_sorters::joint_sorter<>::memory_required<T>(sycl::memory_scope::work_group, N);

    sycl::ext::oneapi::experimental::joint_sort(ext::oneapi::experimental::group_with_scratchpad(cta, sycl::span{(std::byte*)sort_scratch.get_pointer(), bytes}), 
                                                k_indices, 
                                                k_indices + N,
                                                [&](auto x, auto y){ return D[x] < D[y]; });



    

    if(tid < N){
        lambdas[tid] = largest ? D[k_indices[N - 1 - tid]] : D[k_indices[tid]];
    }
    //if (tid < BlockVectors) sycl::ext::oneapi::experimental::printf("Lambdas[%d] = %f\n", tid, lambdas[tid]);

    if(tid < N){
        for(int i = 0; i < BlockVectors; i++){
            auto index = largest ? k_indices[N -1 - i] : k_indices[i];
            eigvects[i*N + tid] = 0;
            for (int j = 0; j < N; j++){
                eigvects[i*N + tid] += LanczosVectors[j*N + tid] * Q[index*N + j];
            }
        }
    }

    //Normalize eigenvectors
    for(int i = 0; i < BlockVectors; i++){
        T val = tid < N ? eigvects[i*N + tid] : 0;
        T rnorm = rsqrt(reduce_over_group(cta, val * val, sycl::plus<T>{}));

        //Invert the sign of the eigenvector if the first element is negative
        T first = tid < N ? eigvects[i*N] : 1;
        group_barrier(cta);
        if(tid<N) eigvects[i*N + tid] = sign(first) * val * rnorm;
    }

    sycl::group_barrier(cta);
} 

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

template <typename T, int SN>
void compute_gram_matrix(sycl::group<1>& cta, T* St, T* S, T* StS, int m){
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
    for(int i = 0; i < SN; i++){
        for(int j = i; j < SN; j++){
            //StS[i*SN + j] = joint_reduce(cta, St + i*m, S + j*m)
        }
    }
    for (int ij = tid; ij < SN*SN; ij += bdim){
        int i = ij / SN;
        int j = ij % SN;
        StS[j*SN + i] = StS[i*SN + j];
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

//M != N = K in place matrix matrix multiplication
template <typename T, int SN>
void inPlaceMatMatMul(sycl::group<1>& cta, Span<T> A, Span<T> B, int m){
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
    T Alocal[SN];
    for(int j = 0; j < SN; j++){
        T sum = 0;
        for(int k = 0; k < SN; k++){
            sum += A[k*m + tid] * B[j*SN + k];
        }
        Alocal[j] = sum;
    }
    sycl::group_barrier(cta);
    for(int j = 0; j < SN; j++){
        A[j*m + tid] = Alocal[j];
    }   
}


template <typename T, int SN>
void inPlaceMatMatMul(sycl::group<1>& cta, Span<T> A, Span<T> B, int m, bool largest){
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
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


template <typename T, int SN>
void orthonormalize(sycl::group<1>& cta, Span<T> S, int m){
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
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
    auto bdim = cta.get_local_range(0);
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

template <typename T, int SN, int BlockVectors>
void update_vectors(sycl::group<1>& cta, Span<T> X, Span<T> R, Span<T> P, Span<T> AX, Span<T> AR, Span<T> AP, Span<T> eigenvects, int m, bool restart){
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
        T* eig_X = eigenvects.data() + i*(SN);
        T* eig_R = eigenvects.data() + i*(SN) + BlockVectors;
        T* eig_P = eigenvects.data() + i*(SN) + 2*BlockVectors;
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

template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest){
    static SyclVector<T> S(batch_size * BlockVectors*9 * m);
    static SyclVector<T> X0(batch_size * BlockVectors * m);
    static SyclVector<T> LastEigVects(batch_size * BlockVectors * BlockVectors*3);
    static SyclVector<T> LastGram(batch_size * 3 * BlockVectors * 3 * BlockVectors);
    static SyclVector<T> vals(batch_size * BlockVectors*3);
    static SyclVector<int> indices(batch_size * BlockVectors);

    auto S_acc = S.to_span();
    auto X0_acc = X0.to_span();
    auto LastEigVects_acc = LastEigVects.to_span();
    auto LastGram_acc = LastGram.to_span();
    auto keys_acc = indices.to_span();
    auto vals_acc = vals.to_span();
    auto indices_acc = indices.to_span();

    const T tol = std::numeric_limits<T>::epsilon() * sycl::sqrt(float(m));
    //
    ctx -> submit([&](sycl::handler& h){
        using TupleType = typename std::iterator_traits<oneapi::dpl::zip_iterator<T*, int*>>::value_type;
        auto StAS_smem      = sycl::local_accessor<T, 1>(BlockVectors*3 * BlockVectors*3, h);
        auto eigvects_smem  = sycl::local_accessor<T, 1>(BlockVectors * BlockVectors*3,h);
        auto lambdas_smem   = sycl::local_accessor<T, 1>(BlockVectors, h);
        auto working_space  = sycl::local_accessor<T, 1>(
                        BlockVectors*3 * BlockVectors*3 * 2 + //Lanczos Vectors and transformation matrix Q
                        BlockVectors*3*8                      /*Everything else */, h);
        auto Scache         = sycl::local_accessor<T, 1>(m, h);
        
        const auto bytes = sycl::ext::oneapi::experimental::default_sorters::joint_sorter<>::memory_required<TupleType>(sycl::memory_scope::work_group, BlockVectors*3);
        auto sort_scratchpad = local_accessor<std::byte, 1>(bytes, h);

        h.parallel_for<LOBPCGg<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            sycl::group<1> cta = item.get_group();
            auto bid = item.get_group_linear_id();
            constexpr auto SN = BlockVectors*3;
            std::bitset<BlockVectors> converged;
            oneapi::dpl::uniform_real_distribution<T> distr(0.0, 1.0);            
            oneapi::dpl::minstd_rand engine(42, tid);
            auto A_acc = A.subspan(bid * m * NZ, m * NZ);
            auto cols_acc = cols.subspan(bid * m * NZ, m * NZ);
            auto StAS = Span<T>(StAS_smem.get_pointer(), BlockVectors*3 * BlockVectors*3);
            auto eigvects = Span<T>(eigvects_smem.get_pointer(), BlockVectors * BlockVectors*3);
            auto lambdas = Span<T>(lambdas_smem.get_pointer(), BlockVectors*3);


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
            DECLARE_SUBSPACE

            for (int i = tid + m*BlockVectors; i < 9*m*BlockVectors; i+=cta.get_local_range(0)){
                subspace[i] = T(0);
            }
            for (int i = tid; i < m*BlockVectors; i+=cta.get_local_range(0)){
                blockX[i] = distr(engine);
                //X0_acc[i] = blockX[i];
                X0_acc[bid*m*BlockVectors + i] = blockX[i];
            }
             
            group_barrier(cta);
            //Normalize S vectors
            orthonormalize<T, BlockVectors>(cta, blockX, m); //Modified Gram-Schmidt

            //A * X
            //matBlockVector<T, K, 1, N>(cta, A_tid, cols_tid, blockX, blockAX, NZ, m); 
            //X^T * A * X
            
            STAS<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockX, m, NZ, StAS, Scache.get_pointer());
            //Fill in with simple 3x3 matrix
            compute_k_eigenpairs<T, BlockVectors, BlockVectors>(cta, StAS, eigvects, lambdas, working_space.get_pointer(), sort_scratchpad, largest);
            matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockX, blockAX, NZ, m);
            inPlaceMatMatMul<T, BlockVectors>(cta, blockX, eigvects, m);
            inPlaceMatMatMul<T, BlockVectors>(cta, blockAX, eigvects, m);

            //Compute the residual R = A*X - X*Lambda
            int iter = 0;
            bool restart = true;
            while(!converged.all() && iter < maxiters){
                //R = A*X - X*Lambda
                for(int i = 0; i < BlockVectors; i++) blockR[i*m + tid] = blockAX[i*m + tid] - lambdas[i] * blockX[i*m + tid];
                //Convergence Check
                for(int i = 0; i < BlockVectors; i++) {if(converged[i]) continue; converged[i] = sycl::sqrt(std::abs(reduce_over_group(cta, blockR[i*m + tid]*blockR[i*m + tid], sycl::plus<T>{}))) < tol;}
                sycl::group_barrier(cta);


                //R = R - X * (X^T * R)
                applyConstraints<T, BlockVectors>(cta, blockR, blockX, m);
                sycl::group_barrier(cta);
                orthonormalize<T, BlockVectors>(cta, blockR, m);
                sycl::group_barrier(cta);
                matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockR, blockAR, NZ, m);
                /* 
                 */

                if (!restart) {
                     orthonormalize<T, BlockVectors>(cta, blockP, m);
                    for(int i = 0; i < BlockVectors; i++){ 
                        blockXtemp[i*m + tid] = blockX[i*m + tid];
                        blockRtemp[i*m + tid] = blockR[i*m + tid];
                        blockPtemp[i*m + tid] = blockP[i*m + tid];
                    }
                    matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockP, blockAP, NZ, m);
                    orthonormalize<T, BlockVectors*3>(cta, blockX, m);
                    STAS<T, K, 1, BlockVectors*3>(cta, A_tid, cols_tid, subspace, m, NZ, StAS, Scache.get_pointer());
                    compute_k_eigenpairs<T, BlockVectors*3, BlockVectors>(cta, StAS, eigvects, lambdas, working_space.get_pointer(), sort_scratchpad, largest);
                    
                    for (int i = 0; i < BlockVectors; i++){ 
                        blockX[i*m + tid] = blockXtemp[i*m + tid];
                        blockR[i*m + tid] = blockRtemp[i*m + tid];
                        blockP[i*m + tid] = blockPtemp[i*m + tid];
                    } 
                    group_barrier(cta);
                    update_vectors<T, BlockVectors*3, BlockVectors>(cta, blockX, blockR, blockP, blockAX, blockAR, blockAP, eigvects, m, false);
                    group_barrier(cta);

                } else {
                    for (int i = 0; i < BlockVectors; i++){
                        blockXtemp[i*m + tid] = blockX[i*m + tid];
                        blockRtemp[i*m + tid] = blockR[i*m + tid];
                    }
                    orthonormalize<T, BlockVectors*2>(cta, blockX, m);
                    STAS<T, K, 1, BlockVectors*2>(cta, A_tid, cols_tid, subspace, m, NZ, StAS, Scache.get_pointer());
                    compute_k_eigenpairs<T, BlockVectors*2, BlockVectors>(cta, StAS, eigvects, lambdas, working_space.get_pointer(), sort_scratchpad, largest);
                    for (int i = 0; i < BlockVectors; i++){
                        blockX[i*m + tid] = blockXtemp[i*m + tid];
                        blockR[i*m + tid] = blockRtemp[i*m + tid];
                    }
                    group_barrier(cta);
                    update_vectors<T, BlockVectors*2, BlockVectors>(cta, blockX, blockR, blockP, blockAX, blockAR, blockAP, eigvects, m, true);
                    restart = false;
                }
                sycl::group_barrier(cta);
                iter++;
            }
            for (int i = tid; i < BlockVectors*BlockVectors*3; i+=cta.get_local_range(0)){
                LastEigVects_acc[i] = eigvects[i];
            }
            for (int i = tid; i < 3*BlockVectors*3*BlockVectors; i+=cta.get_local_range(0)){
                LastGram_acc[i] = StAS[i];
            }
            if(tid < BlockVectors*3)   vals_acc[tid] = lambdas[tid];
        });
    });
    std::string dtype = std::is_same<T, float>::value ? ".float32" : ".float64";
    ctx.wait();
    /* int SD = std::min(int(maxiters+1), 3);
    std::ofstream file("Subspace_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        file.write(reinterpret_cast<const char*>(S_acc.data()), SD * BlockVectors * m * sizeof(T));
    }

    std::ofstream file2("AXARAP_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        file2.write(reinterpret_cast<const char*>(S_acc.data() + 6 * BlockVectors * m), 3 * BlockVectors * m * sizeof(T));
    }

    std::ofstream file3("Initial_X_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + dtype, ios::binary);
    {
        file3.write(reinterpret_cast<const char*>(X0_acc.data()), BlockVectors * m * sizeof(T));
    }

    std::ofstream file4("LastEigVects_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        file4.write(reinterpret_cast<const char*>(LastEigVects_acc.data()), BlockVectors * BlockVectors * SD * sizeof(T));
    }

    std::ofstream file5("LastGram_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        file5.write(reinterpret_cast<const char*>(LastGram_acc.data()), SD * BlockVectors * SD * BlockVectors * sizeof(T));
    }

    std::ofstream file6("Lambdas_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        file6.write(reinterpret_cast<const char*>(vals_acc.data()), BlockVectors * sizeof(T));
    }
    */
    std::cout << vals_acc << std::endl; 
}


template <typename T, typename K>
struct LOBPCG_Helper{
    cusolverDnHandle_t handle;
    cusolverStatus_t solver_status;
    cusparseStatus_t sparse_status;
    cublasStatus_t cublas_status;

    size_t SpMM_buffer_size;
    SyclVector<T> SpMM_buffer;
    size_t h_SYEV_buffer_size;
    SyclVector<T> h_SYEV_buffer;
    size_t d_SYEV_buffer_size;
    SyclVector<T> d_SYEV_buffer;
    
    props_t props;
    params_t params;

    cusparseHandle_t SpMM_AX_handle;
    cusparseHandle_t SpMM_AXR_handle;
    cusparseHandle_t SpMM_AS_handle;
    cusparseSpMatDescr_t descrA;    //Left matrix
    cusparseDnMatDescr_t descrX;    //Right matrix
    cusparseDnMatDescr_t descrXR;   //Concatination matrix [X, R]
    cusparseDnMatDescr_t descrS;    //The full subspace matrix [X, R, P]
    cusparseDnMatDescr_t descrAX;  //The result of the matrix matrix multiplication
    cusparseDnMatDescr_t descrAXR; //The result of the matrix matrix multiplication
    cusparseDnMatDescr_t descrAS;  //The result of the matrix matrix multiplication

    cublasHandle_t GEMM_XAX_handle;
    cublasHandle_t GEMM_XRAXR_handle;
    cublasHandle_t GEMM_SAS_handle;

    LOBPCG_Helper(){
        create_handle(&handle);
        create_params(&params);
        sparse_status = cusparseCreate(&SpMM_AX_handle);
        sparse_status = cusparseCreate(&SpMM_AXR_handle);
        sparse_status = cusparseCreate(&SpMM_AS_handle);

        cublas_status = cublasCreate(&GEMM_XAX_handle);
        cublas_status = cublasCreate(&GEMM_XRAXR_handle);
        cublas_status = cublasCreate(&GEMM_SAS_handle);
    }

    LOBPCG_Helper(cudaStream_t& stream, int k, int n, int nz, int batch_size, T* A, K* cols, K* csr_offsets,
                    T* S, T* SAS, T* AS, T* lambda){
        
        T alpha = 1.0;
        T beta = 0.0;
        set_stream(handle, stream);

        sparse_status = cusparseCreateCsr(&descrA, n, n, nz, csr_offsets, cols, A, CUSPARSE_INDEX_16U, CUSPARSE_INDEX_16U, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
        sparse_status = cusparseCsrSetStridedBatch(descrA, batch_size, n*nz, 0);
        sparse_status = cusparseCreateDnMat(&descrX, n, k, n, S, CUDA_R_32F, CUSPARSE_ORDER_COL);
        sparse_status = cusparseDnMatSetStridedBatch(descrX, batch_size, n*k);
        sparse_status = cusparseCreateDnMat(&descrXR, n, 2*k, n, S, CUDA_R_32F, CUSPARSE_ORDER_COL);
        sparse_status = cusparseDnMatSetStridedBatch(descrXR, batch_size, n*k);
        sparse_status = cusparseCreateDnMat(&descrS, n, 3*k, n, S, CUDA_R_32F, CUSPARSE_ORDER_COL);
        sparse_status = cusparseDnMatSetStridedBatch(descrS, batch_size, n*k);
        sparse_status = cusparseCreateDnMat(&descrAX, k, k, k, SAS, CUDA_R_32F, CUSPARSE_ORDER_COL);
        sparse_status = cusparseDnMatSetStridedBatch(descrAX, batch_size, k*k);
        sparse_status = cusparseCreateDnMat(&descrAXR, 2*k, 2*k, 2*k, SAS, CUDA_R_32F, CUSPARSE_ORDER_COL);
        sparse_status = cusparseDnMatSetStridedBatch(descrAXR, batch_size, k*k);
        sparse_status = cusparseCreateDnMat(&descrAS, 3*k, 3*k, 3*k, SAS, CUDA_R_32F, CUSPARSE_ORDER_COL);
        sparse_status = cusparseDnMatSetStridedBatch(descrAS, batch_size, k*k);

        solver_status = cusolverDnXsyevBatched_bufferSize(  handle,
                                                            params,
                                                            CUSOLVER_EIG_MODE_VECTOR, 
                                                            CUBLAS_FILL_MODE_LOWER,
                                                            k*3,
                                                            CUDA_R_32F,
                                                            SAS,
                                                            k*3,
                                                            CUDA_R_32F,
                                                            lambda,
                                                            CUDA_R_32F,
                                                            &d_SYEV_buffer_size,
                                                            &h_SYEV_buffer_size,
                                                            batch_size);
        
        sparse_status = cusparseSpMM_bufferSize(handle, 
                                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                &alpha,
                                                descrA,
                                                descrS,
                                                &beta,
                                                descrAS,
                                                CUDA_R_32F,
                                                CUSPARSE_SPMM_ALG_DEFAULT,
                                                &SpMM_buffer_size);
                                                
        h_SYEV_buffer.resize(h_SYEV_buffer_size);
        d_SYEV_buffer.resize(d_SYEV_buffer_size);
        SpMM_buffer.resize(SpMM_buffer_size);
    }

    ~LOBPCG_Helper(){
        cusolverDnDestroy(handle);
        cusparseDestroy(SpMM_AX_handle);
        cusparseDestroy(SpMM_AXR_handle);
        cusparseDestroy(SpMM_AS_handle);
        cusparseDestroyDnMat(descrAS);
        cusparseDestroyDnMat(descrAX);
        cusparseDestroyDnMat(descrXR);
        cusparseDestroyDnMat(descrX);
        cusparseDestroySpMat(descrA);

        cublasDestroy(GEMM_XAX_handle);
        cublasDestroy(GEMM_XRAXR_handle);
        cublasDestroy(GEMM_SAS_handle);
    }
};


template <typename T, typename K, int BlockVectors, int NZ>
void LOBPCG_V1(SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t maxiters, bool largest){
    SyclVector<T> S(batch_size * BlockVectors*9 * m);
    SyclVector<T> X0(batch_size * BlockVectors * m);
    SyclVector<T> LastEigVects(batch_size * BlockVectors * BlockVectors*3);
    SyclVector<T> LastGram(batch_size * 3 * BlockVectors * 3 * BlockVectors);
    SyclVector<T> vals(batch_size * BlockVectors);
    SyclVector<T> StAS(batch_size * BlockVectors*3 * BlockVectors*3);
    auto S_acc = S.to_span(); 
    auto X0_acc = X0.to_span();
    auto StAS_acc = StAS.to_span();


    //sycl::buffer<T, 1> lambdas(batch_size * BlockVectors*3);
    //sycl::buffer<int, 1> indices(batch_size * BlockVectors);
    
    size_t syevd_scratchpad_size = 0;
    size_t syevd_scratchpad_host_size = 0;
    SyclVector<int> syevd_info(batch_size);
    SyclVector<T> lambdas(batch_size * BlockVectors*3);
    SyclVector<int> csr_row_offsets(m);
    
    //csr_row_offsets = [0, NZ, 2*NZ, 3*NZ, ...]
    //primitives::iota(ctx, csr_row_offsets, 0);
    //primitives::transform(ctx, csr_row_offsets, csr_row_offsets, [&](int i){return i * NZ;}); 

    auto lambda_acc = lambdas.to_span();

    handle_t handle{};
    stream_t stream{};
    props_t props{};
    params_t params{};
    create_handle(&handle);
    create_stream(&stream);
    set_stream(handle, stream);
    create_params(&params);

    #ifdef USE_CUSOLVER
        status_t status = cusolverDnXsyevBatched_bufferSize(handle, 
                                params, 
                                EIG_SOLVE_MODE, 
                                FILL_MODE_LOWER,
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
    #elif defined(USE_ROCSOLVER)
        status_t status = rocsolver_ssyevd_bufferSize(handle, params, EIG_SOLVE_MODE, FILL_MODE_LOWER, m, A.data(), m, lambdas.data(), &syevd_scratchpad_size, &syevd_scratchpad_host_size, batch_size);
    #endif

    SyclVector<T> syevd_scratchpad (syevd_scratchpad_size);
    SyclVector<T> syevd_scratchpad_host (syevd_scratchpad_host_size);

    #ifdef USE_CUSOLVER
        auto compute_eigenpairs = [&](T* A, T* lambdas, int m, int batch_size){
            status_t status = cusolverDnXsyevBatched(handle, 
                                params, 
                                EIG_SOLVE_MODE, 
                                FILL_MODE_LOWER,
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
        };
    #elif defined(USE_ROCSOLVER)
        auto compute_eigenpairs = [&](T* A, T* lambdas, T* eigenvectors, int m, int batch_size){
            status_t status = rocsolver_ssyevd(handle, params, EIG_SOLVE_MODE, FILL_MODE_LOWER, m, A, m, lambdas, eigenvectors, m, syevd_scratchpad, syevd_scratchpad_size, syevd_scratchpad_host_size, syevd_info.data(), batch_size);
        };
    #endif
    const T tol = sycl::sqrt(std::numeric_limits<T>::epsilon()) * T(m);

    //Setup
    ctx -> submit([&](sycl::handler& h ){
        auto Scache = sycl::local_accessor<T, 1>(m, h);

        h.parallel_for<LOBPCG_V1_1<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            sycl::group<1> cta = item.get_group();
            auto bid = item.get_group_linear_id();
            constexpr auto SN = BlockVectors*3;
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

            for (int i = tid + m*BlockVectors; i < 9*m*BlockVectors; i+=cta.get_local_range(0)){
                subspace[i] = T(0);
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
            

            auto StASspan = StAS_acc.subspan(bid * BlockVectors * BlockVectors, BlockVectors * BlockVectors);
            //A * X
            matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockX, blockAX, NZ, m); 
            //X^T * A * X
            STAS<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockX, m, NZ, StASspan, Scache.get_pointer());
        });
    });
    ctx.wait();
    //std::cout << "StAS: \n" << StAS << std::endl;
    //Solve the eigenvalue problem
    compute_eigenpairs(StAS.data(), lambdas.data(), BlockVectors, batch_size);
    sync_stream(stream);
    //std::cout << "Lambdas: \n" << lambdas << std::endl;
    //std::cout << "Eigenvectors: \n" << StAS << std::endl;



    //
    ctx -> submit([&](sycl::handler& h){
        auto Scache = sycl::local_accessor<T, 1>(m, h);

        h.parallel_for<LOBPCG_V1_2<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            sycl::group<1> cta = item.get_group();
            auto bid = item.get_group_linear_id();
            constexpr auto SN = BlockVectors*3;

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
    auto print_blockvectors = [&](Span<T> block, int num_vects = BlockVectors){
        std::cout << "[";
        for (int j = 0; j < m; j++){
            std::cout << "[";
            for (int i = 0; i < num_vects; i++){
                std::cout << block[i*m + j] << " ";
            }
            std::cout << "]" << std::endl;
        }
        std::cout << "]" << std::endl;
    };

    SyclVector<std::bitset<BlockVectors>> converged(batch_size);
    SyclVector<T> residuals(batch_size * BlockVectors);
    auto converged_acc = converged.to_span();
    auto residuals_acc = residuals.to_span();

    auto iteration = [&](SyclQueue &ctx, Span<T> A, Span<K> cols, int batch_size, int m, size_t iter, bool restart){
        
        //Setup
        ctx -> submit([&](sycl::handler& h){
            auto Scache = sycl::local_accessor<T, 1>(m, h);

            h.parallel_for<LOBPCG_V1_3<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
                auto tid = item.get_local_linear_id();
                auto cta = item.get_group();
                auto bid = item.get_group_linear_id();
                constexpr auto SN = BlockVectors*3;
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
                    if(converged[i]) continue; 
                    residuals[i] = sycl::sqrt(reduce_over_group(cta, blockR[i*m + tid]*blockR[i*m + tid], sycl::plus<T>{}));
                    converged[i] = residuals[i] < tol;
                }
                sycl::group_barrier(cta);
                //R = R - X * (X^T * R)
                applyConstraints<T, BlockVectors>(cta, blockR, blockX, m);
                sycl::group_barrier(cta);
                orthonormalize<T, BlockVectors>(cta, blockR, m);
                matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockR, blockAR, NZ, m);
                auto StAS_offset = restart ? bid * BlockVectors*2 * BlockVectors*2 : bid * BlockVectors*3 * BlockVectors*3;
                auto StAS_size = restart ? BlockVectors*2 * BlockVectors*2 : BlockVectors*3 * BlockVectors*3;
                auto StAS = StAS_acc.subspan(StAS_offset, StAS_size);
                if (!restart) {
                    orthonormalize<T, BlockVectors>(cta, blockP, m);
                    for(int i = 0; i < BlockVectors; i++){ 
                        blockXtemp[i*m + tid] = blockX[i*m + tid];
                        blockRtemp[i*m + tid] = blockR[i*m + tid];
                        blockPtemp[i*m + tid] = blockP[i*m + tid];
                    }
                    matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockP, blockAP, NZ, m);
                    orthonormalize<T, BlockVectors*3>(cta, blockX, m);
                    STAS<T, K, 1, BlockVectors*3>(cta, A_tid, cols_tid, subspace, m, NZ, StAS, Scache.get_pointer());
                } else {
                    for (int i = 0; i < BlockVectors; i++){
                        blockXtemp[i*m + tid] = blockX[i*m + tid];
                        blockRtemp[i*m + tid] = blockR[i*m + tid];
                    }
                    orthonormalize<T, BlockVectors*2>(cta, blockX, m);
                    STAS<T, K, 1, BlockVectors*2>(cta, A_tid, cols_tid, subspace, m, NZ, StAS, Scache.get_pointer());
                }
            });
        });
        ctx -> wait();


        /* std::cout << "Subspace Iteration " << iter << ": \n";
        print_blockvectors(S_acc); */

        //std::cout << "StAS: \n" << StAS << std::endl;
        //Solve the eigenvalue problem
        compute_eigenpairs(StAS.data(), lambdas.data(), restart ? BlockVectors*2 : BlockVectors*3, batch_size);
        sync_stream(stream);

        //std::cout << "EigenVectors: \n" << StAS_acc << std::endl;
        //std::cout << "Lambdas: \n" << lambda_acc << std::endl;


        ctx -> submit([&](sycl::handler& h){
            auto Scache = sycl::local_accessor<T, 1>(m, h);
            h.parallel_for<LOBPCG_V1_4<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(batch_size*m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
                auto tid = item.get_local_linear_id();
                sycl::group<1> cta = item.get_group();
                auto bid = item.get_group_linear_id();
                constexpr auto SN = BlockVectors*3;
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
                auto eigvects_offset = restart ? bid * BlockVectors * 2 * BlockVectors*2 : bid * BlockVectors * 3 * BlockVectors*3;
                auto eigvects_size = restart ? BlockVectors * 2 * BlockVectors*2 : BlockVectors * 3 * BlockVectors*3;
                auto eigvects = StAS_acc.subspan(eigvects_offset, eigvects_size);

                if (!restart) {
                    for (int i = 0; i < BlockVectors; i++){ 
                        blockX[i*m + tid] = blockXtemp[i*m + tid];
                        blockR[i*m + tid] = blockRtemp[i*m + tid];
                        blockP[i*m + tid] = blockPtemp[i*m + tid];
                    } 
                    group_barrier(cta);
                    update_vectors<T, BlockVectors*3, BlockVectors>(cta, blockX, blockR, blockP, blockAX, blockAR, blockAP, eigvects, m, false, largest);
                } else {
                    for (int i = 0; i < BlockVectors; i++){
                        blockX[i*m + tid] = blockXtemp[i*m + tid];
                        blockR[i*m + tid] = blockRtemp[i*m + tid];
                    }
                    group_barrier(cta);
                    update_vectors<T, BlockVectors*2, BlockVectors>(cta, blockX, blockR, blockP, blockAX, blockAR, blockAP, eigvects, m, true, largest);
                }
            });
        });
    };

    for(size_t i = 0; i < maxiters; i++){
        iteration(ctx, A, cols, batch_size, m, i, i < 1 ? true : false);
    }
    ctx.wait();
     

    //std::string dtype = std::is_same<T, float>::value ? ".float32" : ".float64";
   /*  int SD = std::min(int(maxiters+1), 3);
    
    
    Span<T> lambdas_vec(lambdas.begin(), lambdas.end());
    if (largest) std::sort(lambdas_vec.begin(), lambdas_vec.end(), std::greater<T>());
    else std::sort(lambdas_vec.begin(), lambdas_vec.end()); */
    
    

    std::cout << "Eigenvalues: " << lambdas << std::endl;
}


template void LOBPCG<float, uint16_t, 3, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest);
template void LOBPCG<float, uint16_t, 21, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest);
template void LOBPCG_V1<float, uint16_t, 3, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest);
template void LOBPCG_V1<float, uint16_t, 21, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters, bool largest);
//template void LOBPCG<float, uint16_t, 9, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters);
//template void LOBPCG<double, uint16_t, 21, 30>(SyclQueue &ctx, Span<double> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters);
//template void LOBPCG<double, uint16_t, 3, 30>(SyclQueue &ctx, Span<double> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters);
//template void LOBPCG<double, uint16_t, 9, 30>(SyclQueue &ctx, Span<double> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters);

//template void LOBPCG_V1<float, uint16_t, 21, 30>(SyclQueue &ctx, Span<float> A, Span<uint16_t> cols, int batch_size, int m, size_t maxiters);