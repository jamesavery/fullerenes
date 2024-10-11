#include <array>
#include <algorithm>
#include <fstream>
#include "sycl/sycl.hpp"
#include <oneapi/dpl/random>
#include <oneapi/dpl/iterator>
#include <fullerenes/graph.hh>
#include <fullerenes/sycl-headers/all-kernels.hh>
#include <../src/sycl/queue-impl.cc>
#include <iostream>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/argparser.hh>
#include <fullerenes/sycl-headers/fill.hh>
#define CEIL_DIV(x, y) (((x) + (y) - 1) / (y))
template <typename T, typename K, int BlockVectors, int NZ> class LOBPCGg {};

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
        local_max = std::max(local_max, sycl::abs(D[i]) + 2*sycl::abs(L[i]));
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
            if (sycl::abs(L[k]) > numerical_zero){
            a[0] = D[k]; a[1] = L[k];       // a = T[k:k+2,k] is the vector of nonzeros in kth subdiagonal column.
            
            real_t anorm = sqrt(a[0]*a[0] + a[1]*a[1]); 

            // // Udrullet
            // //    reflection_vector(a,anorm,v);
            v[0] = D[k]; v[1] = L[k];
            real_t alpha = -sycl::copysign(anorm,a[0]); // Koster ingenting
            v[0] -= alpha;

            real_t vnorm = sqrt(v[0]*v[0]+v[1]*v[1]);
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
    T D = sqrt(4*b*c+(a-d)*(a-d));
    return {(a+d-D)/2, (a+d+D)/2};
}
//Assumes A is symmetric
template <typename T, int N> 
void lanczos(sycl::group<1>& cta, const local_accessor<T,1>& A, T* X, T* alphas, T* betas, T* V){
    auto tid = cta.get_local_linear_id();
    V = V + tid;
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
        result /= sqrt(reduce_over_group(cta, result * result, sycl::plus<T>{}));
        return result;
    };

    auto v0 = (tid < N) ? distr(engine) : T(0);
    auto norm = sqrt(reduce_over_group(cta, v0 * v0, sycl::plus<T>{}));
    if (tid < N)  V[0*N] = v0 / norm;
    v0 = MGS(0);
    if (tid < N) V[0*N] = v0;
    for (int i = 0; i < N; i++){
        if (i % 2 == 0 && i > 1){
            T v_i = MGS(i-1);
            T v_ip = MGS(i);
            if(tid < N) V[(i-1)*N] = v_i;
            if(tid < N) V[i*N] = v_ip;
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
        T beta = sqrt(reduce_over_group(cta, v * v, sycl::plus<T>{}));
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
        T GR = (k>0?sycl::abs(L[k-1]):0)+sycl::abs(L[k]);
        int not_done = 1;
        while (not_done > 0){
            i++;
            T_QTQ(cta, k+1, D, L, U, V, shift);
            apply_all_reflections(cta, V,k,N,Q);
            GR = (k>0?sycl::abs(L[k-1]):0)+(k+1<N?sycl::abs(L[k]):0);

            if(k>0){
                std::array<T,4> args = {D[k-1], L[k-1], L[k-1], D[k]};
                auto [l0, l1] = eigvalsh2x2(args);
                shift = sycl::abs(l0-d) < sycl::abs(l1-d)? l0:l1;
            } else {shift = D[k];}
        
            if(GR <= std::numeric_limits<T>::epsilon()*T(10.)) not_done--; // Do one (or optionally more) steps after reaching tolerance, to get all off-diagonals below.
                                            // GPU NB: Se GPU NB ovenfor.
            if(i>10){
                //printf("%dth run: Cannot converge eigenvalue %d to tolerance " G " using machine precision %g (d=%g, shift=%g, G=%g)\n" "D[k] = %g, L[k-1] = %g, L[k] = %g\n", nth_time,k,tolerance, std::numeric_limits<T>::epsilon(),d,shift,GR, D[k], (k>0)?L[k-1]:0, (k+1<n)?L[k]:0);
                auto max_error = std::max(std::numeric_limits<T>::epsilon()*T(10.),GR);
                break;
            }
        }
    }
    sycl::group_barrier(cta);
    
    //Eigenvalues now reside in D
}
//Finds the BlockVectors' smallest or largest eigenvalues of the matrix A
template <typename T, int N, int BlockVectors>
void compute_k_eigenpairs(sycl::group<1>& cta, const local_accessor<T,1>& A, const local_accessor<T,1>& eigvects, const local_accessor<T,1>& lambdas, T* working_space, const local_accessor<std::byte,1>& sort_scratch, const bool largest){
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
    lanczos<T, N>(cta, A, X, alphas, betas, LanczosVectors);
    diagonalize<T, N>(cta, U, L, D, V, Q);
    if(tid < N) k_indices[tid] = tid;
    sycl::group_barrier(cta);

    constexpr auto bytes = sycl::ext::oneapi::experimental::default_sorter<>::memory_required<T>(sycl::memory_scope::work_group, N);

    sycl::ext::oneapi::experimental::joint_sort(ext::oneapi::experimental::group_with_scratchpad(cta, sycl::span<std::byte,bytes>(sort_scratch.get_pointer(), N)), 
                                                k_indices, 
                                                k_indices + N,
                                                [&](auto x, auto y){ return D[x] < D[y]; });



    

    if(tid < BlockVectors){
        lambdas[tid] = largest ? D[k_indices[N - 1 - tid]] : D[k_indices[tid]];
    }
    if(tid < BlockVectors) sycl::ext::oneapi::experimental::printf("Lambdas[%d] = %f\n", tid, lambdas[tid]);

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
void matBlockVector(sycl::group<1>& cta ,const private_ptr<T> A, const private_ptr<K> cols, const global_ptr<T> blockX, global_ptr<T> AX, int n, int m){    
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
void STAS(sycl::group<1>& cta, const private_ptr<T> A, const private_ptr<K> cols, global_ptr<T> S, int m, int n, local_ptr<T> StAS, local_ptr<T> Scache){
//S^T * A * S
//To produce the i,j-th element of StAS, we need to compute the dot product of the i-th row of S^T and the j-th column of A*S
//The j-th column of A*S is the matrix vector product of A and the j-th column of S
//For this reason it makes sense to store the j-th column of S in shared memory and tcd .he i-th row of S^T in shared memory (the i-th column of S)

//global_ptr<T> S_global(S);
//local_ptr<T> S_j(Scache, m);
    auto tid = cta.get_local_id(0);
    for(int j = 0; j < SN; j++){
        //Load the j-th column of S into shared memory
        auto event = cta.async_work_group_copy(Scache, S + j*m, m);
        event.wait();
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
void inPlaceMatMatMul(sycl::group<1>& cta, global_ptr<T> A, local_ptr<T> B, int m){
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
void orthonormalize(sycl::group<1>& cta, T* S, int m){
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
    for(int i = 0; i < SN; i++){
        T* S_i = S + i*m;
        T norm = sqrt(reduce_over_group(cta, S_i[tid]*S_i[tid], sycl::plus<T>{}));
        S_i[tid] /= norm;
        sycl::group_barrier(cta);
        for(int j = i+1; j < SN; j++){
            T* S_j = S + j*m;
            T projection = reduce_over_group(cta, S_i[tid]*S_j[tid], sycl::plus<T>{}) * S_i[tid];
            sycl::group_barrier(cta);
            S_j[tid] -= projection;
        }
        sycl::group_barrier(cta);
    }
}

template<typename T, int SN>
void applyConstraints(sycl::group<1>& cta, global_ptr<T> S, const global_ptr<T> C, int m){
    //S = S - C @ (C^T @ S)
    auto tid = cta.get_local_id(0);
    auto bdim = cta.get_local_range(0);
    for(int i = 0; i < SN; i++){
        T* S_i = S + i*m;
        T projection = 0;
        for(int j = 0; j < SN; j++){
            T* C_j = C + i*m;
            projection += sycl::reduce_over_group(cta, S_i[tid]*C_j[tid], sycl::plus<T>{}) * C_j[tid];
        }
        sycl::group_barrier(cta);
        S_i[tid] -= projection;
    }
}

template <typename T, int SN, int BlockVectors>
void update_vectors(sycl::group<1>& cta, global_ptr<T> X, global_ptr<T> R, global_ptr<T> P, global_ptr<T> AX, global_ptr<T> AR, global_ptr<T> AP, local_ptr<T> eigenvects, int m, bool restart){
    auto tid = cta.get_local_id(0);
    group_barrier(cta);
    T Ps[BlockVectors], APs[BlockVectors], Xs[BlockVectors], AXs[BlockVectors];
    for(int i = 0; i < BlockVectors; i++){
        T* X_i = X + i*m; T Xl = X_i[tid];
        T* R_i = R + i*m; T Rl = R_i[tid];
        T* P_i = P + i*m; T Pl = restart ? 0 : P_i[tid];
        T* AX_i = AX + i*m; T AXl = AX_i[tid];
        T* AR_i = AR + i*m; T ARl = AR_i[tid];
        T* AP_i = AP + i*m; T APl = restart ? 0 : AP_i[tid];
        T* eig_X = eigenvects + i*(SN);
        T* eig_R = eigenvects + i*(SN) + BlockVectors;
        T* eig_P = eigenvects + i*(SN) + 2*BlockVectors;
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
void LOBPCG(SyclQueue &ctx, Span<T> A, Span<K> cols, int m, size_t maxiters){
    sycl::buffer<T, 1> S(BlockVectors*9 * m);
    sycl::buffer<T, 1> X0(BlockVectors * m);
    sycl::buffer<T, 1> LastEigVects(BlockVectors * BlockVectors*3);
    sycl::buffer<T, 1> LastGram(3 * BlockVectors * 3 * BlockVectors);
    sycl::buffer<T, 1> vals(BlockVectors);
    sycl::buffer<int, 1> indices(BlockVectors);
    constexpr T tol = 1e-6;
    //
    ctx -> submit([&](sycl::handler& h){
        using TupleType = typename std::iterator_traits<oneapi::dpl::zip_iterator<T*, int*>>::value_type;
        constexpr auto bytes = sycl::ext::oneapi::experimental::default_sorter<>::memory_required<TupleType>(sycl::memory_scope::work_group, BlockVectors*3);

        auto A_acc = A;
        auto cols_acc = cols;
        auto S_acc = sycl::accessor<T, 1, sycl::access::mode::write>(S, h);
        auto X0_acc = sycl::accessor<T, 1, sycl::access::mode::write>(X0, h);
        auto LastEigVects_acc = sycl::accessor<T, 1, sycl::access::mode::write>(LastEigVects, h);
        auto LastGram_acc = sycl::accessor<T, 1, sycl::access::mode::write>(LastGram, h);

        auto StAS = sycl::local_accessor<T, 1>(BlockVectors*3 * BlockVectors*3, h);
        auto eigvects = sycl::local_accessor<T, 1>(BlockVectors * BlockVectors*3,h);
        auto lambdas = sycl::local_accessor<T, 1>(BlockVectors, h);
        auto working_space = sycl::local_accessor<T, 1>(BlockVectors*3 * BlockVectors*3 * 3, h);
        auto Scache = sycl::local_accessor<T, 1>(m, h);
        auto sort_scratchpad = local_accessor<std::byte, 1>(bytes, h);

        auto keys_acc = sycl::accessor<int, 1, sycl::access::mode::read_write>(indices, h);
        auto vals_acc = sycl::accessor<T, 1, sycl::access::mode::read_write>(vals, h);

        h.parallel_for<LOBPCGg<T,K,BlockVectors,NZ>>(nd_range<1>(sycl::range{size_t(m)}, sycl::range{size_t(m)}), [=](sycl::nd_item<1> item){
            auto tid = item.get_local_linear_id();
            sycl::group<1> cta = item.get_group();
            constexpr auto SN = BlockVectors*3;
            std::bitset<BlockVectors> converged;
            oneapi::dpl::uniform_real_distribution<T> distr(0.0, 1.0);            
            oneapi::dpl::minstd_rand engine(42, tid);

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

            global_ptr<T> blockX = S_acc.get_pointer();
            global_ptr<T> blockR = blockX + BlockVectors*m;
            global_ptr<T> blockP = blockR + BlockVectors*m;

            global_ptr<T> blockXtemp = blockP + BlockVectors*m;
            global_ptr<T> blockRtemp = blockXtemp + BlockVectors*m;
            global_ptr<T> blockPtemp = blockRtemp + BlockVectors*m;

            global_ptr<T> blockAX = blockPtemp + BlockVectors*m;
            global_ptr<T> blockAR = blockAX + BlockVectors*m;
            global_ptr<T> blockAP = blockAR + BlockVectors*m;
            for (int i = tid + m*BlockVectors; i < 9*m*BlockVectors; i+=cta.get_local_range(0)){
                blockX[i] = T(0);
            }
            for (int i = tid; i < m*BlockVectors; i+=cta.get_local_range(0)){
                blockX[i] = distr(engine);
                X0_acc[i] = blockX[i];
            }
            group_barrier(cta);
            //Normalize S vectors
            orthonormalize<T, BlockVectors>(cta, blockX, m); //Modified Gram-Schmidt
            

            //A * X
            //matBlockVector<T, K, 1, N>(cta, A_tid, cols_tid, blockX, blockAX, NZ, m); 
            //X^T * A * X
            STAS<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockX, m, NZ, StAS.get_pointer(), Scache.get_pointer());
            //Fill in with simple 3x3 matrix

            compute_k_eigenpairs<T, BlockVectors, BlockVectors>(cta, StAS, eigvects, lambdas, working_space.get_pointer(), sort_scratchpad, true);
            //Print Eigenvectors
            for(int i = 0; i < BlockVectors; i++){
                for(int j = 0; j < BlockVectors; j++){
                    if(tid == 0) sycl::ext::oneapi::experimental::printf("Eigenvector[%d][%d] = %f\n", i, j, eigvects[i*BlockVectors + j]);
                }
            }

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
                for(int i = 0; i < BlockVectors; i++) {if(converged[i]) continue; converged[i] = sqrt(sycl::abs(reduce_over_group(cta, blockR[i*m + tid]*blockR[i*m + tid], sycl::plus<T>{}))) < tol;}
                if(tid == 0) sycl::ext::oneapi::experimental::printf("Iteration %d\n", iter);
                for(int i = 0; i < BlockVectors; i++){if(tid == 0) sycl::ext::oneapi::experimental::printf("Unorthogonalized BlockR[%d][%d] = %f\n", i, 0, blockR[i*m + 0]);}
                //Print ResidualNorms
                /* for(int i = 0; i < BlockVectors; i++){
                    T residual = sqrt(abs(reduce_over_group(cta, blockR[i*m + tid]*blockR[i*m + tid], sycl::plus<T>{})));
                    if(tid == 0) sycl::ext::oneapi::experimental::printf("Residual Norm %d = %f\n", i, residual);
                } */
                
                
                //Print blockR
                sycl::group_barrier(cta);
                //R = R - X * (X^T * R)
                applyConstraints<T, BlockVectors>(cta, blockR, blockX, m);
                sycl::group_barrier(cta);
                orthonormalize<T, BlockVectors>(cta, blockR, m);
                if (tid == 0) sycl::ext::oneapi::experimental::printf("Orthogonalized BlockR After Iteration %d\n", iter);
                for(int i = 0; i < BlockVectors; i++){if(tid == 0) sycl::ext::oneapi::experimental::printf("BlockR[%d][%d] = %f\n", i, 0, blockR[i*m + 0]);}
                
                matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockR, blockAR, NZ, m);
                

                if (!restart) {
                    orthonormalize<T, BlockVectors>(cta, blockP, m);
                    for(int i = 0; i < BlockVectors; i++){ 
                        blockXtemp[i*m + tid] = blockX[i*m + tid];
                        blockRtemp[i*m + tid] = blockR[i*m + tid];
                        blockPtemp[i*m + tid] = blockP[i*m + tid];
                    }

                    matBlockVector<T, K, 1, BlockVectors>(cta, A_tid, cols_tid, blockP, blockAP, NZ, m);

                    orthonormalize<T, BlockVectors*3>(cta, blockX, m);
                    STAS<T, K, 1, BlockVectors*3>(cta, A_tid, cols_tid, blockX, m, NZ, StAS.get_pointer(), Scache.get_pointer());
                    /* if (tid == 0) sycl::ext::oneapi::experimental::printf("Subspace Used In Iteration %d\n", iter);
                    for(int i = 0; i < BlockVectors*3; i++){if(tid == 0) sycl::ext::oneapi::experimental::printf("Subspace[%d][%d] = %f\n", i, 0, blockX[i*m + 0]);}
                     */
                    compute_k_eigenpairs<T, BlockVectors*3, BlockVectors>(cta, StAS, eigvects, lambdas, working_space.get_pointer(), sort_scratchpad, true);
                    
                    for (int i = 0; i < BlockVectors; i++){ 
                        blockX[i*m + tid] = blockXtemp[i*m + tid];
                        blockR[i*m + tid] = blockRtemp[i*m + tid];
                        blockP[i*m + tid] = blockPtemp[i*m + tid];
                    } 
                    group_barrier(cta);
                    update_vectors<T, BlockVectors*3, BlockVectors>(cta, blockX, blockR, blockP, blockAX, blockAR, blockAP, eigvects, m, false);
                    /* if (tid == 0) sycl::ext::oneapi::experimental::printf("BlockAX After Iteration %d\n", iter);
                    for(int i = 0; i < BlockVectors; i++){if(tid == 0) sycl::ext::oneapi::experimental::printf("blockAX[%d][%d] = %f\n", i, 0, blockAX[i*m + 0]);}
                    if (tid == 0) sycl::ext::oneapi::experimental::printf("BlockAP After Iteration %d\n", iter);
                    for(int i = 0; i < BlockVectors; i++){if(tid == 0) sycl::ext::oneapi::experimental::printf("blockAP[%d][%d] = %f\n", i, 0, blockAP[i*m + 0]);} */
                    group_barrier(cta);

                } else {
                    for (int i = 0; i < BlockVectors; i++){
                        blockXtemp[i*m + tid] = blockX[i*m + tid];
                        blockRtemp[i*m + tid] = blockR[i*m + tid];
                    }
                    orthonormalize<T, BlockVectors*2>(cta, blockX, m);
                    if (tid == 0) sycl::ext::oneapi::experimental::printf("Subspace Used In Iteration %d\n", iter);
                    for(int i = 0; i < 2*BlockVectors; i++){if(tid == 0) sycl::ext::oneapi::experimental::printf("Subspace[%d][%d] = %f\n", i, 0, blockX[i*m + 0]);}

                    STAS<T, K, 1, BlockVectors*2>(cta, A_tid, cols_tid, blockX, m, NZ, StAS.get_pointer(), Scache.get_pointer());
                    compute_k_eigenpairs<T, BlockVectors*2, BlockVectors>(cta, StAS, eigvects, lambdas, working_space.get_pointer(), sort_scratchpad, true);
                    /* if (tid == 0) sycl::ext::oneapi::experimental::printf("Eigenvectors Iteration %d\n", iter);
                    for(int i = 0; i < BlockVectors; i++){
                        for(int j = 0; j < BlockVectors*2; j++){
                            if(tid == 0) sycl::ext::oneapi::experimental::printf("eigvects[%d][%d] = %f\n", i, j, eigvects[i*BlockVectors*2 + j]);
                        }
                    } */
                    for (int i = 0; i < BlockVectors; i++){
                        blockX[i*m + tid] = blockXtemp[i*m + tid];
                        blockR[i*m + tid] = blockRtemp[i*m + tid];
                    }
                    group_barrier(cta);
                    update_vectors<T, BlockVectors*2, BlockVectors>(cta, blockX, blockR, blockP, blockAX, blockAR, blockAP, eigvects, m, true);
                    restart = false;
                }
                //Print eigenvectors
                sycl::group_barrier(cta);

                iter++;
            }
            for (int i = tid; i < BlockVectors*BlockVectors*3; i+=cta.get_local_range(0)){
                LastEigVects_acc[i] = eigvects[i];
            }
            for (int i = tid; i < 3*BlockVectors*3*BlockVectors; i+=cta.get_local_range(0)){
                LastGram_acc[i] = StAS[i];
            }

        });
    });
    std::string dtype = std::is_same<T, float>::value ? ".float32" : ".float64";
    ctx.wait();
    int SD = std::min(int(maxiters+1), 3);
    std::ofstream file("Subspace_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        auto S_acc = host_accessor(S, read_only);
        file.write(reinterpret_cast<const char*>(S_acc.get_pointer()), SD * BlockVectors * m * sizeof(T));
    }

    std::ofstream file2("AXARAP_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        auto S_acc = host_accessor(S, read_only);
        file2.write(reinterpret_cast<const char*>(S_acc.get_pointer() + 6 * BlockVectors * m), 3 * BlockVectors * m * sizeof(T));
    }

    std::ofstream file3("Initial_X_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + dtype, ios::binary);
    {
        auto X0_acc = host_accessor(X0, read_only);
        file3.write(reinterpret_cast<const char*>(X0_acc.get_pointer()), BlockVectors * m * sizeof(T));
    }

    std::ofstream file4("LastEigVects_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        auto LastEigVects_acc = host_accessor(LastEigVects, read_only);
        file4.write(reinterpret_cast<const char*>(LastEigVects_acc.get_pointer()), BlockVectors * BlockVectors * SD * sizeof(T));
    }

    std::ofstream file5("LastGram_BV=" + to_string(BlockVectors) + "_N=" + to_string(m) + "_iters=" + to_string(maxiters) + dtype, ios::binary);
    {
        auto LastGram_acc = host_accessor(LastGram, read_only);
        file5.write(reinterpret_cast<const char*>(LastGram_acc.get_pointer()), SD * BlockVectors * SD * BlockVectors * sizeof(T));
    }
}


int main(int argc, char** argv){

    CmdArgs args;
    parseArguments(argc, argv, args);
    size_t N = args.natoms;
    size_t BatchSize = args.nisomers;
    std::string device_type = args.device_type;
    size_t maxiter = args.nlanczos;
    auto queue = SyclQueue(device_type);

    FullereneBatch<float, uint16_t> batch(N, BatchSize);
    DualizeFunctor<float, uint16_t> dualize;
    TutteFunctor<float, uint16_t> tutte_layout;
    SphericalProjectionFunctor<float, uint16_t> spherical_projection;
    ForcefieldOptimizeFunctor<PEDERSEN, float, uint16_t> forcefield_optimize;
    ForcefieldOptimizeFunctor<PEDERSEN, double, uint16_t> forcefield_optimize_double;

    HessianFunctor<PEDERSEN, float, uint16_t> compute_hessians;
    HessianFunctor<PEDERSEN, double, uint16_t> compute_hessians_double;
    fill(batch);
    dualize(queue, batch, LaunchPolicy::SYNC);
    tutte_layout(queue, batch, LaunchPolicy::SYNC);
    spherical_projection(queue, batch, LaunchPolicy::SYNC);
    FullereneBatch<double, uint16_t> batch_double(N, BatchSize);
    {
        auto batch_acc_X = batch.d_.X_cubic_;
        auto batch_double_acc_X = batch_double.d_.X_cubic_;
        auto batch_acc_cubic_neighbours = batch.d_.A_cubic_;
        auto batch_double_acc_cubic_neighbours = batch_double.d_.A_cubic_;

        for(int i = 0; i < N*BatchSize; i++){
            for(int j = 0; j < 3; j++){
                batch_double_acc_X[i][j] = batch_acc_X[i][j];
                batch_double_acc_cubic_neighbours[i][j] = batch_acc_cubic_neighbours[i][j];
            }
        }
    }
    forcefield_optimize_double(queue, batch_double, LaunchPolicy::SYNC, 5*N, 5*N);
    forcefield_optimize(queue, batch, LaunchPolicy::SYNC, 5*N, 5*N);
    

    SyclVector<float> hessians((N*90*BatchSize));
    SyclVector<double> hessians_double((N*90*BatchSize));
    SyclVector<uint16_t> cols((N*90*BatchSize));

    compute_hessians_double(queue, batch_double, LaunchPolicy::SYNC, hessians_double, cols);
    compute_hessians(queue, batch, LaunchPolicy::SYNC, hessians, cols);



    int m = 3;
    int n = 3;
    //std::vector<float> A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    //std::vector<int> cols = {0, 1, 2, 1, 2, 0, 2, 0, 1};

    
    LOBPCG<float, uint16_t, 3, 30>(queue, hessians, cols, N*3, maxiter);
    LOBPCG<double, uint16_t, 3, 30>(queue, hessians_double, cols, N*3, maxiter);
    return 0;
}

