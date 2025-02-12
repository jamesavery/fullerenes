#include <sycl/sycl.hpp>
#include <fullerenes/sycl-headers/sycl-span.hh>

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