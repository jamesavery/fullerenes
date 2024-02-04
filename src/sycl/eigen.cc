using namespace sycl;
#define SQRT sycl::sqrt
#include <fullerenes/sycl-isomer-batch.hh>
#include "forcefield-includes.cc"
#include <array>
enum class EigensolveMode {NO_VECTORS, VECTORS, ENDS, FULL_SPECTRUM, FULL_SPECTRUM_MOLECULE};

template <typename T>
void apply_all_reflections(const T *V, const int n, const int m, T* Q)
{   
    static_assert(std::is_floating_point<T>::value, "T must be floating point");
    for(int k=0;k<n;k++){
        const T &v0 = V[2*k], &v1 = V[2*k+1];      
        // Udrullet:
        //       apply_reflection(Q({k,k+2},{0,m}), v);
        for(int l=threadIdx.x;l<m; l+=blockDim.x){
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
    FLOAT_TYPEDEFS(T);
    //  QTQ_calls ++;
    // Unrolled
    //  real_t numerical_zero = T.max_norm()*10*std::numeric_limits<real_t>::epsilon();
    // specialized max_norm = max(sum(abs(A),axis=1)) for tridiagonal matrix. 
    real_t local_max = real_t(0.);
    for (int i = tix; i < n; i += blockDim.x){
        local_max = std::max(local_max, abs(D[i]) + 2*abs(L[i]));
    }
    real_t max_norm = reduce_over_group(cta, local_max, sycl::greater<real_t>());
    real_t numerical_zero = 10*std::numeric_limits<real_t>::epsilon();
    device_real_t d_n, l_n, l_nm1;
    d_n = D[n]; l_n = L[n]; l_nm1 = L[n-1];
    sycl::group_barrier(cta);
    //real_t a[2], v[2], D[n+1], L[n+1], U[2*(n+1)];
    real_t a[2], v[2];//, D[n+1], L[n+1], U[2*(n+1)];
    for(int k = tix; k < n + 1; k += blockDim.x){
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
            
            real_t anorm = SQRT(a[0]*a[0] + a[1]*a[1]); 

            // // Udrullet
            // //    reflection_vector(a,anorm,v);
            v[0] = D[k]; v[1] = L[k];
            real_t alpha = -copysign(anorm,a[0]); // Koster ingenting
            v[0] -= alpha;

            real_t vnorm = SQRT(v[0]*v[0]+v[1]*v[1]);
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
    for (int k = tix; k<n; k+=blockDim.x){  // Copy working diagonals to output
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
std::array<T,2> INLINE eigvalsh2x2(const std::array<T,4> &A){
    auto [a,b,c,d] = A;
    T D = SQRT(4*b*c+(a-d)*(a-d));
    return {(a+d-D)/2, (a+d+D)/2};
}

template <EigensolveMode mode, typename T, typename K, int M>
void lanczos_(sycl::queue& Q, const IsomerBatch<T,K>& B, buffer<T>& V_, buffer<T>& U, buffer<T>& D, const buffer<T>& H, const buffer<K>& cols, int m){
    TEMPLATE_TYPEDEFS(T,K);
    auto Natoms = B.N();
    auto capacity = B.capacity();

    
}


template <EigensolveMode mode, Device DEV,typename T, typename K>
void compute_eigenvectors_(sycl::queue& ctx, const IsomerBatch<T,K>& B, const sycl::buffer<T>& Q, sycl::buffer<T>& V, sycl::buffer<T>& E, int m){
    int atom_idx = threadIdx.x/3; //Atom index (Integer division so the result is rounded down)
    TEMPLATE_TYPEDEFS(T,K);
    SMEM(T);
    int n = B.n_atoms * 3;
    int offset = 0;
    if (mode == EigensolveMode::ENDS || mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
        offset = 6;
    }
    for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x){
        real_t* v = V.data + I * m * n;
        real_t* e = E.data + I * n * n;
        real_t* q = Q.data + I * m * m;
        if (mode == EigensolveMode::FULL_SPECTRUM_MOLECULE) {
            for(int i = 0; i < 3; i++){
                e[i*n + threadIdx.x] = real_t(threadIdx.x%3 == i)/SQRT(B.n_atoms); 
            }
            float* X_ptr = B.X + n*I; //WARNING float hardcoded here to match the type of B.X (IsomerBatch)
            if (threadIdx.x%3 == 0) {
                e[3*n + threadIdx.x] = real_t(0.);
                e[4*n + threadIdx.x] = -X_ptr[atom_idx*3 + 2];
                e[5*n + threadIdx.x] = -X_ptr[atom_idx*3 + 1];
            } else if (threadIdx.x%3 == 1) {
                e[3*n + threadIdx.x] = -X_ptr[atom_idx*3 + 2];
                e[4*n + threadIdx.x] = real_t(0.);
                e[5*n + threadIdx.x] = X_ptr[atom_idx*3 + 0];
            } else {
                e[3*n + threadIdx.x] = X_ptr[atom_idx*3 + 1];
                e[4*n + threadIdx.x] = X_ptr[atom_idx*3 + 0];
                e[5*n + threadIdx.x] = real_t(0.);
            }
                clear_cache(smem,  n + warpSize);
                e[3*n + threadIdx.x] /= SQRT(reduction(smem, e[3*n + threadIdx.x]*e[3*n + threadIdx.x]));
                clear_cache(smem,  n + warpSize);
                e[4*n + threadIdx.x] /= SQRT(reduction(smem, e[4*n + threadIdx.x]*e[4*n + threadIdx.x]));
                clear_cache(smem,  n + warpSize);
                e[5*n + threadIdx.x] /= SQRT(reduction(smem, e[5*n + threadIdx.x]*e[5*n + threadIdx.x]));
        }
        
        if(threadIdx.x < n)
        for (int k = offset; k < n; k++){
            e = E.data + I * n * n + k * n;
            e[threadIdx.x] = real_t(0.);
            q = Q.data + I * m * m + (k-offset) * m;
            for (int i = 0; i < m; i++){
                e[threadIdx.x] += v[i*n + threadIdx.x] * q[i];
            }
        }
    }
}


//template <typename T, typename K>
//void eigensolve(const IsomerBatch<T,K>& B, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy){
std::vector<sycl::device> get_devices(){
    auto platforms = sycl::platform::get_platforms();
    std::vector<sycl::device> devices;
    for(auto& platform : platforms){
        auto platform_devices = platform.get_devices();
        for(auto& device : platform_devices){
            devices.push_back(device);
        }
    }
    return devices;
}

int get_device_index(std::vector<sycl::device>& devices, sycl::device& device){
    for(int i = 0; i < devices.size(); i++){
        if(devices[i] == device){
            return i;
        }
    }
    return -1;
}


template <EigensolveMode mode, typename T, typename K>
void eigensolve(sycl::queue& ctx, const IsomerBatch<T,K> B, sycl::buffer<T>& Q, const sycl::buffer<T>& hessians, const sycl::buffer<K>& cols, const sycl::buffer<T>& eigenvalues){
    TEMPLATE_TYPEDEFS(T,K);
    static size_t Natoms = B.N();
    static size_t capacity = B.capacity();
    
    static auto devices = std::vector<sycl::device>             {ctx.get_device()};
    static std::vector<sycl::buffer<T>> offDiagBuffers          {B.N()*B.capacity()*3};         //Off-diagonals
    static std::vector<sycl::buffer<T>> QmatBuffers             {B.N()*B.N()*B.capacity()*3*3}; //Q matrices
    static std::vector<sycl::buffer<T>> lanczosVectorBuffers    {B.N()*B.N()*B.capacity()*3*3}; //Lanczos vectors
    auto index = get_device_index(devices, ctx.get_device());

    if(index == -1){
        devices.push_back(ctx.get_device());
        index = devices.size() - 1;
        offDiagBuffers.         push_back(sycl::buffer<T>(B.N()*B.capacity()*3));
        QmatBuffers.            push_back(sycl::buffer<T>(B.N()*B.N()*B.capacity()*3*3));
        lanczosVectorBuffers.   push_back(sycl::buffer<T>(B.N()*B.N()*B.capacity()*3*3));
    }
    if(Natoms != B.N() || capacity != B.capacity()){
        Natoms = B.N();
        capacity = B.capacity();
        for (int i = 0; i < devices.size(); i++){
            offDiagBuffers[i]       = sycl::buffer<T>(B.N()*B.capacity()*3);
            QmatBuffers[i]          = sycl::buffer<T>(B.N()*B.N()*B.capacity()*3*3);
            lanczosVectorBuffers[i] = sycl::buffer<T>(B.N()*B.N()*B.capacity()*3*3);
        }
    }

    //Buffers required: buffer<T>& D_, buffer<T>& L_, buffer<T>& U_, buffer<T>& Q_, 

    //Lanczos
    //CUDA Lanczos Call:        void* kargs[]{(void*)&B,    (void*)&Vs[dev],                (void*)&Ls[dev],        (void*)&eigenvalues,    (void*)&hessians,   (void*)&cols,           (void*)&Nlanczos};
    //CUDA Lanczos Signature:   (const IsomerBatch<DEV> B,  CuArray<T> V_,                  CuArray<T> U,           CuArray<T> D,           const CuArray<T> H, const CuArray<K> cols,  int m)
    //Should be:                B,                          lanczosVectorBuffers[index],    offDiagBuffers[index],  eigenvalues, hessians, cols, n

    //Us is never so we can just delete it
    //Diagonalization (QR) Call:        void* kargs_qr[]{(void*)&B, (void*)&eigenvalues,    (void*)&Ls[dev],    (void*)&Us[dev],    (void*)&Qs[dev],    (void*)&Nlanczos};
    //Diagonalization (QR) Signature:   (const IsomerBatch<DEV> B,  CuArray<T> D_,          CuArray<T> L_,      CuArray<T> U_,      CuArray<T> Q_,      int n)
    //Should be:                        B,                          eigenvalues,            L_,                 NOTHING,                Q_,                n

    Q.submit([&](sycl::handler& h){
        accessor V_acc(lanczosVectorBuffers[index], h, write);
        accessor D_acc(eigenvalues, h, read_write);
        accessor U_acc(offDiagBuffers[index], h, read_write);
        accessor H_acc(hessians, h, read_only);
        accessor cols_acc(cols, h, read_only);
        accessor X_acc(B.X, h, read_only);
        local_accessor smem(range<1>(Natoms*3), h);
        local_accessor betas(range<1>(Natoms*3), h);
        local_accessor alphas(range<1>(Natoms*3), h);
        h.parallel_for<class lanczos_kernel>(sycl::nd_range(sycl::range{Natoms*3*capacity}, sycl::range{Natoms*3}), [=](nd_item<1> nditem){
            auto cta = nditem.get_group();
            auto tid = nditem.get_local_linear_id();
            auto bid = nditem.get_group_linear_id();
            auto N = nditem.get_global_range()[0];
            auto atom_idx = tid/3; //Atom index
            real_t A[M]; //Hessian matrix, threadIdx.x'th row
            node_t C[M]; //Column indices of the threadIdx.x'th row 3-fold degenerate
            real_t* V = V_acc.get_pointer() + bid * m * N + tid;
            real_t* X_ptr = X_acc.get_pointer() + N*bid; 

            auto LCG = [&](const unsigned long seed){
                //Parameters from Numerical Recipes, Knuth and H. W. Lewis
                unsigned long a = 1664525;
                unsigned long c = 1013904223;
                unsigned long m = 4294967296;
                unsigned long x = seed;
                return (a*x + c) % m;
            };

            auto mat_vect = [&](const real_t x){
                real_t result = real_t(0);
                smem[tid] = x;
                sycl::group_barrier(cta);
                #pragma unroll
                for (int j = 0; j < M; j++){
                    int col = C[j];
                    result += A[j] * smem[col];
                }
                sycl::group_barrier(cta);
                return result;
            };

            //Modified Gram-Schmidt, Also orthogonalizes against the deflation space
            auto MGS = [&](int index){
                sycl::group_barrier(cta);
                real_t result = V[index*N];
                if (mode == EigensolveMode::ENDS || mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
                    #pragma unroll
                    for (int j = 0; j < 6; j++){
                        auto proj = reduce_over_group(cta, result * Z[j]) * Z[j];
                        result -= proj; //Remove the component along Z[j] from result
                    }
                }

                #pragma unroll
                for (int j = 0; j < index; j++){
                    auto proj = reduce_over_group(cta, result * V[j*N]) * V[j*N];
                    result -= proj; //Remove the component along V[j*N] from result
                }
                result /= sqrt(reduce_over_group(cta, result * result));
                return result;
            };

            real_t Z[6]; //Eigenvectors spanning the kernel of the hessian (Rotations, Translations)
            if (mode == EigensolveMode::ENDS || mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
                Z[0] = real_t(tid%3 == 0)/SQRT(Natoms); Z[1] = real_t(tid%3 == 1)/SQRT(Natoms); Z[2] = real_t(tid%3 == 2)/SQRT(Natoms); // Translation eigenvectors
                if(tid%3 == 0){
                    Z[3] = real_t(0.);
                    Z[4] = -X_ptr[atom_idx*3 + 2];
                    Z[5] = -X_ptr[atom_idx*3 + 1];
                }
                if(tid%3 == 1){
                    Z[3] = -X_ptr[atom_idx*3 + 2];
                    Z[4] = real_t(0.);
                    Z[5] = X_ptr[atom_idx*3 + 0];
                }
                if(tid%3 == 2){
                    Z[3] = X_ptr[atom_idx*3 + 1];
                    Z[4] = X_ptr[atom_idx*3 + 0];
                    Z[5] = real_t(0.);
                }
                Z[3] /= sycl::sqrt(reduce_over_group(cta, Z[3]*Z[3]));
                Z[4] /= sycl::sqrt(reduce_over_group(cta, Z[4]*Z[4]));
                Z[5] /= sycl::sqrt(reduce_over_group(cta, Z[5]*Z[5]));
            }

            sycl::group_barrier(cta);

            for(int j = 0; j < M; j++){
                A[j] = H_acc[bid*N*M + tid*M + j];
                C[j] = cols_acc[bid*N*M + tid*M + j];
            }
            sycl::group_barrier(cta);

            #include <chrono>

            // Generate a random number between 0 and 99
            V[0*N] = real_t(LCG(tid)); //Seed the random number generator with the thread id
            V[0*N] /= sycl::sqrt(reduce_over_group(cta, V[0*N] * V[0*N]));
            V[0*N] = MGS(0);
            for (int i = 0; i < m; i++){
                if (i % 2 == 0 && i > 1){
                    V[(i-1)*N] = MGS(i-1);
                    V[i*N] = MGS(i);
                }
                real_t v = mat_vect(V[i*N]);
                real_t alpha = reduce_over_group(cta, v * V[i*N]);
                if (tid == i) alphas[i] = alpha;
                if (i == 0){
                    v -= alpha * V[i*N];
                } else {
                    v -= betas[i- 1] * V[(i-1)*N] + alpha * V[i*N];
                }
                real_t beta = SQRT(reduce_over_group(cta, v * v));
                if (tid == i) betas[i] = beta;
                if (i < m-1) V[(i+1)*N] = v / beta;
            }
            if (tid < m){
                D_acc[bid*m + tid] = alphas[tid];
                U_acc[bid*m + tid] = betas[tid];
            }
        });
    });
    
    //Diagonalization (QR) Call:        void* kargs_qr[]{(void*)&B, (void*)&eigenvalues,    (void*)&Ls[dev],    (void*)&Us[dev],    (void*)&Qs[dev],    (void*)&Nlanczos};
    //Diagonalization (QR) Signature:   (const IsomerBatch<DEV> B,  CuArray<T> D_,          CuArray<T> L_,      CuArray<T> U_,      CuArray<T> Q_,      int n)
    //Should be:                        B,                          eigenvalues,            offDiagBuffers[index],  NOTHING,                QmatBuffers[index],                n

    Q.submit([&](sycl::handler& h){
        local_accessor D(range<1>(Natoms*3), h); //Diagonal
        local_accessor L(range<1>(Natoms*3), h); //Lower subdiagonal
        local_accessor U(range<1>(Natoms*3), h); //Upper subdiagonal
        local_accessor V(range<1>(Natoms*3), h); 
        
        accessor D_acc(eigenvalues, h, read_write);
        accessor L_acc(offDiagBuffers[index], h, read_write);
        accessor Q_acc(QmatBuffers[index], h, read_write);

        h.parallel_for<class eigensolve>(sycl::nd_range(sycl::range{capacity*64}, sycl::range{64}), [=](sycl::nd_item<1> nditem){
            auto tid = nditem.get_local_linear_id();
            auto bid = nditem.get_group_linear_id();
            auto cta = nditem.get_group();
            auto bdim = nditem.get_local_range();

            for (int i = tid; i < n; i += bdim){
                D[i]        = D_acc[bid*n + i];
                L[i]        = L_acc[bid*n + i];
                U[i]        = L_acc[bid*n + i];
                Q_acc[bid*n*n + i*(n+1)] = 1; //Identity matrix
            }
            sycl::group_barrier(cta);
            for(int k=n-1;k>=0;k--){
                real_t d = D[k];
                real_t shift = d;

                int i = 0;
                real_t GR = (k>0 ? sycl::abs(L[k-1]) : 0) + sycl::abs(L[k]);
                int not_done = 1;
                while (not_done > 0){
                    i++;
                    T_QTQ(cta, k+1, D.get_pointer(), L.get_pointer(), U.get_pointer(), V.get_pointer(), shift);
                    if(mode == EigensolveMode::VECTORS || mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
                        apply_all_reflections(V.get_pointer(),k,n,Q_acc.get_pointer() + n*n*bid);
                    }
                    GR = (k>0? sycl::abs(L[k-1]) : 0) + (k+1<n? sycl::abs(L[k]): 0);

                    if(k>0){
                        std::array<T,4> args = {D[k-1], L[k-1], L[k-1], D[k]};
                        auto [l0, l1] = eigvalsh2x2(args);
                        shift = sycl::abs(l0-d) < sycl::abs(l1-d)? l0:l1;
                    } else shift = D[k];

                    shift    = D[k];
                
                    if(GR <= std::numeric_limits<real_t>::epsilon()*real_t(10.)) not_done--; // Do one (or optionally more) steps after reaching tolerance, to get all off-diagonals below.
                                                    // GPU NB: Se GPU NB ovenfor.
                    if(i>5){
                        //printf("%dth run: Cannot converge eigenvalue %d to tolerance " G " using machine precision %g (d=%g, shift=%g, G=%g)\n" "D[k] = %g, L[k-1] = %g, L[k] = %g\n", nth_time,k,tolerance, std::numeric_limits<real_t>::epsilon(),d,shift,GR, D[k], (k>0)?L[k-1]:0, (k+1<n)?L[k]:0);
                        auto max_error = std::max(std::numeric_limits<real_t>::epsilon()*real_t(10.),GR);
                        break;
                    }
                }
                sycl::group_barrier(cta);
                //Copy back to global memory
                for (int i = tid; i < n; i += bdim){
                    if( mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
                        if(i < 6) {D_acc[(n+6)*bid + i] = 0;}
                        D_acc[(n+6)*bid + 6 + i] = D[i];
                    } else {
                        D_acc[n*bid + i] = D[i];
                    }
                }
            }
        });
    });

}

 
