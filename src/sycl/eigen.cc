#include <fullerenes/sycl-isomer-batch.hh>
#include "forcefield-includes.cc"
#include <array>
using namespace sycl;
#define SQRT sycl::sqrt

template <typename T>
void apply_all_reflections(const sycl::group<1> &cta, const sycl::local_accessor<T,1>& V, const int n, const int m, multi_ptr<float, sycl::access::address_space::global_space, sycl::access::decorated::legacy> Q)
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
void T_QTQ(sycl::group<1>& cta, const int n, const sycl::local_accessor<T,1>& D, const sycl::local_accessor<T,1>& L, const sycl::local_accessor<T,1>& U, const sycl::local_accessor<T,1>& Vout, T shift=0){
    int tix = cta.get_local_linear_id();
    int bdim = cta.get_local_linear_range();
    FLOAT_TYPEDEFS(T);
    //  QTQ_calls ++;
    // Unrolled
    //  real_t numerical_zero = T.max_norm()*10*std::numeric_limits<real_t>::epsilon();
    // specialized max_norm = max(sum(abs(A),axis=1)) for tridiagonal matrix. 
    real_t local_max = real_t(0.);
    for (int i = tix; i < n; i += bdim){
        local_max = std::max(local_max, abs(D[i]) + 2*abs(L[i]));
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
            if (ABS(L[k]) > numerical_zero){
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
    T D = SQRT(4*b*c+(a-d)*(a-d));
    return {(a+d-D)/2, (a+d+D)/2};
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


template <typename T>
struct EigenBuffers{
    std::vector<sycl::device> devices;

    std::vector<sycl::buffer<T,1>> offDiagBuffers;
    std::vector<sycl::buffer<T,1>> QmatBuffers;
    std::vector<sycl::buffer<T,1>> lanczosBuffers;
    std::vector<std::array<size_t, 3>> dims;

    //SYCL lacks a find for retrieving a unique identifier for a given queue, platform, or device so we have to do this manually
    int get_device_index(const sycl::device& device){
        for(int i = 0; i < devices.size(); i++){
            if(devices[i] == device){
                return i;
            }
        }
        return -1;
    }

    EigenBuffers(size_t N, size_t capacity, size_t nLanczos, const sycl::device& device){
        devices = {device};
        offDiagBuffers = std::vector<sycl::buffer<T,1>>(1, sycl::buffer<T,1>(nLanczos*capacity));
        QmatBuffers = std::vector<sycl::buffer<T,1>>(1, sycl::buffer<T,1>(nLanczos*nLanczos*capacity));
        lanczosBuffers = std::vector<sycl::buffer<T,1>>(1, sycl::buffer<T,1>(nLanczos*N*3*capacity));
        dims = std::vector<std::array<size_t, 3>>(1, {N, capacity, nLanczos});
    }

    void resize(size_t N, size_t capacity, size_t nLanczos, size_t idx){
        offDiagBuffers[idx]   = sycl::buffer<T,1>(nLanczos*capacity);
        QmatBuffers[idx]      = sycl::buffer<T,1>(nLanczos*nLanczos*capacity);
        lanczosBuffers[idx]          = sycl::buffer<T,1>(nLanczos*N*3*capacity);
        dims[idx]             = {N, capacity, nLanczos};
    }
    
    void append_buffers(size_t N, size_t capacity, size_t nLanczos, const sycl::device& device){
        offDiagBuffers.push_back(sycl::buffer<T,1>(nLanczos*capacity));
        QmatBuffers.push_back(sycl::buffer<T,1>(nLanczos*nLanczos*capacity));
        lanczosBuffers.push_back(sycl::buffer<T,1>(nLanczos*N*3*capacity));
        dims.push_back({N, capacity, nLanczos});
        devices.push_back(device);
    }

    void update(size_t N, size_t capacity, size_t nLanczos, const sycl::device& device){
        auto idx = get_device_index(device);
        if(idx == -1){
            append_buffers(N, capacity, nLanczos, device);
            std::cout << "Appended buffers for device " << device.get_info<sycl::info::device::name>() << " ID: " << devices.size() - 1 << "\n";
        }
        else if (N != dims[idx][0] || capacity != dims[idx][1] || nLanczos != dims[idx][2]){
            resize(N, capacity, nLanczos, idx);
            std::cout << "Resized buffers for device " << device.get_info<sycl::info::device::name>() << " ID: " << idx << "\n";
        }
    }

};

template <EigensolveMode mode, typename T, typename K>
void eigensolve(sycl::queue& ctx, IsomerBatch<T,K> B, sycl::buffer<T,1>& hessians, sycl::buffer<K,1>& cols, sycl::buffer<T,1>& eigenvalues){
    TEMPLATE_TYPEDEFS(T,K);
    size_t nLanczos = B.N()*3, capacity = B.capacity(), Natoms = B.N();
    static EigenBuffers<T> buffers(Natoms, capacity, nLanczos, ctx.get_device());
    buffers.update(Natoms, capacity, nLanczos, ctx.get_device());
    auto index = buffers.get_device_index(ctx.get_device());

    //Buffers required: buffer<T>& D_, buffer<T>& L_, buffer<T>& U_, buffer<T>& Q_, 

    //Lanczos
    //CUDA Lanczos Call:        void* kargs[]{(void*)&B,    (void*)&Vs[dev],                (void*)&Ls[dev],        (void*)&eigenvalues,    (void*)&hessians,   (void*)&cols,           (void*)&Nlanczos};
    //CUDA Lanczos Signature:   (const IsomerBatch<DEV> B,  CuArray<T> V_,                  CuArray<T> U,           CuArray<T> D,           const CuArray<T> H, const CuArray<K> cols,  int m)
    //Should be:                B,                          lanczosBuffers[index],    offDiagBuffers[index],  eigenvalues, hessians, cols, n

    //Us is never so we can just delete it
    //Diagonalization (QR) Call:        void* kargs_qr[]{(void*)&B, (void*)&eigenvalues,    (void*)&Ls[dev],    (void*)&Us[dev],    (void*)&Qs[dev],    (void*)&Nlanczos};
    //Diagonalization (QR) Signature:   (const IsomerBatch<DEV> B,  CuArray<T> D_,          CuArray<T> L_,      CuArray<T> U_,      CuArray<T> Q_,      int n)
    //Should be:                        B,                          eigenvalues,            L_,                 NOTHING,                Q_,                n

    ctx.submit([&](sycl::handler& h){
        accessor V_acc(buffers.lanczosBuffers[index], h, write_only);
        accessor D_acc(eigenvalues, h, read_write);
        accessor U_acc(buffers.offDiagBuffers[index], h, read_write);
        accessor H_acc(hessians, h, read_only);
        accessor cols_acc(cols, h, read_only);
        accessor X_acc(B.X, h, read_only);
        local_accessor<T,1> smem(range<1>(Natoms*3), h);
        local_accessor<T,1> betas(range<1>(Natoms*3), h);
        local_accessor<T,1> alphas(range<1>(Natoms*3), h);
        h.parallel_for<class lanczos_kernel>(sycl::nd_range(sycl::range{Natoms*3*capacity}, sycl::range{Natoms*3}), [=](nd_item<1> nditem){
            auto cta = nditem.get_group();
            auto tid = nditem.get_local_linear_id();
            auto bid = nditem.get_group_linear_id();
            auto atom_idx = tid/3; //Atom index
            auto N = Natoms*3;
            constexpr int M = 30;
            real_t A[M]; //Hessian matrix, threadIdx.x'th row
            node_t C[M]; //Column indices of the threadIdx.x'th row 3-fold degenerate
            real_t* V = V_acc.get_pointer() + bid * nLanczos * N + tid;
            coord3d* X_ptr = X_acc.get_pointer() + Natoms*bid; 

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

            real_t Z[6]; //Eigenvectors spanning the kernel of the hessian (Rotations, Translations)
            if (mode == EigensolveMode::ENDS || mode == EigensolveMode::SPECIAL){
                Z[0] = real_t(tid%3 == 0)/SQRT(T(Natoms)); Z[1] = real_t(tid%3 == 1)/SQRT(T(Natoms)); Z[2] = real_t(tid%3 == 2)/SQRT(T(Natoms)); // Translation eigenvectors
                if(tid%3 == 0){
                    Z[3] = real_t(0.);
                    Z[4] = -X_ptr[atom_idx][2];
                    Z[5] = -X_ptr[atom_idx][1];
                }
                if(tid%3 == 1){
                    Z[3] = -X_ptr[atom_idx][2];
                    Z[4] = real_t(0.);
                    Z[5] = X_ptr[atom_idx][0];
                }
                if(tid%3 == 2){
                    Z[3] = X_ptr[atom_idx][1];
                    Z[4] = X_ptr[atom_idx][0];
                    Z[5] = real_t(0.);
                }
                Z[3] /= sycl::sqrt(reduce_over_group(cta, Z[3]*Z[3], sycl::plus<real_t>{}));
                Z[4] /= sycl::sqrt(reduce_over_group(cta, Z[4]*Z[4], sycl::plus<real_t>{}));
                Z[5] /= sycl::sqrt(reduce_over_group(cta, Z[5]*Z[5], sycl::plus<real_t>{}));
            }

            //Modified Gram-Schmidt, Also orthogonalizes against the deflation space
            auto MGS = [&](int index){
                sycl::group_barrier(cta);
                real_t result = V[index*N];
                if (mode == EigensolveMode::ENDS || mode == EigensolveMode::SPECIAL){
                    #pragma unroll
                    for (int j = 0; j < 6; j++){
                        auto proj = reduce_over_group(cta, result * Z[j], sycl::plus<real_t>{}) * Z[j];
                        result -= proj; //Remove the component along Z[j] from result
                    }
                }

                #pragma unroll
                for (int j = 0; j < index; j++){
                    auto proj = reduce_over_group(cta, result * V[j*N], sycl::plus<real_t>{}) * V[j*N];
                    result -= proj; //Remove the component along V[j*N] from result
                }
                result /= sqrt(reduce_over_group(cta, result * result, sycl::plus<real_t>{}));
                return result;
            };

            sycl::group_barrier(cta);

            for(int j = 0; j < M; j++){
                A[j] = H_acc[bid*N*M + tid*M + j];
                C[j] = cols_acc[bid*N*M + tid*M + j];
            }
            sycl::group_barrier(cta);

            // Generate a random number between 0 and 99
            V[0*N] = real_t(LCG(tid)); //Seed the random number generator with the thread id
            V[0*N] /= sycl::sqrt(reduce_over_group(cta, V[0*N] * V[0*N], sycl::plus<real_t>{}));
            V[0*N] = MGS(0);
            for (int i = 0; i < nLanczos; i++){
                if (i % 2 == 0 && i > 1){
                    V[(i-1)*N] = MGS(i-1);
                    V[i*N] = MGS(i);
                }
                real_t v = mat_vect(V[i*N]);
                real_t alpha = reduce_over_group(cta, v * V[i*N],sycl::plus<real_t>{});
                if (tid == i) alphas[i] = alpha;
                if (i == 0){
                    v -= alpha * V[i*N];
                } else {
                    v -= betas[i- 1] * V[(i-1)*N] + alpha * V[i*N];
                }
                real_t beta = SQRT(reduce_over_group(cta, v * v, sycl::plus<real_t>{}));
                if (tid == i) betas[i] = beta;
                if (i < nLanczos-1) V[(i+1)*N] = v / beta;
            }
            if (tid < nLanczos){
                D_acc[bid*nLanczos + tid] = alphas[tid];
                U_acc[bid*nLanczos + tid] = betas[tid];
            }
        });
    });
    
    //Diagonalization (QR) Call:        void* kargs_qr[]{(void*)&B, (void*)&eigenvalues,    (void*)&Ls[dev],    (void*)&Us[dev],    (void*)&Qs[dev],    (void*)&Nlanczos};
    //Diagonalization (QR) Signature:   (const IsomerBatch<DEV> B,  CuArray<T> D_,          CuArray<T> L_,      CuArray<T> U_,      CuArray<T> Q_,      int n)
    //Should be:                        B,                          eigenvalues,            offDiagBuffers[index],  NOTHING,                QmatBuffers[index],                n

    ctx.submit([&](sycl::handler& h){
        local_accessor<T,1> D(range<1>(Natoms*3), h); //Diagonal
        local_accessor<T,1> L(range<1>(Natoms*3), h); //Lower subdiagonal
        local_accessor<T,1> U(range<1>(Natoms*3), h); //Upper subdiagonal
        local_accessor<T,1> V(range<1>(Natoms*7), h); 
        
        accessor D_acc(eigenvalues, h, read_write);
        accessor L_acc(buffers.offDiagBuffers[index], h, read_write);
        accessor Q_acc(buffers.QmatBuffers[index], h, read_write);

        h.parallel_for<class eigensolve>(sycl::nd_range(sycl::range{capacity*64}, sycl::range{64}), [=](sycl::nd_item<1> nditem){
            auto tid = nditem.get_local_linear_id();
            auto bid = nditem.get_group_linear_id();
            auto cta = nditem.get_group();
            auto bdim = cta.get_local_linear_range();
            auto n = nLanczos;

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
                real_t GR = (k>0?abs(L[k-1]):0)+abs(L[k]);
                int not_done = 1;
                while (not_done > 0){
                    i++;
                    T_QTQ(cta, k+1, D, L, U, V, shift);
                    if(mode == EigensolveMode::FULL_SPECTRUM_VECTORS || mode == EigensolveMode::ENDS_VECTORS){
                        apply_all_reflections(cta, V,k,n,Q_acc.get_pointer() + n*n*bid);
                    }
                    GR = (k>0?abs(L[k-1]):0)+(k+1<n?abs(L[k]):0);

                    if(k>0){
                        std::array<T,4> args = {D[k-1], L[k-1], L[k-1], D[k]};
                        auto [l0, l1] = eigvalsh2x2(args);
                        shift = sycl::abs(l0-d) < sycl::abs(l1-d)? l0:l1;
                    } else {shift = D[k];}
                
                    if(GR <= std::numeric_limits<real_t>::epsilon()*real_t(10.)) not_done--; // Do one (or optionally more) steps after reaching tolerance, to get all off-diagonals below.
                                                    // GPU NB: Se GPU NB ovenfor.
                    if(i>5){
                        //printf("%dth run: Cannot converge eigenvalue %d to tolerance " G " using machine precision %g (d=%g, shift=%g, G=%g)\n" "D[k] = %g, L[k-1] = %g, L[k] = %g\n", nth_time,k,tolerance, std::numeric_limits<real_t>::epsilon(),d,shift,GR, D[k], (k>0)?L[k-1]:0, (k+1<n)?L[k]:0);
                        auto max_error = std::max(std::numeric_limits<real_t>::epsilon()*real_t(10.),GR);
                        break;
                    }
                }
            }
            sycl::group_barrier(cta);
           
            //Copy back to global memory
            for (int i = tid; i < n; i += bdim){
                if( mode == EigensolveMode::SPECIAL){
                    if(i < 6) {D_acc[(n+6)*bid + i] = 0;}
                    D_acc[(n+6)*bid + 6 + i] = D[i];
                } else {
                    D_acc[n*bid + i] = D[i];
                }
            } 
        });
    });
    ctx.wait();
    std::cout << "Eigensolve complete\n";
}
template void eigensolve<EigensolveMode::FULL_SPECTRUM, float, uint16_t>(sycl::queue& ctx, const IsomerBatch<float,uint16_t> B, sycl::buffer<float,1>& hessians, sycl::buffer<uint16_t,1>& cols, sycl::buffer<float,1>& eigenvalues); 