#include <fullerenes/sycl-headers/eigen-kernels.hh>
#include "forcefield-includes.cc"
#include <array>
#include "sycl/sycl.hpp"
#include "queue-impl.cc"
#include <oneapi/dpl/random>
using namespace sycl;

template <typename T>
void apply_all_reflections(const sycl::group<1> &cta, const sycl::local_accessor<T,1>& V, const int n, const int m, T* Q)
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

template <EigensolveMode mode, typename T, typename K> class Lanczos {};
template <EigensolveMode mode, typename T, typename K> class QR {};
template <EigensolveMode mode, typename T, typename K> class Vectors {};
template <EigensolveMode mode, typename T, typename K> class Reduction {};

template <typename T, typename K>
struct EigenBuffers{
    std::vector<sycl::device> devices;

    std::vector<sycl::buffer<T,1>> offDiagBuffers;
    std::vector<sycl::buffer<T,1>> diagBuffers;
    std::vector<sycl::buffer<T,1>> QmatBuffers;
    std::vector<sycl::buffer<T,1>> lanczosBuffers;
    std::vector<sycl::buffer<K,1>> endsIdxBuffers; //Only used in ENDS and ENDS_VECTORS mode. Purpose is to store the indices of the maximum and minimum eigenvalues. 
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
        diagBuffers = std::vector<sycl::buffer<T,1>>(1, sycl::buffer<T,1>(nLanczos*capacity));
        endsIdxBuffers = std::vector<sycl::buffer<K,1>>(1, sycl::buffer<K,1>(2*capacity));
        dims = std::vector<std::array<size_t, 3>>(1, {N, capacity, nLanczos});
    }

    void resize(size_t N, size_t capacity, size_t nLanczos, size_t idx){
        offDiagBuffers[idx]   = sycl::buffer<T,1>(nLanczos*capacity);
        diagBuffers[idx]      = sycl::buffer<T,1>(nLanczos*capacity);
        QmatBuffers[idx]      = sycl::buffer<T,1>(nLanczos*nLanczos*capacity);
        lanczosBuffers[idx]   = sycl::buffer<T,1>(nLanczos*N*3*capacity);
        endsIdxBuffers[idx]   = sycl::buffer<K,1>(2*capacity);
        dims[idx]             = {N, capacity, nLanczos};
    }
    
    void append_buffers(size_t N, size_t capacity, size_t nLanczos, const sycl::device& device){
        offDiagBuffers.push_back(sycl::buffer<T,1>(nLanczos*capacity));
        diagBuffers.push_back(sycl::buffer<T,1>(nLanczos*capacity));
        QmatBuffers.push_back(sycl::buffer<T,1>(nLanczos*nLanczos*capacity));
        lanczosBuffers.push_back(sycl::buffer<T,1>(nLanczos*N*3*capacity));
        endsIdxBuffers.push_back(sycl::buffer<K,1>(2*capacity));
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


/* 
    * @brief: Multi-purpose eigensolver for the batch of Hessian matrices derived from the batch of isomers.
    * @param mode: EigensolveMode, the mode of the eigensolver.
    * @param B: IsomerBatch<T,K>, the batch of isomers.
    * @param hessians: sycl::buffer<T,1>, the buffer containing the hessians.
    * @param cols: sycl::buffer<K,1>, the buffer containing the column indices of the hessians.
    * @param eigenvalues: sycl::buffer<T,1>, the buffer to store the eigenvalues.
    * @param _nLanczos: size_t, the number of lanczos iterations.
    * @param eigenvectors: sycl::buffer<T,1>, the buffer to store the eigenvectors.
    * @return: void
    * @note: If mode is ENDS_VECTORS is specified, the implementation stores the eigenvalues as [min_0, max_0, min_1, max_1, ...] and
    * the eigenvectors as:
    * [rotx_0, roty_0, rotz_0, transx_0, transy_0, transz_0, min_0, max_0, rotx_1, roty_1, rotz_1, transx_1, transy_1, transz_1, min_1, max_1, ...]
 */
template <EigensolveMode mode, typename T, typename K>
SyclEvent eigensolve(SyclQueue& Q, FullereneBatchView<T,K> B, 
                            Span<T> hessians, 
                            Span<K> cols, 
                            size_t _nLanczos, 
                            Span<T> eigenvalues, 
                            Span<T> eigenvectors,
                            Span<K> indices_,
                            Span<T> off_diagonal,
                            Span<T> qmat,
                            Span<T> lanczos,
                            Span<T> diag,
                            Span<K> ends_idx){
    TEMPLATE_TYPEDEFS(T,K);
    //If mode is ENDS or ENDS_VECTORS, we can make do with fewer lanczos iterations, default is 50.
    size_t nLanczos = (mode == EigensolveMode::ENDS || mode == EigensolveMode::ENDS_VECTORS) ? _nLanczos : B.N_*3 - 6; //-6 for 6 degrees of freedom
    if (nLanczos  > B.N_*3) {
        throw std::runtime_error("Number of lanczos iterations ("+ std::to_string(nLanczos) +") exceeds the number of rank of the hessian matrix");
    }
    size_t capacity = B.capacity(), Natoms = B.N_, batch_size = B.size();
 
    //Buffers required: buffer<T>& D_, buffer<T>& L_, buffer<T>& U_, buffer<T>& Q_, 

    //Lanczos
    //CUDA Lanczos Call:        void* kargs[]{(void*)&B,    (void*)&Vs[dev],                (void*)&Ls[dev],        (void*)&eigenvalues,    (void*)&hessians,   (void*)&cols,           (void*)&Nlanczos};
    //CUDA Lanczos Signature:   (const IsomerBatch<DEV> B,  CuArray<T> V_,                  CuArray<T> U,           CuArray<T> D,           const CuArray<T> H, const CuArray<K> cols,  int m)
    //Should be:                B,                          lanczosBuffers[index],    offDiagBuffers[index],  eigenvalues, hessians, cols, n

    //Us is never so we can just delete it
    //Diagonalization (QR) Call:        void* kargs_qr[]{(void*)&B, (void*)&eigenvalues,    (void*)&Ls[dev],    (void*)&Us[dev],    (void*)&Qs[dev],    (void*)&Nlanczos};
    //Diagonalization (QR) Signature:   (const IsomerBatch<DEV> B,  CuArray<T> D_,          CuArray<T> L_,      CuArray<T> U_,      CuArray<T> Q_,      int n)
    //Should be:                        B,                          eigenvalues,            L_,                 NOTHING,                Q_,                n

    auto V_acc = lanczos;
    auto D_acc = diag;
    auto U_acc = off_diagonal;
    auto H_acc = hessians;
    auto cols_acc = cols;
    auto X_acc = B.d_.X_cubic_;
    if (lanczos.size() < Natoms*3*nLanczos*batch_size) throw std::runtime_error("Lanczos buffer (of size " + std::to_string(lanczos.size()) + ") is too small for the batch size " + std::to_string(batch_size) + " and number of lanczos iterations " + std::to_string(nLanczos) + " and number of atoms " + std::to_string(Natoms));
    Q -> submit([&](sycl::handler& h){
        //accessor V_acc(buffers.lanczosBuffers[index], h, write_only);
        //accessor D_acc(buffers.diagBuffers[index], h, read_write);
        //accessor U_acc(buffers.offDiagBuffers[index], h, read_write);
        //accessor H_acc(hessians, h, read_only);
        //accessor cols_acc(cols, h, read_only);
        //accessor X_acc(B.X, h, read_only);
        local_accessor<T,1> smem(range<1>(Natoms*3), h);
        local_accessor<T,1> betas(range<1>(Natoms*3), h);
        local_accessor<T,1> alphas(range<1>(Natoms*3), h);
        h.parallel_for<Lanczos<mode,T,K>>(sycl::nd_range(sycl::range{Natoms*3*batch_size}, sycl::range{Natoms*3}), [=](nd_item<1> nditem){
            auto cta = nditem.get_group();
            auto tid = nditem.get_local_linear_id();
            auto bid = nditem.get_group_linear_id();
            auto atom_idx = tid/3; //Atom index
            auto N = Natoms*3;
            constexpr int M = 30;
            real_t A[M]; //Hessian matrix, threadIdx.x'th row
            node_t C[M]; //Column indices of the threadIdx.x'th row 3-fold degenerate
            real_t* V = V_acc.data() + bid * nLanczos * N + tid;
            coord3d* X_ptr = X_acc.template as_span<coord3d>().data() + Natoms*bid; 

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
            Z[0] = real_t(tid%3 == 0)/sqrt(T(Natoms)); Z[1] = real_t(tid%3 == 1)/sqrt(T(Natoms)); Z[2] = real_t(tid%3 == 2)/sqrt(T(Natoms)); // Translation eigenvectors
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


            //Modified Gram-Schmidt, Also orthogonalizes against the deflation space
            auto MGS = [&](int index){
                sycl::group_barrier(cta);
                real_t result = V[index*N];
                #pragma unroll
                for (int j = 0; j < 6; j++){
                    auto proj = reduce_over_group(cta, result * Z[j], sycl::plus<real_t>{}) * Z[j];
                    result -= proj; //Remove the component along Z[j] from result
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

            // Generate random starting vector
            oneapi::dpl::uniform_real_distribution<T> distr(0.0, 1.0);            
            oneapi::dpl::minstd_rand engine(42, tid);

            V[0*N] =  distr(engine); //Seed the random number generator with the thread id
            V[0*N] /= sqrt(reduce_over_group(cta, V[0*N] * V[0*N], sycl::plus<real_t>{}));
            V[0*N] = MGS(0);
            for (int i = 0; i < nLanczos; i++){
                if (i > 1){
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
                real_t beta = sqrt(reduce_over_group(cta, v * v, sycl::plus<real_t>{}));
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
    auto L_acc = off_diagonal;
    auto Q_acc = qmat;
    auto Eig_acc = eigenvalues;
    auto Idx_acc = ends_idx;
    SyclEventImpl final_compute = Q -> submit([&](sycl::handler& h){
        local_accessor<T,1> D(range<1>(nLanczos + 1), h); //Diagonal
        local_accessor<T,1> L(range<1>(nLanczos + 1), h); //Lower subdiagonal
        local_accessor<T,1> U(range<1>(nLanczos*2 + 1), h); //Upper subdiagonal
        local_accessor<T,1> V(range<1>(nLanczos*2), h); 
        
        //accessor D_acc(buffers.diagBuffers[index], h, read_only);
        //accessor L_acc(buffers.offDiagBuffers[index], h, read_write);
        //accessor Q_acc(buffers.QmatBuffers[index], h, read_write);
        //accessor Eig_acc(eigenvalues, h, write_only);
        //accessor Idx_acc(buffers.endsIdxBuffers[index], h, write_only);

        h.parallel_for<QR<mode,T,K>>(sycl::nd_range(sycl::range{batch_size*64}, sycl::range{64}), [=](sycl::nd_item<1> nditem){
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
                real_t GR = (k>0?sycl::abs(L[k-1]):0)+sycl::abs(L[k]);
                int not_done = 1;
                while (not_done > 0){
                    i++;
                    T_QTQ(cta, k+1, D, L, U, V, shift);
                    if(mode == EigensolveMode::FULL_SPECTRUM_VECTORS || mode == EigensolveMode::ENDS_VECTORS){
                        apply_all_reflections(cta, V,k,n,Q_acc.data() + n*n*bid);
                    }
                    GR = (k>0?sycl::abs(L[k-1]):0)+(k+1<n?sycl::abs(L[k]):0);

                    if(k>0){
                        std::array<T,4> args = {D[k-1], L[k-1], L[k-1], D[k]};
                        auto [l0, l1] = eigvalsh2x2(args);
                        shift = sycl::abs(l0-d) < sycl::abs(l1-d)? l0:l1;
                    } else {shift = D[k];}
                
                    if(GR <= std::numeric_limits<real_t>::epsilon()*real_t(10.)) not_done--; // Do one (or optionally more) steps after reaching tolerance, to get all off-diagonals below.
                                                    // GPU NB: Se GPU NB ovenfor.
                    if(i>10){
                        //printf("%dth run: Cannot converge eigenvalue %d to tolerance " G " using machine precision %g (d=%g, shift=%g, G=%g)\n" "D[k] = %g, L[k-1] = %g, L[k] = %g\n", nth_time,k,tolerance, std::numeric_limits<real_t>::epsilon(),d,shift,GR, D[k], (k>0)?L[k-1]:0, (k+1<n)?L[k]:0);
                        auto max_error = std::max(std::numeric_limits<real_t>::epsilon()*real_t(10.),GR);
                        break;
                    }
                }
            }
            sycl::group_barrier(cta);
           
            //Copy back to global memory
            if (mode == EigensolveMode::FULL_SPECTRUM || mode == EigensolveMode::FULL_SPECTRUM_VECTORS)
            for (int i = tid; i < n; i += bdim){
                if(i < 6) {Eig_acc[(n+6)*bid + i] = 0;}
                Eig_acc[(n+6)*bid + 6 + i] = D[i];
            } 

            if (mode == EigensolveMode::ENDS || mode == EigensolveMode::ENDS_VECTORS){
                real_t local_max = real_t(0.);
                real_t local_min = std::numeric_limits<real_t>::max();
                int local_min_idx = 0;
                int local_max_idx = 0;
                for (int i = tid; i < n; i += bdim){
                    local_max = isnan(D[i]) ? NAN : std::max(local_max, sycl::abs(D[i]));
                    local_min = isnan(D[i]) ? NAN : std::min(local_min, sycl::abs(D[i]));
                    if (mode == EigensolveMode::ENDS_VECTORS){
                    local_min_idx = isnan(D[i]) ? NAN : (local_min == sycl::abs(D[i]) ? i : local_min_idx);
                    local_max_idx = isnan(D[i]) ? NAN : (local_max == sycl::abs(D[i]) ? i : local_max_idx);
                    }
                }
                real_t max_eig = reduce_over_group(cta, local_max, sycl::maximum<real_t>{});
                real_t min_eig = reduce_over_group(cta, local_min > 1e-1 ? local_min : std::numeric_limits<real_t>::max(), sycl::minimum<real_t>{});
                if (inclusive_scan_over_group(cta, K(min_eig == local_min), sycl::plus<K>{}) == 1) {Idx_acc[2*bid] = local_min_idx; Eig_acc[2*bid] = min_eig;}   //If the thread is the lowest thread-id to find the minimum eigenvalue, store the index
                if (inclusive_scan_over_group(cta, K(max_eig == local_max), sycl::plus<K>{}) == 1) {Idx_acc[2*bid + 1] = local_max_idx; Eig_acc[2*bid + 1] = max_eig;} //If the thread is the lowest thread-id to find the maximum eigenvalue, store the index
            }
        });
    });
    //Eigenvector calculation call: void* kargs_vector[]{(void*)&B, (void*)&Qs[dev], (void*)&Vs[dev], (void*)&Q, (void*)&Nlanczos};
    //Eigenvector calculation signature: (const IsomerBatch<DEV> B, CuArray<T> Q, CuArray<T> V, CuArray<T> E, int m)
    auto E_acc = eigenvectors;
    if (mode == EigensolveMode::FULL_SPECTRUM_VECTORS || mode == EigensolveMode::ENDS_VECTORS){
    final_compute = Q -> submit([&](sycl::handler& h){
        //accessor V_acc(buffers.lanczosBuffers[index], h, read_only);
        //accessor Q_acc(buffers.QmatBuffers[index], h, read_write);
        //accessor E_acc(eigenvectors, h, write_only);
        //accessor X_acc(B.X, h, read_only);
        //accessor Idx_acc(buffers.endsIdxBuffers[index], h, read_only);
        h.parallel_for<Vectors<mode,T,K>>(sycl::nd_range(sycl::range{batch_size*Natoms*3}, sycl::range{Natoms*3}), [=](sycl::nd_item<1> nditem){
            auto tid = nditem.get_local_linear_id();
            auto bid = nditem.get_group_linear_id();
            auto cta = nditem.get_group();
            auto bdim = cta.get_local_linear_range();
            int atom_idx = tid/3; //Atom index (Integer division so the result is rounded down)
            int n = Natoms*3;
            int m = nLanczos;
            int offset = 6;
            real_t* v = V_acc.data() + bid * m * n;
            real_t* e = E_acc.data() + bid * n * n;
            real_t* q = Q_acc.data() + bid * m * m;

            for(int i = 0; i < 3; i++){
                e[i*n + tid] = real_t(tid%3 == i)/sqrt(Natoms); 
            }
            coord3d* X_ptr = X_acc.template as_span<coord3d>().data() + Natoms*bid;
            if (tid%3 == 0) {
                e[3*n + tid] = real_t(0.);
                e[4*n + tid] = -X_ptr[atom_idx][2];
                e[5*n + tid] = -X_ptr[atom_idx][1];
            } else if (tid%3 == 1) {
                e[3*n + tid] = -X_ptr[atom_idx][2];
                e[4*n + tid] = real_t(0.);
                e[5*n + tid] = X_ptr[atom_idx][0];
            } else {
                e[3*n + tid] = X_ptr[atom_idx][1];
                e[4*n + tid] = X_ptr[atom_idx][0];
                e[5*n + tid] = real_t(0.);
            }
                e[3*n + tid] /= sqrt(reduce_over_group(cta, e[3*n + tid]*e[3*n + tid], sycl::plus<real_t>{}));
                e[4*n + tid] /= sqrt(reduce_over_group(cta, e[4*n + tid]*e[4*n + tid], sycl::plus<real_t>{}));
                e[5*n + tid] /= sqrt(reduce_over_group(cta, e[5*n + tid]*e[5*n + tid], sycl::plus<real_t>{}));
            
            if(mode == EigensolveMode::FULL_SPECTRUM_VECTORS)
            for (int k = offset; k < n; k++){
                e = E_acc.data() + bid * n * n + k * n;
                e[tid] = real_t(0.);
                q = Q_acc.data() + bid * m * m + (k-offset) * m;
                for (int i = 0; i < m; i++){
                    e[tid] += v[i*n + tid] * q[i];
                }
            }
            if(mode == EigensolveMode::ENDS_VECTORS){
                int min_idx = Idx_acc[2*bid];
                int max_idx = Idx_acc[2*bid + 1];
                real_t* emin = E_acc.data() + bid * n * 8;
                real_t* emax = E_acc.data() + bid * n * 8 + n;
                real_t* qmin = Q_acc.data() + bid * m * m + min_idx * m;
                real_t* qmax = Q_acc.data() + bid * m * m + max_idx * m;
                emin[tid] = real_t(0.);
                emax[tid] = real_t(0.);
                for (int i = 0; i < m; i++){
                    emin[tid] += v[i*n + tid] * qmin[i];
                    emax[tid] += v[i*n + tid] * qmax[i];
                }

            }
            
        });
    });
    }
    return final_compute;
}

template <EigensolveMode mode, typename T, typename K>
template <typename... Data>
SyclEvent EigenFunctor<mode, T, K>::compute(SyclQueue& Q, FullereneBatchView<T,K> B, 
                            Span<T> hessians,
                            Span<K> cols,
                            size_t _nLanczos,
                            Span<T> eigenvalues,
                            Span<T> eigenvectors,
                            Data&&... data){
    return eigensolve<mode>(Q, B, hessians, cols, _nLanczos, eigenvalues, eigenvectors, std::forward<Data>(data)...);
}

template <EigensolveMode mode, typename T, typename K>
template <typename... Data>
SyclEvent EigenFunctor<mode, T, K>::compute(SyclQueue& Q, Fullerene<T,K> B, 
                            Span<T> hessians,
                            Span<K> cols,
                            size_t _nLanczos,
                            Span<T> eigenvalues,
                            Span<T> eigenvectors,
                            Data&&... data){
    throw std::logic_error("Not implemented");
}
    
template struct EigenFunctor<EigensolveMode::FULL_SPECTRUM, float, uint16_t>;
template struct EigenFunctor<EigensolveMode::ENDS, float, uint16_t>;
template struct EigenFunctor<EigensolveMode::ENDS_VECTORS, float, uint16_t>;
template struct EigenFunctor<EigensolveMode::FULL_SPECTRUM_VECTORS, float, uint16_t>;

template SyclEvent EigenFunctor<EigensolveMode::FULL_SPECTRUM, float, uint16_t>::compute(SyclQueue& Q, FullereneBatchView<float,uint16_t> B, 
                            Span<float> hessians,
                            Span<uint16_t> cols,
                            size_t _nLanczos,
                            Span<float> eigenvalues,
                            Span<float> eigenvectors,
                            Span<uint16_t>& indices,
                            Span<float>& off_diagonal,
                            Span<float>& diag,
                            Span<float>& qmat,
                            Span<float>& lanczos,
                            Span<uint16_t>& ends_idx);

template SyclEvent EigenFunctor<EigensolveMode::FULL_SPECTRUM, float, uint16_t>::compute(SyclQueue& Q, Fullerene<float,uint16_t> B, 
                            Span<float> hessians,
                            Span<uint16_t> cols,
                            size_t _nLanczos,
                            Span<float> eigenvalues,
                            Span<float> eigenvectors,
                            Span<uint16_t>& indices,
                            Span<float>& off_diagonal,
                            Span<float>& diag,
                            Span<float>& qmat,
                            Span<float>& lanczos,
                            Span<uint16_t>& ends_idx);

template SyclEvent EigenFunctor<EigensolveMode::ENDS, float, uint16_t>::compute(SyclQueue& Q, FullereneBatchView<float,uint16_t> B,
                            Span<float> hessians,
                            Span<uint16_t> cols,
                            size_t _nLanczos,
                            Span<float> eigenvalues,
                            Span<float> eigenvectors,
                            Span<uint16_t>& indices,
                            Span<float>& off_diagonal,
                            Span<float>& diag,
                            Span<float>& qmat,
                            Span<float>& lanczos,
                            Span<uint16_t>& ends_idx);

template SyclEvent EigenFunctor<EigensolveMode::ENDS, float, uint16_t>::compute(SyclQueue& Q, Fullerene<float,uint16_t> B,
                            Span<float> hessians,
                            Span<uint16_t> cols,
                            size_t _nLanczos,
                            Span<float> eigenvalues,
                            Span<float> eigenvectors,
                            Span<uint16_t>& indices,
                            Span<float>& off_diagonal,
                            Span<float>& diag,
                            Span<float>& qmat,
                            Span<float>& lanczos,
                            Span<uint16_t>& ends_idx);

template SyclEvent EigenFunctor<EigensolveMode::ENDS_VECTORS, float, uint16_t>::compute(SyclQueue& Q, FullereneBatchView<float,uint16_t> B,
                            Span<float> hessians,
                            Span<uint16_t> cols,
                            size_t _nLanczos,
                            Span<float> eigenvalues,
                            Span<float> eigenvectors,
                            Span<uint16_t>& indices,
                            Span<float>& off_diagonal,
                            Span<float>& diag,
                            Span<float>& qmat,
                            Span<float>& lanczos,
                            Span<uint16_t>& ends_idx);

template SyclEvent EigenFunctor<EigensolveMode::ENDS_VECTORS, float, uint16_t>::compute(SyclQueue& Q, Fullerene<float,uint16_t> B,
                            Span<float> hessians,
                            Span<uint16_t> cols,
                            size_t _nLanczos,
                            Span<float> eigenvalues,
                            Span<float> eigenvectors,
                            Span<uint16_t>& indices,
                            Span<float>& off_diagonal,
                            Span<float>& diag,
                            Span<float>& qmat,
                            Span<float>& lanczos,
                            Span<uint16_t>& ends_idx);

template SyclEvent EigenFunctor<EigensolveMode::FULL_SPECTRUM_VECTORS, float, uint16_t>::compute(SyclQueue& Q, FullereneBatchView<float,uint16_t> B,
                            Span<float> hessians,
                            Span<uint16_t> cols,
                            size_t _nLanczos,
                            Span<float> eigenvalues,
                            Span<float> eigenvectors,
                            Span<uint16_t>& indices,
                            Span<float>& off_diagonal,
                            Span<float>& diag,
                            Span<float>& qmat,
                            Span<float>& lanczos,
                            Span<uint16_t>& ends_idx);

template SyclEvent EigenFunctor<EigensolveMode::FULL_SPECTRUM_VECTORS, float, uint16_t>::compute(SyclQueue& Q, Fullerene<float,uint16_t> B,
                            Span<float> hessians,
                            Span<uint16_t> cols,
                            size_t _nLanczos,
                            Span<float> eigenvalues,
                            Span<float> eigenvectors,
                            Span<uint16_t>& indices,
                            Span<float>& off_diagonal,
                            Span<float>& diag,
                            Span<float>& qmat,
                            Span<float>& lanczos,
                            Span<uint16_t>& ends_idx);

/* template void eigensolve<EigensolveMode::FULL_SPECTRUM, float, uint16_t>(sycl::queue& ctx, const IsomerBatch<float,uint16_t> B, sycl::buffer<float,1>& hessians, sycl::buffer<uint16_t,1>& cols, sycl::buffer<float,1>& eigenvalues, const LaunchPolicy policy, size_t _nLanczos, sycl::buffer<float,1>& eigenvectors);
template void eigensolve<EigensolveMode::ENDS, float, uint16_t>(sycl::queue& ctx, const IsomerBatch<float,uint16_t> B, sycl::buffer<float,1>& hessians, sycl::buffer<uint16_t,1>& cols, sycl::buffer<float,1>& eigenvalues, const LaunchPolicy policy, size_t _nLanczos, sycl::buffer<float,1>& eigenvectors);
template void eigensolve<EigensolveMode::ENDS_VECTORS, float, uint16_t>(sycl::queue& ctx, const IsomerBatch<float,uint16_t> B, sycl::buffer<float,1>& hessians, sycl::buffer<uint16_t,1>& cols, sycl::buffer<float,1>& eigenvalues, const LaunchPolicy policy, size_t _nLanczos, sycl::buffer<float,1>& eigenvectors);
template void eigensolve<EigensolveMode::FULL_SPECTRUM_VECTORS, float, uint16_t>(sycl::queue& ctx, const IsomerBatch<float,uint16_t> B, sycl::buffer<float,1>& hessians, sycl::buffer<uint16_t,1>& cols, sycl::buffer<float,1>& eigenvalues, const LaunchPolicy policy, size_t _nLanczos, sycl::buffer<float,1>& eigenvectors);

template void eigensolve<EigensolveMode::FULL_SPECTRUM, double, uint16_t>(sycl::queue& ctx, const IsomerBatch<double,uint16_t> B, sycl::buffer<double,1>& hessians, sycl::buffer<uint16_t,1>& cols, sycl::buffer<double,1>& eigenvalues, const LaunchPolicy policy, size_t _nLanczos, sycl::buffer<double,1>& eigenvectors);
template void eigensolve<EigensolveMode::ENDS, double, uint16_t>(sycl::queue& ctx, const IsomerBatch<double,uint16_t> B, sycl::buffer<double,1>& hessians, sycl::buffer<uint16_t,1>& cols, sycl::buffer<double,1>& eigenvalues, const LaunchPolicy policy, size_t _nLanczos, sycl::buffer<double,1>& eigenvectors);
template void eigensolve<EigensolveMode::ENDS_VECTORS, double, uint16_t>(sycl::queue& ctx, const IsomerBatch<double,uint16_t> B, sycl::buffer<double,1>& hessians, sycl::buffer<uint16_t,1>& cols, sycl::buffer<double,1>& eigenvalues, const LaunchPolicy policy, size_t _nLanczos, sycl::buffer<double,1>& eigenvectors);
template void eigensolve<EigensolveMode::FULL_SPECTRUM_VECTORS, double, uint16_t>(sycl::queue& ctx, const IsomerBatch<double,uint16_t> B, sycl::buffer<double,1>& hessians, sycl::buffer<uint16_t,1>& cols, sycl::buffer<double,1>& eigenvalues, const LaunchPolicy policy, size_t _nLanczos, sycl::buffer<double,1>& eigenvectors); */
//template void eigensolve<EigensolveMode::FULL_SPECTRUM, double, uint16_t>(sycl::queue& ctx, const IsomerBatch<double,uint16_t> B, sycl::buffer<double,1>& hessians, sycl::buffer<uint16_t,1>& cols, sycl::buffer<double,1>& eigenvalues, const LaunchPolicy policy, size_t _nLanczos, sycl::buffer<double,1>& eigenvectors);