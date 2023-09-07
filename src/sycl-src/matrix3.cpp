#ifndef DEVICE_MAT3
#define DEVICE_MAT3
#include <array>
struct identity3
{
    identity3(){}
};



template <typename T>
struct mat3
{   
    FLOAT_TYPEDEFS(T);
    T A[9];

    mat3()
    {
        #pragma unroll
        for (int i = 0; i < 9; i++)
            A[i] = T(0.f);
    }

  mat3(const std::array<coord3d,3> &V){
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	A[i*3+j] = V[i][j];
  }
  
    mat3(T a, T b, T c, T d, T e, T f, T g, T h, T i)
    {
        A[0] = a;
        A[1] = b;
        A[2] = c;
        A[3] = d;
        A[4] = e;
        A[5] = f;
        A[6] = g;
        A[7] = h;
        A[8] = i;
    }

    T* operator[](int i) const
    {
        return const_cast<T*>(&A[i * 3]);
    }

    //Assignment operators
    mat3(const mat3& a)
    {
        #pragma unroll
        for (int i = 0; i < 9; i++)
            A[i] = a.A[i];
    }
    static mat3 identity(){ return mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);}
};



template <typename T> mat3<T> operator*(const mat3<T>& a, const mat3<T>& b)
{
    mat3<T> result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++){
            #pragma unroll
            for (int k = 0; k < 3; k++)
                result[i][j] += a[i][k] * b[k][j];
        }
    return result;
}


template <typename T> mat3<T> operator*(const mat3<T>& a, T b)
{
    mat3<T> result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[i][j] * b;
    return result;
}

template <typename T> mat3<T> operator/(const mat3<T>& a, T b){ return a * (T(1.f) / b); }

template <typename T> mat3<T> operator*(const mat3<T>& a, const identity3& b){ return a; }
template <typename T> mat3<T> operator*(const identity3& a, const mat3<T>& b){ return b; }
template <typename T> mat3<T> operator+(const mat3<T>& a, const identity3& b){ 
    mat3<T> result = a;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        result[i][i] += T(1.f);
    return result;
}
template <typename T> mat3<T> operator+(const identity3& a, const mat3<T>& b){ return b + a; }
template <typename T> mat3<T> operator-(const mat3<T>& a, const identity3& b){ 
    mat3<T> result = a;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        result[i][i] -= T(1.f);
    return result;
}
template <typename T> mat3<T> operator-(const identity3& a, const mat3<T>& b){ 
    mat3<T> result = b;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++){
            if (i == j) result[i][j] = T(1.f) - result[i][j];
            else result[i][j] = -result[i][j];
        }
    return result;
}

template <typename T> mat3<T> operator-(const mat3<T>& a){ 
    mat3<T> result;
    #pragma unroll
    for (int i = 0; i < 9; i++)
        result.A[i] = -a.A[i];
    return result;
}
/* __device__ mat3<T> opeartor-(const mat3<T>& a){ 
    return mat3<T>(-a[0][0], -a[0][1], -a[0][2],
                -a[1][0], -a[1][1], -a[1][2],
                -a[2][0], -a[2][1], -a[2][2]);
} */

template <typename T> mat3<T> operator*(T a, const mat3<T>& b){ return b * a; }

template <typename T> mat3<T> operator+(const mat3<T>& a, const mat3<T>& b)
{
    mat3<T> result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[i][j] + b[i][j];
    return result;
}

template <typename T> void operator+=(mat3<T>& a, const mat3<T>& b)
{
    #pragma unroll
    for (int i = 0; i < 9; i++)
        a.A[i] += b.A[i];
}

template <typename T> mat3<T> operator-(const mat3<T>& a, const mat3<T>& b)
{
    mat3<T> result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[i][j] - b[i][j];
    return result;
}

//Computes the L2 norm of a matrix
template <typename T> T norm(const mat3<T>& a){
    T result = T(0.f);
    #pragma unroll
    for (int i = 0; i < 9; i++)
        result += a.A[i] * a.A[i];
    return SQRT(result);
} 

template <typename T> mat3<T> transpose(const mat3<T>& a)
{
    mat3<T> result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[j][i];
    return result;
}

template <typename T> std::array<T,3> dot(const mat3<T>& a, const std::array<T,3>& b)
{
    std::array<T,3> result;
    #pragma unroll
    for (int i = 0; i < 3; i++){
        ((T*)&result)[i] = T(0.f);
        #pragma unroll
        for (int j = 0; j < 3; j++)
            ((T*)&result)[i] += a[i][j] * ((T*)&b)[j];
    }
    return result;
}

template <typename T> std::array<T,3> dot(const std::array<T,3>& a, const mat3<T>& b)
{
    std::array<T,3> result;
    #pragma unroll
    for (int i = 0; i < 3; i++){
        ((T*)&result)[i] = T(0.f);
        #pragma unroll
        for (int j = 0; j < 3; j++)
            ((T*)&result)[i] += ((T*)&a)[j] * b[j][i];
    }
    return result;
}

template <typename T> mat3<T> tensor_product(const std::array<T,3>& a, const std::array<T,3>& b)
{
    mat3<T> result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[i] * b[j];
    return result;
}

//Cross product between a vector and a dyadic tensor (c X A) https://stemandmusic.in/maths/mvt-algebra/dvCP.php
template <typename T> mat3<T> cross(const std::array<T,3>& c, const mat3<T>& A)
{
    return mat3<T>(A[2][0]*c[1] - A[1][0]*c[2], A[2][1]*c[1] - A[1][1]*c[2], A[2][2]*c[1] - A[1][2]*c[2],
                A[0][0]*c[2] - A[2][0]*c[0], A[0][1]*c[2] - A[2][1]*c[0], A[0][2]*c[2] - A[2][2]*c[0],
                A[1][0]*c[0] - A[0][0]*c[1], A[1][1]*c[0] - A[0][1]*c[1], A[1][2]*c[0] - A[0][2]*c[1]);
}

//Cross product between a dyadic tensor and a vector (A X c) https://stemandmusic.in/maths/mvt-algebra/dvCP.php
template <typename T> mat3<T> cross(const mat3<T>& A, const std::array<T,3>& c)
{
    return mat3<T>(A[0][1]*c[2] - A[0][2]*c[1], A[0][2]*c[0] - A[0][0]*c[2], A[0][1]*c[1] - A[0][2]*c[0],
                A[1][1]*c[2] - A[1][2]*c[1], A[1][2]*c[0] - A[1][0]*c[2], A[1][1]*c[1] - A[1][2]*c[0],
                A[2][1]*c[2] - A[2][2]*c[1], A[2][2]*c[0] - A[2][0]*c[2], A[2][1]*c[1] - A[2][2]*c[0]);
}

/* template <int elements_per_thread>
struct sym_tridiag_t{
    T alpha[elements_per_thread];
    T beta[elements_per_thread];

    template <typename T> sym_tridiag_t()
    {
        #pragma unroll
        for (int i = 0; i < elements_per_thread; i++){
            alpha[i] = T(0.f);
            beta[i] = T(0.f);
        }
    }
};

template <int elements_per_thread>
template <typename T> void mat_mul(const sym_tridiag_t<elements_per_thread>& a, const std::array<T,3> b){
    std::array<T,3> result;
    #pragma unroll
    for (int i = 0; i < elements_per_thread; i++){
        result[i] = a.alpha[i] * b[i];
        if(i > 0) result[i] += a.beta[i-1] * b[i-1];
        if(i < elements_per_thread-1) result[i] += a.beta[i] * b[i+1];
    }
} */

template <typename T> void print(const mat3<T>& a, int thread_id = 0)
{   /* 
    if(threadIdx.x != thread_id)
        return;
    printf("{");
    #pragma unroll
    for (int i = 0; i < 3; i++){
        #pragma unroll
        for (int j = 0; j < 3; j++)
            printf("%f ", a[i][j]);
        
        if(i != 2) printf("\n");
        else printf("}\n");
    } */
}

template <typename T, typename K>
struct hessian_t{
    TEMPLATE_TYPEDEFS(T,K);
    //Sparse hessian_t matrix stored in such a way that each thread owns 1 row of the matrix.
    mat3<T> A[10]; // Each thread has 10*3*3 entries in the hessian_t matrix. Itself and 9 neighbours
    node_t indices[10]; 
    // Each atom has 3 neighbours and 6 outer neighbours
    // We store them in the following order:
    // a, b, c, d, b_m, c_m, d_m, b_p, c_p, d_p

    hessian_t(const sycl::group<1>& cta, const NodeNeighbours<K>& G){
        indices[0] = cta.get_local_linear_id();
        indices[1] = d_get(G.cubic_neighbours,0);
        indices[2] = d_get(G.cubic_neighbours,1);
        indices[3] = d_get(G.cubic_neighbours,2);
        indices[4] = d_get(G.prev_on_face,0);
        indices[5] = d_get(G.prev_on_face,1);
        indices[6] = d_get(G.prev_on_face,2);
        indices[7] = d_get(G.next_on_face,0);
        indices[8] = d_get(G.next_on_face,1);
        indices[9] = d_get(G.next_on_face,2);

        #pragma unroll
        for (int i = 0; i < 10; i++)
            A[i] = mat3<T>();
    }

    hessian_t(const hessian_t& other){
        #pragma unroll
        for (int i = 0; i < 10; i++){
            A[i] = other.A[i];
            indices[i] = other.indices[i];
        }
    }

    const mat3<T>& operator[](const int i) const { return A[i]; } 
    mat3<T>& operator[](const int i) { return A[i]; }
/*     T power_iteration(std::array<T,3>* smem);
    template <int eigenvalues_per_thread>
    sym_tridiag_t<eigenvalues_per_thread> lanczos_iteration(T* smem);
 */
};
/* 
device_coord3d mat_vect_mult(const hessian_t& A, const device_coord3d& x,device_coord3d* smem)
{
    BLOCK_SYNC;
    smem[threadIdx.x] = x;
    BLOCK_SYNC;
    device_coord3d result = {0.f,0.f,0.f};
    #pragma unroll
    for (int i = 0; i < 10; i++)
        result += dot(A[i], smem[A.indices[i]]);
    return result;
}

T hessian_t::power_iteration(device_coord3d* smem)
{
    clear_cache(reinterpret_cast<T*>(smem), 3*blockDim.x);
    device_coord3d b_k;
    device_coord3d b_kp1;
    curandState state;
    curand_init(clock() + threadIdx.x, 0, 0, &state);
    b_k[0] = curand_uniform(&state);
    b_k[1] = curand_uniform(&state);
    b_k[2] = curand_uniform(&state);
    b_k = {T(threadIdx.x*3), T(threadIdx.x*3 + 1), T(threadIdx.x*3 + 2)};
    T norm;
    for (size_t i = 0; i < 500; i++){ 
        b_kp1 =  mat_vect_mult(*this, b_k, smem + blockDim.x);
        norm = SQRT(reduction(reinterpret_cast<T*>(smem), dot(b_kp1, b_kp1)));
        b_k = b_kp1 /norm;
    }
    //Return rayleigh quotient
    T lambda = reduction(reinterpret_cast<T*>(smem), dot(b_k, mat_vect_mult(*this, b_k, smem + blockDim.x))) / reduction(reinterpret_cast<T*>(smem), dot(b_k, b_k));
    auto bk_new = mat_vect_mult(*this, b_k, smem);
    auto numerator = reduction(reinterpret_cast<T*>(smem), dot(b_k, mat_vect_mult(*this, b_k, smem + blockDim.x)));
    clear_cache(reinterpret_cast<T*>(smem), 3*blockDim.x);
    auto denominator = reduction(reinterpret_cast<T*>(smem), dot(b_k, b_k));
    return numerator/denominator;
}
T dot(const device_coord3d& v, const device_coord3d& u, T* smem){
    return reduction(smem, dot(v,v));
}

device_coord3d orthogonalize(const device_coord3d& v, const device_coord3d& u, T* smem){
    return v - dot(v,u,smem) / dot(u,u,smem) * u;   
}

//Lanczos algorithm: http://www.cs.cmu.edu/afs/cs/academic/class/15859n-f16/Handouts/TrefethenBau/LanczosIteration-36.pdf
//This method produces a symmetric tridiagonal matrix T_m of dimension m x m
//The matrix is represented by 2 vectors: alpha and beta
//Each thread stores and returns up to 3 values of alpha and beta
template <int eigenvalues_per_thread>
sym_tridiag_t<eigenvalues_per_thread> hessian_t::lanczos_iteration(T* smem){
  DEVICE_TYPEDEFS;
  
    int N = blockDim.x;
    coord3d* mat_mem = reinterpret_cast<coord3d*>(smem + N*2);
    int m = eigenvalues_per_thread * N;
    assert(m <= N*3); //A is a 3Nx3N matrix, so m must be <= 3N
    real_t beta  = real_t(0.);
    real_t alpha = real_t(0.);
    curandState state;
    curand_init(42 + threadIdx.x, 0, 0, &state); //( ͡° ͜ʖ ͡°)
    coord3d b = (coord3d){curand_uniform(&state), curand_uniform(&state), curand_uniform(&state)};
    coord3d q_km1 = (coord3d){0.,0.,0.};
    coord3d q_k = b / SQRT(reduction(smem, dot(b, b)));
    sym_tridiag_t<eigenvalues_per_thread> result;
    for (int i = 0; i < m ; i++){
        coord3d v = mat_vect_mult(*this, q_k, mat_mem);
        alpha = reduction(smem, dot(q_k, v));
        v -= beta * q_km1 + alpha * q_k;
        beta = SQRT(reduction(smem, dot(v, v)));
        q_k = v / beta;
        if (i  >= threadIdx.x * 3 && i < (threadIdx.x + 1)*3){
            result.alpha[i % 3] = alpha;
            result.beta[i % 3] = beta;
        }
    }
    sequential_print((coord3d){result.alpha[0], result.alpha[1], result.alpha[2]}, 0);
    return result; 
}

 */


#endif
