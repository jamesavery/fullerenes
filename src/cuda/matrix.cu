#ifndef DEVICE_MAT3
#define DEVICE_MAT3

struct identity3
{
    INLINE identity3(){}
};



struct mat3
{
    device_real_t A[9];

    INLINE mat3()
    {
        #pragma unroll
        for (int i = 0; i < 9; i++)
            A[i] = device_real_t(0.f);
    }

    INLINE mat3(device_real_t a, device_real_t b, device_real_t c, device_real_t d, device_real_t e, device_real_t f, device_real_t g, device_real_t h, device_real_t i)
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

    INLINE device_real_t* operator[](int i) const
    {
        return const_cast<device_real_t*>(&A[i * 3]);
    }

    //Assignment operators
    INLINE mat3(const mat3& a)
    {
        #pragma unroll
        for (int i = 0; i < 9; i++)
            A[i] = a.A[i];
    }
};


INLINE mat3 operator*(const mat3& a, const mat3& b)
{
    mat3 result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++){
            result[i][j] = device_real_t(0.f);
            #pragma unroll
            for (int k = 0; k < 3; k++)
                result[i][j] += a[i][k] * b[k][j];
        }
    return result;
}


INLINE mat3 operator*(const mat3& a, device_real_t b)
{
    mat3 result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[i][j] * b;
    return result;
}

INLINE mat3 operator/(const mat3& a, device_real_t b){ return a * (device_real_t(1.f) / b); }

INLINE mat3 operator*(const mat3& a, const identity3& b){ return a; }
INLINE mat3 operator*(const identity3& a, const mat3& b){ return b; }
INLINE mat3 operator+(const mat3& a, const identity3& b){ 
    mat3 result = a;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        result[i][i] += device_real_t(1.f);
    return result;
}
INLINE mat3 operator+(const identity3& a, const mat3& b){ return b + a; }
INLINE mat3 operator-(const mat3& a, const identity3& b){ 
    mat3 result = a;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        result[i][i] -= device_real_t(1.f);
    return result;
}
INLINE mat3 operator-(const identity3& a, const mat3& b){ 
    mat3 result = b;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++){
            if (i == j) result[i][j] = device_real_t(1.f) - result[i][j];
            else result[i][j] = -result[i][j];
        }
    return result;
}

INLINE mat3 operator-(const mat3& a){ 
    mat3 result;
    #pragma unroll
    for (int i = 0; i < 9; i++)
        result.A[i] = -a.A[i];
    return result;
}
/* __device__ mat3 opeartor-(const mat3& a){ 
    return mat3(-a[0][0], -a[0][1], -a[0][2],
                -a[1][0], -a[1][1], -a[1][2],
                -a[2][0], -a[2][1], -a[2][2]);
} */

INLINE mat3 operator*(device_real_t a, const mat3& b){ return b * a; }

INLINE mat3 operator+(const mat3& a, const mat3& b)
{
    mat3 result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[i][j] + b[i][j];
    return result;
}

INLINE void operator+=(mat3& a, const mat3& b)
{
    #pragma unroll
    for (int i = 0; i < 9; i++)
        a.A[i] += b.A[i];
}

INLINE mat3 operator-(const mat3& a, const mat3& b)
{
    mat3 result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[i][j] - b[i][j];
    return result;
}

INLINE mat3 transpose(const mat3& a)
{
    mat3 result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[j][i];
    return result;
}

INLINE device_coord3d dot(const mat3& a, const device_coord3d& b)
{
    device_coord3d result;
    #pragma unroll
    for (int i = 0; i < 3; i++){
        ((device_real_t*)&result)[i] = device_real_t(0.f);
        #pragma unroll
        for (int j = 0; j < 3; j++)
            ((device_real_t*)&result)[i] += a[i][j] * ((device_real_t*)&b)[j];
    }
    return result;
}

INLINE device_coord3d dot(const device_coord3d& a, const mat3& b)
{
    device_coord3d result;
    #pragma unroll
    for (int i = 0; i < 3; i++){
        ((device_real_t*)&result)[i] = device_real_t(0.f);
        #pragma unroll
        for (int j = 0; j < 3; j++)
            ((device_real_t*)&result)[i] += ((device_real_t*)&a)[j] * b[j][i];
    }
    return result;
}

INLINE mat3 tensor_product(const device_coord3d& a, const device_coord3d& b)
{
    mat3 result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[i] * b[j];
    return result;
}

//Cross product between a vector and a dyadic tensor (c X A) https://stemandmusic.in/maths/mvt-algebra/dvCP.php
INLINE mat3 cross(const device_coord3d& c, const mat3& A)
{
    return mat3(A[2][0]*c[1] - A[1][0]*c[2], A[2][1]*c[1] - A[1][1]*c[2], A[2][2]*c[1] - A[1][2]*c[2],
                A[0][0]*c[2] - A[2][0]*c[0], A[0][1]*c[2] - A[2][1]*c[0], A[0][2]*c[2] - A[2][2]*c[0],
                A[1][0]*c[0] - A[0][0]*c[1], A[1][1]*c[0] - A[0][1]*c[1], A[1][2]*c[0] - A[0][2]*c[1]);
}

//Cross product between a dyadic tensor and a vector (A X c) https://stemandmusic.in/maths/mvt-algebra/dvCP.php
INLINE mat3 cross(const mat3& A, const device_coord3d& c)
{
    return mat3(A[0][1]*c[2] - A[0][2]*c[1], A[0][2]*c[0] - A[0][0]*c[2], A[0][1]*c[1] - A[0][2]*c[0],
                A[1][1]*c[2] - A[1][2]*c[1], A[1][2]*c[0] - A[1][0]*c[2], A[1][1]*c[1] - A[1][2]*c[0],
                A[2][1]*c[2] - A[2][2]*c[1], A[2][2]*c[0] - A[2][0]*c[2], A[2][1]*c[1] - A[2][2]*c[0]);
}

INLINE void print(const mat3& a, int thread_id = 0)
{   
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
    }
}


struct hessian_t{
    //Sparse hessian_t matrix stored in such a way that each thread owns 1 row of the matrix.
    mat3 A[10]; // Each thread has 10*3*3 entries in the hessian_t matrix. Itself and 9 neighbours
    device_node_t indices[9]; 
    // Each atom has 3 neighbours and 6 outer neighbours
    // We store them in the following order:
    // a, b, c, d, b_m, c_m, d_m, b_p, c_p, d_p

    INLINE hessian_t(const NodeNeighbours& G){
        indices[0] = threadIdx.x;
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
        for (int i = 0; i < 9; i++)
            A[i] = mat3();
    }

    INLINE hessian_t(const hessian_t& other){
        #pragma unroll
        for (int i = 0; i < 10; i++){
            A[i] = other.A[i];
        }
        #pragma unroll
        for (int i = 0; i < 10; i++){
            indices[i] = other.indices[i];
        }
    }

    INLINE const mat3& operator[](const int i) const { return A[i]; } 
    INLINE mat3& operator[](const int i) { return A[i]; }
};

INLINE device_coord3d mat_vect_mult(const hessian_t& A, const device_coord3d& x,device_coord3d* smem)
{
    BLOCK_SYNC;
    smem[threadIdx.x] = x;
    device_coord3d result;
    #pragma unroll
    for (int i = 0; i < 10; i++)
        result += dot(A[i], smem[A.indices[i]]);
    return result;
}

#endif