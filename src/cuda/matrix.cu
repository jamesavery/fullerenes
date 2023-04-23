#ifndef DEVICE_MAT3
#define DEVICE_MAT3

struct sym3
{
    float A[6];

    INLINE sym3()
    {
        #pragma unroll
        for (int i = 0; i < 6; i++)
            A[i] = device_real_t(0.f);
    }

    INLINE sym3(float a, float b, float c, float d, float e, float f)
    {
        A[0] = a;
        A[1] = b;
        A[2] = c;
        A[3] = d;
        A[4] = e;
        A[5] = f;
    }

    INLINE float operator[](int i) const
    {
        return A[i];
    }
};

struct identity3
{
    INLINE identity3(){}
};


struct mat3
{
    float A[9];

    INLINE mat3()
    {
        #pragma unroll
        for (int i = 0; i < 9; i++)
            A[i] = device_real_t(0.f);
    }

    INLINE mat3(float a, float b, float c, float d, float e, float f, float g, float h, float i)
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

    INLINE float* operator[](int i) const
    {
        return const_cast<float*>(&A[i * 3]);
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


INLINE mat3 operator*(const mat3& a, float b)
{
    mat3 result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = a[i][j] * b;
    return result;
}

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
        result[i][i] = device_real_t(1.f) - result[i][i];
    return result;
}


INLINE mat3 operator*(float a, const mat3& b){ return b * a; }

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
        ((float*)&result)[i] = device_real_t(0.f);
        #pragma unroll
        for (int j = 0; j < 3; j++)
            ((float*)&result)[i] += a[i][j] * ((float*)&b)[j];
    }
    return result;
}

INLINE device_coord3d dot(const device_coord3d& a, const mat3& b)
{
    device_coord3d result;
    #pragma unroll
    for (int i = 0; i < 3; i++){
        ((float*)&result)[i] = device_real_t(0.f);
        #pragma unroll
        for (int j = 0; j < 3; j++)
            ((float*)&result)[i] += ((float*)&a)[j] * b[j][i];
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
            result[i][j] = ((float*)&a)[i] * ((float*)&b)[j];
    return result;
}

INLINE mat3 identity()
{
    mat3 result;
    result[0][0] = device_real_t(1.f);
    result[1][1] = device_real_t(1.f);
    result[2][2] = device_real_t(1.f);
    return result;
}

//TODO: check if this is correct
INLINE mat3 cross(const device_coord3d& a, const mat3& b)
{
    mat3 result;
    #pragma unroll
    for (int i = 0; i < 3; i++)
        #pragma unroll
        for (int j = 0; j < 3; j++)
            result[i][j] = ((float*)&a)[i] * b[j][i];
    return result;
}





#endif