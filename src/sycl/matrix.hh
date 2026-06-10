#pragma once
#include <array>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <climits>
#include <assert.h>
#include <string.h>

using std::array;
using std::string;

//typedef double real_t;
#define PRINTF_G "%g"
//typedef long double real_t;
//#define PRINTF_G "%Lg"
typedef std::complex<real_t> complex_t;
//constexpr real_t machine_precision = pow(std::numeric_limits<real_t>::radix,-std::numeric_limits<real_t>::digits);




#define INLINE inline __attribute__((always_inline))
INLINE real_t    conj(real_t x)   { return x; }
INLINE complex_t conj(complex_t x){ return std::conj(x); }
INLINE real_t    Re(complex_t x)  { return x.real(); }
INLINE real_t    Re(real_t    x)  { return x; }
INLINE real_t    Im(complex_t x)  { return x.imag(); }
INLINE real_t    Im(real_t    x)  { return 0; }

// A simple matrix class with no dynamic allocation, suitable for use with fast stack memory
// Because everything operates on pre-allocated memory, many operations require a pointer to a
// result_data pointer.
template <typename scalar> struct matrix {
  scalar *data;
  array<int,2> shape, stride;
  size_t offset;

  matrix(scalar *data=0, array<int,2> shape={0,0}, array<int,2> stride_ = {INT_MIN,INT_MIN}, size_t offset=0) : 
  data(data), shape(shape), offset(offset)
  {
    if(stride_[0]==INT_MIN){ stride[0] = shape[1]; stride[1] = 1; }
    else stride = stride_;
  }

  INLINE matrix transpose() const { return matrix(data,{shape[1],shape[0]}, {stride[1],stride[0]},offset); }
  INLINE matrix T()         const { return transpose(); }
  
  matrix copy(scalar *result_data) const {
    auto [m,n] = shape;
    matrix result(result_data, shape);
    const auto &A(*this);    
    
    for(int i=0;i<m;i++)
      for(int j=0;j<n; j++)
	result(i,j) = A(i,j);

    return result;
  }

  matrix conj(scalar *result_data) const {
    auto [m,n] = shape;
    matrix result(result_data, shape);
    const auto &A(*this);    
    
    for(int i=0;i<m;i++)
      for(int j=0;j<n; j++)
	result(i,j) = ::conj(A(i,j));

    return result;    
  }

  /* Matrix element access pattern: A(i,j) is i,j-element of A, which may be a sub-view of a larger matrix.
   */
  INLINE scalar& operator()(size_t i, size_t j)
  {
    return data[offset + i*stride[0] + j*stride[1]];
  }

  INLINE scalar operator()(size_t i, size_t j) const 
  {
    return data[offset + i*stride[0] + j*stride[1]];
  }

  /* Flat base-array memory access pattern. A[i] is the ith element in the flat base array memory. 
     This won't work on a view: as it is the ith element in flat memory, it does not use the stride or offset.
  */
  INLINE scalar& operator[](ssize_t i)       { return data[i];  }
  INLINE scalar  operator[](ssize_t i) const { return data[i];  }  
  
  // Calculates the matrix max-norm / infinity-norm: max(sum(abs(A),axis=1))
  real_t max_norm() const { 
    const auto &A = *this;
    auto  [m,n]   = A.shape;
    
    real_t mx = 0;
    for(int i=0;i<m;i++){ 	
      real_t row_norm = 0;
      for(int j=0;j<n;j++) row_norm += std::abs(A(i,j));
      mx = std::max(mx,row_norm);
    }
    return mx;
  }

  // Shift diagonal by a scalar amount
  matrix &operator+=(const scalar shift)
  {
    auto &A(*this);
    int n = std::min(shape[0],shape[1]);
    
    for(int i=0;i<n;i++) A(i,i) += shift;

    return A;
  }

  matrix &operator+=(const matrix B)
  {
    auto &A(*this);
    auto [m,n] = shape;
    
    for(int i=0;i<m;i++)
      for(int j=0;j<n;j++)
	A(i,j) += B(i,j);

    return A;
  }


  matrix &operator-=(const matrix B)
  {
    auto &A(*this);
    auto [m,n] = shape;
    
    for(int i=0;i<m;i++)
      for(int j=0;j<n;j++)
	A(i,j) -= B(i,j);

    return A;
  }  
  
  bool is_hermitian(real_t tolerance=0) const {
    auto [m,n] = shape;
    const auto &A(*this);

    if(m!=n) return false;
    
    for(int i=0;i<n;i++)
      for(int j=i;j<n;j++) if(std::abs(A(i,j) - ::conj(A(j,i)))>tolerance) return false;
    return true;
  }
  
  matrix diagonal(scalar *result_data) const
  {
    int n = std::max(shape[0],shape[1]);    
    const auto& A(*this);
    matrix<scalar> B(result_data,{n,n});

    if(shape[0] == 1){ // from vector to diagonal matrix
      memset(result_data,0,n*n*sizeof(scalar));
      for(int i=0;i<n;i++) B(i,i) = A(0,i);
    } else if(shape[1] == 1) {
      memset(result_data,0,n*n*sizeof(scalar));      
      for(int i=0;i<n;i++) B(i,i) = A(i,0); 
    } else {				// extract diagonal into vector
      B = matrix<scalar>(result_data,{1,n});
      for(int i=0;i<n;i++) B(0,i) = A(i,i);
    }

    return B;
  }
  
  
  INLINE void swap_column(int k, int p) {
    auto [m,n] = shape;
    auto &A(*this);

    for(int i=0;i<m;i++) std::swap(A(i,k), A(i,p));
  }

  INLINE void swap_row(int k, int p) {
    auto [m,n] = shape;
    auto &A(*this);

    for(int j=0;j<n;j++) std::swap(A(k,j), A(p,j));
  }  
  
  INLINE matrix operator()(array<int,2> I, array<int,2> J) const { return slice(I,J); }
  INLINE matrix operator()(array<int,2> I, int j) const    { return slice(I,{j,j+1}); }
  INLINE matrix operator()(int i, array<int,2> J) const    { return slice({i,i+1},J); }  


  INLINE ssize_t size() const { return shape[0]*shape[1]; }
 
  INLINE matrix& operator *=(const scalar s)
  {
    for(size_t i=0;i<size();i++) data[i] *= s;
    return *this;
  }

  template <typename T>
  INLINE matrix& operator +=(const matrix<T> &B)
  {
    auto &A(*this);
    auto [m,n] = A.shape;
    auto [p,q] = B.shape;
    assert(m==p);
    assert(n==q);
 
    for(int i=0;i<m;i++)
      for(int j=0;j<q;j++) A(i,j) += B(i,j);

    return *this;
  }

  template <typename T>
  INLINE matrix& operator -=(const matrix<T> &B)
  {
    auto &A(*this);
    auto [m,n] = A.shape;
    auto [p,q] = B.shape;
    assert(m==p);
    assert(n==q);
 
    for(int i=0;i<m;i++)
      for(int j=0;j<q;j++) A(i,j) -= B(i,j);

    return *this;
  }

  
  template <typename T>
  INLINE matrix& operator *=(const matrix<T> &B)
  {
    auto &A(*this);
    auto [m,n] = A.shape;
    auto [p,q] = B.shape;
    assert(n==p);
    assert(n==q);
 
   
    for(int i=0;i<m;i++){
      scalar row[q];
      for(int j=0;j<q;j++){
	scalar sum = 0;
	for(int k=0;k<n;k++) sum += A(i,k)*B(k,j);
	row[j] = sum;
      }
      for(int j=0;j<q;j++) A(i,j) = row[j];
    }
    return *this;
  }

  INLINE matrix slice(array<int,2> I, array<int,2> J) const {
    I[1] = std::min(I[1],shape[0]);
    J[1] = std::min(J[1],shape[1]);
    
    int m = I[1]-I[0], n = J[1]-J[0];
    int offset = I[0]*stride[0] + J[0]*stride[1];

    matrix result(data, {m,n}, stride, offset);    
    return result;
  }

  string str(real_t trunc_tolerance=100*machine_precision,
	     array<string,2> separators={",",",\n "},
	     array<string,2> brackets={"[","]"},
	     string null="Null") const
  {
    std::stringstream ss;
    const auto &A(*this);
    auto [m,n] = shape;
    auto [brack_open,brack_close] = brackets;

    if(trunc_tolerance < 0) trunc_tolerance = 100*machine_precision;
    real_t numerical_zero = trunc_tolerance*max_norm();
    
    if(data == 0) ss << null;
    else {
      ss << brack_open;
      for(int i=0;i<m;i++){
	ss << brack_open;
	for(int j=0;j<n;j++){
	  auto Aij = A(i,j);
	  ss << (std::abs(Aij)<=numerical_zero? 0 : Aij) << (j+1<n? separators[0] : "");
	}
	ss << brack_close + (i+1<m? separators[1] : "");
      } 
      ss << brack_close;
    }
    return ss.str();
  }

  operator string() const { return str(); }
};


typedef matrix<real_t>    real_matrix;
typedef matrix<complex_t> complex_matrix;



template <typename scalar>
matrix<scalar> identity(scalar *data, array<int,2> shape)
{
  matrix<scalar> A(data,shape);
  
  memset(data,0,shape[0]*shape[1]*sizeof(scalar));

  for(int i=0;i<std::min(shape[0],shape[1]);i++) A(i,i) = 1;

  return A;
}


template <typename scalar>
matrix<scalar>& mul(const matrix<scalar>& A, const matrix<scalar> &B, matrix<scalar> &C, const scalar scale=1){
  if(!(A.shape[1] == B.shape[0] &&
       C.shape[0] == A.shape[0] && C.shape[1] == B.shape[1])){
    fprintf(stderr,"Incompatible dimensions for matrix multiplication. Multiplicands: A.shape=(%d,%d), B.shape=(%d,%d); Result storage: C.shape=(%d,%d)\n",A.shape[0],A.shape[1],B.shape[0],B.shape[1], C.shape[0],C.shape[1]);
    abort();    
  }
    size_t m = A.shape[0], p = A.shape[1], n = B.shape[1];

    for(size_t i=0;i<m;i++)
      for(size_t j=0;j<n;j++){
	scalar sum=0;
	for(size_t k=0;k<p;k++) sum += scale*A(i,k)*B(k,j);
	C(i,j) = sum;
      }
  return C;    
}


template <typename scalar> matrix<real_t>& Re(const matrix<scalar> A, matrix<real_t> &C)
{
  auto [m,n] = A.shape;
  for(int i=0;i<m;i++)
    for(int j=0;j<n;j++)
      C(i,j) = Re(A(i,j));
  return C;
}

template <typename scalar> matrix<real_t>& Im(const matrix<scalar> A, matrix<real_t> &C)
{
  auto [m,n] = A.shape;
  for(int i=0;i<m;i++)
    for(int j=0;j<n;j++)
      C(i,j) = Im(A(i,j));
  return C;
}

// S*A*S.H()
template <typename S1, typename S2, typename scalar>
matrix<scalar> transform(const matrix<S1> A, const matrix<S2> S, matrix<scalar> &C)
{
  auto [m,n] = A.shape;
  auto [p,q] = C.shape;
  assert(m==n);
  assert(p==n && q==n);

  scalar SA[n*n];
  
  for(size_t i=0;i<n;i++)
    for(size_t j=0;j<m;j++){
      scalar sum = 0;
      for(size_t k=0;k<n;k++) sum += S(i,k)*A(k,j);
      SA[i*n+j] = sum;
    }

  for(size_t i=0;i<n;i++)
    for(size_t j=0;j<n;j++){
      scalar sum = 0;
      for(size_t k=0;k<n;k++) sum += SA[i*n+k]*conj(S(j,k));
      C(i,j) = sum;
    } 
  
  return C;
}






// NB: These all allocate new memore -- responsibility of caller to free it after use.
template <typename scalar>
matrix<scalar> operator-(const matrix<scalar>& A, const matrix<scalar>&B)
{
  int m = A.shape[0], p = A.shape[1], n = B.shape[1];
  
  scalar *C_data = (scalar*)calloc(m*n,sizeof(scalar));
  matrix<scalar> C = A.copy(C_data);
  C -= B;
  return C;
}

template <typename scalar>
matrix<scalar> operator*(const matrix<scalar>& A, const matrix<scalar>&B)
{
  int m = A.shape[0], p = A.shape[1], n = B.shape[1];
  
  scalar *C_data = (scalar*)calloc(m*n,sizeof(scalar));
  matrix<scalar> C(C_data,{m,n});
  mul(A,B,C);
  return C;
}



template <typename scalar>
INLINE scalar dot(const matrix<scalar>& a, const matrix<scalar> &b)
{
  // Check that a and b are both column vectors
  auto [m,n] = a.shape;
  auto [p,q] = b.shape;
  assert(n==q==1);
  assert(m==p);
  
  scalar sum = 0;
  
  for(size_t i=0;i<m;i++) sum += conj(a(i,0))*b(i,0);

  return sum;
}


//TODO: Explicit instantiation in matrix.cc
// extern template class matrix<real_t>;
// extern template class matrix<complex_t>;
