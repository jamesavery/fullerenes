#ifndef MATRIX_HH
#define MATRIX_HH

#include <vector>
#include <iostream>
#include <climits>
#include <cassert>
#include <cmath>

using namespace std;

template <typename T> class matrix : public vector<T> {
public:
  int m,n; 

  matrix(int m, int n, const vector<T>& A) : vector<T>(A.begin(), A.end()), m(m), n(n) { assert(A.size()==m*n); }
  matrix(int m, int n, const T& zero = 0) : vector<T>(m*n,zero), m(m), n(n) {}

  // Convert from compatible type
  template <typename S> matrix(const matrix<S>& A) : vector<T>(A.begin(),A.end()), m(A.m), n(A.n) {}

  T& operator()       (int i, int j){ return (*this)[i*n+j]; }
  T  operator()(int i, int j) const { return (*this)[i*n+j]; }

  static matrix minplus_multiply(const matrix& A, const matrix& B, const T& infty_value = USHRT_MAX)
  {
    assert(A.n == B.m);
    matrix C(A.m,B.n);

    for(int i=0;i<A.m;i++)
      for(int j=0;j<B.n;j++){
	T x = infty_value;
	for(int k=0;k<A.n;k++) x = min(x, T(A[i*A.n+k]+B[k*B.n+j]));
	x = min(x,infty_value);
	C[i*C.n+j] = x;
      }
    return C;    
  }

  matrix APSP() const {
    int count = ceil(log2(m));
    matrix A(*this);
    for(int i=0;i<count;i++) A = minplus_multiply(A,A);
    
    return A;
  }

  matrix operator*(const matrix& B)
  {
    assert(n == B.m);
    matrix C(m,B.n);

    for(int i=0;i<m;i++)
      for(int j=0;j<B.n;j++){
	T x = 0;
	for(int k=0;k<n;k++) x += (*this)[i*n+k]*B[k*B.n+j];
	C(i,j) = x;
      }
    return C;    
  }

  matrix operator+(const matrix& B){ 
    assert(n == B.n && m == B.m);
    matrix C(*this);
    for(int i=0;i<C.size();i++) C[i] += B[i];
	return C;
  }

  friend ostream& operator<<(ostream& S, const matrix& A)
  {
    vector< vector<T> > VV(A.m, vector<T>(A.n));
    for(int i=0;i<A.m;i++) 
      for(int j=0;j<A.n;j++)
	VV[i][j] = A[i*A.n+j];

    S << VV;
    return S;
  }
};

#endif
