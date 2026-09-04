#pragma once

#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>
#include <limits.h>
#include <math.h>
#include <assert.h>

#include "auxiliary.hh"
#include "fullerenes/dense_linalg_view.hh"   // LinAlg::matmul (operator* delegates)

using namespace std;

// Forward decls so matrix::APSP_with_paths can declare APSPResult<T> as
// its return type; APSPResult is defined after matrix because its
// members are matrix<T> / matrix<int> (which need the complete type).
template <typename T> class matrix;
template <typename T> struct APSPResult;

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
	for(int k=0;k<A.n;k++) x = std::min(x, T(A[i*A.n+k]+B[k*B.n+j]));
	x = std::min(x,infty_value);
	C[i*C.n+j] = x;
      }
    return C;    
  }

  matrix APSP(bool zero_diagonal=true) const {

    if(zero_diagonal){
      // When A(i,i) = 0, any path of length < m is included in the set
      // of paths of length m (as we can prefix with i->i->...->i).
      // Hence, we only need to calculate A^(2^n), with log2(m) matrix muls.
      int count = ceil(log2(m));
      matrix A(*this);
      for(int i=0;i<count;i++) A = minplus_multiply(A,A);

      return A;
    } else {
      // When we want non-trivial self-paths (A(i,i) > 0), we need all terms
      // in the Kleene sum, except for the identity.
      matrix A(*this), Ak(*this), S(*this);
      for(int i=0;i<m;i++){
	Ak = minplus_multiply(Ak,A);
	S  = min(S,Ak);
      }
      return S;
    }
  }

  // Floyd-Warshall with predecessor recording.  Reads (*this) as an edge-
  // weight matrix; entries >= `infty` are treated as "no direct edge".
  // Returns the all-pairs shortest-path distances and a "second-to-last"
  // predecessor matrix suitable for reconstruct_path().
  //
  // Self-distance is preserved as-is in dist; preds(u, u) is u when the
  // self-edge weight is finite, else -1.
  APSPResult<T> APSP_with_paths(const T& infty = std::numeric_limits<T>::max() / 2) const;

  // Element-wise sqrt; result type same as input. Naming this primitive
  // (rather than open-coding the loop) makes pipeline code read as math.
  matrix sqrt_elementwise() const {
    matrix R(*this);
    for (int i = 0; i < (int)this->size(); i++) R[i] = std::sqrt(R[i]);
    return R;
  }

  // Element-wise square (x -> x*x).
  matrix square_elementwise() const {
    matrix R(*this);
    for (int i = 0; i < (int)this->size(); i++) R[i] = R[i] * R[i];
    return R;
  }

  static matrix min(const matrix& A, const matrix& B)
  {
    matrix C(A.m,A.n);
    for(int i=0;i<A.size();i++) C[i] = std::min(A[i],B[i]);
    return C;
  }

  matrix operator*(const matrix& B) const
  {
    assert(n == B.m);
    matrix C(m,B.n);

    if constexpr (std::is_same_v<T, double>) {
      // @ref matmul-ijk-order (dense_linalg_view.hh): the double product IS
      // LinAlg::matmul -- ONE bit-pinned body shared with the solver family.
      // Do not reintroduce a local loop here.
      if (m > 0 && B.n > 0)
        LinAlg::matmul({std::span<const double>(this->data(), (size_t)m * n), m, n, n},
                       {std::span<const double>(B.data(), (size_t)B.m * B.n), B.m, B.n, B.n},
                       std::span<double>(C.data(), (size_t)m * B.n));
    } else {
      for(int i=0;i<m;i++)
        for(int j=0;j<B.n;j++){
          T x = 0;
          for(int k=0;k<n;k++) x += (*this)[i*n+k]*B[k*B.n+j];
          C(i,j) = x;
        }
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

  vector<vector<pair<int, T> > > sparse_representation() const {
    const matrix &A(*this);
    vector<vector<pair<int, T> > > result(m);

    for(int i=0;i<m;i++)
      for(int j=0;j<n;j++) if(A(i,j) != 0) result[i].push_back(make_pair(j,A(i,j)));

    return result;
  }
};


// All-pairs shortest-path result with predecessor recording. preds is
// indexed and interpreted by reconstruct_path() below.
template <typename T>
struct APSPResult {
  matrix<T>   dist;    // dist(u, v) = shortest u->v path length
  matrix<int> preds;   // preds(u, v) = second-to-last node K on shortest
                       //   u->v path (so the path reconstructs as
                       //   path(u, K) ++ [v]); preds(u, v) = u for a
                       //   direct edge; = -1 if v is unreachable from u.
};


template <typename T>
APSPResult<T> matrix<T>::APSP_with_paths(const T& infty) const {
  APSPResult<T> R{ matrix<T>(*this), matrix<int>(m, n, -1) };
  for (int u = 0; u < m; u++)
    for (int v = 0; v < n; v++)
      if (R.dist(u, v) < infty) R.preds(u, v) = u;
  for (int k = 0; k < n; k++)
    for (int u = 0; u < m; u++) {
      if (!(R.dist(u, k) < infty)) continue;
      for (int v = 0; v < n; v++) {
        if (!(R.dist(k, v) < infty)) continue;
        T via = R.dist(u, k) + R.dist(k, v);
        if (via < R.dist(u, v)) {
          R.dist(u, v)  = via;
          R.preds(u, v) = R.preds(k, v);
        }
      }
    }
  return R;
}


// Reconstructs the node sequence of the shortest U->V path from a
// "second-to-last" predecessor matrix (the second member of APSPResult).
//
//   preds(u, v) = K  means the shortest u->v path is u -...- K - v.
//   preds(u, v) = u  for a direct edge u - v.
//   preds(u, v) = -1 means v is unreachable from u.
//
// Returns:
//   [U]                           if U == V (trivial path).
//   [U, K_1, K_2, ..., V]         normal case.
//   {} (empty)                    if V is unreachable from U.
inline vector<int> reconstruct_path(const matrix<int>& preds, int U, int V) {
  if (U == V) return {U};
  if (preds(U, V) < 0) return {};
  vector<int> tail;
  int cur = V;
  while (cur != U) {
    tail.push_back(cur);
    int p = preds(U, cur);
    if (p < 0) return {};   // chain broken (defensive: shouldn't fire if preds is well-formed)
    cur = p;
  }
  tail.push_back(U);
  std::reverse(tail.begin(), tail.end());
  return tail;
}
