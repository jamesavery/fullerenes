#ifndef MATRIX_HH
#define MATRIX_HH

#include "ring.hh"
#include <vector>
#include <algorithm>

int det(const std::vector< std::vector<int> >& A, int n);
int det_dense_clow_int(const std::vector< std::vector<int> >& A, int n);
template <class T>
int det_clow(const std::vector< std::vector<T> >& A, int n, ring<T> &r, bool sparse);
template <class T>
void sparse_representation(const std::vector< std::vector<T> >& A, 
                           std::vector< std::vector< std::pair<int, T> > >& sparse,
                           T zero);



//Returns a representation s.t. each column has a list.
template <class T>
void sparse_representation(const std::vector< std::vector<T> >& A, 
                           std::vector< std::vector< std::pair<int, T> > >& sparse,
                           T zero) {
    int N = (int) A.size();
    sparse.resize(N);
    std::fill(sparse.begin(), sparse.end(), std::vector< std::pair<int, T> >(0));
    for(int i = 0;i < N; ++i) for(int j = 0;j < N; ++j) {
        if(A[i][j] != zero) {
            sparse[j].push_back(std::make_pair(i, A[i][j]));
        }
    }
}

//sparse=false: Time: O(n^4), space: O(n^2)
//sparse=true:  Time: O(n^2*(n+m)), space: O(n^2), where m is the number of non-zero entries.
template <class T>
int det_clow(const std::vector< std::vector<T> >& A, int n, ring<T> &r, bool sparse) {
  std::vector< std::vector< std::vector<T> > > s(2, std::vector< std::vector<T> >(n, std::vector<T>(n,r.add_id())));
    
    //init
  T t = n%2 == 0 ? r.add_inv(r.mul_id()) : r.mul_id();
  for(int c_0 = 0; c_0 < n; ++c_0) for(int c = 0; c < n; ++c) s[1][c_0][c] = c == c_0 ? t : r.add_id();
  
  if(sparse) {
    std::vector< std::vector< std::pair<int, T> > > B;
    sparse_representation<T>(A, B, r.add_id());
    for(int l = 2; l <= n; ++l) {
      for(int c_0 = 0; c_0 < n; ++c_0) for(int c = 0; c < n; ++c) {
	  s[l%2][c_0][c] = r.add_id();
	  if(c >= c_0) {
	    for(size_t i = 0; i < B[c].size(); ++i) {
	      int c1 = B[c][i].first;
	      T v    = B[c][i].second; 
	      s[l%2][c_0][c] = r.add(s[l%2][c_0][c], r.mul(s[(l-1)%2][c_0][c1], v));
	    }
	    if(c_0 == c) {
	      for(int c1=0; c1 < n; ++c1) {
		for(size_t i = 0; i < B[c1].size(); ++i) {
		  int c2 = B[c1][i].first;
		  T v    = B[c1][i].second;
		  if(c1 >= c_0) {
		    s[l%2][c_0][c] = r.sub(s[l%2][c_0][c], r.mul(s[(l-1)%2][c1][c2], v));
		  }
		}
	      }
	    } 
	  }
	}
    }
  } else {
    for(int l = 2; l <= n; ++l) {
      for(int c_0 = 0; c_0 < n; ++c_0) for(int c = 0; c < n; ++c) {
	  s[l%2][c_0][c] = r.add_id();
	  if(c >= c_0) {
	    for(int c1 = 0; c1 < n; ++c1) {
	      s[l%2][c_0][c] = r.add(s[l%2][c_0][c], r.mul(s[(l-1)%2][c_0][c1], A[c1][c]));
	    }
	    if(c_0 == c) {
	      for(int c1=0; c1 < n; ++c1) for(int c2=0; c2 < n; ++c2) {
		  if(c1 >= c_0) {
		    s[l%2][c_0][c] = r.sub(s[l%2][c_0][c], r.mul(s[(l-1)%2][c1][c2], A[c2][c1]));
		  }
		}
	    } 
	  }
	}
    }
  }
  T det = r.add_id();
  for(int c_0 = 0; c_0 < n; ++c_0) for(int c = 0; c < n; ++c) {
      det = r.add(det, r.mul(s[n%2][c_0][c], A[c][c_0]));
    } 
  return det; 
}



#endif
