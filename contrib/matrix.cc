#include <cstdio>
#include <algorithm>
#include <vector>
#include <complex>
#include <iostream>
#include <cassert>
#include "ring.hh"


//Memory use: O(n^2), time: O(n^4)
int det_dense_clow_int(const std::vector< std::vector<int> >& A, int n) {
  std::vector< std::vector< std::vector<int> > > s(2, std::vector< std::vector<int> >(n, std::vector<int>(n, 0)));
    
  //init
  int t = n%2 == 0 ? -1 : 1;
  for(int c_0 = 0; c_0 < n; ++c_0) for(int c = 0; c < n; ++c) s[1][c_0][c] = c == c_0 ? t : 0;
    
  for(int l = 2; l <= n; ++l) {
    for(int c_0 = 0; c_0 < n; ++c_0) for(int c = 0; c < n; ++c) {
	s[l%2][c_0][c] = 0;
	if(c >= c_0) {
	  for(int c1 = 0; c1 < n; ++c1) {
	    s[l%2][c_0][c] += s[(l-1)%2][c_0][c1] * A[c1][c];
	  }
	  if(c_0 == c) {
	    for(int c1=0; c1 < n; ++c1) for(int c2=0; c2 < n; ++c2) {
		if(c1 >= c_0) {
		  s[l%2][c_0][c] -= s[(l-1)%2][c1][c2] * A[c2][c1];
		}
	      }
	  } 
	}
      }
  }
    
  int det = 0;
  for(int c_0 = 0; c_0 < n; ++c_0) for(int c = 0; c < n; ++c) {
      det += s[n%2][c_0][c] * A[c][c_0];
    } 
  return det; 
}


int det(const std::vector< std::vector<int> >& A, int n) {
	assert(A.size() == n && A[0].size() == n);
	//a[k][i][j] represents a^{(k)}_{ij}
	//and k = k mod 3
	std::vector< std::vector< std::vector<int> > > a(3, std::vector< std::vector<int> >(n, std::vector<int>(n, 0)));
	std::vector<bool> usedRow(n, false), usedCol(n, false);
	int i2, j2, i1 = -1, j1 = -1, i0 = -1, j0 = -1;
	for(int i = 0;i < n; ++i) for(int j = 0;j < n; ++j) {
		a[0][i][j] = A[i][j];
		if(i0 < 0 && a[0][i][j]) {
			i0 = i; j0 = j;
		}
	}
	usedRow[i0] = 1; usedCol[j0] = 1;
	for(int k = 1;k <= n-1; ++k) {
		i2 = i1; j2 = j1; //(i_2, j_2) = (i_1, j_1)
		i1 = i0; j1 = j0; //(i_1, j_1) = (i_0, j_0)
		i0 = -1; j0 = -1; //(i_0, j_0) = reset
		for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) if(!usedRow[i] && !usedCol[j]) {
			a[k%3][i][j] = a[(k-1)%3][i1][j1] * a[(k-1)%3][i][j] - a[(k-1)%3][i][j1] * a[(k-1)%3][i1][j];
			a[k%3][i][j] /= k > 1 ? a[(k-2)%3][i2][j2] : 1;
			if(i0 == -1 && a[k%3][i][j]) {
				i0 = i; j0 = j;
			}
		}
        if(i0 == -1) return 0; //Singular matrix
		usedRow[i0] = 1; usedCol[j0] = 1;
	}
	return a[(n-1)%3][i0][j0];
}
