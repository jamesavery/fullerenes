#include <cstdio>
#include <algorithm>
#include <vector>
#include <complex>
#include <iostream>
#include <cassert>

int det(std::vector< std::vector<int> > A, int n) {
	assert(A.size() == n && A[0].size() == n);
	//a[k][i][j] represents a^{(k)}_{ij}
	//and k = k mod 3
	std::vector< std::vector< std::vector<int> > > a(3, std::vector< std::vector<int> >(n, std::vector<int>(n, 0)));
	std::vector<bool> usedRow(n, false), usedCol(n, false);
	int i2, j2, i1, j1, i0 = -1, j0 = -1;
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
