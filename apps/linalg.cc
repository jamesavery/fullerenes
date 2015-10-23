#include <stdlib.h>
#include <vector>

using namespace std;

extern "C" {
#if defined LAPACK_IS_MKL
# include <mkl_lapack.h>
#else 
void dsyevx_(char *jobz, char *range,char *uplo,int *n, double *a, int *lda,  double *vl,   
	    double *vu, int *il, int  *iu, double *abstol,  int *m,  double *w,  double *z,  int *ldz, 
	     double *work, int *lwork, int *iwork, int *ifail,  int *info);
#endif

  void dsyevx(const char jobz, const char range,const char uplo,int n, double *a, int lda,  double vl,   
	      double vu, int il, int  iu, double abstol,  int *m,  double *w,  double *z,  int ldz, 
	      int *ifail,  int *info)
  {
    char JOBZ = jobz, RANGE=range,UPLO=uplo;
    int lwork = 10*n;
    double *work = (double*)malloc(10*n*sizeof(double));
    int   *iwork = (int*)malloc(5*n*sizeof(int));

    dsyevx_(&JOBZ,&RANGE,&UPLO,&n,a,&lda, &vl,&vu,&il,&iu,&abstol,m,w,z,&ldz,work,&lwork,iwork,ifail,info);


    free(work);
    free(iwork);
  }
}

class DenseSqrMatrix: vector<double> {
public:
  int n;

  DenseSqrMatrix(int n=1) : vector<double>(n*n),n(n) {}

  vector<double> eigenvalues(int start=0, int end=-1){
    if(end<0) end = n-1;
    const int nlambda = end-start+1;
    int M, INFO;

    vector<double> lambda(nlambda), ap(n*n);
    vector<int> ifail(n);

    for(int i=0;i<n*n;i++) ap[i] = (*this)[i];

    dsyevx('N',end==n-1?'A':'I','L',n,&ap[0],n,-1,-1,start+1,end+1,2*1e-200,&M,&lambda[0],NULL/*Z*/,n,&ifail[0],&INFO);

    return lambda;
  }

  double& operator()(int i, int j){ return (*this)[i*n+j]; }
  double  operator() (int i, int j) const { return (*this)[i*n+j]; }
};



