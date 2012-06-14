#ifndef FORTRAN_HH
# define FORTRAN_HH

extern "C" { 
  extern void windup_(const int* m,const int *ipr, int *ier, const int *s, int *d); 
  extern void unwind2_(const Graph **dual, const int *m, int *ier,
		       const int *s, int *nmr, char *group, const int length); 
  extern void dualanalyze2_(const int *n, const int *m, const Graph **dual, const Graph **g,
			   int *IRhag5,int *IRhag6,int *np,int *nelec, int *ndeg, double *sigmah, double *gap);
}


#endif
